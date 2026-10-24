#!/usr/bin/env bash
# ============================================================================
# 单样本 cellranger count：rsync fastq→md5→重命名前缀→count→验证→BAM归档→删暂存
# 用法: 01_cellranger_count.sh <sample_id> <fastq_source_dir> <fastq_prefix> <bam_archive_dir> <cell_or_nuclei>
# 依赖 run_config.sh（PROJECT_NAME/TRANSCRIPTOME/FASTQ_STAGING/OUTPUT_ROOT；
#                    LOCALCORES/LOCALMEM 为可选覆盖位，默认空=不显式传线程/内存参数）
# ============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/run_config.sh"

SID="$1"; SRC="$2"; PREFIX="$3"; BAMDIR="$4"; CELL="$5"

STAGE="$FASTQ_STAGING/$PROJECT_NAME/$SID"
OUT_ROOT="$OUTPUT_ROOT/$PROJECT_NAME"
TRANSCRIPT="$TRANSCRIPTOME"
CR="cellranger"; [ -n "$CELLRANGER_DIR" ] && CR="$CELLRANGER_DIR/cellranger"

echo "==== [$(date '+%F %T')] cellranger $SID 开始 ===="

# 0. 残留检查（严格流式前提：staging 与样本输出目录必须为空，防中断重跑/并行混样）
[ -d "$STAGE" ] && [ "$(ls -A "$STAGE" 2>/dev/null | wc -l)" -gt 0 ] && { echo "[err] 残留: $STAGE 非空（上一个样本未清空？）"; exit 1; }
[ -d "$OUT_ROOT/$SID" ] && [ "$(ls -A "$OUT_ROOT/$SID" 2>/dev/null | wc -l)" -gt 0 ] && { echo "[err] 残留: $OUT_ROOT/$SID 非空"; exit 1; }

mkdir -p "$STAGE" "$OUT_ROOT/$SID/outs"

# 1. rsync fastq 到本地暂存（排除 macOS 资源叉 ._ 与 md5）
rsync -a --exclude='._*' --include='*.fastq.gz' --exclude='*' "$SRC/" "$STAGE/"

# 2. 重命名前缀为 sample_id，使 --sample 命中
for suf in R1 R2; do
  old="$STAGE/${PREFIX}_S1_L001_${suf}_001.fastq.gz"
  new="$STAGE/${SID}_S1_L001_${suf}_001.fastq.gz"
  [ -f "$old" ] || { echo "[err] 缺失 $old"; exit 1; }
  mv "$old" "$new"
done

# 3. md5 校验（用源目录原始 .md5 的 hash 比对）
for suf in R1 R2; do
  md5file="$SRC/${PREFIX}_S1_L001_${suf}_001.fastq.gz.md5"
  if [ -f "$md5file" ]; then
    exp=$(awk '{print $1}' "$md5file")
    act=$(md5sum "$STAGE/${SID}_S1_L001_${suf}_001.fastq.gz" | awk '{print $1}')
    [ "$exp" = "$act" ] || { echo "[err] md5 校验失败 $SID $suf: $exp != $act"; exit 1; }
    echo "  md5 $suf OK"
  else
    echo "  [warn] 无 $md5file，跳过校验"
  fi
done

# 4. cellranger count（--nosecondary 跳过二级分析，--create-bam 保留 BAM 供 velocyto）
#    线程/内存默认由软件与系统级资源配置决定；LOCALCORES/LOCALMEM 非空时为用户显式覆盖
CR_RES=()
[ -n "${LOCALCORES:-}" ] && CR_RES+=(--localcores="$LOCALCORES")
[ -n "${LOCALMEM:-}" ] && CR_RES+=(--localmem="$LOCALMEM")
cd "$OUT_ROOT"
$CR count \
  --id="$SID" \
  --transcriptome="$TRANSCRIPT" \
  --fastqs="$STAGE" \
  --sample="$SID" \
  --nosecondary \
  --create-bam=true \
  "${CR_RES[@]}"

# 5. 验证产出
OUTS="$OUT_ROOT/$SID/outs"
for item in filtered_feature_bc_matrix filtered_feature_bc_matrix.h5 raw_feature_bc_matrix.h5 \
            possorted_genome_bam.bam possorted_genome_bam.bam.bai metrics_summary.csv web_summary.html; do
  [ -e "$OUTS/$item" ] || { echo "[err] 缺失产出 $OUTS/$item"; exit 1; }
done
echo "  产出验证通过"

# 6. BAM 归档到源目录（velocyto 阶段会 rsync 回本地；大文件用 rsync 不用 cp——
#    cp 写 drvfs/9p 挂载盘大文件会报"无法分配内存"）
if [ -n "$BAMDIR" ]; then
  rsync -a --append-verify "$OUTS/possorted_genome_bam.bam" "$BAMDIR/"
  rsync -a --append-verify "$OUTS/possorted_genome_bam.bam.bai" "$BAMDIR/"
  local_sz=$(stat -c%s "$OUTS/possorted_genome_bam.bam")
  zbam_sz=$(stat -c%s "$BAMDIR/possorted_genome_bam.bam")
  [ "$local_sz" = "$zbam_sz" ] || { echo "[err] BAM 归档大小不一致"; exit 1; }
  rm -f "$OUTS/possorted_genome_bam.bam" "$OUTS/possorted_genome_bam.bam.bai"
  echo "  BAM 已归档至 $BAMDIR，本地副本已删"
fi

# 7. 清理 cellranger 临时文件（白名单：仅保留 outs 与 velocyto）
#    删除 SC_RNA_COUNTER_CS/、extras/、.mri.tgz 及所有 _* 元数据，释放空间
SAMPLEDIR="$OUT_ROOT/$SID"
if [ -d "$SAMPLEDIR/outs" ]; then
  echo ">>> 清理 cellranger 临时文件（保留 outs, velocyto）"
  removed=0
  shopt -s dotglob nullglob
  for item in "$SAMPLEDIR"/*; do
    case "$(basename "$item")" in outs|velocyto|_log) continue ;; esac
    echo "  删除 $(basename "$item")"
    rm -rf "$item"
    removed=$((removed+1))
  done
  shopt -u dotglob nullglob
  echo "  共清理 $removed 项"
else
  echo "[warn] $SAMPLEDIR/outs 不存在，跳过清理以免误删"
fi

# 8. 删除本地 fastq 暂存
rm -rf "$STAGE"
echo "==== [$(date '+%F %T')] cellranger $SID 完成 ===="