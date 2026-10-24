#!/usr/bin/env bash
# ============================================================================
# 单样本 velocyto：阶段边界闸→rsync 回 BAM→校验比对成功→run10x/run→验证 loom→删中间BAM→删本地BAM
# 用法: 02_velocyto.sh <sample_id> <samplefolder> <bam_archive_dir>
#   cellranger 分支（ALIGNER=cellranger）：velocyto run10x，输入 cellranger 样本目录
#   mobivision 分支（ALIGNER=mobivision）：velocyto run -b -e，输入 mobivision outs 目录
#     （-e <sid> 保证 loom 命名为 {sid}.loom；交付盘 BAM 会 rsync 回本地）
# 依赖 run_config.sh（ALIGNER/VELO_REF/CONDA_SH/ENV_VELOCYTO/SAMPLE_SHEET/OUTPUT_ROOT/
#                    PROJECT_NAME；VELO_THREADS 为可选线程覆盖位，默认空=不显式传 -@；
#                    VELO_PER_SAMPLE_HANDOFF=1 显式豁免阶段边界闸，默认空=强制）
# 调度：默认样本表全部样本 quantify 完成才放行（开头闸口逐样本核验）；放行后可 4 并发
#       + 15 分钟发车间隔（外层编排，见 SKILL.md 关键阶段边界）
# ============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/run_config.sh"

SID="$1"; SAMPLEDIR="$2"; ZBAMDIR="$3"

# 0a. 阶段边界闸（默认强制）：样本表全部样本 quantify 完成标记确认后才放行 velocyto。
#     动机：quantify 高峰期独占整机线程/内存（mobivision 另有数百 GB 级中间产物磁盘压力），
#     与 velocyto 重叠既拖垮在跑的 quantify，也破坏阶段1 串行铁律的隔离前提——
#     阶段边界不是建议，未获显式豁免前不得越界发车（见 SKILL.md 关键阶段边界）。
#     标记语义（逐样本；本地工作目录为主，样本表 bam_archive_dir 列/第三参数指向目录为回退，
#     覆盖两种既有布局——交付根 `<dir>/` 与已含样本层的 `<dir>/<sid>`；marker 再找 outs/ 下）：
#       mobivision：_flagdone（本地 $OUTPUT_ROOT/$PROJECT_NAME/<sid>/；回退交付盘）
#       cellranger：outs/web_summary.html（count 成功产物；回退归档目录）
#       交付盘候选目录 = bam_archive_dir 列（非空时）∪ 第三参数（非空时），逐候选测
#       <cand>/<sid>/、<cand>/ 两种层级，marker 路径 <层>/_flagdone（mobivision）
#       或 <层>/outs/web_summary.html（cellranger）。
#     豁免：用户在 run_config 显式填 VELO_PER_SAMPLE_HANDOFF=1（声明"每样本 quantify
#           完成即可直接接续该样本 velocyto"）——未声明一律拒跑并列出缺项。
if [ -z "${VELO_PER_SAMPLE_HANDOFF:-}" ]; then
  echo "==== [$(date '+%F %T')] 阶段边界闸：核验样本表全部样本 quantify 完成 ===="
  [ -f "$SAMPLE_SHEET" ] || { echo "[err] SAMPLE_SHEET 不存在，无法核验阶段边界: $SAMPLE_SHEET"; exit 1; }
  mapfile -t ALL_SIDS < <(sed '1s/^﻿//' "$SAMPLE_SHEET" | tr -d '\r' | awk -F, 'NR>1 && $1!~/^#/ && NF>=2 {print $1}')
  [ "${#ALL_SIDS[@]}" -gt 0 ] || { echo "[err] 样本表无有效样本行，无法核验阶段边界"; exit 1; }
  # col4（bam_archive_dir）提列用 awk：read（无论 -a 还是逐变量）在 IFS=, 下丢/并不空字段，
  # 无法保证第 4 列定位。awk -F, 的 $4 与 NF>=2 样本集同行同序（同一文件、同过滤条件）。
  declare -A ARC_COL=()
  mapfile -t _arc_lines < <(sed '1s/^﻿//' "$SAMPLE_SHEET" | tr -d '\r' | awk -F, 'NR>1 && $1!~/^#/ && NF>=2 {print $4}')
  for i in "${!ALL_SIDS[@]}"; do
    ARC_COL["${ALL_SIDS[$i]}"]="${_arc_lines[$i]:-}"
  done
  NOT_DONE=()
  if [ "$ALIGNER" = "mobivision" ]; then MARKER=_flagdone; else MARKER=outs/web_summary.html; fi
  for s in "${ALL_SIDS[@]}"; do
    ok=0
    [ -f "$OUTPUT_ROOT/$PROJECT_NAME/$s/$MARKER" ] && ok=1
    if [ "$ok" = 0 ]; then
      # 回退：候选目录（col4 非空 ∪ 第三参数非空）× 层级（<cand>/<sid>、<cand>）
      for cand in "${ARC_COL[$s]:-}" "$ZBAMDIR"; do
        [ "$ok" = 1 ] && break
        [ -n "$cand" ] || continue
        [ -f "$cand/$s/$MARKER" ] && ok=1
        [ -f "$cand/$MARKER" ] && ok=1
      done
    fi
    [ "$ok" = 1 ] || NOT_DONE+=("$s")
  done
  if [ "${#NOT_DONE[@]}" -gt 0 ]; then
    echo "[err] 阶段边界闸拒绝 velocyto $SID：以下样本 quantify 完成标记未确认："
    for s in "${NOT_DONE[@]}"; do echo "        - $s"; done
    echo "      须样本表全部样本 quantify 完成后才可发车（防止与在跑 quantify 争抢资源/磁盘）。"
    echo "      如用户显式声明每样本可直接接续，在 run_config 置 VELO_PER_SAMPLE_HANDOFF=1。"
    exit 1
  fi
  echo "  通过：${#ALL_SIDS[@]} 个样本 quantify 完成标记全部确认"
else
  echo "[warn] VELO_PER_SAMPLE_HANDOFF=1：阶段边界闸已按用户显式声明豁免，仅跑本样本完成校验"
fi

RMSK="$VELO_REF/rmsk.gtf"
GENES="$VELO_REF/genes.gtf"
OUTS="$SAMPLEDIR/outs"

# 0. VELO_REF 必须与比对阶段使用同一参考 GTF；无法确认时中止，不以基因数相等替代同源性校验
[ -f "$GENES" ] || { echo "[err] VELO_REF 缺少 genes.gtf: $GENES"; exit 1; }
if [ "$ALIGNER" = "mobivision" ]; then
  [ -n "${MOBI_INDEX:-}" ] || { echo "[err] MOBI_INDEX 未配置"; exit 1; }
  [ -f "$MOBI_INDEX/genes/genes.gtf" ] || { echo "[err] MOBI_INDEX 缺少 genes/genes.gtf"; exit 1; }
  ref_hash=$(sha256sum "$MOBI_INDEX/genes/genes.gtf" | cut -d' ' -f1)
  velo_hash=$(sha256sum "$GENES" | cut -d' ' -f1)
  [ "$ref_hash" = "$velo_hash" ] || { echo "[err] VELO_REF/genes.gtf 与 MOBI_INDEX/genes/genes.gtf 不一致；必须使用同一参考 GTF"; exit 1; }
elif [ "$ALIGNER" = "cellranger" ]; then
  [ -n "${TRANSCRIPTOME:-}" ] || { echo "[err] TRANSCRIPTOME 未配置"; exit 1; }
  [ -f "$TRANSCRIPTOME/genes/genes.gtf" ] || { echo "[err] TRANSCRIPTOME 缺少 genes/genes.gtf，无法校验参考一致性"; exit 1; }
  ref_hash=$(sha256sum "$TRANSCRIPTOME/genes/genes.gtf" | cut -d' ' -f1)
  velo_hash=$(sha256sum "$GENES" | cut -d' ' -f1)
  [ "$ref_hash" = "$velo_hash" ] || { echo "[err] VELO_REF/genes.gtf 与 TRANSCRIPTOME/genes/genes.gtf 不一致；必须使用同一参考 GTF"; exit 1; }
else
  echo "[err] ALIGNER 仅支持 cellranger|mobivision，当前: $ALIGNER"; exit 1
fi

# 线程默认由 velocyto 决定，不显式传 -@/线程参数；VELO_THREADS 非空时为用户显式覆盖。
# 大文件读取由外层调度保持串行/间隔。
VELO_T=()
[ -n "${VELO_THREADS:-}" ] && VELO_T=(-@"$VELO_THREADS")

echo "==== [$(date '+%F %T')] velocyto $SID 开始（ALIGNER=$ALIGNER）===="

if [ "$ALIGNER" = "mobivision" ]; then
  # ---- mobivision 分支：BAM 在 ${SID}_outs/ 下，velocyto run -b -e ----
  M_OUTS="$SAMPLEDIR/${SID}_outs"
  BAM="$M_OUTS/${SID}_Aligned.sort.bam"
  BARCODES="$M_OUTS/filtered_cell_gene_matrix/barcodes.tsv.gz"
  CELLBAM="$M_OUTS/cellsorted_${SID}_Aligned.sort.bam"

  # 1. BAM 可能已随交付清空，从归档/交付盘 rsync 回本地 outs
  if [ ! -f "$BAM" ]; then
    [ -n "$ZBAMDIR" ] || { echo "[err] BAM 缺失且未给归档目录"; exit 1; }
    rsync -a --exclude='._*' --include="${SID}_Aligned.sort.bam" --exclude='*' "$ZBAMDIR/" "$M_OUTS/"
  fi
  [ -f "$BAM" ] || { echo "[err] rsync 后 BAM 仍缺失"; exit 1; }

  # 2. 校验 mobivision 完成标记
  [ -f "$SAMPLEDIR/_flagdone" ] || { echo "[err] 无 _flagdone（quantify 未完成）"; exit 1; }
  [ -f "$BARCODES" ] || { echo "[err] 缺 filtered barcodes"; exit 1; }

  # 3. velocyto run（-e <sid> 命名 loom；-b 显式 barcodes；BAM index 缺失仅是警告）
  source "$CONDA_SH"; conda activate "$ENV_VELOCYTO"
  mkdir -p "$SAMPLEDIR/velocyto"
  velocyto run \
    -b "$BARCODES" \
    -o "$SAMPLEDIR/velocyto" \
    -e "$SID" \
    -m "$RMSK" \
    "${VELO_T[@]}" \
    "$BAM" \
    "$GENES"

  # 4. 验证 loom（必须 {sid}.loom 命名）
  LOOM="$SAMPLEDIR/velocyto/$SID.loom"
  [ -f "$LOOM" ] && [ -s "$LOOM" ] || { echo "[err] loom 未生成或为空（检查是否漏传 -e $SID）"; exit 1; }
  echo "  loom $(stat -c%s "$LOOM") B"

  # 5-6. 删除中间 cellsorted BAM 与本地 BAM 副本（交付盘原件保留）
  rm -f "$CELLBAM" "$BAM"
else
  # ---- cellranger 分支（原逻辑）----
  BAM="$OUTS/possorted_genome_bam.bam"
  CELLBAM="$OUTS/cellsorted_possorted_genome_bam.bam"

  # 1. rsync 原 BAM 回本地 outs（velocyto 需读 BAM；直读挂载盘极慢）
  if [ ! -f "$BAM" ]; then
    rsync -a --exclude='._*' --include='possorted_genome_bam.bam' --exclude='*' "$ZBAMDIR/" "$OUTS/"
  fi
  [ -f "$BAM" ] || { echo "[err] rsync 后 BAM 仍缺失"; exit 1; }

  # 2. 校验 cellranger 完成标记
  if [ -f "$SAMPLEDIR/_log" ]; then
    grep -q "Pipestance completed successfully!" "$SAMPLEDIR/_log" || { echo "[err] cellranger _log 无成功标记"; exit 1; }
  else
    # _log 可能已在阶段1 清理中移除；以官方的 web_summary 作为等价的完成标记
    [ -f "$OUTS/web_summary.html" ] || { echo "[err] 缺 _log 且无 outs/web_summary.html，无法确认比对完成"; exit 1; }
  fi
  [ -f "$OUTS/filtered_feature_bc_matrix/barcodes.tsv.gz" ] || { echo "[err] 缺 filtered barcodes"; exit 1; }

  # 3. velocyto run10x（-m mask 重复区；loom 产出到 SAMPLEDIR/velocyto/）
  source "$CONDA_SH"; conda activate "$ENV_VELOCYTO"
  velocyto run10x -m "$RMSK" "${VELO_T[@]}" "$SAMPLEDIR" "$GENES"

  # 4. 验证 loom
  LOOM="$SAMPLEDIR/velocyto/$SID.loom"
  [ -f "$LOOM" ] && [ -s "$LOOM" ] || { echo "[err] loom 未生成或为空"; exit 1; }
  echo "  loom $(stat -c%s "$LOOM") B"

  # 5-6. 删除中间 BAM 与本地原 BAM 副本（Z 盘原件保留）
  rm -f "$CELLBAM" "$BAM"
fi

echo "==== [$(date '+%F %T')] velocyto $SID 完成 ===="
