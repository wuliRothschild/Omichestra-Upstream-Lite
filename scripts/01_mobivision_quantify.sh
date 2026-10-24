#!/usr/bin/env bash
# ============================================================================
# 单样本 mobivision quantify：平铺复制 fastq→md5→quantify（严格串行）→
#   阶段监听清理中间产物→验证→rsync 交付→清理
# 用法: 01_mobivision_quantify.sh <sample_id> <fastq_source_dir> <delivery_dir>
# 依赖 run_config.sh（ALIGNER=mobivision 时必须；MOBIVISION_DIR/MOBI_INDEX/
#                       FASTQ_STAGING/OUTPUT_ROOT；LOCALCORES 为可选线程覆盖位，默认空=不显式传 -t）
# 铁律（详见 references/mobivision-platform.md §3-4）：
#   - quantify 严格串行：上一个样本 _flagdone 且全部文件清空后才能跑本脚本
#   - 中间产物按阶段监听删除（convert_fastq 体积可达数百 GB 以上，是磁盘杀手）
#   - 交付用 rsync（cp 写挂载盘大文件会报内存错误）
# ============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/run_config.sh"

SID="$1"; SRC="$2"; DELIVERY="$3"
[ "$ALIGNER" = "mobivision" ] || { echo "[err] 本脚本仅用于 ALIGNER=mobivision"; exit 1; }

STAGE="$FASTQ_STAGING/$PROJECT_NAME/$SID"
OUT_ROOT="$OUTPUT_ROOT/$PROJECT_NAME/$SID"
MOBI_LOG="$OUT_ROOT/logs/mobivision.log"

echo "==== [$(date '+%F %T')] mobivision quantify $SID 开始 ===="

# 0. 残留检查（严格流式前提：staging 与输出目录必须为空）
[ -d "$STAGE" ] && [ "$(ls -A "$STAGE" 2>/dev/null | wc -l)" -gt 0 ] && { echo "[err] 残留: $STAGE 非空（上一个样本未清空？）"; exit 1; }
[ -d "$OUT_ROOT" ] && [ "$(ls -A "$OUT_ROOT" 2>/dev/null | wc -l)" -gt 0 ] && { echo "[err] 残留: $OUT_ROOT 非空"; exit 1; }

# 1. rsync 递归复制 fastq + md5 到本地暂存（墨卓样本目录可能含 S 编号子目录；
#    只读源，绝不 mv/删源文件——原始数据只允许复制）
mkdir -p "$STAGE" "$OUT_ROOT"
rsync -a --exclude='._*' --include='*.fastq.gz' --include='*/' --exclude='*' "$SRC/" "$STAGE/"
rsync -a --exclude='._*' --include='*md5*' --exclude='*' "$SRC/" "$STAGE/" || true
# staging 内平铺：子目录中的 fastq 移到暂存根（只动本地副本，不碰源；mobivision 用 os.listdir 不递归子目录）
while IFS= read -r -d '' f; do
  base="$(basename "$f")"
  [ -e "$STAGE/$base" ] && { echo "[err] 平铺文件名冲突: $base"; exit 1; }
  mv "$f" "$STAGE/$base"
done < <(find "$STAGE" -mindepth 2 -name '*.fastq.gz' -print0 2>/dev/null || true)
find "$STAGE" -mindepth 1 -type d -empty -delete 2>/dev/null || true

NFQ=$(find "$STAGE" -name '*.fastq.gz' | wc -l)
[ "$NFQ" -ge 2 ] || { echo "[err] 暂存目录 fastq 不足: $STAGE"; exit 1; }
echo "  fastq 平铺完成: $NFQ 个文件 -> $STAGE"

# 2. md5 校验（墨卓交付的 {sid}_Filename_only_md5.xls 为 md5sum -c 兼容平铺格式；
#    xls 可能含已删的其他样本条目，--ignore-missing 跳过缺失项，现存文件严格校验）
MD5XLS=$(find "$STAGE" -name '*_Filename_only_md5.xls' | head -1)
if [ -n "$MD5XLS" ]; then
  ( cd "$STAGE" && md5sum -c --ignore-missing "$(basename "$MD5XLS")" 2>/dev/null ) || { echo "[err] md5 校验失败"; exit 1; }
  echo "  md5 校验通过"
else
  echo "  [warn] 未找到 _Filename_only_md5.xls，跳过校验"
fi

# 3. quantify（后台运行 + 阶段监听清理）
#    线程/内存默认由软件决定，不显式传 -t（显式限制曾致软件崩溃；系统级配置负责资源上限）；
#    LOCALCORES 非空时为用户显式覆盖
MOBI_T=()
[ -n "${LOCALCORES:-}" ] && MOBI_T=(-t "$LOCALCORES")
source "$MOBIVISION_DIR/source.sh"
cd "$OUT_ROOT"
mobivision quantify \
  -i "$MOBI_INDEX" \
  -f "$STAGE" \
  -o "$OUT_ROOT" \
  -s "$SID" \
  "${MOBI_T[@]}" \
  --intron included &
QPID=$!

# 4. 阶段监听清理（核心经验脚本化：盯 mobivision.log 阶段关键字 + 磁盘告警）
echo ">>> 阶段监听清理启动（PID $QPID）"
src_cleaned=0; mid_cleaned=0
while kill -0 "$QPID" 2>/dev/null; do
  # merge 完成（ANALYSIS STARTED = gzip 校验开始，merged 已就绪）→ 删源 fastq
  if [ "$src_cleaned" -eq 0 ] && grep -q "MOBIVISION v4.1 ANALYSIS STARTED" "$MOBI_LOG" 2>/dev/null; then
    rm -f "$STAGE"/*.fastq.gz "$STAGE"/*md5*
    echo "  [$(date '+%T')] merge 完成，已删源 fastq（释放 $(du -sh "$STAGE" 2>/dev/null | cut -f1) 占用）"
    src_cleaned=1
  fi
  # reads-trimming 完成后删除 convert/split/merged（cutadapt 在此阶段前仍读取 convert_fastq）
  if [ "$mid_cleaned" -eq 0 ] && grep -q "READS-TRIMMING FINISHED" "$MOBI_LOG" 2>/dev/null; then
    rm -rf "$OUT_ROOT/convert_fastq" "$OUT_ROOT/split_fastq"
    rm -rf "$STAGE/merged"
    echo "  [$(date '+%T')] READS-TRIMMING FINISHED，已删 convert/split/merged"
    mid_cleaned=1
  fi
  # 磁盘告警
  avail_gb=$(df -BG / | awk 'NR==2 {gsub(/G/,"",$4); print $4}')
  if [ "${avail_gb:-0}" -lt 50 ] 2>/dev/null; then
    echo "  [warn] 磁盘可用仅 ${avail_gb}G，注意 convert_fastq/map_result 增长"
  fi
  sleep 300
done
wait "$QPID" || { echo "[err] mobivision quantify 失败，见 $MOBI_LOG"; exit 1; }
[ -f "$OUT_ROOT/_flagdone" ] || { echo "[err] 无 _flagdone，quantify 未正常完成"; exit 1; }

# 5. 验证产出
OUTS="$OUT_ROOT/${SID}_outs"
for item in "${SID}_Aligned.sort.bam" "${SID}_Aligned.sort.bam.bai" \
            filtered_cell_gene_matrix raw_cell_gene_matrix \
            "${SID}_Report.html" "${SID}_Report.json" "${SID}_summary.csv"; do
  [ -e "$OUTS/$item" ] || { echo "[err] 缺失产出 $OUTS/$item"; exit 1; }
done
echo "  产出验证通过（$(du -sh "$OUTS" | cut -f1)）"

# 6. rsync 交付（outs 重命名为 outs/；大文件必须 rsync，cp 写挂载盘会内存错误）
mkdir -p "$DELIVERY/$SID"
rsync -a --append-verify "$OUTS/" "$DELIVERY/$SID/outs/"
cp -r "$OUT_ROOT/_flagdone" "$OUT_ROOT/logs" "$OUT_ROOT/run_analysis_cmds.txt" "$DELIVERY/$SID/"
# 交付校验：无 cellsorted 临时 BAM、BAM 字节一致
[ -z "$(ls "$DELIVERY/$SID/outs/"*cellsorted* 2>/dev/null || true)" ] || { echo "[err] 交付目录含 cellsorted 临时 BAM"; exit 1; }
src_sz=$(stat -c%s "$OUTS/${SID}_Aligned.sort.bam")
dst_sz=$(stat -c%s "$DELIVERY/$SID/outs/${SID}_Aligned.sort.bam")
[ "$src_sz" = "$dst_sz" ] || { echo "[err] BAM 交付大小不一致: $src_sz != $dst_sz"; exit 1; }
echo "  交付完成: $DELIVERY/$SID（BAM ${src_sz} B 校验一致）"

# 7. 清理（交付验证后删除 staging 与全部输出）
rm -rf "$STAGE" "$OUT_ROOT"
echo "==== [$(date '+%F %T')] mobivision quantify $SID 完成（WSL 已清空）===="
