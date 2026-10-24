#!/usr/bin/env bash
# ============================================================================
# 串行 cellbender remove-background（GPU）：对全部样本去背景
# 输入：各样本 outs/raw_feature_bc_matrix.h5（由本脚本内 format 步骤保证存在）
# 输出：<sample>/outs/cellbender/cellbender_matrix.h5
# 用法: 03_cellbender.sh [run_config.sh]（遍历样本表）
# 串行：GPU 独占资源，并行不提速反慢
# 幂等（指纹判据）：skip 判据为输入指纹（raw h5 size+mtime、参数、cellbender 版本、脚本版本）
#   逐项一致（_runmeta/<sid>.cellbender.fingerprint.json）；产物在但指纹不符/无记录即重跑
# ============================================================================
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CFG="${1:-run_config.sh}"; [ -f "$CFG" ] || CFG="$SCRIPT_DIR/$CFG"
# shellcheck disable=SC1090
source "$CFG"
source "$CONDA_SH"; conda activate "$ENV_CELLBENDER"

OUT_ROOT="$OUTPUT_ROOT/$PROJECT_NAME"
# 读样本表（跳过注释与空行），取 sample_id 列
mapfile -t SIDS < <(awk -F, 'NR>1 && $1!~/^#/ && NF>=2 {print $1}' "$SAMPLE_SHEET")

# cellbender 前置步骤一：format——确保 raw/filtered mtx 目录与 h5 四件齐全
# （MobiVision 分支在此复制标准命名 mtx 目录并构建 h5；cellbender 输入
#  raw_feature_bc_matrix.h5 由此保证；CellRanger 原生齐全时为 no-op）
VALIDATE_PYTHON="${VALIDATE_PYTHON:-$(conda info --base 2>/dev/null)/envs/${ENV_SCAUTOQC:-scautoqc}/bin/python}"
[ -x "$VALIDATE_PYTHON" ] || { echo "[err] format/校验 Python 不存在: $VALIDATE_PYTHON"; exit 1; }
"$VALIDATE_PYTHON" "$SCRIPT_DIR/02b_format_matrices.py" --cr-root "$OUT_ROOT" \
  --samples "$(IFS=,; echo "${SIDS[*]}")" --species "$SPECIES" \
  || { echo "[err] format 失败，禁止运行 cellbender"; exit 1; }

# cellbender 前置步骤二：三位一体校验（比对矩阵/loom/formatted h5）
# 任一样本不一致即停止，不跳过、不继续。
for SID in "${SIDS[@]}"; do
  "$VALIDATE_PYTHON" "$SCRIPT_DIR/03b_validate_triple_consistency.py" --sid "$SID" --cr-root "$OUT_ROOT" \
    || { echo "[err] $SID 三位一体校验失败，禁止运行 cellbender"; exit 1; }
done

# cellbender 版本（进指纹：版本变化即重跑，指纹判据）
CB_VER=$(cellbender --version 2>/dev/null | head -n 1); CB_VER="${CB_VER:-unknown}"

for SID in "${SIDS[@]}"; do
  INPUT="$OUT_ROOT/$SID/outs/raw_feature_bc_matrix.h5"
  OUTDIR="$OUT_ROOT/$SID/outs/cellbender"
  OUTFILE="$OUTDIR/cellbender_matrix.h5"
  FP="$OUT_ROOT/_runmeta/$SID.cellbender.fingerprint.json"
  FP_ARGS=(--fp "$FP" --script-file "$SCRIPT_DIR/03_cellbender.sh" --set "cuda=1" \
    --input "$INPUT" --version "cellbender=$CB_VER" --base "$OUT_ROOT")
  if "$VALIDATE_PYTHON" "$SCRIPT_DIR/runmeta.py" check "${FP_ARGS[@]}" --require "$OUTFILE"; then
    echo "[$(date '+%F %T')] skip $SID（指纹一致，cellbender 已完成）"
    continue
  fi
  echo "########## [$(date '+%F %T')] cellbender $SID ##########"
  mkdir -p "$OUTDIR"
  # 子 shell 内 cd 到样本输出目录再跑：cellbender 默认检查点 ckpt.tar.gz 固定落在 $OUTDIR，
  # 就地删除——与执行者运行本脚本时的 cwd 无关（否则检查点散落调用目录，每样本残留 1GB+）
  ( cd "$OUTDIR" && cellbender remove-background --cuda --input "$INPUT" --output "$OUTFILE" 2>&1 | tee "$OUTDIR/cellbender.log" )
  [ -f "$OUTFILE" ] || { echo "[err] $SID cellbender 无产出"; continue; }
  rm -f "$OUTDIR/ckpt.tar.gz"
  "$VALIDATE_PYTHON" "$SCRIPT_DIR/runmeta.py" write "${FP_ARGS[@]}"
  echo "########## [$(date '+%F %T')] $SID done ##########"
done
echo "ALL_CELLBENDER_DONE"