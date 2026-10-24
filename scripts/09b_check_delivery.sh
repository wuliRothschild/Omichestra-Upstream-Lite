#!/usr/bin/env bash
# ============================================================================
# 交付完整性校验（全部跑完后执行；任一必需项缺失即非零退出，打印缺失明细）
# 用法: 09b_check_delivery.sh [run_config.sh]
# 校验项（对应交付契约）：
#  1. 每样本 outs 标准四件：{raw,filtered}_feature_bc_matrix/（mtx 三件套）+ 对应 .h5
#  2. 每样本 velocyto/<sid>.loom
#  3. 每样本 outs/cellbender/ 标准输出（3 h5 + cell_barcodes.csv + metrics.csv + pdf；log 仅警告）
#  4. 每样本 outs/metadata.json（十字段，含 species）
#  5. 质控产出：scautoqc/ 官方命名全套（单样本数据集按单样本止步 止步合并前：
#       根部合并件不检查、warn 提示，逐样本产出仍校验；scDblFinder 双运行四列校验：
#       逐样本 csv 表头四列 + pooled/pooled_filtered obs 四列；簇级检验校验：逐样本
#       3_doublet_scores/<sid>_scdblfinder_cluster_qc.csv 存在且表头含四列
#       （scdblfinder_bh_val/scdblfinder_homotypic_bh_val/scdblfinder_stringent/
#       scdblfinder_homotypic_stringent）+ 交付对象 obs 含该四列与改名后
#       scrublet_bh_val/scrublet_stringent 且不含旧名 bh_pval/stringent_doublet；
#       阈值闸校验：逐样本 2_qc_objects/<sid>_threshold_gate.csv 存在且含
#       pass_threshold_gate 列 + pooled/pooled_filtered obs 含该列）；
#  6. _runmeta 指纹记录（指纹判据）：每样本阶段指纹齐备、字段完整、合并级成员一致
# ============================================================================
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CFG="run_config.sh"; BRANCH="scautoqc"
while [ $# -gt 0 ]; do
  case "$1" in
    *) [ -f "$1" ] && CFG="$1"; shift ;;
  esac
done
[ -f "$CFG" ] || CFG="$SCRIPT_DIR/$CFG"
# shellcheck disable=SC1090
source "$CFG"

CR_ROOT="$OUTPUT_ROOT/$PROJECT_NAME"
DS="${DATASET_NAME:-$PROJECT_NAME}"
SC_DIR="${SC_ROOT:-$CR_ROOT/scautoqc}"
[ -n "${SCANPY_PYTHON:-}" ] || { echo "[config 错误] SCANPY_PYTHON 未填（交付校验需用指定 python 实例解析 metadata.json，禁止 PATH 抽签）" >&2; exit 1; }
mapfile -t SIDS < <(sed '1s/^﻿//' "$SAMPLE_SHEET" | tr -d '\r' | awk -F, 'NR>1 && $1!~/^#/ && NF>=2 {print $1}')
[ "${#SIDS[@]}" -gt 0 ] || { echo "[err] 样本表为空: $SAMPLE_SHEET"; exit 1; }

ERRS=()
WARN=()

# ---- 1-4. per-sample 契约 ----
for SID in "${SIDS[@]}"; do
  OUTS="$CR_ROOT/$SID/outs"
  for kind in raw filtered; do
    for f in barcodes.tsv.gz features.tsv.gz matrix.mtx.gz; do
      [ -f "$OUTS/${kind}_feature_bc_matrix/$f" ] || ERRS+=("$SID 缺 outs/${kind}_feature_bc_matrix/$f")
    done
    [ -f "$OUTS/${kind}_feature_bc_matrix.h5" ] || ERRS+=("$SID 缺 outs/${kind}_feature_bc_matrix.h5")
  done
  [ -f "$CR_ROOT/$SID/velocyto/$SID.loom" ] || ERRS+=("$SID 缺 velocyto/$SID.loom")
  CB="$OUTS/cellbender"
  for f in cellbender_matrix.h5 cellbender_matrix_filtered.h5 cellbender_matrix_posterior.h5 \
           cellbender_matrix_cell_barcodes.csv cellbender_matrix_metrics.csv cellbender_matrix.pdf; do
    [ -f "$CB/$f" ] || ERRS+=("$SID 缺 outs/cellbender/$f")
  done
  { [ -f "$CB/cellbender_matrix.log" ] || [ -f "$CB/cellbender.log" ]; } || WARN+=("$SID cellbender 无日志文件")
  if [ -f "$OUTS/metadata.json" ]; then
    "$SCANPY_PYTHON" - "$OUTS/metadata.json" <<'EOF' || ERRS+=("$SID metadata.json 十字段不全/解析失败")
import json, sys
m = json.load(open(sys.argv[1], encoding="utf-8"))
need = ("group","age","batch","datasets","sample","tech","gender","tissue","chemistry","species")
missing = [k for k in need if k not in m]
sys.exit(1 if missing else 0)
EOF
  else
    ERRS+=("$SID 缺 outs/metadata.json")
  fi
done

# ---- 5. 分支产出 ----
case "$BRANCH" in
  scautoqc)
    SC="$SC_DIR"
    for d in 1_gathered_objects 2_qc_objects 2_qc_plots_individual 3_doublet_scores 5_qc_plots_overall; do
      [ -d "$SC/$d" ] || ERRS+=("scautoqc 缺目录 $d/")
    done
    for SID in "${SIDS[@]}"; do
      ls "$SC/1_gathered_objects/$SID".*.h5ad >/dev/null 2>&1 || ERRS+=("$SID 缺 1_gathered_objects gathered h5ad")
      [ -f "$SC/2_qc_objects/${SID}_postqc.h5ad" ] || ERRS+=("$SID 缺 2_qc_objects/${SID}_postqc.h5ad")
      [ -f "$SC/3_doublet_scores/${SID}_scrublet.csv" ] || ERRS+=("$SID 缺 3_doublet_scores/${SID}_scrublet.csv")
      [ -f "$SC/3_doublet_scores/${SID}_scdblfinder.csv" ] || ERRS+=("$SID 缺 3_doublet_scores/${SID}_scdblfinder.csv")
      [ -f "$SC/2_qc/$SID/${SID}_qc_thresholds.csv" ] || ERRS+=("$SID 缺 2_qc/$SID/${SID}_qc_thresholds.csv")
      # 阈值闸（06d，Science 2026 口径硬过滤）：产物存在 + 表头含 pass_threshold_gate
      [ -f "$SC/2_qc_objects/${SID}_threshold_gate.csv" ] || ERRS+=("$SID 缺 2_qc_objects/${SID}_threshold_gate.csv")
      GHDR=$(head -n 1 "$SC/2_qc_objects/${SID}_threshold_gate.csv" 2>/dev/null)
      case "$GHDR" in *pass_threshold_gate*) ;; *) ERRS+=("$SID threshold_gate.csv 缺列 pass_threshold_gate");; esac
      # scDblFinder 双运行四列（默认 + 同型，官方 FAQ 口径）
      HDR=$(head -n 1 "$SC/3_doublet_scores/${SID}_scdblfinder.csv" 2>/dev/null)
      for c in score_scdblfinder doublet_scdblfinder score_scdblfinder_homotypic doublet_scdblfinder_homotypic; do
        case "$HDR" in *"$c"*) ;; *) ERRS+=("$SID scdblfinder.csv 缺列 $c（双运行四列契约）");; esac
      done
      # 簇级检验（06e，sctk 同构）：产物存在 + 表头四列
      [ -f "$SC/3_doublet_scores/${SID}_scdblfinder_cluster_qc.csv" ] || ERRS+=("$SID 缺 3_doublet_scores/${SID}_scdblfinder_cluster_qc.csv")
      CHDR=$(head -n 1 "$SC/3_doublet_scores/${SID}_scdblfinder_cluster_qc.csv" 2>/dev/null)
      for c in scdblfinder_bh_val scdblfinder_homotypic_bh_val scdblfinder_stringent scdblfinder_homotypic_stringent; do
        case "$CHDR" in *"$c"*) ;; *) ERRS+=("$SID scdblfinder_cluster_qc.csv 缺列 $c（簇级检验四列契约）");; esac
      done
    done
    # 单样本数据集（关键单样本止步）：06 在 pool_all/finalize_qc 前有意止步——根部合并
    # 产物（pooled 对象/合并阈值表/通过率表）不产出是契约设计，不是交付缺失；
    # 单样本时只 warn 不 FAIL，校验逐样本产出。
    if [ "${#SIDS[@]}" -eq 1 ]; then
      WARN+=("数据集仅包含 1 个样本（${SIDS[*]}）：scautoqc 按单样本契约止步于合并前，根部 pooled 三件套/合并阈值表不产出（见 SKILL.md 关键单样本止步）")
    else
      for f in qc_thresholds.csv scautoqc_pooled0.h5ad scautoqc_pooled.h5ad scautoqc_pooled_filtered.h5ad sample_passqc_df.csv scdblfinder_scores.csv; do
        [ -f "$SC/$f" ] || ERRS+=("scautoqc 缺根部 $f")
      done
      # pooled/pooled_filtered obs 必须含 scDblFinder 双运行四列 + 簇级检验四列 + 阈值闸列
      # + 改名后 scrublet 判定列，且不得残留旧名（backed 只读 obs，不载矩阵）
      "$SCANPY_PYTHON" - "$SC/scautoqc_pooled.h5ad" "$SC/scautoqc_pooled_filtered.h5ad" <<'EOF' \
        || ERRS+=("pooled obs 双胞列契约不满足（缺必需列/残留旧名/不可读，见上方明细）")
import sys
import anndata as ad
need = ("score_scdblfinder", "doublet_scdblfinder",
        "score_scdblfinder_homotypic", "doublet_scdblfinder_homotypic",
        "scdblfinder_bh_val", "scdblfinder_homotypic_bh_val",
        "scdblfinder_stringent", "scdblfinder_homotypic_stringent",
        "pass_threshold_gate", "scrublet_bh_val", "scrublet_stringent")
forbidden = ("bh_pval", "stringent_doublet")
bad = False
for p in sys.argv[1:]:
    a = ad.read_h5ad(p, backed="r")
    miss = [c for c in need if c not in a.obs.columns]
    stale = [c for c in forbidden if c in a.obs.columns]
    if miss:
        bad = True
        print("  [err] %s obs 缺列: %s" % (p, miss))
    if stale:
        bad = True
        print("  [err] %s obs 残留旧名（应已更名 scrublet_bh_val/scrublet_stringent）: %s" % (p, stale))
    try:
        a.file.close()
    except Exception:
        pass
sys.exit(1 if bad else 0)
EOF
    fi
    ;;
esac

# ---- 6. _runmeta 指纹记录（指纹判据：阶段指纹齐备 + 字段完整 + 合并级成员一致）----
if [ -n "$BRANCH" ]; then
  "$SCANPY_PYTHON" "$SCRIPT_DIR/runmeta.py" verify --dir "$CR_ROOT/_runmeta" \
    --samples "$(IFS=,; echo "${SIDS[*]}")" --branch "$BRANCH" \
    || ERRS+=("_runmeta 指纹记录不完整（见上方 [fingerprint] 明细）")
fi

# ---- 汇总 ----
for w in "${WARN[@]:-}"; do [ -n "$w" ] && echo "[warn] $w"; done
if [ "${#ERRS[@]}" -gt 0 ]; then
  echo "[FAIL] 交付校验未通过（${#ERRS[@]} 项缺失）："
  for e in "${ERRS[@]}"; do echo "  - $e"; done
  exit 1
fi
echo "[OK] 交付完整：${#SIDS[@]} 样本 | 分支=${BRANCH:-未知} | 根=$CR_ROOT"
