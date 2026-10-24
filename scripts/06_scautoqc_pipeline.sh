#!/usr/bin/env bash
# ============================================================================
# scautoqc 全流程（original，per-sample）：bridge→gather→qc→threshold_gate→flag_doublet/scDblFinder/scDblFinder 簇级检验 + pool_all + finalize_qc
# 用法: 06_scautoqc_pipeline.sh [run_config.sh]
# 产出: $OUTPUT_ROOT/$PROJECT_NAME/scautoqc/（可用 run_config SC_ROOT 覆盖落盘），目录命名完全对齐 nf-scautoqc 官方：
#       1_gathered_objects/2_qc_objects/2_qc_plots_individual/3_doublet_scores/
#       5_qc_plots_overall/ + 根部 qc_thresholds.csv/scautoqc_pooled{,_filtered}.h5ad/
#       sample_passqc_df.csv（+ skill 附加 scdblfinder_scores.csv）
#       中间产物保留在 1_bridge/<sid>/ 与 2_qc/<sid>/（可断点重跑，非官方交付项）
#       最终 scautoqc_pooled_filtered.h5ad（质控好的对象，供下游去批次/注释）
# 幂等（指纹判据）：每样本与合并级 skip 均以输入指纹为准——_runmeta/<sid>.scautoqc.fingerprint.json、
#       _runmeta/merged_scautoqc_pooled_filtered.fingerprint.json（成员指纹集合 + 合并参数）；
#       产物存在但指纹不符/无记录即重跑该范围（marker/参数/脚本/输入任一变化触发）
# 单样本数据集：在多样本处理（pool_all/finalize_qc）前止步并显式警告——单样本做"合并"
#       没有意义（scVI/harmony 等多样本整合对单样本是空跑，属 0/0 式无意义形式化）；
#       这不是软件缺陷，而是有意设防。仅保留单样本质控结果，不产出 pooled 三件套、
#       不写 HANDOFF 契约段落；请从根本上避免（合并多个数据集再跑 / 全程单样本分析路径）
# ============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CFG="${1:-run_config.sh}"; [ -f "$CFG" ] || CFG="$SCRIPT_DIR/$CFG"
# shellcheck disable=SC1090
source "$CFG"
source "$CONDA_SH"; conda activate "$ENV_SCAUTOQC"
export BASE_DIR="$NF_SCAUTOQC_DIR"   # finalize_qc 据此定位 genes_list

CR_ROOT="$OUTPUT_ROOT/$PROJECT_NAME"
DS="${DATASET_NAME:-$PROJECT_NAME}"
SC_ROOT="${SC_ROOT:-$CR_ROOT/scautoqc}"
mkdir -p "$SC_ROOT/1_gathered_objects" "$SC_ROOT/2_qc_objects" "$SC_ROOT/2_qc_plots_individual" \
  "$SC_ROOT/3_doublet_scores" "$SC_ROOT/5_qc_plots_overall"

# tech/chemistry：run_config 显式值优先，否则按 ALIGNER 推导（写 metadata.json 与 HANDOFF 用）
META_TECH="${TECH:-}"; [ -z "$META_TECH" ] && { [ "${ALIGNER:-cellranger}" = "mobivision" ] && META_TECH="mobivision" || META_TECH="sc"; }
META_CHEM="${CHEMISTRY:-}"; [ -z "$META_CHEM" ] && { [ "${ALIGNER:-cellranger}" = "mobivision" ] && META_CHEM="mobi_3prime_v1" || META_CHEM="10x_3prime_v3"; }

# ---- 0. 准备 genes_list（物种适配，写入 nf-scautoqc/genes_list/）----
echo ">>> [0] 准备 genes_list (species=$SPECIES)"
python "$SCRIPT_DIR/05_prepare_genes_list.py" --species "$SPECIES" \
  --out-dir "$NF_SCAUTOQC_DIR/genes_list"

# ---- 读样本表 ----
mapfile -t SIDS < <(sed '1s/^﻿//' "$SAMPLE_SHEET" | tr -d '\r' | awk -F, 'NR>1 && $1!~/^#/ && NF>=2 {print $1}')

# 数据集一致性硬校验：一次运行=一个数据集（tissue×tech×chemistry 唯一组合）；
# tissue/chemistry 为运行级单值，此处校验样本表 cell_or_nuclei（sc/sn）全一致——混批即拒跑
mapfile -t CONS < <(sed '1s/^﻿//' "$SAMPLE_SHEET" | tr -d '\r' | awk -F, 'NR>1 && $1!~/^#/ && NF>=2 {v=$5; gsub(/[ ]/,"",v); print (v==""?"cell":v)}')
[ "${#CONS[@]}" -eq "${#SIDS[@]}" ] || { echo "[err] 样本表 cell_or_nuclei 列与样本数不符: $SAMPLE_SHEET"; exit 1; }
UNIQ_CON=$(printf '%s\n' "${CONS[@]}" | sort -u | paste -sd,)
[ "$(printf '%s\n' "${CONS[@]}" | sort -u | wc -l)" -eq 1 ] \
  || { echo "[err] 样本表 cell_or_nuclei 混杂（$UNIQ_CON）：单细胞(sc)与单细胞核(sn)不得归入同一数据集，请拆分样本表分批运行"; exit 1; }
echo ">>> 数据集一致性: cell_or_nuclei=$UNIQ_CON（${#SIDS[@]} 样本）"

# 本阶段相关工具版本（进指纹：版本变化即重跑，指纹判据）
VER_STACK=$(python -c "import importlib.metadata as m
def v(p):
    try: return m.version(p)
    except Exception: return 'unknown'
print(','.join(f'{p}={v(p)}' for p in ('scanpy','sctk','scrublet')))" 2>/dev/null || echo "unknown")

# Windows R 中文防护：R 4.5.3 UCRT 下 "C.UTF-8" 名称不可用；非 UTF-8 locale（GBK）会使
# 含中文的 R 源解析与日志乱码。R 调用显式指定 UTF-8 locale（Windows 用 .utf8/en_US.UTF-8，
# WSL/容器可用 C.UTF-8 覆盖：R_LC_CTYPE=C.UTF-8）
R_LC_CTYPE="${R_LC_CTYPE:-en_US.UTF-8}"

# scDblFinder 前阈值闸参数（06d；Science 2026 口径契约默认，run_config GATE_* 可覆盖；
# 实取值入指纹 params——阈值变化即触发受影响样本重跑，指纹判据）
GATE_MIN_COUNTS="${GATE_MIN_COUNTS:-200}"
GATE_MIN_GENES="${GATE_MIN_GENES:-200}"
GATE_MAX_MITO_CELL="${GATE_MAX_MITO_CELL:-30}"
GATE_MAX_MITO_NUCLEI="${GATE_MAX_MITO_NUCLEI:-5}"

# ---- 1-3. per-sample: bridge → gather → qc original → flag_doublet ----
for SID in "${SIDS[@]}"; do
  echo "########## [$SID] scautoqc start $(date '+%F %T') ##########"
  cr_mtx="$CR_ROOT/$SID/outs/filtered_feature_bc_matrix"
  velo_loom="$CR_ROOT/$SID/velocyto/$SID.loom"
  cb_h5="$CR_ROOT/$SID/outs/cellbender/cellbender_matrix.h5"
  for f in "$cr_mtx/barcodes.tsv.gz" "$cr_mtx/features.tsv.gz" "$cr_mtx/matrix.mtx.gz" "$velo_loom" "$cb_h5"; do
    [ -f "$f" ] || { echo "[err] $SID 输入缺失: $f（检查上游 format/velocyto/cellbender 是否完成）"; exit 1; }
  done

  # cell_or_nuclei（gather 消费，也是指纹参数之一）
  CON=$(sed '1s/^﻿//' "$SAMPLE_SHEET" | tr -d '\r' | awk -F, -v s="$SID" 'NR>1 && $1==s {gsub(/[ ]/,"",$5); print $5}')
  CON="${CON:-cell}"

  # 幂等（指纹判据）：输入指纹逐项一致 + 六件产物齐备才 skip
  FP="$CR_ROOT/_runmeta/$SID.scautoqc.fingerprint.json"
  FP_ARGS=(--fp "$FP" --script-file "$SCRIPT_DIR/06_scautoqc_pipeline.sh" \
    --script-file "$SCRIPT_DIR/06d_threshold_gate.py" \
    --script-file "$SCRIPT_DIR/06e_scdblfinder_cluster_qc.py" \
    --set "qc_mode=original" --set "gmm_cutoff=inner" --set "cell_or_nuclei=$CON" --set "species=$SPECIES" \
    --set "gate_min_counts=$GATE_MIN_COUNTS" --set "gate_min_genes=$GATE_MIN_GENES" \
    --set "gate_max_mito_cell=$GATE_MAX_MITO_CELL" --set "gate_max_mito_nuclei=$GATE_MAX_MITO_NUCLEI" \
    --input "$cr_mtx/barcodes.tsv.gz" --input "$cr_mtx/features.tsv.gz" --input "$cr_mtx/matrix.mtx.gz" \
    --input "$velo_loom" --input "$cb_h5" --version "python_stack=$VER_STACK" --base "$CR_ROOT")
  if python "$SCRIPT_DIR/runmeta.py" check "${FP_ARGS[@]}" \
      --require "$SC_ROOT/2_qc_objects/${SID}_postqc.h5ad" \
      --require "$SC_ROOT/2_qc_objects/${SID}_threshold_gate.csv" \
      --require "$SC_ROOT/3_doublet_scores/${SID}_scrublet.csv" \
      --require "$SC_ROOT/3_doublet_scores/${SID}_scdblfinder.csv" \
      --require "$SC_ROOT/3_doublet_scores/${SID}_scdblfinder_cluster_qc.csv"; then
    echo "[skip] $SID 指纹一致（scautoqc 已完成）"
    continue
  fi

  # 0b. cellbender h5 features 约定对齐（make_10x_h5 的 R 点号去重 vs bridge 的 anndata
  #     连字符去重是两套系统性错配约定，不对齐 gather 必炸 ValueError；幂等，
  #     仅写 staging 副本的 features/id+name，矩阵数值/barcodes 零改动）
  python "$SCRIPT_DIR/06a_align_cellbender_features.py" --cr-root "$CR_ROOT" --samples "$SID" \
    || { echo "[err] $SID cellbender features 与 filtered mtx 不可对齐（上方 [err] 有具体原因）"; exit 1; }

  # 1. bridge → STARsolo 目录（mtx 目录 features.tsv 第 1 列 ENSG/ENSMUSG 作 loom 对齐键）
  python "$SCRIPT_DIR/04_bridge_to_starsolo.py" \
    --cr_mtx "$cr_mtx" --velo_loom "$velo_loom" --outdir "$SC_ROOT/1_bridge/$SID" 2>>"$SC_ROOT/1_bridge_$SID.log" \
    || { echo "[err] $SID bridge 失败，$SC_ROOT/1_bridge_$SID.log 尾部:"; tail -n 30 "$SC_ROOT/1_bridge_$SID.log"; exit 1; }

  # 2. gather（官方脚本 starsolo 模式，cell_or_nuclei 从样本表读，默认 cell）
  #    cell: n_counts≥1000/mito 0.1–max/spliced 50–97.5；nuclei: n_counts≥300/mito 0.1–max/spliced 0.01–100
  #    产出归入官方 1_gathered_objects/<sid>.<filename>
  ( cd "$SC_ROOT/1_bridge/$SID" && python "$NF_SCAUTOQC_DIR/bin/gather_matrices.py" \
      --cr_gene ./Gene/filtered --cr_velo ./Velocyto --cb_h5 "$cb_h5" \
      --gather_mode starsolo --cell_or_nuclei "$CON" ) 2>>"$SC_ROOT/1_bridge_$SID.log" \
    || { echo "[err] $SID gather 失败，$SC_ROOT/1_bridge_$SID.log 尾部:"; tail -n 30 "$SC_ROOT/1_bridge_$SID.log"; exit 1; }
  GATH="$SC_ROOT/1_bridge/$SID/gene_velo_cellbender.filtered.h5ad"
  [ -f "$GATH" ] || { echo "[err] $SID gather 无产出: $GATH"; exit 1; }
  mv "$GATH" "$SC_ROOT/1_gathered_objects/$SID.gene_velo_cellbender.filtered.h5ad"

  # 3. qc original（per-sample：mito循环20/50/80，4指标GMM阈值硬编码）
  #    用 skill 的 qc.py（物种无关 flags，兼容人/鼠的线粒体/核糖体/血红蛋白基因识别）
  #    postqc h5ad → 2_qc_objects/；QC 图 → 2_qc_plots_individual/；per-sample 阈值 csv 留 2_qc/（pool_all 消费）
  mkdir -p "$SC_ROOT/2_qc/$SID"
  ( cd "$SC_ROOT/2_qc/$SID" && python "$SCRIPT_DIR/qc.py" \
      --sample_id "$SID" \
      --gath_obj "$SC_ROOT/1_gathered_objects/$SID.gene_velo_cellbender.filtered.h5ad" \
      --qc_mode original --gmm_cutoff inner ) 2>>"$SC_ROOT/2_qc_$SID.log" \
    || { echo "[err] $SID qc 失败，$SC_ROOT/2_qc_$SID.log 尾部:"; tail -n 30 "$SC_ROOT/2_qc_$SID.log"; exit 1; }
  [ -f "$SC_ROOT/2_qc/$SID/${SID}_postqc.h5ad" ] || { echo "[err] $SID qc 无 postqc"; exit 1; }
  mv "$SC_ROOT/2_qc/$SID/${SID}_postqc.h5ad" "$SC_ROOT/2_qc_objects/${SID}_postqc.h5ad"
  find "$SC_ROOT/2_qc/$SID" -maxdepth 1 -name "*.png" -exec mv {} "$SC_ROOT/2_qc_plots_individual/" \;

  # 3b. 阈值闸（06d；Science 2026 口径硬过滤：counts/genes/mito 纯阈值，sc/sn 分流 mito；
  #     不做 GMM——GMM 属 qc.py cellwise_qc 职责；scrublet 不入滤）
  #     产 barcode 级 mask → 2_qc_objects/（08 消费作鉴定输入；09 合并入 obs + 裁剪 filtered）
  python "$SCRIPT_DIR/06d_threshold_gate.py" \
    --postqc "$SC_ROOT/2_qc_objects/${SID}_postqc.h5ad" --sid "$SID" \
    --out "$SC_ROOT/2_qc_objects/${SID}_threshold_gate.csv" \
    --cell-or-nuclei "$CON" \
    --gate-min-counts "$GATE_MIN_COUNTS" --gate-min-genes "$GATE_MIN_GENES" \
    --gate-max-mito-cell "$GATE_MAX_MITO_CELL" --gate-max-mito-nuclei "$GATE_MAX_MITO_NUCLEI" \
    2>>"$SC_ROOT/2_qc_$SID.log" \
    || { echo "[err] $SID 阈值闸失败，$SC_ROOT/2_qc_$SID.log 尾部:"; tail -n 30 "$SC_ROOT/2_qc_$SID.log"; exit 1; }

  # 4. flag_doublet（scrublet；直接产出到官方 3_doublet_scores/）
  ( cd "$SC_ROOT/3_doublet_scores" && python "$NF_SCAUTOQC_DIR/bin/flag_doublet.py" \
      --input "$SC_ROOT/2_qc_objects/${SID}_postqc.h5ad" --samp "$SID" --filter cluster_passed_qc ) 2>>"$SC_ROOT/3_doublets_$SID.log" \
    || { echo "[err] $SID flag_doublet 失败，$SC_ROOT/3_doublets_$SID.log 尾部:"; tail -n 30 "$SC_ROOT/3_doublets_$SID.log"; exit 1; }

  # 4b. scDblFinder 双运行（per-sample postqc，第二条独立 doublet 线；鉴定输入=阈值闸
  #     通过细胞（--keep-csv，被闸细胞 score 默认 1）+ 默认/同型（clusters=FALSE,
  #     removeUnidentifiable=FALSE）两次平行鉴定，输出四列）
  LC_CTYPE="$R_LC_CTYPE" "$RSCRIPT" "$SCRIPT_DIR/08_scdblfinder.R" \
    --postqc "$SC_ROOT/2_qc_objects/${SID}_postqc.h5ad" --sid "$SID" \
    --keep-csv "$SC_ROOT/2_qc_objects/${SID}_threshold_gate.csv" \
    --out "$SC_ROOT/3_doublet_scores/${SID}_scdblfinder.csv" \
    --scanpy-python "$SCANPY_PYTHON" 2>>"$SC_ROOT/3_doublets_$SID.log" \
    || { echo "[err] $SID scDblFinder 失败，$SC_ROOT/3_doublets_$SID.log 尾部:"; tail -n 30 "$SC_ROOT/3_doublets_$SID.log"; exit 1; }

  # 4c. scDblFinder 簇级 BH 检验（06e；sctk run_scrublet 同构：sctk 配方重算簇 →
  #     簇中位数→上侧 MAD z→正态 p→BH；被闸细胞不参与统计、bh_val=1；
  #     stringent=score>0.5|bh_val<0.05；检验不可用兜底 bh_val=1 不中断；只标注不过滤）
  python "$SCRIPT_DIR/06e_scdblfinder_cluster_qc.py" \
    --postqc "$SC_ROOT/2_qc_objects/${SID}_postqc.h5ad" --sid "$SID" \
    --scdblfinder-csv "$SC_ROOT/3_doublet_scores/${SID}_scdblfinder.csv" \
    --gate-csv "$SC_ROOT/2_qc_objects/${SID}_threshold_gate.csv" \
    --out "$SC_ROOT/3_doublet_scores/${SID}_scdblfinder_cluster_qc.csv" \
    2>>"$SC_ROOT/3_doublets_$SID.log" \
    || { echo "[err] $SID scDblFinder 簇级检验失败，$SC_ROOT/3_doublets_$SID.log 尾部:"; tail -n 30 "$SC_ROOT/3_doublets_$SID.log"; exit 1; }

  # 5. 写本样本指纹（指纹判据：六件产物齐备后立据，断点续跑判据）
  python "$SCRIPT_DIR/runmeta.py" write "${FP_ARGS[@]}"
  echo "########## [$SID] done $(date '+%F %T') ##########"
done

# ---- 4-5b. 池化段（≥2 样本）。单样本数据集在此止步（有意设防，非软件缺陷）：
#      合并单个样本为 scVI/harmony 等多样本整合备料没有意义（0/0 式空跑）；
#      显式警告并要求从根本上避免（合并多个数据集 / 全程单样本分析路径）。
#      仅保留单样本质控结果（1_gathered/2_qc_objects/3_doublet_scores/qc_thresholds.csv），
#      不产出 pooled 三件套、不写契约段落（## upstream 阶段 承诺三个 pooled h5ad）----
SAMPLES=$(IFS=,; echo "${SIDS[*]}")
if [ "${#SIDS[@]}" -ge 2 ]; then
# ---- 合并级指纹（指纹判据）：成员指纹集合 + 合并参数逐项一致才 skip 池化段 ----
MERGED_FP="$CR_ROOT/_runmeta/merged_scautoqc_pooled_filtered.fingerprint.json"
MERGED_ARGS=(--fp "$MERGED_FP" --script-file "$SCRIPT_DIR/06_scautoqc_pipeline.sh" \
  --set "qc_mode=original" --set "samples=$SAMPLES" --base "$CR_ROOT")
MEMBER_ARGS=(); for s in "${SIDS[@]}"; do MEMBER_ARGS+=(--member "$CR_ROOT/_runmeta/$s.scautoqc.fingerprint.json"); done
if python "$SCRIPT_DIR/runmeta.py" check "${MERGED_ARGS[@]}" "${MEMBER_ARGS[@]}" \
    --require "$SC_ROOT/scautoqc_pooled.h5ad" --require "$SC_ROOT/scautoqc_pooled_filtered.h5ad" \
    --require "$SC_ROOT/qc_thresholds.csv" --require "$SC_ROOT/sample_passqc_df.csv" \
    --require "$SC_ROOT/scdblfinder_scores.csv"; then
  echo ">>> [4-5b] 合并级指纹一致，跳过池化段"
else
# ---- 4. pool_all（根部产出 scautoqc_pooled0.h5ad + qc_thresholds.csv，官方命名）----
OBJECTS=$(for s in "${SIDS[@]}"; do echo "$SC_ROOT/2_qc_objects/${s}_postqc.h5ad"; done | paste -sd,)
RANGES=$(for s in "${SIDS[@]}"; do echo "$SC_ROOT/2_qc/$s/${s}_qc_thresholds.csv"; done | paste -sd,)
( cd "$SC_ROOT" && python "$NF_SCAUTOQC_DIR/bin/pool_all.py" \
    --samples "$SAMPLES" --objects "$OBJECTS" --ranges "$RANGES" ) 2>>"$SC_ROOT/pool_all.log" \
  || { echo "[err] pool_all 失败，$SC_ROOT/pool_all.log 尾部:"; tail -n 30 "$SC_ROOT/pool_all.log"; exit 1; }

# ---- 5. finalize_qc（根部产出 pooled/pooled_filtered h5ad + sample_passqc_df.csv；overall 图 → 5_qc_plots_overall/）----
# 先校验每样本 scrublet csv 存在（缺失会在 finalize 内部才暴露，提前显式报错）
for s in "${SIDS[@]}"; do
  [ -f "$SC_ROOT/3_doublet_scores/${s}_scrublet.csv" ] || { echo "[err] 缺 scrublet csv: $SC_ROOT/3_doublet_scores/${s}_scrublet.csv"; exit 1; }
done
SCR=$(for s in "${SIDS[@]}"; do echo "$SC_ROOT/3_doublet_scores/${s}_scrublet.csv"; done | paste -sd,)
( cd "$SC_ROOT" && python "$NF_SCAUTOQC_DIR/bin/finalize_qc.py" \
    --obj scautoqc_pooled0.h5ad --scr "$SCR" --qc_mode original ) 2>>"$SC_ROOT/finalize_qc.log" \
  || { echo "[err] finalize_qc 失败，$SC_ROOT/finalize_qc.log 尾部:"; tail -n 30 "$SC_ROOT/finalize_qc.log"; exit 1; }
find "$SC_ROOT" -maxdepth 1 -name "*.png" -exec mv {} "$SC_ROOT/5_qc_plots_overall/" \;

# ---- 5b. 合并 scDblFinder 得分写入 pooled/filtered obs（与 scrublet 并存）----
python "$SCRIPT_DIR/09_merge_scdblfinder.py" --sc-root "$SC_ROOT" --samples "$SAMPLES" \
  2>>"$SC_ROOT/scdblfinder_merge.log"

# 合并级指纹（指纹判据：成员指纹集合 + 合并参数，成员任一失效即池化段失效）
python "$SCRIPT_DIR/runmeta.py" write "${MERGED_ARGS[@]}" "${MEMBER_ARGS[@]}"
fi
else
echo "" >&2
echo "########################################################################" >&2
echo "# [警告] 单样本数据集（${#SIDS[@]} 样本）：scautoqc 止步于多样本处理前" >&2
echo "#   pool_all/finalize_qc 是为下游 scVI/harmony 等多样本整合备料的合并步骤；" >&2
echo "#   对单样本做合并是空跑（0/0 式无意义形式化）——这不是软件缺陷，是有意设防。" >&2
echo "#   请从根本上避免单样本数据集进入本阶段：" >&2
echo "#     · 合并多个数据集为一个数据集（≥2 样本）再跑；或" >&2
echo "#     · 全程维持单样本分析路径（不进入 scautoqc 池化/整合）。" >&2
echo "#   本次仅交付单样本质控结果（1_gathered_objects/2_qc_objects/" >&2
echo "#   2_qc_plots_individual/3_doublet_scores/qc_thresholds.csv）+ metadata.json，" >&2
echo "#   不产出 pooled 三件套、不写 ## upstream 阶段 契约段落。" >&2
echo "########################################################################" >&2
echo "" >&2
fi

# ---- 5c. metadata.json：十类基础字段（已存在则跳过，缺则以 Unknown 兜底建）+ scautoqc 阶段基因记录 ----
python "$SCRIPT_DIR/01c_write_metadata.py" --cr-root "$CR_ROOT" --samples "$SAMPLES" \
  --datasets "$DS" --meta-csv "${META_CSV:-}" \
  --tech "$META_TECH" --chemistry "$META_CHEM" --tissue "${TISSUE:-Unknown}" \
  --species "${SPECIES:-Unknown}"
python "$SCRIPT_DIR/06c_record_stage_genes.py" --cr-root "$CR_ROOT" --samples "$SAMPLES"

# ---- 6. upsert ## upstream 阶段 到项目 HANDOFF.md（段级累积，给下游；单样本无 pooled 产出，跳过契约段落）----
if [ "${#SIDS[@]}" -ge 2 ]; then
# 同一 HANDOFF 承载多数据集时 DATASET_NAME 非空 → 段落标题带数据集后缀，避免互覆盖
SEC_ARGS=(); [ -n "${DATASET_NAME:-}" ] && SEC_ARGS=(--section "$DATASET_NAME")
python "$SCRIPT_DIR/07_write_handoff.py" \
  --project "$PROJECT_NAME" --species "$SPECIES" \
  --tech "$META_TECH" --chemistry "$META_CHEM" \
  --sc-root "$SC_ROOT" --samples "$SAMPLES" \
  --handoff-md "$HANDOFF_MD" ${SEC_ARGS[@]+"${SEC_ARGS[@]}"} \
  --gate-min-counts "$GATE_MIN_COUNTS" --gate-min-genes "$GATE_MIN_GENES" \
  --gate-max-mito-cell "$GATE_MAX_MITO_CELL" --gate-max-mito-nuclei "$GATE_MAX_MITO_NUCLEI"
else
echo ">>> [6] 单样本数据集：无 pooled 产出，不写 ## upstream 阶段 契约段落（见上方警告）"
fi

# ---- 7. 清理临时/缓存/垃圾文件（保留流程产物与标准输出）----
# __pycache__（python 缓存）、R 临时（.Rhistory/.RData）；中间产物（1_bridge/2_qc）保留可重跑
find "$SC_ROOT" "$SCRIPT_DIR" -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
find "$SC_ROOT" -type f \( -name ".Rhistory" -o -name ".RData" \) -delete 2>/dev/null || true
if [ "${#SIDS[@]}" -ge 2 ]; then
  echo "==== scautoqc 完成：$SC_ROOT/scautoqc_pooled_filtered.h5ad ===="
else
  echo "==== scautoqc 完成（单样本，已止步于多样本处理前）：$SC_ROOT ===="
fi
