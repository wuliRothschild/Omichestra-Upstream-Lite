#!/usr/bin/env python
# ========================================================================
# 06e scDblFinder 簇级 BH 检验（scautoqc 阶段 步骤4c，08 之后、pool 之前）
#
# 机制：对 08 双运行的两列 score（score_scdblfinder /
#       score_scdblfinder_homotypic）各做一次簇级 BH 检验，口径与 sctk
#       run_scrublet 严格同构（sctk/_doublets.py 源码实证，检验函数直接
#       import 复用，零重写漂移）：
#         1) 在 postqc 全量细胞上按 sctk 聚类配方重算簇（leiden 默认
#            random_state=0，同机同版本下与 scrublet 运行同簇）；
#         2) 簇 score 中位数 → 上侧 MAD z → 正态 p → BH（fdr_bh）；
#         3) 细胞取其簇的校正值，得 scdblfinder_bh_val /
#            scdblfinder_homotypic_bh_val。
#       stringent 判定（用户裁定）：score>0.5 | bh_val<0.05，两运行各一列。
#       0.5 对齐 integration 既有剔除阈值（score<0.5 保留）。
#
# 与 06d 阈值闸联动：簇中位数只用闸通过细胞的 score（被闸细胞的 score=1
#       是回填惯例而非测量值）；被闸细胞与无闸通过成员的簇 bh_val=1；
#       被闸细胞因 score=1>0.5 其 stringent=TRUE（与回填语义一致）。
#
# 兜底（镜像 sctk ignore_error 哲学）：聚类/检验异常或上侧 MAD 失效
#       （如全部簇中位数相等）→ 对应运行 bh_val 全量记 1 并日志标注，
#       不中断流水线（bh_val=1 表示"无簇级证据"而非出错）。
#
# 用法：
#   python 06e_scdblfinder_cluster_qc.py --postqc <postqc.h5ad> --sid <sid>
#        --scdblfinder-csv <3_doublet_scores/<sid>_scdblfinder.csv>
#        --gate-csv <2_qc_objects/<sid>_threshold_gate.csv>
#        --out <3_doublet_scores/<sid>_scdblfinder_cluster_qc.csv>
# 输出：barcode + scdblfinder_bh_val/scdblfinder_homotypic_bh_val/
#       scdblfinder_stringent/scdblfinder_homotypic_stringent 四列，
#       覆盖 postqc 全量细胞（只标注不过滤）。
# ========================================================================
import argparse
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from sctk._doublets import _calculate_z, _test_outlier  # noqa: F401  （_test_outlier 内部依赖 _calculate_z，同模块复用）

BH_ALPHA = 0.05        # stringent 的 bh_val 阈值（对齐 scrublet stringent 档）
SCORE_CUT = 0.5        # stringent 的 score 阈值（对齐 integration 剔除阈值）

# (输出前缀, score 列名)
RUNS = [
    ("scdblfinder", "score_scdblfinder"),
    ("scdblfinder_homotypic", "score_scdblfinder_homotypic"),
]


def _sctk_recipe_clusters(X, obs_names):
    """按 sctk run_scrublet 内联配方重算簇（逐行对齐 sctk/_doublets.py）。"""
    adata_copy = anndata.AnnData(X=X, obs=pd.DataFrame(index=obs_names))
    sc.pp.filter_genes(adata_copy, min_cells=3)
    sc.pp.normalize_total(adata_copy, target_sum=1e4)
    sc.pp.log1p(adata_copy)
    sc.pp.highly_variable_genes(
        adata_copy, min_mean=0.0125, max_mean=3, min_disp=0.5, subset=True
    )
    sc.pp.scale(adata_copy, zero_center=False)
    sc.pp.pca(adata_copy, svd_solver="arpack", zero_center=False)
    sc.pp.neighbors(adata_copy, n_pcs=30)
    sc.tl.leiden(adata_copy, resolution=1)
    resolution_function = lambda x: np.maximum(np.maximum(np.log10(x) - 1, 0) ** 2, 0.1)  # noqa: E731
    for clst in np.unique(adata_copy.obs["leiden"]):
        clst_size = sum(adata_copy.obs["leiden"] == clst)
        sc.tl.leiden(
            adata_copy,
            restrict_to=("leiden", [clst]),
            resolution=resolution_function(clst_size),
            key_added="leiden_R",
        )
        adata_copy.obs["leiden"] = adata_copy.obs["leiden_R"]
    return adata_copy.obs["leiden"].astype(str).values


def _cluster_bh(score, clusters, fit_mask):
    """单运行的簇级 BH 检验。score/clusters 全量长度；fit_mask=参与统计的细胞
    （闸通过）。返回全量 bh_val（被排除细胞与无有效成员的簇记 1）。"""
    bh_val = np.ones(score.shape[0], dtype=float)
    clst_ids = np.unique(clusters)
    meds, kept = [], []
    for c in clst_ids:
        m = (clusters == c) & fit_mask
        if m.sum() == 0:
            continue  # 无闸通过成员的簇：bh_val 维持 1
        meds.append(np.median(score[m]))
        kept.append(c)
    if len(meds) < 2:
        raise ValueError(f"有效簇数不足（{len(meds)}），簇级检验不可用")
    pvals, bh_pvals, _z = _test_outlier(np.asarray(meds))
    if not np.all(np.isfinite(bh_pvals)):
        raise ValueError("上侧 MAD 失效（簇中位数退化），簇级检验不可用")
    bh_by_clst = dict(zip(kept, bh_pvals))
    for i, c in enumerate(clusters):
        if fit_mask[i] and c in bh_by_clst:
            bh_val[i] = bh_by_clst[c]
    return bh_val


def main():
    ap = argparse.ArgumentParser(description="06e scDblFinder 簇级 BH 检验（sctk 同构）")
    ap.add_argument("--postqc", required=True)
    ap.add_argument("--sid", required=True)
    ap.add_argument("--scdblfinder-csv", required=True)
    ap.add_argument("--gate-csv", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    for p in (args.postqc, args.scdblfinder_csv, args.gate_csv):
        if not os.path.exists(p):
            raise SystemExit(f"[err] 输入缺失: {p}")

    a = ad.read_h5ad(args.postqc)
    n = a.n_obs
    obs_names = a.obs_names.astype(str)

    scdf = pd.read_csv(args.scdblfinder_csv, dtype={"barcode": str})
    need_cols = {"barcode"} | {c for _, c in RUNS}
    missing = need_cols - set(scdf.columns)
    if missing:
        raise SystemExit(f"[err] {args.scdblfinder_csv} 缺列: {sorted(missing)}")
    scdf = scdf.set_index("barcode")

    gate = pd.read_csv(args.gate_csv, dtype={"barcode": str}).set_index("barcode")
    if "pass_threshold_gate" not in gate.columns:
        raise SystemExit(f"[err] {args.gate_csv} 缺 pass_threshold_gate 列")

    # 对齐与硬校验（全量覆盖，缺一即报错，对齐 09 的 NA 硬校验风格）
    if not obs_names.isin(scdf.index).all():
        miss = int((~obs_names.isin(scdf.index)).sum())
        raise SystemExit(f"[err] scdblfinder csv 缺 {miss} 个 barcode（08 输出未覆盖 postqc 全量）")
    if not obs_names.isin(gate.index).all():
        miss = int((~obs_names.isin(gate.index)).sum())
        raise SystemExit(f"[err] threshold_gate csv 缺 {miss} 个 barcode（06d 输出未覆盖 postqc 全量）")
    scdf = scdf.loc[obs_names]
    gate_pass = gate.loc[obs_names, "pass_threshold_gate"].astype(str).str.lower().eq("true").to_numpy()

    scores = {}
    for name, col in RUNS:
        v = pd.to_numeric(scdf[col], errors="coerce").to_numpy(dtype=float)
        if np.isnan(v).any():
            raise SystemExit(f"[err] {col} 存在 NA/非数值（{int(np.isnan(v).sum())} 个），中止")
        scores[name] = v
    print(f"[06e] {args.sid}: 细胞 {n}，闸通过 {int(gate_pass.sum())}，被闸 {int((~gate_pass).sum())}", flush=True)

    # 聚类（全量细胞，与 scrublet 运行同细胞集合）；聚类失败则两运行一起兜底
    try:
        clusters = _sctk_recipe_clusters(a.X, obs_names)
        n_clst = len(np.unique(clusters))
        print(f"[06e] {args.sid}: 聚类完成（sctk 配方），{n_clst} 簇", flush=True)
    except (ValueError, ZeroDivisionError) as e:
        clusters = None
        print(f"[warn] {args.sid}: 聚类失败（{e}），两运行 bh_val 全量回填 1", flush=True)

    out = pd.DataFrame({"barcode": obs_names})
    for name, _col in RUNS:
        bh_col, str_col = f"{name}_bh_val", f"{name}_stringent"
        if clusters is None:
            bh_val = np.ones(n, dtype=float)
        else:
            try:
                bh_val = _cluster_bh(scores[name], clusters, gate_pass)
            except (ValueError, ZeroDivisionError) as e:
                bh_val = np.ones(n, dtype=float)
                print(f"[warn] {args.sid}: {name} 簇级检验不可用（{e}），bh_val 全量回填 1", flush=True)
        stringent = (scores[name] > SCORE_CUT) | (bh_val < BH_ALPHA)
        out[bh_col] = bh_val
        out[str_col] = stringent
        print(
            f"[06e] {args.sid}: {name} bh_val<{BH_ALPHA} 细胞 {int((bh_val < BH_ALPHA).sum())}，"
            f"stringent=TRUE {int(stringent.sum())}",
            flush=True,
        )

    out.to_csv(args.out, index=False)
    print(f"[06e] {args.sid}: 写出 {args.out}（{len(out)} 行）", flush=True)


if __name__ == "__main__":
    sys.exit(main())
