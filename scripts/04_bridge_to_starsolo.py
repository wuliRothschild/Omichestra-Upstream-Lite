#!/usr/bin/env python3
"""桥接脚本：比对矩阵（mtx 三件套） + velocyto(loom) → STARsolo 目录结构
供 nf-scautoqc gather_matrices.py --gather_mode starsolo 使用。

产物：
  outdir/Gene/filtered/{barcodes,features}.tsv.gz + matrix.mtx.gz   (filtered 矩阵)
  outdir/Velocyto/{barcodes,features}.tsv.gz + {spliced,unspliced,ambiguous}.mtx.gz

barcode 统一为 16bp（去 cellranger 的 -1、去 velocyto 的 sample: 前缀和 x 后缀），
与 cellbender 输出对齐；基因以比对矩阵顺序为准，用 ENSG/ENSMUSG（features.tsv
第 1 列 / loom Accession）做对齐以规避重复基因符号。

输入为 mtx 目录而非 h5：交付 h5 的行名语义随来源不同（CellRanger 原生=ENSG id；
make_10x_h5 构建=make.unique 符号 id，见 02b_format_matrices.py），mtx 目录的
features.tsv 第 1 列恒为 ENSG/ENSMUSG，是两种分支下统一的对齐键。
features.tsv 第 2 列（符号）经 anndata make_index_unique 去重写入产物（与原
scanpy read_10x_h5 路径的产物格式一致）。
"""
import argparse
import gzip
import os

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad
from anndata.utils import make_index_unique


def write_mtx_gz(path, mat):
    """mat: cells×genes sparse，流式写出 10x 规范的 genes×cells .mtx.gz
    （直接写 gzip 文件对象，不落未压缩中间文件，避免大矩阵磁盘峰值翻倍）"""
    with gzip.open(path, "wb") as f:
        scipy.io.mmwrite(f, sp.csr_matrix(mat.T))  # 转置为 genes×cells


def write_lines_gz(path, lines):
    with gzip.open(path, "wt") as f:
        f.write("\n".join(lines) + "\n")


def read_mtx_dir(mtx_dir):
    """读 10x mtx 三件套 → (genes×cells CSR, barcodes, gene_ids, unique_symbols, feature_types)"""
    def _read(fn):
        with gzip.open(os.path.join(mtx_dir, fn), "rt") as f:
            return [l.rstrip("\n") for l in f]

    barcodes = _read("barcodes.tsv.gz")
    feats = [l.split("\t") for l in _read("features.tsv.gz")]
    gene_ids = [x[0] for x in feats]
    symbols = [x[1] if len(x) > 1 and x[1] not in ("", "NA") else x[0] for x in feats]
    ftypes = [x[2] if len(x) > 2 else "Gene Expression" for x in feats]
    X = sp.csr_matrix(scipy.io.mmread(os.path.join(mtx_dir, "matrix.mtx.gz")))  # genes×cells
    if X.shape != (len(gene_ids), len(barcodes)):
        raise SystemExit(f"[err] 矩阵维度 {X.shape} 与 features({len(gene_ids)})/barcodes({len(barcodes)}) 不符")
    uniq_symbols = list(make_index_unique(pd.Index(symbols)))  # 重复符号 -1/-2 后缀
    return X, barcodes, gene_ids, uniq_symbols, ftypes


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cr_mtx", required=True, help="filtered_feature_bc_matrix/ 目录（mtx 三件套）")
    ap.add_argument("--velo_loom", required=True, help="velocyto .loom")
    ap.add_argument("--outdir", required=True, help="输出根目录")
    args = ap.parse_args()

    # 1. 比对 filtered 矩阵 → Gene/filtered
    crX, cr_bc, gene_ids, uniq_symbols, ftypes = read_mtx_dir(args.cr_mtx)
    barcodes = pd.Index(cr_bc).str.replace("-1$", "", regex=True)  # 16bp 无后缀
    feats = [f"{gid}\t{name}\t{ft}" for gid, name, ft in zip(gene_ids, uniq_symbols, ftypes)]

    gene_dir = os.path.join(args.outdir, "Gene", "filtered")
    os.makedirs(gene_dir, exist_ok=True)
    write_lines_gz(os.path.join(gene_dir, "barcodes.tsv.gz"), list(barcodes))
    write_lines_gz(os.path.join(gene_dir, "features.tsv.gz"), feats)
    write_mtx_gz(os.path.join(gene_dir, "matrix.mtx.gz"), crX.T)

    # 2. velocyto loom → Velocyto（spliced/unspliced/ambiguous 三 mtx）
    vl = ad.read_loom(args.velo_loom)
    vl.obs_names = vl.obs_names.str.split(":").str[-1].str.replace("x$", "", regex=True)  # 16bp 无前后缀

    # 用 ENSG/ENSMUSG 对齐 velocyto 基因到比对矩阵顺序（规避重复基因符号）
    acc = vl.var["Accession"]
    dup = acc[acc.duplicated()].unique()
    if len(dup):
        raise SystemExit(f"[err] loom Accession 存在重复（{len(dup)} 个，前5: {list(dup[:5])}）——无法唯一对齐基因顺序")
    missing = sorted(set(gene_ids) - set(acc))
    if missing:
        raise SystemExit(f"[err] {len(missing)} 个比对矩阵基因在 loom 中缺失（前5: {missing[:5]}）——检查 GTF/参考是否一致")
    vl_pos = pd.Series(np.arange(vl.n_vars), index=acc)
    order = vl_pos.loc[gene_ids].values  # 比对矩阵基因顺序在 vl 中的位置

    velo_dir = os.path.join(args.outdir, "Velocyto")
    os.makedirs(velo_dir, exist_ok=True)
    write_lines_gz(os.path.join(velo_dir, "barcodes.tsv.gz"), list(vl.obs_names))
    write_lines_gz(os.path.join(velo_dir, "features.tsv.gz"), feats)  # 同比对矩阵 features（已对齐顺序）
    for layer_name in ("spliced", "unspliced", "ambiguous"):
        write_mtx_gz(os.path.join(velo_dir, f"{layer_name}.mtx.gz"), vl.layers[layer_name][:, order])

    print(f"[bridge] {crX.shape[1]} cells × {crX.shape[0]} genes → {args.outdir}")


if __name__ == "__main__":
    main()
