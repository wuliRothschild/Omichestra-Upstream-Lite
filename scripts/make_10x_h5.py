#!/usr/bin/env python3
"""mtx 三件套 → 10x v3 HDF5，逐语义复刻 Seurat::Read10X(gene.column=2) +
DropletUtils::write10xCounts(version="3", type="HDF5") 标准组合的输出，
保证行名/列名/数值与 R 侧完全一致：

- 行名 = features.tsv 第 --gene-column 列（默认 2 = 基因符号，对齐 Read10X gene.column=2），
  该行名为空/NA 时回退第 1 列（ENSG；Read10X 对 NA 符号的既定回退行为），
  重复行名按 R make.unique 规则追加 .1/.2 后缀（Read10X 默认 unique.names=TRUE）
- features/id = features/name = 行名（write10xCounts 的 gene.id=gene.symbol=rownames 语义）
- 列名 = barcodes.tsv 原样（不改名、不去后缀）
- 数值 = matrix.mtx.gz 原值，CSC 布局（genes×cells）；data float64、indices/indptr/shape int32
- 定长字节串（宽度=该列最长元素）；features/_all_tag_keys=["genome"]；
  每个 dataset 带 rhdf5-NA.OK=[1] 属性（rhdf5 写出产物的固有属性）；
  全部 dataset 为 gzip(level=6)+shuffle、单 chunk=全长的存储布局（rhdf5 默认写出方式）
- root attrs（write10xCounts 固定写出）：chemistry_description="Single Cell 3' v3"、
  filetype="matrix"、library_ids="custom"、original_gem_groups=[1]、version=[2]（int32 标量组）

用法: make_10x_h5.py --mtx-dir <dir> --out <path.h5> --genome <GRCh38|GRCm39|...>
                     [--gene-column 2] [--feature-type "Gene Expression"]
依赖: h5py/scipy/numpy（ENV_SCAUTOQC）
"""
import argparse
import gzip
import os

import h5py
import numpy as np
import scipy.io
import scipy.sparse as sp


def make_unique(names):
    """R make.unique：首次出现原样，第 k 次重复追加 .{k-1}（按出现顺序）"""
    seen = {}
    out = []
    for n in names:
        if n not in seen:
            seen[n] = 0
            out.append(n)
        else:
            seen[n] += 1
            out.append(f"{n}.{seen[n]}")
    return out


def read_mtx_dir(mtx_dir, gene_column=2):
    """读 10x mtx 三件套，返回 (CSC genes×cells 矩阵, 行名 list, 列名 list)"""
    with gzip.open(os.path.join(mtx_dir, "barcodes.tsv.gz"), "rt") as f:
        barcodes = [l.rstrip("\n") for l in f]
    with gzip.open(os.path.join(mtx_dir, "features.tsv.gz"), "rt") as f:
        feats = [l.rstrip("\n").split("\t") for l in f]
    gi = gene_column - 1
    if any(len(x) <= gi for x in feats):
        raise SystemExit(f"[err] features.tsv.gz 无第 {gene_column} 列（首行列数={len(feats[0])}）")
    # 空/NA 行名回退第 1 列（Read10X 既定行为），再 make.unique
    rownames = make_unique([x[gi] if x[gi] not in ("", "NA") else x[0] for x in feats])
    X = scipy.io.mmread(os.path.join(mtx_dir, "matrix.mtx.gz"))  # genes×cells
    X = sp.csc_matrix(X)
    if X.shape != (len(rownames), len(barcodes)):
        raise SystemExit(f"[err] 矩阵维度 {X.shape} 与 features({len(rownames)})/barcodes({len(barcodes)}) 不符")
    return X, rownames, barcodes


def _create(grp, name, arr):
    """gzip6+shuffle、单 chunk=全长、rhdf5-NA.OK 属性（逐项对齐 rhdf5 写出产物）"""
    d = grp.create_dataset(name, data=arr, compression="gzip", compression_opts=6,
                           shuffle=True, chunks=(len(arr),))
    d.attrs["rhdf5-NA.OK"] = np.array([1], dtype="int32")
    return d


def _wstr(grp, name, values):
    """定长字节串 dataset（宽度=最长元素）"""
    return _create(grp, name, np.array(values, dtype="S"))


def write_10x_h5(path, X, rownames, barcodes, genome, feature_type="Gene Expression"):
    """X: CSC genes×cells。结构/ dtype 逐项对齐 write10xCounts(version=3, type=HDF5)"""
    n_genes, n_cells = X.shape
    with h5py.File(path, "w") as f:
        f.attrs["chemistry_description"] = "Single Cell 3' v3"
        f.attrs["filetype"] = "matrix"
        f.attrs["library_ids"] = "custom"
        f.attrs["original_gem_groups"] = np.array([1], dtype="int32")
        f.attrs["version"] = np.array([2], dtype="int32")
        grp = f.create_group("matrix")
        _create(grp, "barcodes", np.array(barcodes, dtype="S"))
        for name, arr in (("data", X.data.astype("float64")),
                          ("indices", X.indices.astype("int32")),
                          ("indptr", X.indptr.astype("int32")),
                          ("shape", np.array([n_genes, n_cells], dtype="int32"))):
            _create(grp, name, arr)
        feat = grp.create_group("features")
        _wstr(feat, "_all_tag_keys", ["genome"])
        _wstr(feat, "feature_type", [feature_type] * n_genes)
        _wstr(feat, "genome", [genome] * n_genes)
        _wstr(feat, "id", rownames)
        _wstr(feat, "name", rownames)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mtx-dir", required=True, help="10x mtx 目录（barcodes.tsv.gz/features.tsv.gz/matrix.mtx.gz）")
    ap.add_argument("--out", required=True, help="输出 h5 路径")
    ap.add_argument("--genome", required=True, help="写入 features/genome 的基因组名（如 GRCh38/GRCm39）")
    ap.add_argument("--gene-column", type=int, default=2, help="features.tsv 作行名的列（默认 2=符号，对齐 Read10X gene.column=2）")
    ap.add_argument("--feature-type", default="Gene Expression", help="features/feature_type 统一值（默认 Gene Expression）")
    args = ap.parse_args()

    X, rownames, barcodes = read_mtx_dir(args.mtx_dir, args.gene_column)
    write_10x_h5(args.out, X, rownames, barcodes, args.genome, args.feature_type)
    print(f"[make_10x_h5] {args.out}: {X.shape[0]} genes × {X.shape[1]} cells, nnz={X.nnz}, genome={args.genome}")


if __name__ == "__main__":
    main()
