#!/usr/bin/env python3
"""交付格式规范化：确保每样本 outs/ 下 MobiVision/CellRanger 标准四件齐全——
  raw_feature_bc_matrix/  与 filtered_feature_bc_matrix/  （mtx 三件套目录）
  raw_feature_bc_matrix.h5 与 filtered_feature_bc_matrix.h5

CellRanger 原生四件齐全 → no-op；MobiVision 分支（raw_cell_gene_matrix/
filtered_cell_gene_matrix 命名、无 h5）：mtx 目录复制为标准命名，h5 用
make_10x_h5.py 基于 mtx 目录构建（行名=基因符号经 make.unique，与
Seurat::Read10X(gene.column=2) + DropletUtils::write10xCounts(version="3")
标准组合的输出逐项等价）。

任一必需源缺失即非零退出（不跳过样本）。幂等：已存在且不缺的产物不重建。
本步骤必须在 cellbender（03）之前完成——03 的输入 raw_feature_bc_matrix.h5 由此保证。

用法: 02b_format_matrices.py --cr-root <根> --samples sid1,sid2,... --species <human|mouse>
依赖: h5py/scipy/numpy（ENV_SCAUTOQC）
"""
import argparse
import os
import shutil

from make_10x_h5 import read_mtx_dir, write_10x_h5

GENOME_BY_SPECIES = {"human": "GRCh38", "mouse": "GRCm39"}
MTX_FILES = ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")


def mtx_complete(d):
    return os.path.isdir(d) and all(os.path.isfile(os.path.join(d, f)) for f in MTX_FILES)


def copy_mtx_dir(src, dst, retry=10):
    # 源与目标同处本地 OUTPUT_ROOT（非跨挂载），本地拷贝用 shutil；
    # 跨挂载/机械硬盘大文件统一走编排脚本的 retry_rsync（rsync -a --append-verify）
    os.makedirs(dst, exist_ok=True)
    for fn in MTX_FILES:
        s, d = os.path.join(src, fn), os.path.join(dst, fn)
        if os.path.exists(d):
            continue
        for i in range(retry):
            try:
                shutil.copyfile(s, d)
                break
            except OSError:
                if i == retry - 1:
                    raise
                print(f"[retry {i+1}] copy {fn}")
        if not os.path.exists(d):
            raise OSError(f"拷贝失败: {s} -> {d}")


def ensure_mtx_dir(outs, kind):
    """kind ∈ {raw, filtered}；返回标准命名 mtx 目录路径"""
    dst = os.path.join(outs, f"{kind}_feature_bc_matrix")
    if mtx_complete(dst):
        return dst
    src = os.path.join(outs, f"{kind}_cell_gene_matrix")  # MobiVision 命名
    if mtx_complete(src):
        print(f"[format] 复制 mtx 目录 {kind}_cell_gene_matrix -> {kind}_feature_bc_matrix")
        copy_mtx_dir(src, dst)
        return dst
    raise SystemExit(f"[err] {kind} mtx 三件套缺失：{dst} 与 {src} 均不完整（检查上游比对是否完成）")


def ensure_h5(outs, kind, mtx_dir, genome):
    h5 = os.path.join(outs, f"{kind}_feature_bc_matrix.h5")
    if os.path.isfile(h5):
        print(f"[format] h5 已存在: {h5}")
        return
    X, rownames, barcodes = read_mtx_dir(mtx_dir, gene_column=2)
    write_10x_h5(h5, X, rownames, barcodes, genome)
    print(f"[format] 生成 h5: {h5}（{X.shape[0]} genes × {X.shape[1]} cells, nnz={X.nnz}）")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cr-root", required=True, help="样本根目录（含 <sid>/outs/）")
    ap.add_argument("--samples", required=True, help="逗号分隔 sampleID")
    ap.add_argument("--species", required=True, help="物种（决定写入 h5 的 genome 名：human→GRCh38, mouse→GRCm39）")
    args = ap.parse_args()

    genome = GENOME_BY_SPECIES.get(args.species)
    if genome is None:
        raise SystemExit(f"[err] species={args.species} 无基因组名映射（human→GRCh38, mouse→GRCm39）；"
                         "custom 物种需先扩展 GENOME_BY_SPECIES")

    for sid in [s.strip() for s in args.samples.split(",") if s.strip()]:
        outs = os.path.join(args.cr_root, sid, "outs")
        if not os.path.isdir(outs):
            raise SystemExit(f"[err] {sid} 无 outs 目录: {outs}")
        for kind in ("raw", "filtered"):
            mtx_dir = ensure_mtx_dir(outs, kind)
            ensure_h5(outs, kind, mtx_dir, genome)
        print(f"[format] {sid} 四件齐全")


if __name__ == "__main__":
    main()
