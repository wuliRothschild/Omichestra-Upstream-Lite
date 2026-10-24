#!/usr/bin/env python3
# ============================================================================
# 06a_align_cellbender_features.py — cellbender h5 与 filtered mtx 的 features
# 约定一致性检查 + 修复（scautoqc gather 前置守卫，幂等）
#
# 背景（两次批量运行实证，稳定复现）：
#   make_10x_h5.py 按 write10xCounts 语义把 cellbender 输入 h5 的
#   features/id=features/name=基因符号，重复符号按 R make.unique 追加 .1/.2；
#   cellbender 原样透传该约定到 cellbender_matrix.h5（且行序可与 filtered
#   features.tsv 不同——raw/filtered 参考行序差异）。
#   而 04 bridge 产出的 Gene/filtered 经 anndata make_index_unique 是连字符
#   风格（X-1）。gather_matrices.py（starsolo 模式）要求 cb 与 cr_gene 的
#   var_names 集合严格相等 → 重复符号基因两风格对不上即
#   ValueError: raw features do not match main matrix (N missing, N extra)。
#
# 修复（仅写 cr_root 下的 staging 副本，不动源数据）：
#   模式1（canon 双射，主路径）：在 tsv 符号序列上重现 R make.unique（与
#     make_10x_h5.py 同一算法），与 h5 name 集合做双射配对（X.N ↔ tsv 第 N 次
#     出现）；逐 h5 行写回 features/id←tsv ENSG、features/name←tsv 原始符号。
#     配对按 N 对齐保证 sctk 读取后 var_names_make_unique 与 bridge 的
#     make_index_unique 指向同一基因。
#   模式2（ENSG 直配）：h5 name 本身即 ENSG 时，按 ENSG 直接配对。
#   修复后 sctk.read_cellbender（var_names=name → make_unique）与 bridge 侧
#   集合必然相等（同一符号多重集）。矩阵数值/barcodes/droplet_latents 零改动。
#
# 用法:
#   06a_align_cellbender_features.py --cr-root <root> --samples <sid1,sid2,...>
#                                    [--check-only]
# 退出码: 0=全部一致或修复成功；1=任一样本不可对齐（不猜测、宁报错）
# 依赖: h5py/numpy/pandas（ENV_SCAUTOQC）
# ============================================================================
import argparse
import gzip
import os
import sys

import h5py
import numpy as np


def dash_unique(names):
    """与 anndata var_names_make_unique 等价：X, X-1, X-2 ...（按出现顺序）"""
    counts = {}
    out = []
    for n in names:
        k = counts.get(n, 0)
        out.append(n if k == 0 else f"{n}-{k}")
        counts[n] = k + 1
    return out


def r_make_unique(names):
    """R make.unique（点号风格，同 make_10x_h5.py）：第 k 次重复追加 .{k-1}"""
    counts = {}
    out = []
    for n in names:
        k = counts.get(n, 0)
        out.append(n if k == 0 else f"{n}.{k}")
        counts[n] = k + 1
    return out


def read_tsv_symbols(tsv_path):
    """features.tsv.gz → (ids, symbols)；空/NA 符号回退 id（对齐 make_10x_h5.py:56 与 Read10X 行为）"""
    with gzip.open(tsv_path, "rt") as f:
        rows = [l.rstrip("\n").split("\t") for l in f]
    ids = [r[0] for r in rows]
    syms = [r[1] if len(r) > 1 and r[1] not in ("", "NA") else r[0] for r in rows]
    return ids, syms


def read_h5_features(h5_path):
    with h5py.File(h5_path, "r") as f:
        feat = f["matrix/features"]
        name = [x.decode() for x in feat["name"][:]]
        hid = [x.decode() for x in feat["id"][:]] if "id" in feat else None
    return hid, name


def simulate_gather(tsv_syms, h5_names):
    """模拟 gather 两侧 var_names 集合比对：相等返回 None，否则 (missing, extra)"""
    a = set(dash_unique(tsv_syms))   # bridge 侧（sc.read_10x_mtx 默认唯一化）
    b = set(dash_unique(h5_names))   # sctk.read_cellbender 侧
    if a == b:
        return None
    return sorted(a - b), sorted(b - a)


def canon_bijection(tsv_syms, h5_names):
    """模式1：tsv 符号经 R make.unique 后与 h5 name 集合双射配对。
    成功返回 h5行号→tsv行号 的映射 list；失败返回 None。"""
    canon = r_make_unique(tsv_syms)  # tsv行号 → canon 名
    if sorted(canon) != sorted(h5_names):
        return None
    pos = {}
    for i, c in enumerate(canon):
        pos[c] = i  # canon 唯一（make.unique 保证），直接建索引
    try:
        return [pos[n] for n in h5_names]
    except KeyError:
        return None


def ensg_bijection(tsv_ids, h5_names):
    """模式2：h5 name 即 ENSG，直接配对。成功返回映射 list，失败 None。"""
    pos = {}
    for i, g in enumerate(tsv_ids):
        if g in pos:  # ENSG 应唯一；不唯一则此模式不可靠
            return None
        pos[g] = i
    if any(n not in pos for n in h5_names):
        return None
    return [pos[n] for n in h5_names]


def np_bytes(lst):
    return np.array([x.encode() for x in lst], dtype="S")


def process(cr_root, sid, check_only):
    tsv = os.path.join(cr_root, sid, "outs", "filtered_feature_bc_matrix", "features.tsv.gz")
    h5p = os.path.join(cr_root, sid, "outs", "cellbender", "cellbender_matrix.h5")
    for p in (tsv, h5p):
        if not os.path.isfile(p):
            print(f"[err] {sid}: 输入缺失 {p}")
            return 1
    tsv_ids, tsv_syms = read_tsv_symbols(tsv)
    h5_id, h5_names = read_h5_features(h5p)
    if len(h5_names) != len(tsv_ids):
        print(f"[err] {sid}: 行数不符 h5={len(h5_names)} tsv={len(tsv_ids)}")
        return 1

    diff = simulate_gather(tsv_syms, h5_names)
    if diff is None:
        print(f"[ok] {sid}: features 集合一致（{len(tsv_ids)}），无需修复")
        return 0
    missing, extra = diff
    print(f"[chk] {sid}: {len(missing)} missing / {len(extra)} extra（例: {missing[:3]} vs {extra[:3]}）")

    mapping = canon_bijection(tsv_syms, h5_names)
    mode = "canon(R点号双射)"
    if mapping is None:
        mapping = ensg_bijection(tsv_ids, h5_names)
        mode = "ensg(直配)"
    if mapping is None:
        print(f"[err] {sid}: 两种配对模式均失败，不猜测——请人工核查该样本 features 来源")
        return 1
    if check_only:
        print(f"[chk] {sid}: 可经 {mode} 修复（--check-only 未写入）")
        return 0

    new_id = [tsv_ids[i] for i in mapping]
    new_nm = [tsv_syms[i] for i in mapping]
    with h5py.File(h5p, "r+") as f:
        feat = f["matrix/features"]
        old_attrs = {k: dict(feat[k].attrs) for k in ("id", "name") if k in feat}
        if "id" in feat:
            del feat["id"]
        if "name" in feat:
            del feat["name"]
        for key, vals in (("id", new_id), ("name", new_nm)):
            d = feat.create_dataset(key, data=np_bytes(vals), compression="gzip", compression_opts=6, shuffle=True, chunks=(len(vals),))
            for ak, av in old_attrs.get(key, {}).items():
                d.attrs[ak] = av
    # 回读复验：模拟 gather 必须一致
    rb_id, rb_nm = read_h5_features(h5p)
    diff2 = simulate_gather(tsv_syms, rb_nm)
    if diff2 is not None or rb_id != new_id:
        print(f"[err] {sid}: 修复后复验失败")
        return 1
    n_remap = sum(1 for a, b in zip(h5_names, new_nm) if a != b)
    print(f"[fix] {sid}: {mode} 对齐 {len(new_id)} features（重命名 {n_remap} 行），id→ENSG, name→tsv 符号，复验一致")
    return 0


def main():
    ap = argparse.ArgumentParser(description="cellbender h5 与 filtered mtx features 约定一致性检查+修复")
    ap.add_argument("--cr-root", required=True)
    ap.add_argument("--samples", required=True, help="逗号分隔 sid 列表")
    ap.add_argument("--check-only", action="store_true", help="只检查不写入")
    args = ap.parse_args()
    rc = 0
    for sid in [s.strip() for s in args.samples.split(",") if s.strip()]:
        rc |= process(args.cr_root, sid, args.check_only)
    sys.exit(rc)


if __name__ == "__main__":
    main()
