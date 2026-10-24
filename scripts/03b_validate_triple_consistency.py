#!/usr/bin/env python3
"""严格校验 filtered mtx、10x formatted h5 与 velocyto loom 三位一体一致性。

mtx ↔ formatted h5（同源转换，硬等价）：基因行顺序、barcode 列顺序、
稀疏表达值三元组必须逐值一致。h5 行名兼容两种语义：features/id 全为
ENSG/ENSMUSG（10x 官方/CellRanger 原生）时按 ENSG 与 mtx 第 1 列比对；
否则按符号语义（make_10x_h5 构建 / write10xCounts 语义：id=name=make.unique(符号)）与
mtx 第 2 列经同一 make.unique 规则比对，并要求 h5 的 id 与 name 两列逐元素相同。
loom（velocyto 剪接分层计数，语义不同于 GeneFull）：与 mtx/h5 坐标系一一对应——
基因（ENSG）集合与 barcode（归一后）集合完全一致、各自无重复、shape 一致；
行列顺序可与比对矩阵不同（velocyto 按自身 GTF 解析与 BAM 顺序输出），由下游
bridge 按 ENSG+barcode 严格重排对齐，不允许去重/补零/符号对齐等任何自定义改动。
loom 层契约：spliced/unspliced 层存在且非空（主层为空占位，不参与）。
loom 若另含显式 GeneFull 层，则按坐标映射重排后与 mtx/h5 逐值比较（加严通道）。
关系式 counts = spliced+unspliced+ambiguous 对 GeneFull 不成立（velocyto 计数
模型与比对器不同：-m mask 重复区、分子须可明确分类、multi-gene 丢弃），不做
逐值比较；s+u+a 与 mtx 总量比仅作信息性输出，不参与判定。
输入只读；任一缺失或不一致均以非零退出。
"""
import argparse
import hashlib
import gzip
import os
import re
import sys

import h5py
import numpy as np
import scipy.sparse as sp

from make_10x_h5 import make_unique

ENSG_RE = re.compile(r"^ENS[A-Z]*G\d+")  # ENSG…/ENSMUSG…（人/鼠 Ensembl 基因 ID）


def _text(x):
    return x.decode("utf-8", "replace") if isinstance(x, (bytes, np.bytes_)) else str(x)


def read_lines(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as f:
        return [line.rstrip("\r\n") for line in f]


def read_features(path):
    rows = [line.split("\t") for line in read_lines(path)]
    if not rows or any(len(r) < 1 for r in rows):
        raise ValueError(f"features 文件为空或格式错误: {path}")
    return [_text(r[0]) for r in rows], [_text(r[1]) if len(r) > 1 else "" for r in rows]


def read_mtx(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as f:
        for line in f:
            if not line.startswith("%"):
                nrow, ncol, nnz = map(int, line.split()[:3]); break
        rows = np.empty(nnz, dtype=np.int64); cols = np.empty(nnz, dtype=np.int64); vals = np.empty(nnz, dtype=np.float64)
        i = 0
        for line in f:
            if not line.strip(): continue
            r, c, v = line.split()[:3]
            if i >= nnz: raise ValueError("matrix.mtx 非零条目超过 header 声明")
            rows[i], cols[i], vals[i] = int(r)-1, int(c)-1, float(v); i += 1
    if i != nnz: raise ValueError(f"matrix.mtx 非零条目不足: {i}/{nnz}")
    return sp.coo_matrix((vals, (rows, cols)), shape=(nrow, ncol))


def read_h5(path):
    with h5py.File(path, "r") as h:
        g = h["matrix"]
        bc = [_text(x) for x in g["barcodes"][:]]
        fg = g["features"]
        ids = [_text(x) for x in fg["id"][:]]
        name_ds = fg["name"] if "name" in fg else fg["gene_names"]
        names = [_text(x) for x in name_ds[:]]
        m = sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(g["shape"][:]))
    return bc, ids, names, m


def norm_barcode(x):
    b = _text(x).strip()
    if ":" in b: b = b.rsplit(":", 1)[-1]
    if b.lower().endswith("x"): b = b[:-1]
    if b.endswith("-1"): b = b[:-2]
    return b.upper()


def ensure_unique(label, values, errors):
    dup = len(values) - len(set(values))
    if dup: errors.append(f"{label}: 归一化后有 {dup} 个重复值")


def _finish_hash(shape, nnz, stream_digest):
    h = hashlib.blake2b(digest_size=32)
    h.update(np.asarray([shape[0], shape[1], nnz], dtype="<i8").tobytes())
    h.update(stream_digest)
    return h.hexdigest()


def hash_sparse_coo(m):
    # Canonicalize duplicate coordinates before hashing (Matrix Market may contain them).
    c = m.tocoo(copy=True); c.sum_duplicates()
    order = np.lexsort((c.row, c.col))  # column-major canonical order
    stream = hashlib.blake2b(digest_size=32)
    for idx in order:
        stream.update(np.asarray([c.row[idx], c.col[idx]], dtype="<i8").tobytes())
        stream.update(np.asarray(c.data[idx], dtype="<f8").tobytes())
    return _finish_hash(m.shape, len(order), stream.digest()), len(order)


def hash_loom_layer(layer, shape, row_to_mtx=None, col_to_mtx=None, chunk=512):
    # Bounded column chunks; collect triplets, map to mtx coordinates when the loom
    # ordering differs, then hash in global column-major canonical order.
    rs, cs, vs = [], [], []
    for start in range(0, shape[1], chunk):
        a = np.asarray(layer[:, start:min(start + chunk, shape[1])])
        rr, cc = np.nonzero(a)
        vv = a[rr, cc].astype(np.float64)
        mr = row_to_mtx[rr] if row_to_mtx is not None else rr.astype(np.int64)
        mc = col_to_mtx[start + cc] if col_to_mtx is not None else (start + cc).astype(np.int64)
        rs.append(mr.astype(np.int64)); cs.append(mc.astype(np.int64)); vs.append(vv)
    rr = np.concatenate(rs); cc = np.concatenate(cs); vv = np.concatenate(vs)
    order = np.lexsort((rr, cc))
    stream = hashlib.blake2b(digest_size=32)
    for idx in order:
        stream.update(np.asarray([rr[idx], cc[idx]], dtype="<i8").tobytes())
        stream.update(np.asarray(vv[idx], dtype="<f8").tobytes())
    return _finish_hash(shape, len(order), stream.digest()), len(order)


def loom_layer_stats(layer, shape, chunk=512):
    # Bounded column chunks; count nonzeros and total without materializing the layer.
    nnz = 0; total = 0.0
    for start in range(0, shape[1], chunk):
        a = np.asarray(layer[:, start:min(start + chunk, shape[1])], dtype=np.float64)
        nnz += int(np.count_nonzero(a)); total += float(a.sum())
    return nnz, total


def read_loom_meta(path):
    try:
        import loompy
    except ImportError as e:
        raise RuntimeError("缺少 loompy，无法只读检查 loom") from e
    ds = loompy.connect(path, mode="r", validate=False)
    try:
        acc = [_text(x) for x in ds.ra["Accession"]]
        cells = [_text(x) for x in ds.ca["CellID"]]
        layers = list(ds.layers.keys())
        stats = {}
        for name in ("spliced", "unspliced", "ambiguous"):
            if name in layers:
                stats[name] = loom_layer_stats(ds.layers[name], ds.shape)
        genefull = [k for k in layers if str(k).strip().lower() == "genefull"]
        return acc, cells, layers, stats, (genefull[0] if genefull else None), tuple(ds.shape)
    finally:
        ds.close()


def hash_loom_genefull(path, layer_name, shape, row_to_mtx, col_to_mtx):
    import loompy
    ds = loompy.connect(path, mode="r", validate=False)
    try:
        return hash_loom_layer(ds.layers[layer_name], shape, row_to_mtx, col_to_mtx)
    finally:
        ds.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sid", required=True); ap.add_argument("--cr-root", required=True)
    args = ap.parse_args(); root = os.path.abspath(args.cr_root); sid = args.sid
    outs = os.path.join(root, sid, "outs"); errors = []
    mtx_dir = next((os.path.join(outs, x) for x in ("filtered_cell_gene_matrix", "filtered_feature_bc_matrix") if os.path.isdir(os.path.join(outs, x))), None)
    if not mtx_dir: errors.append(f"缺少 mtx 目录: {outs}")
    velo = os.path.join(root, sid, "velocyto", f"{sid}.loom")
    h5_path = os.path.join(outs, "filtered_feature_bc_matrix.h5")
    for p, label in ((velo, "loom"), (h5_path, "formatted h5")):
        if not os.path.isfile(p): errors.append(f"缺少 {label}: {p}")
    if errors:
        print("[FAIL] 三位一体输入缺失"); [print("  -", e) for e in errors]; return 1
    try:
        mtx_bc = read_lines(os.path.join(mtx_dir, "barcodes.tsv.gz")); mtx_id, mtx_name = read_features(os.path.join(mtx_dir, "features.tsv.gz")); mtx = read_mtx(os.path.join(mtx_dir, "matrix.mtx.gz"))
        h5_bc, h5_id, h5_name, h5m = read_h5(h5_path)
        loom_id, loom_cell, layers, layer_stats, genefull_layer, loom_shape = read_loom_meta(velo)
    except Exception as e:
        print(f"[FAIL] 读取输入失败: {e}"); return 1
    for label, ids in (("mtx", mtx_id), ("loom", loom_id)):
        bad = [x for x in ids if not ENSG_RE.match(x)]
        if bad: errors.append(f"{label} gene ID 非 ENSG/ENSMUSG: {len(bad)} 个（前5: {bad[:5]}）")

    # ---- mtx ↔ h5（同源转换）：行序/列序/表达值逐值硬等价 ----
    if all(ENSG_RE.match(x) for x in h5_id):
        # ENSG 语义 h5（CellRanger 原生/10x 官方）：features/id 与 mtx 第 1 列逐序一致
        if mtx_id != h5_id: errors.append("mtx 与 h5 基因 ID 顺序不一致（features 第1列 vs features/id）")
    else:
        # 符号语义 h5（make_10x_h5 构建 / write10xCounts 语义）：id=name=make.unique(符号)，空符号回退第 1 列
        if h5_id != h5_name: errors.append("h5 features/id 与 features/name 不一致（符号语义要求两列逐元素相同）")
        sym = [n if n not in ("", "NA") else i for i, n in zip(mtx_id, mtx_name)]
        if make_unique(sym) != h5_id: errors.append("h5 行名与 mtx 符号列（make.unique 去重后）不一致")
    mb, hb = [norm_barcode(x) for x in mtx_bc], [norm_barcode(x) for x in h5_bc]
    lb = [norm_barcode(x) for x in loom_cell]
    for label, vals in (("mtx", mb), ("h5", hb), ("loom", lb)): ensure_unique(label, vals, errors)
    if mb != hb: errors.append("mtx 与 h5 barcode 顺序不一致")
    if mtx.shape != h5m.shape: errors.append(f"mtx 与 h5 shape 不一致: {mtx.shape} vs {h5m.shape}")
    md, mn = hash_sparse_coo(mtx); hd, hn = hash_sparse_coo(h5m)
    if (mn, md) != (hn, hd): errors.append(f"mtx 与 h5 稀疏表达值不一致（nnz {mn} vs {hn}）")

    # ---- loom：坐标系集合级一一对应（顺序可不同，由下游 bridge 按 ENSG+barcode 重排）----
    if len(loom_id) != len(set(loom_id)):
        errors.append(f"loom 基因 ID 有重复: {len(loom_id) - len(set(loom_id))} 个")
    elif set(loom_id) != set(mtx_id):
        diff = set(mtx_id) ^ set(loom_id)
        errors.append(f"loom 与 mtx 基因集合不一致: {len(diff)} 个差异（前5: {sorted(diff)[:5]}）")
    if set(lb) != set(mb):
        diff = set(mb) ^ set(lb)
        errors.append(f"loom 与 mtx barcode 集合不一致: {len(diff)} 个差异（前5: {sorted(diff)[:5]}）")
    if loom_shape != mtx.shape:
        errors.append(f"loom 与 mtx shape 不一致: loom={loom_shape} vs mtx={mtx.shape}")

    # ---- loom 层契约：velocyto 标准层存在且非空（主层为空占位，不参与）----
    missing_layers = [x for x in ("spliced", "unspliced") if x not in layers]
    if missing_layers:
        errors.append(f"loom 缺少 velocyto 标准层 {missing_layers}（现有层: {layers}）")
    else:
        for name in ("spliced", "unspliced"):
            if layer_stats.get(name, (0, 0.0))[0] == 0:
                errors.append(f"loom {name} 层全零（nnz=0），不是有效的 velocyto 产出")

    # ---- 加严通道：loom 含显式 GeneFull 层时，按坐标映射重排后与 mtx/h5 逐值比较 ----
    genefull_note = ""
    if genefull_layer is not None and not errors:
        gpos = {g: i for i, g in enumerate(mtx_id)}
        cpos = {b: j for j, b in enumerate(mb)}
        row_map = np.array([gpos[g] for g in loom_id]) if loom_id != mtx_id else None
        col_map = np.array([cpos[b] for b in lb]) if lb != mb else None
        loom_digest, loom_nnz = hash_loom_genefull(velo, genefull_layer, loom_shape, row_map, col_map)
        if (mn, md) != (loom_nnz, loom_digest):
            errors.append(f"loom GeneFull 层与 mtx/h5 表达值不一致（nnz {loom_nnz} vs {mn}）")
        else:
            genefull_note = f", GeneFull nnz={loom_nnz}（与 mtx 逐值一致）"

    if errors:
        print(f"[FAIL] 三位一体校验不一致（{sid}）:"); [print("  -", e) for e in errors]; return 1

    # ---- 信息性输出（不参与判定）----
    same_order = loom_id == mtx_id and lb == mb
    sua_total = sum(s[1] for s in layer_stats.values())
    sua_ratio = sua_total / mtx.data.sum() if mtx.data.sum() else float("nan")
    print(f"[OK] 三位一体一致: {sid} | genes={len(mtx_id)} cells={len(mb)} mtx/h5 nnz={mn}")
    print(f"     loom 行列序与 mtx 相同: {same_order}（不同则由下游 bridge 按 ENSG+barcode 重排）")
    print(f"     loom 层: spliced nnz={layer_stats.get('spliced', (0,0))[0]}, unspliced nnz={layer_stats.get('unspliced', (0,0))[0]}{genefull_note}")
    print(f"     [INFO] sum(spliced+unspliced+ambiguous)/sum(mtx) = {sua_ratio:.4f}（velocyto 计数模型系统性低于 GeneFull，比值仅作健康度参考，不参与判定）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
