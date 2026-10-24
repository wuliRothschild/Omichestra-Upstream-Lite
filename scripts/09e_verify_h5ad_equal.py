#!/usr/bin/env python3
# ========================================================================
# h5ad 逐格等价校验器：结构（组 / 数据集 / 属性）+ 全 dataset 值流式比对；只读
# 用法: 09e_verify_h5ad_equal.py <src.h5ad> <dst.h5ad>
# 退出码: 0 = 逐格全等；非 0 = 存在差异或参数错误
# ========================================================================
import sys
import numpy as np
import h5py

CHUNK = 8_000_000  # elements per compare step


def arr_equal(a, b):
    if not isinstance(a, np.ndarray):
        a = np.asarray(a)
    if not isinstance(b, np.ndarray):
        b = np.asarray(b)
    if a.shape != b.shape:
        return False
    if a.dtype.kind == "O" or b.dtype.kind == "O":
        try:
            for x, y in zip(a.ravel().tolist(), b.ravel().tolist()):
                if x == y:
                    continue
                try:
                    if np.isnan(x) and np.isnan(y):
                        continue
                except Exception:
                    pass
                return False
            return True
        except Exception:
            return False
    try:
        if a.dtype.kind == "f" and b.dtype.kind == "f":
            return bool(np.array_equal(a, b, equal_nan=True))
        return bool(np.array_equal(a, b))
    except Exception:
        return False


def attrs_equal(ga, gb, path, problems):
    ka, kb = set(ga.keys()), set(gb.keys())
    if ka != kb:
        problems.append("ATTR-KEYS %s: only-src=%s only-dst=%s" % (path, sorted(ka - kb), sorted(kb - ka)))
        return
    for k in ka:
        va, vb = ga[k], gb[k]
        va = va.tolist() if isinstance(va, np.ndarray) else va
        vb = vb.tolist() if isinstance(vb, np.ndarray) else vb
        if isinstance(va, float) and isinstance(vb, float):
            ok = (va == vb) or (np.isnan(va) and np.isnan(vb))
        else:
            try:
                ok = (va == vb)
                if isinstance(ok, np.ndarray):
                    ok = bool(ok.all())
            except Exception:
                ok = False
        if not ok:
            problems.append("ATTR %s/%s: %r vs %r" % (path, k, va, vb))


def compare_group(ga, gb, name, problems, ds_count):
    attrs_equal(ga.attrs, gb.attrs, name, problems)
    ka, kb = set(ga.keys()), set(gb.keys())
    if ka != kb:
        problems.append("KEYS %s: only-src=%s only-dst=%s" % (name, sorted(ka - kb), sorted(kb - ka)))
    for k in sorted(ka & kb):
        ca, cb = ga[k], gb[k]
        sub = name + "/" + k
        if isinstance(ca, h5py.Group):
            if not isinstance(cb, h5py.Group):
                problems.append("KIND %s: group vs dataset" % sub)
                continue
            ds_count = compare_group(ca, cb, sub, problems, ds_count)
            continue
        if isinstance(cb, h5py.Group):
            problems.append("KIND %s: dataset vs group" % sub)
            continue
        ds_count += 1
        if ca.dtype != cb.dtype:
            problems.append("DTYPE %s: %s vs %s" % (sub, ca.dtype, cb.dtype))
            continue
        if ca.shape != cb.shape:
            problems.append("SHAPE %s: %s vs %s" % (sub, ca.shape, cb.shape))
            continue
        attrs_equal(ca.attrs, cb.attrs, sub, problems)
        if len(ca.shape) == 0:
            if not arr_equal(ca[()], cb[()]):
                problems.append("VALUE %s: scalar differs" % sub)
            continue
        per_row = int(np.prod(ca.shape[1:])) if ca.ndim > 1 else 1
        step = max(1, CHUNK // max(1, per_row))
        rows = ca.shape[0]
        i = 0
        while i < rows:
            j = min(rows, i + step)
            if not arr_equal(ca[i:j], cb[i:j]):
                problems.append("VALUE %s: mismatch in rows [%d:%d]" % (sub, i, j))
                break
            i = j
    return ds_count


def main():
    if len(sys.argv) != 3:
        print("usage: 09e_verify_h5ad_equal.py <src.h5ad> <dst.h5ad>")
        sys.exit(2)
    problems = []
    with h5py.File(sys.argv[1], "r") as fa, h5py.File(sys.argv[2], "r") as fb:
        attrs_equal(fa.attrs, fb.attrs, "/", problems)
        if set(fa.keys()) != set(fb.keys()):
            problems.append("ROOT KEYS: %s vs %s" % (sorted(fa.keys()), sorted(fb.keys())))
        ds_count = 0
        for k in sorted(set(fa.keys()) & set(fb.keys())):
            oa, ob = fa[k], fb[k]
            if isinstance(oa, h5py.Group):
                if not isinstance(ob, h5py.Group):
                    problems.append("ROOT KIND %s" % k)
                    continue
                ds_count = compare_group(oa, ob, k, problems, ds_count)
            else:
                ds_count += 1
                if oa.dtype != ob.dtype or oa.shape != ob.shape:
                    problems.append("ROOT DS %s dtype/shape" % k)
                    continue
                if not arr_equal(oa[()], ob[()]):
                    problems.append("ROOT DS %s values" % k)
    print("VERIFY datasets=%d problems=%d" % (ds_count, len(problems)))
    for p in problems[:60]:
        print("  " + p)
    sys.exit(0 if not problems else 3)


if __name__ == "__main__":
    main()
