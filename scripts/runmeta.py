#!/usr/bin/env python3
# ============================================================================
# _runmeta 指纹 manifest —— 断点续跑的 skip 判据（upstream 指纹判据 / formatted F1）
#
# 判据是"输入指纹逐项一致"，而非"产物文件存在"：
#   <CR_ROOT>/_runmeta/<sid>.<stage>.fingerprint.json   每样本阶段（stage=cellbender|scautoqc）
#   <CR_ROOT>/_runmeta/merged_<对象>.fingerprint.json    合并级产物（成员指纹集合 + 合并参数）
#   字段：script（脚本名@sha256[:12]，可多个）/ params（运行参数）/ marker_db（marker 库 sha256）
#         / inputs（相对路径 → size + mtime_ns）/ versions（本阶段相关工具版本）
#         / members（合并级：成员指纹文件 sha256）/ produced_at
# 任一字段不等 → 该范围失效须重跑（marker 库、质控参数、参考、脚本版本、输入文件、样本集合任一变化；
# 语义等价的中断续传不受影响）。已存在的产物若无指纹记录，视为"无法确证来源"→ 按需重跑。
#
# CLI（bash 直调；python 脚本可 sys.path 引入 import）：
#   runmeta.py build --fp <path> --script-file <脚本> [--script-file <脚本>…]
#                    [--set k=v …] [--marker-sha <hex>] [--input <文件> …]
#                    [--version k=v …] [--member <指纹json> …]
#                    [--base <相对化根>]
#   runmeta.py check <同 build 参数> [--require <产物> …]
#                    退出码 0=指纹一致可 skip；1=失效/无记录/必需产物缺失（须重跑）
#   runmeta.py verify --dir <_runmeta> --samples s1,s2 [--branch scautoqc]
#                    校验指纹记录齐备/字段完整/合并级成员一致；0=通过；1=缺失明细（打印）
# 依赖：仅 Python 标准库（任何 conda env 的 python 均可运行）。
# 本文件随仓库发布；若被其他流程复用，以本副本为准。
# ============================================================================
import argparse
import hashlib
import json
import os
import sys
import time

FIELDS = ("script", "params", "marker_db", "inputs", "versions", "produced_at")


def sha256_file(path, n=12):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:n]


def canonical(obj):
    return json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":"))


def make_fingerprint(script_files=None, sets=None, marker_sha="", inputs=None,
                     versions=None, members=None, base=None):
    """构造指纹 dict（inputs 缺失即抛 FileNotFoundError，调用方转成"须重跑"）"""
    scripts = ["%s@%s" % (os.path.basename(p), sha256_file(p)) for p in (script_files or [])]
    fp = {
        "script": scripts,
        "params": {str(k): str(v) for k, v in (sets or {}).items()},
        "marker_db": marker_sha or "",
        "inputs": {},
        "versions": {str(k): str(v) for k, v in (versions or {}).items()},
        "produced_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    for p in (inputs or []):
        st = os.stat(p)
        key = os.path.relpath(p, base).replace(os.sep, "/") if base else os.path.abspath(p).replace(os.sep, "/")
        fp["inputs"][key] = {"size": st.st_size, "mtime_ns": st.st_mtime_ns}
    if members:
        fp["members"] = {}
        for m in members:
            name = os.path.basename(m)
            fp["members"][name] = sha256_file(m, 64) if os.path.isfile(m) else "MISSING"
    return fp


def compare(old, new):
    """返回差异说明列表（忽略 produced_at）；空列表=逐项一致"""
    diffs = []
    ok = dict(old); ok.pop("produced_at", None)
    nk = dict(new); nk.pop("produced_at", None)
    for k in sorted(set(ok) | set(nk)):
        ov, nv = ok.get(k), nk.get(k)
        if ov == nv:
            continue
        if k in ("inputs", "members"):
            for f in sorted(set(ov or {}) | set(nv or {})):
                a, b = (ov or {}).get(f), (nv or {}).get(f)
                if a != b:
                    diffs.append("%s[%s]: %s -> %s" % (k, f, a, b))
        else:
            diffs.append("%s: %r -> %r" % (k, ov, nv))
    return diffs


def check(fp_path, require=None, **kw):
    """返回 (rc, msg)：rc=0 指纹一致可 skip；1 须重跑"""
    missing = [p for p in (require or []) if not os.path.exists(p)]
    if missing:
        msg = "[fingerprint] 必需产物缺失 -> 重跑: %s" % missing[0]
        print(msg)
        return 1, msg
    if not os.path.isfile(fp_path):
        msg = "[fingerprint] 无指纹记录（首次运行/无来源可确证）-> 重跑: %s" % fp_path
        print(msg)
        return 1, msg
    try:
        with open(fp_path, encoding="utf-8") as f:
            old = json.load(f)
    except Exception as e:
        msg = "[fingerprint] 记录读取失败 -> 重跑: %s" % e
        print(msg)
        return 1, msg
    try:
        new = make_fingerprint(**kw)
    except FileNotFoundError as e:
        msg = "[fingerprint] 输入缺失 -> 重跑: %s" % e
        print(msg)
        return 1, msg
    diffs = compare(old, new)
    if diffs:
        msg = "[fingerprint] 指纹变化 -> 重跑: " + "; ".join(diffs[:8])
        print(msg)
        return 1, msg
    print("[fingerprint] 一致（skip）: %s" % fp_path)
    return 0, ""


def write(fp_path, **kw):
    fp = make_fingerprint(**kw)
    os.makedirs(os.path.dirname(os.path.abspath(fp_path)), exist_ok=True)
    with open(fp_path, "w", encoding="utf-8") as f:
        json.dump(fp, f, ensure_ascii=False, indent=2, sort_keys=True)
        f.write("\n")
    print("[fingerprint] 写入: %s" % fp_path)
    return fp


def verify(rm_dir, samples, branch=""):
    """交付校验：指纹记录齐备 + 字段完整 + 合并级成员一致。返回 0/1（明细打印）"""
    errs = []

    def load_checked(path):
        try:
            with open(path, encoding="utf-8") as f:
                fp = json.load(f)
        except Exception as e:
            errs.append("%s 解析失败: %s" % (os.path.basename(path), e))
            return None
        miss = [k for k in FIELDS if k not in fp]
        if miss:
            errs.append("%s 缺字段 %s" % (os.path.basename(path), miss))
        return fp

    stage = "scautoqc" if branch == "scautoqc" else ""
    expect = ["%s.cellbender.fingerprint.json" % s for s in samples]
    if stage:
        expect += ["%s.%s.fingerprint.json" % (s, stage) for s in samples]
    merged = []
    if branch == "scautoqc" and len(samples) > 1:
        merged.append("merged_scautoqc_pooled_filtered.fingerprint.json")

    for name in expect + merged:
        p = os.path.join(rm_dir, name)
        if not os.path.isfile(p):
            errs.append("缺指纹记录 _runmeta/%s" % name)
            continue
        fp = load_checked(p)
        if fp is None or name not in merged:
            continue
        for m, digest in sorted((fp.get("members") or {}).items()):
            if digest == "MISSING":
                errs.append("_runmeta/%s: 成员 %s 构建时无指纹记录" % (name, m))
                continue
            mp = os.path.join(rm_dir, m)
            if not os.path.isfile(mp):
                errs.append("_runmeta/%s: 成员指纹 %s 缺失" % (name, m))
            elif sha256_file(mp, 64) != digest:
                errs.append("_runmeta/%s: 成员 %s 指纹已变化（合并级产物未随成员更新）" % (name, m))
    if errs:
        for e in errs:
            print("[fingerprint] " + e)
        return 1
    print("[fingerprint] 指纹记录完整（%d 样本，branch=%s）" % (len(samples), branch or "auto"))
    return 0


def _kv(items):
    out = {}
    for it in (items or []):
        if "=" not in it:
            raise SystemExit("[err] --set/--version 需要 k=v 形式: %s" % it)
        k, v = it.split("=", 1)
        out[k] = v
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)
    for m in ("build", "check"):
        p = sub.add_parser(m)
        p.add_argument("--fp", required=True)
        p.add_argument("--script-file", action="append", default=[])
        p.add_argument("--set", action="append", default=[])
        p.add_argument("--marker-sha", default="")
        p.add_argument("--input", action="append", default=[])
        p.add_argument("--version", action="append", default=[])
        p.add_argument("--member", action="append", default=[])
        p.add_argument("--base", default="")
        p.add_argument("--require", action="append", default=[])
    pv = sub.add_parser("verify")
    pv.add_argument("--dir", required=True)
    pv.add_argument("--samples", required=True)
    pv.add_argument("--branch", default="")
    a = ap.parse_args(argv)

    if a.mode == "verify":
        return verify(a.dir, [s.strip() for s in a.samples.split(",") if s.strip()], a.branch)

    kw = dict(script_files=a.script_file, sets=_kv(a.set), marker_sha=a.marker_sha,
              inputs=a.input, versions=_kv(a.version), members=a.member,
              base=(a.base or None))
    if a.mode == "build":
        write(a.fp, **kw)
        return 0
    rc, _ = check(a.fp, require=a.require, **kw)
    return rc


if __name__ == "__main__":
    sys.exit(main())
