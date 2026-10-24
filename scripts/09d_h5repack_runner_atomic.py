#!/usr/bin/env python3
# ========================================================================
# 存量 h5ad 原子压缩 runner：h5repack -f GZIP=4 生成 ASCII 暂存临时件 → 逐格
# 校验 → 硬链接保全原件 → 单次 os.replace 原子换名；清单驱动、断点续跑，
# 原名全程不缺失、无半成品窗口
# 用法: 09d_h5repack_runner_atomic.py <manifest.csv> [backup_root] [staging_dir]
#   manifest 列: path,status；仅 status=pending 参与，按 OK / OK_resume /
#   FAIL_* / MISSING 原地写回
# 依赖: 同目录 09e_verify_h5ad_equal.py（须先经阳性/阴性回归；判据见
#   references/delivery-verification.md）；h5repack 可执行文件（环境变量
#   H5REPACK 或 PATH 之一）
# 注意: 传入 h5repack 的路径必须纯 ASCII（源目录与暂存目录均保持 ASCII）
# ========================================================================
import sys
import os
import csv
import subprocess
import shutil
import datetime

H5REPACK = os.environ.get("H5REPACK") or shutil.which("h5repack")
HERE = os.path.dirname(os.path.abspath(__file__))
VERIFY = os.path.join(HERE, "09e_verify_h5ad_equal.py")
PY = sys.executable


def log(logcsv, row):
    new = not os.path.exists(logcsv)
    with open(logcsv, "a", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        if new:
            w.writerow(["time", "path", "status", "old_bytes", "new_bytes", "stage", "note"])
        w.writerow(row)


def run(cmd):
    p = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
    return p.returncode, ((p.stdout or "") + (p.stderr or ""))


def write_manifest(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    if not H5REPACK:
        sys.exit("h5repack not found: set env H5REPACK to the executable path")
    man = sys.argv[1]
    workdir = os.path.dirname(os.path.abspath(man))
    backup_root = sys.argv[2] if len(sys.argv) > 2 else os.path.join(workdir, "h5ad_orig_backup")
    staging = sys.argv[3] if len(sys.argv) > 3 else os.path.join(workdir, "_tmp_gz")
    os.makedirs(staging, exist_ok=True)
    logcsv = os.path.join(workdir, "h5ad_convert_log.csv")
    rows = list(csv.DictReader(open(man, encoding="utf-8")))

    def now():
        return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    for i, row in enumerate(rows):
        if row.get("status") != "pending":
            continue
        src = row["path"]
        ds = os.path.basename(os.path.dirname(src))
        bk = os.path.join(backup_root, ds, os.path.basename(src))
        tmp = os.path.join(staging, ds + "__" + os.path.basename(src) + ".gztmp")
        if not os.path.exists(src):
            row["status"] = "MISSING"
            log(logcsv, [now(), src, "MISSING", 0, 0, "check", "file not found"])
            write_manifest(man, rows)
            continue
        if os.path.exists(bk):
            if os.path.exists(tmp):
                # crashed between hardlink and replace: tmp ready, bk links the still-intact original
                rc, out = run([PY, VERIFY, src, tmp])
                if rc == 0:
                    os.replace(tmp, src)
                    row["status"] = "OK_resume"
                    log(logcsv, [now(), src, "OK_resume", os.path.getsize(bk), os.path.getsize(src), "resume", "pre-replace crash recovered | bk(hardlink)=" + bk])
                else:
                    row["status"] = "FAIL_verify"
                    log(logcsv, [now(), src, "FAIL_verify", os.path.getsize(src), os.path.getsize(tmp), "resume", "re-verify failed | bk=" + bk])
            else:
                # crashed after replace (src already new); bk links the original
                rc, out = run([PY, VERIFY, bk, src])
                if rc == 0:
                    row["status"] = "OK_resume"
                    log(logcsv, [now(), src, "OK_resume", os.path.getsize(bk), os.path.getsize(src), "resume", "post-replace state confirmed | bk(hardlink)=" + bk])
                else:
                    row["status"] = "FAIL_precheck"
                    log(logcsv, [now(), src, "FAIL_precheck", os.path.getsize(src), 0, "check", "bk exists but src/bk verify failed | bk=" + bk])
            write_manifest(man, rows)
            continue
        old = os.path.getsize(src)
        if os.path.exists(tmp):
            os.remove(tmp)
        t0 = now()
        rc, out = run([H5REPACK, "-f", "GZIP=4", src, tmp])
        if rc != 0 or not os.path.exists(tmp):
            row["status"] = "FAIL_repack"
            log(logcsv, [t0, src, "FAIL_repack", old, 0, "h5repack", out[-300:]])
            if os.path.exists(tmp):
                os.remove(tmp)
            write_manifest(man, rows)
            continue
        rc, out = run([PY, VERIFY, src, tmp])
        if rc != 0:
            row["status"] = "FAIL_verify"
            log(logcsv, [t0, src, "FAIL_verify", old, os.path.getsize(tmp), "verify", out[-300:].replace("\n", " | ")])
            write_manifest(man, rows)
            continue
        new = os.path.getsize(tmp)
        os.makedirs(os.path.dirname(bk), exist_ok=True)
        os.link(src, bk)
        os.replace(tmp, src)
        row["status"] = "OK"
        log(logcsv, [t0, src, "OK", old, new, "convert", out.strip().replace("\n", " | ") + " | bk(hardlink)=" + bk])
        write_manifest(man, rows)
        print("[%d/%d] OK %s | %.2fGB -> %.2fGB | bk=%s" % (i + 1, len(rows), os.path.basename(src), old / 1e9, new / 1e9, bk), flush=True)
    print("H5AD-ATOMIC-RUNNER-DONE", flush=True)


if __name__ == "__main__":
    main()
