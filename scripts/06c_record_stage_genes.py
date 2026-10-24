#!/usr/bin/env python3
"""把质控阶段实际使用的基因信息写入每样本 outs/metadata.json 的 stages 键。

本流程（scautoqc）不启用 celltypist → stages.scautoqc.celltypist_model = null（如实记录，
便于下游追溯这一阶段有没有用到参考基因集）。

幂等：只更新 stages 键，十类基础字段（01c 写入）不动。
用法:
  06c_record_stage_genes.py --cr-root <根> --samples sid1,sid2
"""
import argparse
import json
import os


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cr-root", required=True)
    ap.add_argument("--samples", required=True)
    args = ap.parse_args()

    stages = {"scautoqc": {"celltypist_model": None}}
    for sid in [s.strip() for s in args.samples.split(",") if s.strip()]:
        path = os.path.join(args.cr_root, sid, "outs", "metadata.json")
        if not os.path.isfile(path):
            print(f"[warn] {sid}: 缺 outs/metadata.json，跳过")
            continue
        with open(path, encoding="utf-8") as f:
            meta = json.load(f)
        meta.setdefault("stages", {}).update(stages)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(meta, f, ensure_ascii=False, indent=2)
            f.write("\n")
        print(f"[metadata] {sid} stages 更新（scautoqc）")


if __name__ == "__main__":
    main()
