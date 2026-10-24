#!/usr/bin/env python3
"""写每样本 outs/metadata.json 的十类基础信息（交付契约第 5 项）。

十字段：group/age/batch/datasets/sample/tech/gender/tissue/chemistry/species。
来源约定（不可从下游 sample_metadata.csv 取——upstream 运行时它尚不存在）：
- 运行 upstream 前由用户主动交代，或由 agent 阻塞式向用户询问后填写；
  用户明确不提供的一律取值 "Unknown"（不猜测、不编造）
- per-sample 字段（group/age/batch/gender/sample）经 --meta-csv 提供
  （列：sample_id + group/age/gender/batch/sample 任意子集，表头必需；缺文件/缺列/缺行/空值 → Unknown）
- gender 取值规范 F/M：常见写法（男/女/Male/Female，大小写不敏感）自动归一；
  其余值原样保留（不猜测），未知仍落 Unknown
- 数据集级字段经 CLI 提供：--datasets（默认=cr-root 末级目录名，即 PROJECT_NAME）、
  --tech/--chemistry/--tissue/--species（缺省 → Unknown；物种通常直接传 run_config SPECIES，
  不从样本名推断；custom 物种走scautoqc 阶段 scautoqc 时按运行级值原样落）
- sample 默认=样本名（sample_id）；meta-csv 提供 sample 列时按列值（如业务样本名）

幂等：已存在的 metadata.json 默认跳过（防编排脚本的兜底调用覆盖已收集的真实值）；
--force 时覆盖十个基础键并保留 stages 等其余键（stages 由 06c_record_stage_genes.py
在质控分支完成后写入）。

用法: 01c_write_metadata.py --cr-root <根> --samples sid1,sid2,...
       [--meta-csv <csv>] [--datasets <名>] [--tech <值>] [--chemistry <值>] [--tissue <值>]
       [--species <human|mouse|Unknown>] [--force]
"""
import argparse
import csv
import json
import os

BASE_KEYS = ("group", "age", "batch", "datasets", "sample", "tech", "gender", "tissue", "chemistry", "species")
PER_SAMPLE_KEYS = ("group", "age", "gender", "batch", "sample")
UNKNOWN = "Unknown"

# gender 规范 F/M 的同义写法（大小写不敏感）；不在表内的值原样保留，不猜测
GENDER_MAP = {
    "f": "F", "female": "F", "女": "F",
    "m": "M", "male": "M", "男": "M",
}


def load_meta_csv(path):
    """返回 {sample_id: {key: value}}；文件不存在/无表头 → 空表（调用方已允许缺省）"""
    if not path or not os.path.isfile(path):
        return {}
    with open(path, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    out = {}
    for r in rows:
        sid = (r.get("sample_id") or "").strip()
        if not sid:
            continue
        rec = {}
        for k in PER_SAMPLE_KEYS:
            v = (r.get(k) or "").strip()
            if k == "gender" and v:
                nv = GENDER_MAP.get(v.lower(), v)
                if nv != v:
                    print(f"[metadata] gender 归一化: {v} -> {nv}（{sid}）")
                v = nv
            rec[k] = v
        out[sid] = rec
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cr-root", required=True, help="样本根目录（含 <sid>/outs/）")
    ap.add_argument("--samples", required=True, help="逗号分隔 sampleID")
    ap.add_argument("--meta-csv", default="", help="per-sample 元数据 csv（可选；列 sample_id,group,age,gender,batch[,sample]）")
    ap.add_argument("--datasets", default="", help="数据集名（默认=cr-root 末级目录名，即 PROJECT_NAME）")
    ap.add_argument("--tech", default=UNKNOWN)
    ap.add_argument("--chemistry", default=UNKNOWN)
    ap.add_argument("--tissue", default=UNKNOWN)
    ap.add_argument("--species", default=UNKNOWN, help="物种（human|mouse|custom|Unknown；通常传 run_config SPECIES，不从样本名推断）")
    ap.add_argument("--force", action="store_true", help="覆盖已存在的 metadata.json 基础键（保留 stages 等其余键）")
    args = ap.parse_args()

    per_sample = load_meta_csv(args.meta_csv)
    datasets = args.datasets.strip() or os.path.basename(os.path.normpath(args.cr_root))

    for sid in [s.strip() for s in args.samples.split(",") if s.strip()]:
        outs = os.path.join(args.cr_root, sid, "outs")
        if not os.path.isdir(outs):
            raise SystemExit(f"[err] {sid} 无 outs 目录: {outs}（01 比对完成后才能写 metadata.json）")
        row = per_sample.get(sid, {})
        meta = {
            "group": row.get("group") or UNKNOWN,
            "age": row.get("age") or UNKNOWN,
            "batch": row.get("batch") or UNKNOWN,
            "datasets": datasets,
            "sample": row.get("sample") or sid,
            "tech": args.tech.strip() or UNKNOWN,
            "gender": row.get("gender") or UNKNOWN,
            "tissue": args.tissue.strip() or UNKNOWN,
            "chemistry": args.chemistry.strip() or UNKNOWN,
            "species": args.species.strip() or UNKNOWN,
        }
        path = os.path.join(outs, "metadata.json")
        if os.path.isfile(path):
            if not args.force:
                print(f"[metadata] {sid} 已存在，跳过（--force 覆盖）")
                continue
            with open(path, encoding="utf-8") as f:
                existing = json.load(f)
            existing.update({k: meta[k] for k in BASE_KEYS})
            meta = existing
        with open(path, "w", encoding="utf-8") as f:
            json.dump(meta, f, ensure_ascii=False, indent=2)
            f.write("\n")
        print(f"[metadata] {sid} -> {path}")


if __name__ == "__main__":
    main()
