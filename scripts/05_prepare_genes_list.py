#!/usr/bin/env python3
"""从 skill 内置的 assets/ref/{species} 提取 scautoqc finalize_qc 所需的基因列表。

产出三个文件（写入 --out-dir，默认当前目录；06 流水线显式传 nf-scautoqc/genes_list/）：
  JP_cycle_genes.list  细胞周期基因
  ig_genes.list        免疫球蛋白基因
  tcr_genes.list       T 细胞受体基因

物种来源（内置 assets/ref/）：
  mouse  周期←cell_cycle_genes_updated_2019.csv(94)；IG/TCR←mart_export.txt 按 Gene type
  human  周期←JP_cycle_genes.list(559, 直接复制)；IG/TCR←mart_export.txt 按 Gene type
  custom 用户用 --src-cc 与 --src-mart 提供来源

说明：finalize_qc.py 用正则 ^IGH[ADEMG][0-9]*$ 拆 IG 恒定区（igc），该正则为人类格式，
对小鼠恒定区不匹配——仅使 var['ig'] 含恒定区基因，不影响核心 QC 过滤。
"""
import argparse
import csv
import os
import shutil


def from_mart(mart_path, prefix):
    """从 biomaRt mart_export.txt 按 Gene type 提取 IG_/TR_ 非 pseudogene 的 Gene name"""
    out = set()
    with open(mart_path) as f:
        r = csv.reader(f); next(r)
        for row in r:
            if len(row) < 5:
                continue
            gname, gtype = row[2], row[3]
            if gtype.startswith(prefix) and "pseudo" not in gtype.lower():
                out.add(gname)
    return sorted(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", required=True, choices=["mouse", "human", "custom"])
    ap.add_argument("--out-dir", default=None, help="输出目录（默认当前目录）")
    ap.add_argument("--src-cc", help="[custom] 细胞周期基因 csv（首行当表头跳过；取第 2 列，单列无序号时取第 1 列）")
    ap.add_argument("--src-mart", help="[custom] biomaRt mart_export.txt（含 Gene type 列）")
    ap.add_argument("--assets-dir", default=None, help="内置 assets/ref 根目录（默认脚本相对定位）")
    args = ap.parse_args()

    if args.assets_dir is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        args.assets_dir = os.path.join(os.path.dirname(script_dir), "assets", "ref")
    if args.out_dir is None:
        args.out_dir = os.getcwd()

    if args.species == "custom":
        cc_csv, mart = args.src_cc, args.src_mart
        if not cc_csv or not mart:
            ap.error("custom 模式需 --src-cc 与 --src-mart")
    else:
        sp_dir = os.path.join(args.assets_dir, args.species)
        cc_csv = os.path.join(sp_dir, "cell_cycle_genes_updated_2019.csv") if args.species == "mouse" else None
        mart = os.path.join(sp_dir, "mart_export.txt")

    os.makedirs(args.out_dir, exist_ok=True)

    # 周期基因
    cc_out = os.path.join(args.out_dir, "JP_cycle_genes.list")
    if args.species == "human":
        shutil.copy(os.path.join(args.assets_dir, "human", "JP_cycle_genes.list"), cc_out)
    else:
        cc = []
        with open(cc_csv) as f:
            r = csv.reader(f)
            if next(r, None) is None:
                pass
            for row in r:
                if len(row) >= 2 and row[1].strip():
                    cc.append(row[1].strip())
                elif len(row) == 1 and row[0].strip() and not row[0].strip().isdigit():
                    cc.append(row[0].strip())
        with open(cc_out, "w") as f:
            f.write("\n".join(cc) + "\n")

    # IG / TCR（从 mart 按 Gene type 提取）
    ig = from_mart(mart, "IG_")
    tr = from_mart(mart, "TR_")
    with open(os.path.join(args.out_dir, "ig_genes.list"), "w") as f:
        f.write("\n".join(ig) + "\n")
    with open(os.path.join(args.out_dir, "tcr_genes.list"), "w") as f:
        f.write("\n".join(tr) + "\n")

    print(f"[genes_list] species={args.species} → {args.out_dir}")
    print(f"  JP_cycle_genes.list: {sum(1 for _ in open(cc_out))} 基因")
    print(f"  ig_genes.list:       {len(ig)} 基因")
    print(f"  tcr_genes.list:      {len(tr)} 基因")


if __name__ == "__main__":
    main()
