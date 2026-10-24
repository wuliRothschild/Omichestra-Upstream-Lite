#!/usr/bin/env python3
# ============================================================================
# upsert ## upstream 阶段 段落到项目 HANDOFF.md（交接给下游）
# 用法: 07_write_handoff.py --project <P> --species <S> --sc-root <dir> --samples <sid1,sid2,...>
#                            [--handoff-md Document/HANDOFF.md] [--tech sc] [--chemistry 10x_3prime_v3]
#                            [--pooled0 <path> --pooled <path> --filtered <path>]
# 产出: 段级 upsert 到 --handoff-md（默认 Document/HANDOFF.md），段落标题 ## upstream 阶段
#       （段落内子标题用 ### 避免与段落边界 ## 冲突）
# 依赖: conda scautoqc（anndata，仅用于读 cell counts）
# ============================================================================
import argparse
import os
import re
import datetime
import anndata as ad


def ncells(path):
    if path and os.path.exists(path):
        try:
            return ad.read_h5ad(path, backed='r').shape[0]
        except Exception as e:
            print(f'[warn] 读 {path} 失败: {e}')
            return '?'
    return '?'


def upsert_handoff_section(md_path, title, body):
    """段级 upsert：## title 段落存在则替换（到下一个 ## 或文件末），不存在则追加；文件不存在则创建。"""
    header = f"## {title}"
    section = header + "\n" + body.rstrip()
    if os.path.exists(md_path):
        with open(md_path, encoding='utf-8') as f:
            content = f.read()
        pat = re.compile(r'^## ' + re.escape(title) + r'\n.*?(?=\n^## |\Z)', re.MULTILINE | re.DOTALL)
        if pat.search(content):
            new_content = pat.sub(lambda m: section, content)
        else:
            new_content = content.rstrip() + "\n\n" + section
    else:
        new_content = "# HANDOFF — 项目全流程记录\n\n" + section
    os.makedirs(os.path.dirname(os.path.abspath(md_path)), exist_ok=True)
    with open(md_path, 'w', encoding='utf-8') as f:
        f.write(new_content.rstrip() + "\n")


def main():
    ap = argparse.ArgumentParser(description='upsert ## upstream 阶段 到项目 HANDOFF.md')
    ap.add_argument('--project', required=True)
    ap.add_argument('--species', required=True)
    ap.add_argument('--tech', default='sc')
    ap.add_argument('--chemistry', default='10x_3prime_v3')
    ap.add_argument('--sc-root', required=True, help='scautoqc 产出根目录（官方布局，pooled h5ad 在根部）')
    ap.add_argument('--samples', required=True, help='逗号分隔的 sampleID')
    ap.add_argument('--handoff-md', default='Document/HANDOFF.md', help='项目 HANDOFF.md 路径（默认 Document/HANDOFF.md）')
    ap.add_argument('--pooled0', default='')
    ap.add_argument('--pooled', default='')
    ap.add_argument('--filtered', default='')
    ap.add_argument('--section', default='', help='段落标题后缀（同一 HANDOFF 承载多数据集时填数据集名；非空则标题为 ## upstream 阶段（<值>），缺省为契约精确标题 ## upstream 阶段）')
    ap.add_argument('--gate-min-counts', default='200')
    ap.add_argument('--gate-min-genes', default='200')
    ap.add_argument('--gate-max-mito-cell', default='30')
    ap.add_argument('--gate-max-mito-nuclei', default='5')
    args = ap.parse_args()

    samples = [s.strip() for s in args.samples.split(',') if s.strip()]
    pooled0 = args.pooled0 or os.path.join(args.sc_root, 'scautoqc_pooled0.h5ad')
    pooled = args.pooled or os.path.join(args.sc_root, 'scautoqc_pooled.h5ad')
    filtered = args.filtered or os.path.join(args.sc_root, 'scautoqc_pooled_filtered.h5ad')

    n0, n1, nf = ncells(pooled0), ncells(pooled), ncells(filtered)
    # 阈值闸（06d）通过细胞数：per-sample csv 汇总（filtered 细胞数已含闸裁剪，此处单列闸层）
    n_gate = '?'
    try:
        import pandas as pd
        n_gate = 0
        for s in samples:
            g = pd.read_csv(os.path.join(args.sc_root, '2_qc_objects', f'{s}_threshold_gate.csv'))
            n_gate += int(g['pass_threshold_gate'].sum())
    except Exception as e:
        print(f'[warn] 读阈值闸 csv 失败: {e}')
        n_gate = '?'
    sample_rows = '\n'.join(f"| {s} |" for s in samples)

    body = f"""生成时间: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | skill: Omichestra-upstream
项目: {args.project}

### 技术信息
- species: {args.species}
- tech: {args.tech}
- chemistry: {args.chemistry}
- barcode_format: barcode-sampleID

### 样本
| sample_id |
|---|
{sample_rows}

### scautoqc 产出
根目录: {args.sc_root}
| 文件 | 细胞数 | .X | layers | QC列 |
|---|---|---|---|---|
| scautoqc_pooled0.h5ad | {n0} | cellbender | raw/spliced/unspliced/ambiguous | - |
| scautoqc_pooled.h5ad | {n1} | cellbender | (已删) | scrublet系7列(scrublet_score/scrublet_score_z/cluster_scrublet_score/scrublet_bh_val/scrublet_done/doublet/scrublet_stringent)/scDblFinder系8列(score_scdblfinder/doublet_scdblfinder/score_scdblfinder_homotypic/doublet_scdblfinder_homotypic/scdblfinder_bh_val/scdblfinder_homotypic_bh_val/scdblfinder_stringent/scdblfinder_homotypic_stringent)/pass_threshold_gate/good_qc_cluster/pass_auto_filter |
| scautoqc_pooled_filtered.h5ad | {nf} | cellbender | (已删) | 同pooled(子集；簇级 QC + 阈值闸双过滤) |

### 关键说明
- pooled0 与 pooled 的 obs_names 顺序完全一致（可按位置合并 layers+QC列）
- .X = cellbender 去背景整数 counts（未 normalize/log）
- barcode 格式 = barcode-sampleID（如 AAACCCAAGTGTGTTC-<sid>）
- scdblfinder_scores.csv：scDblFinder 合并表（双运行 4 列 + 簇级检验 4 列），已写入 pooled/filtered obs；pooled0 不含（finalize 前阶段）
- 簇级 BH 检验（06e，sctk 同构）：scDblFinder 双运行各配 `*_bh_val`（簇得分中位数→z→BH-FDR，簇中位数仅取阈值闸通过细胞）与 `*_stringent`（score>0.5 或 bh_val<0.05）；被闸细胞 bh_val=1、stringent=TRUE（score=1 回填所致）；聚类在 postqc 全量上按 sctk 配方重算（leiden 分辨率1，确定性随机种子），失败兜底 bh_val=1 仅警告不中断
- 交付层字段命名消歧：交付对象中 scrublet 的 bh_pval/stringent_doublet 已更名 scrublet_bh_val/scrublet_stringent；per-sample scrublet csv 保留工具原生列名（溯源）
- 阈值闸（06d，Science 2026 口径硬过滤）：n_counts>={args.gate_min_counts}、n_genes>={args.gate_min_genes}、percent_mito<{args.gate_max_mito_cell}（sc）/<{args.gate_max_mito_nuclei}（sn）；三层细胞数 = pooled 全量 {n1} → 闸通过 {n_gate} → filtered 交付 {nf}（簇级 QC ∩ 阈值闸）；被闸细胞 scDblFinder score=1 回填，pooled obs 含 pass_threshold_gate 列
"""
    sec = args.section.strip()
    title = f'upstream 阶段（{sec}）' if sec else 'upstream 阶段'
    upsert_handoff_section(args.handoff_md, title, body)
    print(f'[handoff] upsert ## {title} → {args.handoff_md} | samples={len(samples)} | cells pooled0/pooled/filtered = {n0}/{n1}/{nf}')


if __name__ == '__main__':
    main()
