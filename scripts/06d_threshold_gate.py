#!/usr/bin/env python3
# ============================================================================
# 06d_threshold_gate.py — scDblFinder 前阈值硬过滤（per-sample，Science 2026 口径）
# 位点: qc.py（scAutoQC 细胞级标注）之后、08_scdblfinder.R 之前（06 内步骤 4a）。
# 机制: 对 postqc obs 已有的 n_counts / n_genes / percent_mito 做纯阈值判定
#       （三项全过为通过），不做 GMM——GMM+initial bounds 是 scAutoQC cellwise_qc
#       （qc.py）的既有职责，不在此重复。scrublet 等双胞指标不参与过滤。
# 契约默认阈值（Science 2026, doi:10.1126/science.aeb2672，不含其 scrublet 项）:
#   n_counts>=200、n_genes>=200、percent_mito<30（sc）/ <5（sn）
#   run_config GATE_* 变量可覆盖（06 传参），实取值入指纹 params（指纹判据）。
# 用法: python 06d_threshold_gate.py --postqc <postqc.h5ad> --sid <sampleID>
#        --out <2_qc_objects/<sid>_threshold_gate.csv> --cell-or-nuclei <cell|nuclei>
#        [--gate-min-counts 200] [--gate-min-genes 200]
#        [--gate-max-mito-cell 30] [--gate-max-mito-nuclei 5]
# 输出: csv（barcode, pass_n_counts, pass_n_genes, pass_percent_mito,
#       pass_threshold_gate）——供 08_scdblfinder.R（--keep-csv 消费，未通过细胞
#       回填 score=1）与 09_merge_scdblfinder.py（合并入 obs + 裁剪 filtered）；
#       兼作断点续跑的 per-sample 完成标记之一
# 依赖: conda scautoqc（anndata）
# ============================================================================
import argparse
import anndata as ad
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description='scDblFinder 前阈值硬过滤（Science 2026 口径，纯阈值无 GMM）')
    ap.add_argument('--postqc', required=True, help='2_qc_objects/<sid>_postqc.h5ad')
    ap.add_argument('--sid', required=True, help='sampleID（日志用）')
    ap.add_argument('--out', required=True, help='输出 csv 路径')
    ap.add_argument('--cell-or-nuclei', default='cell', choices=['cell', 'nuclei'],
                    help='测序类型：mito 阈值分流（cell<30 / nuclei<5）')
    ap.add_argument('--gate-min-counts', type=float, default=200)
    ap.add_argument('--gate-min-genes', type=float, default=200)
    ap.add_argument('--gate-max-mito-cell', type=float, default=30)
    ap.add_argument('--gate-max-mito-nuclei', type=float, default=5)
    args = ap.parse_args()

    max_mito = args.gate_max_mito_nuclei if args.cell_or_nuclei == 'nuclei' else args.gate_max_mito_cell

    # backed 模式只读 obs（X/layers 留盘）；NaN 比较恒 False = fail-closed（保守判负）
    a = ad.read_h5ad(args.postqc, backed='r')
    obs = a.obs
    for col in ('n_counts', 'n_genes', 'percent_mito'):
        if col not in obs.columns:
            a.file.close()
            raise SystemExit(f'[err] {args.sid}: postqc obs 缺 {col} 列（qc.py 应先算好 QC 指标）')

    gate = pd.DataFrame({
        'barcode': pd.Index(a.obs_names.astype(str)),
        'pass_n_counts': (obs['n_counts'] >= args.gate_min_counts).to_numpy(),
        'pass_n_genes': (obs['n_genes'] >= args.gate_min_genes).to_numpy(),
        'pass_percent_mito': (obs['percent_mito'] < max_mito).to_numpy(),
    })
    n_obs = len(gate)
    a.file.close()
    gate['pass_threshold_gate'] = (
        gate['pass_n_counts'] & gate['pass_n_genes'] & gate['pass_percent_mito'])
    gate.to_csv(args.out, index=False)

    n_pass = int(gate['pass_threshold_gate'].sum())
    print(f'[threshold_gate] {args.sid}: n={n_obs} | '
          f'counts>={args.gate_min_counts:g} 剔除 {(~gate["pass_n_counts"]).sum()} | '
          f'genes>={args.gate_min_genes:g} 剔除 {(~gate["pass_n_genes"]).sum()} | '
          f'mito<{max_mito:g}({args.cell_or_nuclei}) 剔除 {(~gate["pass_percent_mito"]).sum()} | '
          f'通过 {n_pass}/{n_obs}（{n_pass / max(n_obs, 1):.1%}）-> {args.out}')
    if n_pass < 1:
        raise SystemExit(f'[err] {args.sid}: 阈值闸后无细胞剩余，样本质量过低或阈值不适配'
                         '（run_config GATE_* 可校准）')


if __name__ == '__main__':
    main()
