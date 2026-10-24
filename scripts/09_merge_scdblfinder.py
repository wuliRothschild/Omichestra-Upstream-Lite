#!/usr/bin/env python3
# ============================================================================
# 合并 per-sample scDblFinder csv + 簇级检验 csv + 阈值闸 mask → 写入
# pooled/filtered obs + 输出标准 csv + 交付层字段重命名
# 用法: 09_merge_scdblfinder.py --sc-root <scautoqc根> --samples <sid1,sid2,...>
# 输入: <sc-root>/3_doublet_scores/<sid>_scdblfinder.csv（08 产出，per-sample；
#       双运行四列；被阈值闸硬过滤的细胞 score=1 / class=singlet）
#       + <sc-root>/3_doublet_scores/<sid>_scdblfinder_cluster_qc.csv（06e 产出：
#       scdblfinder_bh_val/scdblfinder_homotypic_bh_val/
#       scdblfinder_stringent/scdblfinder_homotypic_stringent，sctk 同构簇级检验）
#       + <sc-root>/2_qc_objects/<sid>_threshold_gate.csv（06d 产出，barcode 级
#       pass_threshold_gate mask，Science 2026 口径硬编码阈值）
# 输出: <sc-root>/scdblfinder_scores.csv（合并，标准产出；双运行 4 列 + 簇级检验 4 列）
#       + scautoqc_pooled.h5ad / scautoqc_pooled_filtered.h5ad 的 obs 新增九列
#         score_scdblfinder(float) + doublet_scdblfinder(bool)
#         + score_scdblfinder_homotypic(float) + doublet_scdblfinder_homotypic(bool)
#         + scdblfinder_bh_val(float) + scdblfinder_homotypic_bh_val(float)
#         + scdblfinder_stringent(bool) + scdblfinder_homotypic_stringent(bool)
#         + pass_threshold_gate(bool)，与 scrublet 并存
#       交付层重命名（仅交付对象；per-sample csv 保留工具原生名溯源）：
#         bh_pval → scrublet_bh_val、stringent_doublet → scrublet_stringent
#       pooled_filtered 随后按 pass_threshold_gate 裁剪（硬过滤生效于交付对象）
#       pooled0 不写（finalize 前，无 scrublet，保持阶段一致）
# 依赖: conda scautoqc（anndata）
# ============================================================================
import argparse
import os
import anndata as ad
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description='合并 scDblFinder 得分写入 pooled/filtered obs')
    ap.add_argument('--sc-root', required=True, help='scautoqc 产出根目录（含 3_doublet_scores/ 与根部 pooled h5ad）')
    ap.add_argument('--samples', required=True, help='逗号分隔的 sampleID')
    args = ap.parse_args()

    samples = [s.strip() for s in args.samples.split(',') if s.strip()]
    pooled_dir = args.sc_root  # 官方布局：pooled h5ad 与合并 csv 均在根部

    # ---- 合并 per-sample csv ----
    parts = []
    for sid in samples:
        csv = os.path.join(args.sc_root, '3_doublet_scores', f'{sid}_scdblfinder.csv')
        if os.path.exists(csv):
            parts.append(pd.read_csv(csv))
        else:
            print(f'[warn] 缺失 {csv}（该样本 scDblFinder 未跑或失败）')
    if not parts:
        raise SystemExit('[err] 无 scDblFinder csv 可合并，检查 08_scdblfinder.R 是否跑成功')
    df = pd.concat(parts, ignore_index=True)

    # 构建 barcode-sampleID（对齐 pooled obs_names）+ class 转 bool
    df['barcode_sampleID'] = df['barcode'].astype(str) + '-' + df['sampleID'].astype(str)
    for dbl_col in ['doublet_scdblfinder', 'doublet_scdblfinder_homotypic']:
        df[dbl_col] = df[dbl_col].astype(str) == 'doublet'
    df = df.set_index('barcode_sampleID')

    # ---- 合并 per-sample 簇级检验 csv（06e 产出；缺即硬报错，同阈值闸口径）----
    CLUST_COLS = ['scdblfinder_bh_val', 'scdblfinder_homotypic_bh_val',
                  'scdblfinder_stringent', 'scdblfinder_homotypic_stringent']
    clust_parts = []
    for sid in samples:
        ccsv = os.path.join(args.sc_root, '3_doublet_scores', f'{sid}_scdblfinder_cluster_qc.csv')
        if not os.path.exists(ccsv):
            raise SystemExit(f'[err] 缺簇级检验产物 {ccsv}（06e_scdblfinder_cluster_qc.py 未跑或失败）')
        c = pd.read_csv(ccsv)
        miss = [c0 for c0 in ['barcode'] + CLUST_COLS if c0 not in c.columns]
        if miss:
            raise SystemExit(f'[err] {ccsv} 缺列: {miss}')
        clust_parts.append(c[['barcode'] + CLUST_COLS].assign(sampleID=sid))
    clust_df = pd.concat(clust_parts, ignore_index=True)
    clust_df['barcode_sampleID'] = clust_df['barcode'].astype(str) + '-' + clust_df['sampleID'].astype(str)
    clust_df = clust_df.set_index('barcode_sampleID')[CLUST_COLS]
    for c0 in CLUST_COLS[:2]:
        clust_df[c0] = pd.to_numeric(clust_df[c0], errors='raise').astype('float64')
    for c0 in CLUST_COLS[2:]:
        clust_df[c0] = clust_df[c0].astype(str).str.lower() == 'true'
    df = df.join(clust_df, how='left')
    if df[CLUST_COLS].isna().any().any():
        raise SystemExit('[err] 簇级检验合并存在未命中 barcode（08 与 06e 细胞集合不一致）')

    # ---- 输出合并 csv（标准产出）----
    out_csv = os.path.join(pooled_dir, 'scdblfinder_scores.csv')
    df.reset_index()[['barcode_sampleID', 'score_scdblfinder', 'doublet_scdblfinder',
                      'score_scdblfinder_homotypic', 'doublet_scdblfinder_homotypic'] + CLUST_COLS].to_csv(
        out_csv, index=False)
    print(f'[scdblfinder] 合并 {len(df)} 细胞 -> {out_csv} '
          f'| doublet rate={df["doublet_scdblfinder"].mean():.4f} '
          f'| homotypic rate={df["doublet_scdblfinder_homotypic"].mean():.4f}')

    # ---- 合并 per-sample 阈值闸 mask（06d 产出，2_qc_objects/）----
    gate_parts = []
    for sid in samples:
        gcsv = os.path.join(args.sc_root, '2_qc_objects', f'{sid}_threshold_gate.csv')
        if not os.path.exists(gcsv):
            raise SystemExit(f'[err] 缺阈值闸产物 {gcsv}（06d_threshold_gate.py 未跑或失败）')
        g = pd.read_csv(gcsv)
        gate_parts.append(g[['barcode', 'pass_threshold_gate']].assign(sampleID=sid))
    gate_df = pd.concat(gate_parts, ignore_index=True)
    gate_df['barcode_sampleID'] = gate_df['barcode'].astype(str) + '-' + gate_df['sampleID'].astype(str)
    gate_pass = gate_df.set_index('barcode_sampleID')['pass_threshold_gate'].astype(bool)

    # ---- 写入 pooled/filtered obs（pooled0 不写：finalize 前无 scrublet，保持一致）----
    pairs = [('score_scdblfinder', 'doublet_scdblfinder'),
             ('score_scdblfinder_homotypic', 'doublet_scdblfinder_homotypic')]
    for f in ['scautoqc_pooled.h5ad', 'scautoqc_pooled_filtered.h5ad']:
        path = os.path.join(pooled_dir, f)
        a = ad.read_h5ad(path)
        for score_col, dbl_col in pairs:
            score = df[score_col].reindex(a.obs_names)
            n_na = int(score.isna().sum())
            if n_na:
                missing = a.obs_names[score.isna().to_numpy()][:5].tolist()
                raise SystemExit(
                    f'[err] {f}: {n_na}/{a.n_obs} 细胞无 {score_col} 记录（前5: {missing}）；'
                    '未命中细胞无合法值可填（填 False/0 属编造）——检查该样本 08_scdblfinder.R 是否成功、'
                    'barcode-sampleID 格式是否一致')
            a.obs[score_col] = score.astype('float64').values
            a.obs[dbl_col] = df[dbl_col].reindex(a.obs_names).astype(bool).values
        gate = gate_pass.reindex(a.obs_names)
        n_na = int(gate.isna().sum())
        if n_na:
            missing = a.obs_names[gate.isna().to_numpy()][:5].tolist()
            raise SystemExit(
                f'[err] {f}: {n_na}/{a.n_obs} 细胞无阈值闸记录（前5: {missing}）——'
                '检查 06d_threshold_gate.py 是否全样本跑过、barcode-sampleID 格式是否一致')
        a.obs['pass_threshold_gate'] = gate.astype(bool).values
        # 簇级检验 4 列（06e；NA 硬校验同上）
        for c0 in CLUST_COLS:
            v = clust_df[c0].reindex(a.obs_names)
            n_na = int(v.isna().sum())
            if n_na:
                missing = a.obs_names[v.isna().to_numpy()][:5].tolist()
                raise SystemExit(
                    f'[err] {f}: {n_na}/{a.n_obs} 细胞无 {c0} 记录（前5: {missing}）——'
                    '检查 06e_scdblfinder_cluster_qc.py 是否全样本跑过')
            a.obs[c0] = (v.astype('float64').values if c0.endswith('_bh_val')
                         else v.astype(bool).values)
        # 交付层重命名（幂等：旧名在→改；仅新名在→已改过重跑跳过；皆无→finalize 未写入，报错）
        for old, new in [('bh_pval', 'scrublet_bh_val'),
                         ('stringent_doublet', 'scrublet_stringent')]:
            if old in a.obs.columns:
                if new in a.obs.columns:
                    raise SystemExit(f'[err] {f}: obs 同时存在 {old} 与 {new}，更名冲突')
                a.obs.rename(columns={old: new}, inplace=True)
                print(f'[rename] {f}: {old} -> {new}')
            elif new not in a.obs.columns:
                raise SystemExit(f'[err] {f}: obs 既无 {old} 也无 {new}（finalize_qc 未写入 scrublet 判定列）')
        # 阈值闸硬过滤生效于交付对象：filtered 按 mask 裁剪（pooled 保留全量+标注）
        if f == 'scautoqc_pooled_filtered.h5ad':
            n0 = a.n_obs
            a = a[a.obs['pass_threshold_gate'].to_numpy()].copy()
            print(f'[threshold_gate] {f}: {n0} -> {a.n_obs}（阈值闸剔除 {n0 - a.n_obs}）')
        a.write_h5ad(path, compression="gzip")
        print(f'[scdblfinder] {f}: n={a.n_obs} | score NA=0 '
              f'| doublet rate={a.obs["doublet_scdblfinder"].mean():.4f} '
              f'| homotypic rate={a.obs["doublet_scdblfinder_homotypic"].mean():.4f} '
              f'| stringent={int(a.obs["scdblfinder_stringent"].sum())} '
              f'| homotypic_stringent={int(a.obs["scdblfinder_homotypic_stringent"].sum())}')


if __name__ == '__main__':
    main()
