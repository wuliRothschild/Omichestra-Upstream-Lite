# 流程总览与数据流

> 参考构建默认官方过滤参考（`REF_SOURCE` 四分支：`refdata` 官方包默认 / `download` / `build` / `custom`），白名单与路由见 [reference-setup.md](reference-setup.md)。

## 四阶段数据流

```
              阶段1 比对定量                     阶段2 velocyto
fastq(.gz+md5) ─────────────────▶ outs/ ─────────────────▶ <sid>/velocyto/<sid>.loom
 rsync 到本地      filtered/raw 矩阵    rsync BAM 回本地        (spliced/unspliced/ambiguous)
                  possorted_genome_bam.bam ──归档──▶（中间产物与本地 BAM 用完即删）

              阶段3 cellbender（GPU，串行）
              02b_format_matrices.py ──▶ 标准四件：raw/filtered_feature_bc_matrix{/, .h5}
              （MobiVision 的 *_cell_gene_matrix 复制为标准命名；缺 h5 由 make_10x_h5.py 构建，
                与 DropletUtils::write10xCounts(version="3") 输出逐项等价）
              03b_validate_triple_consistency.py ──▶ 比对矩阵/loom/h5 三位一体硬校验
raw_*.h5 ──────────────────────────────▶ outs/cellbender/cellbender_matrix.h5（去背景矩阵）

              阶段4 质控（scautoqc）
filtered_feature_bc_matrix/ + <sid>.loom + cellbender.h5
   ├─ 04_bridge_to_starsolo.py ──▶ STARsolo 目录（Gene/filtered + Velocyto/spliced|unspliced|ambiguous）
   ├─ gather_matrices.py（starsolo）──▶ gene_velo_cellbender.filtered.h5ad
   ├─ qc.py（original）──▶ <sid>_postqc.h5ad + 阈值 csv + 图
   ├─ 06d_threshold_gate.py ──▶ 每样本阈值闸判定表（pass_threshold_gate）
   ├─ flag_doublet.py ──▶ 每样本 scrublet 分数
   ├─ 08_scdblfinder.R ──▶ scDblFinder 双运行（默认 + 同型）分数
   ├─ 06e_scdblfinder_cluster_qc.py ──▶ 每样本簇级 BH 检验表
   ├─ pool_all.py ──▶ scautoqc_pooled0.h5ad
   ├─ finalize_qc.py ──▶ scautoqc_pooled.h5ad / scautoqc_pooled_filtered.h5ad + sample_passqc_df.csv
   └─ 09_merge_scdblfinder.py ──▶ 合并 scDblFinder 列（双运行 4 列 + 簇级 4 列）
```

## 每阶段契约

| 阶段 | 输入 | 输出 | 关键产物 |
|---|---|---|---|
| 比对 | fastq.gz + md5 | filtered/raw 矩阵 + BAM | `outs/filtered_feature_bc_matrix.h5` |
| velocyto | BAM + filtered barcodes | spliced/unspliced loom | `<sid>/velocyto/<sid>.loom` |
| format | 比对矩阵（含 MobiVision 命名） | 标准四件（mtx 目录 + h5） | `outs/{raw,filtered}_feature_bc_matrix{,.h5}` |
| cellbender | raw_feature_bc_matrix.h5 | 去背景矩阵 | `outs/cellbender/cellbender_matrix.h5` |
| 质控（scautoqc） | filtered mtx + loom + cellbender.h5 | 质控池化对象 | `scautoqc_pooled_filtered.h5ad` |

## 产物目录

```
$OUTPUT_ROOT/$PROJECT_NAME/
├── <sid>/
│   ├── outs/
│   │   ├── filtered_feature_bc_matrix{,.h5}
│   │   ├── raw_feature_bc_matrix{,.h5}
│   │   ├── metrics_summary.csv, web_summary.html
│   │   └── cellbender/cellbender_matrix.h5
│   ├── velocyto/<sid>.loom
│   └── _log                              # 比对成功标记
├── scautoqc/                             # 目录命名对齐 nf-scautoqc 官方
│   ├── 1_gathered_objects/<sid>.gene_velo_cellbender.filtered.h5ad
│   ├── 2_qc_objects/<sid>_postqc.h5ad
│   ├── 2_qc_objects/<sid>_threshold_gate.csv
│   ├── 2_qc_plots_individual/
│   ├── 3_doublet_scores/<sid>_{scrublet,scdblfinder}.csv
│   ├── 3_doublet_scores/<sid>_scdblfinder_cluster_qc.csv
│   ├── 5_qc_plots_overall/
│   ├── qc_thresholds.csv
│   ├── scautoqc_pooled0.h5ad
│   ├── scautoqc_pooled.h5ad
│   ├── scautoqc_pooled_filtered.h5ad     # ★ 下游输入
│   ├── sample_passqc_df.csv
│   ├── scdblfinder_scores.csv            # 合并 scDblFinder 表（默认 + 同型）
│   └── 1_bridge/, 2_qc/                  # 断点续跑所需中间产物
```

## 顺序理由

- **先比对**：产出矩阵与 BAM，是后续所有阶段的基础。
- **紧跟 velocyto**：趁样本还在本地消费 BAM，出 loom 后即可清理 BAM。
- **再做 cellbender**：只需 raw 矩阵，可与 velocyto 并行（GPU vs CPU）。
- **最后质控**：整合 filtered 计数、loom 与 cellbender 矩阵；三者共用同一文件契约。
