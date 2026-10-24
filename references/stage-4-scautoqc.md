# scAutoQC 算法手册（nf-scautoqc / sctk）

> 本文档是 Omichestra-upstream 阶段4 默认分支（`06_scautoqc_pipeline.sh`）的完整算法手册：scAutoQC 算法机制（四步工作流、GMM 自适应阈值、QC 簇投票、mito 三档循环）、官方 nf-scautoqc 能力全景、本 skill 的固定执行层与参数全表、产出契约与故障排除。内容来源：sctk 官方教程页（teichlab.github.io/sctk/notebooks/automatic_qc.html）逐字引用 + nf-scautoqc 仓库（github.com/cellgeni/nf-scautoqc：README/CHANGELOG/main.nf/nextflow.config/bin 脚本）逐行 + 本环境 sctk 安装包源码实证 + Oliver et al. 2024（Nature，doi:10.1038/s41586-024-07571-1）方法段。
> 适用前提：`run_config.sh` 已完成配置（`ENV_SCAUTOQC`/`NF_SCAUTOQC_DIR`/`SCANPY_PYTHON`/`RSCRIPT` 等，见 [environment-setup.md](environment-setup.md)）。

## 目录

1. [算法总览与两步判定哲学](#1-算法总览与两步判定哲学)
2. [八项 QC 指标与物种基因标记（calculate_qc）](#2-八项-qc-指标与物种基因标记calculate_qc)
3. [GMM 自适应阈值机制（fit_gaussian / cellwise_qc）](#3-gmm-自适应阈值机制fit_gaussian--cellwise_qc)
4. [QC 嵌入与簇级投票（generate_qc_clusters / clusterwise_qc）](#4-qc-嵌入与簇级投票generate_qc_clusters--clusterwise_qc)
5. [original 模式：mito 三档循环与编码汇总](#5-original-模式mito-三档循环与编码汇总)
6. [官方能力全景：三种 qc_mode 与六种 run mode](#6-官方能力全景三种-qc_mode-与六种-run-mode)
7. [双胞检测双线：scrublet 与 scDblFinder](#7-双胞检测双线scrublet-与-scdblfinder)
8. [池化与终判（pool_all / finalize_qc）](#8-池化与终判pool_all--finalize_qc)
9. [本 skill 执行层（06_scautoqc_pipeline.sh）](#9-本-skill-执行层06_scautoqc_pipelinesh)
10. [输入与产出契约](#10-输入与产出契约)
11. [与 Oliver et al. 2024 方法的对应关系](#11-与-oliver-et-al-2024-方法的对应关系)
12. [故障排除](#12-故障排除)
13. [下游衔接](#13-下游衔接)

---

## 1. 算法总览与两步判定哲学

scAutoQC 是 Teich 实验室 sctk 包（Single Cell analysis Tool Kit）的自动 QC 工作流；nf-scautoqc 是它的 Nextflow 实现（仓库 README 原文："nf-scautoqc is the Nextflow implementation of scAutoQC pipeline used in Oliver et al, 2024"），即 Oliver et al. 2024 Pan-GI atlas 的质控流程。

官方教程给出的四个核心步骤（逐字）：

- `sctk.calculate_qc()` computes QC measures
- `sctk.cellwise_qc()` performs cell-level QC calling
- `sctk.generate_qc_clusters()` creates QC-space clusters of cells
- `sctk.clusterwise_qc()` computes cluster-wide statistics of cell quality

并明确建议（逐字）："It is recommended to run scAutoQC on a per-sample basis."

**两步判定哲学**（Oliver 2024 Methods 原文）："Classification of cells/droplets as passing or failing QC was then performed in a two-step process, first by classifying each cell as passing or failing QC based on four-metric parameters and thresholds set by a Gaussian mixture model (GMM). ... Then, whole clusters were classified as passing QC if 50% or more of individual cells within the cluster passed QC."——即先用 GMM 在**细胞级**逐指标定界判定，再以 **QC 簇多数表决**吸收单细胞级的误判波动。簇级投票是强平滑器：零星的单细胞判定翻转改变不了簇的多数派归属（这也是论文强调该流程"对起始参数不敏感"的机制基础）。

**本 skill 的取舍**：不经 Nextflow 壳（官方 nextflow.config 绑定 Sanger LSF/singularity 绝对路径），由 `06_scautoqc_pipeline.sh` 直接串联官方 `bin/` 脚本等价执行，产出目录命名完全对齐官方；所有参数固定（§9 全表），不经命令行临时改参。

## 2. 八项 QC 指标与物种基因标记（calculate_qc）

### 2.1 指标计算机制（qc.py `calculate_qc`）

对 gather 产出的 h5ad（`.X`=cellbender 去背景计数，layers=raw/spliced/unspliced/ambiguous）：

1. `sk.calculate_qc(ad, log1p=False)`：内部调 `scanpy.pp.calculate_qc_metrics(qc_vars=[mito, ribo, hb], percent_top=[50])`，产出 n_counts/log1p_n_counts/n_genes/log1p_n_genes/percent_mito/percent_ribo/percent_hb/percent_top50（及对应 n_counts_*）。
2. 若有 `raw` 层：对 raw 层再算一遍（suffix=`_raw`），并派生 `percent_soup = (1 − n_counts/n_counts_raw) × 100`——cellbender 从 droplet 计数中扣除背景 RNA 的比例（观测派生量，0–100）。
3. 若有 `spliced`+`unspliced` 层：各算一遍（suffix=`_spliced`/`_unspliced`），派生 `percent_spliced = n_counts_spliced / (n_counts_spliced + n_counts_unspliced) × 100`（两层皆 0 的细胞 `fillna(0)`）。**缺剪接层时 percent_spliced 自动从指标列表剔除**（`qc_metrics.pop()`，官方行为），不报错。

### 2.2 八项指标定义（`DEFAULT_QC_METRICS`，顺序固定）

| 指标 | 定义 | 依赖 |
|---|---|---|
| `log1p_n_counts` | log1p（细胞总计数），`.X`（cellbender 去背景后） | 基础 |
| `log1p_n_genes` | log1p（检出基因数） | 基础 |
| `percent_mito` | 线粒体基因计数占比（%） | var['mito'] 标记 |
| `percent_ribo` | 核糖体基因计数占比（%） | var['ribo'] 标记 |
| `percent_hb` | 血红蛋白基因计数占比（%） | var['hb'] 标记 |
| `percent_top50` | 前 50 高表达基因计数占比（scanpy `pct_counts_in_top_50_genes`） | 基础 |
| `percent_soup` | (1 − n_counts/n_counts_raw)×100，cellbender 扣除比例 | raw 层 |
| `percent_spliced` | 剪接计数占比（%） | spliced+unspliced 层 |

**观测与过滤的区分**：8 项全部写入 obs 并作为特征参与 QC 嵌入聚类（§4）；但只有 4 个核心指标（n_counts/n_genes/percent_mito/percent_spliced）进入 GMM 过滤（§5），其余 4 项只观测不过滤。

### 2.3 基因标记正则（物种适配，skill 修补点）

sctk `calculate_qc` 的默认 flags 是**人类全大写**正则（源码逐字）：`{"mito": r"^MT-", "ribo": r"^RP[LS]", "hb": r"^HB"}`——在小鼠（`mt-`/`Rpl`/`Hbb`）上近乎全失配。skill 的 qc.py 改用物种无关正则并显式传入（`sk.calculate_qc(..., flags=SPECIES_AGNOSTIC_FLAGS)`）：

| 标记 | skill 正则 | 说明 |
|---|---|---|
| mito | `^[mM][tT]-` | mt-/MT-/Mt-/mT- 任意大小写前缀 |
| ribo | `^[rR][pP][lLsS]` | Rpl/Rps/RPL/RPS/... |
| hb | `^[Hh][Bb](?![pP])` | Hb/HB/hb，lookahead 排除假基因 Hbp/HBP |

对人源数据两版等价（零实质影响）；对鼠源数据 skill 版是必要修正。

## 3. GMM 自适应阈值机制（fit_gaussian / cellwise_qc）

### 3.1 机制（官方教程逐字）

"Under the hood, `sctk.cellwise_qc()` models each of the specified QC metrics as a Gaussian mixture model, with the model proposing possible low and high value cutoffs at points where the probability density function of the mixture falls below a threshold (0.05 by default)."

即：每个过滤指标的一维分布用高斯混合模型（GMM）拟合，取混合概率密度函数（PDF）从峰值向两侧跌破 0.05 的位置作为该指标的数据驱动上下界。

### 3.2 指标参数五列（官方教程逐字解释）

- `min`：the minimum value to use in modelling the distribution; values lower than this will be skipped; set to np.nan to use the encountered minimum
- `max`：the maximum value to use in modelling the distribution; values higher than this will be skipped; set to np.nan to use the encountered maximum
- `scale`："log" to log-transform the QC metric values prior to modelling, "linear" to keep the values un-transformed
- `side`："min_only" to only perform filtering of low values, "max_only" to only perform filtering of high values, "both" to perform filtering from both sides
- `min_pass_rate`：the minimum proportion of cells that need to pass the filtering on this feature; if the fraction is too low, the GMM is discarded and the QC values are thresholded based on the provided minimum/maximum

**要点**：min/max 不是硬过滤线，而是 **GMM 拟合域**（先把指标截断到 [min, max] 再拟合）；`scale="log"` 时对值与 min/max 同步 `log1p` 变换（源码实证）；side 决定最终保留 GMM 界的那一侧（`min_only` 则上界放开为数据最大值，`max_only` 反之，`both` 双侧都用 GMM 界）。

### 3.3 fit_gaussian 源码行为链（sctk `_pipeline.py` 实证）

1. 截断：`x_fit = x[(x >= xmin) & (x <= xmax)]`；
2. 数值缩放：`_scale_factor(x_fit)` 把拟合域量级归一（防极端值数值问题）；
3. 逐组件数拟合：`for n in np.arange(10)+1`（1–10 个高斯分量），每个 `GaussianMixture(n_components=n, random_state=0)`，未收敛则 `warm_start` 续拟合至收敛；
4. **BIC 选模型**：`bics.append(gmm.bic(x_fit))`，取 BIC 最小者（`np.argmin`）——源码 docstring 逐字："The best GMM will be selected based on BIC"，注释逐字："# pick best one based on BIC (the lower the better)"；
5. 定界：在 500 点网格上算混合 PDF，从 PDF 最大值位置向两侧找跌破 `threshold=0.05` 的点；`cutoff="inner"` 取离峰值**最近**的跌破点（默认），`"outer"` 取**最远**的（docstring 逐字："the closest ("inner") or furthest ("outer") point relative to the GMM PDF maximum where the PDF falls below `threshold`"）。

**BIC 说明**：BIC（Bayesian Information Criterion，贝叶斯信息准则）是模型选择的定量指标，`BIC = k·ln(n) − 2·ln(L̂)`（k=自由参数个数，n=样本数，L̂=最大似然），拟合优度与复杂度权衡，**越小越好**。分量越多似然越高（拉低 BIC）但参数惩罚越重（抬高 BIC），最低点即"够用"的组件数。出处诚实标注：BIC 见于 sctk 源码（docstring/注释/代码）与 Oliver 2024 Methods（"scAutoQC was subsequently improved to automate the best model fit between 1 and 10 components based on the Bayesian information criterion"）；**sctk 官方教程页与 nf-scautoqc 仓库文档均未提及 BIC**——文档滞后于代码，软件行为以源码为准。论文 atlas 实际跑法是固定 10 组件（"the number of GMM components was set to 10 for an overfit model"），现行 sctk 的 BIC 1–10 自动选择是论文自认的后续改进版（两法差异的实测影响见 §11）。

### 3.4 cellwise_qc 的兜底链（源码实证，按序生效）

1. `fit_gaussian` 抛 `ValueError`（拟合失败）→ 回退 xmin/xmax（未给则数据极值）；
2. GMM 界判定后通过细胞数 `< n_obs × min_pass_rate` → **丢弃 GMM**，按 side 回退硬界 min/max（教程逐字："the GMM is discarded and the QC values are thresholded based on the provided minimum/maximum"）；
3. 各指标通过布尔向量按 **AND** 合成 `cell_qc_key`（四指标全过才算细胞过）；
4. 若最终 0 细胞通过 → 打印 "No cells passed. Performing simple filtering on counts, genes and mito%"，回退简单过滤 `n_counts≥min & n_genes≥min & percent_mito<max`；
5. 每指标的实际界与通过数实时打印（`{metric}: [{low}, {high}], {n_pass}/{n_obs} passed`），并写入 `uns['scautoqc_ranges']`（qc.py 据此产出 per-sample `{sid}_qc_thresholds.csv` tidy 长表：sample/metric/low/high/n_pass/n_total/pass_rate + `all_metrics` 汇总行）。

## 4. QC 嵌入与簇级投票（generate_qc_clusters / clusterwise_qc）

### 4.1 QC 簇构建（sctk 源码参数链，全部固定）

以 obs 中 8 项指标为特征建辅助对象 `aux_ad`，跑一条迷你 scanpy 流程：

| 步骤 | 调用 | 固定参数 |
|---|---|---|
| 标准化 | `sc.pp.scale(aux_ad)` | 默认 |
| PCA | `sc.pp.pca(aux_ad, n_comps=n_pcs)` | `n_pcs = max(2, 指标数−2)`（8 指标→6） |
| 近邻图 | `sc.pp.neighbors(aux_ad, n_neighbors=n_neighbors)` | `n_neighbors = min(max(5, n_obs/500), 10)` |
| UMAP | `sc.tl.umap(aux_ad, min_dist=0.1)` | `min_dist=0.1`，写入 `obsm['X_umap_qc']` |
| Leiden 聚类 | `sc.tl.leiden(aux_ad, resolution=res)` | `res=0.2`，写入 `obs['qc_cluster']` |

产出一个 QC 空间（相似 QC 画像的细胞聚在一起），与生物学表达嵌入无关——这是刻意的：QC 判定只应由质量指标驱动。

### 4.2 簇级投票（clusterwise_qc）

对每个 `qc_cluster` 算细胞级通过比例（`pd.crosstab` 行归一化），**通过比例 ≥ threshold（默认 0.5，由 `--min_frac` 控制）的簇判为 good QC cluster**，簇内全部细胞（含个别 GMM 未过的）继承簇判定，写入 `key_added`（original 模式下为 `good_qc_cluster_mito{X}`）。

官方教程对平滑效果的描述（逐字）："We seem to get a nontrivial number of cells individually flagged as poor quality ... but only a small subset of them remain after cluster membership is leveraged."

## 5. original 模式：mito 三档循环与编码汇总

### 5.1 机制（qc.py `run_qc_mito_loop_original`）

1. `generate_qc_clusters(ad, metrics=qc_metrics)`：单次聚类（8 指标、res=0.2），建 QC UMAP + qc_cluster；
2. 循环 mito 上限 **`MITO_THRESHOLDS = (20, 50, 80)`**：每档把 percent_mito 的 GMM 拟合域上限设为当档值，重跑 `cellwise_qc`（4 指标 GMM）+ `clusterwise_qc`（threshold=0.5），产出 `good_qc_cell_mito{X}` / `good_qc_cluster_mito{X}` / `pass_auto_filter_mito{X}`；
3. 编码汇总：`pass_auto_filter` 与 `good_qc_cluster` 初值 100，按 `reversed((20,50,80))`（即 80→50→20）依次覆盖——**值 = 还能通过的最严 mito 档**：20=三档全过（质量最好）、50=过 80+50 但不过 20、80=仅过 80、100=三档全败（哨兵值）。

### 5.2 4 指标 GMM 拟合域（硬编码于 qc.py，cell/nuclei 双分支）

不传 `--metrics_csv` 即按 `ad.uns['cell_or_nuclei']` 选分支（gather 阶段从样本表第 5 列读入，默认 cell）：

| 指标 | min | max | scale | side | min_pass_rate |
|------|-----|-----|-------|------|---------------|
| **cell 分支** ||||||
| n_counts | 1000 | None | log | min_only | 0.1 |
| n_genes | 100 | None | log | min_only | 0.1 |
| percent_mito | 0.1 | max_mito（循环 20/50/80） | log | max_only | 0.1 |
| percent_spliced | 50 | 97.5 | log | both | 0.1 |
| **nuclei 分支** ||||||
| n_counts | 300 | None | log | min_only | 0.1 |
| n_genes | 100 | None | log | min_only | 0.1 |
| percent_mito | 0.1 | max_mito（循环 20/50/80） | log | max_only | 0.1 |
| percent_spliced | 0.01 | 100 | log | max_only | 0.1 |

> nuclei 分支差异（官方自带，CHANGELOG 24-088 "Added support for single-nuclei samples"）：n_counts 下限放宽到 300（细胞核 RNA 总量低）；percent_spliced 改 (0.01, 100) max_only——细胞核富 pre-mRNA，spliced 天然偏低，故不卡下限只卡上限。解离组织用 cell，单核测序用 nuclei。
> 自定义阈值传 `--metrics_csv`（CSV 格式：index=metric，列=min/max/scale/side/min_pass_rate，官方样例 `example_cutoffs_qc.csv`）；注意 percent_mito 的 max 列会被 mito 循环当档值覆盖。**本 skill 不传 metrics_csv，固定用内置分支。**

### 5.3 为何选 original 而非 multires

multires 固定 mito=20，会把 mito 20–25% 的正常实质细胞（肝/心/肌肉，代谢活跃、线粒体丰富）误判为低质。original 的 mito 循环 + 终判用 ≤80 档对实体组织友好。若数据以免疫/血液为主（mito 普遍 <15%），两者差异不大；含实质组织时 original 更合理。

## 6. 官方能力全景：三种 qc_mode 与六种 run mode

> 本节为官方仓库（README/nextflow.config/qc.py 源码）的能力清单。**本 skill 固定：qc_mode=original、gmm_cutoff=inner、min_frac=0.5、不传 metrics_csv、不启用 celltypist**（§9），其余模式仅作知识备查。

### 6.1 三种 qc_mode（`--qc_mode`，默认 original）

- **original**（默认，Pan-GI atlas 用）：§5 的 mito 三档循环 + 单分辨率簇投票（res=0.2）。
- **multires**：同样的指标计算与 QC 嵌入，但**固定 mito=20**（无循环）；GMM 细胞级判定一次（`cell_passed_qc`），然后多分辨率聚类（0.1–1.0，步长 0.1）各做一次簇投票，记录每细胞的 `consensus_fraction`（通过的聚类分辨率占比）；显式搜索使"失败细胞集合"与细胞级判定的 Jaccard 最大的共识阈值 → `consensus_passed_qc`；退化守卫（consensus 全 ≈0 → 强制全 False）；`keep_multires` = 簇通过 或（共识通过且 consensus_fraction ≥ 0.90 营救线）。产出对单一聚类分辨率不敏感的判定。
- **combined**：把 multires 的多分辨率共识嵌进 mito 三档循环（每档重拟合 GMM + 多分辨率 + Jaccard 阈值 + 0.90 营救 + 退化守卫），并向宽松档传播通过（20→50→80），汇总 `consensus_pass_auto_filter`（同样的 ≤80/100 编码）。

### 6.2 六种 run mode（nextflow `-entry`，官方 Nextflow 流程编排）

| run mode | 步骤链 | 说明 |
|---|---|---|
| `all` | 1-2-3-4-5-6 | 全流程（含 scVI 整合） |
| `only_qc` | 1-2-3 | 到池化前（含双胞） |
| `after_qc` | 4-5-6 | 从池化起（复用既往 QC 产物） |
| `until_integrate` | 1-2-3-4-5 | 到整合前 |
| `only_integrate` | 6 | 仅 scVI 整合 |
| `subset` | 2a-3-4-5a | 不跑自动 QC，按 `--limits_csv` 硬阈值子集（`subset_object` 替换 run_qc；`finalize_qc_basic` 不删任何细胞/样本） |

输入两种格式：Legacy（`--SAMPLEFILE` 样本列表 + `--ss_prefix`/`--cb_prefix`/`--cr_prefix` 路径前缀）或新版 CSV（列 `sampleID` + `path_to_starsolo`/`path_to_cellranger` + 可选 `path_to_cellbender`/`cell_or_nuclei`，官方样例 `example_input.csv`）。nextflow.config 全局默认：`min_frac=0.5`、`ss_out=GeneFull`、`gather_mode=cellbender`、`gmm_cutoff=inner`、`qc_mode=original`、`n_top_genes=7500`、`cellbender_input=cellbender.h5`。

`run_qc` 还可选 CellTypist 注释（`--celltypist_model gut` 预设 5 模型：cecilia22_predH/predL（免疫，32/98 型）、elmentaite21_pred（肠道 134 型）、suo22_pred（胎儿基质+免疫 138 型）、megagut_pred（Pan-GI 89 型），或自定义 .pkl；min_prob=0.3 标 Uncertain）。**本 skill 不启用**（metadata.json 如实记 `celltypist_model=null`）。

## 7. 双胞检测双线：scrublet 与 scDblFinder

### 7.1 scrublet 机制（sctk `run_scrublet` 源码实证）

官方 `flag_doublet.py` 调 `sk.run_scrublet(ad, inplace=False)`，在 per-sample postqc 对象（cellbender 去背景 `.X`）的**全量细胞**上跑：

1. `Scrublet(X).scrub_doublets()`：模拟双胞得 `scrublet_score`（逐细胞 0–1）；
2. 独立标准流程建簇：filter_genes(min_cells=3) → normalize_total(1e4) → log1p → HVG(min_mean=0.0125, max_mean=3, min_disp=0.5, subset) → scale(zero_center=False) → PCA(svd_solver="arpack") → neighbors(n_pcs=30) → leiden(resolution=1) → 按簇大小逐簇再细分（`resolution_function = max(max(log10(簇大小)−1, 0)², 0.1)`）；
3. 簇级检验：每簇取 scrublet_score 中位数 → MAD z 分数（upper_mad_only，只算上侧离群）→ 正态 p 值 → **BH 校正**得 `bh_pval`；
4. 输出 5 列：`scrublet_score`、`scrublet_score_z`、`cluster_scrublet_score`、`bh_pval`、`scrublet_done`；拟合异常（ValueError/ZeroDivisionError）时兜底全 0/p=1 并标 `scrublet_done=False`。
5. `flag_doublet.py` 守卫：若**非零计数基因数 < 20% 总基因数**（`(ad.X.sum(0) > 0).sum() < ad.shape[1] * 0.2`）直接 exit，不产 csv——06 在 finalize 前对每样本 scrublet csv 做存在性硬校验，缺失即显式报错。

> 忠实记录：`flag_doublet.py` 的 CLI 接受 `--filter`/`--groupby` 参数，但当前脚本主体未消费它们——scrublet 始终在全量细胞上运行（06 传的 `--filter cluster_passed_qc` 为无效参数，保留仅为对齐官方调用形）。

### 7.2 finalize 的两列标注（只标注、不过滤）

finalize_qc.py 写入 pooled obs：

- `doublet`（宽松）：`scrublet_score > 0.3 | bh_pval < 0.65`——疑似双胞，假阳性较高
- `stringent_doublet`（严格）：`scrublet_score > 0.3 | bh_pval < 0.05`——FDR 5% 统计显著，确信双胞，假阳性极低

本流程两者都标注、都不删除，保留完整数据。**交付层命名消歧（09 执行）**：scrublet 的 `bh_pval`/`stringent_doublet` 与 scDblFinder 簇级检验新列（§7.4）同义易混，09 合并时把交付对象（pooled/pooled_filtered）obs 中的 `bh_pval`→`scrublet_bh_val`、`stringent_doublet`→`scrublet_stringent`（幂等：旧名在则改、仅新名在则跳过、皆无则报错）；per-sample `<sid>_scrublet.csv` 保留工具原生列名（溯源）。下游整合流程若需去双胞，用保守的 `scrublet_stringent` 删除：`adata = adata[~adata.obs.scrublet_stringent]`（与官方 step6 integrate 的行为一致，见 §13）。

### 7.3 scDblFinder 第二线（skill 附加，与 scrublet 并存）

除 scrublet 外，本流程在相同位置（per-sample postqc 对象）跑 **scDblFinder**（Bioconductor R 包）作为第二条独立 doublet 检测线，输出覆盖与 scrublet 相同的细胞集合（postqc 的 cellbender `.X` 全量细胞），只换内核。

**鉴定输入：阈值闸（06d，Science 2026 口径硬过滤）**：08 运行前由 `06d_threshold_gate.py` 对 postqc obs 的 `n_counts`/`n_genes`/`percent_mito` 做**纯阈值判定**（三项全过为通过；不做 GMM——GMM+initial bounds 是 cellwise_qc（§3）的既有职责，不在此重复；scrublet 不参与任何过滤）。契约默认阈值取自 Science 2026 Hormone Cell Atlas（doi:10.1126/science.aeb2672）的 post-scAutoQC 硬过滤包（≥200 UMI / ≥200 genes / mito<30%，不含其 scrublet<0.3 项），mito 按 sc/sn 分流：

| 指标 | sc | sn |
|---|---|---|
| n_counts | ≥200 | ≥200 |
| n_genes | ≥200 | ≥200 |
| percent_mito | <30 | <5 |

run_config `GATE_MIN_COUNTS`/`GATE_MIN_GENES`/`GATE_MAX_MITO_CELL`/`GATE_MAX_MITO_NUCLEI` 可覆盖（留空=契约默认），实取值入指纹 params（改动即触发重跑）。产出 `2_qc_objects/<sid>_threshold_gate.csv`（barcode + 逐项判定 + `pass_threshold_gate`），兼作断点续跑完成标记。被闸细胞不参与鉴定，其 `score_*` 默认 **1**（供下游按分数剔除——integration 默认 `score_scdblfinder < 0.5` 会顺带移除它们）、`doublet_*` 记 `singlet`（算法未对其判定，类别列不编造 doublet 结论）。闸同时生效于交付对象：09 合并时 `pass_threshold_gate` 写入 pooled obs，pooled_filtered 按闸裁剪（§10.2）。

**双运行（官方 FAQ 同型口径）**：官方 FAQ 原文——"However, should you for some reason try to identify also homotypic doublets with scDblFinder, be sure to not to use the cluster-based approach, and to set removeUnidentifiable=FALSE..."。默认 `removeUnidentifiable=TRUE` 会把"普遍不可识别"组合（近似同型）的人工双胞从训练集剔除，使分类器聚焦异型双胞。每样本平行跑两次：

| 运行 | 调用 | 输出列 |
|---|---|---|
| 默认（异型聚焦） | `scDblFinder(sce, dbr=0.1)` | `score_scdblfinder` / `doublet_scdblfinder` |
| 同型 | `scDblFinder(sce, clusters=FALSE, removeUnidentifiable=FALSE, dbr=0.1)` | `score_scdblfinder_homotypic` / `doublet_scdblfinder_homotypic` |

- **运行**：`08_scdblfinder.R`（per-sample，与 `flag_doublet.py` 同级），`schard` 读 postqc → `as.SingleCellExperiment` → 按阈值闸 mask 子集（`--keep-csv` 消费 06d 产物）→ 两次 `scDblFinder`（固定 `set.seed(1234)`，人工双胞模拟含随机抽样，种子保证可复算）；`09_merge_scdblfinder.py` 合并 per-sample 得分写入 pooled/filtered obs。
- **标注（不过滤）**：score 两列 float（0-1）+ doublet 两列 bool（class==doublet），类型对齐 scrublet 的 `scrublet_score`(float) + `doublet`(bool)。pooled0 不写（finalize 前阶段，与 scrublet 一致）。
- **为何并存两条线**：scrublet（最近邻+模拟）与 scDblFinder（共表达 cxds + 模拟双胞 + xgb 分类器，`score` 参数默认 `"xgb"`）算法独立，得分中度相关、可交叉验证，任一单独方法都有盲区。下游整合时可选其一或取共识（两方法都判 doublet）；同型运行列用于需要审查同型双胞的场景。
- **标准产出**：`scdblfinder_scores.csv`（scautoqc 根部）（barcode-sampleID + 双运行 4 列 + §7.4 簇级检验 4 列）。

> pooled/filtered 的 obs 同时保留 scrublet（5+2 列，交付层改名后 `bh_pval`→`scrublet_bh_val`、`stringent_doublet`→`scrublet_stringent`）+ scDblFinder（双运行 4 列 + 簇级检验 4 列）两套标注，都不参与 `good_qc_cluster` 过滤，删除留给下游整合。

### 7.4 scDblFinder 簇级 BH 检验（06e，sctk 同构）

scDblFinder 原生不产出 per-cell FDR（其 `doubletThresholding` 只给集合级阈值与 FDR 曲线，无逐细胞 p 值列），与 scrublet 的 `bh_pval` 不对等。`06e_scdblfinder_cluster_qc.py` 复用 sctk 的簇级检验补平这一口径——直接 `from sctk._doublets import _calculate_z, _test_outlier`（与 scrublet 侧同一实现，零漂移），对 08 双运行的两个 score 各做一次：

1. **聚类**：在 postqc **全量细胞**上按 sctk `run_scrublet` 的配方逐字重算（filter_genes(min_cells=3) → normalize_total(1e4) → log1p → HVG(0.0125/3/0.5, subset) → scale(zero_center=False) → PCA(svd_solver=arpack) → neighbors(n_pcs=30) → leiden(resolution=1，默认 random_state=0，确定性)）——同一细胞集、同一配方，簇划分与 scrublet 侧一致；该检验与 score 来源无关（score 仅在簇中位数处进入），对 scDblFinder 的 xgb 概率同样成立。06e 的簇与 scDblFinder 内部聚类不同源（独立配方独立重算），检验统计量（簇中位数→MAD z→BH）与 xgb 逐细胞分类机制不同轴——是与分类器误差部分去相关的独立审查轴；`bh_val<0.05` 而 score 低的细胞值得人工复核。
2. **簇级检验**：每簇取 score 中位数 → `_calculate_z`（upper_mad_only，×1.4826）→ 正态 p（1−cdf）→ BH 校正（`_test_outlier`，statsmodels fdr_bh），得 `scdblfinder_bh_val`/`scdblfinder_homotypic_bh_val`。
3. **闸联动**：簇中位数仅取阈值闸通过的成员——被闸细胞 score=1 是 08 的回填约定值而非实测，入中位数会污染统计；被闸细胞与零通过成员簇的 bh_val 记 1。
4. **stringent**：`scdblfinder_stringent = (score_scdblfinder > 0.5) | (scdblfinder_bh_val < 0.05)`（同型列同理，`scdblfinder_homotypic_stringent`）。不用 class 列合成——被闸细胞 class=singlet，按 class 合成会让被闸细胞逃逸；score 阈值 0.5 对齐下游整合的 scDblFinder 去除线（§13），bh<0.05 对齐 scrublet stringent 的 FDR 口径（§7.2）。
5. **兜底**：聚类或检验失败（ValueError/ZeroDivisionError/MAD 失效）→ 该运行 bh_val 全 1 + 日志 `[warn]`，不中断流程（scrublet 侧同设计，§7.1 第 4 条）；此时 stringent 退化为纯 score>0.5 判定。

产出 `3_doublet_scores/<sid>_scdblfinder_cluster_qc.csv`（barcode + 上述 4 列，覆盖 postqc 全量），纳入断点续跑完成标记与交付校验；09 合并写入 pooled/pooled_filtered obs（§10.2）。只标注不过滤，与全线 doublet 列一致。

## 8. 池化与终判（pool_all / finalize_qc）

### 8.1 pool_all（≥2 样本才执行）

- `sc.AnnData.concatenate(batch_key="sampleID", batch_categories=samples)` 合并全部 postqc 对象 → `scautoqc_pooled0.h5ad`；
- var 跨样本汇总：各层 `n_counts{suffix}-{sid}`/`n_cells{suffix}-{sid}` 列求和写入 `var`（并留 `varm` 矩阵），逐样本后缀列删除；
- per-sample `{sid}_qc_thresholds.csv` 拼接为根部 **`qc_thresholds.csv`**（tidy 长表，全样本全指标的实际 GMM 界与通过率）。

### 8.2 finalize_qc（original 模式过滤逻辑）

- **基因列表标记**（`BASE_DIR/genes_list/`）：`var['cc'/'ig'/'tcr']`（细胞周期/免疫球蛋白/TCR，供下游 scVI 移除周期基因；不参与核心过滤）。注意 finalize 内部用正则 `^IGH[ADEMG][0-9]*$` 拆 IG 恒定区（igc），此为人类格式，小鼠（Ighg1/Ighm）不匹配——仅使 `var['ig']` 含恒定区（细微差异），不影响 QC 过滤。
- **样本级汇总** `sample_passqc_df.csv`：per-sample 三档 frac/count（20/50/80）+ total_cell_count。
- **样本级剔除**：`frac80<0.1 且 count80<100` 的样本归入 `samples_to_check2`（通过率<10% 且通过细胞<100 的低质样本）；另打印两组提示名单：check1（count20>10000）与 check3（count80==0，零通过），仅打印不剔除。
- **细胞级保留**：`good_qc_cluster≤80`（至少过了最宽的 mito80 档）且 sampleID 不在 samples_to_check2 → `scautoqc_pooled_filtered.h5ad`。
- **产出**：`scautoqc_pooled.h5ad`（过滤前，含全部标注）、`scautoqc_pooled_filtered.h5ad`（过滤后）、`sample_passqc_df.csv`、两张 overall 图（`scatterplot_sample_passqc.png`：total vs count80 symlog 散点；`barplots_sample_passqc.png`：三档 count 的 log10 直方图）→ 归入 `5_qc_plots_overall/`。
- 其他模式过滤口径（备查）：combined → `consensus_pass_auto_filter≤80`；multires → `keep_multires`。

> pooled_filtered **只按 good_qc_cluster 过滤，不删 doublet**。doublet/stringent_doublet 仅作 obs 列标注，删除留给下游整合步骤。

### 8.3 genes_list 物种适配（05_prepare_genes_list.py）

finalize_qc.py 需要三个基因列表；nf-scautoqc 自带的是**人类全大写符号**，对小鼠（首字母大写）几乎全失配。05 脚本从 skill 内置 `assets/ref/{species}` 提取并写入 `$NF_SCAUTOQC_DIR/genes_list/`：

- mouse：周期←`cell_cycle_genes_updated_2019.csv`(94)；IG/TCR←`mart_export.txt` 按 `Gene type` 筛 IG_/TR_ 非 pseudogene
- human：周期←`JP_cycle_genes.list`(559，直接复制)；IG/TCR←`mart_export.txt`
- custom：`--src-cc` 与 `--src-mart` 自备来源

### 8.4 元数据记录（metadata.json，06 的 5c 段）

finalize 后、写 HANDOFF 前，06 对每个样本落定 `outs/metadata.json`（交付契约第 5 项）：

1. `01c_write_metadata.py`：十类基础字段（group/age/batch/datasets/sample/tech/gender/tissue/chemistry/species）。**已存在则整体跳过**（不覆盖）；缺失则以 `Unknown` 兜底创建。tech/chemistry 按 ALIGNER 推导（cellranger→`sc`/`10x_3prime_v3`；mobivision→`mobivision`/`mobi_3prime_v1`），可被 run_config `TECH`/`CHEMISTRY` 覆盖；tissue 缺省 `Unknown`（run_config `TISSUE` 可填）；species 由 `--species`（通常传 run_config `SPECIES`）写入，不从样本名推断。per-sample 字段（group/age/gender/batch）来自配置阶段整理的 `--meta-csv`（列 `sample_id,group,age,gender,batch[,sample][,species]`，见 SKILL.md 配置步骤 3）；gender 规范 `F`/`M`，01c 读表时对男/女/Male/Female 等常见写法自动归一（大小写不敏感），其余值原样保留不猜测。
2. `06c_record_stage_genes.py`：写入 `stages.scautoqc.celltypist_model = null`——本流程**不启用** celltypist 注释，如实记录 null 占位（交付契约"阶段所用基因（如有）"：无则记 null，不编造）。

## 9. 本 skill 执行层（06_scautoqc_pipeline.sh）

### 9.1 为什么需要桥接（04_bridge_to_starsolo.py）

nf-scautoqc 的 `gather_matrices.py` 原生支持 STARsolo 输出 + CellBender + Cell Ranger，但有四处不兼容我们的数据组合（cellranger/mobivision + velocyto.py loom + cellbender）：

1. **`sk.read_velocyto` 只认 STARsolo velocyto 目录**（`spliced.mtx.gz`/`unspliced.mtx.gz`/`ambiguous.mtx.gz` + barcodes + features），**不认标准 velocyto.py 的单个 `.loom`**
2. cellranger+cellbender 模式会读 `raw_feature_bc_matrix`（~1.9M droplet）去索引 velocyto（几千细胞），barcode 集合不匹配
3. cellranger 与 velocyto 的 barcode 后缀/前缀不一致（`-1`、`sample:`、`x`）
4. **cellbender h5 与 bridge 产物的 features 命名约定系统性错配**：cellbender 输入 h5 由 `make_10x_h5.py` 按 write10xCounts 语义生成（features/id=name=符号，重复符号按 R make.unique 追加 `.1`/`.2`），cellbender 原样透传到 `cellbender_matrix.h5`；而 bridge 的 `Gene/filtered` 经 anndata 唯一化是连字符风格（`X-1`）。gather（starsolo 模式）要求 cb 与 cr_gene 的 var_names 集合严格相等——人类参考的 22 个重复符号基因两风格对不上，直接 `ValueError: raw features do not match main matrix (22 missing, 22 extra)`。前 3 处由 04 桥接解决，**第 4 处由 06 在 bridge 前调 `06a_align_cellbender_features.py` 自动检查+修复**（幂等：在 tsv 符号序列上重现 R make.unique 与 h5 name 双射配对，逐行写回 features/id←ENSG、features/name←tsv 原始符号；只写 staging 副本，矩阵数值/barcodes/droplet_latents 零改动，修复后 sctk.read_cellbender 的 var_names 与 bridge 侧集合必然相等）

桥接脚本把 cellranger filtered + velocyto loom 转成 STARsolo 目录结构：

- `Gene/filtered/`：cellranger filtered 矩阵（barcode 去 `-1`）
- `Velocyto/`：loom 拆成 spliced/unspliced/ambiguous 三个 mtx（barcode 去 `sample:` 前缀和 `x` 后缀）
- 基因用 **Ensembl gene_id**（cellranger 的 `gene_ids` ↔ velocyto 的 `Accession`）对齐，规避重复基因符号（小鼠 Gm* 有 37 个重复）

然后用 `gather_matrices.py --gather_mode starsolo` 消费，产出 `.X`=cellbender 去背景、layers=raw/spliced/unspliced/ambiguous 的 h5ad。

### 9.2 步骤链与固定命令（全部固定，不经命令行改参）

```
[0]  05_prepare_genes_list.py --species $SPECIES --out-dir $NF_SCAUTOQC_DIR/genes_list
     ↓ 样本表读取 + cell_or_nuclei 一致性硬校验（混杂即 exit 1）
[1-4] per-sample 循环：
  0. 输入硬校验（mtx 三件套 + loom + cellbender_matrix.h5，缺一即 exit 1）
  0b. features 对齐：06a_align_cellbender_features.py --cr-root $CR_ROOT --samples <sid>
      （§9.1 第 4 处错配的守卫；已一致的样本秒过，错配的自动修复 staging 副本，不可修复即 exit 1）
  1. bridge：04_bridge_to_starsolo.py --cr_mtx <filtered_feature_bc_matrix> --velo_loom <sid>.loom --outdir 1_bridge/<sid>
  2. gather（官方 bin/gather_matrices.py，工作目录 1_bridge/<sid>）：
     --cr_gene ./Gene/filtered --cr_velo ./Velocyto --cb_h5 <cellbender_matrix.h5>
     --gather_mode starsolo --cell_or_nuclei <样本表第5列，默认 cell>
     → 产出移归 1_gathered_objects/<sid>.gene_velo_cellbender.filtered.h5ad
  3. qc（skill 的 qc.py = 官方版 + 物种无关 flags，另含绘图细节与 pandas 2.x 兼容修补，人源数据行为与官方版等价；工作目录 2_qc/<sid>）：
     --sample_id <sid> --gath_obj <gathered.h5ad> --qc_mode original --gmm_cutoff inner
     （不传 metrics_csv/celltypist_model/min_frac → min_frac=0.5 默认）
     → <sid>_postqc.h5ad 移归 2_qc_objects/；6 张 qc_plot PNG 移归 2_qc_plots_individual/；
       <sid>_qc_thresholds.csv 留 2_qc/<sid>/（pool_all 消费）
  3b. 阈值闸（06d_threshold_gate.py；Science 2026 口径硬过滤，纯阈值无 GMM，§7.3）：
     --postqc 2_qc_objects/<sid>_postqc.h5ad --sid <sid>
     --out 2_qc_objects/<sid>_threshold_gate.csv --cell-or-nuclei <样本表第5列>
     （阈值默认 counts>=200/genes>=200/mito<30(sc)/<5(sn)；run_config GATE_* 覆盖，入指纹）
  4. flag_doublet（官方 bin/flag_doublet.py，工作目录 3_doublet_scores/）：
     --input <postqc.h5ad> --samp <sid> --filter cluster_passed_qc → <sid>_scrublet.csv
     （--filter 被 CLI 接受但脚本主体未消费，实际在全量细胞上跑，见 §7.1）
  4b. scDblFinder 双运行（08_scdblfinder.R，LC_CTYPE=UTF-8 前缀 + RSCRIPT 指定实例）：
     --postqc <postqc.h5ad> --sid <sid> --keep-csv 2_qc_objects/<sid>_threshold_gate.csv
     --out 3_doublet_scores/<sid>_scdblfinder.csv --scanpy-python $SCANPY_PYTHON（dbr=0.1）
     鉴定输入=阈值闸通过细胞（被闸细胞 score=1/class=singlet）：默认运行 +
     同型运行（clusters=FALSE, removeUnidentifiable=FALSE），输出 4 列（§7.3）
  4c. scDblFinder 簇级 BH 检验（06e_scdblfinder_cluster_qc.py；sctk 同构，§7.4）：
     --postqc <postqc.h5ad> --sid <sid>
     --scdblfinder-csv 3_doublet_scores/<sid>_scdblfinder.csv
     --gate-csv 2_qc_objects/<sid>_threshold_gate.csv
     --out 3_doublet_scores/<sid>_scdblfinder_cluster_qc.csv
     （sctk 配方重算簇 → 簇中位数 z→BH；被闸细胞不入中位数、bh_val=1；
       stringent=score>0.5 或 bh_val<0.05；失败兜底 bh_val=1 仅警告不中断）
[幂等] 每样本指纹 _runmeta/<sid>.scautoqc.fingerprint.json：脚本 sha（06 + 06d + 06e）/参数
     （qc_mode/gmm_cutoff/cell_or_nuclei/species/GATE_* 四阈值）/输入五件
     （filtered mtx 三件 + loom + cellbender_matrix.h5）
     size+mtime/工具版本（scanpy/sctk/scrublet）逐项一致 + 产物齐备才 skip
    （--require 核验五件：postqc/threshold_gate csv/scrublet csv/scdblfinder csv/
     scdblfinder_cluster_qc csv；编排层完成标记六件另含 1_gathered_objects/ 的
     gathered h5ad）；否则整样本重跑（marker 库、参数、参考、脚本版本任一变化触发）
[4-5b] 池化段（≥2 样本才执行；合并级指纹 _runmeta/merged_scautoqc_pooled_filtered.fingerprint.json
     = 各样本 scautoqc 指纹成员集合 + 合并参数，逐项一致则整段跳过）：
  4. pool_all（官方）：--samples <逗号列表> --objects <postqc 列表> --ranges <qc_thresholds 列表>
     → 根部 scautoqc_pooled0.h5ad + qc_thresholds.csv
  5. finalize_qc（官方）：--obj scautoqc_pooled0.h5ad --scr <scrublet csv 列表> --qc_mode original
     → scautoqc_pooled.h5ad / scautoqc_pooled_filtered.h5ad / sample_passqc_df.csv；
       overall PNG 移归 5_qc_plots_overall/（finalize 前 06 先硬校验每样本 scrublet csv 存在）
  5b. 09_merge_scdblfinder.py --sc-root $SC_ROOT --samples <列表>
[5c] metadata.json：01c_write_metadata.py（十字段，含 species；已存在跳过）+ 06c_record_stage_genes.py
[6]  07_write_handoff.py 写 HANDOFF 契约段（多数据集时 --section $DATASET_NAME）
[7]  清理 __pycache__/.Rhistory/.RData（中间产物 1_bridge/2_qc 保留可断点重跑）
```

### 9.3 目录布局与日志

`SC_ROOT` 默认 `$OUTPUT_ROOT/$PROJECT_NAME/scautoqc`（run_config `SC_ROOT` 可覆盖），命名完全对齐官方：

```
$SC_ROOT/
├── 1_gathered_objects/                 # 官方 output 1（per-sample gathered h5ad）
├── 2_qc_objects/                       # 官方 output 2（per-sample postqc h5ad + skill 附加 <sid>_threshold_gate.csv）
├── 2_qc_plots_individual/              # 官方 output 3（per-sample 6 张 QC PNG）
├── 3_doublet_scores/                   # 官方 output 4（scrublet csv + skill 附加 scdblfinder csv + scdblfinder_cluster_qc csv）
├── 5_qc_plots_overall/                 # 官方 output 8（scatterplot/barplots sample_passqc）
├── qc_thresholds.csv                   # 官方 output 5（全样本 GMM 阈值 tidy 长表）
├── sample_passqc_df.csv                # 官方 output 9（样本级三档通过率）
├── scautoqc_pooled.h5ad                # 官方 output 6（过滤前全标注）
├── scautoqc_pooled_filtered.h5ad       # 官方 output 7（质控后对象，供下游）
├── scdblfinder_scores.csv              # skill 附加（scDblFinder 合并得分）
├── 1_bridge/<sid>/                     # 中间产物（STARsolo 结构，保留可重跑）
├── 2_qc/<sid>/                         # 中间产物（per-sample qc_thresholds.csv 等）
└── 日志：1_bridge_<sid>.log / 2_qc_<sid>.log / 3_doublets_<sid>.log /
    pool_all.log / finalize_qc.log / scdblfinder_merge.log
```

指纹记录（指纹判据）落 `$CR_ROOT/_runmeta/`（`<sid>.scautoqc.fingerprint.json`、`merged_scautoqc_pooled_filtered.fingerprint.json`），与 SC_ROOT 平级、不随 SC_ROOT 覆盖而漂移；`09b_check_delivery.sh` 交付收口时经 `runmeta.py verify` 校验齐备与成员一致。

### 9.4 单样本数据集（止步于多样本处理前，有意设防）

scautoqc 的 4-5b 段（pool_all 合并 → finalize_qc 池化过滤）是**多样本处理**——为下游 scVI/harmony 等多样本整合备料。数据集只有 1 个样本时，"合并"没有意义：对单样本跑整合是空跑（0/0 式无意义形式化）。**这不是软件缺陷，而是有意设防**：06 在多样本处理前止步，打印显眼警告，仅交付单样本质控结果（`1_gathered_objects/`、`2_qc_objects/`、`2_qc_plots_individual/`、`3_doublet_scores/`、`qc_thresholds.csv`）+ metadata.json（5c 段照跑），不产出 `scautoqc_pooled*.h5ad` 三件套、不写 `## upstream 阶段` 契约段落（该段落承诺三个 pooled h5ad，写了就是谎报）。

**请从根本上避免单样本数据集进入本阶段**：

- 合并多个数据集为一个数据集（≥2 样本）再跑；或
- 全程维持单样本分析路径（不进入 scautoqc 池化/整合）。

机制注记（供理解守卫位置，非缺陷陈述）：官方 `pool_all.py` 依赖 `AnnData.concatenate` 产生的 `-<sid>` 后缀 var 列，而 concatenate 单批次不加后缀，强行跑会在求和处 `KeyError`——06 的守卫让止步发生在烧算力之前，并给出明确指引而非栈追踪。

## 10. 输入与产出契约

### 10.1 输入（每样本 3 组，06 开头硬校验，缺一即 exit 1）

| 输入 | 路径 | 来源阶段 |
|---|---|---|
| 表达矩阵三件套 | `<sid>/outs/filtered_feature_bc_matrix/{barcodes,features,matrix}.mtx.gz`（mobivision 分支经 02b 规范化） | 阶段1 比对 |
| 剪接 loom | `<sid>/velocyto/<sid>.loom` | 阶段2 velocyto |
| cellbender 输出 | `<sid>/outs/cellbender/cellbender_matrix.h5`（gather 的 `--cb_h5` 输入：sctk.read_cellbender 读其矩阵与 cell_probability/droplet_efficiency） | 阶段3 cellbender |

另需：样本表 `SAMPLE_SHEET`（第 5 列 `cell_or_nuclei`，空=cell；全表必须同值，sc/sn 混批拒跑）、`NF_SCAUTOQC_DIR`（官方仓库 clone）、`ENV_SCAUTOQC`（conda 环境，scanpy 1.10.4 + sctk + scrublet + celltypist）、`RSCRIPT`/`SCANPY_PYTHON`（scDblFinder 线；R 调用带 `LC_CTYPE` UTF-8 前缀，见环境手册）、run_config 十字段元数据来源（`META_CSV`/`TECH`/`CHEMISTRY`/`TISSUE`/`SPECIES`）。

### 10.2 产出语义（两份 pooled 对象的分工）

- **`scautoqc_pooled.h5ad`**（过滤前）：obs 全标注——8 项 QC 指标、三档 mito 的细胞/簇判定（`good_qc_cell_mito{20,50,80}` 等）、`pass_auto_filter`/`good_qc_cluster` 编码、scrublet 5 列 + doublet/scrublet_stringent（交付层改名后，§7.2）、scDblFinder 8 列（双运行 4 列 + 簇级检验 4 列，§7.3/§7.4）、`pass_threshold_gate`（阈值闸判定）、cell_probability/droplet_efficiency；layers 已删（finalize_qc 以 X/obs/var 重建对象、不携带 layers——原始 raw/spliced/unspliced/ambiguous 四层保留在 `scautoqc_pooled0.h5ad`）。**含流程全部处理信息，可自定义阈值做下游处理**。
- **`scautoqc_pooled_filtered.h5ad`**（过滤后）：完成**指标 QC**（good_qc_cluster≤80 + 样本级剔除）**+ 阈值闸硬过滤**（`pass_threshold_gate=TRUE`；Science 2026 口径，§7.3），**保留全部双胞标注未删**——可直接手动一步 `adata[~adata.obs.scrublet_stringent]` 后进入下游分析。

## 11. 与 Oliver et al. 2024 方法的对应关系

论文 Methods 的两个关键句（逐字）：

> "For the atlas, the number of GMM components was set to 10 for an overfit model. scAutoQC was subsequently improved to automate the best model fit between 1 and 10 components based on the Bayesian information criterion."

> "Then, whole clusters were classified as passing QC if 50% or more of individual cells within the cluster passed QC."

**逐值一致项**（本流程 = 论文 scAutoQC 方法）：8 指标集合、QC 嵌入参数链（PCA n_comps=max(2,m−2)、n_neighbors=min(max(5,n/500),10)、UMAP min_dist=0.1、Leiden res=0.2）、4 指标 GMM 拟合域（cell 分支）、mito 三档 [20,50,80]、GMM threshold=0.05、cutoff=inner、簇 ≥50% 多数表决（min_frac=0.5）、编码汇总规则、pool 后 `good_qc_cluster≤80` 过滤、`frac80<0.1 且 count80<100` 样本级剔除。

**差异项**（均可溯源，非偏离）：

| 差异 | 论文 atlas | 本流程 | 性质 |
|---|---|---|---|
| GMM 组件数 | 固定 10（overfit model） | sctk 现行 BIC 自动选 1–10 | 论文自认的后续改进版；论文明确参数"只起引导作用、对起始点不敏感"，簇级 50% 表决吸收组件数级别的拟合差异 |
| 双胞处理 | 整合步骤删除 stringent 双胞 | 只标注不删；`scrublet_stringent`（交付层改名自 `stringent_doublet`；`scrublet_score>0.3 \| bh_pval<0.05`）列等价对应，删除留给下游整合（官方 step6 同样用 stringent 删） | 数据保留更全，删与不删由下游决定 |
| scDblFinder 第二线 | 无 | 有（默认+同型双运行，标注不过滤） | skill 附加的独立验证线 |
| 细胞级硬过滤 | 无（细胞级只标注，删除走簇级/样本级） | 有——阈值闸（06d，Science 2026 口径：counts≥200/genes≥200/mito<30% sc/<5% sn，纯阈值无 GMM），闸失败细胞从双胞鉴定输入与 filtered 交付对象中排除 | skill 附加层（§7.3）；与 Oliver 2024 的标注-only 策略差异明确可溯源 |
| 基因标记正则 | 人类大写 `^MT-` 等 | 物种无关 `^[mM][tT]-` 等 | skill 修补；人源数据零实质影响，鼠源必要修正 |
| celltypist | atlas 用 gut 模型组 | 不启用（metadata 记 null） | 注释非质控环节，留给下游 |
| 平台 | 10X | 10X / MobiNova 双平台 | 非方法差异 |

## 12. 故障排除

- **布尔列绘图崩 / numpy2 import 崩**：scanpy 须 1.10.4（非官方 Dockerfile 的 1.9.3），原因与环境组装见 [environment-setup.md](environment-setup.md)。
- **gather 报 `raw features do not match main matrix (N missing, N extra)`**：cellbender h5 的重复基因符号是 R 点号去重（`X.1`）而 bridge 侧是 anndata 连字符去重（`X-1`）——§9.1 第 4 处的系统性约定错配（重复符号基因数量随参考注释版本而定，人类 GRCh38 为 22 个）。06 已在 bridge 前经 `06a_align_cellbender_features.py` 自动检查+修复（幂等、只写 staging 副本、矩阵数值零改动）；手工流程绕开 06 时需先手动跑该脚本（`--check-only` 为只读预检）。两种配对模式（canon 点号双射 / ENSG 直配）都失败时脚本拒修并报错——此时不要强行改名，人工核查该样本 features 来源。
- **样本表 cell_or_nuclei 混杂**：06 开头硬校验，sc/sn 混批直接 exit 1——拆分样本表分批运行（数据集唯一性硬规则，见 SKILL.md）。样本表为 5 列 `sample_id,fastq_source_dir,fastq_prefix,bam_archive_dir,cell_or_nuclei`，06 只消费第 1、5 列并剥离 `\r`/空格；**自建校验/读取脚本常见坑是漏剥 CRLF/BOM 导致误报**，直接用 06 自带校验为准。
- **输入缺失**：06 逐样本校验 3 组输入（mtx 三件套/loom/cellbender_matrix.h5），缺失报 `[err] <sid> 输入缺失: <路径>`，先补齐上游 format/velocyto/cellbender。
- **scDblFinder 报 size factor/拟合错误**：官方 FAQ——输入不得含空 droplet 或极低覆盖细胞（<200/500 reads 可致报错）。鉴定输入已由 06d 阈值闸把关（counts≥200/genes≥200/mito 分流，§7.3）；若仍报错，查 3_doublets_<sid>.log 中 `[scdblfinder]` 行的阈值闸剔除数与 2_qc_objects/<sid>_threshold_gate.csv 的逐项判定，评估该样本低覆盖占比是否异常；确属阈值不适配时用 run_config `GATE_*` 校准（改动入指纹即触发重跑），质量过低的样本由 finalize 的样本级剔除处理。
- **06e 簇级检验打 `[warn]` 兜底**：聚类或 MAD 检验失败时该运行 bh_val 全 1 并继续（设计行为，§7.4 第 5 条；scrublet 侧同设计）——bh_val 全 1 即"该运行无簇级证据"，stringent 退化为纯 score>0.5 判定。若某样本持续兜底，查 3_doublets_<sid>.log 的 `[warn]` 行：常见原因是细胞数过少致聚类/中位数退化，或 sctk 版本漂移（`_calculate_z`/`_test_outlier` 签名变化会被 import 守卫显式报出）。
- **06 失败即打印日志尾部**：per-sample 各步与 pool/finalize 失败时自动打印对应日志尾部 30 行并 exit 1——不再需要事后翻 per-sample 日志定位 ValueError 被静默吞掉的位置。
- **scrublet csv 缺失**：flag_doublet 的 20% 基因守卫可能跳过产出；06 在 finalize 前显式校验并报 `[err] 缺 scrublet csv`，不会放任 finalize 内部才暴错。
- **GMM 拟合失败/通过率过低**：cellwise_qc 三级兜底（§3.4）自动回退硬界，日志可见 "WARNING: Failed to find lower bound" / "No cells passed" 字样——属官方设计行为，核查对应指标分布即可。
- **单样本数据集**：06 止步于多样本处理前并打印警告（§9.4），这不是失败；按指引合并数据集或走单样本路径。
- **pool_all KeyError（var 后缀列）**：单样本强跑 concatenate 的典型表现——06 守卫已防止；若绕开 06 手跑官方脚本遇到，即此原因。
- **Windows Git Bash 调 wsl.exe 的参数被改写**：MSYS 会把以 `/` 开头的参数（`/mnt`、`/home`、`~`）重写成 Windows 路径（典型报错：`bash: 行 1: C:/Program: 没有那个文件或目录`、`~` 被展开为 Windows 主目录）。规则：`export MSYS2_ARG_CONV_EXCL="*"`；wsl 参数一律用绝对 `/home/...`、`/mnt/...` 路径（不用 `~`）；复杂命令写成脚本文件 `cp` 进 `/tmp`、`sed -i "s/\r$//"` 后执行，**不写经多层 shell 的 awk/sed 嵌套引号 one-liner**（引号必坏，表现为脚本静默跑 0 个样本）；Windows 侧带空格路径（`C:/Program Files/...`）必须整体加引号。
- **rsync 归档校验误报**：`rsync -nci src/ dst/` 对顶层目录固定打印 `skipping directory .`——这是常态输出不是差异；自建校验逻辑须排除该行再判定（`rsync -nci ... | grep -v 'skipping directory'`），否则会把全量一致误判为 ARC_FAIL。
- **批量运行时长预期**：per-sample 串行（bridge→gather→qc→阈值闸→scrublet→scDblFinder→簇级检验），每样本数分钟至数十分钟（视细胞数），整批可达数小时；编排日志归 `$OUTPUT_ROOT/$PROJECT_NAME/logs/orchestration/`。

## 13. 下游衔接

本流程止于官方 step5 finalize_qc，**不执行 step6 scVI 整合**；下游去批次由下游整合流程承接。官方 step6（integration.py，备查口径）：删 `stringent_doublet` → 移除 `var['cc']` 周期基因 → `sc.pp.highly_variable_genes(flavor="seurat_v3", n_top_genes=7500, subset=True)` → `SCVI.setup_anndata(batch_key=sampleID, continuous_covariate_keys=["log1p_n_counts", "percent_mito"])` → scVI(n_latent=20, n_layers=2, use_layer_norm="both", use_batch_norm="none", encode_covariates=True, dropout_rate=0.2) → train(train_size=0.9, early_stopping_patience=45, max_epochs=400, batch_size=256, limit_train_batches=100) → neighbors(use_rep="X_scvi") → UMAP(min_dist=0.3) → `scautoqc_integrated.h5ad` + `elbo_training.png` + `scvi_model.pkl`。

交付对象：`scautoqc_pooled_filtered.h5ad`（质控后、双胞未删）+ `scautoqc_pooled.h5ad`（全信息备查）；下游去双胞用 `scrublet_stringent`（交付层改名自 `stringent_doublet`，§7.2）或 scDblFinder 线的 `scdblfinder_stringent`（§7.4）。HANDOFF 契约段落与 formatted 读侧见 [handoff-contract.md](handoff-contract.md)。
