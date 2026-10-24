# HANDOFF.md 字段约定

Omichestra-upstream 在 scautoqc 阶段末步往项目单一 `Document/HANDOFF.md` 段级 upsert `## upstream 阶段` 段落，下游流程读该段落。本文档定义该段落的字段约定（段落 body 子标题用 `###`，避免与段落边界 `##` 冲突）。

## 段落字段约定（## upstream 阶段 body）

以下为 `07_write_handoff.py` 的实际产出格式（逐字对齐）。HANDOFF.md 文件不存在时脚本先写文件级标题 `# HANDOFF — 项目全流程记录`，再追加本段落；已存在则只 upsert 本段落。

```markdown
## upstream 阶段
生成时间: <YYYY-MM-DD HH:MM:SS> | skill: Omichestra-upstream
项目: <PROJECT_NAME>

### 技术信息
- species: <SPECIES>
- tech: sc
- chemistry: 10x_3prime_v3
- barcode_format: barcode-sampleID

### 样本
| sample_id |
|---|
| <sid1> |
| <sid2> |

### scautoqc 产出
根目录: <绝对路径>/scautoqc
| 文件 | 细胞数 | .X | layers | QC列 |
|---|---|---|---|---|
| scautoqc_pooled0.h5ad | <全部> | cellbender | raw/spliced/unspliced/ambiguous | - |
| scautoqc_pooled.h5ad | <全部> | cellbender | (已删) | scrublet系7列(scrublet_score/scrublet_score_z/cluster_scrublet_score/scrublet_bh_val/scrublet_done/doublet/scrublet_stringent)/scDblFinder系8列(score_scdblfinder/doublet_scdblfinder/score_scdblfinder_homotypic/doublet_scdblfinder_homotypic/scdblfinder_bh_val/scdblfinder_homotypic_bh_val/scdblfinder_stringent/scdblfinder_homotypic_stringent)/pass_threshold_gate/good_qc_cluster/pass_auto_filter |
| scautoqc_pooled_filtered.h5ad | <过滤子集（簇级 QC ∩ 阈值闸）> | cellbender | (已删) | 同pooled(子集；簇级 QC + 阈值闸双过滤) |

### 关键说明
- pooled0 与 pooled 的 obs_names 顺序完全一致（可按位置合并 layers+QC列）
- .X = cellbender 去背景整数 counts（未 normalize/log）
- barcode 格式 = barcode-sampleID（如 AAACCCAAGTGTGTTC-<sid>）
- scdblfinder_scores.csv：scDblFinder 合并表（双运行 4 列 + 簇级检验 4 列），已写入 pooled/filtered obs；pooled0 不含（finalize 前阶段）
- 簇级 BH 检验（06e，sctk 同构）：scDblFinder 双运行各配 `*_bh_val`（簇得分中位数→z→BH-FDR，簇中位数仅取阈值闸通过细胞）与 `*_stringent`（score>0.5 或 bh_val<0.05）；被闸细胞 bh_val=1、stringent=TRUE（score=1 回填所致）；聚类在 postqc 全量上按 sctk 配方重算（leiden 分辨率1，确定性随机种子），失败兜底 bh_val=1 仅警告不中断
- 交付层字段命名消歧：交付对象中 scrublet 的 bh_pval/stringent_doublet 已更名 scrublet_bh_val/scrublet_stringent；per-sample scrublet csv 保留工具原生列名（溯源）
- 阈值闸（06d，Science 2026 口径硬过滤）：n_counts>=<实取值>、n_genes>=<实取值>、percent_mito<<实取值>（sc）/<<实取值>（sn）；三层细胞数 = pooled 全量 <n1> → 闸通过 <n_gate> → filtered 交付 <nf>（簇级 QC ∩ 阈值闸）；被闸细胞 scDblFinder score=1 回填，pooled obs 含 pass_threshold_gate 列
```

## 字段说明（formatted 依赖）

| 字段 | formatted 如何用 |
|------|-----------------|
| `项目` | 下游解析脚本记录到 handoff.json（下游脚本不消费该值） |
| `species/tech/chemistry` | 写入 obs 对应列（与样本元数据表交叉校验） |
| `barcode_format` | 校验为 `barcode-sampleID`（改名正则的前提） |
| `### 样本` 表的 `sample_id` | 改名正则的 sampleID 集合 + 拆分单样本的循环 |
| `### scautoqc 产出` 表 | 三个 h5ad 的绝对路径（下游解析脚本以 `根目录` + 文件名拼接） |
| `### 关键说明` | 人读，提醒 pooled0/pooled 对齐、.X 是 cellbender |

## 技术平台值（tech/chemistry）

`06_scautoqc_pipeline.sh` 调 `07_write_handoff.py`，**按 `ALIGNER` 自动推导并传参**，无需手工干预：

- cellranger 分支（默认）：`tech: sc`、`chemistry: 10x_3prime_v3`
- **MobiVision 分支**（`ALIGNER=mobivision`）：`tech: mobivision`、`chemistry: mobi_3prime_v1`（5' 试剂盒则把 run_config `CHEMISTRY` 填为 `mobi_5prime_v1` 覆盖；`TECH` 同理可覆盖 tech）。barcode 格式仍为 `barcode-sampleID`（loom 内 barcode 带 `x` 后缀、mobivision 矩阵 barcode 无 `-1` 后缀，桥接脚本统一处理）。MobiVision 平台说明见 [mobivision-platform.md](mobivision-platform.md)

run_config `TECH`/`CHEMISTRY` 留空走上述推导，填写则覆盖；两个 HANDOFF 段落与每样本 `outs/metadata.json` 用同一组推导值，三处一致。

## 上游非 Omichestra-upstream 时

若 scautoqc 是手工跑的（无 upstream），照本模板手写一份 HANDOFF.md，填实际项目信息即可。formatted 不关心上游是谁，只解析此文档。
