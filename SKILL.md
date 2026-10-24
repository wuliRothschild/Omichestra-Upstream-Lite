---
name: omichestra-upstream
description: 单细胞 RNA-seq 上游模块（scautoqc 版）。10x / MobiNova 单细胞 RNA-seq 上游 + 自动质控：fastq → 比对定量（cellranger count 或 MobiVision quantify）→ velocyto（剪接）→ cellbender（去背景）→ nf-scautoqc（GMM 自动质控 + 阈值闸 + scrublet/scDblFinder 双线双胞 + 簇级 BH 检验 + 池化），产出可直接下游整合/注释的 scautoqc_pooled_filtered.h5ad。当用户要做单细胞上游分析、fastq 转矩阵、cellranger 比对、mobivision/墨卓/MobiNova 数据处理（quantify/fmt_gtf/mkindex）、velocyto、cellbender、scautoqc/scAutoQC 自动质控、双胞过滤、或把一批 fastq 端到端转成下游可用对象时使用。即使只提到其中一个环节（"帮我跑 cellranger"、"跑 mobivision quantify"、"做下质控"），只要处于这条流水线语境也应触发。
---

# Upstream pipeline（scautoqc 版）

把一批单细胞 RNA-seq 样本从 fastq 端到端转成下游可直接用的质控 h5ad。双平台由 `run_config.sh` 的 `ALIGNER` 选择：`cellranger`（10x）或 `mobivision`（墨卓 MobiNova）。平台差异集中在阶段 0–2，阶段 3–4 共用；质控统一走 `06_scautoqc_pipeline.sh`。

## 何时用

- 有一批 fastq（10x 或墨卓），要得到质控好的表达矩阵做下游分析
- 提到 cellranger / mobivision / 墨卓 / MobiNova / MobiCube / quantify / fmt_gtf / mkindex / velocyto / cellbender / scautoqc / 双胞过滤，且处于完整上游流程语境
- 要把 fastq "直接转成可分析的 h5ad"

## 流程总览

```
fastq ──比对定量──▶ matrix + BAM ──velocyto──▶ spliced/unspliced loom
                                                │
      format（02b）→ 三位一体校验（03b）→ cellbender ──▶ 去背景矩阵
                                                │
   06_scautoqc_pipeline.sh：bridge → gather → qc(original) → 阈值闸(06d)
     → scrublet + scDblFinder 双运行 → 簇级 BH(06e) → pool_all → finalize_qc
     → scautoqc_pooled_filtered.h5ad
全程完成 → 09b_check_delivery.sh（交付校验）→ 09c_cleanup.sh（归档+核验+清理）
```

四阶段独立、可断点续跑；幂等以输入指纹为准（`_runmeta/<sid>.<stage>.fingerprint.json`）。**执行模式：分阶段确认**——每阶段前与用户核对参数（样本表、路径、资源），产出后汇报验证再进入下一阶段；不要无人值守跑完全程。

## 首次部署（仅首次）

1. **环境**：`bash scripts/00_setup_envs.sh` — 建 `scautoqc`（scanpy + sctk + scrublet + celltypist）与 `velocyto` 两个 conda 环境，clone nf-scautoqc；cellbender 复用现有 GPU 环境（含 torch），缺失则提示。scDblFinder 复用任一装齐 `scDblFinder` + `schard` 的 R 环境，`run_config` 填 `RSCRIPT`/`SCANPY_PYTHON`。详见 [references/environment-setup.md](references/environment-setup.md)。
2. **参考**：`bash scripts/00_setup_references.sh` — 按 `REF_SOURCE` 分支处理：cellranger `TRANSCRIPTOME`（官方包直用，或 mkgtf 白名单 + mkref）或 mobivision `MOBI_INDEX`（**先 `fmt_gtf` 规范化，再 `mkindex`**，过滤按官方白名单）；缺参考给下载/构建指引；velocyto `VELO_REF`（genes.gtf + rmsk.gtf）。详见 [references/reference-setup.md](references/reference-setup.md)。

平台手册：[cellranger-platform.md](references/cellranger-platform.md)、[mobivision-platform.md](references/mobivision-platform.md)。

## 配置（每批数据必填）

1. `cp scripts/run_config.template.sh scripts/run_config.sh`，填：`PROJECT_NAME`、`SPECIES`、`SAMPLE_SHEET`、`HANDOFF_MD`、`ALIGNER`（+ `TRANSCRIPTOME` 或 `MOBIVISION_DIR`/`MOBI_INDEX`）、`VELO_REF`、`NF_SCAUTOQC_DIR`、三个环境名、`OUTPUT_ROOT`、`FASTQ_STAGING`、`RSCRIPT`/`SCANPY_PYTHON` 与资源参数。
   - **WSL 下 `OUTPUT_ROOT`/`FASTQ_STAGING` 禁止 `/mnt/*`**（处理期必须原生 ext4；`run_config` source 时硬校验拒跑）。
   - 可选：`DATASET_NAME`（同一 CR_ROOT 跑多数据集时必填）、`SC_ROOT`、`META_CSV`。
2. **写样本表**：仿 [assets/sample_sheet_template.csv](assets/sample_sheet_template.csv)，列 `sample_id,fastq_source_dir,fastq_prefix,bam_archive_dir,cell_or_nuclei`；`sample_id` 仅 `[A-Za-z0-9_-]`。格式见 [references/sample-sheet.md](references/sample-sheet.md)。
3. **收集十类元数据**（硬步骤）：`group/age/batch/datasets/sample/tech/gender/tissue/chemistry/species`，由用户交代或阻塞询问；不提供的记 `Unknown`，不猜测。整理为 csv 后用 `01c_write_metadata.py` 写入每样本 `outs/metadata.json`。

## 分阶段执行

### 阶段1 · 比对定量（每样本，严格串行）

- cellranger：`bash scripts/01_cellranger_count.sh <sid> <fastq_dir> <prefix> <bam_archive> <cell_or_nuclei>`
- mobivision：`bash scripts/01_mobivision_quantify.sh <sid> <fastq_dir> <delivery_dir>`

铁律：比对严格串行；样本完全隔离（上一样本清空后才复制下一个）；原始 fastq 只复制不移动；大文件串行复制、用 rsync 而非 cp。详见 [references/cellranger-platform.md](references/cellranger-platform.md) §4.1 与 [references/mobivision-platform.md](references/mobivision-platform.md) §3–5。

### 阶段2 · velocyto（可并发）

`bash scripts/02_velocyto.sh <sid> <samplefolder> <bam_archive>`

先校验 `VELO_REF` 与比对参考 GTF 的 SHA-256 一致，再跑 `velocyto run10x`（cellranger）或 `velocyto run`（mobivision，`-e <sid>` 保证 loom 命名）。**阶段边界**：样本表全部样本比对完成才放行 velocyto（02 开头闸口强制），除非用户显式声明逐样本接续（`VELO_PER_SAMPLE_HANDOFF=1`）。

### 阶段3 · cellbender（串行，GPU）

```bash
"$SCANPY_PYTHON" scripts/02b_format_matrices.py --cr-root <根> --samples <sids> --species <human|mouse>
"$SCANPY_PYTHON" scripts/03b_validate_triple_consistency.py --sid <sid> --cr-root <根>
bash scripts/03_cellbender.sh
```

`03_cellbender.sh` 内置先 format 再逐样本三位一体校验，任一失败即中止。校验通过后对每样本 `raw_feature_bc_matrix.h5` 跑 `cellbender remove-background --cuda`；GPU 独占，串行。

### 阶段4 · 质控（scautoqc）

```bash
bash scripts/06_scautoqc_pipeline.sh
```

一次跑完：准备物种 genes_list → 每样本输入校验 + features 对齐（06a）→ bridge → gather → qc（original）→ **阈值闸（06d，counts/genes/mito 硬阈值，纯阈值无 GMM）**→ scrublet → scDblFinder 双运行（默认 + 同型）→ **簇级 BH 检验（06e，sctk 同构）**→ pool_all → finalize_qc → 合并 scDblFinder 列 + 阈值闸列 → 清理临时文件。产出 `scautoqc_pooled_filtered.h5ad`。细节与阈值口径见 [references/stage-4-scautoqc.md](references/stage-4-scautoqc.md)。

### 阶段5 · 交付校验

```bash
bash scripts/09b_check_delivery.sh
```

逐项核验交付契约（标准四件/loom/cellbender 全套/metadata.json 十字段/scDblFinder 双运行四列 + 簇级四列/阈值闸产物/指纹记录/质控产出），缺失非零退出。缺失即交付失败——补齐复跑，不得绕过。

### 阶段6 · 归档与清理（先核验后删除）

```bash
bash scripts/09c_cleanup.sh --deliver-cr <交付根> [--deliver-sc <根>] [--apply]
```

先 rsync 归档补齐并对每对本地根→交付根做树级核验（逐文件存在且字节数一致），核验通过才允许删除本地工作根；默认 dry-run，`--apply` 才实际删除；删除目标带危险路径防护。

## 产出

- 每样本 `outs/`：标准四件（raw/filtered mtx 目录 + h5）、metrics、BAM（mobivision 另有报告）
- `<sid>/velocyto/<sid>.loom`；`<sid>/outs/cellbender/cellbender_matrix.h5` 全套
- `<sid>/outs/metadata.json`：十类信息 + `stages` 键（本阶段所用基因/版本快照）
- `scautoqc/`：三个 pooled h5ad、`qc_thresholds.csv`、`sample_passqc_df.csv`、`scdblfinder_scores.csv`、各目录产物
- `_runmeta/*.fingerprint.json`：断点续跑判据

## 交接下游

阶段4 末步向 `Document/HANDOFF.md` 段级 upsert `## upstream 阶段`：三个 pooled h5ad 路径/细胞数/layers/QC 列差异、sampleID、species/tech/chemistry、barcode 格式。字段约定见 [references/handoff-contract.md](references/handoff-contract.md)。

## 关键决策（勿改，除非有据）

1. **scanpy 1.10.4**：更早版本在 numpy 2 下 `np.float_` 缺失、布尔列绘图崩溃。
2. **qc_mode=original**：mito 循环 20/50/80（finalize 用 80），对含正常高 mito 细胞的实体组织友好；4 指标 GMM 阈值硬编码于 `qc.py`。
3. **gather_mode=starsolo + 桥接脚本**：nf-scautoqc 的 `read_velocyto` 只认 STARsolo mtx 目录，用 `04_bridge_to_starsolo.py` 适配（ENSG 对齐，规避重复符号）。
4. **禁止跨文件系统处理读写**：工作目录必须原生文件系统；跨文件系统只允许 rsync 进源数据、rsync 出交付。`run_config` 对 `/mnt/*` 硬校验。
5. **BAM 生命周期**：比对产出 → 归档源盘 → velocyto 用完即删中间与本地副本。
6. **scrublet + scDblFinder 双线 + 阈值闸**：两套独立算法都标注不过滤；阈值闸（06d）决定进入鉴定与交付的细胞；簇级 BH（06e）提供正交轴；交付层列名 `scrublet_bh_val`/`scrublet_stringent` 防歧义。
7. **三位一体校验**：比对矩阵 / loom / formatted h5 的基因与 barcode 集合一致、无重复；loom 为剪接分层计数，不与 counts 逐值比较；velocyto 必须用与比对相同的参考 GTF（SHA-256 校验）。
8. **数据集唯一性**：一次运行 = 一个数据集（一种 tissue×tech×chemistry）；sc 与 sn 不得混批（`cell_or_nuclei` 全一致硬校验）；单样本数据集在合并前止步并显式警告。
9. **交付校验 + 先核验后删除**：交付契约缺失即失败；清理必须先通过树级核验。
10. **幂等指纹**：skip 判据是输入指纹（脚本版本/参数/参考/输入文件/成员集合），不是文件存在；产物在但无指纹 = 按须重跑。
11. **源树只读**：产物只落 `OUTPUT_ROOT` 系，禁止回写源数据树。
12. **符号链接纪律**：会被原位改写的对象必须是运行树内真实副本，禁止写 symlink 指向的源对象。

## 清理

阶段内联清理（脚本自动）、收口归档清理（`09c_cleanup.sh`）、执行者侧临时件（日志重定向进项目 `logs/`，一次性脚本放临时目录并删除）三层都不得留残留。

## 排错

见 [references/troubleshooting.md](references/troubleshooting.md)：scanpy 布尔绘图崩、`np.float_`、flag_doublet 静默退出、cellranger `--id` 非法字符、9p 挂载慢、实体组织 QC 通过率、MobiVision 磁盘/静默损坏 BAM 等。

## 参考文档索引

| 文档 | 内容 |
|---|---|
| [pipeline-overview.md](references/pipeline-overview.md) | 四阶段数据流、目录布局、阶段契约 |
| [stage-1-3-upstream.md](references/stage-1-3-upstream.md) | 阶段1–3 细节：比对、velocyto、cellbender 与三位一体校验 |
| [stage-4-scautoqc.md](references/stage-4-scautoqc.md) | 质控阶段细节、阈值口径、双胞与簇级检验 |
| [cellranger-platform.md](references/cellranger-platform.md) | cellranger 平台手册 |
| [mobivision-platform.md](references/mobivision-platform.md) | MobiVision 平台手册（含调度铁律与磁盘预算） |
| [reference-setup.md](references/reference-setup.md) | 参考基因组构建与白名单 |
| [sample-sheet.md](references/sample-sheet.md) | 样本表格式与 sample_id 规则 |
| [environment-setup.md](references/environment-setup.md) | conda 环境、R 环境、依赖版本 |
| [handoff-contract.md](references/handoff-contract.md) | HANDOFF 段落字段契约 |
| [delivery-verification.md](references/delivery-verification.md) | 交付对象只读核验：obs 契约、压缩重建、断言 |
| [troubleshooting.md](references/troubleshooting.md) | 常见故障与排错 |
| [CONTRIBUTING.md](CONTRIBUTING.md) | 仓库约定：零项目残留、编码/换行、跨平台、数据流原则、提交前检查 |

## 给执行者的提醒

- 先读 [references/pipeline-overview.md](references/pipeline-overview.md)；`ALIGNER=mobivision` 时先读平台手册
- 每阶段产出都要验证（文件存在、shape 合理、无错误日志）才进下一步；cellbender 前必须过三位一体校验
- 线程/内存默认交给软件决定，不显式限制
- 长任务写成脚本文件后 `setsid nohup … &` 脱离启动，定期轮询日志；跨调用持久物放常驻目录
- 规范未覆盖的情况（非 10x、双端 index、多 lane fastq）先查官方文档再决策，不臆断；阻塞时要求用户介入
