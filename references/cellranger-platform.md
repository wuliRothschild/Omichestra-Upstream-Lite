# Cell Ranger 9.0 平台手册（10x Genomics）

> 本文档是 Omichestra-upstream 的 Cell Ranger 分支完整手册：参考基因组构建（cellranger mkref / mkgtf）、官方预构建参考（refdata-gex）的使用、cellranger count 比对最佳实践。内容来自 10x Genomics 官方文档（www.10xgenomics.com）逐字引用 + cellranger-9.0.1 安装目录源码/二进制实证。
> 适用前提：`run_config.sh` 中 `ALIGNER=cellranger`。参考策略（REF_SOURCE 四分支）与官方过滤白名单见 [reference-setup.md](reference-setup.md)；mobivision 分支见 [mobivision-platform.md](mobivision-platform.md)。

## 目录

1. [平台总览与软件环境](#1-平台总览与软件环境)
2. [参考基因组构建（mkref + mkgtf）](#2-参考基因组构建mkref--mkgtf)
3. [官方预构建参考 refdata-gex 的使用](#3-官方预构建参考-refdata-gex-的使用)
4. [cellranger count 比对最佳实践](#4-cellranger-count-比对最佳实践)
5. [交付规范](#5-交付规范)
6. [故障排除](#6-故障排除)
7. [下游衔接（阶段2-4）](#7-下游衔接阶段2-4)

---

## 1. 平台总览与软件环境

- **软件**：cellranger 9.0.1（独立二进制，非 conda；从 10x Genomics 下载 tar.gz 解压，`export PATH=<cellranger-9.0.1>`）。`run_config.sh` 的 `CELLRANGER_DIR` 指向安装目录（含 `cellranger` 可执行），留空则用 PATH。
- **子命令**：count（比对定量）、mkref（构建参考）、mkgtf（GTF 过滤）、vdj/atac/multi 等。
- **参考要求**：`cellranger count --transcriptome` 指向 mkref 或官方 refdata 产出的目录（含 `reference.json`/`fasta/`/`genes/`/`star/`）。
- **run_config.sh 相关变量**：`CELLRANGER_DIR`（安装根）、`TRANSCRIPTOME`（参考目录，REF_SOURCE=refdata 时自动定位官方包；build/custom 时填构建输出）。

## 2. 参考基因组构建（mkref + mkgtf）

> 10x 官方文档原文："To create custom references, use the `cellranger mkref` command, passing it one or more matching sets of FASTA and GTF files. This utility copies your FASTA and GTF, indexes these in several formats, and outputs a folder with the name you pass to `--genome`. Input GTF files are typically filtered with `mkgtf` prior to `mkref`."

### 2.1 mkref 完整用法（官方文档 cr-3p-references + 9.0.1 源码 mkref.rs）

```
cellranger mkref --genome=<输出目录名> --fasta=<genome.fa> --genes=<genes.gtf> [--memgb=16] [--nthreads=1] [--ref-version=<串>]
```

| 参数 | 说明（官方逐字/源码注释） |
|---|---|
| `--genome` | 必需。输出文件夹名 `[a-zA-Z0-9_-]+`；多基因组重复指定，输出 `<name1>_and_<name2>` |
| `--fasta` | 必需。基因组 FASTA 路径 |
| `--genes` | 必需。基因注释 GTF 路径 |
| `--memgb` | 可选。STAR 索引最大内存 GB，**默认 16**（须大于参考 FASTA 的 GB 数） |
| `--nthreads` | 可选。STAR 索引线程数，**默认 1** |
| `--ref-version` | 可选。写入 reference.json 的版本串 |

**约束（官方与源码实证）**：
- FASTA **不可 gzip**（必须纯文本 `.fa`/`.fasta`；官方 refdata 内为 `genome.fa`）。
- GTF 须含 `exon` 行、`transcript_id`/`gene_id` 属性（reference_builder.py 校验：9 列、start<end、strand、contig、跨染色体转录本去重、必须有 exon）。
- **mkref 不做 biotype 过滤**——过滤由 `mkgtf` 单独完成（见 2.2）。
- 系统需求（官方）：人 3Gb FASTA 建索引约 8 核时、需 32 GB 内存。
- `reference.json` 的 `mkref_version` 与 cellranger 版本：**无硬校验**——2024-A 官方包标注 `mkref_version: "8.0.0"` 可被 cellranger 9.0.1 正常使用（9.0.1 与 8.0.1 mkref 逻辑一致）；但 10x 官方建议参考与 cellranger 版本匹配。

### 2.2 mkgtf —— 官方 GTF 过滤工具（过滤的唯一入口）

官方文档原文："Cell Ranger provides `mkgtf`, a simple utility to filter genes based on their key-value pairs in the GTF attribute column. The command syntax requires input and output GTF file names and `--attribute` values specifying gene biotypes to filter from the GTF file:"

```
cellranger mkgtf input.gtf output.gtf --attribute=key:allowable_value
```

官方构建 GRCh38 参考包时实际执行的过滤命令（逐字，17 项白名单；**这是官方参考基因集合的唯一权威来源**）：

```
cellranger mkgtf Homo_sapiens.GRCh38.ensembl.gtf Homo_sapiens.GRCh38.ensembl.filtered.gtf \
  --attribute=gene_biotype:protein_coding \
  --attribute=gene_biotype:lncRNA \
  --attribute=gene_biotype:antisense \
  --attribute=gene_biotype:IG_LV_gene \
  --attribute=gene_biotype:IG_V_gene \
  --attribute=gene_biotype:IG_V_pseudogene \
  --attribute=gene_biotype:IG_D_gene \
  --attribute=gene_biotype:IG_J_gene \
  --attribute=gene_biotype:IG_J_pseudogene \
  --attribute=gene_biotype:IG_C_gene \
  --attribute=gene_biotype:IG_C_pseudogene \
  --attribute=gene_biotype:TR_V_gene \
  --attribute=gene_biotype:TR_V_pseudogene \
  --attribute=gene_biotype:TR_D_gene \
  --attribute=gene_biotype:TR_J_gene \
  --attribute=gene_biotype:TR_J_pseudogene \
  --attribute=gene_biotype:TR_C_gene
```

- 过滤逻辑（源码 `GtfBuilder.build_gtf()`）：GTF 属性中 `key` 存在且值不在白名单即剔除该行。
- 官方文档注：GENCODE v31/M22（Ensembl 97）起 `lincRNA`/`antisense` 并入 `lncRNA` biotype——故 v32+ 输入下 17 项与 10x 构建脚本的 `BIOTYPE_PATTERN`（16 项）等价。
- 白名单值表与两物种实测基因数见 [reference-setup.md](reference-setup.md)（唯一源，此处不重复）。

### 2.3 从头构建流程（REF_SOURCE=build）

1. 下载 GENCODE primary_assembly 注释 + primary_assembly FASTA（官方构建笔记用 `ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/` 的 `gencode.v44.primary_assembly.annotation.gtf.gz`）。
2. `cellranger mkgtf` 按 2.2 白名单过滤 → `{genome}.filtered.gtf`。
3. `cellranger mkref --genome=... --fasta=... --genes=...filtered.gtf --ref-version=...`（不主动传线程/内存参数，由软件默认与系统级资源配置决定；2024-A 官方构建记录使用 `--nthreads=16`，仅作事实记录，非默认命令参数）。
4. 产出目录作 `TRANSCRIPTOME`。

## 3. 官方预构建参考 refdata-gex 的使用

### 3.1 官方包与下载（REF_SOURCE=refdata / download）

官方构建笔记原文："10x Genomics offers pre-built Cell Ranger reference packages from the downloads page. For purposes of reproducibility, the exact build steps are provided here."

| 官方包 | 注释来源 | 下载直链（cf.10xgenomics.com） |
|---|---|---|
| `refdata-gex-GRCh38-2024-A.tar.gz` | GENCODE v44 / Ensembl 110（human） | `https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz` |
| `refdata-gex-GRCm39-2024-A.tar.gz` | GENCODE vM33 / Ensembl 110（mouse） | `https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz` |

- 版本对应（官方发布说明）：2020-A = GENCODE v32 (human)/vM23 (mouse) = Ensembl 98；2024-A = GENCODE v44 (human)/vM33 (mouse) = Ensembl 110，2024-A 增加 `protein_coding_LoF`（多态假基因）并做 PAR masking。
- 解压到 `REFDATA_DIR`（run_config）：`tar -xzf refdata-gex-GRCh38-2024-A.tar.gz -C $REFDATA_DIR` → 目录 `$REFDATA_DIR/refdata-gex-GRCh38-2024-A/`。

### 3.2 refdata 作为 TRANSCRIPTOME（cellranger 直用）

官方包目录结构与 mkref 输出同构（`fasta/genome.fa` + `genes/genes.gtf.gz` + `star/` + `reference.json`），**cellranger 直接以解压目录作 `TRANSCRIPTOME`**，无需重建。

### 3.3 reference.json 字段解读（唯一事实源处）

2024-A 官方包实测 `reference.json`：

```json
{
    "fasta_hash": "...",
    "genomes": ["GRCh38"],
    "gtf_hash.gz": "...",
    "input_fasta_files": ["Homo_sapiens.GRCh38.dna.primary_assembly.fa.modified"],
    "input_gtf_files": ["gencode.v44.primary_assembly.annotation.gtf.filtered"],
    "mem_gb": 16,
    "mkref_version": "8.0.0",
    "threads": 2,
    "version": "2024-A"
}
```

| 字段 | 含义 |
|---|---|
| `input_gtf_files` | 构建输入的 GTF（**`.filtered` = 已按 2.2 白名单过滤**）；`version` = 参考版本（2024-A） |
| `mkref_version` | 构建所用 cellranger 版本；9.0.1 接受 8.0.0 构建的参考（无硬校验，见 2.1） |
| `mem_gb`/`threads` | 构建时资源参数（记录用，运行不依赖） |

## 4. cellranger count 比对最佳实践

```
cellranger count --id=<sid> --transcriptome=$TRANSCRIPTOME --fastqs=<fastq_dir> \
  --sample=<sid> --expect-cells=<N> \
  --nosecondary --create-bam=true
```

- **fastq 先 rsync 到本地 NVMe**（直读 9p/网络挂载 I/O 极差），md5 校验后再 count；源数据只复制，不移动或剪切。复制大文件时保持单文件串行，避免机械硬盘并发 I/O。
- `--nosecondary`（跳过降维/聚类，加速）、`--create-bam=true`（velocyto 需要 BAM）。
- 产出验证：`outs/filtered_feature_bc_matrix.h5` + `raw_feature_bc_matrix.h5` + `possorted_genome_bam.bam` + `web_summary.html`（无错误）、`Estimated Number of Cells` 合理。
- BAM 生命周期：产出 → 归档源盘 → velocyto 用完即删中间 `cellsorted_*` 与本地副本。

### 4.1 调度铁律（与 mobivision 分支同一标准）

> 平台差异只在中间产物规模（cellranger 的 `SC_RNA_COUNTER_CS/` 为 GB 级，远小于 mobivision 的 convert_fastq），调度约束相同。

1. **count 严格串行**：单个 count 默认占满机器线程与数十 GB 内存；并行多个 count 互相争抢 CPU/内存/磁盘 I/O，反而更慢且增加失败面。一个样本完全结束（脚本第 7 步白名单清理完成）后才启动下一个。
2. **样本完全隔离**：上一个样本的 staging 与输出目录清空后，才能复制下一个样本的 fastq；01 脚本开头有残留检查（staging/样本输出目录非空即中止），防止中断重跑或并行造成的混样。
3. **velocyto 间隔发车**：所有样本 count 完成后才进入 velocyto 阶段；并发上限与 15 分钟发车间隔见 SKILL.md 阶段2（`MAX_JOBS`/`INTERVAL`，外层编排）。
4. **禁止跨文件系统处理读写**：staging 与输出根必须在 WSL 原生 ext4，禁止 `/mnt/*`（drvfs/9p 吞吐降 1–2 个数量级）；跨文件系统仅允许 rsync 复制源 fastq 进本地、rsync 交付/BAM 归档回挂载盘。run_config source 时硬校验，违规拒跑（见 SKILL.md 关键原生文件系统约束）。

## 5. 交付规范

- 样本目录 `<sid>/outs/`（cellranger 原生即 `outs/`）：`filtered/raw_feature_bc_matrix{,.h5}`、`metrics_summary.csv`、`web_summary.html`、`possorted_genome_bam.bam(+.bai)`。
- 交付验证：h5/mtx 可读、web_summary 无错误、BAM 非空。
- 大文件传输用 `rsync -a --append-verify` 而非 cp（drvfs/9p 挂载写大文件会"无法分配内存"）。

## 6. 故障排除

- **mkref 报 GTF 无 exon / transcript_id 缺失**：GTF 不合规，换 GENCODE primary_assembly 注释或先 mkgtf。
- **mkref 内存不足**：`--memgb` 须大于参考 FASTA GB 数（官方默认 16 对 3Gb 人基因组需调大，官方建议 32）。
- **STAR 索引太慢**：`--nthreads` 提高（2024-A 官方构建用 16）。
- **参考版本不匹配 cellranger**：10x 建议匹配，但 8.0.0 构建的 2024-A 参考在 9.0.1 实测可用（无硬校验）。

## 7. 下游衔接（阶段2-4）

cellranger count 产出的 `filtered_feature_bc_matrix.h5` + BAM 是后续标准输入：
1. 阶段2 velocyto：`velocyto run10x -m <rmsk.gtf> <cellranger_samplefolder> <genes.gtf>`（线程使用 velocyto 默认值；`VELO_REF` 见 reference-setup.md）。
2. 阶段3 cellbender：`raw_feature_bc_matrix.h5`。
3. 阶段4 质控：`06_scautoqc_pipeline.sh`，文件契约见 SKILL.md。
