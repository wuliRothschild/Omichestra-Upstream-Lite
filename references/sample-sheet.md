# 样本表格式

CSV 格式，驱动所有阶段。首行表头，后续每行一个样本。

## 列定义

| 列 | 必填 | 含义 |
|----|------|------|
| `sample_id` | 是 | 全程标识，作 cellranger `--id`/`--sample`、velocyto/cellbender/scautoqc 目录名。仅 `[A-Za-z0-9_-]`（cellranger 由 martian 层硬校验：数字/字母/连字符/下划线，无首字符限制） |
| `fastq_source_dir` | 是 | 该样本 fastq 所在目录（含 `<prefix>_S1_L001_R{1,2}_001.fastq.gz` 与 `.md5`） |
| `fastq_prefix` | 是 | fastq 文件名前缀（`<prefix>_S1_L001_R1_001.fastq.gz` 的 `<prefix>`） |
| `bam_archive_dir` | 否 | cellranger BAM 归档目录（通常=fastq_source_dir）。留空则不归档（velocyto 需 BAM，一般要填） |
| `cell_or_nuclei` | 是 | `cell` 或 `nuclei`，**必填**（默认 cell）。scautoqc 据此切换 4 指标 GMM 阈值：cell→`n_counts≥1000, mito 0.1–max, spliced 50–97.5`；nuclei→`n_counts≥300, mito 0.1–max, spliced 0.01–100`。解离组织用 cell，单核测序用 nuclei。填错会导致阈值严重偏差。**全表必须同值**：单细胞（sc）与单细胞核（sn）不得归入同一数据集（一个数据集只允许一种 tissue×tech×chemistry 组合）；06 开头硬校验，混杂即拒跑——sc/sn 混批请拆成两张表分批运行 |

## 注释与空行

`#` 开头的行与空行被脚本跳过（awk 过滤），可放说明。

## 示例

```csv
sample_id,fastq_source_dir,fastq_prefix,bam_archive_dir,cell_or_nuclei
sample1,/path/to/fastq/sample1,S1,/path/to/fastq/sample1,cell
sample2,/path/to/fastq/sample2,S2,/path/to/fastq/sample2,cell
sample3,/path/to/fastq/sample3,S3,/path/to/fastq/sample3,cell
```

（sn 数据集另用一张全 `nuclei` 的表——sc/sn 不得同表混跑，见上表 `cell_or_nuclei` 行）

## 命名建议

`sample_id` 建议编码实验设计（组织/处理/重复），便于下游：
- `Fe` = iron / `De` = deficiency / `Ov` = overload
- `Sp` = spleen / `Li` = liver
- `E1`/`C1` = experiment replicate 1 / control replicate 1

避免用 illumina 原始 sample name（A0/N1 等无语义）直接作 sample_id。

## sample_id 重命名的意义

脚本会把 fastq 文件从 `<prefix>_S1_L001_R1_001.fastq.gz` 重命名为 `<sample_id>_S1_L001_R1_001.fastq.gz`，使 cellranger `--sample=<sample_id>` 命中。这样：
- 多个样本即使原始 prefix 相同（如两组都有 N1）也不会混淆（各自在独立暂存目录）
- 下游所有产出以 sample_id 命名，可追溯

## 注意

- 多 lane（L001+L002）或多套 fastq：当前脚本只处理单 lane L001 的 R1/R2 各一个。若有多 lane，需把所有 lane 的 fastq 放同一目录（cellranger 会自动合并同名不同 lane），或调整 01 脚本的 rsync/重命名逻辑
- 双 index（如 cellranger multi 的抗体/CRISPR）：本流程只处理 gene expression，多模态需扩展

## MobiVision 分支（`ALIGNER=mobivision`）

- `fastq_prefix` 可留空：墨卓命名（`[SampleID]_[时间戳].R1/R2.fastq.gz`）由 mobivision 自动识别，`01_mobivision_quantify.sh` 不做前缀重命名
- `fastq_source_dir` 可能含多级子目录（墨卓按测序批次分 `<S编号>_<hm|rg|yg|zf>/` 子目录），脚本会平铺复制到暂存目录
- `bam_archive_dir` 在 mobivision 分支语义为**交付根目录**（`01` 脚本交付 `{sid}/outs/`，`02` 脚本从该目录 rsync 回 BAM 供 velocyto）
- 多批次 fastq 由 mobivision quantify 自动 merge，无需预合并；样本目录内 fastq 文件名必须唯一（平铺冲突会报错）
- 完整平台手册见 [mobivision-platform.md](mobivision-platform.md)
