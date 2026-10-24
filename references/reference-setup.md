# 参考基因组构建与参考策略

本文件是跨平台（cellranger/mobivision）的**参考构建权威文档**：REF_SOURCE 四分支路由、10X 官方过滤白名单（**唯一源**）、官方预构建参考包（refdata-gex）的下载/解压/结构、velocyto VELO_REF 与参考版本一致性。平台专用命令细节见 [cellranger-platform.md](cellranger-platform.md) 与 [mobivision-platform.md](mobivision-platform.md)。

## 参考策略（REF_SOURCE 四分支）

`run_config.sh` 的 `REF_SOURCE` 决定参考来源（默认 `refdata`）：

| REF_SOURCE | 行为 | 何时用 |
|---|---|---|
| `refdata`（默认） | 用 `$REFDATA_DIR/refdata-gex-{GRCh38\|GRCm39}-$REFDATA_VERSION/` 解压目录。cellranger 直接作 `TRANSCRIPTOME`；mobivision 从该目录 `fasta/genome.fa` + `genes/genes.gtf.gz` 跑 `mkindex` 重建专用索引 | 默认；本地有官方包 |
| `download` | 本地无官方包 → 按 `REFDATA_URL`（默认按物种生成官方直链）下载 tar.gz 并解压到 `REFDATA_DIR` | 本地无包、需联网获取 |
| `build` | 从头构建，用**官方过滤白名单**（cellranger `mkgtf` 或 mobivision `fmt_gtf -t`）后 mkref/mkindex | 本地无包且需自建 |
| `custom` | 用户显式指定自定义方案（如 mobivision `fmt_gtf` 不传 `-t` 全量保留、或任意自定义 GTF），脚本仅校验 `TRANSCRIPTOME`/`MOBI_INDEX` 完整性 | 用户显式指定 |

- `REFDATA_VERSION` 默认 `2024-A`，可改（如 `2020-A`，GENCODE v32 对应，见下）。
- `SPECIES` → refdata 名映射：`human`→`GRCh38`、`mouse`→`GRCm39`。
- 所有分支均"检查 + 指引"，不自动下载/构建（10G+ 下载与数小时 mkindex 由用户确认后手动执行）。

## 10X 官方过滤白名单（唯一事实源）

> 其他文档（cellranger-platform.md / mobivision-platform.md）只写"按此白名单"指针，**本表仅在此一处**。来源：10x 官方文档 `cr-3p-references`（mkgtf 可接受 biotype 列表逐字）+ 官方构建笔记 `cr-ref-build-steps`（BIOTYPE_PATTERN）。

官方参考 = 保留 **protein_coding + lncRNA + 全部 IG/TR 基因与假基因**，剔除其余全部 biotype（miRNA/snRNA/snoRNA/rRNA/常规假基因/反义等）。白名单 17 项：

```
protein_coding, lncRNA, antisense,
IG_C_gene, IG_D_gene, IG_J_gene, IG_LV_gene, IG_V_gene,
IG_V_pseudogene, IG_J_pseudogene, IG_C_pseudogene,
TR_C_gene, TR_D_gene, TR_J_gene, TR_V_gene,
TR_V_pseudogene, TR_J_pseudogene
```

- 注：GENCODE v31/M22（Ensembl 97）起 `lincRNA`/`antisense` 并入 `lncRNA`，故 v32+ 输入下 17 项与官方构建脚本的 `BIOTYPE_PATTERN`（16 项）等价。
- 2024-A 官方包额外纳入 `protein_coding_LoF`（多态假基因，发布说明原文："Polymorphic pseudogenes were included by adding protein_coding_LoF to the list of accepted biotypes"）。

### 两物种实测基因数（标注：实测，非官方公布）

官方文档页不公布基因数（5 个官方页面 grep 0 命中），以下为对官方包 `genes/genes.gtf.gz` 的**实测统计**：

| 官方包 | 基因总数 | 构成 |
|---|---|---|
| GRCh38-2024-A | **38,606** | protein_coding 19,408 + lncRNA 18,551 + IG/TR 647 |
| GRCm39-2024-A | **33,696** | protein_coding 22,134 + lncRNA 11,157 + IG/TR 405 |

## 官方预构建参考包（REF_SOURCE=refdata / download）

官方构建笔记原文："10x Genomics offers pre-built Cell Ranger reference packages from the downloads page. For purposes of reproducibility, the exact build steps are provided here."

| 官方包 | 注释 | 官方下载直链 |
|---|---|---|
| `refdata-gex-GRCh38-2024-A.tar.gz`（~10.7GB） | GENCODE v44 / Ensembl 110 | `https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz` |
| `refdata-gex-GRCm39-2024-A.tar.gz` | GENCODE vM33 / Ensembl 110 | `https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz` |

- 版本对应（官方发布说明 `cr-reference-release-notes`）：2020-A = "Transcriptome annotations updated from Ensembl 93 to GENCODE v32 (human) and vM23 (mouse), which are equivalent to Ensembl 98… chromosome names now follow the GENCODE/UCSC convention (e.g., chr1 and chrM)"；2024-A = "updated from GENCODE v32 to GENCODE v44… PAR masked"。
- 解压：`tar -xzf refdata-gex-GRCh38-2024-A.tar.gz -C $REFDATA_DIR` → `$REFDATA_DIR/refdata-gex-GRCh38-2024-A/`（含 `fasta/genome.fa` + `genes/genes.gtf.gz` + `star/` + `reference.json`）。
- 结构解读（reference.json 字段）见 [cellranger-platform.md §3.3](cellranger-platform.md)（唯一处）。

### 两平台消费方式

- **cellranger**：解压目录直接作 `TRANSCRIPTOME`（结构同构，无需重建）。
- **mobivision**：从解压目录取 `fasta/genome.fa` + `genes/genes.gtf.gz`（已是官方过滤集合）作 `mobivision mkindex` 输入，重建 mobivision 专用索引（命令见 [mobivision-platform.md §2.3](mobivision-platform.md)）。**不直接复用 cellranger 的 `star/` 索引**——两软件 STAR 参数/校验字段（reference.json 的 mkref_version）未经实证互通，重建最稳。

## 文献表述纠正

文献/数据仓库常见的"Transcriptome reference **exactly matching Cell Ranger 2020-A**"（如 Oliver 2024 Methods）**不是 10X 官方文档原文**——经 5 个官方页面全文 grep 均 0 命中。官方对应的表述是发布说明的"Transcriptome annotations updated from … to GENCODE v32…"。写文档/论文时引用官方原句；skill 默认用 2024-A（GENCODE v44）。

## 从头构建（REF_SOURCE=build）

按官方过滤白名单从 GENCODE/Ensembl 注释过滤后构建：

- **cellranger**：`cellranger mkgtf` 17 项白名单 → `cellranger mkref`（完整命令与约束见 [cellranger-platform.md §2](cellranger-platform.md)）。
- **mobivision**：`mobivision fmt_gtf -t` 同白名单 → `mobivision mkindex`（见 [mobivision-platform.md §2](mobivision-platform.md)）。

## velocyto VELO_REF

velocyto `run10x`/`run` 需两个文件：
- `genes.gtf`：含 exon 的 GTF（velocyto 据此生成 intron 模型）。**必须与比对阶段所用 GTF 逐字节一致**（02_velocyto.sh 开头做 SHA-256 校验）：
  - `ALIGNER=cellranger`：用 `TRANSCRIPTOME/genes/genes.gtf`（官方包内 `genes/genes.gtf.gz` 解压）——与比对参考同构建，天然一致。
  - `ALIGNER=mobivision`：用 `MOBI_INDEX/genes/genes.gtf`——mkindex 输出内的 GTF 已经过 `fmt_gtf` 规范化，与官方包原始 GTF 并非逐字节相同；**不要**用官方包原始 `genes.gtf.gz` 解压充当，否则 SHA-256 校验中止。
- `rmsk.gtf`：RepeatMasker 重复区注释（`-m` 参数，mask 重复区避免假性 unspliced）。

### rmsk.gtf 获取

- UCSC Table Browser：group=Repeats，track=RepeatMasker，output=GTF，下载对应基因组（mm39/hg38）
- 或用 velocyto 仓库的 `get_rmsk.py` 从 RepeatMasker `.out` 文件转换

两个文件放同一目录，该目录填入 `run_config.sh` 的 `VELO_REF`。

## 参考版本一致性（三位一体前提）

- cellranger `TRANSCRIPTOME`（或 mobivision `MOBI_INDEX`）须与 velocyto `VELO_REF` **同一参考 GTF**——`VELO_REF/genes.gtf` 与比对参考的 `genes/genes.gtf` 必须逐字节一致，否则 loom 基因注释错位与剪接量漂移。02_velocyto.sh 开头校验两者基因 GTF 的 **SHA-256**，不一致中止；无法取得比对参考 GTF 时同样中止（不以基因数相等替代同源性校验）。
- REF_SOURCE=refdata 时：cellranger 分支比对参考与 `VELO_REF/genes.gtf` 均来自同一官方包（`genes/genes.gtf.gz` 解压），天然一致；mobivision 分支须取 `MOBI_INDEX/genes/genes.gtf`（fmt_gtf 规范化产物，见上）。
- **三位一体**：比对矩阵（filtered mtx）/ velocyto loom / formatted h5 必须坐标系严格对齐、任何一方未经自定义改动（去重/补零/符号对齐）——mtx 与 h5 为同源转换：基因行（ENSG）顺序、barcode 列序、表达值逐值硬等价；loom 与二者：基因（ENSG）集合与 barcode（归一后）集合完全一致、无重复、shape 一致（行列顺序可与比对矩阵不同——velocyto 按自身 GTF 解析/BAM 顺序输出，由下游 bridge 按 ENSG+barcode 严格重排）；loom 为剪接分层计数（spliced/unspliced/ambiguous，主层空占位），`s+u+a ≠ GeneFull counts`（velocyto 计数模型与比对器不同），不做逐值比较，但层结构必须完整（spliced/unspliced 存在、非空）；loom 若另含显式 `GeneFull` 层，则按坐标映射重排后与 mtx/h5 逐值比较。cellbender 前运行 `scripts/03b_validate_triple_consistency.py` 校验（不一致中止）。
- 基因符号命名规则：小鼠首字母大写（Mki67），人类全大写（MKI67）——genes_list 物种适配脚本 05 据此处理。
