# MobiVision V4.1 平台手册（墨卓 MobiNova）

> 本文档是 Omichestra-upstream 的 MobiVision 分支完整手册：参考基因组构建（fmt_gtf / mkindex）、原始测序数据比对（quantify）、严格流式调度与交付规范。
> 适用前提：`run_config.sh` 中 `ALIGNER=mobivision`。参考策略（REF_SOURCE 四分支）与官方过滤白名单见 [reference-setup.md](reference-setup.md)；cellranger 分支见 [cellranger-platform.md](cellranger-platform.md)。

## 目录

1. [平台总览与软件环境](#1-平台总览与软件环境)
2. [参考基因组构建（fmt_gtf + mkindex）](#2-参考基因组构建fmt_gtf--mkindex)
3. [quantify 全流程与中间产物生命周期](#3-quantify-全流程与中间产物生命周期)
4. [严格流式调度铁律](#4-严格流式调度铁律)
5. [资源与时长参考](#5-资源与时长参考)
6. [交付规范](#6-交付规范)
7. [故障排除](#7-故障排除)
8. [下游衔接（阶段3-4）](#8-下游衔接阶段3-4)

---

## 1. 平台总览与软件环境

- **软件**：`mobivision-v4.1`（conda 环境分发，自带 `bin/` 下 STAR 2.7.10b、samtools、seqkit、cutadapt 等全部依赖）
- **激活**：每个新 shell 必须 `source <mobivision-v4.1>/source.sh`（设置 PATH + ulimit）
- **子命令**：quantify（比对定量）、fmt_gtf（注释过滤）、mkindex（构建索引）、rename/mtx2csv/re_call_cell/integrate、vdj/cite/chip/atac 等
- **输出命名**：样本产物为 `{sample_id}_outs/` 目录（交付时按本 skill 规范改名为 `outs/`）
- **阶段日志**：`logs/mobivision.log` 逐阶段记录；`_flagdone` 表示成功结束
- **run_config.sh 相关变量**：`MOBIVISION_DIR`（安装根，含 source.sh）、`MOBI_INDEX`（mkindex 输出）、`MOBI_MEM_GB`（mkindex `-m` 覆盖位，默认空=不显式传）

## 2. 参考基因组构建（fmt_gtf + mkindex）

### 2.1 fmt_gtf —— 注释 GTF 过滤

```
mobivision fmt_gtf -i <input.gtf> -o <out.gtf> [-t <biotype>]... [--mt <mt基因列表>]
```

- **`fmt_gtf` 是 `mkindex` 前的必经规范化步骤**：即使需要保留全部基因、实质上不按 biotype 过滤，也必须运行 `fmt_gtf`，但不传任何 `-t` 参数；不得直接把未经 `fmt_gtf` 规范化的原始 GTF 交给 `mkindex`。
- **`-t` 匹配的是 biotype 值**（`gene_type`/`transcript_type`），**不是 feature 类型**。默认（REF_SOURCE=build）按 **10X 官方过滤白名单**过滤（17 项：protein_coding + lncRNA + 全部 IG/TR 基因与假基因，见 [reference-setup.md](reference-setup.md) 白名单唯一源）。
- **无 `-t` 的全量模式仍会规范化 GTF**：它保留输入中的全部 gene、transcript、exon、CDS、UTR 等注释行，但会按 MobiVision 要求规范化属性和格式；因此输出文件不应被视为输入文件的字节级副本。 
- **说明书官方样例**（MobiVision v4.1 使用说明书 §1.2 fmt_gtf，逐字，17 个 `-t` 含 `antisense`，与 10X 官方 mkgtf 白名单等价）：
  ```
  $ mobivision fmt_gtf -i /media/db/gencode.v38.primary_assembly.annotation.gtf.gz \
    -o /media/db_out/gencode.v38.primary_assembly.annotation.filtered.gtf \
    -t protein_coding \
    -t lncRNA \
    -t antisense \
    -t IG_LV_gene -t IG_V_gene -t IG_V_pseudogene -t IG_D_gene -t IG_J_gene \
    -t IG_J_pseudogene -t IG_C_gene -t IG_C_pseudogene \
    -t TR_V_gene -t TR_V_pseudogene -t TR_D_gene -t TR_J_gene \
    -t TR_J_pseudogene -t TR_C_gene
  ```
- **全量保留所有基因 = 不传 `-t`**（REF_SOURCE=custom 选项）：源码 `bool_biotype` 在过滤列表为空时全部放行，gene 行全量保留，exon/CDS/transcript 等行数与输入完全一致（零丢失）。默认走官方白名单（对齐 10X 官方参考基因集合）；需全量时才显式 custom。
- **注意 `-t gene` 是错误写法**：`gene` 不是任何 biotype 值，会匹配不到任何基因导致输出为空。
- **数据源推荐**（说明书 §1.2 与（五）参考序列信息）：人 = GRCh38 的 **GENCODE v38**（MobiNova 平台提供版本，备选 Ensembl v104）；鼠 = GRCm39 的 **GENCODE vM27**（备选 Ensembl v104）。输入支持 `.gtf` / `.gtf.gz`。
- **Ensembl GTF 兼容**：自动把 `gene_biotype`/`transcript_biotype` 替换为 `gene_type`/`transcript_type` 再匹配；`gene_id "ENSG...17"` 拆分为 `gene_id` + `gene_version` 两个属性（transcript_id/exon_id 同理）。
- **硬编码排除**：带 `tag "readthrough_transcript"` 或 `tag "PAR"` 的转录本永远被排除（Ensembl 人类注释实测无基因因此丢失）。
- **`--mt`**：给列表命中的基因加 `mt-` 前缀（匹配 gene_id 值或 gene_name 值；若列表是基因符号如 `MT-ND1`，走 gene_name 改名分支）。
- **已知无害现象**：对带显式 `gene_version` 属性的 Ensembl GTF，输出会同时含 `gene_version "01"`（拆分默认值）与原值，STAR 解析不受影响。
- 输出 GTF 建议命名为 `{genome}.filtered.gtf` 保留备查。

### 2.2 mkindex —— 构建 MobiVision 兼容索引

```
mobivision mkindex -n <GENOME> -f <genome.fa> -g <filtered.gtf> [-r <版本串>] [-m <内存GB>]
```

- **参数表**（MobiVision v4.1 使用说明书 §1.3 mkindex）：`-n` 基因组名称（`[a-zA-Z0-9_]`，多基因组重复指定）、`-f` 基因组 fasta（支持 `.fasta/.fa.gz/.fna.gz`）、`-g` 注释 GTF（`.gtf/.gtf.gz`，多基因组按顺序对应）、`-r` 参考版本号（可选，如 `"ref_v1.0"`）、`-m` **STAR 构建参考最大内存 GB**（可选）。**无 `-t` 线程参数**（固定单线程）。
- **说明书样例**（逐字，单物种）：
  ```
  $ mobivision mkindex -n GRCh38_gencode \
    -f /media/db_out/GRCh38.primary_assembly.genome.fa.gz \
    -g /media/db_out/gencode.v38.primary_assembly.annotation.filtered.gtf
  ```
  多物种：`-n/-f/-g` 按顺序重复指定，输出文件夹 `<genome1>_and_<genome2>`。
- **输出**：运行命令所在目录下 `<GENOME>/` 文件夹，结构（与 cellranger refdata 同构）：
  ```
  <GENOME>/
  ├── fasta/genome.fa (+ .fai via samtools faidx)
  ├── genes/genes.gtf
  ├── star/          # STAR genomeGenerate 索引（Genome/SA/SAindex/sjdbList.*）
  └── reference.json # fasta/gtf 的 SHA1、输入文件、mem_gb、mkref_version="MobiVision vv4.1"、version
  ```
- **支持输入**：fasta 支持 `.fasta/.fa.gz/.fna.gz`；gtf 支持 `.gtf/.gtf.gz`（自动解压，用完删除）。
- **染色体命名**：单基因组时 fasta/GTF 原样拷贝、染色体名不变（Ensembl 的 `1`/`MT` 与 GTF 天然匹配，直接可用）；多基因组（`-n` 重复指定）时给染色体加 `{genome名}_` 前缀防冲突，输出目录 `<g1>_and_<g2>`。
- **内存校验**：`-m` 低于最低需求（基因组 GB + SA 索引 GB + 3）时拒绝执行并提示 `--memgb=`；人类 3.1G 基因组最低约 8G，实测 `-m 64` 充裕。
- **STAR 调用**（mkindex 无线程参数，**固定单线程**）：人类基因组构建约 **2–3 小时**（STARtmp 体积可达数百 GB，写盘后自动清理）；GTF 解析（sjdb）成功以 Log.out 出现 junction 统计为准。
- **quantify 用法**：`-i <GENOME>/` 即索引路径。
- 官方参考是 Gencode（`chr1` 命名 + `gene_type`）；Ensembl 114（无 chr 前缀 + `gene_biotype`）完全兼容，两者均可。

### 2.3 用 10X 官方 refdata 重建 mobivision 索引（REF_SOURCE=refdata 的 mobivision 默认路径）

默认方案下消费官方 `refdata-gex-{GRCh38|GRCm39}-2024-A`（`REFDATA_DIR` 解压目录）——其 `fasta/genome.fa` + `genes/genes.gtf.gz` 已是 10X 官方过滤集合（白名单见 [reference-setup.md](reference-setup.md)）。即使不指定 `-t`、实质上保留全部注释，也必须先运行 `fmt_gtf` 完成 MobiVision 规范化，再将规范化输出交给 `mkindex`；不得直接把原始 GTF 交给 `mkindex`。

```
source $MOBIVISION_DIR/source.sh
mobivision fmt_gtf -i $REFDATA_DIR/refdata-gex-<REF>-$REFDATA_VERSION/genes/genes.gtf.gz \
  -o <normalized.gtf>
mobivision mkindex -n <GENOME> \
  -f $REFDATA_DIR/refdata-gex-<REF>-$REFDATA_VERSION/fasta/genome.fa \
  -g <normalized.gtf> \
  -r <REFDATA_VERSION>
```

命令不主动传入线程或内存参数，由 MobiVision 默认值与系统级资源配置决定。产出目录作 `MOBI_INDEX`。`REF` 由 `SPECIES` 映射（human→GRCh38、mouse→GRCm39），详见 [reference-setup.md](reference-setup.md)。

## 3. quantify 全流程与中间产物生命周期

### 3.1 命令与输入

```
mobivision quantify -i $MOBI_INDEX -f <fastq_dir> -o <out_dir> -s <sample_id> --intron included
```

- **线程/内存默认由软件决定，不显式传 `-t`**（显式传 `-t` 曾致软件崩溃；系统级资源配置负责上限）。
- **fastq 命名**（墨卓平台标准）：`[SampleID]_[时间戳].R1.fastq.gz` / `.R2.fastq.gz`；SampleID 可含 `-`。正则 `(.+)[-_\.](R1|R2)(_0\d+)*.(fastq|fq)(.gz)?` 自动识别，**无需重命名**。
- **多 lane/多批次**：同一样本多次测序（如 `<S编号>_hm/`、`<S编号>_rg/` 等子目录，每个 1–2 对 fastq）必须**平铺**到同一个 `-f` 目录（`os.listdir` 不递归子目录）。不同 SampleID 的多对 fastq 官方支持自动 merge 合并为一对。
- **只复制不移动**：源 fastq 由脚本 rsync 复制到本地暂存（子目录→staging 内平铺），**绝不用 mv/剪切移动源文件**（原始数据只允许复制）。
- **大文件串行**：复制 >2G 大文件时不同时复制其他内容（单线程串行，否则机械硬盘崩溃）。
- **`--intron included`**：默认即 included（基因定量含内含子），显式写出以明确意图。

### 3.2 阶段链（logs/mobivision.log 关键字）

```
MOBIVISION v4.1 ANALYSIS STARTED
  → BARCODE-CORRECTING        # merge → gzip -t 校验(全部 merged) → seqkit split → barcode 校正(解压输出)
  → READS-TRIMMING            # cutadapt 修剪，输出 cellranger 风格 S0_L001
  → FILE-REARRANGING
  → READS-MAPPING             # STAR（比对 → BAMsort → 写 BAM → STARtmp 清理）
  → SATURATION-CALCULATING    # 最慢阶段：samtools view -q 255 -s 0.1/0.2/.../1.0 六点抽样
  → RESULTS-SUMMARIZING
  → REPORT-GENERATING         # outs 目录生成（BAM 拷贝 + 报告）
  → (samtools reheader 精简 BAM header → Aligned.sort.bam)
  → MOBIVISION v4.1 ANALYSIS FINISHED   # 同时写 _flagdone
```

### 3.3 中间产物与删除时机（磁盘管理核心）

| 产物 | 位置 | 大小参考 | 删除时机 |
|---|---|---|---|
| 源 fastq | `-f` 目录 | 与输入同 | **merge 完成后立即删**（merged 是唯一分析输入） |
| `merged/` | `-f` 目录下 | ≈ 源压缩总量 | **split 完成后立即删** |
| `split_fastq/` | `-o` 目录 | ~90–190G | `READS-TRIMMING FINISHED` 后删（cutadapt 完成后） |
| `convert_fastq/` | `-o` 目录 | **最大中间产物，可达数百 GB 至 TB 级**（解压版 starsolo.fastq） | `READS-TRIMMING FINISHED` 后立即删（cutadapt 以其为输入；BARCODE-CORRECTING 后删会致 cutadapt FileNotFoundError） |
| `merge_fastq/` | `-o` 目录 | ~300–470G（trimming 输出） | STAR 启动后 mobivision 自动删 |
| `map_result/`（含 `*__STARtmp`） | `-o` 目录 | STARtmp 体积可达数百 GB | STAR 完成（BAM 写出）后 STARtmp 自动删；BAM+Solo.out 保留 |
| `{sid}_outs/` | `-o` 目录 | BAM（数十至数百 GB）+ 矩阵 + 报告 | 交付验证后删 |

**经验**：convert_fastq 是"磁盘杀手"——它增长到 TB 级时若与其他任务叠加必然爆盘。但删除时机必须是 **READS-TRIMMING FINISHED**（cutadapt 以 convert_fastq 为输入，过早删会致 `FileNotFoundError` → 比对空 BAM 假成功）；在该阶段前盯着磁盘余量（低于安全余量时告警）。BAM 大小与测序深度成正比。

### 3.4 验证产出（阶段1 脚本已内置）

- `_flagdone` 存在
- `outs/{sid}_Aligned.sort.bam` + `.bai`（最终 BAM；中间 `sortedByCoord.out.bam` 会被 mobivision 清理）
- `outs/filtered_cell_gene_matrix/{barcodes,features,matrix}.mtx.gz`（过滤空液滴后）
- `outs/raw_cell_gene_matrix/` 同构（过滤前）
- `outs/{sid}_Report.html/.json`、`{sid}_summary.csv`、`{sid}_filtered.h5ad`、`{sid}.qlentille`

## 4. 严格流式调度铁律

> 以下四条是墨卓数据跑批的硬约束，违反会直接导致数据损坏或交付失败（有真实事故记录）。

1. **quantify 严格串行**：一个样本的 quantify 必须完全结束（`_flagdone`）才能启动下一个。**禁止并行**——两个 quantify 的中间产物（convert 可达数百 GB 至 TB 级）叠加必然磁盘爆满，STAR 在写盘失败时产出的 BAM 是**静默损坏**的（`Invalid BGZF header`，samtools 读时才发现），整轮必须重跑。
2. **样本完全隔离**：上一个样本的**所有**文件（fastq、merged、全部中间产物、outs、交付后残留）在 WSL/工作区完全清空后，才能开始复制下一个样本的 fastq。禁止提前复制。
3. **velocyto 阶段边界 + 并行有间隔**：样本表**全部**样本的 quantify 完成并验证后才进入 velocyto 阶段（02 脚本开头闸口逐样本核验 `_flagdone` 等完成标记，缺任一即拒跑）——quantify 独占整机资源且有 数百 GB 级中间产物磁盘压力，与 velocyto 重叠会破坏第 1/2 条的串行隔离前提。唯一例外：用户显式声明"每样本 quantify 完成即可直接接续该样本 velocyto"，此时 run_config 置 `VELO_PER_SAMPLE_HANDOFF=1`。velocyto 可并行（上限 4），**每 15 分钟发车一个**（`INTERVAL=900`），避免同时启动抢内存。
4. **禁止跨文件系统处理读写**：staging 与工作根（`FASTQ_STAGING`/`OUTPUT_ROOT`）必须在 WSL 原生 ext4，禁止 `/mnt/*`——WSL 经 drvfs/9p 访问 Windows 挂载盘经协议翻译，吞吐会掉一个数量级以上，长任务被拖到不可接受的时长（量级估计见 §5）。跨文件系统仅允许两类一次性传输：① 从源盘 rsync 复制 fastq 进 staging；② 完成后 rsync 交付 outs 回挂载盘。run_config source 时硬校验，违规拒跑。

## 5. 资源与时长参考

（线程数由软件默认决定，系统级配置负责资源上限）

| 阶段 | 相对开销 | 说明 |
|---|---|---|
| merge + gzip -t 校验 | 中 | 读全量 merged 做完整性校验 |
| split + convert | 高 | convert 输出为解压 fastq，占用最大 |
| trimming | 低 | |
| STAR 比对 + 排序 | 高 | 排序阶段 CPU 打满、内存占用大 |
| SATURATION 抽样 | 中–高 | 1.0 点最慢（全量 reads） |
| summarize + report + reheader | 低 | |
| **quantify 单样本合计** | 数小时至十余小时（随数据量与机型浮动） | |

> 上述时长基于工作根位于 WSL 原生 ext4；若误放 `/mnt/*`（drvfs/9p 挂载盘），单样本可超 20 h（见 §4 铁律4）。

- **velocyto**（run 模式，`-b -e`）：扫描 + counting + loom，每样本数小时级（随 BAM 体积与细胞数浮动）。
- **磁盘峰值预算**：单样本 quantify 峰值 ≈ 源 fastq + merged + split + convert + STARtmp + BAM，最坏可达 TB 级。严格串行并预留余量；交付（rsync 到挂载盘）期间 WSL 侧保留 outs 副本，注意双份占用。
- **velocyto 阶段磁盘**：两个 BAM 副本 + 两个 cellsorted 中间（各 ~70–80% BAM 大小）+ loom，峰值可达数百 GB。

## 6. 交付规范

交付目录结构（`<交付根>/{sample_id}/`，交付根由项目 run_config 决定，可为挂载盘）：

```
{sample_id}/
├── outs/                       # 必须叫 outs，不是 {sample_id}_outs
│   ├── {sample_id}_Aligned.sort.bam (+ .bai)
│   ├── filtered_cell_gene_matrix/   # 过滤空液滴后 counts
│   ├── raw_cell_gene_matrix/        # 过滤空液滴前 counts
│   ├── {sample_id}_Report.html/.json、{sample_id}_summary.csv
│   ├── {sample_id}_filtered.h5ad、{sample_id}.qlentille、cell_metrics.csv 等
├── velocyto/
│   └── {sample_id}.loom        # 必须 -e {sample_id} 才能得到此命名
├── _flagdone
├── logs/                       # mobivision.log / stdout.log / stderr.log
└── run_analysis_cmds.txt
```

- **loom 命名**：`velocyto run` **必须加 `-e {sample_id}`**，否则输出为 `<bam名>_<随机后缀>.loom`（不合规）。
- **禁止交付 cellsorted BAM**：velocyto 会在 BAM 同目录生成 `cellsorted_*.bam`（~70–80% BAM 大小）的临时文件，交付前必须删除。
- **大文件传输用 rsync 不用 cp**：`cp` 向 drvfs/9p 挂载盘（Windows 盘）写超大文件会报"无法分配内存"中断；用 `rsync -av --append-verify`（可断点续传，大文件亦可稳定续传）。
- **交付验证**：BAM 字节数与 WSL 源一致（`stat -c%s` 比对）；outs 内无 `cellsorted*`；loom 非空且名为 `{sid}.loom`。
- **交付教训**：`cp -r` 不带 `-r` 会静默跳过目录（logs 目录丢失且无法补回）——交付命令统一 `rsync -a` 目录 + 显式文件。

## 7. 故障排除

### 磁盘满导致 BAM 静默损坏（必须重跑）

- **现象**：mobivision.log 出现 `parallel: Error: Is the disk full?` / `cat: 写入错误: 设备上没有空间`；之后 `SATURATION-CALCULATING FAILED : could not convert string to float: ''`（读损坏统计文件）；或 samtools 报 `[E::bgzf_read_block] Invalid BGZF header`。
- **根因**：并行 quantify 或中间产物未及时清理导致磁盘 100% 满，STAR 排序分片写盘失败但不报错退出，最终 BAM 在损坏 offset 处截断。
- **处置**：删除全部损坏产物（convert/split/map_result/outs），**释放磁盘后重跑**该样本。BAM 无法局部修复。
- **预防**：遵守第 4 节串行铁律 + 第 3.3 节删除时机表；监控 `df -h`，<50G 即告警。

### saturation 阶段异常慢

- 1.0 抽样点读全量 BAM，属正常；各抽样点耗时逐步递增。不要中途 kill。

### cp 大文件到 Windows 盘报"无法分配内存"

- drvfs 写缓冲问题，用 `rsync -av --append-verify` 替代（见第 6 节）。

### velocyto 警告 BAM index 缺失

- `[E::idx_find_and_load] Could not retrieve index file` 仅警告，velocyto 不依赖 .bai，可忽略。

### fmt_gtf 输出重复 gene_version 属性

- Ensembl GTF 特有（`gene_version "01"` + 原值并存），STAR/quantify 不受影响，无需处理。

### mkindex 用错 `-t gene`

- `fmt_gtf -t gene` 匹配不到任何 biotype → 输出空 GTF → mkindex 报 GTF 异常。全量保留应**不传 `-t`**（见 2.1）。

## 8. 下游衔接（阶段3-4）

mobivision quantify 产出的矩阵是 **mtx.gz 三件套 + scanpy h5ad**，**不是 10x 的 `raw_feature_bc_matrix.h5`**。进入 cellbender / 质控（阶段 3–4）前的衔接由 `03_cellbender.sh` 内置自动完成：

1. **格式规范化（02b）**：`scripts/02b_format_matrices.py` 把 `raw/filtered_cell_gene_matrix/` 复制为标准命名 `raw/filtered_feature_bc_matrix/`，并用 `make_10x_h5.py` 基于 mtx 构建对应 `.h5`（10x v3 格式；features/id = features/name = make.unique(基因符号)、空符号回退 ENSG——与 Seurat::Read10X(gene.column=2) + DropletUtils::write10xCounts(version="3") 标准组合输出逐项等价）。CellRanger 分支原生四件齐全，此步为 no-op。
2. **三位一体校验（03b，强制）**：`scripts/03b_validate_triple_consistency.py --sid <sid> --cr-root <root>`——mtx 与 h5 行/列（barcode）/表达值逐值一致（h5 行名兼容两语义：CellRanger 原生 ENSG id，或 02b 构建的符号 id=name）；loom 与二者基因（ENSG）集合、barcode（归一后）集合一一对应（无重复/缺失，顺序可不同——velocyto 按自身顺序输出，bridge 按 ENSG+barcode 严格重排）、shape 一致且 spliced/unspliced 层完整非空（loom 剪接分层语义不同于 GeneFull：`s+u+a ≠ counts`，不做逐值比较；含显式 GeneFull 层时按坐标映射重排后逐值比较）；不一致中止并报告（禁止通过/跳过）。
3. **bridge 适配（04，仅 scautoqc 分支）**：`04_bridge_to_starsolo.py` 用 loom 的 `Accession`（ENSG）对齐 mtx 目录 features.tsv 第 1 列（ENSG/ENSMUSG，两比对分支恒定），不依赖 h5 行名语义；三位一体坐标对齐是前提。
4. **handoff**：tech/chemistry 由 `06_scautoqc_pipeline.sh` 按 `ALIGNER=mobivision` 自动推导（`tech=mobivision`、`chemistry=mobi_3prime_v1`）并传给写 HANDOFF 的脚本；5' 试剂盒把 run_config `CHEMISTRY` 填为 `mobi_5prime_v1` 覆盖（`TECH` 同理可覆盖）。

> 遇阻塞要求用户介入，不自行绕过。
