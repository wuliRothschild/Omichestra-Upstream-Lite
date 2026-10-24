# 阶段 1-3 上游细节（cellranger / velocyto / cellbender）

## cellranger count（01_cellranger_count.sh）

### 关键参数
- `--nosecondary`：跳过二级分析（降维/聚类），下游用 Seurat/Scanpy自己做
- `--create-bam=true`：保留 BAM（velocyto 必需，否则无法做 RNA velocity）
- 线程与内存：默认由 Cell Ranger 与系统级资源配置决定；仅在用户明确要求时覆盖

### I/O 规避（重要）
fastq 若在 9p 挂载（WSL `/mnt/z`、`/mnt/c`）或网络盘，**必须先 rsync 到本地 NVMe/SSD** 再跑。直读挂载盘会让比对慢一个数量级以上。脚本已内置 rsync 到 `$FASTQ_STAGING`。

### fastq 命名约定
脚本假定 10x bcl2fastq 标准命名：`<prefix>_S1_L001_R{1,2}_001.fastq.gz` + 同名 `.md5`。
- `<prefix>` 是 illumina sample name（如 A0、Fe_0、N1）
- 脚本重命名前缀为 `sample_id`，使 cellranger `--sample=<sample_id>` 能命中
- 多 lane（L001/L002...）或多套 fastq：需调整脚本的 fastq 匹配逻辑（当前只处理单 lane L001）

### BAM 生命周期
cellranger 产出 BAM → 归档到 `$BAMDIR`（通常=源 fastq 目录）→ 删本地副本 → velocyto 阶段 rsync 回本地。这样 BAM 始终在源盘有一份原件，本地只在用时临时存在。

### `--id` 合法字符
cellranger `--id` 仅允许 `[A-Za-z0-9_-]`（无首字符限制）。`+`/`.`/空格 等非法。sample_id 校验在样本表填写时把关。

### 跑完自动清理（保留 outs / velocyto）
cellranger 跑完后，样本目录里会残留大量中间文件（`SC_RNA_COUNTER_CS/`、`extras/telemetry/`、`<id>.mri.tgz`、`_cmdline`/`_log`/`_finalstate`/`_perf`/`_filelist`/`_mrosource`/`_vdrkill` 等 `_*` 元数据）。01 脚本末尾用**白名单清理**：仅保留 `outs/` 与 `velocyto/`，其余顶层条目一律删除。保护逻辑：清理前确认 `outs/` 存在，否则跳过以免误删空目录。每项打印，可追溯。

> 对已跑完但未清理的历史样本，可对该样本目录单独执行同样的白名单清理（bash 或 python `os.listdir` 跳过 outs/velocyto，其余 `shutil.rmtree`/`os.remove`）。

## velocyto（02_velocyto.sh）

### 命令
`velocyto run10x -m rmsk.gtf <cellranger_samplefolder> <genes.gtf>`
- `-m`：mask 重复区（RepeatMasker rmsk.gtf），避免重复区产生假性 unspliced 信号
- 线程参数不主动指定，使用 velocyto 默认值
- 输入是 cellranger 样本目录（含 `outs/`），velocyto 自动找 BAM + filtered barcodes
- 产出 `<samplefolder>/velocyto/<samplefolder_name>.loom`

### 注意
- velocyto loom 的 cell 名格式为 `<sample>:<barcode>x`（如 `sample1:AACACACA...x`），桥接脚本会清洗为 16bp
- loom 标准结构：主层为空（全零占位），计数在 `spliced`/`unspliced`/`ambiguous` 三个命名层；行属性 `Accession`=ENSG、`Gene`=符号/ENSG 混合显示，列属性 `CellID`
- loom 行列顺序按 velocyto 自身 GTF 解析/BAM 顺序输出，**与比对矩阵行列序通常不同**；下游桥接按 ENSG+barcode 严格重排（集合必须一一对应，03b 校验强制）
- **`spliced+unspliced+ambiguous ≠ 比对矩阵 GeneFull counts`**：velocyto 计数模型不同（-m mask 重复区、分子须可明确分类为 exon/intron 来源、multi-gene 分子丢弃），总量系统性偏小、逐格不等。禁止用该关系式做逐值等价校验（会 FAIL 一切正常数据）
- 不主动指定线程或内存参数，使用软件默认值与系统级资源配置

### 并发调度
velocyto 单样本约 3 小时，可并发。建议脚本外层调度：4 并发 + 15 分钟发车间隔（同时启动会因 samtools sort 抢内存 OOM）。调度编排按"MobiVision/CellRanger 阶段全部完成后才开始、间隔 15 分钟发车"原则。

## cellbender（03_cellbender.sh）

### 命令
`cellbender remove-background --cuda --input raw_feature_bc_matrix.h5 --output cellbender_matrix.h5`
- **必须用 raw matrix**（不是 filtered），cellbender 需要空 droplet 估计背景
- `--cuda`：GPU 加速（必须）
- 串行：GPU 独占，并行不提速反慢（显存争抢）

### 检查点清理
cellbender 默认每 10 分钟存 `ckpt.tar.gz` 到当前工作目录。03 在子 shell 内 `cd` 到样本 `outs/cellbender/` 再跑 cellbender，使检查点固定落该目录并由脚本就地 `rm -f`（与执行者 cwd 无关；否则检查点散落调用目录，每样本残留 1GB+ 垃圾）。

### 已知尾部错误
cellbender 0.3.2 运行结束时会尝试生成 HTML 报告（`cellbender_matrix_report.html`），常因 nbconvert 临时文件路径失败抛 Traceback。**不影响 h5 矩阵产出**，核心结果已写出。

### 显存
12 GB 级消费级 GPU 足够默认参数。若加 `--expected-cells` 或提高 `--total-droplets` 需注意 OOM。
