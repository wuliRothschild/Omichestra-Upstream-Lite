# 排错（已知坑与解法）

> MobiVision 分支（`ALIGNER=mobivision`）的坑详见 [mobivision-platform.md](mobivision-platform.md) §7；本文档列 cellranger 分支与共用坑。

## 磁盘满导致 mobivision BAM 静默损坏（必须重跑）

**现象**：mobivision.log 出现 `parallel: Error: Is the disk full?` / `cat: 写入错误: 设备上没有空间`；之后 `SATURATION-CALCULATING FAILED : could not convert string to float: ''`；或 samtools 报 `[E::bgzf_read_block] Invalid BGZF header`。

**原因**：并行跑多个 quantify 或中间产物（convert_fastq 可达数百 GB 至 TB 级）未及时清理导致磁盘 100% 满，STAR 排序分片写盘失败但不报错，最终 BAM 静默损坏。

**解法**：删除全部损坏产物（convert/split/map_result/outs）后重跑。**预防**：quantify 严格串行 + 阶段监听清理（01_mobivision_quantify.sh 内置）+ 磁盘 <50G 告警。详见 mobivision-platform.md §3–4。

## mobivision 交付盘大文件 cp 报"无法分配内存"

**现象**：`cp` 大体积 BAM 到 Windows 挂载盘（drvfs/9p）时报 `无法分配内存` 中断。

**解法**：用 `rsync -av --append-verify`（可断点续传）。`01_mobivision_quantify.sh` 交付段已用 rsync。

## mobivision loom 命名不合规（`{sid}_Aligned_<随机>.loom`）

**原因**：`velocyto run` 未传 `-e <sid>`，loom 以 BAM 名 + 随机后缀命名。

**解法**：必须 `-e <sid>`（`02_velocyto.sh` mobivision 分支已内置）。交付前验证 `velocyto/{sid}.loom` 存在。

## fmt_gtf 用错 `-t gene`

**现象**：`fmt_gtf -t gene` 输出空 GTF（gene 不是 biotype 值），mkindex 报 GTF 异常。

**解法**：全量保留应**不传 `-t`**；`-t` 只接受 biotype 值（protein_coding/lncRNA 等）。详见 mobivision-platform.md §2.1。

## scanpy 布尔列绘图崩溃

**现象**：`qc.py` 跑完 QC 计算（"N/M pass"）后，在 `generate_qc_plots` 崩：`TypeError: The numpy boolean negative, the '-' operator, is not supported`。

**原因**：scanpy ≤1.9.x 的 `scatterplots.py` 对布尔型 color（`consensus_passed_qc` 等）用 `-color_vector` 排序，numpy 拒绝布尔取负。

**解法**：用 scanpy **1.10.4**（见 environment-setup.md）。sctk 文档明确 1.10 修复。不要用 nf-scautoqc Dockerfile 的 1.9.3。

## `np.float_` AttributeError

**现象**：`import scanpy` 即 `AttributeError: np.float_ was removed in the NumPy 2.0 release`。

**原因**：scanpy 1.9.3 的 `_morans_i.py` 用 `np.float_`，numpy 2.0 移除。

**解法**：numpy 降到 1.26.4 + scanpy 升到 1.10.4（1.10 已修此 bug，且兼容 numpy 1.26）。

## flag_doublet.py 静默 exit（无 scrublet.csv 产出）

**现象**：`flag_doublet.py` 退出码 0 但目录里没有 `<sid>_scrublet.csv`。

**原因**：脚本内部判断 `(ad.X.sum(0) > 0).sum() < ad.shape[1] * 0.2` 时（>80% 基因零表达）直接 `exit()` 不产文件。

**解法**：检查该样本细胞数是否过低或矩阵异常。正常 10x 样本不会触发。若触发，复查 cellbender 输出是否正常。`06_scautoqc_pipeline.sh` 已对每个样本检查 csv 存在性。

## cellranger `--id` 非法字符

**现象**：cellranger 报 `Invalid name: <id> (only numbers, letters, dash, and underscore allowed)`。

**原因**：`--id` 经 martian 层硬校验，仅允许 `[A-Za-z0-9_-]`（无首字符限制）。原始 illumina name 含 `+`/`.`/空格 会失败。

**解法**：sample_id 用合法命名（见 sample-sheet.md）。脚本会用 sample_id 作 --id。

## 9p/挂载盘 I/O 极慢（含工作根误放挂载盘）

**现象**：cellranger/mobivision/velocyto 跑得极慢，或长时间无进度。实测 mobivision quantify 工作根放挂载盘（drvfs），单样本可被拖到十几小时以上（原生盘快一个数量级）。

**原因**：WSL 进程读写 `/mnt/*`（drvfs/9p 挂载的 Windows 盘）经协议翻译，随机读写性能差 10–50 倍——不仅是输入输出文件，处理期工作根（`OUTPUT_ROOT`/`FASTQ_STAGING`）位于 `/mnt/*` 时，全部中间产物 I/O 都被拖慢。

**解法（硬规则）**：处理期工作根必须位于 WSL 原生 ext4（如 `$HOME/...`）；run_config source 时对 `OUTPUT_ROOT`/`FASTQ_STAGING` 硬校验，位于 `/mnt/*` 直接拒跑。跨文件系统仅允许两类一次性传输：① 从源盘 rsync 复制 fastq/BAM 进本地暂存；② 完成后 rsync 交付产出回挂载盘。各阶段已内置 rsync 到 `$FASTQ_STAGING`。

## cellbender ckpt.tar.gz 残留

**现象**：`ckpt.tar.gz`（1GB+）出现在工作目录。

**原因**：cellbender 默认每 10 分钟存检查点到 cwd。

**解法**：`03_cellbender.sh` 末尾已 `rm -f ckpt.tar.gz`。若手动跑，记得清理。

## cellbender HTML 报告 Traceback

**现象**：cellbender 日志末尾 `FileNotFoundError: tmp.report.nbconvert.html`。

**说明**：这是 cellbender 0.3.2 生成 HTML 报告的尾部 bug（nbconvert 临时文件路径），**不影响 h5 矩阵产出**。日志应有 `Completed remove-background.` 在 Traceback 前后。

## 肝/实体组织 QC 通过率

**现象**：肝/心脏等实体组织样本 QC 通过率常低于脾/血样。

**说明**：实体组织（肝/心/肌肉）细胞正常有线粒体基因占比 5–25%（代谢活跃）。**original 模式（mito循环20/50/80，终用mito80）对实体组织友好**——只有 >80% mito 的濒死细胞才判 fail，正常实质细胞保留。multires 固定 mito20 会把代谢活跃的正常实质细胞误杀，实体组织通过率显著偏低；遇此情况改 original。若 original 下仍异常偏低，复查 QC 图（`scautoqc/2_qc_plots_individual/*.png`）确认是真实低质还是阈值偏严；必要时调 `--gmm_cutoff outer` 或提供 `--metrics_csv` 自定义阈值。

## nf-scautoqc genes_list 被覆盖

**现象**：finalize_qc 的 `var['cc'/'ig'/'tcr']` 全 False。

**原因**：nf-scautoqc 自带的是人类符号，对小鼠失配；或 genes_list 目录被其它操作覆盖。

**解法**：每次跑 `06_scautoqc_pipeline.sh` 前会先调 `05_prepare_genes_list.py` 重新生成物种对应列表到 `$NF_SCAUTOQC_DIR/genes_list/`。确认 `SPECIES` 配置正确。

## var_names 不唯一警告

**现象**：scanpy 警告 `Variable names are not unique`。

**说明**：cellranger GRCm39 参考有 37 个重复基因符号（Gm*）。桥接脚本用 Ensembl gene_id 对齐规避。pool_all 会 `var_names_make_unique`。不影响结果。

## scDblFinder 相关

**schard 读 postqc 的 assay 处理**：postqc 含 .X(cellbender) + layers(raw/spliced/unspliced/ambiguous)。`schard::h5ad2seurat(use.raw=FALSE, assay="RNA")` 只把 .X 读进 RNA counts；`as.SingleCellExperiment` 把 RNA→主 counts（scDblFinder 只读主 counts），其余 layers 进 altExps 被忽略——符合预期，scDblFinder 只需 counts。

**scDblFinder 输出列名**：1.18+ 产出 `scDblFinder.score`/`.class`/`.sample`/`.weighted`/`.cxds_score` 五列，`08_scdblfinder.R` 对默认/同型两次运行各取 `.score`/`.class`。

**per-sample csv 合并对齐**：`08` 输出 barcode + sampleID + 双运行 4 列（`score_scdblfinder`/`doublet_scdblfinder`/`score_scdblfinder_homotypic`/`doublet_scdblfinder_homotypic`，被 06d 阈值闸排除的细胞 score=1/class=singlet），`09_merge_scdblfinder.py` 拼 `barcode-sampleID` 对齐 pooled obs_names（AnnData.concatenate 的 `index_unique='-'` 格式），并合入 `06e` 的簇级检验 4 列（`scdblfinder_bh_val`/`scdblfinder_homotypic_bh_val`/`scdblfinder_stringent`/`scdblfinder_homotypic_stringent`；缺 `<sid>_scdblfinder_cluster_qc.csv` 或缺列即硬报错）+ 阈值闸列，同时在交付层把 scrublet 的 `bh_pval`/`stringent_doublet` 更名为 `scrublet_bh_val`/`scrublet_stringent`（幂等，重跑跳过）。若对齐失败（score NA），检查 barcode 是否含 `-1` 等后缀残留。

**06e 簇级检验的兜底与真失败区分**：聚类/MAD 检验异常时 06e 只降级不中断（该运行 bh_val 全 1 + 日志 `[warn]`，csv 正常落盘，退出码 0）——这是设计行为（bh_val=1 即"无簇级证据"，stringent 退化为纯 score>0.5）；输入文件缺失/列缺失/barcode 不覆盖则退出码非零（真失败，06 编排层中止）。09/10 报"缺 `<sid>_scdblfinder_cluster_qc.csv`"即该 csv 未产出或未满四列——补跑 06e（命令见 stage-4-scautoqc.md §9.2 步骤 4c）后重跑 09/10 即可，无需重跑 08。

**低覆盖细胞报错（size factor/拟合错误）**：官方 FAQ——输入不得含空 droplet 或极低覆盖细胞（<200/500 reads 可致报错）。鉴定输入已由 06d 阈值闸把关（counts≥200/genes≥200/mito<30% sc/<5% sn，Science 2026 口径）；若仍报错，看日志 `[scdblfinder]` 行的阈值闸剔除数与 `2_qc_objects/<sid>_threshold_gate.csv` 的逐项判定，评估该样本低覆盖占比是否异常；确属阈值不适配时用 run_config `GATE_*` 校准（改动入指纹即触发重跑）。

**pooled0 不写 scDblFinder**：pooled0 是 finalize 前对象（pool_all 产），scrublet 在 finalize 才 merge，故 pooled0 无 scrublet 也无 scDblFinder——阶段一致。pooled/filtered（finalize 后）才有两套 doublet 标注。

**scrublet_score 不在 postqc**：postqc.h5ad 无 `scrublet_score`（qc.py 的 `DEFAULT_QC_METRICS` 不含它，不跑 scrublet）；scrublet 由 `flag_doublet.py` 单独跑，结果在 `scautoqc/3_doublet_scores/<sid>_scrublet.csv`，finalize 才 merge 进 pooled。

