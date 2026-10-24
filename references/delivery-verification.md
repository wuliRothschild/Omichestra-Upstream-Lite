# 交付对象核验（scautoqc 产出）

> 用途：对 `scautoqc_pooled0 / scautoqc_pooled / scautoqc_pooled_filtered.h5ad` 做只读核验——回答"obs 应包含什么""体积为何如此""是否被改动"；改动前先建立基线。核验只用 h5py 读 obs 键值与 dataset 属性，**不加载表达矩阵**。列语义见 [handoff-contract.md](handoff-contract.md)、[stage-4-scautoqc.md](stage-4-scautoqc.md) §7。

## 1. obs 契约（直接当断言用）

终态 `scautoqc_pooled_filtered.h5ad` 的 obs = **79 个数据列 + `_index`（barcode-sampleID）**：

| 组 | 列数 | 说明 |
|---|---|---|
| 基础计量/占比 | 46 | 11 族（`n_counts`、`n_genes`、`n_counts_mito/ribo/hb`、`percent_mito/ribo/hb`、`percent_top50`、`log1p_n_counts`、`log1p_n_genes`）各含 `_raw`/`_spliced`/`_unspliced` 变体（44）+ `percent_spliced`、`percent_soup` |
| QC 判定 | 13 | `pass_auto_filter(_mito20/50/80)`、`pass_default`、`good_qc_cell_mito20/50/80`、`good_qc_cluster(_mito20/50/80)`、`qc_cluster` |
| 样本/技术 | 4 | `sampleID`、`cell_or_nuclei`、`cell_probability`、`droplet_efficiency` |
| doublet + 阈值闸 | 16 | scrublet 7 + scDblFinder 8 + `pass_threshold_gate` |

不含（别误查）：十类元数据（在 `outs/metadata.json`）；剪接派生量（filtered 无 layers，需要时从 pooled0 回算）。

- 对照 pooled0（63 数据列）：`filtered − pooled0` 应恰为上述 16 列，`pooled0 − filtered` 为空；pooled 与 filtered 列集合相同。
- 交付层改名：`bh_pval → scrublet_bh_val`、`stringent_doublet → scrublet_stringent`，旧名不得残留。
- 值级断言（逐格重算应 0 失配）：
  - `pass_threshold_gate` = `n_counts ≥ 200 & n_genes ≥ 200 & percent_mito < 30`（sn：< 5），以 `GATE_*` 实取值为准；
  - scrublet：`doublet` = `score > 0.3 | bh_val < 0.65`；`scrublet_stringent` = `score > 0.3 | bh_val < 0.05`；
  - scDblFinder（同型同理）：`scdblfinder_stringent` = `score > 0.5 | bh_val < 0.05`；
  - 被闸细胞回填：`score_scdblfinder*` = 1、`doublet_scdblfinder*` = False；被闸数 = pooled − 闸通过层（与 HANDOFF 三层计数对账）。
- filtered 是 pooled 的闸通过子集：`pass_threshold_gate` 应全真，且同名列与 pooled 一致（同一次写入）。

## 2. 存储与压缩重建

- 压缩由最后写盘者决定：全链 h5ad 写点统一 `compression="gzip"`（含 `09_merge_scdblfinder.py`）。压缩只改存储、不改内容；pooled0 含四层反而可能更小，别用体积反推内容。
- 压缩重建属受控写操作，须另行授权。做法：`h5repack -f GZIP=4` 生成临时件（ASCII 暂存目录）→ **逐格校验**（结构/属性/全 dataset 值；判据是值全等，不是字节全等）→ 硬链接保全原件 → 单次 `os.replace` 原子换名（原名不缺失、无半成品窗口；禁止先移走或先删）。校验器先经阳性 + 阴性（蓄意篡改须被捕获）回归；批量用清单 + 逐件日志，可断点续跑。`h5repack` 只接受 ASCII 路径。参考实现：`scripts/09d_h5repack_runner_atomic.py`；校验器：`scripts/09e_verify_h5ad_equal.py`（退出码 0 = 逐格全等）。
- h5py 比对坑：vlen 字符串标量 `[()]` 返回 `bytes`（非 ndarray），比较前 `np.asarray()`；根组属性也要比对；压缩后 contiguous 小数据集会被切成单块（读透明）。
- 查压缩参数：读各 dataset 的 `.compression`/`.dtype`；mtime 可佐证"最后写入后未再改动"。

## 3. 只读核验步骤

1. 任一具备 h5py + anndata 的环境均可；给原生绝对路径（不要把 MSYS 风格的 `/c/...` 传给原生 Windows 程序）。
2. h5py 坑：`obs.keys()` 含 `_index`；分类列（`sampleID`、`qc_cluster`、`cell_or_nuclei`、`scrublet_done`、`pass_default` 等）以 Group（`codes` + `categories`）存储，须按此解包，否则 `too many values to unpack`；判等要 codes 与 categories 双查。
3. 断言顺序：列集合差（16 列）→ 旧名零残留 → 三层行数对账 → 公式逐格重算 + 被闸回填 →（可选）与基线快照逐列对照。
4. 报告只输出数量/布尔/差异清单，不打印细胞级数值。
5. 转换/重建后的读取级终验：`anndata.read_h5ad(path, backed='r')` 核对 shape 与 obs 列数。
