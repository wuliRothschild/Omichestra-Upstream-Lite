# 环境配置详解

本流程用三个独立 conda 环境，互不污染。

## scautoqc 环境（阶段4 + 桥接脚本 04/05）

scanpy **1.10.4** 是关键——非 nf-scautoqc Dockerfile 的 1.9.3。原因：
- scanpy 1.9.3 的 `_morans_i.py` 用 `np.float_`，numpy 2.0 已移除，import 即崩
- scanpy 1.9.3 的 `scatterplots.py` 对布尔型 color（如 `consensus_passed_qc`）用 `-color_vector` 排序，numpy 拒绝布尔取负 → 绘图崩溃
- sctk 官方 notebook 明确：`# this won't be necessary in scanpy 1.10.0, booleans will become directly plottable`

```bash
conda create -n scautoqc python=3.10 -y
conda activate scautoqc
pip install "scanpy==1.10.4" "numpy==1.26.4" "scikit-learn==1.5.0" \
  scrublet celltypist statsmodels adjusttext numpy-groupies harmonypy \
  "leidenalg>=0.10,<0.11" matplotlib-venn seaborn umap-learn loompy
pip install "git+https://github.com/Teichlab/sctk.git"
```

验证：`python -c "import sctk,scanpy,scrublet,celltypist;from sctk import _pipeline;print('ok')"`，确认 `_pipeline.multi_resolution_cluster_qc` 存在。

- sctk 0.2.4（master 分支，Poetry 管理，依赖 scanpy ^1.9.3 即 >=1.9.3,<2.0，1.10.4 满足）
- sctk 的 `run_scrublet` 隐式 `import scrublet`（独立包，非 scanpy.external），必须装 scrublet 包

## velocyto 环境（阶段2）

```bash
conda create -n velocyto python=3.10 -y
conda activate velocyto
conda install -y -c bioconda samtools
pip install velocyto.py loompy
```

验证：`velocyto --version`（0.17.17）。velocyto 依赖 samtools，必须 conda 装。

## cellbender 环境（阶段3，需 GPU）

可复用任一已装齐 cellbender + CUDA 版 torch 的现有环境（run_config 的 `ENV_CELLBENDER` 指向它）。若新建：

```bash
conda create -n cellbender python=3.10 -y
conda activate cellbender
pip install cellbender
# 按 GPU CUDA 版本装 torch（示例 cu126）
pip install torch --index-url https://download.pytorch.org/whl/cu126
```

验证：`python -c "import cellbender,torch;assert torch.cuda.is_available()"`。cellbender 必须 GPU，单样本约 1 小时（12GB 级显存消费级 GPU 即可跑默认参数）。

## scDblFinder 环境（阶段4 步骤4b，需 R）

scDblFinder 是 Bioconductor R 包，作为 scrublet 之外的第二条独立 doublet 检测线（per-sample postqc，鉴定输入由 06d 阈值闸把关后默认+同型双运行，标注不过滤）。Bioconductor 发行版由 BiocManager 按 R 版本自动匹配（R 4.5.x → Bioc 3.21+）。

**不必新建专用 env**——任一已装齐 `scDblFinder` + `schard` + `SingleCellExperiment` 的 R env 即可（如 `seurat` env）。先检查现有 R env：
```bash
# 替换 <env> 为现有 R env（如 seurat/r_env）
<env>/bin/Rscript -e 'cat(requireNamespace("scDblFinder",quietly=T),
                       requireNamespace("schard",quietly=T),
                       requireNamespace("SingleCellExperiment",quietly=T),"\n")'
# 输出 TRUE TRUE TRUE 即可用，run_config 的 RSCRIPT 指向该 env 的 Rscript
```

若现有 env 缺库，补装（不新建 env）：
```bash
conda activate <env>   # R 4.5.3 的 env
R -e '
  install.packages(c("reticulate","remotes"))
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
  BiocManager::install(c("scDblFinder","SingleCellExperiment"))
  remotes::install_github("cellgeni/schard")   # schard：h5ad→Seurat 读取（reticulate-free，依赖 rhdf5）
'
```

或新建专用 env（可选，仅当无合适 R env 时）：
```bash
conda create -n scdblfinder r-base=4.5.3 -y
conda activate scdblfinder
R -e '  install.packages(c("reticulate","remotes"))
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
  BiocManager::install(c("scDblFinder","SingleCellExperiment"))
  remotes::install_github("cellgeni/schard")
'
```

**schard 源与安装注意**：现行 schard = GitHub 独源 `cellgeni/schard`（1.0.x，定位 "Reticulate-free"：直接依赖 Bioconductor `rhdf5`，不依赖 reticulate）；早年配方里的 `haoeric/scHarD`（reticulate 桥）属旧线，不要混装。Windows 侧安装推荐（纯 R、无编译）：从该仓库取目标 commit 的 codeload tarball（目标 SHA 可读任一既有安装 `DESCRIPTION` 的 `RemoteSha`，保持与既有环境同源），再 `install.packages(<tarball>, repos=NULL, type='source')` 装入目标库；装后以加载测试确认。

验证：`Rscript -e "cat(requireNamespace('scDblFinder',quietly=T), requireNamespace('schard',quietly=T), requireNamespace('SingleCellExperiment',quietly=T))"` → `TRUE TRUE TRUE`。

- run_config 的 `RSCRIPT` 指向该 env 的 Rscript；`SCANPY_PYTHON` 指向 scautoqc env 的 python（视 schard 版本：reticulate 版经 SCANPY_PYTHON 读 h5ad；reticulate-free 的 1.0.x 直读，不需要）。
- scDblFinder 输出 `scDblFinder.score`/`.class`/`.sample`/`.weighted`/`.cxds_score`（08 脚本对默认/同型两次运行各取 score/class，产出 4 列）。
- **R 调用必须带 UTF-8 locale**：以 `LC_CTYPE` 前缀调用（06 已内置 `R_LC_CTYPE`，默认 `en_US.UTF-8`；WSL/容器可用 `C.UTF-8` 覆盖）。Windows R 4.5.3（UCRT）下 `C.UTF-8` 名称无效，若 locale 落到 GBK/C（如 GBK 中文系统），含中文的 R 源解析与日志输出会乱码/报错；08 脚本内置 locale 守卫（首个语句切 UTF-8）与纯 ASCII 回退路径兜底。

## R environment (for scDblFinder)

scDblFinder runs in R. Use any R installation with `scDblFinder` and `schard` available; point
`RSCRIPT` at it in `run_config.sh`. `08_scdblfinder.R` is plain ASCII and reads paths from the
command line, so it works on any OS/R version that satisfies the package requirements.

If you keep R on a different filesystem or host, expose the same work tree to it and pass
absolute paths; keep the script ASCII-only to avoid locale parsing issues.
