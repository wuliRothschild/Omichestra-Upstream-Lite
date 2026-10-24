# ============================================================================
# scDblFinder 双细胞鉴定（per-sample postqc，与 scrublet 同级，第二条独立 doublet 线）
# 鉴定输入 = 06d_threshold_gate.py 阈值闸通过细胞（Science 2026 口径硬编码阈值：
#   counts>=200、genes>=200、mito<30% sc / <5% sn；纯阈值判定，无 GMM——GMM 属
#   qc.py cellwise_qc 职责）。被阈值闸硬过滤的细胞不参与鉴定，其 score 默认 1
#   （供下游按分数剔除），class 记 singlet
#   （算法未对其判定，不编造 doublet 结论）。
# 官方 FAQ 同型双胞规则："should you try to identify also homotypic doublets with
#   scDblFinder, be sure to not to use the cluster-based approach, and to set
#   removeUnidentifiable=FALSE" —— 每样本平行两次运行：默认运行
#   （removeUnidentifiable=TRUE，聚焦异型双胞）+ 同型运行（clusters=FALSE,
#   removeUnidentifiable=FALSE），输出四列：
#   score_scdblfinder / doublet_scdblfinder /
#   score_scdblfinder_homotypic / doublet_scdblfinder_homotypic
# 用法: Rscript 08_scdblfinder.R --postqc <postqc.h5ad> --sid <sampleID>
#                               --keep-csv <2_qc_objects/<sid>_threshold_gate.csv>
#                               --out <per-sample csv> --scanpy-python <path>
#                               [--dbr 0.1]
# 依赖: R 4.5.3 + scDblFinder + schard + reticulate + scanpy python（anndata）
# 输入: 2_qc_objects/<sid>_postqc.h5ad（gather+qc 产，cellbender .X，全量细胞）
#       + <sid>_threshold_gate.csv（06d 产，barcode 级 pass_threshold_gate mask）
# 输出: <out> csv（barcode, sampleID + 上述四列；覆盖 postqc 全量细胞），
#       供 09_merge_scdblfinder.py 合并写入 pooled/filtered obs
# ============================================================================
# --- LC_CTYPE → UTF-8（须先于一切；Windows R 非 UTF-8 locale 下 GBK 可能吞引号致解析错，
#     日志中文亦依赖之；按序探测，失败时告警但不中断）---
{
  .utf8_locale <- NULL
  for (.loc in c("C.UTF-8", "Chinese_China.utf8", "en_US.UTF-8")) {
    .cur <- suppressWarnings(tryCatch(Sys.setlocale("LC_CTYPE", .loc), error = function(e) ""))
    if (grepl("utf-?8", .cur, ignore.case = TRUE)) { .utf8_locale <- .cur; break }
  }
  if (is.null(.utf8_locale))
    warning("no UTF-8 LC_CTYPE available (current: ", Sys.getlocale("LC_CTYPE"),
            "): Chinese log text may be garbled", call. = FALSE)
}
suppressMessages({
  library(schard); library(SingleCellExperiment); library(scDblFinder); library(reticulate)
})


Sys.setenv("VROOM_CONNECTION_SIZE" = 9999999999 * 1024^2)
options(future.globals.maxSize = 9999999999 * 1024^2)



# ---- 命令行参数 ----
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NA) {
  i <- which(args == key)
  if (length(i) == 0) default else args[i + 1]
}
POSTQC        <- get_arg("--postqc")
SID           <- get_arg("--sid")
KEEP_CSV      <- get_arg("--keep-csv")
OUT           <- get_arg("--out")
SCANPY_PYTHON <- get_arg("--scanpy-python")
DBR           <- as.numeric(get_arg("--dbr", "0.1"))
stopifnot(!is.na(POSTQC), !is.na(SID), !is.na(KEEP_CSV), !is.na(OUT), !is.na(SCANPY_PYTHON),
          file.exists(POSTQC), file.exists(KEEP_CSV))

use_python(SCANPY_PYTHON, required = TRUE)
set.seed(1234)   # 人工双胞模拟含随机抽样，固定种子保证可复算



# ---- schard 读 postqc → SCE ----
# postqc 含 .X(cellbender) + layers(raw/spliced/unspliced/ambiguous)；
# schard 读 .X 进 RNA counts，as.SingleCellExperiment 把 RNA→主 counts（scDblFinder 只读主 counts）
ad <- schard::h5ad2seurat(POSTQC, use.raw = FALSE, assay = "RNA"); gc()
sce <- as.SingleCellExperiment(ad)

# ---- 阈值闸 mask（06d 产：Science 2026 口径硬过滤；鉴定输入=pass 细胞）----
gate <- read.csv(KEEP_CSV, colClasses = c("barcode" = "character"))
m <- match(colnames(sce), gate$barcode)
if (any(is.na(m))) {
  stop(sprintf("[scdblfinder] %s: 阈值闸 csv 缺 %d 个 barcode（06d 与 postqc 不同源？重跑 06d）",
               SID, sum(is.na(m))))
}
keep <- tolower(as.character(gate$pass_threshold_gate[m])) == "true"
if (any(is.na(keep))) {
  stop(sprintf("[scdblfinder] %s: 阈值闸 csv 的 pass_threshold_gate 列含无法解析值", SID))
}
if (sum(keep) < 1) {
  stop(sprintf("[scdblfinder] %s: 阈值闸（06d）后无细胞剩余，样本质量过低或阈值不适配（run_config GATE_* 可校准）", SID))
}
cat(sprintf("[scdblfinder] %s: n=%d | 阈值闸剔除 %d 细胞，%d 进入鉴定\n",
            SID, ncol(sce), sum(!keep), sum(keep)))
sce_f <- sce[, keep]

# ---- 两次平行运行：默认（异型聚焦）+ 同型（官方 FAQ 口径）----
sce_d <- scDblFinder(sce_f, dbr = DBR, verbose = FALSE)
sce_h <- scDblFinder(sce_f, clusters = FALSE, removeUnidentifiable = FALSE,
                     dbr = DBR, verbose = FALSE)
cd_d <- colData(sce_d); cd_h <- colData(sce_h)



# ---- 输出 per-sample csv（覆盖 postqc 全量细胞；阈值闸过滤细胞 score=1 / class=singlet）----
# barcode 用 postqc 原始 obs_names，09 合并时拼 -sampleID
out <- data.frame(
  barcode                       = colnames(sce),                # 纯 barcode
  sampleID                      = SID,
  score_scdblfinder             = 1,                            # 未鉴定默认 1
  doublet_scdblfinder           = "singlet",                    # 未判定不编造 doublet
  score_scdblfinder_homotypic   = 1,
  doublet_scdblfinder_homotypic = "singlet",
  stringsAsFactors              = FALSE
)
m <- match(colnames(sce_f), out$barcode)
out$score_scdblfinder[m]             <- cd_d$scDblFinder.score    # 0-1 doublet 得分
out$doublet_scdblfinder[m]           <- as.character(cd_d$scDblFinder.class)
out$score_scdblfinder_homotypic[m]   <- cd_h$scDblFinder.score
out$doublet_scdblfinder_homotypic[m] <- as.character(cd_h$scDblFinder.class)
write.csv(out, OUT, row.names = FALSE)
cat(sprintf("[scdblfinder] %s: 鉴定 n=%d na=%d/%d | dbl=%.4f dbl_homotypic=%.4f | 阈值闸 %d 细胞 score=1 -> %s\n",
            SID, ncol(sce_f), sum(is.na(cd_d$scDblFinder.score)), sum(is.na(cd_h$scDblFinder.score)),
            mean(cd_d$scDblFinder.class == "doublet", na.rm = TRUE),
            mean(cd_h$scDblFinder.class == "doublet", na.rm = TRUE),
            sum(!keep), OUT))
