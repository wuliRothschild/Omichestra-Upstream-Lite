#!/usr/bin/env bash
# ============================================================================
# Omichestra-upstream 流程配置模板
# 用法：复制为 run_config.sh，填入实际值，所有阶段脚本 source 它
# 原则：所有路径在此集中定义，脚本零硬编码，便于跨机器/账号迁移
# ============================================================================

# ---- 项目 ----
PROJECT_NAME=""                           # 输出目录名（输出落在 $OUTPUT_ROOT/$PROJECT_NAME/），填项目名如 scRNA_project1
SPECIES="mouse"                           # mouse | human | custom（决定 genes_list 来源；02b 的 h5 基因组名映射目前只覆盖 human/mouse）
# 注：相对路径（如默认 Document/...）按脚本运行时的 cwd 解析——请在项目工作目录下运行各阶段脚本，或改填绝对路径
SAMPLE_SHEET="Document/sample_sheet.csv"  # 样本表 CSV（默认 Document/；见 assets/sample_sheet_template.csv）
HANDOFF_MD="Document/HANDOFF.md"          # 项目单一 HANDOFF（07 upsert ## upstream 阶段；默认 Document/，可改）

# ---- 比对平台（阶段1 分流依据；cellranger 默认，mobivision 需另填下组）----
ALIGNER="cellranger"                      # cellranger | mobivision

# ---- 参考基因组（00_setup_references.sh 可自动检查/构建）----
# 参考策略 REF_SOURCE 四分支（默认 refdata 官方包；详见 references/reference-setup.md）：
#   refdata（默认）：用 $REFDATA_DIR/refdata-gex-<GRCh38|GRCm39>-$REFDATA_VERSION/ 解压目录；
#                    cellranger 直接作 TRANSCRIPTOME，mobivision 从其中 fasta+genes.gtf 重建 mkindex
#   download：本地无官方包 → 按 REFDATA_URL（默认官方直链模板）下载并解压到 REFDATA_DIR
#   build：从头构建（cellranger mkgtf 白名单 + mkref；mobivision fmt_gtf 白名单 + mkindex），输出目录填 TRANSCRIPTOME/MOBI_INDEX
#   custom：用户自定义方案（如 mobivision fmt_gtf 不传 -t 全量保留），TRANSCRIPTOME/MOBI_INDEX 必填
REF_SOURCE="refdata"                      # refdata | download | build | custom
REFDATA_DIR="${HOME}/reference/refdata"   # 官方包存放/解压根（refdata/download 时用）
REFDATA_VERSION="2024-A"                  # 官方包版本（2024-A 默认；可改 2020-A 等）
REFDATA_URL=""                            # 官方包下载 URL；留空按物种生成 cf.10xgenomics.com 直链模板（human→GRCh38、mouse→GRCm39），可覆盖为镜像/内网
TRANSCRIPTOME=""                          # [cellranger] 参考目录（含 reference.json/fasta/genes/star）；REF_SOURCE=refdata 时自动定位官方包可不填，build/custom 时填 mkref 输出
MOBIVISION_DIR=""                         # [mobivision] mobivision-v4.1 安装根目录（含 source.sh）
MOBI_INDEX=""                             # [mobivision] mkindex 输出目录（含 fasta/genome.fa + genes/genes.gtf + star/ + reference.json）；REF_SOURCE=refdata 时自动定位重建可不填，build/custom 时填 mkindex 输出
MOBI_MEM_GB=""                            # [mobivision] mkindex -m 覆盖位（GB）；默认空=不显式传，由软件默认与系统级资源配置决定；显式填写后 00 指引命令附加 -m
VELO_REF=""                               # velocyto 参考目录（含 genes.gtf 与 rmsk.gtf；必须与比对阶段同一官方参考，见 reference-setup.md）

# ---- 工具与环境 ----
CELLRANGER_DIR=""                         # cellranger 安装目录（含 cellranger 可执行），留空则用 PATH
CONDA_SH="${HOME}/miniconda3/etc/profile.d/conda.sh"  # conda 初始化脚本
ENV_SCAUTOQC="scautoqc"                   # scautoqc 阶段 conda 环境
ENV_VELOCYTO="velocyto"                   # velocyto 阶段 conda 环境
ENV_CELLBENDER="scvi"                     # cellbender 阶段 conda 环境（需含 GPU torch）
NF_SCAUTOQC_DIR="${HOME}/nf-scautoqc"     # nf-scautoqc clone 位置（scautoqc 阶段依赖其 bin/）
RSCRIPT=""                                # 指向任一装齐 scDblFinder/schard 的 Rscript 绝对路径
SCANPY_PYTHON=""                          # scanpy 环境的 python（schard 经 reticulate 读 h5ad 用，如 scautoqc env 的 python）

# ---- 输入输出 ----
FASTQ_SOURCE_ROOT=""                      # 源 fastq 根目录（如 9p/网络挂载 /mnt/z/...；仅作 rsync 复制源，绝不做处理工作目录）
OUTPUT_ROOT="${HOME}/sc_analysis"         # 输出根目录（cellranger/velocyto/cellbender/scautoqc 产出；处理期工作根，WSL 下禁止 /mnt/*——见文件末尾跨文件系统防护；交付完成后由 09c_cleanup.sh 归档核验并清理（先核验后删除），不留本地副本）
FASTQ_STAGING="${HOME}/sc_staging"        # 本地高速暂存（rsync 中转，须在本地 NVMe/SSD；WSL 下禁止 /mnt/*；每样本跑完即清，全程结束由 09c_cleanup.sh 收尾）
# 编排日志约定：后台任务日志重定向到 $OUTPUT_ROOT/$PROJECT_NAME/logs/orchestration/（随 09c_cleanup 归档进交付根，不散落 $HOME）

# ---- 元数据（可选；写 outs/metadata.json 与 HANDOFF 的 tech/chemistry/tissue/species）----
# 十类信息（group/age/batch/datasets/sample/tech/gender/tissue/chemistry/species）由用户主动交代或 agent 阻塞询问；
# 用户明确不提供则一律 Unknown。per-sample 字段（group/age/gender/batch）不进本文件——
# 由 agent 收集后调 01c_write_metadata.py --meta-csv 写入；此处仅运行级单值（留空走推导/Unknown）
TECH=""                                   # 技术平台：留空按 ALIGNER 推导（cellranger→sc，mobivision→mobivision）
CHEMISTRY=""                              # 建库化学：留空按 ALIGNER 推导（10x_3prime_v3 / mobi_3prime_v1；按实际试剂盒显式填写更准）
TISSUE=""                                 # 组织（数据集级单值——一个数据集只允许一种 tissue×tech×chemistry 组合）；留空=Unknown
META_CSV=""                               # per-sample 元数据 csv（列 sample_id,group,age,gender,batch[,sample][,species]；species 列可选、缺省取运行级 SPECIES；gender 规范 F/M，常见写法自动归一）；留空则 per-sample 字段全 Unknown
DATASET_NAME=""                           # 数据集名：metadata datasets 字段/HANDOFF 段落标题后缀；留空=$PROJECT_NAME（同一 CR_ROOT 跑多个数据集时必填，防产物互覆盖）
SC_ROOT=""                                # scautoqc 产出根覆盖；留空=$OUTPUT_ROOT/$PROJECT_NAME/scautoqc（契约默认，与样本同级）

VALIDATE_PYTHON=""                        # 阶段3 三位一体校验用 python；留空=由 ENV_SCAUTOQC 推导

# ---- 资源（线程/内存默认由软件决定，不在命令中显式限制——显式限制曾致软件崩溃；
#      下列变量默认空=不传入任何命令；确需覆盖时由用户显式填写，脚本仅在非空时传入对应参数）----
LOCALCORES=""                             # 线程覆盖位：cellranger --localcores / mobivision quantify -t（默认空=软件默认线程）
LOCALMEM=""                               # 内存覆盖位（GB）：cellranger --localmem（默认空=软件默认内存）
VELO_THREADS=""                           # velocyto 每样本线程覆盖位（-@；默认空=velocyto 默认线程）
MAX_JOBS=4                                # velocyto 并发上限（外层编排，非进程内资源限制）
INTERVAL=900                              # velocyto 发车间隔秒（避免同时启动抢内存）
VELO_PER_SAMPLE_HANDOFF=""                # 阶段边界覆盖位：默认空=样本表全部样本 quantify 完成才放行 velocyto（02 开头闸口强制）；
                                          # 仅当用户显式声明"每样本 quantify 完成即可直接接续该样本 velocyto"时填 1

# ---- scDblFinder 前阈值闸（scautoqc 阶段，06d_threshold_gate.py；硬阈值口径）----
# 纯阈值硬过滤（不做 GMM——GMM 属 qc.py cellwise_qc 职责）：counts>=200、genes>=200、
# mito<30%（sc）/ <5%（sn）；scrublet 不参与过滤。留空=契约默认；实取值入指纹 params，
# 改动即触发受影响样本重跑（指纹判据）
# GATE_MIN_COUNTS=""                       # n_counts 下限（默认 200）
# GATE_MIN_GENES=""                        # n_genes 下限（默认 200）
# GATE_MAX_MITO_CELL=""                    # percent_mito 上限·sc（默认 30）
# GATE_MAX_MITO_NUCLEI=""                  # percent_mito 上限·sn（默认 5）

                                          
# ---- 校验配置完整性（脚本 source 时自动跑；缺失即中止 source 方脚本）----
_cfg_err() { echo "[config 错误] $1" >&2; exit 1; }
[ -n "$PROJECT_NAME" ] || _cfg_err "PROJECT_NAME 未填"
case "$REF_SOURCE" in
  refdata|download)
    [ -n "$REFDATA_DIR" ] || _cfg_err "REFDATA_DIR 未填（refdata/download 时：官方包存放/解压根）"
    ;;
  build|custom) ;;
  *) _cfg_err "REF_SOURCE 仅支持 refdata|download|build|custom，当前: $REF_SOURCE" ;;
esac
case "$ALIGNER" in
  cellranger)
    [ "$REF_SOURCE" = "refdata" ] || [ -n "$TRANSCRIPTOME" ] \
      || _cfg_err "TRANSCRIPTOME 未填（REF_SOURCE=build/custom 时填 mkref 输出目录；或先跑 00_setup_references.sh）"
    ;;
  mobivision)
    [ -n "$MOBIVISION_DIR" ] || _cfg_err "MOBIVISION_DIR 未填（mobivision 安装根，含 source.sh）"
    [ "$REF_SOURCE" = "refdata" ] || [ -n "$MOBI_INDEX" ] \
      || _cfg_err "MOBI_INDEX 未填（REF_SOURCE=build/custom 时填 mkindex 输出目录；或先跑 00_setup_references.sh）"
    ;;
  *) _cfg_err "ALIGNER 仅支持 cellranger|mobivision，当前: $ALIGNER" ;;
esac
[ -n "$VELO_REF" ] || _cfg_err "VELO_REF 未填"
[ -n "$SAMPLE_SHEET" ] || _cfg_err "SAMPLE_SHEET 未填"
[ -f "$SAMPLE_SHEET" ] || _cfg_err "样本表不存在: $SAMPLE_SHEET"
[ -n "$RSCRIPT" ] || _cfg_err "RSCRIPT 未填（scDblFinder 需要可用 R 环境）"
[ -n "$SCANPY_PYTHON" ] || _cfg_err "SCANPY_PYTHON 未填（schard 读 h5ad 需要）"

# ---- 跨文件系统防护（WSL 硬约束；source 时自动跑，违规即拒跑）----
# 原因：WSL 进程读写 /mnt/*（drvfs/9p 挂载的 Windows 盘）经协议翻译，吞吐会掉一个数量级以上，
# 长任务会被拖到不可接受的时长（量化实测见 references/mobivision-platform.md §5）。
# 处理期工作根（OUTPUT_ROOT/FASTQ_STAGING）必须位于执行环境原生文件系统（WSL=ext4，如 $HOME/...）。
# 跨文件系统仅允许两类一次性传输：① 从源盘 rsync 复制 fastq/BAM 进本地暂存；② 完成后 rsync 交付出库。
if grep -qiE 'microsoft|wsl' /proc/version 2>/dev/null; then
  for _v in OUTPUT_ROOT FASTQ_STAGING; do
    case "${!_v}" in
      /mnt/*) _cfg_err "$_v 禁止指向 Windows 挂载盘（/mnt/*）：${!_v}
  WSL 经 drvfs/9p 访问 Windows 盘吞吐极低。
  请将 $_v 改为 WSL 原生 ext4 路径（如 \$HOME/...）；
  跨文件系统仅允许 rsync 复制源数据进本地、rsync 交付产出回挂载盘。" ;;
    esac
  done
fi
