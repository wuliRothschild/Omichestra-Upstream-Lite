#!/usr/bin/env bash
# ============================================================================
# 一次性部署：检查/构建比对参考（cellranger TRANSCRIPTOME 或 mobivision MOBI_INDEX）
#             与 velocyto 参考（VELO_REF）
# 参考策略按 run_config 的 REF_SOURCE 四分支（默认 refdata，见 references/reference-setup.md）：
#   refdata（默认）：定位 $REFDATA_DIR/refdata-gex-<GRCh38|GRCm39>-$REFDATA_VERSION/ 官方包
#                    cellranger 直接作 TRANSCRIPTOME；mobivision 从包内 fasta+genes.gtf 重建 mkindex
#   download：本地无官方包 → 打印官方下载 URL（REFDATA_URL 可覆盖）与解压指引
#   build：从头构建 → cellranger 用 mkgtf 白名单 + mkref；mobivision 先 fmt_gtf 规范化再 mkindex
#          （命令模板指针 references/cellranger-platform.md 与 references/mobivision-platform.md；无 -t 表示全量保留，不跳过规范化）
#   custom：用户自定义方案 → 仅校验 TRANSCRIPTOME/MOBI_INDEX 完整性
# 资源参数默认不主动传入（由软件默认值与系统级资源配置决定）；MOBI_MEM_GB 非空时 mkindex 指引附加 -m
# velocyto 参考（VELO_REF）：需要 genes.gtf + rmsk.gtf；缺失则给出下载指引
# 用法：bash 00_setup_references.sh [run_config.sh]
# 风格：只检查 + 打印指引，不自动下载/构建（10G+ 下载与数小时 mkindex 由用户确认后手动执行）
# ============================================================================
set -euo pipefail
CFG="${1:-run_config.sh}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[ -f "$CFG" ] || CFG="$SCRIPT_DIR/$CFG"
# shellcheck disable=SC1090
source "$CFG"

# SPECIES → refdata 名映射（human→GRCh38, mouse→GRCm39）
case "$SPECIES" in
  human) REFDATA_GENOME="GRCh38" ;;
  mouse) REFDATA_GENOME="GRCm39" ;;
  *) REFDATA_GENOME="" ;;
esac
# 官方包解压目录（REF_SOURCE=refdata/download 用）
REFDATA_PKG="${REFDATA_DIR}/refdata-gex-${REFDATA_GENOME}-${REFDATA_VERSION}"
# 官方下载 URL（REFDATA_URL 留空则按物种生成直链模板；可覆盖为镜像/内网）
if [ -z "${REFDATA_URL:-}" ]; then
  case "$REFDATA_GENOME" in
    GRCh38) REFDATA_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-${REFDATA_VERSION}.tar.gz" ;;
    GRCm39) REFDATA_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-${REFDATA_VERSION}.tar.gz" ;;
  esac
fi

# mkindex 内存覆盖位：默认空=不显式传 -m；用户显式填写 MOBI_MEM_GB 时指引附加
MOBI_M=""
[ -n "${MOBI_MEM_GB:-}" ] && MOBI_M=" -m $MOBI_MEM_GB"

echo "===== 比对参考检查（ALIGNER=$ALIGNER, REF_SOURCE=${REF_SOURCE:-refdata}）====="

ref_check_ready() {  # $1=参考目录
  [ -f "$1/reference.json" ] && [ -d "$1/star" ]
}

if [ "$ALIGNER" = "cellranger" ]; then
  # cellranger 可执行
  CR="cellranger"
  [ -n "$CELLRANGER_DIR" ] && CR="$CELLRANGER_DIR/cellranger"

  case "${REF_SOURCE:-refdata}" in
    refdata)
      if [ -n "${REFDATA_GENOME:-}" ] && [ -d "$REFDATA_PKG" ] && ref_check_ready "$REFDATA_PKG"; then
        echo "[ok] 官方参考就绪：$REFDATA_PKG"
        echo "  cellranger 直接用其路径作 TRANSCRIPTOME（无需重建）："
        echo "    TRANSCRIPTOME=$REFDATA_PKG"
      else
        echo "[missing] 官方参考 $REFDATA_PKG 不完整（缺 reference.json 或 star/）"
        echo "  将 REF_SOURCE 改为 download 走下载指引，或把已有官方包解压到 REFDATA_DIR；详见 references/reference-setup.md"
      fi
      ;;
    download)
      if [ -n "${REFDATA_GENOME:-}" ] && [ -d "$REFDATA_PKG" ] && ref_check_ready "$REFDATA_PKG"; then
        echo "[ok] 官方参考已存在：$REFDATA_PKG（REF_SOURCE=download 但本地已有，直接用）"
        echo "    TRANSCRIPTOME=$REFDATA_PKG"
      else
        echo "[missing] 官方参考 $REFDATA_PKG 未就绪，请下载并解压（~10GB）："
        echo "    mkdir -p $REFDATA_DIR"
        echo "    wget $REFDATA_URL -O $REFDATA_DIR/refdata-gex-${REFDATA_GENOME}-${REFDATA_VERSION}.tar.gz"
        echo "    tar -xzf $REFDATA_DIR/refdata-gex-${REFDATA_GENOME}-${REFDATA_VERSION}.tar.gz -C $REFDATA_DIR"
        echo "  下载直链可改 run_config 的 REFDATA_URL（覆盖为镜像/内网）"
      fi
      ;;
    build)
      echo "  cellranger 从头构建（官方过滤白名单，详见 references/cellranger-platform.md §2）："
      echo "    # 1. 下载 GENCODE primary_assembly 注释 + FASTA（如 gencode.v44.primary_assembly.annotation.gtf.gz）"
      echo "    # 2. mkgtf 过滤（17 项白名单示例，全部 biotype 见 references/reference-setup.md）"
      echo "    $CR mkgtf <gencode>.gtf <gencode>.filtered.gtf \\"
      echo "      --attribute=gene_biotype:protein_coding --attribute=gene_biotype:lncRNA --attribute=gene_biotype:antisense \\"
      echo "      --attribute=gene_biotype:IG_LV_gene --attribute=gene_biotype:IG_V_gene --attribute=gene_biotype:IG_V_pseudogene \\"
      echo "      --attribute=gene_biotype:IG_D_gene --attribute=gene_biotype:IG_J_gene --attribute=gene_biotype:IG_J_pseudogene \\"
      echo "      --attribute=gene_biotype:IG_C_gene --attribute=gene_biotype:IG_C_pseudogene \\"
      echo "      --attribute=gene_biotype:TR_V_gene --attribute=gene_biotype:TR_V_pseudogene --attribute=gene_biotype:TR_D_gene \\"
      echo "      --attribute=gene_biotype:TR_J_gene --attribute=gene_biotype:TR_J_pseudogene --attribute=gene_biotype:TR_C_gene"
      echo "    # 3. mkref（FASTA 不可 gzip，GTF 须含 exon；--memgb 须大于 FASTA GB 数）"
      echo "    $CR mkref --genome=<GENOME> --fasta=<genome>.fa --genes=<genome>.filtered.gtf --ref-version=$REFDATA_VERSION"
      echo "  构建后将输出目录填入 run_config.sh 的 TRANSCRIPTOME"
      ;;
    custom)
      if [ -n "$TRANSCRIPTOME" ] && ref_check_ready "$TRANSCRIPTOME"; then
        echo "[ok] cellranger 参考就绪：$TRANSCRIPTOME"
      else
        echo "[missing] cellranger 参考 $TRANSCRIPTOME 不完整（缺 reference.json 或 star/）"
        echo "  自定义方案需自行保证参考完整；若需全量/特殊 biotype，确保 GTF 过滤与构建正确"
      fi
      ;;
    *) echo "[err] 未知 REF_SOURCE: ${REF_SOURCE:-refdata}"; exit 1 ;;
  esac
else
  # mobivision 分支（fmt_gtf + mkindex；详见 references/mobivision-platform.md §2）
  case "${REF_SOURCE:-refdata}" in
    refdata)
      if [ -n "${REFDATA_GENOME:-}" ] && [ -d "$REFDATA_PKG" ] && ref_check_ready "$REFDATA_PKG"; then
        echo "[ok] 官方参考就绪：$REFDATA_PKG"
        echo "  mobivision 从官方包重建专用索引（不直接复用 cellranger star/，详见 mobivision-platform.md §2.3）："
        echo "    source $MOBIVISION_DIR/source.sh"
        echo "    mobivision fmt_gtf -i $REFDATA_PKG/genes/genes.gtf.gz -o <normalized.gtf>"
        echo "    mobivision mkindex -n <GENOME> \\"
        echo "      -f $REFDATA_PKG/fasta/genome.fa \\"
        echo "      -g <normalized.gtf> \\"
        echo "      -r $REFDATA_VERSION$MOBI_M"
        echo "  构建后将输出目录填入 run_config.sh 的 MOBI_INDEX"
      else
        echo "[missing] 官方参考 $REFDATA_PKG 不完整（缺 reference.json 或 star/）"
        echo "  将 REF_SOURCE 改为 download 走下载指引，或把已有官方包解压到 REFDATA_DIR；详见 references/reference-setup.md"
      fi
      ;;
    download)
      if [ -n "${REFDATA_GENOME:-}" ] && [ -d "$REFDATA_PKG" ] && ref_check_ready "$REFDATA_PKG"; then
        echo "[ok] 官方参考已存在：$REFDATA_PKG（REF_SOURCE=download 但本地已有，直接用）"
        echo "    按 refdata 分支重建 mkindex（见上）"
      else
        echo "[missing] 官方参考 $REFDATA_PKG 未就绪，请下载并解压（~10GB）："
        echo "    mkdir -p $REFDATA_DIR"
        echo "    wget $REFDATA_URL -O $REFDATA_DIR/refdata-gex-${REFDATA_GENOME}-${REFDATA_VERSION}.tar.gz"
        echo "    tar -xzf $REFDATA_DIR/refdata-gex-${REFDATA_GENOME}-${REFDATA_VERSION}.tar.gz -C $REFDATA_DIR"
        echo "  下载直链可改 run_config 的 REFDATA_URL（覆盖为镜像/内网）"
      fi
      ;;
    build)
      echo "  mobivision 从头构建（先 fmt_gtf 规范化，再 mkindex；详见 references/mobivision-platform.md §2.1）："
      echo "    source $MOBIVISION_DIR/source.sh"
      echo "    # 过滤模式：按官方白名单传入 -t；全量保留模式：省略 -t，但仍必须运行 fmt_gtf"
      echo "    mobivision fmt_gtf -i <gencode>.gtf -o <gencode>.normalized.gtf"
      echo "    # mkindex 不主动传入线程/内存参数（MOBI_MEM_GB 非空时为用户显式覆盖，附加 -m）"
      echo "    mobivision mkindex -n <GENOME> -f <genome>.fa -g <genome>.normalized.gtf -r $REFDATA_VERSION$MOBI_M"
      echo "  构建后将输出目录填入 run_config.sh 的 MOBI_INDEX"
      ;;
    custom)
      if [ -n "$MOBI_INDEX" ] && ref_check_ready "$MOBI_INDEX"; then
        echo "[ok] mobivision 参考就绪：$MOBI_INDEX"
      else
        echo "[missing] mobivision 参考 $MOBI_INDEX 不完整（缺 reference.json 或 star/）"
        echo "  自定义方案（如 fmt_gtf 不传 -t 全量保留）需自行保证参考完整；详见 references/mobivision-platform.md §2.1"
      fi
      ;;
    *) echo "[err] 未知 REF_SOURCE: ${REF_SOURCE:-refdata}"; exit 1 ;;
  esac
fi

echo
echo "===== velocyto 参考检查 ====="
VELO_OK=true
[ -f "$VELO_REF/genes.gtf" ] || { echo "[missing] $VELO_REF/genes.gtf"; VELO_OK=false; }
[ -f "$VELO_REF/rmsk.gtf" ] || { echo "[missing] $VELO_REF/rmsk.gtf"; VELO_OK=false; }
if $VELO_OK; then
  echo "[ok] velocyto 参考就绪：$VELO_REF"
else
  echo "  velocyto 参考需要："
  if [ "$ALIGNER" = "mobivision" ]; then
    echo "    genes.gtf — 复制 MOBI_INDEX 输出内的 genes/genes.gtf（fmt_gtf 规范化产物；与比对参考逐字节一致，02 脚本 SHA-256 校验依赖此）"
  else
    echo "    genes.gtf — 用 cellranger 参考 TRANSCRIPTOME/genes/genes.gtf（官方包 genes/genes.gtf.gz 解压；与比对参考逐字节一致）"
  fi
  echo "    rmsk.gtf  — RepeatMasker 重复区注释（UCSC Table Browser 下载，或 velocyto 的 get_rmsk.py 从 repeatmasker 文件生成）"
  echo "  下载与生成步骤见 references/reference-setup.md"
fi
