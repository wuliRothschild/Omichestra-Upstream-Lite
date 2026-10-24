#!/usr/bin/env bash
# ============================================================================
# 一次性部署：scautoqc / velocyto conda 环境 + clone nf-scautoqc
# 首次使用时运行。已存在的环境/仓库会跳过。
# 用法：bash 00_setup_envs.sh [run_config.sh]
# ============================================================================
set -euo pipefail
CFG="${1:-run_config.sh}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[ -f "$CFG" ] || CFG="$SCRIPT_DIR/$CFG"
# shellcheck disable=SC1090
source "$CFG"
source "$CONDA_SH"

# ---- scautoqc 环境（scanpy1.10 + sctk + scrublet + celltypist）----
if conda env list | awk '{print $1}' | grep -qx "$ENV_SCAUTOQC"; then
  echo "[skip] conda 环境 $ENV_SCAUTOQC 已存在"
else
  echo "[build] 创建 $ENV_SCAUTOQC（scanpy 1.10.4 + sctk）"
  conda create -n "$ENV_SCAUTOQC" python=3.10 -y
  conda activate "$ENV_SCAUTOQC"
  # scanpy 1.10 解决布尔列绘图崩溃与 np.float_；numpy 1.26 为 scanpy1.10+sctk 经验证组合
  pip install "scanpy==1.10.4" "numpy==1.26.4" "scikit-learn==1.5.0" \
    scrublet celltypist statsmodels adjusttext numpy-groupies harmonypy \
    "leidenalg>=0.10,<0.11" matplotlib-venn seaborn umap-learn loompy
  pip install "git+https://github.com/Teichlab/sctk.git"
  python -c "import sctk,scanpy,scrublet,celltypist;print('scautoqc env ok')" || { echo "[err] $ENV_SCAUTOQC import 失败"; exit 1; }
  conda deactivate
fi

# ---- velocyto 环境 ----
if conda env list | awk '{print $1}' | grep -qx "$ENV_VELOCYTO"; then
  echo "[skip] conda 环境 $ENV_VELOCYTO 已存在"
else
  echo "[build] 创建 $ENV_VELOCYTO（velocyto）"
  conda create -n "$ENV_VELOCYTO" python=3.10 -y
  conda activate "$ENV_VELOCYTO"
  # velocyto 依赖 samtools（conda 装）+ velocyto.py（pip）
  conda install -y -c bioconda samtools
  pip install velocyto.py loompy
  python -c "import velocyto;print('velocyto env ok')" || { echo "[err] $ENV_VELOCYTO import 失败"; exit 1; }
  conda deactivate
fi

# ---- cellbender 环境（复用现有，仅校验）----
if conda env list | awk '{print $1}' | grep -qx "$ENV_CELLBENDER"; then
  conda activate "$ENV_CELLBENDER"
  if python -c "import cellbender,torch;assert torch.cuda.is_available()" 2>/dev/null; then
    echo "[ok] $ENV_CELLBENDER 含 cellbender + GPU torch"
  else
    echo "[warn] $ENV_CELLBENDER 缺 cellbender 或 GPU torch，cellbender 阶段前请补装"
  fi
  conda deactivate
else
  echo "[warn] $ENV_CELLBENDER 不存在，cellbender 阶段需先建（参考 references/environment-setup.md）"
fi

# ---- clone nf-scautoqc ----
if [ -d "$NF_SCAUTOQC_DIR/.git" ]; then
  echo "[skip] nf-scautoqc 已 clone 于 $NF_SCAUTOQC_DIR"
else
  echo "[clone] nf-scautoqc → $NF_SCAUTOQC_DIR"
  git clone https://github.com/cellgeni/nf-scautoqc.git "$NF_SCAUTOQC_DIR"
fi
echo "[done] 环境部署完成"