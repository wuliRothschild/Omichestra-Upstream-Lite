#!/usr/bin/env bash
# ============================================================================
# 交付归档 + 本地工作根清理（收口末步：先全量核验，后删除；默认 dry-run）
# 用法: 09c_cleanup.sh [run_config.sh] --deliver-cr <交付根> [--deliver-sc <根>]
#        [--verify-only] [--purge-staging] [--exclude <正则>]... [--apply]
# 契约（关键先核验后删除）：
#  1) 归档：每对 本地根→交付根 先 rsync -a 补齐（幂等增量；瞬时件不带过去；
#     --verify-only 跳过本步只核验，用于交付已另行完成的场景）
#  2) 核验：本地每个文件在交付端存在且字节数一致——任一缺失/不符即打印明细
#     并非零退出：交付不完整，禁止清理（先补齐交付再复跑）
#  3) 清理（仅 --apply 执行）：已核验的本地根、FASTQ_STAGING/<PROJECT>（非空需
#     --purge-staging）、两侧瞬时件（__pycache__/*.pyc/
#     .Rhistory/.RData/ckpt.tar.gz）
# 删除防护：本地根禁止 空/"/"/$HOME/层级<3 级/位于 FASTQ_SOURCE_ROOT 内（原始数据只复制不移动）
# 分工：09b=交付契约的语义校验（该有的有没有）；09c=树级归档核验+清理闸口（一件不漏才准删）
# 注：--exclude 仅在 --verify-only 场景有意义（本地有意保留、不要求交付端覆盖的文件，
#     正则按相对路径匹配）；默认归档模式为全树 rsync，本地一切产出都会进交付根。
# ============================================================================
set -uo pipefail
export LC_ALL=C
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CFG="run_config.sh"
DEL_CR=""; DEL_SC=""
APPLY=0; VERIFY_ONLY=0; PURGE_STAGING=0
EXCLUDES=()
while [ $# -gt 0 ]; do
  case "$1" in
    --deliver-cr) DEL_CR="$2"; shift 2 ;;
    --deliver-sc) DEL_SC="$2"; shift 2 ;;
    --verify-only) VERIFY_ONLY=1; shift ;;
    --purge-staging) PURGE_STAGING=1; shift ;;
    --exclude) EXCLUDES+=("$2"); shift 2 ;;
    --apply) APPLY=1; shift ;;
    *) [ -f "$1" ] && CFG="$1"; shift ;;
  esac
done
[ -f "$CFG" ] || CFG="$SCRIPT_DIR/$CFG"
# shellcheck disable=SC1090
source "$CFG"

CR_ROOT="${OUTPUT_ROOT%/}/$PROJECT_NAME"
SC_DIR="${SC_ROOT:-$CR_ROOT/scautoqc}"; SC_DIR="${SC_DIR%/}"
MODE=$([ "$APPLY" = 1 ] && echo "apply" || echo "dry-run")
echo "==== [$(date '+%F %T')] 09c 交付归档与清理（$MODE）===="

# 瞬时件（可再生缓存/检查点）：不归档、不参与核验，--apply 时两侧扫除
TRANSIENT_RE='(^|/)(__pycache__)(/|$)|\.pyc$|(^|/)\.Rhistory$|(^|/)\.RData$|(^|/)ckpt\.tar\.gz$'
EX_RE=""
[ "${#EXCLUDES[@]}" -gt 0 ] && EX_RE="$(IFS='|'; echo "${EXCLUDES[*]}")"
filter_ex() { if [ -n "$EX_RE" ]; then grep -vE "$EX_RE"; else cat; fi; }

guard_local() {  # $1=路径 $2=标签；删除目标硬防护
  local p="${1%/}"
  case "$p" in
    ""|"/"|"$HOME"|"$HOME/"|"/mnt"|"/mnt/"|"${OUTPUT_ROOT%/}")
      echo "[err] $2 指向危险路径: '$1'，拒绝作为删除目标" >&2; exit 1 ;;
  esac
  [ -d "$p" ] || { echo "[err] $2 本地根不存在: $p" >&2; exit 1; }
  [ "$(awk -F/ '{print NF}' <<<"$p")" -ge 4 ] \
    || { echo "[err] $2 路径层级过浅（<3 级），拒绝作为删除目标: $p" >&2; exit 1; }
  if [ -n "${FASTQ_SOURCE_ROOT:-}" ]; then
    case "$p" in
      "${FASTQ_SOURCE_ROOT%/}"|"${FASTQ_SOURCE_ROOT%/}"/*)
        echo "[err] $2 位于 FASTQ_SOURCE_ROOT 内（原始数据只复制不移动），拒绝: $p" >&2; exit 1 ;;
    esac
  fi
}

# ---- 配对：本地根 -> 交付根 ----
LOCALS=(); REMOTES=(); LABELS=()
[ -n "$DEL_CR" ] || { echo "[err] 缺 --deliver-cr <交付根>（CR_ROOT 的交付映射必填）" >&2; exit 1; }
LOCALS+=("$CR_ROOT"); REMOTES+=("${DEL_CR%/}"); LABELS+=("CR_ROOT")
SC_HELD=""
case "${SC_ROOT:-}" in
  ""|"$CR_ROOT"/*) : ;;  # 默认位于 CR_ROOT 内，随 --deliver-cr 一并归档
  /mnt/*) echo "[note] SC_ROOT 位于挂载盘（本身即交付侧），不归档不清理: $SC_ROOT" ;;
  *) if [ -n "$DEL_SC" ]; then
       LOCALS+=("${SC_ROOT%/}"); REMOTES+=("${DEL_SC%/}"); LABELS+=("SC_ROOT")
     else
       echo "[note] SC_ROOT 独立于 CR_ROOT 且未给 --deliver-sc：保留不清理: $SC_ROOT"
       SC_HELD="$SC_ROOT"
     fi ;;
esac

# ---- 归档 + 全量核验 ----
FAIL=0
LMF="$(mktemp)"; RMF="$(mktemp)"
for i in "${!LOCALS[@]}"; do
  L="${LOCALS[$i]}"; R="${REMOTES[$i]}"
  guard_local "$L" "${LABELS[$i]}"
  case "$R" in "$L"|"$L"/*) echo "[err] 交付根不得位于本地根之内: $R" >&2; exit 1 ;; esac
  case "$L" in "$R"|"$R"/*) echo "[err] 本地根不得位于交付根之内: $L" >&2; exit 1 ;; esac
  if [ "$VERIFY_ONLY" = 0 ]; then
    echo ">>> 归档 ${LABELS[$i]}: $L/ -> $R/（rsync 幂等增量）"
    mkdir -p "$R"
    rsync -a --exclude='__pycache__/' --exclude='*.pyc' --exclude='.Rhistory' \
          --exclude='.RData' --exclude='ckpt.tar.gz' "$L/" "$R/" \
      || { echo "[err] rsync 归档失败: $L -> $R" >&2; exit 1; }
  fi
  if [ ! -d "$R" ]; then echo "  [FAIL] 交付根不存在: $R"; FAIL=1; continue; fi
  echo ">>> 核验 ${LABELS[$i]}: 本地每文件须在交付端存在且字节数一致"
  (cd "$L" && find . -type f -printf '%P\t%s\n') | grep -vE "$TRANSIENT_RE" | filter_ex | sort > "$LMF"
  (cd "$R" && find . -type f -printf '%P\t%s\n') | sort > "$RMF"
  MISS="$(comm -23 <(cut -f1 "$LMF") <(cut -f1 "$RMF"))"
  DIFFSZ="$(join -t $'\t' -1 1 -2 1 "$LMF" "$RMF" | awk -F'\t' '$2!=$3{print $1}')"
  NMISS=$(printf '%s' "$MISS" | grep -c .); NDIFF=$(printf '%s' "$DIFFSZ" | grep -c .)
  if [ "$NMISS" -gt 0 ] || [ "$NDIFF" -gt 0 ]; then
    FAIL=1
    echo "  [FAIL] ${LABELS[$i]} 交付不完整：交付端缺失 $NMISS 项、字节数不符 $NDIFF 项（前 50 条）："
    [ -n "$MISS" ] && printf '%s\n' "$MISS" | head -50 | sed 's/^/    缺: /'
    [ -n "$DIFFSZ" ] && printf '%s\n' "$DIFFSZ" | head -50 | sed 's/^/    异: /'
  else
    echo "  [OK] ${LABELS[$i]} 全量一致（$(wc -l < "$LMF") 个文件）"
  fi
done
rm -f "$LMF" "$RMF"
if [ "$FAIL" != 0 ]; then
  echo "[FAIL] 交付不完整——禁止清理。先补齐交付（复跑本脚本、去掉 --verify-only 可自动归档补齐），再复跑核验" >&2
  exit 1
fi
echo "[OK] 全量核验通过：交付端已覆盖本地全部产出"

# ---- 删除清单 ----
STG=""
[ -n "${FASTQ_STAGING:-}" ] && STG="${FASTQ_STAGING%/}/$PROJECT_NAME"
echo ">>> 清理清单（$MODE）："
for d in "${LOCALS[@]}"; do echo "    - $d"; done
[ -n "$STG" ] && [ -d "$STG" ] && echo "    - $STG（staging）"
[ -n "$SC_HELD" ] && echo "    保留（SC_ROOT 未给交付映射）: $SC_HELD"

if [ "$APPLY" != 1 ]; then
  echo "[dry-run] 未删除任何内容；确认清单无误后加 --apply 执行删除"
else
  for i in "${!LOCALS[@]}"; do
    guard_local "${LOCALS[$i]}" "${LABELS[$i]}"  # 删除前再防护一次
    rm -rf "${LOCALS[$i]}" && echo "  [del] ${LOCALS[$i]}"
  done
  if [ -n "$STG" ] && [ -d "$STG" ]; then
    if [ -z "$(ls -A "$STG" 2>/dev/null)" ] || [ "$PURGE_STAGING" = 1 ]; then
      rm -rf "$STG" && echo "  [del] $STG"
    else
      echo "  [warn] staging 非空（可能有中断样本的暂存副本），跳过: $STG（确需删除加 --purge-staging）"
      ls -A "$STG" | head -20 | sed 's/^/         /'
    fi
  fi
  # 空壳收尾：OUTPUT_ROOT/FASTQ_STAGING 若已空则一并移除（有其他项目内容则保留）
  rmdir --ignore-fail-on-non-empty "${OUTPUT_ROOT%/}" 2>/dev/null
  [ -n "${FASTQ_STAGING:-}" ] && rmdir --ignore-fail-on-non-empty "${FASTQ_STAGING%/}" 2>/dev/null
  # 交付侧瞬时件扫除（早前归档可能带过去的缓存）+ skill 脚本目录缓存
  for R in "${REMOTES[@]}"; do
    [ -d "$R" ] || continue
    find "$R" -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null
    find "$R" -type f \( -name '*.pyc' -o -name '.Rhistory' -o -name '.RData' -o -name 'ckpt.tar.gz' \) -delete 2>/dev/null
  done
  find "$SCRIPT_DIR" -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null
  echo "[OK] 清理完成"
fi

# ---- 执行者侧复核清单（skill 管不到的位置，agent 责任）----
echo "[执行者侧复核] 以下临时件由执行者负责，任务结束前确认无残留（SKILL.md「清理临时文件」第 3 层）："
echo "  - 后台编排日志（nohup/setsid 重定向）：今后一律落 \$OUTPUT_ROOT/\$PROJECT_NAME/logs/orchestration/ 随本脚本归档；"
echo "    本次若已在 \$HOME 等处产生，请手动删除"
echo "  - 一次性辅助脚本（mktemp 文件、dry-run 副本）：当会话用完即删"
echo "  - 系统临时目录（/tmp、Windows %TEMP%）中本次运行产生的文件"
