# 仓库约定（CONTRIBUTING）

本文件是仓库级工程约定：改脚本、加脚本、写文档前读一遍。它约束"怎么改"，不描述"流程做什么"——流程看 `SKILL.md` 与 `references/`。

## 1. 零项目残留

仓库里只允许出现**方法、参数、占位符**，不允许出现任何具体项目的内容：

- 项目标识：样本名、生物学分组名、项目名、处理/造模名
- 项目数值：具体细胞数、基因数、维度、样本数、簇数
- 机器路径：用户名、盘符、conda 安装位置、内网地址
- 过程痕迹：修改记录、"更新版/已修正/v2"、会话编号
- 凭据：任何 key、token、base_url

示例一律用占位符（`<sid>`、`<path>`、`<N>`、`Sample1`、`TissueA`）。物种级通用参考数据（黑名单、基因注释导出、周期基因等）属资产，放 `assets/ref/<species>/`。

## 2. 自包含与缓存

- 每个仓库独立可用，不引用仓库外的文件。
- `scripts/` 只放源码；`__pycache__/`、`*.pyc`、`.Rhistory`、`.RData`、运行日志一律不进版本库（见 `.gitignore`）。
- 中间产物（临时目录、调试脚本、探测输出）用完即清，提交前 `git status` 必须干净。

## 3. 跨平台与编码

- 所有文本 UTF-8 **无 BOM**；脚本（`.sh` / `.py` / `.R`）一律 **LF**（`.gitattributes` 已强制）。
- 读取用户提供的表格（csv/tsv）必须容忍 **BOM 与 CRLF**——Excel/Windows 编辑过的文件几乎必然带这两样，这不是异常：
  - shell：`sed $'1s/^\xEF\xBB\xBF//' "$SHEET" | tr -d '\r' | awk ...`（bash ANSI-C 引号，字面写 BOM 字节即可）
  - python：`open(path, newline="", encoding="utf-8-sig")`
- 变量展开一律加引号（`"$VAR"`）：路径含空格时不加引号会断词。
- shell 只用 bash/POSIX 通用语法；`sed -i` 无备份后缀、`date -d`、`readlink -f` 等 GNU 专有写法不用，确需时在 `references/environment-setup.md` 声明依赖。
- 路径拼接用语言内建方法（`os.path.join`、`file.path`），不手拼 `/`。
- Windows 只作宿主：不引入 PowerShell / 计划任务 / 注册表等仅 Windows 的机制。

## 4. 路径与解释器

- 所有路径、解释器、环境名集中在 `scripts/run_config.template.sh`；脚本内零硬编码绝对路径。
- 脚本不裸调 `python` / `Rscript` 走 PATH；一律 `"$SCANPY_PYTHON"` / `"$RSCRIPT"`，或先 `conda activate "$ENV_X"` 再调用。
- 模板里解释器变量留空由用户填写（示例写在注释里），并配非空守卫。

## 5. 文档与产物落位

- 文档类（HANDOFF、sample sheet、metadata、报告、日志）默认 `<工作目录>/Document/`。
- 数据类（h5ad / loom / 中间矩阵）默认 `OUTPUT_ROOT`。
- 落位由 run_config 决定，脚本不自行决定。

## 6. 数据流原则

**大到小可以裁，小合大不能合。**

- 从上级对象（基因维度完整）按 barcode 取子集 → 安全，基因维度仍是上级的子集。
- 把多个下级子集 merge 回上级继续用 → 不安全：merge 后的基因是各子集裁剪后的并集，重新算 feature threshold 会让分布失真。
- 下游交接一律从上一轮完整基因维度的对象 subset，不从已裁剪的子集重跑。

## 7. 单一 HANDOFF，段级 upsert

一个项目只维护一份 `<工作目录>/Document/HANDOFF.md`；每个阶段跑完追加/替换自己那一段（默认替换，用户可显式要求追加），不动其他阶段段落。多数据集项目由 `DATASET_NAME` 决定段落标题后缀。字段契约见 `references/handoff-contract.md`。

## 8. 脚本编号与风格

- 两位数字前缀（`NN_名称.sh` / `.py` / `.R`），本仓库占用 `00`–`09`。
- 脚本头用框说明块，写清输入、输出、依赖、失败即停条件。
- 决策段用"**粗体关键词**：解释 — 理由"三段式，方便后来者判断能不能改。
- 文档中文为主，命令、路径、术语保留英文原文。

## 9. 提交前检查

```bash
bash -n scripts/*.sh
python -m compileall -q scripts
git grep -nIE 'C:/Users/|/home/[a-z]|K:/|sk-[A-Za-z0-9]{16,}|ghp_|github_pat_' -- . ':!.git' ':!.github' ':!CONTRIBUTING.md'
git grep -lI $'\r' -- . ':!.git' || true
```

前两条必须通过；第三条应无匹配（匹配即残留）；第四条应无输出（有输出即 CRLF 混入）。
第 3、4 条已在 CI 中自动执行（`.github/workflows/ci.yml` 的 `hygiene` 作业）。
