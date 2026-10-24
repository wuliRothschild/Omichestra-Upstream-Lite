# Omichestra Upstream (Lite)

A reproducible upstream pipeline for droplet-based single-cell RNA-seq: **FASTQ → alignment/quantification → RNA velocity (velocyto) → ambient-RNA removal (CellBender) → automated QC (scAutoQC) → analysis-ready `h5ad`**.

This is the "Lite" edition: the upstream module with a **single QC path** (nf-scautoqc + threshold gate + scrublet/scDblFinder + cluster-level BH). Platform-specific extensions are intentionally out of scope, so the codebase stays small and auditable.

## What it does

| Stage | Tool | Output |
|---|---|---|
| 1. Alignment / quantification | `cellranger count` (10x) or `mobivision quantify` (MobiNova) | filtered/raw matrices + BAM |
| 2. Splicing | `velocyto` | `<sid>.loom` (spliced/unspliced/ambiguous) |
| 3. Ambient RNA | `cellbender remove-background` | background-removed matrix |
| 4. QC | nf-scautoqc (`qc_mode=original`, GMM thresholds), threshold gate, scrublet + scDblFinder (default + homotypic), cluster-level BH | `scautoqc_pooled_filtered.h5ad` |
| 5. Delivery check | `09b_check_delivery.sh` | contract verification |
| 6. Archive & cleanup | `09c_cleanup.sh` | verify-then-delete local work tree |

Key properties:

- **Two alignment backends** behind one config switch (`ALIGNER=cellranger|mobivision`).
- **Strict serial scheduling** for alignment and CellBender (GPU exclusive), sample-isolated staging.
- **Triple consistency gate** (matrix / loom / h5) before CellBender — no silent coordinate drift.
- **Two independent doublet detectors** plus a hard threshold gate and a cluster-level BH test.
- **Fingerprint-based idempotency**: reruns skip only when inputs, parameters, references and script versions are unchanged.
- **Verify-then-delete cleanup** so the working tree never lingers after delivery.

## Requirements

- Linux (native filesystem for work directories; WSL2 works if work dirs are on ext4 — `/mnt/*` is rejected)
- conda / mamba
- `cellranger` (10x) and/or `mobivision` (MobiNova), `velocyto`
- CUDA GPU for CellBender
- R with `scDblFinder` + `schard` (single-cell doublet detection)
- Python 3.10+ with scanpy, anndata, h5py, scipy

See `references/environment-setup.md` and `references/reference-setup.md`.

## Quick start

```bash
git clone https://github.com/wuliRothschild/Omichestra-Upstream-Lite.git && cd Omichestra-Upstream-Lite

# 1) per project: config + sample sheet (the 00_* scripts read run_config.sh)
cp scripts/run_config.template.sh scripts/run_config.sh
$EDITOR scripts/run_config.sh
mkdir -p Document
cp assets/sample_sheet_template.csv Document/sample_sheet.csv
$EDITOR Document/sample_sheet.csv

# 2) once: environments + references
bash scripts/00_setup_envs.sh
bash scripts/00_setup_references.sh

# 3) stages (confirm parameters between stages)
bash scripts/01_cellranger_count.sh <sid> <fastq_dir> <prefix> <bam_archive> <cell_or_nuclei>
bash scripts/02_velocyto.sh <sid> <samplefolder> <bam_archive>
bash scripts/03_cellbender.sh
bash scripts/06_scautoqc_pipeline.sh

# 4) close out
bash scripts/09b_check_delivery.sh
bash scripts/09c_cleanup.sh --deliver-cr <delivery_root>
```

The skill-style instructions (`SKILL.md`) describe the same workflow for agent-driven execution, including the decision log and pitfalls.

## Configuration

`scripts/run_config.template.sh` is the single source of run-level settings (project name, species, aligner, reference, sample sheet, environments, output roots, optional overrides). Per-sample metadata goes to `outs/metadata.json` via `01c_write_metadata.py`:

```bash
"$SCANPY_PYTHON" scripts/01c_write_metadata.py --cr-root "$OUTPUT_ROOT/$PROJECT_NAME" \
  --samples <sid1,sid2,...> --meta-csv <metadata.csv> \
  --tech <value> --chemistry <value> --tissue <value> --species <human|mouse>
```

## Repository layout

```
SKILL.md                     # workflow + decision log (agent-facing)
README.md                    # this file
LICENSE                      # Apache-2.0
scripts/                     # numbered stages + helpers
references/                  # platform manuals, QC details, troubleshooting
assets/                      # sample sheet template + bundled reference gene lists (human/mouse, ~14 MB)
```

## Quality gates

- `02b_format_matrices.py` normalizes matrices to the standard 10x quartet.
- `03b_validate_triple_consistency.py` checks gene/barcode coordinate consistency before CellBender.
- `06d_threshold_gate.py` applies hard QC thresholds (counts/genes/mito) before doublet calling.
- `06e_scdblfinder_cluster_qc.py` adds a cluster-level BH axis to scDblFinder scores.
- `09b_check_delivery.sh` verifies the full delivery contract; missing items fail the run.
- `09c_cleanup.sh` verifies every local file exists in the delivery copy before deleting the work tree.

## License

Apache License 2.0 — see [LICENSE](LICENSE).

## Acknowledgements

This pipeline orchestrates and builds on excellent open-source work: [nf-scautoqc](https://github.com/cellgeni/nf-scautoqc), [scanpy](https://scanpy.readthedocs.io/), [anndata](https://anndata.readthedocs.io/), [velocyto](http://velocyto.org/), [CellBender](https://github.com/broadinstitute/CellBender), [scrublet](https://github.com/swolock/scrublet), [scDblFinder](https://bioconductor.org/packages/scDblFinder/), [STARsolo](https://github.com/alexdobin/STAR), [Cell Ranger](https://www.10xgenomics.com/), and MobiVision. Please cite the respective tools when you use them.

---

## 中文简介

本仓库是单细胞 RNA-seq 上游流水线的精简版：fastq → 比对定量（cellranger / MobiVision）→ velocyto 剪接 → CellBender 去背景 → nf-scautoqc 自动质控（GMM 阈值 + 硬阈值闸 + scrublet/scDblFinder 双线 + 簇级 BH + 池化），产出可直接用于下游整合与注释的 `scautoqc_pooled_filtered.h5ad`。执行手册与决策记录见 `SKILL.md`，平台细节与排错见 `references/`。许可证 Apache-2.0。
