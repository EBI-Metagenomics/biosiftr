# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

BioSIFTR is a Nextflow (DSL2) / nf-core-style pipeline that produces taxonomic and functional
profiles for **shallow-shotgun** short reads (< ~10M reads) by mapping them against MGnify
**biome-specific genome catalogues**. It is containerised (Docker/Apptainer/Singularity); Conda is
not supported.

## Commands

```bash
# Run the full nf-test suite (uses the bundled trimmed test dataset under tests/)
nf-test test --profile test,docker          # CI runs: nf-test test --ci --profile test,docker

# Run a single test case
nf-test test tests/default.nf.test --profile test,docker

# Run the pipeline directly on the test data (alternative to nf-test, see tests/README.md)
cd tests && nextflow run ../main.nf -profile test,docker --outdir results

# Real run (databases auto-download on first use into --reference_dbs)
nextflow run main.nf -profile docker \
  --biome <CATALOGUE_ID> --input samplesheet.csv --outdir <DIR> \
  --reference_dbs </path/to/dbs> --decontamination_indexes </path/to/bwamem2/indexes>

# Lint python in bin/ (config in pyproject.toml) and format files (prettier, see .pre-commit-config.yaml)
ruff check bin/
pre-commit run --all-files
```

Snapshot tests assert against `tests/default.nf.test.snap`. When output legitimately changes,
regenerate snapshots with `nf-test test --update-snapshot ...`. `tests/.nftignore` lists paths
whose _content_ is volatile and excluded from content snapshots (only their names are checked).

## Architecture

The pipeline entrypoint chain is `main.nf` → `workflows/biosiftr.nf` (the `BIOSIFTR` workflow).
`workflows/biosiftr.nf` is the orchestration hub; almost all logic lives there.

Conceptual stages (all driven from `BIOSIFTR`):

1. **Preprocessing** (skippable via `--skip_decont`): `FASTP` trim → `BWAMEM2DECONTNOBAMS`
   decontamination against human+phiX, then a second pass against the **biome's canonical host**
   when the biome is not human. The host genome name is derived by `params.biome.split('-')[0]`
   (e.g. `cow-rumen-v1-0-1` → `cow`) and its index is expected under `--decontamination_indexes`.
   Then `FASTQC` on clean reads.
2. **Reference resolution**: `subworkflows/download_references.nf` (`DOWNLOAD_REFERENCES`) checks
   whether all required per-biome DB files already exist under `--reference_dbs/<biome>/`; if not it
   downloads them. Same logic for the human_phiX bwa-mem2 index and the DRAM dbs. This subworkflow
   is the single source of truth for DB layout (sourmash sbt, `genomes-all_metadata.tsv`,
   `functional_profiles_DB/`, `kegg_completeness_DB/`, optional `bwa_reps.*`).
3. **Mapping** — two parallel branches over the same clean reads:
   - **Sourmash** (always): `SOURMASH_SKETCH` → `SOURMASH_GATHER` → `POSTPROC_SOURMASHTAXO`.
   - **bwa-mem2** (only when `--run_bwa true`): `ALIGN_BWAMEM2` → `POSTPROC_BWATAXO`. The bwa
     coverage threshold adapts to species richness reported by sourmash (>150 species → 0.01, else 0.1).
4. **Functional inference** per branch: `POSTPROC_FUNCTIONSPRED` (pangenome vs core controlled by
   `--core_mode`), `KEGG_SPECIES`/`KEGG_COMPLETENESS`, and optional `DRAM_DISTILL` (`--run_dram`).
5. **Integration**: `POSTPROC_INTEGRATOR` aliases combine per-sample taxonomy / KO / Pfam / KEGG
   module tables into project-wide matrices. `MULTIQC` aggregates QC.

The sourmash and bwa branches reuse the **same local modules** via Nextflow `include ... as`
aliases (`SM_*` vs `BWA_*`, `INTEGRA_*` vs `BWA_INT_*`). When editing a module's interface, check
every alias in the include block at the top of `workflows/biosiftr.nf`.

### Module layout

- `modules/local/` — pipeline-specific processes. Each wraps a Python script in `bin/` (e.g.
  `postproc/functionspred.nf` runs `bin/species2functions.py`). Business logic lives in `bin/*.py`;
  the `.nf` files are thin wrappers. To change profiling behaviour, edit the corresponding `bin/` script.
- `modules/nf-core/` and `modules/ebi-metagenomics/` — vendored modules tracked in `modules.json`;
  treat as upstream (avoid hand-editing).
- `conf/modules.config` — per-process `publishDir` and `ext.args`; output directory structure is
  defined here, not in the modules.
- `lib/*.groovy` — nf-core template helpers (`WorkflowBiosiftr`, `WorkflowMain`, etc.).

### Key conventions

- Channel items are `[ meta, files ]` tuples where `meta = [id:..., single_end:bool]`, built by the
  `groupReads` closure from the samplesheet (`sample,fastq_1,fastq_2`).
- Params are declared in `nextflow.config` and validated against `nextflow_schema.json` via the
  `nf-validation` plugin. Add new params in **both** places.
- Notable flags: `--run_bwa`, `--core_mode` (use for large catalogues like human-gut to avoid
  accessory-gene over-prediction), `--run_dram`, `--skip_decont`.
- `-profile` selects compute env (`docker`, `singularity`, `local`, `ebi`, `test`, ...). `test`
  pulls config from `conf/test.config` pointing at the bundled `tests/reference_dbs/` and
  `tests/bwamem2/`.

Version is set in `manifest.version` in `nextflow.config`.
