# Test dataset

This test dataset consists of several components that have been manually curated for testing purposes, aggressively trimmed down (with the exception of the DRAM Distill db).

The dataset includes a trimmed version of the human-vaginal-v1.0 dataset alongside a reduced human_phiX BWAMem2 index. For the sourmash sketch index, we created a custom for the genome [MGYG000303700](https://www.ebi.ac.uk/metagenomics/api/v1/genomes/MGYG000303700/downloads/MGYG000303700.fna). Additionally, we generated simulated reads using [samtools wgsim](https://sarahpenir.github.io/bioinformatics/Simulating-Sequence-Reads-with-wgsim/).

## Database build test (`--new_mgnify_db`)

`build_mgnify_db.nf.test` uses two more components:

- `mgnify_catalogue/`: a tiny catalogue with the layout of a catalogue folder on the cluster, made by `mgnify_catalogue/make_fixture.py`. It has a singleton, a pangenome with consistent files, and a pangenome with the duplicate-name bug of old mmseqs catalogues. The accessory genes are SARS-CoV-2 genes, so that eggNOG-mapper finds hits.
- `eggnog_db/`: a tiny eggNOG database. `eggnog_proteins.dmnd` is the DIAMOND database of the SARS-CoV-2 proteome. `eggnog.db` is the nf-core test `eggnog.db` with annotation rows added for the S and N proteins. Both files come from [nf-core/test-datasets](https://github.com/nf-core/test-datasets/tree/modules/data) (`genomics/sarscov2/genome/` and `delete_me/eggnogmapper/eggnog.db`).

To run the database build on the test data directly:

```bash
nextflow run main.nf -profile test_build_db,docker --outdir results
```
