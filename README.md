[![GitHub Actions CI Status](https://github.com/ebi-metagenomics/shallowmapping/workflows/nf-core%20CI/badge.svg)](https://github.com/ebi-metagenomics/shallowmapping/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/ebi-metagenomics/shallowmapping/workflows/nf-core%20linting/badge.svg)](https://github.com/ebi-metagenomics/shallowmapping/actions?query=workflow%3A%22nf-core+linting%22)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/ebi-metagenomics/shallowmapping)

## Introduction

**ebi-metagenomics/shallowmapping** is a bioinformatics pipeline that generates taxonomic and functional profiles for low-yield (shallow shotgun: < 10 M reads) short raw-reads using [`MGnify biome-specific genome catalogues`](https://www.ebi.ac.uk/metagenomics/browse/genomes) as a reference. 

At the moment, the biome selection is limited to the precomputed databases available to downloading (chicken-gut-v1-0-1 and mouse-gut-v1-0). Other databases can be build for any of the [`MGnify genome catalogues`](https://www.ebi.ac.uk/metagenomics/browse/genomes) under request by opening an issue in this repo.

The main sections of the pipeline includes the following steps:
1. Raw-reads quality control ([`fastp`](https://github.com/OpenGene/fastp))
2. HQ reads decontamination versus human, phyX, and host ([`bwa-mem2`](https://github.com/bwa-mem2/bwa-mem2))
3. QC report of decontaminated reads ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. Integrated quality report of reads before and after decontamination ([`MultiQC`](http://multiqc.info/))
5. Mapping HQ clean reads using [`Sourmash`](https://github.com/sourmash-bio/sourmash) and bwa-mem2 (optional)
6. Taxonomic profiles generation
7. Functional profiles inference

The final output includes a species relative abundance table, Pfam and KEGG Orthologs (KO) count tables, a KEGG modules completeness table, and DRAM-style visuals. In addition, the shallow-mapping pipeline will integrate the taxonomic and functional tables of all the samples in the input samplesheet.

<p align="center" width="100%">
   <img src="docs/images/workflow.png" width="90%"/>
</p>


## Install and databases setup

### Required reference databases for decontamination

The first time you run the pipeline you need to put available the genome reference of human+phiX and the host (chicken or mouse) indexed for bwa-mem2. Already indexed reference genomes are hosted in the MGnify FTP site. To download, clone the repo and run the following commands:

```bash
cd shallowmapping
mkdir -p databases/reference_genomes && cd databases/reference_genomes

# Downloading human+phiX reference genomes

wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/references/human_phiX/human_phix_ref_bwamem2.tar.gz
tar -xvf human_phix_ref_bwamem2.tar.gz
mv bwamem2/* .
rm -r bwamem2

# Downloading the host genome. Replace $HOST by the name of the reference genome you are intend to download (e.g. chocken)

wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/references/$HOST/$HOST_ref_bwamem2.tar.gz
tar -xvf $HOST_ref_bwamem2.tar.gz
mv bwamem2/* .
rm -r bwamem2
```

You can move the location of the databases and update the config files accordingly.


#### MGnify genome catalogue databases
In addition, if this is the first time you are running the pipeline with a specific biome, you need to download the precomputed databases from the corresponding MGnify genomes catalogue.

Inside the databases create a directory named exactly as the Catalogue ID:
`chicken-gut-v1-0-1`
`mouse-gut-v1-0`

```bash
cd shallowmapping/databases/
mkdir $CATALOGUE_ID && cd $CATALOGUE_ID

# Downloading the catalogue metadata file. Replace $HOST for the name of the catalogue and $VERSION for the version you are intend to use as a reference

wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/$HOST/$VERSION/genomes-all_metadata.tsv 

# Downloading the pangenome function tables

wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/$HOST/$VERSION/pangenome_functions/functional_profiles.tar.gz
tar -xvf functional_profiles.tar.gz

wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/$HOST/$VERSION/pangenome_functions/kegg_completeness.tar.gz
tar -xvf kegg_completeness.tar.gz

# Downloading the representative genomes indexed for sourmash

wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/$HOST/$VERSION/sourmash_db_$HOST_$VERSION/sourmash_species_representatives_k51.sbt.zip

```

Running the pipeline using bwamem2 is optional. If you want to run the pipeline with this option you need to download the bwamem2 database for representative genomes as well:

```bash
wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/references/$BIOME_reps/$BIOME-$VERSION_bwamem2.tar.gz
tar -cvf $BIOME-$VERSION_bwamem2.tar.gz
mv $BIOME-$VERSION_bwamem2/* .
rm -r $BIOME-$VERSION_bwamem2
```

#### External databases for DRAM

```bash
cd shallowmapping/databases/
mkdir -p external_dbs/dram_distill_dbs && cd external_dbs/dram_distill_dbs
wget https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/amg_database.tsv
wget https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/etc_module_database.tsv
wget https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/function_heatmap_form.tsv
wget https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/genome_summary_form.tsv
wget https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/module_step_form.tsv
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
```


### Usage

Prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
paired_sample,/PATH/test_R1.fq.gz,/PATH/test_R2.fq.gz
single_sample,/PATH/test.fq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

Now, you can run the pipeline using the minumum mandatory arguments:

```bash
nextflow run /PATH/shallowmapping/main.nf \
   --biome <chicken-gut-v1-0-1 or mouse-gut-v1-0> \
   --input samplesheet.csv \
   --outdir <PROJECT_NAME> default = `results`
```

At the moment, the biome selection is limited to the precomputed databases available to downloading (chicken-gut-v1-0-1 and mouse-gut-v1-0). Other databases can be build for any of the [`MGnify genome catalogues`](https://www.ebi.ac.uk/metagenomics/browse/genomes) under request by opening an issue in this repo.


Optional arguments includes:
```bash
   --run_bwa <true or false> default = `false`    # To generate results using bwamem2 besides sourmash
   --core_mode <true or false> default = `false`  # To use core functions instead of pangenome functions
```

Use `--core_mode true` for large catalogues like the mouse-gut to avoid over-prediction due to an extremely large number of accessory genes in the pangenome.
Nextflow option `-profile` can be use to select a suitable config for your computational resources.
Nextflow option `-resume` can be use to re-run the pipeline from the last successfully finished step. 


## Credits

ebi-metagenomics/shallowmapping pipeline was originally written by @Ales-ibt.

We thank the following people for their extensive assistance in the development of this pipeline:

@mberacochea


## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

If you use  ebi-metagenomics/shallowmapping for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
