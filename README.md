[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.04.2-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.3.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.3.1)

## Introduction

**nf-ontgeno** is a bioinformatics pipeline that processes Oxford Nanopore Technologies (ONT) sequencing data for variant detection and genotyping. The pipeline is built using [Nextflow](https://www.nextflow.io/) and follows the [nf-core](https://nf-co.re/) community guidelines to ensure high-quality, reproducible, and portable analyses.

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))
3. Map reads to reference genome with minimap2
4. Produce mapping statistics with SAMtools, Mosdepth, and PanDepth
4. Perform variant calling with Clair3 and DeepVariant
5. Joint genotype calling with GLNexus
6. Phasing with Eagle2 and Whatshap
7. Variant annotation with Ensembl VEP
8. Structural variant calling with Sniffles2 and cuteSV

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

### Samplesheet

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:
```csv
sample,runid,library,fastq
patient01,run1,A,/path_to_data/patient01.fastq.gz
```

Each row represents a fastq file. Files with the same sample id are merged during the analysis. The combination of sample, runid and library should be unique for each file.

### Parameter file

Next, prepare a YAML file for the pipeline parameters, referring to the samplesheet with the `input` parameter.

`params.yaml`:
```yaml
# General settings:
#-----------------------
project: "p2024-0032"
outdir: "results_sup"

# Sample details:
#-----------------------
input: "test/samplelist.sup.csv"

# References and intervals:
#-----------------------
fasta: "/data/references/Homo_sapiens/Ensembl/GRCh38/Sequence/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
bed: "test/region.cftr.bed"
panel: "assets/1kGP_high_coverage_Illumina.7_117470000_117680000.filtered.SNV_INDEL_SV_phased_panel.bcf"

# Filter settings:
#-----------------------
trimmer: 'fastp'
min_length: 100
trim_length: 0
store_trimmed: false
min_mapq: 20
min_coverage: 10

# Workflow settings:
#-----------------------
basecalling_model: 'sup'
genotype_model: "assets/r1041_e82_400bps_sup_v410"
glnexus_config: "assets/clair3.yml"

# Annotation settings:
#-----------------------
genome: 'GRCh38'
species: 'homo_sapiens'
vep_cache_version: 112
vep_cache: "/data/databases/vep"
```

### Running on the IBU cluster

To run the pipeline, start an interactive session on the IBU cluster:

```bash
srun --partition pibu_el8 --account p2024-0032 --cpus-per-task=1 --mem=8000 --time=144:00:00 --pty bash
module load Java
export APPTAINER_CACHEDIR=$SCRATCH
export NXF_SINGULARITY_CACHEDIR=/data/projects/p2024-0032_comprehensive_cftr_gene_sequencing_using_long_read_nanopore_technology/pipelines/singularity_cache
export NXF_TEMP=$SCRATCH
```

Now, you can run the pipeline using:

```bash
nextflow run main.nf \
  -profile unibe_ibu \
  -params-file test/params.yaml
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

### Local execution

To run the pipeline locally on your computer, use the provided `local.config` file and adjust settings as needed.

```bash
nextflow run main.nf \
  -c local.config \
  -params-file test/params.yaml
```

Make sure to set the environment variable `NXF_SINGULARITY_CACHEDIR` to avoid having to download Singularity containers repeatedly:
```bash
echo $NXF_SINGULARITY_CACHEDIR
```

If not defined, set it before running the pipeline to a directory with sufficient free disk space:
```bash
export NXF_SINGULARITY_CACHEDIR=/path/to/singularity_cache
```

### Clean-up

To completely clean up all previous pipeline runs, the following files and folders need to be deleted:

```bash
rm -r work            # The working cache of the pipeline
rm -r .nextflow       # The Nextflow database containing information for nextflow log
rm -r <outdir>        # The folder with published results as defined in --outdir
rm .nextflow.log*     # The log files of previous runs
```
After this, it is no longer possible to resume a previous run, so be careful what you delete. To just delete a specific pipeline run, use `nextflow log` to get the ids of all runs and remove the run with `nextflow clean <RUN_NAME> -f`.

## Pipeline output

All result files will be published to the directory defined by the `--outdir` parameter.

## Credits

nf-ontgeno was originally written by Alexander Nater.

We thank the following people for their extensive assistance in the development of this pipeline:
- Nada El Makhzen

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
