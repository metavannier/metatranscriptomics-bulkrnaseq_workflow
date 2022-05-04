# Differential Expression Workflow: RNA-seq analysis

## Author

Thomas Vannier (@metavannier), https://centuri-livingsystems.org/t-vannier/

## About

This workflow performs an RNA-seq analysis from the sequencing output data to the differential expression analyses.

You need to install [Singularity](https://github.com/hpcng/singularity/blob/master/INSTALL.md#install-golang) on your computer. This workflow also work in a slurm environment.

Each snakemake rules call a specific conda environment. In this way you can easily change/add tools for each step if necessary. 

3 steps for the analysis:
- clean.smk: The quality of the raw reads are assessed using [FastQC v0.11.9 toolkit](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Adapters and low quality reads are trimmed using [Trimmomatic v0.39](https://academic.oup.com/bioinformatics/article/30/15/2114/2390096).
- count.smk: [HiSat2 v2.2.1](https://www.nature.com/articles/nmeth.3317) is used for mapping the raw reads to the reference genome. The expression for each gene is evaluated using featureCounts from the [Subread v2.0.1 package](https://pubmed.ncbi.nlm.nih.gov/30783653/).
- differential_exp.smk: The low expressed genes are removed from further analysis. The raw counts are normalized and used for differential expression testing using [DESeq2 v1.28.0](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8).

## Usage

### Step 1: Install workflow

You can use this workflow by downloading and extracting the latest release. If you intend to modify and further extend this workflow or want to work under version control, you can fork this repository.

We would be pleased if you use this workflow and participate in its improvement. If you use it in a paper, don't forget to give credits to the author by citing the URL of this repository and, if available, its DOI (see above).

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files and repositories:
- 00_RawData need the single or pair-end fastq file of each run to analyse
- 01_Reference the fasta file and gff/gtf of your reference genome for the mapping step
- [sample.tsv](/sample.tsv), [coldata.tsv](/coldata.tsv) and [condition.tsv](/condition.tsv) to indicate the samples, run, condition, etc. for the analyse.
- [config.yaml](/config.yaml) indicating the parameters to use.

### Step 3: Execute workflow

- You need [Singularity v3.5.3](https://github.com/hpcng/singularity/blob/master/INSTALL.md#install-golang) installed on your computer or cluster.

- Load snakemake from a docker container and run the workflow from the root by using these commands:

`singularity run docker://snakemake/snakemake:v6.3.0`

- Then execute the workflow locally via

`snakemake  --use-conda --use-singularity --cores 10`

### Step 4: Investigate results 

After successful execution, you can create a self-contained interactive HTML report with all results via:

`snakemake --report report.html`
