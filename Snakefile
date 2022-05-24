# Docker container based on a minimal Ubuntu installation that includes conda-forge's mambaforge installer.
container: "docker://condaforge/mambaforge"

import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")

##### load config and sample sheets #####

configfile: "config.yaml"
validate(config, schema="06_Schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index(["project", "condition", "sample"], drop=False)
validate(samples, schema="06_Schemas/samples.schema.yaml")

coldata = pd.read_table(config["coldata"]).set_index(["project", "condition", "type"], drop=False)
validate(coldata, schema="06_Schemas/coldata.schema.yaml")

condition = pd.read_table(config["condition"]).set_index(["condition"], drop=False)
validate(coldata, schema="06_Schemas/condition.schema.yaml")

##### Set variables ####
ROOTDIR = os.getcwd()
RAWDATA = srcdir("00_RawData/")
REF = srcdir("01_Reference/")
CONTAINER = srcdir("02_Container/")
SCRIPTDIR = srcdir("03_Script/")
ENVDIR = srcdir("04_Workflow/")
OUTPUTDIR = srcdir("05_Output/")
REPORT = srcdir("07_Report/")

# ----------------------------------------------
# Target rules
# ----------------------------------------------
SAMPLES = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())
RUN =  config["run"]["type"].split(',')
EXT = config["run"]["ext"]
WAY = config["run"]["way"].split(',')
PROJ = config["qiime"]["name"]
# Database information to assign taxonomy
DB_classifier = config["taxonomy"]["database_classified"]
ref_level = config["diffexp"]["ref_level"]
genome = config["ref"]["genome"]
index = config["ref"]["index"]
annotation = config["ref"]["annotation"]
microbeannotatordb = config["microbeannotator"]["folder"]
annotationfolder = config["microbeannotator"]["annotationfolder"]
RUN_ID = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())

rule all:
  input:
    # expand( "05_Output/01_fastqc/{samples}_{run}_fastqc.html", samples=SAMPLES, run=RUN),
    # expand( "05_Output/01_fastqc/{samples}_{run}_fastqc.zip", samples=SAMPLES, run=RUN),
    # sample_trimmed=expand( "05_Output/01_trimmomatic/{samples}_{run}.trimmed.fastq", samples=SAMPLES, run=RUN),
    # sample_untrimmed=expand( "05_Output/01_trimmomatic/{samples}_{run}un.trimmed.fastq", samples=SAMPLES, run=RUN),
    # nonrrna=expand( "05_Output/02_sortmerna/{samples}_non-rRNA-reads_{way}.fq", samples=SAMPLES, way=WAY),
    # rrna=expand( "05_Output/02_sortmerna/{samples}_rRNA-reads_{way}.fq", samples=SAMPLES, way=WAY),
    # expand( "05_Output/03_fastqc/{samples}_non-rRNA-reads_{way}_fastqc.html", samples=SAMPLES, way=WAY),
    # expand( "05_Output/03_fastqc/{samples}_non-rRNA-reads_{way}_fastqc.zip", samples=SAMPLES, way=WAY),
    # OUTPUTDIR + "03_fastqc/non-rRNA-reads_trimmed_multiqc.html",
    #### Functional annotation with microbeannotator
    ### Make the database (To do on your local machine with aspera connect)
    # microbeannotatordb = microbeannotatordb + "/microbeannotator.db",
    # conversiondb = microbeannotatordb + "/conversion.db",
    ### Annotation
    # mockannotationfolder = annotationfolder + "/mockfile.txt",
    #### Qiime import
    # q2_import = OUTPUTDIR + "05_qiime_import/" + PROJ + "-demux-paired-end.qza",
    # primer = OUTPUTDIR + "05_qiime_import/" + PROJ + "-PE-demux-noprimer.qzv",
    # table = OUTPUTDIR + "06_dada2/" + PROJ + "-table-dada2.qza",
    # rep = OUTPUTDIR + "06_dada2/" + PROJ + "-rep-seqs-dada2.qza",
    # stats = OUTPUTDIR + "06_dada2/" + PROJ + "-dada2-stats.qza",
    # filterrep = OUTPUTDIR + "06_dada2/" + PROJ + "-rep-filtered-seqs-dada2.qza",
    # filtertable = OUTPUTDIR + "06_dada2/" + PROJ + "-table-filtered-dada2.qza",
    # sklearn = OUTPUTDIR + "07_taxonomy/" + PROJ + "-tax_sklearn.qza",
    # taxafiltertable = OUTPUTDIR + "07_taxonomy/" + PROJ + "-taxa-table-filtered-dada2.qza",
    # q2_repseq_filtered = OUTPUTDIR + "07_taxonomy/" + PROJ + "-rep-filtered-seqs-taxa-dada2.qza",
    # table_tax = OUTPUTDIR + "07_taxonomy/taxonomy.tsv",
    # table_tax_filtered = report(OUTPUTDIR + "07_taxonomy/taxonomy_filtered.tsv", caption = ROOTDIR + "/07_Report/tax.rst", category="04 taxonomy"),
    # rep_viz = report(OUTPUTDIR + "07_taxonomy/" + PROJ + "-rep-filtered-seqs-taxa-dada2.qzv", caption = ROOTDIR + "/07_Report/dada2seq.rst", category="03 dada2"),
    # stats_viz = report(OUTPUTDIR + "06_dada2/" + PROJ + "-dada2-stats.qzv", caption = ROOTDIR + "/07_Report/dada2summary.rst", category="03 dada2"),
    # featurestat = report(OUTPUTDIR + "07_taxonomy/" + PROJ + "-taxa-table-filtered-dada2.qzv", caption = ROOTDIR + "/07_Report/dada2summary.rst", category="04 taxonomy"),
    # rarefactionfiltertable = OUTPUTDIR + "07_taxonomy/" + PROJ + "-rarefaction-table-filtered-dada2.qza",
    # relativefreqtable = OUTPUTDIR + "07_taxonomy/" + PROJ + "-relative-frequency-dada2.qza",
    # table_biom = OUTPUTDIR + "07_taxonomy/feature-table.biom",
    # taxo_table_biom = OUTPUTDIR + "07_taxonomy/" + PROJ + "-asv-table-with-taxonomy.biom",
    # taxo_table_tsv = OUTPUTDIR + "07_taxonomy/" + PROJ + "-asv-table-with-taxonomy.tsv",
    # taxabarplots = OUTPUTDIR + "07_taxonomy/" + PROJ + "-taxa-bar-plots.qzv",
    # index1 = expand( OUTPUTDIR + "{index}.1.ht2", index=index),
    # index2 = expand( OUTPUTDIR + "{index}.2.ht2", index=index),
    # index3 = expand( OUTPUTDIR + "{index}.3.ht2", index=index),
    # index4 = expand( OUTPUTDIR + "{index}.4.ht2", index=index),
    # index5 = expand( OUTPUTDIR + "{index}.5.ht2", index=index),
    # index6 = expand( OUTPUTDIR + "{index}.6.ht2", index=index),
    # index7 = expand( OUTPUTDIR + "{index}.7.ht2", index=index),
    # index8 = expand( OUTPUTDIR + "{index}.8.ht2", index=index),
    # bam = expand( OUTPUTDIR + "08_hisat/{samples}.bam", samples=SAMPLES),
    #### Samtools coverage: genome coverage
    # coverage = expand( OUTPUTDIR + "08_hisat/{samples}_coverage.txt", samples=SAMPLES),
    # avcoverage = expand( OUTPUTDIR + "08_hisat/average_coverage.txt", samples=SAMPLES),
    countmatrices = expand( OUTPUTDIR + "09_featurecounts/{samples}_count.txt", samples=SAMPLES),
    # count_df = OUTPUTDIR + "07_cpm/count.txt",
    # output_filter_count = OUTPUTDIR + "07_cpm/count_filtered.txt",
    # cpm = OUTPUTDIR + "07_cpm/cpm_filtered.txt",
    # rds = "05_Output/08_deseq2_init/all.rds",
    # normalized_counts_file = "05_Output/08_deseq2_init/normalized_counts.tsv",
    # table=expand(OUTPUTDIR + "09_differential_expression/{condition.condition}_vs_{ref_level}_all_genes_stats.tsv", condition=condition.itertuples(), ref_level=ref_level),
    # sur=expand(OUTPUTDIR + "09_differential_expression/{condition.condition}_vs_{ref_level}_signif-up-regulated.txt", condition=condition.itertuples(), ref_level=ref_level),
    # sdr=expand(OUTPUTDIR + "09_differential_expression/{condition.condition}_vs_{ref_level}_signif-down-regulated.txt", condition=condition.itertuples(), ref_level=ref_level),
    # html_report = OUTPUTDIR + "09_differential_expression/diffexp.html",

# ----------------------------------------------
# setup singularity 
# ----------------------------------------------

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
#singularity: "docker://continuumio/miniconda3"

# ----------------------------------------------
# setup report
# ----------------------------------------------

report: "report/workflow.rst"

# ----------------------------------------------
# Impose rule order for the execution of the workflow 
# ----------------------------------------------

# ruleorder: trimmomatic > fastqc_trimmed > hisat_build > hisat > featureCounts > cpm_filtering > deseq2_init > diffexp

# ----------------------------------------------
# Load rules 
# ----------------------------------------------

include: ENVDIR + "clean.smk"
include: ENVDIR + "functional_annotation.smk"
include: ENVDIR + "taxonomic_annotation.smk"
include: ENVDIR + "count.smk"
include: ENVDIR + "diffexp.smk"
