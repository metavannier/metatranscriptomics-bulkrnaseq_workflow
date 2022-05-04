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
ref_level = config["diffexp"]["ref_level"]
genome = config["ref"]["genome"]
index = config["ref"]["index"]
annotation = config["ref"]["annotation"]
microbeannotatordb = config["microbeannotator"]["folder"]
annotationfolder = config["microbeannotator"]["annotationfolder"]
RUN_ID = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())

rule all:
  input:
    expand( "05_Output/02_trimmomatic/{samples}_{run}.trimmed.fastq", samples=SAMPLES, run=RUN),
    expand( "05_Output/02_trimmomatic/{samples}_{run}un.trimmed.fastq", samples=SAMPLES, run=RUN),
    # expand( "05_Output/02_sortmerna/{samples}_non-rRNA-reads_{way}.fq", samples=SAMPLES, way=WAY),
    # expand( "05_Output/02_sortmerna/{samples}_rRNA-reads_{way}.fq", samples=SAMPLES, way=WAY),
    # expand( "05_Output/03_fastqc/{samples}_non-rRNA-reads_{way}_fastqc.html", samples=SAMPLES, way=WAY),
    # expand( "05_Output/03_fastqc/{samples}_non-rRNA-reads_{way}_fastqc.zip", samples=SAMPLES, way=WAY),
    # OUTPUTDIR + "03_fastqc/non-rRNA-reads_trimmed_multiqc.html",
    #### Functional annotation with microbeannotator
    ### Make the database (To do on your local machine with aspera connect)
    # microbeannotatordb = microbeannotatordb + "/microbeannotator.db",
    # conversiondb = microbeannotatordb + "/conversion.db",
    ### Annotation
    # mockannotationfolder = annotationfolder + "/mockfile.txt",
    # index1 = expand( OUTPUTDIR + "{index}.1.ht2", index=index),
    # index2 = expand( OUTPUTDIR + "{index}.2.ht2", index=index),
    # index3 = expand( OUTPUTDIR + "{index}.3.ht2", index=index),
    # index4 = expand( OUTPUTDIR + "{index}.4.ht2", index=index),
    # index5 = expand( OUTPUTDIR + "{index}.5.ht2", index=index),
    # index6 = expand( OUTPUTDIR + "{index}.6.ht2", index=index),
    # index7 = expand( OUTPUTDIR + "{index}.7.ht2", index=index),
    # index8 = expand( OUTPUTDIR + "{index}.8.ht2", index=index),
    # bam = expand( OUTPUTDIR + "05_hisat/{samples}.bam", samples=SAMPLES),
    # countmatrices = expand( OUTPUTDIR + "06_featurecounts/{samples}_count.txt", samples=SAMPLES),
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

ruleorder: trimmomatic > fastqc_trimmed > hisat_build > hisat > featureCounts > cpm_filtering > deseq2_init > diffexp

# ----------------------------------------------
# Load rules 
# ----------------------------------------------

include: ENVDIR + "clean.smk"
include: ENVDIR + "functional_annotation.smk"
include: ENVDIR + "count.smk"
include: ENVDIR + "diffexp.smk"
