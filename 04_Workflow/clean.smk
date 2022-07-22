SAMPLES = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())
RUN =  config["run"]["type"].split(',')
WAY = config["run"]["way"].split(',')
EXT = config["run"]["ext"]

# ----------------------------------------------
# FastQC to check the reads quality
# ----------------------------------------------
# rule fastqc:
#   output:
#     expand( "05_Output/01_fastqc/{samples}_{run}_fastqc.html", samples=SAMPLES, run=RUN),
#     expand( "05_Output/01_fastqc/{samples}_{run}_fastqc.zip", samples=SAMPLES, run=RUN)

#   input:
#     expand( "00_RawData/{samples}_{run}.{ext}", samples=SAMPLES, run=RUN, ext=EXT)

#   conda: 
#     CONTAINER + "fastqc.yaml"

#   shell:
#     "fastqc --outdir 05_Output/01_fastqc/ {input}"

# # ----------------------------------------------
# # Trimmomatic: trimming reads and removing adapter sequences
# # ----------------------------------------------

rule trimmomatic:
  input:
    sample=expand("00_RawData/{samples}_{run}.{ext}", samples=SAMPLES, run=RUN, ext=EXT)

  output:
    sample_trimmed=expand( "05_Output/01_trimmomatic/{samples}_{run}.trimmed.fastq", samples=SAMPLES, run=RUN),
    sample_untrimmed=expand( "05_Output/01_trimmomatic/{samples}_{run}un.trimmed.fastq", samples=SAMPLES, run=RUN)

  params:
    adaptaters=config["trimmomatic"]["adaptaters"],
    leading=config["trimmomatic"]["leading"],
    trailing=config["trimmomatic"]["trailing"],
    minlen=config["trimmomatic"]["minlen"]

  conda: 
    CONTAINER + "trimmomatic.yaml"

  shell:
    """
    sample=({input.sample})
    adaptaters=({params.adaptaters})
    leading=({params.leading})
    trailing=({params.trailing})
    minlen=({params.minlen})
    sample_trimmed=({output.sample_trimmed})
    sample_untrimmed=({output.sample_untrimmed})
    len=${{#sample[@]}}
    for (( i=0; i<$len; i=i+2 ))
    do trimmomatic PE -threads 4 ${{sample[$i]}} ${{sample[$i+1]}} ${{sample_trimmed[$i]}} ${{sample_untrimmed[$i]}} ${{sample_trimmed[$i+1]}} ${{sample_untrimmed[$i+1]}} ILLUMINACLIP:${{adaptaters}}:2:30:10 LEADING:${{leading}} TRAILING:${{trailing}} SLIDINGWINDOW:4:15 MINLEN:${{minlen}}
    rm ${{sample[$i]}} ${{sample[$i+1]}}
    done
    """

# ----------------------------------------------
# SortmeRNA: Remove RNA
# ----------------------------------------------

rule sortmerna:
  input:
    sample_trimmed=expand( "05_Output/01_trimmomatic/{samples}_{run}.trimmed.fastq", samples=SAMPLES, run=RUN),
    ddb = config["ref"]["ddb"],
  output:
    nonrrna=expand( "05_Output/02_sortmerna/{samples}_non-rRNA-reads_{way}.fq", samples=SAMPLES, way=WAY),
    rrna=expand( "05_Output/02_sortmerna/{samples}_rRNA-reads_{way}.fq", samples=SAMPLES, way=WAY),
  params:
    workdir = "05_Output/02_sortmerna",
    rmkvdb = "05_Output/02_sortmerna/kvdb",
    rmreadb = "05_Output/02_sortmerna/readb",
    nonrrnapath = expand( "05_Output/02_sortmerna/{samples}_non-rRNA-reads", samples=SAMPLES),
    rrnapath = expand( "05_Output/02_sortmerna/{samples}_rRNA-reads", samples=SAMPLES)
  conda:
    CONTAINER + "sortmerna.yaml"
  shell:
    """
    sample_trimmed=({input.sample_trimmed})
    nonrrnapath=({params.nonrrnapath})
    rrnapath=({params.rrnapath})
    len=${{#sample_trimmed[@]}}
    flag=0
    for (( i=0; i<$len; i=i+2 ))
    do sortmerna --reads ${{sample_trimmed[$i]}} --reads ${{sample_trimmed[$i+1]}} --workdir {params.workdir} --ref {input.ddb} --fastx --paired_out --out2 --aligned ${{rrnapath[$flag]}} --other ${{nonrrnapath[$flag]}} -threads 6
    flag=$((flag+1))
    rm -r {params.rmkvdb}
    rm -r {params.rmreadb}
    rm ${{sample_trimmed[$i]}} ${{sample_trimmed[$i+1]}}
    done
    """

# ----------------------------------------------
# FastQC to check the reads trimmed quality
# ----------------------------------------------
rule fastqc_trimmed:
  output:
    expand( "05_Output/03_fastqc/{samples}_non-rRNA-reads_{way}_fastqc.html", samples=SAMPLES, way=WAY),
    expand( "05_Output/03_fastqc/{samples}_non-rRNA-reads_{way}_fastqc.zip", samples=SAMPLES, way=WAY)

  input:
    nonrrna=expand( "05_Output/02_sortmerna/{samples}_non-rRNA-reads_{way}.fq", samples=SAMPLES, way=WAY),

  conda: 
    CONTAINER + "fastqc.yaml"

  shell:
    "fastqc --outdir 05_Output/03_fastqc/ {input}"

# ----------------------------------------------
# MultiQC to check the reads trimmed quality
# ----------------------------------------------

rule multiqc_trimmed:
  input:
    trim_qc = expand( "05_Output/03_fastqc/{samples}_non-rRNA-reads_{way}_fastqc.zip", samples=SAMPLES, way=WAY)
  output:
    trim_multi_html = report(OUTPUTDIR + "03_fastqc/non-rRNA-reads_trimmed_multiqc.html", caption = ROOTDIR + "/07_Report/multiqc.rst", category="01 Quality report"), 
  params:
    multiqc_output_trim = OUTPUTDIR + "03_fastqc/non-rRNA-reads_trimmed_multiqc_data"
  conda:
    CONTAINER + "multiqc.yaml"
  shell: 
    """
    multiqc -n {output.trim_multi_html} {input.trim_qc} --force #run multiqc
    rm -rf {params.multiqc_output_trim} #clean-up
    """
