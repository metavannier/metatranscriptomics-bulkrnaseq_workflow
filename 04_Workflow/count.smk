SAMPLES = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())

# ----------------------------------------------
# HISAT2-build: Indexing a reference genome
# ----------------------------------------------
rule hisat_build:
  output:
    expand( OUTPUTDIR + "{index}.1.ht2", index=index),
    expand( OUTPUTDIR + "{index}.2.ht2", index=index),
    expand( OUTPUTDIR + "{index}.3.ht2", index=index),
    expand( OUTPUTDIR + "{index}.4.ht2", index=index),
    expand( OUTPUTDIR + "{index}.5.ht2", index=index),
    expand( OUTPUTDIR + "{index}.6.ht2", index=index),
    expand( OUTPUTDIR + "{index}.7.ht2", index=index),
    expand( OUTPUTDIR + "{index}.8.ht2", index=index)

  input:
    ref = REF + config["ref"]["genome"]

  params:
    index = OUTPUTDIR + config["ref"]["index"]

  conda: 
    CONTAINER + "hisat2.yaml"

  shell:
    """
    input={input.ref}
    output={params.index}
    if [ ${{input}} == '*.gz' ];then
      gunzip ${{input}}
      input="${{input%.*}}"
    fi
    hisat2-build ${{input}} ${{output}}
    """

# ----------------------------------------------
# HISAT2: alignment of NGS reads to a population of genomes
# ----------------------------------------------
rule hisat:
  output:
    bam = expand( OUTPUTDIR + "08_hisat/{samples}.bam", samples=SAMPLES)
    
  input:
    nonrrna = expand( "05_Output/02_sortmerna/{samples}_non-rRNA-reads_{way}.fq", samples=SAMPLES, way=WAY),

  conda: 
    CONTAINER + "hisat2.yaml"

  params:
    sam = expand( OUTPUTDIR + "08_hisat/{samples}.sam", samples=SAMPLES),
    reads = config["run"]["reads"],
    index = OUTPUTDIR + config["ref"]["index"]

  shell:
    """
    index={params.index}
    nonrrna=({input.nonrrna})
    len=${{#nonrrna[@]}}
    reads=({params.reads})
    sam=({params.sam})
    bam=({output.bam})
    flag=0
    if [ ${{reads}} == 'paired' ];then
      for (( i=0; i<$len; i=i+2 ))
        do hisat2 -p 12 -x ${{index}} -1 ${{nonrrna[$i]}} -2 ${{nonrrna[$i+1]}} -S ${{sam[$flag]}}
          samtools sort ${{sam[$flag]}} > ${{bam[$flag]}}
          rm ${{sam[$flag]}}
          flag=$((${{flag}}+1))
      done
    elif [ ${{reads}} == 'unpaired' ];then
      for (( i=0; i<$len; i++ ))
        do hisat2 -p 12 -x ${{index}} -U ${{nonrrna[$i]}} -S ${{sam[$i]}}
          samtools sort ${{sam[$i]}} > ${{bam[$i]}}
          rm ${{sam[$i]}}
      done
    else
      echo "Your fastq files have to be paired or unpaired. Please check the config.yaml"
    fi
    """

# ----------------------------------------------
# samtools coverage: genome coverage
# ----------------------------------------------

rule coverage:
  input:
    bam = expand( OUTPUTDIR + "08_hisat/{samples}.bam", samples=SAMPLES) 
  output:
    coverage = expand( OUTPUTDIR + "08_hisat/{samples}_coverage.txt", samples=SAMPLES) 
  conda:
    CONTAINER + "hisat2.yaml"
  shell:
    """
    bam=({input.bam})
    coverage=({output.coverage})
    len=${{#bam[@]}}
    for (( i=0; i<$len; i++ ))
      do samtools coverage ${{bam[$i]}} -o ${{coverage[$i]}}
    done
    """

rule average_coverage:
  input:
    coverage = expand( OUTPUTDIR + "08_hisat/{samples}_coverage.txt", samples=SAMPLES) 
  output:
    avcoverage = expand( OUTPUTDIR + "08_hisat/average_coverage.txt", samples=SAMPLES)
  container: None
  shell:
    """
      coverage_file=({input.coverage})
      avcoverage_file=({output.avcoverage})
      len=${{#coverage_file[@]}}
      for (( i=0; i<$len; i++ ))
          do file=${{coverage_file[$i]}}
          echo $file >> $avcoverage_file
          awk -F " " "{{SUM+=\$6}} END {{print SUM/(NR-1)}}" $file >> $avcoverage_file
      done
    """
# ----------------------------------------------
# featureCounts: read count/summarization program
# ----------------------------------------------

rule featureCounts:
  output:
    countmatrices = expand( OUTPUTDIR + "09_featurecounts/{samples}_count.txt", samples=SAMPLES)
    
  input:
    annotation = REF + config["ref"]["annotation"],
    bam = expand( OUTPUTDIR + "08_hisat/{samples}.bam", samples=SAMPLES)

  conda: 
    CONTAINER + "featureCounts.yaml"

  params:
    reads = config["run"]["reads"],
    geneid = config["ref"]["geneid"],
    geneidentifier = config["ref"]["geneidentifier"]

  shell:
    """
    reads=({params.reads})
    geneid=({params.geneid})
    annotation={input.annotation}
    bam=({input.bam})
    countmatrices=({output.countmatrices})
    geneidentifier=({params.geneidentifier})
    len=${{#bam[@]}}
    if [ ${{reads}} == 'paired' ];then
      for (( i=0; i<$len; i++ ))
        do featureCounts -T 12 -p -t ${{geneidentifier}} -g ${{geneid}} -a ${{annotation}} -o ${{countmatrices[$i]}} ${{bam[$i]}}
      done
    elif [ ${{reads}} == 'unpaired' ];then
      for (( i=0; i<$len; i++ ))
        do featureCounts -T 12 -t ${{geneidentifier}} -g ${{geneid}} -a ${{annotation}} -o ${{countmatrices[$i]}} ${{bam[$i]}}
      done
    fi
    """

# ----------------------------------------------
# cpm: Counts Per Million and filtering
# ----------------------------------------------

rule cpm_filtering:
  output:
    count_df = report(OUTPUTDIR + "10_cpm/count.txt", caption="../07_Report/count.rst", category="02 Count matrices"),
    output_filter_count = report(OUTPUTDIR + "10_cpm/count_filtered.txt", caption="../07_Report/count_filtered.rst", category="02 Count matrices"),
    cpm = report(OUTPUTDIR + "10_cpm/cpm_filtered.txt", caption="../07_Report/cpm_filtered.rst", category="02 Count matrices")

  input:
    expand( OUTPUTDIR + "09_featurecounts/{samples}_count.txt", samples=SAMPLES)
    
  params:
    thresh_cpm = config["filtering"]["thresh_cpm"],
    thresh_sample = config["filtering"]["thresh_sample"],
    rmrun_list = config["filtering"]["rmrun"],
    path = OUTPUTDIR + "09_featurecounts"

  conda: 
    CONTAINER + "cpm.yaml"

  script:
    SCRIPTDIR + "cpm_filtering.R"
