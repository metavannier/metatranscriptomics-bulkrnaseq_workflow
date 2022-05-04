RUN_ID = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())
ref_level = config["diffexp"]["ref_level"]

# ----------------------------------------------
# DESeq2: Differential expression analysis
# ----------------------------------------------

def get_deseq2_threads(wildcards=None):
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6

rule deseq2_init:
  input:
    cts = "05_Output/07_cpm/count_filtered.txt",
    coldata = config["coldata"]

  output:
    rds = "05_Output/08_deseq2_init/all.rds",
    normalized_counts_file = report("05_Output/08_deseq2_init/normalized_counts.tsv", caption="../report/normalized_counts.rst", category="02 Count matrices")
  conda:
    "../02_Container/deseq2.yaml"

  params:
    project = expand("{samples.project}",samples=samples.itertuples()),
    samples = expand("{samples.sample}",samples=samples.itertuples()),
    ref_level = config["diffexp"]["ref_level"],
    rmproj_list = config["filtering"]["rmproj"]

  threads: get_deseq2_threads()

  script:
    "../03_Script/deseq2.R"


rule diffexp:
  input:
    rds = "05_Output/08_deseq2_init/all.rds",
    rmd = "03_Script/diffexp.Rmd",
    pca = "03_Script/data_quality.R",
    coldata = config["coldata"]
  output:
    html_report=report(OUTPUTDIR + "09_differential_expression/diffexp.html", caption="../report/diffexp.rst", category="03 Report differential expression"),
    table=report(expand(OUTPUTDIR + "09_differential_expression/{condition.condition}_vs_{ref_level}_all_genes_stats.tsv", condition=condition.itertuples(), ref_level=ref_level), caption="../report/stat.rst", category="04 Table differential expression"),
    sur=report(expand(OUTPUTDIR + "09_differential_expression/{condition.condition}_vs_{ref_level}_signif-up-regulated.txt", condition=condition.itertuples(), ref_level=ref_level), caption="../report/stat.rst", category="04 Table differential expression"),
    sdr=report(expand(OUTPUTDIR + "09_differential_expression/{condition.condition}_vs_{ref_level}_signif-down-regulated.txt", condition=condition.itertuples(), ref_level=ref_level), caption="../report/stat.rst", category="04 Table differential expression")
  conda:
    CONTAINER + "diffexp.yaml"
  params:
    pca_labels=config["pca"]["labels"],
    ref_level = config["diffexp"]["ref_level"],
    lfcshrink_type = config["diffexp"]["lfcshrink_type"],
    gene_name = config["diffexp"]["gene_name"],
    mutant_level = config["diffexp"]["mutant_level"],
    nbpval = config["diffexp"]["nbpval"],
    FCcutoff=config["diffexp"]["FCcutoff"],
    pCutoff=config["diffexp"]["pCutoff"],
  message: 
    "Run the differential expression analyses"  
  script:
    SCRIPTDIR + "diffexp_reports_compilation.R"
