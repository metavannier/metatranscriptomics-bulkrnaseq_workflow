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
    cts = "05_Output/10_cpm/count_filtered.txt",
    coldata = config["coldata"]

  output:
    rds = "05_Output/11_deseq2_init/all.rds",
    normalized_counts_file = report("05_Output/11_deseq2_init/normalized_counts.tsv", caption="../07_Report/normalized_counts.rst", category="05 Count matrices")
  conda:
    "../02_Container/deseq2.yaml"

  params:
    project = expand("{samples.project}",samples=samples.itertuples()),
    samples = expand("{samples.sample}",samples=samples.itertuples()),
    condition = expand("{samples.condition}",samples=samples.itertuples()),
    ref_level = config["diffexp"]["ref_level"],
    rmproj_list = config["filtering"]["rmproj"]

  threads: get_deseq2_threads()

  script:
    "../03_Script/deseq2.R"

# ----------------------------------------------
# Change the gene id by the annotation id 
# ----------------------------------------------

# rule annotation_id:
#   input:
#     normalized_counts_file = "05_Output/11_deseq2_init/normalized_counts.tsv",
#     annotationfile = REF + config["microbeannotator"]["annotationfolder"] + "/annotation_results/" + config["microbeannotator"]["protfiles"] + ".annot"


#   output:
#     normalized_counts_annotation_file = report("05_Output/11_deseq2_init/normalized_annotation_counts.tsv", caption="../07_Report/normalized_counts.rst", category="05 Count matrices")

#   shell:
#     """
#     chmod +x 03_Script/gene_id_correspondence.sh
#     03_Script/gene_id_correspondence.sh {input.normalized_counts_file} {input.annotationfile} {output.normalized_counts_annotation_file}
#     """

# ----------------------------------------------
# Statistacal analyses and production of the report with Rmarkdown
# ----------------------------------------------

rule diffexp:
  input:
    rds = "05_Output/11_deseq2_init/all.rds",
    rmd = "03_Script/diffexp.Rmd",
    pca = "03_Script/data_quality.R",
    coldata = config["coldata"],
    annotationfile = REF + config["microbeannotator"]["annotationfolder"] + "/annotation_results/" + config["microbeannotator"]["protfiles"] + ".annot",
    ko_list = config["ko_list"]

  output:
    html_report=report(expand(OUTPUTDIR + config["diffexp"]["html_report"]), caption="../07_Report/diffexp.rst", category="07 Report differential expression"),
    table=report(expand(OUTPUTDIR + "12_differential_expression/{condition.condition}_vs_{ref_level}_all_genes_stats.tsv", condition=condition.itertuples(), ref_level=ref_level), caption="../07_Report/stat.rst", category="06 Table differential expression"),
    sur=report(expand(OUTPUTDIR + "12_differential_expression/{condition.condition}_vs_{ref_level}_signif-up-regulated.txt", condition=condition.itertuples(), ref_level=ref_level), caption="../07_Report/stat.rst", category="06 Table differential expression"),
    sdr=report(expand(OUTPUTDIR + "12_differential_expression/{condition.condition}_vs_{ref_level}_signif-down-regulated.txt", condition=condition.itertuples(), ref_level=ref_level), caption="../07_Report/stat.rst", category="06 Table differential expression"),
    metabolo=report(expand(OUTPUTDIR + "12_differential_expression/{condition.condition}_vs_{ref_level}_metabolic_marker_genes.txt", condition=condition.itertuples(), ref_level=ref_level), caption="../07_Report/stat.rst", category="06 Table differential expression")
  
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
