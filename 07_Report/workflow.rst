This workflow performs metatranscriptome differential expression analysis on single-end RNA-seq data of the diploid strains (NE/CA) and wild type strain (Sc/Sc).
It also does a barcoding analyses to know the taxonomic composition of your samples.

Materials and methods
---------------------

The quality of the raw reads was assessed using FastQC v0.11.9 toolkit (Andrews, 2010). NexteraPE Adapters and low-quality reads were trimmed using Trimmomatic v0.39 (Bolger et al., 2014). SortMeRNA v4.3.4 (Kopylova et al., 2012) is used to filter rRNA fragments from the reads. 

The non rRNA reads were aligned with HiSat2 v2.2.1 (Kim et al., 2015) using default options. 

Gene counts were quantified using feature Counts from the Subread v2.0.1 package (Liao et al., 2019). 

Alignment and gene counts were generated against the reference genome Crocosphaera watsonii WH 8501 (downloaded for JGI Genome Portal). 

The low expressed genes which did not have more than one count per million reads (1CPM) in at least three samples within each dataset were removed from further analysis. 

Gene counts were then normalized and used for differential expression testing using DESeq2 v1.28.0 (Love et al., 2014). 

For the differential RNA-seq analyses, the genes with a Fold change between -2 to 2 and an adjusted p-value (padj) > 0.01 were not considered as significantly expressed. 

Volcanoplot where generated with the software EnhancedVolcano v1.8.0. 

The parameters used for each tool are described on the config.yaml file of the workflow folder.