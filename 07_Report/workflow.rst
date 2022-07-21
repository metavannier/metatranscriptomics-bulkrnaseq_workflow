This workflow performs metatranscriptome analysis on paired-end RNA-seq data from cultures of the diazotroph Crocosphaera to three climate change scenarios with & without organic matter additions.

To Download the files of this report : Ctrl + left clic on the blue link.


Materials and methods
---------------------

The main steps in this workflow are :

.. _clean.smk: https://github.com/metavannier/mio_filella_crocosphaera/blob/master/04_Workflow/clean.smk`_
.. _functional_annotation.smk: https://github.com/metavannier/mio_filella_crocosphaera/blob/master/04_Workflow/functional_annotation.smk`_
.. _taxonomic_annotation.smk: https://github.com/metavannier/mio_filella_crocosphaera/blob/master/04_Workflow/taxonomic_annotation.smk`_
.. _count.smk: https://github.com/metavannier/mio_filella_crocosphaera/blob/master/04_Workflow/count.smk`_
.. _diffexp.smk: https://github.com/metavannier/mio_filella_crocosphaera/blob/master/04_Workflow/diffexp.smk`_

- Quality control and cleaning (`clean.smk`_).

The quality of the raw reads was assessed using FastQC v0.11.9 toolkit (Andrews, 2010). NexteraPE Adapters and low-quality reads were trimmed using Trimmomatic v0.39 (Bolger et al., 2014).

SortMeRNA v4.3.4 (Kopylova et al., 2012) was filtering rRNA gene fragments from the metatranscriptomic data using the 16S reference sequences from the SILVA database (Pruesse et al., 2007) (release 132: SSU Ref NR 90). 

- Taxonomic composition analyses to know the taxonomic composition of the differents samples (`taxonomic_annotation.smk`_).

The taxonomic composition metrics for the cultures composition correspond to the abundance profile of Amplicon Sequence Variants (ASVs).

ASVs were produce using DADA2 (Callahan et al., 2016) in QIIME 2 (Bolyen et al., 2019) with the rRNA gene fragments previously detected. A Naive Bayes classifiers trained on the Silva (release 138) 99% OTUs full-length sequences were applied to obtain the pre-trained taxonomy classifiers used for the taxonomic assignation of the features with the classify-sklearn method (Pedregosa et al. 2011). The table was rarefied to 333 sequences and filtered to exclude samples with less than 333 sequences.

- Fonctional annotation of Crocosphaera watsonii WH 8501 genome (`functional_annotation.smk`_).

We used MicrobeAnnotator v2.0.5 (Ruiz-Perez et al.,2021) to update the fonctional annotation of Crocosphaera watsonii WH 8501 genes (downloaded from JGI Genome Portal, GenBank GI #67858163).

Kofam (date: April 27, 2022) (Aramaki T et al. 2020), UniProt’s Swissprot, trEMBL (date:January 2022) (UniProt C. 2019) and NCBI's RefSeq (release 211) (O’Leary NA et al., 2016) databases were used to annotate the proteins.

We used Diamond (Buchfink B et al., 2015) as parameter for the searching step in the microbeannotator_db_builder script.

- Alignment and gene counts (`functional_annotation.smk`_).

The non rRNA reads were aligned with HiSat2 v2.2.1 (Kim et al., 2015) using default options.

Gene counts were quantified using feature Counts from the Subread v2.0.1 package (Liao et al., 2019). 

Alignment and gene counts were generated against the reference genome Crocosphaera watsonii WH 8501. 

The low expressed genes which did not have more than one count per million reads (1CPM) in at least three samples within each dataset were removed from further analysis. 

Gene counts were then normalized and used for differential expression testing using DESeq2 v1.28.0 (Love et al., 2014). 

- Differential RNA-seq analyses (`diffexp.smk`_).

For the differential RNA-seq analyses, the genes with a Fold change between -1.5 to 1.5 and an adjusted p-value (padj) > 0.01 were not considered significantly regulated. 

The statistical representation of data where generated with the functions EnhancedVolcano v1.8.0 and ggplot2 v3.3.2 from the R software v4.0. 

The parameters used for each tool are described on the config.yaml file of the workflow folder.