# ----------------------------------------------
# Qiime2 import multiplexed paired-end FASTQ
# ----------------------------------------------

rule qiime_import:
  input:
    rowdatapath = ROOTDIR + "/manifest.tsv"

  output:
    outputimport = OUTPUTDIR + "05_qiime_import/" + PROJ + "-demux-paired-end.qza"

  conda: 
    ROOTDIR + "/02_Container/qiime2.yaml"

  shell:
    """
    qiime tools import \
      --type {config[qiime][type]} \
      --input-path {input.rowdatapath} \
      --input-format {config[qiime][format]} \
      --output-path {output.outputimport}
    qiime tools validate {output.outputimport}
    """

rule get_stats:
  input:
    q2_primerRM = OUTPUTDIR + "05_qiime_import/" + PROJ + "-demux-paired-end.qza",
  output:
    primer = report(OUTPUTDIR + "05_qiime_import/" + PROJ + "-PE-demux-noprimer.qzv", caption = ROOTDIR + "/07_Report/rawsum.rst", category="02 reads report"),
  conda: 
    ROOTDIR + "/02_Container/qiime2.yaml"
  shell: 
    """
    qiime demux summarize --i-data {input.q2_primerRM} --o-visualization {output.primer}
    """

#########################
# Denoising with dada2  #
#########################

rule dada2:
    input:
        q2_primerRM = OUTPUTDIR + "05_qiime_import/" + PROJ + "-demux-paired-end.qza"
    output:
        table = OUTPUTDIR + "06_dada2/" + PROJ + "-table-dada2.qza",
        rep = OUTPUTDIR + "06_dada2/" + PROJ + "-rep-seqs-dada2.qza",
        stats = OUTPUTDIR + "06_dada2/" + PROJ + "-dada2-stats.qza"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime dada2 denoise-paired \
        --i-demultiplexed-seqs {input.q2_primerRM} \
        --p-trunc-q {config[dada2][truncation_err]} \
        --p-trunc-len-f {config[dada2][truncation_len-f]} \
        --p-trunc-len-r {config[dada2][truncation_len-r]} \
        --p-trim-left-f {config[dada2][trim-left-f]} \
        --p-trim-left-r {config[dada2][trim-left-r]} \
        --p-max-ee-f {config[dada2][quality_err-f]} \
        --p-max-ee-r {config[dada2][quality_err-r]} \
        --p-n-reads-learn {config[dada2][training]} \
        --p-n-threads {config[dada2][threads]} \
        --p-chimera-method {config[dada2][chimera]} \
        --o-table {output.table} \
        --o-representative-sequences {output.rep} \
        --o-denoising-stats {output.stats}
        """

rule filterseq:
    input:
        rep = OUTPUTDIR + "06_dada2/" + PROJ + "-rep-seqs-dada2.qza",
    output:
        filterrep = OUTPUTDIR + "06_dada2/" + PROJ + "-rep-filtered-seqs-dada2.qza",
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime feature-table filter-seqs \
        --i-data {input.rep} \
        --m-metadata-file {input.rep} \
        --p-where 'length(sequence) >= {config[dada2][seqlenth]}' \
        --o-filtered-data {output.filterrep} 
        """

rule filterfeature:
    input:
        table = OUTPUTDIR + "06_dada2/" + PROJ + "-table-dada2.qza",
    output:
        filtertable = OUTPUTDIR + "06_dada2/" + PROJ + "-table-filtered-dada2.qza",
    params:
        filterrep = OUTPUTDIR + "06_dada2/" + PROJ + "-rep-filtered-seqs-dada2.qza",
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime feature-table filter-features \
        --i-table {input.table} \
        --m-metadata-file {params.filterrep} \
        --o-filtered-table {output.filtertable}        
        """

# ----------------------------------------------
# Assigning taxonomy
# ----------------------------------------------

rule assign_tax:
    input:
        q2_repseq = OUTPUTDIR + "06_dada2/" + PROJ + "-rep-filtered-seqs-dada2.qza",
        db_classified = REF + DB_classifier
    output:
        sklearn = OUTPUTDIR + "07_taxonomy/" + PROJ + "-tax_sklearn.qza"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"  
    shell:
        """
        qiime feature-classifier classify-sklearn \
        --i-classifier {input.db_classified} \
        --i-reads {input.q2_repseq} \
        --o-classification {output.sklearn} \
        --p-confidence {config[taxonomy][confidence]} \
        --p-n-jobs {config[taxonomy][njobs]}
        """

rule taxa_filter_table:
    input:
        filtertable = OUTPUTDIR + "06_dada2/" + PROJ + "-table-filtered-dada2.qza",
        sklearn = OUTPUTDIR + "07_taxonomy/" + PROJ + "-tax_sklearn.qza"
    output:
        taxafiltertable = OUTPUTDIR + "07_taxonomy/" + PROJ + "-taxa-table-filtered-dada2.qza"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell: 
        """
        qiime taxa filter-table \
        --i-table {input.filtertable} \
        --i-taxonomy {input.sklearn} \
        --p-exclude {config[taxonomy][filter]} \
        --o-filtered-table {output.taxafiltertable}
        """

rule filter_seq:
    input:
        q2_repseq = OUTPUTDIR + "06_dada2/" + PROJ + "-rep-filtered-seqs-dada2.qza",
        taxafiltertable = OUTPUTDIR + "07_taxonomy/" + PROJ + "-taxa-table-filtered-dada2.qza"
    output:
        q2_repseq_filtered = OUTPUTDIR + "07_taxonomy/" + PROJ + "-rep-filtered-seqs-taxa-dada2.qza"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell: 
        """
        qiime feature-table filter-seqs \
        --i-data {input.q2_repseq} \
        --i-table {input.taxafiltertable} \
        --o-filtered-data {output.q2_repseq_filtered}
        """

rule gen_tax:
    input:
        sklearn = OUTPUTDIR + "07_taxonomy/" + PROJ + "-tax_sklearn.qza"
    output:
        table_tax = OUTPUTDIR + "07_taxonomy/taxonomy.tsv",
        table_tax_filtered = report(OUTPUTDIR + "07_taxonomy/taxonomy_filtered.tsv", caption = ROOTDIR + "/07_Report/tax.rst", category="04 taxonomy"),
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    params:
        directory(OUTPUTDIR + "07_taxonomy/")
    shell:
        """
        qiime tools export --input-path {input.sklearn} --output-path {params}
        filter=({config[taxonomy][filter]})
        words="${{filter//,/|}}"
        grep -Ev $words {output.table_tax} > {output.table_tax_filtered}
        """

rule stats:
    input:
        filterrep = OUTPUTDIR + "07_taxonomy/" + PROJ + "-rep-filtered-seqs-taxa-dada2.qza",
        stats = OUTPUTDIR + "06_dada2/" + PROJ + "-dada2-stats.qza",
        taxafiltertable = OUTPUTDIR + "07_taxonomy/" + PROJ + "-taxa-table-filtered-dada2.qza",
    output:
        rep_viz = report(OUTPUTDIR + "07_taxonomy/" + PROJ + "-rep-filtered-seqs-taxa-dada2.qzv", caption = ROOTDIR + "07_Report/dada2seq.rst", category="03 dada2"),
        stats_viz = report(OUTPUTDIR + "06_dada2/" + PROJ + "-dada2-stats.qzv", caption = ROOTDIR + "07_Report/dada2summary.rst", category="03 dada2"),
        featurestat = report(OUTPUTDIR + "07_taxonomy/" + PROJ + "-taxa-table-filtered-dada2.qzv", caption = ROOTDIR + "07_Report/dada2summary.rst", category="04 taxonomy"),
    params:
        metadata = ROOTDIR + "/sample-metadata.tsv"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime metadata tabulate --m-input-file {input.stats} --o-visualization {output.stats_viz}
        qiime feature-table tabulate-seqs --i-data {input.filterrep} --o-visualization {output.rep_viz}
        qiime feature-table summarize --i-table {input.taxafiltertable} --m-sample-metadata-file {params.metadata} --o-visualization {output.featurestat} 
        """

rule filter_rarefaction:
    input:
        taxafiltertable = OUTPUTDIR + "07_taxonomy/" + PROJ + "-taxa-table-filtered-dada2.qza"
    output:
        rarefactionfiltertable = OUTPUTDIR + "07_taxonomy/" + PROJ + "-rarefaction-table-filtered-dada2.qza"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime feature-table filter-samples \
        --i-table {input.taxafiltertable} \
        --p-min-frequency {config[diversity][samplingdepth]} \
        --o-filtered-table {output.rarefactionfiltertable}
        """

rule relative_frequency:
    input:
        rarefactionfiltertable = OUTPUTDIR + "07_taxonomy/" + PROJ + "-rarefaction-table-filtered-dada2.qza"
    output:
        relativefreqtable = OUTPUTDIR + "07_taxonomy/" + PROJ + "-relative-frequency-dada2.qza"
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        qiime feature-table relative-frequency \
        --i-table {input.rarefactionfiltertable} \
        --o-relative-frequency-table {output.relativefreqtable}
        """

rule gen_table:
    input:
        rarefactionfiltertable = OUTPUTDIR + "07_taxonomy/" + PROJ + "-rarefaction-table-filtered-dada2.qza"
    output:
        table_biom = OUTPUTDIR + "07_taxonomy/feature-table.biom",
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    params:
        directory(OUTPUTDIR + "07_taxonomy/")
    shell:
        """
        qiime tools export --input-path {input.rarefactionfiltertable} --output-path {params}
        """

rule convert:
    input:
        table_biom = OUTPUTDIR + "07_taxonomy/feature-table.biom",
        taxo_table = OUTPUTDIR + "07_taxonomy/taxonomy.tsv",
    output:
        taxo_table_biom = OUTPUTDIR + "07_taxonomy/" + PROJ + "-asv-table-with-taxonomy.biom",
        taxo_table_tsv = report(OUTPUTDIR + "07_taxonomy/" + PROJ + "-asv-table-with-taxonomy.tsv", caption = ROOTDIR + "/07_Report/asv_table.rst", category="03 dada2"),
    params:
        directory(OUTPUTDIR + "07_taxonomy/")
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    shell:
        """
        taxo_table={input.taxo_table}
        var="#OTUID\ttaxonomy\tconfidence"
        sed -i "1s/.*/$var/" ${{taxo_table}}
        biom add-metadata -i {input.table_biom} -o {output.taxo_table_biom} --observation-metadata-fp {input.taxo_table} --sc-separated taxonomy
        biom convert -i {output.taxo_table_biom} -o {output.taxo_table_tsv} --to-tsv --header-key taxonomy
        """

rule taxa_barplot:
    input:
        rarefactionfiltertable = OUTPUTDIR + "07_taxonomy/" + PROJ + "-rarefaction-table-filtered-dada2.qza",
        sklearn = OUTPUTDIR + "07_taxonomy/" + PROJ + "-tax_sklearn.qza",
    output:
        taxabarplots = report(OUTPUTDIR + "07_taxonomy/" + PROJ + "-taxa-bar-plots.qzv", caption = ROOTDIR + "/07_Report/taxbarplot.rst", category="04 taxonomy")
    conda:
        ROOTDIR + "/02_Container/qiime2.yaml"
    params:
        metadata = ROOTDIR + "/sample-metadata.tsv"
    shell:
        """
        qiime taxa barplot \
        --i-table {input.rarefactionfiltertable} \
        --i-taxonomy {input.sklearn} \
        --m-metadata-file {params.metadata} \
        --o-visualization {output.taxabarplots}
        """