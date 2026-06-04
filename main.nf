#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// log pipeline info upon start
log.info """\
    =============================================
    S I N G L E - C E L L   P I P E L I N E
    =============================================
    Dataset      : ${params.input_url}
    Output Dir   : ${params.outdir}
    Min Genes    : ${params.min_genes}
    Max Mito %   : ${params.max_mito}
    =============================================
    """

// processes
process FETCH_DATA {
    tag "Downloading 10x data"
    publishDir "${params.outdir}/raw_data", mode: 'copy'

    input:
    val url

    output:
    path "filtered_gene_bc_matrices", emit: matrix_dir

    script:
    """
    wget -qO pbmc3k.tar.gz ${url}
    tar -xzf pbmc3k.tar.gz
    """
}

process RUN_SCANPY_QC {
    tag "QC and filtering"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path matrix_dir

    output:
    path "filtered_adata.h5ad", emit: h5ad
    path "qc_violins.png", emit: plots

    script:
    // run_scanpy_qc.py --input ${matrix_dir} --min_genes ${params.min_genes}
    """
    python -c "
import scanpy as sc
adata = sc.read_10x_mtx('${matrix_dir}/hg19/', var_names='gene_symbols')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, show=False)
import matplotlib.pyplot as plt
plt.savefig('qc_violins.png')
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < ${params.max_mito}, :]
adata.write('filtered_adata.h5ad')
    "
    """
}

// workflow
workflow {
    // fetch data
    matrix_ch = FETCH_DATA(params.input_url)
    
    // pass fetched data to QC process
    RUN_SCANPY_QC(matrix_ch.matrix_dir)
}
