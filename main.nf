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

  process RUN_QC {
    tag "QC and filtering"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path matrix_dir
    val sample_id

    output:
    path "*_filtered.*", emit: filtered_data
    path "*_pre_filter_violins.png", emit: plots

    script:
    if (params.backend == 'scanpy') {
        """
        run_scanpy_qc.py --matrix_dir ${matrix_dir} \\
                         --min_genes ${params.min_genes} \\
                         --max_mito ${params.max_mito} \\
                         --out_prefix ${sample_id}
        """
    } else if (params.backend == 'seurat') {
        """
        run_seurat_qc.R --matrix_dir ${matrix_dir} \\
                        --min_genes ${params.min_genes} \\
                        --max_mito ${params.max_mito} \\
                        --out_prefix ${sample_id}
        """
    } else {
        error "Unrecognised backend: ${params.backend}. Please choose 'scanpy' or 'seurat'."
    }
}

process RUN_DIMRED {
    tag "PCA and UMAP"
    publishDir "${params.outdir}/dimred", mode: 'copy'

    input:
    path filtered_data
    val sample_id

    output:
    path "*_dimred.*", emit: final_data
    path "*_umap.png", emit: plots

    script:
    if (params.backend == 'scanpy') {
        """
        run_scanpy_dimred.py --input_h5ad ${filtered_data} \\
                             --out_prefix ${sample_id}
        """
    } else if (params.backend == 'seurat') {
        """
        run_seurat_dimred.R --input_rds ${filtered_data} \\
                            --out_prefix ${sample_id}
        """
    } else {
        error "Unrecognised backend."
    }
}
// workflow
workflow {
    matrix_ch = FETCH_DATA(params.input_url)
    RUN_QC(matrix_ch.matrix_dir, 'pbmc3k')
    RUN_DIMRED(qc_ch.filtered_data, 'pbmc3k')
}
