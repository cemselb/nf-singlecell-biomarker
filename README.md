# 🧬 scRNA-seq Biomarker Discovery Pipeline
[![Pipeline CI](https://github.com/cemselb/nf-singlecell-biomarker/actions/workflows/ci.yml/badge.svg)](https://github.com/cemselb/nf-singlecell-biomarker/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.10.1-brightgreen.svg)](https://www.nextflow.io/)

A containerised, scalable **Nextflow** pipeline for scRNA-seq processing.

## 🧬 Pipeline Architecture

This pipeline is modular and supports multiple backends (**Scanpy** for Python enthusiasts and **Seurat** for R users).

```mermaid
graph LR
    A[Raw 10x Data] --> B(FETCH_DATA)
    B --> C(RUN_QC)
    C --> D(RUN_DIMRED)
    D --> E[Clustered UMAP Output]
```

# Repo structure

```
nf-singlecell-biomarker/
├── main.nf
├── nextflow.config
├── bin/                    
│   └── run_scanpy_qc.py 
├── assets/                 
│   └── multiqc_config.yml ### not in use
├── .github/workflows/      
│   └── ci.yml   
└── README.md 
```
