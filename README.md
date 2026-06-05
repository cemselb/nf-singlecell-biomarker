# 🧬 scRNA-seq Automated Nextflow Processing Pipeline
[![Pipeline CI](https://github.com/cemselb/nf-singlecell-biomarker/actions/workflows/ci.yml/badge.svg)](https://github.com/cemselb/nf-singlecell-biomarker/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.10.1-brightgreen.svg)](https://www.nextflow.io/)

A containerised, scalable **Nextflow** pipeline for scRNA-seq processing, machine learning-based annotation, and biomarker discovery.

## 🧬 Pipeline Architecture

This pipeline is modular and supports multiple backends (**Scanpy** for Python enthusiasts and **Seurat** for R users).

```mermaid
flowchart TD
    %% Styling
    classDef process fill:#1768AC,stroke:#333,stroke-width:2px,color:#fff,rx:5px,ry:5px;
    classDef input fill:#f4f4f4,stroke:#666,stroke-width:1px,stroke-dasharray: 5 5;
    classDef output fill:#44A833,stroke:#333,stroke-width:1px,color:#fff;

    %% Nodes
    Input[("☁️ 10x Genomics<br/>Raw Data URL")]
    FETCH["FETCH_DATA<br/>Download & extract"]
    QC["RUN_QC<br/>Filtering & QC metrics"]
    DIMRED["RUN_DIMRED<br/>PCA, UMAP & neighbours"]
    
    %% Parallel ML and DGE Nodes
    BIOMARKERS["FIND_BIOMARKERS<br/>DGE"]
    ANNOTATE["ANNOTATE_CELLS<br/>ML Annotation (CellTypist / SingleR)"]

    %% Output Nodes
    OutBio[\"📊 Top markers CSV<br/>🔥 Heatmaps"/]
    OutAnno[\"🧬 Annotated h5ad/rds<br/>📈 Cell Type Counts"/]

    %% Connections
    Input --> FETCH
    FETCH -->|Raw Matrices| QC
    QC -->|Filtered Data| DIMRED
    
    %% Forking the workflow
    DIMRED -->|Processed obj| BIOMARKERS
    DIMRED -->|Processed obj| ANNOTATE
    
    BIOMARKERS --> OutBio
    ANNOTATE --> OutAnno

    %% Apply Classes Safely
    class Input input;
    class FETCH,QC,DIMRED,BIOMARKERS,ANNOTATE process;
    class OutBio,OutAnno output;
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

# How to
```markdown
curl -s [https://get.nextflow.io](https://get.nextflow.io) | bash
```

To run with **scanpy**
```markdown
nextflow run main.nf -profile conda --backend scanpy
```

To run with **Seurat**
```markdown
nextflow run main.nf -profile conda --backend seurat
```
