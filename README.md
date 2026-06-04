# 🧬 scRNA-seq Biomarker Discovery Pipeline

![Nextflow](https://img.shields.io/badge/Nextflow-1DB182?style=for-the-badge&logo=nextflow&logoColor=white)
![Docker](https://img.shields.io/badge/Docker-2496ED?style=for-the-badge&logo=docker&logoColor=white)
![Python](https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white)

A containerised, scalable Nextflow pipeline for scRNA-seq processing.

## 🏗️ Architecture
1. **Data:** Automatically fetches public 10x Genomics `.mtx` matrices (defaults to 3k PBMCs).
2. **QC:** Filters ambient RNA, dead cells (mitochondrial thresholds), and doublets using `Scanpy`.
3. **Dimensionality reduction:** Performs highly variable gene (HVG) selection, PCA, and UMAP projection.

# Repo structure

```
nf-singlecell-biomarker/
├── main.nf
├── nextflow.config
├── modules/                
│   ├── fetch_data.nf 
│   ├── qc_filter.nf
│   └── dim_reduction.nf
├── bin/                    
│   └── run_scanpy_qc.py 
├── assets/                 
│   └── multiqc_config.yml 
├── .github/workflows/      
│   └── ci.yml   
└── README.md 
```


## How to
```bash
# Clone the repository
git clone [https://github.com/cemselb/nf-singlecell-biomarker.git](https://github.com/cemselb/nf-singlecell-biomarker.git)
cd nf-singlecell-biomarker

# Run the pipeline with the Docker profile
nextflow run main.nf -profile docker
