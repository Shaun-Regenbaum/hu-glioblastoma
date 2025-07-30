# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a single-cell RNA sequencing analysis project studying human glioblastoma (brain tumor) differentiation patterns. The research focuses on identifying and characterizing different tumor cell states and their transitions.

## Tech Stack and Dependencies

- **Python** with Jupyter notebooks for interactive analysis
- **Key packages**: scanpy, numpy, pandas, seaborn, matplotlib, scikit-learn, leidenalg, louvain, umap-learn, celltypist
- **Package manager**: Uses `uv` for fast Python package installation

## Common Commands

```bash
# Install dependencies (if requirements.txt exists)
uv pip install -r requirements.txt

# Start Jupyter notebook
jupyter notebook

# Start Jupyter lab
jupyter lab
```

## Code Architecture

### Data Organization
- `data/raw/`: Original 10x Genomics h5 files (samples 1914, 1919)
- `data/processed/`: Analysis outputs and intermediate files
- `data/external/`: Reference datasets (SRX samples, Well plates)
- `cache/`: Cached computation results
- `figures/`: Generated analysis plots
- `Yuval/`: Project-specific data and analyses

### Analysis Workflow
1. **Data Loading**: Load 10x h5 files using scanpy
2. **Quality Control**: Filter cells based on mitochondrial content, gene counts
3. **Normalization**: Log-normalization and scaling
4. **Dimensionality Reduction**: PCA → UMAP/tSNE
5. **Clustering**: Leiden/Louvain algorithms
6. **Cell Type Annotation**: Based on marker gene expression

### Key Cell States Being Studied
- **Mesenchymal (MES)**: VIM, CD44, CHI3L1, ANXA1, ANXA2
- **Astrocyte-like (AC)**: GFAP, S100B, CST3, HOPX
- **Oligodendrocyte-like (OPC)**: OLIG1, MBP, PLP1, ALCAM
- **Neural Progenitor (NPC)**: SOX4, DCX, CD24, DLL3

### Important Notebooks
- `finalv2.ipynb`, `finalv3.ipynb`: Main analysis pipelines
- `integration.ipynb`, `integration2.ipynb`: Data integration workflows
- `ivygap.ipynb`, `ivygap_integration.ipynb`: IvyGAP dataset integration
- `heatmap.ipynb`: Expression heatmap generation

## Analysis Considerations

When working with this codebase:
- Cell type classification focuses on four main trajectories: NPCs, OPCs, Mesenchymal, and Astrocyte
- Hybrid states exist: AC/MES, NPC/OPC, AC/OPC
- No immune cells expected (tumors grown outside immune system)
- Use established marker genes from README.md for cell type validation
- Cache intermediate results to speed up re-analysis