import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Create output directories if they don't exist
os.makedirs('figures/sox2_enrichment', exist_ok=True)
os.makedirs('data/sox2_enrichment', exist_ok=True)

# Read the annotated data
adata = sc.read_h5ad('write/human_combined_4.h5ad')

# Set up visualization parameters
sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=300, facecolor='white', frameon=True, fontsize=12)

# Set the scanpy figure saving directory
sc.settings.figdir = 'figures/sox2_enrichment'

# Create a binary classification for SOX2-enriched vs non-enriched cells
# We'll define SOX2-enriched cells as those with expression above the 75th percentile
sox2_expression = adata[:, 'SOX2'].X.toarray().flatten()
sox2_threshold = np.percentile(sox2_expression, 75)
adata.obs['sox2_enriched'] = (sox2_expression > sox2_threshold).astype(str)

# Visualize SOX2 expression on UMAP
with plt.rc_context({'figure.figsize': (10, 10)}):
    sc.pl.umap(adata, 
               color='SOX2', 
               size=120, 
               alpha=0.8, 
               cmap='RdBu_r', 
               title='SOX2 Expression Level',
               frameon=False,
               show=False,
               save='sox2_expression_umap.png')

# Visualize SOX2 enrichment status on UMAP
with plt.rc_context({'figure.figsize': (10, 10)}):
    sc.pl.umap(adata, 
               color='sox2_enriched', 
               size=120, 
               alpha=0.8, 
               title='SOX2 Enrichment Status',
               frameon=False,
               show=False,
               save='sox2_enrichment_status_umap.png')

# Perform differential expression analysis
de_results = sc.tl.rank_genes_groups(
    adata,
    groupby='sox2_enriched',
    groups=['True'],  # Compare SOX2-enriched cells against all others
    reference='False',  # Non-enriched cells are the reference
    method='wilcoxon',
    key_added='sox2_enriched_vs_rest'
)

# Extract results into a DataFrame
def get_de_results(adata, key='sox2_enriched_vs_rest', group='True'):
    results = pd.DataFrame({
        'names': adata.uns[key]['names'][group],
        'scores': adata.uns[key]['scores'][group],
        'pvals': adata.uns[key]['pvals'][group],
        'pvals_adj': adata.uns[key]['pvals_adj'][group],
        'logfoldchanges': adata.uns[key]['logfoldchanges'][group]
    })
    results.index = results['names']
    return results

de_df = get_de_results(adata)

# Filter for significantly differentially expressed genes
significant_genes = de_df[
    (de_df['pvals_adj'] < 0.01) & 
    (abs(de_df['logfoldchanges']) > 1)
].sort_values('pvals_adj')

# Save results to CSV
significant_genes.to_csv('data/sox2_enrichment/sox2_enriched_de_genes.csv')

# Create expression matrix for significant genes using leiden clusters
cell_types = adata.obs['leiden'].unique()
expression_matrix = pd.DataFrame(index=significant_genes.index)

for cell_type in cell_types:
    mask = adata.obs['leiden'] == cell_type
    cell_type_mean = np.array(adata[:, significant_genes.index].X[mask.values].mean(axis=0)).flatten()
    expression_matrix[cell_type] = cell_type_mean

# Save the expression matrix to CSV
expression_matrix.to_csv('data/sox2_enrichment/sox2_enriched_genes_expression_matrix.csv')

# Create violin plots for top genes
if len(significant_genes) > 0:
    top_genes = significant_genes.index[:min(10, len(significant_genes))]
    with plt.rc_context({'figure.figsize': (10, 6)}):
        sc.pl.violin(adata, 
                     top_genes, 
                     groupby='sox2_enriched', 
                     rotation=45,
                     show=False,
                     save='sox2_enriched_top_genes_violin.png')

    # Create matrix plot for top genes using leiden clusters
    with plt.rc_context({'figure.figsize': (12, 6)}):
        sc.pl.matrixplot(adata, 
                         top_genes, 
                         groupby='leiden', 
                         dendrogram=False,
                         show=False,
                         save='sox2_enriched_markers_matrix.png')

# Print summary statistics
print(f"Total number of significant DE genes: {len(significant_genes)}")
print("\nTop 20 differentially expressed genes in SOX2-enriched cells:")
print(significant_genes.head(20)[['logfoldchanges', 'pvals']])

# Print cell type distribution
print("\nDistribution of cells across SOX2 enrichment status:")
print(adata.obs['sox2_enriched'].value_counts())
print("\nPercentage of cells in each category:")
print(adata.obs['sox2_enriched'].value_counts(normalize=True) * 100)

# Distribution of SOX2-enriched cells across leiden clusters
print("\nDistribution of SOX2-enriched cells across clusters:")
cluster_enrichment = pd.crosstab(adata.obs['leiden'], adata.obs['sox2_enriched'], normalize='index') * 100
print(cluster_enrichment)

# Save cluster enrichment to CSV
cluster_enrichment.to_csv('data/sox2_enrichment/cluster_enrichment_stats.csv')

# Check for available stem cell markers
potential_markers = ['SOX2', 'NANOG', 'MYC', 'KLF4']
available_markers = [gene for gene in potential_markers if gene in adata.var_names]

if available_markers:
    with plt.rc_context({'figure.figsize': (10, 6)}):
        sc.pl.violin(adata, 
                     available_markers, 
                     groupby='sox2_enriched', 
                     rotation=45,
                     show=False,
                     save='available_stem_cell_markers.png') 