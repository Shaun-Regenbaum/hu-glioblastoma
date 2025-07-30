import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# Create output directories if they don't exist
os.makedirs('data/sox2_de', exist_ok=True)
os.makedirs('figures/sox2_de', exist_ok=True)

# Read the annotated data
print("Loading data...")
adata = sc.read_h5ad('write/human_combined_4.h5ad')

# Set up visualization parameters
sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=300, facecolor='white', frameon=True, fontsize=12)
sc.settings.figdir = 'figures/sox2_de'

# Create SOX2 enrichment classification
print("Classifying cells by SOX2 expression...")
sox2_expression = adata[:, 'SOX2'].X.toarray().flatten()
sox2_threshold = np.percentile(sox2_expression, 75)
adata.obs['sox2_enriched'] = (sox2_expression > sox2_threshold).astype(str)

# Perform differential expression analysis using t-test in both directions
print("Performing differential expression analysis...")
# SOX2-high vs SOX2-low (for upregulated genes)
sc.tl.rank_genes_groups(
    adata,
    groupby='sox2_enriched',
    groups=['True'],
    reference='False',
    method='t-test',
    n_genes=2000,
    key_added='sox2_de_up'
)

# SOX2-low vs SOX2-high (for downregulated genes)
sc.tl.rank_genes_groups(
    adata,
    groupby='sox2_enriched',
    groups=['False'],
    reference='True',
    method='t-test',
    n_genes=2000,
    key_added='sox2_de_down'
)

# Extract results into a DataFrame with additional statistics
def get_de_results(adata, key='sox2_de_up', group='True'):
    results = pd.DataFrame({
        'names': adata.uns[key]['names'][group],
        'scores': adata.uns[key]['scores'][group],
        'pvals': adata.uns[key]['pvals'][group],
        'pvals_adj': adata.uns[key]['pvals_adj'][group],
        'logfoldchanges': adata.uns[key]['logfoldchanges'][group]
    })
    results.index = results['names']
    
    # Calculate mean expression in each group more efficiently
    genes = results.index
    exp_matrix = adata[:, genes].X.toarray()
    enriched_mask = adata.obs['sox2_enriched'] == 'True'
    
    results['mean_enriched'] = np.mean(exp_matrix[enriched_mask], axis=0)
    results['mean_non_enriched'] = np.mean(exp_matrix[~enriched_mask], axis=0)
    results['log2mean_enriched'] = np.log2(results['mean_enriched'] + 1)
    results['log2mean_non_enriched'] = np.log2(results['mean_non_enriched'] + 1)
    
    return results

print("Processing differential expression results...")
# Get upregulated genes (SOX2-high vs SOX2-low)
de_results_up = get_de_results(adata, key='sox2_de_up', group='True')

# Get downregulated genes (SOX2-low vs SOX2-high)
# Note: we need to flip the sign of logfoldchanges for consistency
de_results_down = get_de_results(adata, key='sox2_de_down', group='False')
de_results_down['logfoldchanges'] = -de_results_down['logfoldchanges']  # Flip the sign

# Combine results
de_results = pd.concat([de_results_up, de_results_down])

# Print distribution statistics before filtering
print("\nDistribution of fold changes before filtering:")
print("Number of genes with negative fold change:", sum(de_results['logfoldchanges'] < 0))
print("Number of genes with positive fold change:", sum(de_results['logfoldchanges'] > 0))
print("\nMin fold change:", de_results['logfoldchanges'].min())
print("Max fold change:", de_results['logfoldchanges'].max())

# Print number of genes that meet each criterion separately
print("\nNumber of genes meeting individual criteria:")
print("Genes with p-adj < 0.05:", sum(de_results['pvals_adj'] < 0.05))
print("Genes with abs(log2FC) > 0.5:", sum(abs(de_results['logfoldchanges']) > 0.5))
print("Genes with log2FC < -0.5:", sum(de_results['logfoldchanges'] < -0.5))
print("Genes with log2FC > 0.5:", sum(de_results['logfoldchanges'] > 0.5))

# Print number of genes meeting both criteria
print("\nNumber of genes meeting both criteria:")
down_regulated = de_results[
    (de_results['pvals_adj'] < 0.05) & 
    (de_results['logfoldchanges'] < -0.5)
]
up_regulated = de_results[
    (de_results['pvals_adj'] < 0.05) & 
    (de_results['logfoldchanges'] > 0.5)
]
print("Significantly down-regulated:", len(down_regulated))
print("Significantly up-regulated:", len(up_regulated))

# Filter for significantly differentially expressed genes with more lenient thresholds
significant_genes = de_results[
    (de_results['pvals_adj'] < 0.05) &  # Changed from 0.01 to 0.05
    (abs(de_results['logfoldchanges']) > 0.5)  # Changed from 1 to 0.5
].sort_values('pvals_adj')

# Save results
print("Saving results...")
significant_genes.to_csv('data/sox2_de/significant_de_genes.csv')

# Create volcano plot
print("Generating volcano plot...")
plt.figure(figsize=(12, 8))

# Plot all genes
plt.scatter(de_results['logfoldchanges'], 
           -np.log10(de_results['pvals']),
           alpha=0.5,
           color='gray',
           s=20,
           label='Not significant')

# Separate upregulated and downregulated significant genes
sig_up = significant_genes[significant_genes['logfoldchanges'] > 0.5]
sig_down = significant_genes[significant_genes['logfoldchanges'] < -0.5]

# Plot significant upregulated genes
plt.scatter(sig_up['logfoldchanges'],
           -np.log10(sig_up['pvals']),
           alpha=0.8,
           color='red',
           s=30,
           label='Upregulated')

# Plot significant downregulated genes
plt.scatter(sig_down['logfoldchanges'],
           -np.log10(sig_down['pvals']),
           alpha=0.8,
           color='blue',
           s=30,
           label='Downregulated')

# Label top genes (both up and downregulated)
top_up = sig_up.head(10)
top_down = sig_down.head(10)
for idx, row in pd.concat([top_up, top_down]).iterrows():
    plt.annotate(idx, 
                (row['logfoldchanges'], -np.log10(row['pvals'])),
                xytext=(5, 5),
                textcoords='offset points',
                fontsize=8)

# Add threshold lines
plt.axhline(y=-np.log10(0.05), color='gray', linestyle='--', alpha=0.5, label='p-value = 0.05')
plt.axvline(x=-0.5, color='gray', linestyle='--', alpha=0.5, label='log2FC = -0.5')
plt.axvline(x=0.5, color='gray', linestyle='--', alpha=0.5, label='log2FC = 0.5')

plt.xlabel('Log2 Fold Change')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot: SOX2-enriched vs non-enriched cells')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig('figures/sox2_de/volcano_plot.png', dpi=300, bbox_inches='tight')
plt.close()

# Create heatmap of top differentially expressed genes
print("Generating heatmap...")
top_genes = significant_genes.head(50).index.tolist()
if len(top_genes) > 0:
    with plt.rc_context({'figure.figsize': (12, 8)}):
        sc.pl.heatmap(adata, 
                     var_names=top_genes, 
                     groupby='sox2_enriched',
                     show_gene_labels=True,
                     dendrogram=True,
                     show=False,
                     save='top_de_genes_heatmap.png')

# Print summary statistics
print("\nAnalysis Summary:")
print(f"Total number of genes tested: {len(de_results)}")
print(f"Number of significantly differentially expressed genes: {len(significant_genes)}")

# Separate upregulated and downregulated genes
upregulated = significant_genes[significant_genes['logfoldchanges'] > 0].sort_values('logfoldchanges', ascending=False)
downregulated = significant_genes[significant_genes['logfoldchanges'] < 0].sort_values('logfoldchanges')

print(f"\nNumber of upregulated genes: {len(upregulated)}")
print(f"Number of downregulated genes: {len(downregulated)}")

print("\nTop 20 upregulated genes in SOX2-enriched cells:")
print(upregulated.head(20)[['logfoldchanges', 'pvals_adj']])
print("\nTop 20 downregulated genes in SOX2-enriched cells:")
print(downregulated.head(20)[['logfoldchanges', 'pvals_adj']])

# Save summary statistics to file
with open('data/sox2_de/analysis_summary.txt', 'w') as f:
    f.write("SOX2 Differential Expression Analysis Summary\n")
    f.write("===========================================\n\n")
    f.write(f"Total number of genes tested: {len(de_results)}\n")
    f.write(f"Number of significantly differentially expressed genes: {len(significant_genes)}\n")
    f.write(f"Number of upregulated genes: {len(upregulated)}\n")
    f.write(f"Number of downregulated genes: {len(downregulated)}\n")
    f.write(f"Significance criteria: adjusted p-value < 0.05 and |log2 fold change| > 0.5\n\n")
    
    f.write("Top 20 upregulated genes in SOX2-enriched cells:\n")
    f.write(upregulated.head(20)[['logfoldchanges', 'pvals_adj']].to_string())
    f.write("\n\nTop 20 downregulated genes in SOX2-enriched cells:\n")
    f.write(downregulated.head(20)[['logfoldchanges', 'pvals_adj']].to_string())

print("\nAnalysis complete! Results saved in data/sox2_de/ and figures/sox2_de/") 