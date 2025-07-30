    import scanpy as sc
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Read the annotated data
    adata = sc.read_h5ad('write/human_combined_4.h5ad')

    # Annotating types
    adata.obs.loc[adata.obs.leiden == '0', 'cell_type_f'] = 'Mesenchymal'
    adata.obs.loc[adata.obs.leiden == '7', 'cell_type_f'] = 'Mesenchymal'
    adata.obs.loc[adata.obs.leiden == '9', 'cell_type_f'] = 'Mesenchymal'
    adata.obs.loc[adata.obs.leiden == '10', 'cell_type_f'] = 'Mesenchymal'
    adata.obs.loc[adata.obs.leiden == '6', 'cell_type_f'] = 'Mesenchymal'



    adata.obs.loc[adata.obs.leiden == '8', 'cell_type_f'] = 'Astrocyte'
    adata.obs.loc[adata.obs.leiden == '4', 'cell_type_f'] = 'Astrocyte'

    adata.obs.loc[adata.obs.leiden == '3', 'cell_type_f'] = 'Oligodendrocyte and Nueral Progenitor Cells'
    adata.obs.loc[adata.obs.leiden == '5', 'cell_type_f'] = 'Oligodendrocyte and Nueral Progenitor Cells'
    adata.obs.loc[adata.obs.leiden == '2', 'cell_type_f'] = 'Oligodendrocyte and Nueral Progenitor Cells'

    adata.obs.loc[adata.obs.leiden == '3', 'cell_type_f'] = 'Cycling Cells'

    adata.obs.loc[adata.obs.leiden == '12', 'cell_type_f'] = 'Endothelial'

    # remove all cells in leiden group 3 from adata
    # adata = adata[adata.obs.leiden != '3']

    # Create binary classification for OPC vs non-OPC cells

    sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=True, fontsize=20, figsize=(12, 12))
    palette = sns.color_palette("Paired")
    # convert to hex
    palette = palette.as_hex()
    # change the order of colors
    palette = [palette[0], palette[8], palette[2], palette[6], palette[4]]
    # make the last a light gray
    palette.append('#D3D3D3')

    sc.pl.umap(adata, color="leiden", title="Annotated Cell Types", size=120, alpha=0.6, palette=palette, frameon=False, colorbar_loc=None)
    sc.pl.umap(adata, color="cell_type_f", title="Annotated Cell Types", size=120, alpha=0.6, palette=palette, frameon=False, colorbar_loc=None)

    adata.obs['is_opc'] = (adata.obs['leiden'] == '2').astype(str)

    # Perform differential expression analysis
    de_results = sc.tl.rank_genes_groups(
        adata,
        groupby='cell_type_f',
        groups=['Oligodendrocyte and Nueral Progenitor Cells'],  # Compare OPC cells against all others
        reference='Astrocyte',  # Non-OPC cells are the reference
        method='wilcoxon',
        key_added='opc_vs_rest'
    )

    # Extract results into a DataFrame
    def get_de_results(adata, key='opc_vs_rest', group='Oligodendrocyte and Nueral Progenitor Cells'):
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

    significant_genes
    markers = significant_genes.index
    # limit to 10
    markers = markers[:20]
    sc.pl.matrixplot(adata, markers, groupby='cell_type_f', dendrogram=False)
    # Create a copy of adata with only significant genes
    adata_subset = adata[:, significant_genes.index]

    # Save results to CSV
    significant_genes.to_csv('cycling_de_genes.csv')



    # Print summary statistics
    print(f"Total number of significant DE genes: {len(significant_genes)}")
    print("\nTop 20 differentially expressed genes in cycling cells:")
    print(significant_genes.head(20)[['logfoldchanges', 'pvals']])

    # Print cell type distribution
    print("\nDistribution of cells across cell types:")
    print(adata.obs['cell_type_f'].value_counts())

    sc.pl.umap(adata, color='PPP1R9A', size=120, alpha=0.8, cmap='RdBu_r' ,frameon=False, legend_loc=None, colorbar_loc=None)
