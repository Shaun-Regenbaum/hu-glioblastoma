#!/usr/bin/env python3
"""
Create figures for authentic Neftel analysis with SOX2+/SOX2- splits
Matching their Figure 1B approach with immune cells + malignant/non-malignant separation
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

plt.style.use('default')
sns.set_palette("husl")

def load_authentic_neftel_data():
    """Load authentic Neftel data with SOX2 splits"""
    neftel_dir = Path('/Users/shaunie/Desktop/hu-glioblastoma/results/neftel')
    
    authentic_data = pd.read_csv(neftel_dir / 'Authentic-Neftel-SOX2-Split.csv')
    malignant_summary = pd.read_csv(neftel_dir / 'Malignant-Cells-SOX2-Summary.csv')
    nonmalignant_summary = pd.read_csv(neftel_dir / 'NonMalignant-Cells-SOX2-Summary.csv')
    meta_modules = pd.read_csv(neftel_dir / 'Meta-Modules-SOX2-Split.csv')
    
    return authentic_data, malignant_summary, nonmalignant_summary, meta_modules

def create_authentic_neftel_overview(authentic_data):
    """Create overview figure matching Neftel Figure 1B style"""
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Malignant vs Non-malignant by SOX2 status
    malignant_data = authentic_data.groupby(['dataset', 'malignant', 'sox2_status'])['total_cells'].sum().reset_index()
    
    # Pivot for plotting
    pivot_malignant = malignant_data.pivot_table(
        index=['dataset', 'malignant'], 
        columns='sox2_status', 
        values='total_cells', 
        fill_value=0
    ).reset_index()
    
    # Create grouped bar plot
    datasets = pivot_malignant['dataset'].unique()
    x = np.arange(len(datasets))
    width = 0.15
    
    colors = {'Malignant': ['#FF6B6B', '#FF9999'], 'Non-malignant': ['#4ECDC4', '#7FDDDD']}
    
    for i, dataset in enumerate(datasets):
        dataset_data = pivot_malignant[pivot_malignant['dataset'] == dataset]
        
        malignant_sox2_pos = dataset_data[dataset_data['malignant'] == True]['SOX2+'].iloc[0] if len(dataset_data[dataset_data['malignant'] == True]) > 0 else 0
        malignant_sox2_neg = dataset_data[dataset_data['malignant'] == True]['SOX2-'].iloc[0] if len(dataset_data[dataset_data['malignant'] == True]) > 0 else 0
        nonmal_sox2_pos = dataset_data[dataset_data['malignant'] == False]['SOX2+'].iloc[0] if len(dataset_data[dataset_data['malignant'] == False]) > 0 else 0
        nonmal_sox2_neg = dataset_data[dataset_data['malignant'] == False]['SOX2-'].iloc[0] if len(dataset_data[dataset_data['malignant'] == False]) > 0 else 0
        
        ax1.bar(i - width*1.5, malignant_sox2_pos, width, label='Malignant SOX2+' if i == 0 else "", color=colors['Malignant'][0])
        ax1.bar(i - width*0.5, malignant_sox2_neg, width, label='Malignant SOX2-' if i == 0 else "", color=colors['Malignant'][1])
        ax1.bar(i + width*0.5, nonmal_sox2_pos, width, label='Non-malignant SOX2+' if i == 0 else "", color=colors['Non-malignant'][0])
        ax1.bar(i + width*1.5, nonmal_sox2_neg, width, label='Non-malignant SOX2-' if i == 0 else "", color=colors['Non-malignant'][1])
    
    ax1.set_title('Cell Distribution: Malignant vs Non-malignant by SOX2\\n(Authentic Neftel Approach)', fontweight='bold')
    ax1.set_xlabel('Dataset')
    ax1.set_ylabel('Total Cells')
    ax1.set_xticks(x)
    ax1.set_xticklabels(datasets)
    ax1.legend()
    
    # Plot 2: Cell type composition with SOX2 splits
    cell_type_summary = authentic_data.groupby(['dataset', 'base_cell_type', 'sox2_status'])['total_cells'].sum().reset_index()
    
    # Focus on main cell types for clarity
    main_types = ['Astrocytic', 'Mesenchymal', 'Neural_Progenitor', 'Oligodendrocytic']
    main_data = cell_type_summary[cell_type_summary['base_cell_type'].isin(main_types)]
    
    # Create stacked bar for each dataset
    for i, dataset in enumerate(['GSE131928', 'Organoid']):
        dataset_data = main_data[main_data['dataset'] == dataset]
        pivot_data = dataset_data.pivot(index='base_cell_type', columns='sox2_status', values='total_cells').fillna(0)
        
        if len(pivot_data) > 0:
            bottom_pos = np.zeros(len(pivot_data))
            bottom_neg = np.zeros(len(pivot_data))
            
            # Stack SOX2+ and SOX2- for each cell type
            if 'SOX2+' in pivot_data.columns:
                ax2.bar(np.arange(len(pivot_data)) + i*0.4, pivot_data['SOX2+'], 0.35, 
                       label=f'{dataset} SOX2+' if i == 0 else "", alpha=0.8)
                bottom_pos = pivot_data['SOX2+']
            
            if 'SOX2-' in pivot_data.columns:
                ax2.bar(np.arange(len(pivot_data)) + i*0.4, pivot_data['SOX2-'], 0.35, 
                       bottom=bottom_pos, label=f'{dataset} SOX2-' if i == 0 else "", alpha=0.6)
    
    ax2.set_title('Malignant Cell Types by SOX2 Status\\n(4-State Neftel Analysis)', fontweight='bold')
    ax2.set_xlabel('Malignant Cell Type')
    ax2.set_ylabel('Total Cells')
    if len(main_data) > 0:
        ax2.set_xticks(np.arange(len(pivot_data)) + 0.2)
        ax2.set_xticklabels(pivot_data.index, rotation=45)
    ax2.legend()
    
    # Plot 3: Immune cell analysis
    immune_data = authentic_data[authentic_data['base_cell_type'].isin(['Macrophage', 'T_cell'])]
    if len(immune_data) > 0:
        immune_summary = immune_data.groupby(['dataset', 'base_cell_type', 'sox2_status'])['total_cells'].sum().reset_index()
        
        # Create immune cell plot
        immune_pivot = immune_summary.pivot_table(
            index='base_cell_type', 
            columns=['dataset', 'sox2_status'], 
            values='total_cells', 
            fill_value=0
        )
        
        immune_pivot.plot(kind='bar', ax=ax3, alpha=0.7)
        ax3.set_title('Immune Cell Analysis\\n(Macrophages & T Cells)', fontweight='bold')
        ax3.set_xlabel('Immune Cell Type')
        ax3.set_ylabel('Total Cells')
        ax3.legend(title='Dataset & SOX2', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax3.tick_params(axis='x', rotation=45)
    else:
        ax3.text(0.5, 0.5, 'No immune cells detected\\nin current analysis', 
                ha='center', va='center', transform=ax3.transAxes, fontsize=12)
        ax3.set_title('Immune Cell Analysis', fontweight='bold')
    
    # Plot 4: SOX2 enrichment by cell type and malignancy
    sox2_enrichment = authentic_data.groupby(['base_cell_type', 'malignant', 'sox2_status'])['total_cells'].sum().reset_index()
    sox2_totals = sox2_enrichment.groupby(['base_cell_type', 'malignant'])['total_cells'].sum().reset_index()
    sox2_totals = sox2_totals.rename(columns={'total_cells': 'total_cells_type'})
    
    sox2_enrichment = sox2_enrichment.merge(sox2_totals, on=['base_cell_type', 'malignant'])
    sox2_enrichment['percentage'] = (sox2_enrichment['total_cells'] / sox2_enrichment['total_cells_type']) * 100
    
    # Focus on SOX2+ percentage
    sox2_pos = sox2_enrichment[sox2_enrichment['sox2_status'] == 'SOX2+']
    
    if len(sox2_pos) > 0:
        # Separate malignant and non-malignant
        malignant_sox2 = sox2_pos[sox2_pos['malignant'] == True]
        nonmal_sox2 = sox2_pos[sox2_pos['malignant'] == False]
        
        x = np.arange(len(sox2_pos['base_cell_type'].unique()))
        width = 0.35
        
        cell_types = sox2_pos['base_cell_type'].unique()
        malignant_pcts = [malignant_sox2[malignant_sox2['base_cell_type'] == ct]['percentage'].iloc[0] 
                         if len(malignant_sox2[malignant_sox2['base_cell_type'] == ct]) > 0 else 0 
                         for ct in cell_types]
        nonmal_pcts = [nonmal_sox2[nonmal_sox2['base_cell_type'] == ct]['percentage'].iloc[0] 
                      if len(nonmal_sox2[nonmal_sox2['base_cell_type'] == ct]) > 0 else 0 
                      for ct in cell_types]
        
        ax4.bar(x - width/2, malignant_pcts, width, label='Malignant', alpha=0.8, color='#FF6B6B')
        ax4.bar(x + width/2, nonmal_pcts, width, label='Non-malignant', alpha=0.8, color='#4ECDC4')
        
        ax4.set_title('SOX2+ Enrichment by Cell Type\\n(Malignant vs Non-malignant)', fontweight='bold')
        ax4.set_xlabel('Cell Type')
        ax4.set_ylabel('SOX2+ Percentage (%)')
        ax4.set_xticks(x)
        ax4.set_xticklabels(cell_types, rotation=45)
        ax4.legend()
        
        # Add percentage labels
        for i, (mal_pct, non_pct) in enumerate(zip(malignant_pcts, nonmal_pcts)):
            if mal_pct > 0:
                ax4.text(i - width/2, mal_pct + 1, f'{mal_pct:.1f}%', ha='center', va='bottom', fontsize=9)
            if non_pct > 0:
                ax4.text(i + width/2, non_pct + 1, f'{non_pct:.1f}%', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    return fig

def create_meta_module_sox2_analysis(meta_modules):
    """Create meta-module analysis with SOX2 splits"""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Plot 1: Meta-module distribution by SOX2 status
    # Use the actual meta_modules data directly since it already has the combinations
    
    # Create heatmap data using the meta_module_sox2 column directly
    heatmap_data = meta_modules.pivot_table(
        index='meta_module_sox2', 
        columns='dataset', 
        values='cell_count', 
        aggfunc='sum',
        fill_value=0
    )
    
    # Reorder to group by module type
    module_order = []
    for base_module in ['AC', 'OPC', 'NPC1', 'NPC2', 'MES1', 'MES2']:
        for sox2_status in ['SOX2+', 'SOX2-']:
            module_sox2 = f"{base_module}_{sox2_status}"
            if module_sox2 in heatmap_data.index:
                module_order.append(module_sox2)
    
    if module_order:
        heatmap_data = heatmap_data.reindex(module_order)
        
        sns.heatmap(heatmap_data, annot=True, fmt='d', cmap='viridis', ax=ax1,
                   cbar_kws={'label': 'Cell Count'})
        ax1.set_title('Meta-Module Expression by SOX2 Status\\n(Neftel 6-Module Approach)', fontweight='bold')
        ax1.set_xlabel('Dataset')
        ax1.set_ylabel('Meta-Module + SOX2 Status')
    
    # Plot 2: SOX2 enrichment by meta-module
    if len(meta_modules) > 0:
        # Calculate SOX2+ percentage for each meta-module
        meta_totals = meta_modules.groupby(['dataset', 'meta_module'])['cell_count'].sum().reset_index()
        meta_totals = meta_totals.rename(columns={'cell_count': 'total_count'})
        
        meta_with_totals = meta_modules.merge(meta_totals, on=['dataset', 'meta_module'])
        meta_with_totals['percentage'] = (meta_with_totals['cell_count'] / meta_with_totals['total_count']) * 100
        
        # Focus on SOX2+ percentages
        sox2_pos_meta = meta_with_totals[meta_with_totals['sox2_status'] == 'SOX2+']
        
        if len(sox2_pos_meta) > 0:
            # Aggregate first to avoid duplicates
            sox2_agg = sox2_pos_meta.groupby(['meta_module', 'dataset'])['percentage'].mean().reset_index()
            pivot_sox2 = sox2_agg.pivot(index='meta_module', columns='dataset', values='percentage').fillna(0)
            
            x = np.arange(len(pivot_sox2))
            width = 0.35
            
            bars1 = ax2.bar(x - width/2, pivot_sox2['GSE131928'] if 'GSE131928' in pivot_sox2.columns else [0]*len(pivot_sox2), 
                           width, label='Patient', alpha=0.8)
            bars2 = ax2.bar(x + width/2, pivot_sox2['Organoid'] if 'Organoid' in pivot_sox2.columns else [0]*len(pivot_sox2), 
                           width, label='Organoid', alpha=0.8)
            
            ax2.set_title('SOX2+ Enrichment by Meta-Module\\n(Patient vs Organoid)', fontweight='bold')
            ax2.set_xlabel('Meta-Module')
            ax2.set_ylabel('SOX2+ Percentage (%)')
            ax2.set_xticks(x)
            ax2.set_xticklabels(pivot_sox2.index)
            ax2.legend()
            
            # Add percentage labels
            for bars in [bars1, bars2]:
                for bar in bars:
                    height = bar.get_height()
                    if height > 0:
                        ax2.text(bar.get_x() + bar.get_width()/2., height + 1,
                                f'{height:.1f}%', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    return fig

def create_methodology_summary():
    """Create summary of authentic Neftel methodology"""
    
    fig, ax = plt.subplots(figsize=(14, 10))
    ax.axis('off')
    
    summary_text = """
AUTHENTIC NEFTEL METHODOLOGY IMPLEMENTATION

🧬 BASED ON NEFTEL ET AL. CELL 2019 - FIGURE 1B

📊 CELL TYPE CLASSIFICATION:
✅ CNA+ Malignant Cells (4 States):
   • Astrocytic (AC-like) - GFAP, S100B markers
   • Mesenchymal (MES-like) - VIM, CD44 markers  
   • Neural Progenitor (NPC-like) - SOX4, DCX markers
   • Oligodendrocytic (OPC-like) - OLIG1, MBP markers

🛡️ Non-Malignant Cells:
   • Macrophages (cyan in original Figure 1B)
   • T cells (green in original Figure 1B)
   • Endothelial cells
   • Cycling cells (mixed malignant/non-malignant)

🎯 KEY FEATURES:
• SOX2+/SOX2- STRATIFICATION: Every cell type split by SOX2 expression
• IMMUNE CELLS INCLUDED: Macrophages and T cells as in original
• MALIGNANT vs NON-MALIGNANT: Separate analysis tracks
• META-MODULES: 6 overlapping modules (AC, OPC, NPC1, NPC2, MES1, MES2)
• OVERLAPPING STATES: Cells can belong to multiple categories

📈 SCIENTIFIC RATIONALE:
This implementation captures the biological reality described in Neftel et al.:
1. Glioblastoma cells exist in plastic, overlapping states
2. SOX2 is a critical stem cell marker across all cell types
3. Immune infiltration is present despite CD45- sorting
4. Malignant and non-malignant cells coexist in tumors

🔬 METHODOLOGICAL ADVANTAGES:
• Preserves cellular plasticity information
• Maintains biological complexity
• Enables SOX2-stratified analysis
• Separates malignant from immune/stromal components
• Allows organoid validation of malignant states specifically

📊 DATA STRUCTURE:
• 16 cell types (8 base types × 2 SOX2 states)
• Malignant: 4 base types × 2 SOX2 states = 8 categories
• Non-malignant: 4 base types × 2 SOX2 states = 8 categories
• Meta-modules: 6 modules × 2 SOX2 states = 12 categories

✅ This approach enables authentic comparison with the foundational
   Neftel et al. methodology while adding SOX2 stratification
"""
    
    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, 
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle="round,pad=0.5", facecolor='lightgreen', alpha=0.3))
    
    ax.set_title('Authentic Neftel Methodology: Complete Implementation', 
                fontsize=16, fontweight='bold', pad=20)
    
    return fig

def main():
    """Generate all authentic Neftel figures"""
    
    print("🧬 Creating Authentic Neftel Analysis Figures")
    print("=" * 60)
    
    # Load data
    authentic_data, malignant_summary, nonmalignant_summary, meta_modules = load_authentic_neftel_data()
    print(f"Loaded: {len(authentic_data)} cell type entries")
    
    # Create output directory
    output_dir = Path('/Users/shaunie/Desktop/hu-glioblastoma/results/neftel/figures/authentic')
    output_dir.mkdir(exist_ok=True)
    
    # Generate figures
    print("\n1. Creating authentic Neftel overview...")
    fig1 = create_authentic_neftel_overview(authentic_data)
    fig1.savefig(output_dir / 'authentic_neftel_overview.png', dpi=300, bbox_inches='tight')
    plt.close(fig1)
    
    print("2. Creating meta-module SOX2 analysis...")
    fig2 = create_meta_module_sox2_analysis(meta_modules)
    fig2.savefig(output_dir / 'meta_module_sox2_analysis.png', dpi=300, bbox_inches='tight')
    plt.close(fig2)
    
    print("3. Creating methodology summary...")
    fig3 = create_methodology_summary()
    fig3.savefig(output_dir / 'authentic_neftel_methodology.png', dpi=300, bbox_inches='tight')
    plt.close(fig3)
    
    print(f"\n✅ All figures saved to: {output_dir}")
    print("\nGenerated figures:")
    print("  1. authentic_neftel_overview.png - Complete cell type analysis")
    print("  2. meta_module_sox2_analysis.png - Meta-module with SOX2 splits")
    print("  3. authentic_neftel_methodology.png - Methodology summary")
    
    print("\n" + "=" * 60)
    print("🎯 Authentic Neftel implementation complete with SOX2 stratification!")
    print("=" * 60)

if __name__ == "__main__":
    main()