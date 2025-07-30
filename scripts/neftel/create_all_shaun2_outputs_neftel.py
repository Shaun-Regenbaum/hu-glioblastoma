#!/usr/bin/env python3
"""
Create ALL Shaun2-style outputs using authentic Neftel methodology
- FourSample-FourPop CSVs
- CellCounts, Percentages, Methods
- OrganoidVsPatient comparisons
- All figures matching Shaun2 exactly
- NO immune cells in organoids (they're grown in vitro)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set style to match Shaun2
plt.style.use('default')
sns.set_palette("husl")

def load_and_prepare_neftel_data():
    """Load Neftel data and prepare for Shaun2-style outputs"""
    
    # Load the authentic Neftel data
    neftel_data = pd.read_csv('/Users/shaunie/Desktop/hu-glioblastoma/results/neftel/Authentic-Neftel-SOX2-Split.csv')
    
    # CRITICAL: Remove immune cells from organoids (they don't exist in vitro)
    print("Removing immune cells from organoid samples...")
    neftel_data = neftel_data[~(
        (neftel_data['dataset'] == 'Organoid') & 
        (neftel_data['base_cell_type'].isin(['Macrophage', 'T_cell']))
    )]
    
    # For Shaun2 compatibility, we'll focus on the 4 main malignant states + cycling + endothelial
    # Merge SOX2+ and SOX2- back together for main analysis (but keep split data for SOX2 analysis)
    merged_data = neftel_data.groupby(['dataset', 'sample', 'base_cell_type', 'malignant']).agg({
        'total_cells': 'sum',
        'SOX2_expressing': 'sum',
        'percentage': 'sum'  # This will give us overlapping percentages
    }).reset_index()
    
    # Calculate SOX2 percentage
    merged_data['SOX2_percentage'] = (merged_data['SOX2_expressing'] / merged_data['total_cells']) * 100
    
    return neftel_data, merged_data

def create_four_sample_four_pop_csvs(merged_data):
    """Create FourSample-FourPop CSVs exactly like Shaun2"""
    
    print("\nCreating FourSample-FourPop CSVs...")
    
    # Focus on the 4 main malignant populations
    four_pops = ['Astrocytic', 'Mesenchymal', 'Neural_Progenitor', 'Oligodendrocytic']
    four_pop_data = merged_data[merged_data['base_cell_type'].isin(four_pops)]
    
    # Get 4 representative samples from each dataset
    patient_samples = ['BT749', 'BT771', 'BT830', 'BT1160']
    organoid_samples = ['Organoid_D77', 'Organoid_D88-1', 'Organoid_D88-3', 'Organoid_D133-1']
    
    # Filter for these samples
    four_sample_data = four_pop_data[
        four_pop_data['sample'].isin(patient_samples + organoid_samples)
    ]
    
    # Create the main FourSample-FourPop CSV
    output_path = Path('/Users/shaunie/Desktop/hu-glioblastoma/results/neftel')
    four_sample_data.to_csv(output_path / 'FourSample-FourPop-Neftel.csv', index=False)
    
    # Create summary by cell type
    summary = four_sample_data.groupby(['dataset', 'base_cell_type']).agg({
        'total_cells': 'sum',
        'SOX2_expressing': 'sum'
    }).reset_index()
    summary['SOX2_percentage'] = (summary['SOX2_expressing'] / summary['total_cells']) * 100
    
    summary.to_csv(output_path / 'FourPop-Summary-Neftel.csv', index=False)
    
    return four_sample_data

def create_cell_counts_percentages_methods(merged_data):
    """Create CellCounts, Percentages, and Methods CSVs"""
    
    print("Creating CellCounts, Percentages, and Methods CSVs...")
    
    output_path = Path('/Users/shaunie/Desktop/hu-glioblastoma/results/neftel')
    
    # 1. Cell Counts CSV
    cell_counts = merged_data.pivot_table(
        index='sample',
        columns='base_cell_type',
        values='total_cells',
        fill_value=0
    )
    cell_counts.to_csv(output_path / 'CellCounts-Neftel.csv')
    
    # 2. Percentages CSV (note: may exceed 100% due to overlapping states)
    percentages = merged_data.pivot_table(
        index='sample',
        columns='base_cell_type',
        values='percentage',
        fill_value=0
    )
    percentages.to_csv(output_path / 'Percentages-Neftel.csv')
    
    # 3. Methods CSV
    methods_data = {
        'Method': ['Neftel et al. (2019)'],
        'Description': ['Authentic implementation of Neftel methodology with 6 meta-modules'],
        'Cell_Types': ['AC-like, MES-like (MES1/MES2), NPC-like (NPC1/NPC2), OPC-like'],
        'Overlapping': ['Yes - cells can express multiple programs'],
        'Immune_Cells': ['Macrophages and T cells in patients only'],
        'SOX2_Analysis': ['All cell types split into SOX2+ and SOX2-'],
        'Reference': ['Neftel et al., Cell 2019, Figure 1B approach']
    }
    methods_df = pd.DataFrame(methods_data)
    methods_df.to_csv(output_path / 'Methods-Neftel.csv', index=False)
    
    return cell_counts, percentages

def create_organoid_vs_patient_comparison(merged_data, neftel_data):
    """Create OrganoidVsPatient comparison focusing on 4 main populations"""
    
    print("Creating OrganoidVsPatient FourPop comparison...")
    
    # Focus on 4 main malignant populations
    four_pops = ['Astrocytic', 'Mesenchymal', 'Neural_Progenitor', 'Oligodendrocytic']
    four_pop_data = merged_data[merged_data['base_cell_type'].isin(four_pops)]
    
    # Calculate totals by dataset
    comparison = four_pop_data.groupby(['dataset', 'base_cell_type']).agg({
        'total_cells': 'sum',
        'SOX2_expressing': 'sum'
    }).reset_index()
    
    # Calculate dataset totals
    dataset_totals = comparison.groupby('dataset')['total_cells'].sum()
    comparison['percentage_of_dataset'] = comparison.apply(
        lambda x: (x['total_cells'] / dataset_totals[x['dataset']]) * 100, axis=1
    )
    comparison['SOX2_percentage'] = (comparison['SOX2_expressing'] / comparison['total_cells']) * 100
    
    # Save comparison
    output_path = Path('/Users/shaunie/Desktop/hu-glioblastoma/results/neftel')
    comparison.to_csv(output_path / 'OrganoidVsPatient-FourPop-Neftel.csv', index=False)
    
    # Create SOX2 split comparison
    sox2_comparison = neftel_data[neftel_data['base_cell_type'].isin(four_pops)].groupby(
        ['dataset', 'base_cell_type', 'sox2_status']
    )['total_cells'].sum().reset_index()
    
    sox2_comparison.to_csv(output_path / 'OrganoidVsPatient-SOX2Split-Neftel.csv', index=False)
    
    return comparison, sox2_comparison

def create_all_shaun2_style_figures(merged_data, neftel_data, comparison):
    """Create all figures exactly matching Shaun2 style"""
    
    print("Creating all Shaun2-style figures...")
    
    output_dir = Path('/Users/shaunie/Desktop/hu-glioblastoma/results/neftel/figures')
    output_dir.mkdir(exist_ok=True)
    
    # Figure 1: Cell Type Distribution (exactly like Shaun2)
    fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Focus on main cell types
    main_types = ['Astrocytic', 'Mesenchymal', 'Neural_Progenitor', 'Oligodendrocytic', 'Cycling', 'Endothelial']
    plot_data = merged_data[merged_data['base_cell_type'].isin(main_types)]
    
    # Plot 1: Total cells by type and dataset
    cell_summary = plot_data.groupby(['dataset', 'base_cell_type'])['total_cells'].sum().reset_index()
    pivot_cells = cell_summary.pivot(index='base_cell_type', columns='dataset', values='total_cells').fillna(0)
    
    pivot_cells.plot(kind='bar', ax=ax1, alpha=0.8, width=0.7)
    ax1.set_title('Cell Type Distribution Across Datasets\n(Neftel Methodology)', 
                  fontsize=14, fontweight='bold')
    ax1.set_xlabel('Cell Type')
    ax1.set_ylabel('Total Cells')
    ax1.legend(title='Dataset')
    ax1.tick_params(axis='x', rotation=45)
    
    # Add value labels
    for container in ax1.containers:
        ax1.bar_label(container, fmt='%d', rotation=90, fontsize=9)
    
    # Plot 2: SOX2 expression by cell type
    sox2_summary = plot_data.groupby(['dataset', 'base_cell_type']).agg({
        'SOX2_expressing': 'sum',
        'total_cells': 'sum'
    }).reset_index()
    sox2_summary['SOX2_percentage'] = (sox2_summary['SOX2_expressing'] / sox2_summary['total_cells']) * 100
    
    # Plot SOX2 percentages
    for i, dataset in enumerate(['GSE131928', 'Organoid']):
        dataset_data = sox2_summary[sox2_summary['dataset'] == dataset]
        positions = np.arange(len(dataset_data)) + (0.4 if dataset == 'Organoid' else -0.4)
        ax2.bar(positions, dataset_data['SOX2_percentage'], 
                width=0.35, label=dataset, alpha=0.8)
    
    ax2.set_title('SOX2 Expression by Cell Type\n(Neftel Methodology)', 
                  fontsize=14, fontweight='bold')
    ax2.set_xlabel('Cell Type')
    ax2.set_ylabel('SOX2+ Cells (%)')
    ax2.set_xticks(range(len(dataset_data)))
    ax2.set_xticklabels(dataset_data['base_cell_type'], rotation=45)
    ax2.legend(title='Dataset')
    
    plt.tight_layout()
    fig1.savefig(output_dir / 'cell_type_distribution.png', dpi=300, bbox_inches='tight')
    plt.close(fig1)
    
    # Figure 2: Four Population Comparison (exactly like Shaun2)
    fig2, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Use the comparison data for 4 main populations
    four_pops = ['Astrocytic', 'Mesenchymal', 'Neural_Progenitor', 'Oligodendrocytic']
    four_pop_comparison = comparison[comparison['base_cell_type'].isin(four_pops)]
    
    # Plot 1: Side-by-side percentage comparison
    pivot_pct = four_pop_comparison.pivot(index='base_cell_type', columns='dataset', 
                                          values='percentage_of_dataset').fillna(0)
    
    x = np.arange(len(pivot_pct))
    width = 0.35
    
    bars1 = ax1.bar(x - width/2, pivot_pct['GSE131928'], width, label='Patient (GSE131928)', alpha=0.8)
    bars2 = ax1.bar(x + width/2, pivot_pct['Organoid'], width, label='Organoid', alpha=0.8)
    
    ax1.set_title('Four Population Composition: Patient vs Organoid\n(Neftel Overlapping States)', 
                  fontweight='bold')
    ax1.set_xlabel('Cell Type')
    ax1.set_ylabel('Percentage of Dataset (%)')
    ax1.set_xticks(x)
    ax1.set_xticklabels(pivot_pct.index, rotation=45)
    ax1.legend()
    
    # Add percentage labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax1.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                        f'{height:.1f}%', ha='center', va='bottom', fontsize=9)
    
    # Plot 2: Correlation scatter
    patient_pcts = pivot_pct['GSE131928'].values
    organoid_pcts = pivot_pct['Organoid'].values
    
    ax2.scatter(patient_pcts, organoid_pcts, s=200, alpha=0.7, edgecolors='black')
    
    # Add correlation line
    if len(patient_pcts) > 1:
        correlation = np.corrcoef(patient_pcts, organoid_pcts)[0, 1]
        m, b = np.polyfit(patient_pcts, organoid_pcts, 1)
        ax2.plot(patient_pcts, m*patient_pcts + b, 'r--', alpha=0.7, linewidth=2)
        ax2.text(0.05, 0.95, f'R = {correlation:.3f}', transform=ax2.transAxes, 
                fontsize=12, bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
    
    ax2.set_xlabel('Patient Percentage (%)')
    ax2.set_ylabel('Organoid Percentage (%)')
    ax2.set_title('Patient-Organoid Correlation\n(Four Malignant Populations)', fontweight='bold')
    ax2.set_xlim(-5, max(patient_pcts) + 5)
    ax2.set_ylim(-5, max(organoid_pcts) + 5)
    
    # Add cell type labels
    for i, cell_type in enumerate(pivot_pct.index):
        ax2.annotate(cell_type, (patient_pcts[i], organoid_pcts[i]), 
                    xytext=(5, 5), textcoords='offset points', fontsize=10)
    
    # Plot 3: Cell counts comparison
    pivot_counts = four_pop_comparison.pivot(index='base_cell_type', columns='dataset', 
                                            values='total_cells').fillna(0)
    
    pivot_counts.plot(kind='bar', ax=ax3, alpha=0.8)
    ax3.set_title('Cell Counts by Type\n(Neftel Methodology)', fontweight='bold')
    ax3.set_xlabel('Cell Type')
    ax3.set_ylabel('Total Cells')
    ax3.legend(title='Dataset')
    ax3.tick_params(axis='x', rotation=45)
    
    # Plot 4: SOX2 comparison for 4 populations
    pivot_sox2 = four_pop_comparison.pivot(index='base_cell_type', columns='dataset', 
                                          values='SOX2_percentage').fillna(0)
    
    bars1 = ax4.bar(x - width/2, pivot_sox2['GSE131928'], width, label='Patient', alpha=0.8)
    bars2 = ax4.bar(x + width/2, pivot_sox2['Organoid'], width, label='Organoid', alpha=0.8)
    
    ax4.set_title('SOX2 Expression in Four Populations\n(Patient vs Organoid)', fontweight='bold')
    ax4.set_xlabel('Cell Type')
    ax4.set_ylabel('SOX2+ Percentage (%)')
    ax4.set_xticks(x)
    ax4.set_xticklabels(pivot_sox2.index, rotation=45)
    ax4.legend()
    
    plt.tight_layout()
    fig2.savefig(output_dir / 'four_population_comparison.png', dpi=300, bbox_inches='tight')
    plt.close(fig2)
    
    # Figure 3: Sample Heatmap (exactly like Shaun2)
    fig3, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 12))
    
    # Create heatmap data
    heatmap_data = merged_data.pivot_table(
        index='sample',
        columns='base_cell_type',
        values='percentage',
        fill_value=0
    )
    
    # Separate by dataset
    patient_samples = [s for s in heatmap_data.index if not s.startswith('Organoid')]
    organoid_samples = [s for s in heatmap_data.index if s.startswith('Organoid')]
    
    # Patient heatmap
    if patient_samples:
        patient_data = heatmap_data.loc[patient_samples]
        sns.heatmap(patient_data, annot=True, fmt='.1f', cmap='viridis', ax=ax1,
                   cbar_kws={'label': 'Percentage (%)'})
        ax1.set_title('Patient Sample Cell Type Composition\n(Neftel Overlapping States)', 
                     fontweight='bold')
        ax1.set_xlabel('Cell Type')
        ax1.set_ylabel('Patient Sample')
    
    # Organoid heatmap
    if organoid_samples:
        organoid_data = heatmap_data.loc[organoid_samples]
        sns.heatmap(organoid_data, annot=True, fmt='.1f', cmap='plasma', ax=ax2,
                   cbar_kws={'label': 'Percentage (%)'})
        ax2.set_title('Organoid Sample Cell Type Composition\n(Neftel Overlapping States)', 
                     fontweight='bold')
        ax2.set_xlabel('Cell Type')
        ax2.set_ylabel('Organoid Sample')
    
    plt.tight_layout()
    fig3.savefig(output_dir / 'sample_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close(fig3)
    
    print("✅ All Shaun2-style figures created successfully!")

def main():
    """Generate all Shaun2-style outputs using Neftel methodology"""
    
    print("🧬 Creating ALL Shaun2-Style Outputs with Neftel Methodology")
    print("=" * 70)
    
    # Load and prepare data
    print("\n1. LOADING AND PREPARING DATA")
    neftel_data, merged_data = load_and_prepare_neftel_data()
    print(f"Loaded {len(neftel_data)} SOX2-split entries")
    print(f"Merged to {len(merged_data)} cell type entries")
    
    # Create FourSample-FourPop CSVs
    print("\n2. CREATING FOURSAMPLE-FOURPOP CSVs")
    four_sample_data = create_four_sample_four_pop_csvs(merged_data)
    
    # Create CellCounts, Percentages, Methods
    print("\n3. CREATING CELLCOUNTS, PERCENTAGES, AND METHODS")
    cell_counts, percentages = create_cell_counts_percentages_methods(merged_data)
    
    # Create OrganoidVsPatient comparison
    print("\n4. CREATING ORGANOID VS PATIENT COMPARISONS")
    comparison, sox2_comparison = create_organoid_vs_patient_comparison(merged_data, neftel_data)
    
    # Create all figures
    print("\n5. CREATING ALL SHAUN2-STYLE FIGURES")
    create_all_shaun2_style_figures(merged_data, neftel_data, comparison)
    
    print("\n" + "=" * 70)
    print("✅ ALL SHAUN2-STYLE OUTPUTS CREATED WITH NEFTEL METHODOLOGY!")
    print("=" * 70)
    
    print("\n📁 Files created in /results/neftel/:")
    print("  • FourSample-FourPop-Neftel.csv")
    print("  • FourPop-Summary-Neftel.csv")
    print("  • CellCounts-Neftel.csv")
    print("  • Percentages-Neftel.csv")
    print("  • Methods-Neftel.csv")
    print("  • OrganoidVsPatient-FourPop-Neftel.csv")
    print("  • OrganoidVsPatient-SOX2Split-Neftel.csv")
    
    print("\n🎨 Figures created in /results/neftel/figures/:")
    print("  • cell_type_distribution.png")
    print("  • four_population_comparison.png")
    print("  • sample_heatmap.png")
    
    print("\n📌 KEY FEATURES:")
    print("  ✅ Exact same format as Shaun2")
    print("  ✅ NO immune cells in organoids (removed)")
    print("  ✅ Overlapping cell states (Neftel methodology)")
    print("  ✅ SOX2+/SOX2- analysis included")
    print("  ✅ 4 main malignant populations focus")

if __name__ == "__main__":
    main()