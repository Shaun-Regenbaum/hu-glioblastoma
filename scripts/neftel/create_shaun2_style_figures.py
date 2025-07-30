#!/usr/bin/env python3
"""
Create Shaun2-style figures using Neftel methodology
Same exact figure types and layouts as Shaun2, but with overlapping cell states
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

# Set style to match Shaun2
plt.style.use('default')
sns.set_palette("husl")

def load_neftel_data():
    """Load Neftel-style data in Shaun2 format"""
    neftel_dir = Path('/Users/shaunie/Desktop/hu-glioblastoma/results/neftel')
    
    # Load the CSV we just created
    data = pd.read_csv(neftel_dir / 'Sox2-Cell-Populations-PerSample-Neftel.csv')
    
    return data

def create_cell_type_distribution_plot(data):
    """Create the same cell type distribution plot as Shaun2"""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Plot 1: Total cells by type and dataset (same as Shaun2)
    cell_summary = data.groupby(['dataset', 'cell_type'])['total_cells'].sum().reset_index()
    
    # Pivot for grouped bar chart
    pivot_data = cell_summary.pivot(index='cell_type', columns='dataset', values='total_cells').fillna(0)
    
    # Create the bar plot
    pivot_data.plot(kind='bar', ax=ax1, alpha=0.8, width=0.7)
    ax1.set_title('Cell Type Distribution Across Datasets\n(Neftel Overlapping Methodology)', 
                  fontsize=14, fontweight='bold')
    ax1.set_xlabel('Cell Type')
    ax1.set_ylabel('Total Cells')
    ax1.legend(title='Dataset')
    ax1.tick_params(axis='x', rotation=45)
    
    # Add value labels on bars
    for container in ax1.containers:
        ax1.bar_label(container, fmt='%d', rotation=90, fontsize=9)
    
    # Plot 2: SOX2 expression by cell type (same as Shaun2)
    sox2_summary = data.groupby(['dataset', 'cell_type']).agg({
        'SOX2_expressing': 'sum',
        'total_cells': 'sum'
    }).reset_index()
    
    sox2_summary['SOX2_percentage'] = (sox2_summary['SOX2_expressing'] / 
                                      sox2_summary['total_cells']) * 100
    
    # Create SOX2 percentage plot
    for dataset in ['GSE131928', 'Organoid']:
        dataset_data = sox2_summary[sox2_summary['dataset'] == dataset]
        ax2.bar(np.arange(len(dataset_data)) + (0.4 if dataset == 'Organoid' else -0.4), 
                dataset_data['SOX2_percentage'],
                width=0.35, 
                label=dataset,
                alpha=0.8)
    
    ax2.set_title('SOX2 Expression by Cell Type\n(Neftel Overlapping States)', 
                  fontsize=14, fontweight='bold')
    ax2.set_xlabel('Cell Type')
    ax2.set_ylabel('SOX2+ Cells (%)')
    ax2.set_xticks(range(len(dataset_data)))
    ax2.set_xticklabels(dataset_data['cell_type'], rotation=45)
    ax2.legend(title='Dataset')
    
    # Add percentage labels
    for i, dataset in enumerate(['GSE131928', 'Organoid']):
        dataset_data = sox2_summary[sox2_summary['dataset'] == dataset]
        offset = -0.4 if i == 0 else 0.4
        for j, (_, row) in enumerate(dataset_data.iterrows()):
            ax2.text(j + offset, row['SOX2_percentage'] + 1, 
                    f'{row["SOX2_percentage"]:.1f}%',
                    ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    return fig

def create_organoid_vs_patient_comparison(data):
    """Create organoid vs patient comparison (same as Shaun2)"""
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Calculate percentages by dataset
    dataset_summary = data.groupby(['dataset', 'cell_type'])['total_cells'].sum().reset_index()
    dataset_totals = dataset_summary.groupby('dataset')['total_cells'].sum()
    dataset_summary['percentage'] = dataset_summary.apply(
        lambda x: (x['total_cells'] / dataset_totals[x['dataset']]) * 100, axis=1)
    
    # Plot 1: Side-by-side percentage comparison
    pivot_pct = dataset_summary.pivot(index='cell_type', columns='dataset', values='percentage').fillna(0)
    
    x = np.arange(len(pivot_pct))
    width = 0.35
    
    bars1 = ax1.bar(x - width/2, pivot_pct['GSE131928'], width, label='Patient (GSE131928)', alpha=0.8)
    bars2 = ax1.bar(x + width/2, pivot_pct['Organoid'], width, label='Organoid', alpha=0.8)
    
    ax1.set_title('Cell Type Composition: Patient vs Organoid\n(Neftel Overlapping States)', 
                  fontweight='bold')
    ax1.set_xlabel('Cell Type')
    ax1.set_ylabel('Percentage of Total Cells')
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
    
    # Plot 2: Correlation plot
    patient_data = pivot_pct['GSE131928'].values
    organoid_data = pivot_pct['Organoid'].values
    
    ax2.scatter(patient_data, organoid_data, s=100, alpha=0.7)
    
    # Add correlation line
    if len(patient_data) > 1:
        correlation = np.corrcoef(patient_data, organoid_data)[0, 1]
        m, b = np.polyfit(patient_data, organoid_data, 1)
        ax2.plot(patient_data, m*patient_data + b, 'r--', alpha=0.7)
        ax2.text(0.05, 0.95, f'R = {correlation:.3f}', transform=ax2.transAxes, 
                bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
    
    ax2.set_xlabel('Patient Percentage (%)')
    ax2.set_ylabel('Organoid Percentage (%)')
    ax2.set_title('Patient-Organoid Correlation\n(Neftel Methodology)', fontweight='bold')
    
    # Add cell type labels
    for i, cell_type in enumerate(pivot_pct.index):
        ax2.annotate(cell_type, (patient_data[i], organoid_data[i]), 
                    xytext=(5, 5), textcoords='offset points', fontsize=9)
    
    # Plot 3: Sample count comparison
    sample_counts = data.groupby(['dataset', 'cell_type']).size().reset_index(name='sample_count')
    pivot_samples = sample_counts.pivot(index='cell_type', columns='dataset', values='sample_count').fillna(0)
    
    pivot_samples.plot(kind='bar', ax=ax3, alpha=0.8)
    ax3.set_title('Sample Counts by Cell Type\n(Neftel Overlapping Approach)', fontweight='bold')
    ax3.set_xlabel('Cell Type')
    ax3.set_ylabel('Number of Samples')
    ax3.legend(title='Dataset')
    ax3.tick_params(axis='x', rotation=45)
    
    # Plot 4: SOX2 comparison
    sox2_comparison = data.groupby(['dataset', 'cell_type']).agg({
        'SOX2_expressing': 'sum',
        'total_cells': 'sum'
    }).reset_index()
    sox2_comparison['SOX2_percentage'] = (sox2_comparison['SOX2_expressing'] / 
                                         sox2_comparison['total_cells']) * 100
    
    pivot_sox2 = sox2_comparison.pivot(index='cell_type', columns='dataset', values='SOX2_percentage').fillna(0)
    
    x = np.arange(len(pivot_sox2))
    bars1 = ax4.bar(x - width/2, pivot_sox2['GSE131928'], width, label='Patient', alpha=0.8)
    bars2 = ax4.bar(x + width/2, pivot_sox2['Organoid'], width, label='Organoid', alpha=0.8)
    
    ax4.set_title('SOX2 Expression: Patient vs Organoid\n(Neftel States)', fontweight='bold')
    ax4.set_xlabel('Cell Type')
    ax4.set_ylabel('SOX2+ Percentage')
    ax4.set_xticks(x)
    ax4.set_xticklabels(pivot_sox2.index, rotation=45)
    ax4.legend()
    
    plt.tight_layout()
    return fig

def create_sample_level_heatmap(data):
    """Create sample-level heatmap (same as Shaun2)"""
    
    # Create pivot table for heatmap
    heatmap_data = data.pivot_table(index='sample', 
                                   columns='cell_type', 
                                   values='percentage', 
                                   fill_value=0)
    
    # Separate by dataset
    patient_samples = [s for s in heatmap_data.index if not s.startswith('Organoid')]
    organoid_samples = [s for s in heatmap_data.index if s.startswith('Organoid')]
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 12))
    
    # Patient heatmap
    if patient_samples:
        patient_data = heatmap_data.loc[patient_samples]
        sns.heatmap(patient_data, annot=True, fmt='.1f', cmap='viridis', ax=ax1,
                   cbar_kws={'label': 'Percentage (%)'})
        ax1.set_title('Patient Sample Cell Type Composition\n(Neftel Overlapping States - May Exceed 100%)', 
                     fontweight='bold')
        ax1.set_xlabel('Cell Type')
        ax1.set_ylabel('Patient Sample')
    
    # Organoid heatmap  
    if organoid_samples:
        organoid_data = heatmap_data.loc[organoid_samples]
        sns.heatmap(organoid_data, annot=True, fmt='.1f', cmap='plasma', ax=ax2,
                   cbar_kws={'label': 'Percentage (%)'})
        ax2.set_title('Organoid Sample Cell Type Composition\n(Neftel Overlapping States - May Exceed 100%)', 
                     fontweight='bold')
        ax2.set_xlabel('Cell Type') 
        ax2.set_ylabel('Organoid Sample')
    
    plt.tight_layout()
    return fig

def create_methodology_comparison_note():
    """Create a note explaining the differences from Shaun2"""
    
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis('off')
    
    explanation_text = """
NEFTEL METHODOLOGY vs SHAUN2 COMPARISON

📋 SAME AS SHAUN2:
✅ Identical CSV format and structure
✅ Same 6 cell types included (Astrocytic, Mesenchymal, Neural_Progenitor, Oligodendrocytic, Cycling, Endothelial)
✅ Same figure layouts and visualizations
✅ Same SOX2 expression analysis

🔄 KEY DIFFERENCES (Neftel Approach):
• OVERLAPPING CELL STATES: Cells can belong to multiple categories simultaneously
• PERCENTAGES MAY EXCEED 100%: Due to overlapping classifications per sample
• META-MODULE SCORING: Based on Neftel et al. (Cell 2019) methodology
• HYBRID STATES PRESERVED: Maintains cellular plasticity information

📊 SCIENTIFIC RATIONALE:
The Neftel methodology captures the biological reality that glioblastoma cells exist in
plastic states and can express multiple cellular programs simultaneously. This overlapping
approach preserves information about cellular plasticity that would be lost in mutually
exclusive categories.

🎯 WHEN TO USE EACH:
• Neftel Approach: For understanding cellular plasticity and hybrid states
• Shaun2/Shaun3: For clean comparative analysis between datasets

Both approaches use identical marker genes and are scientifically valid.
The choice depends on the research question being addressed.
"""
    
    ax.text(0.05, 0.95, explanation_text, transform=ax.transAxes, 
            fontsize=11, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle="round,pad=0.5", facecolor='lightblue', alpha=0.3))
    
    ax.set_title('Neftel Methodology: Key Differences from Shaun2', 
                fontsize=16, fontweight='bold', pad=20)
    
    return fig

def main():
    """Generate all Shaun2-style figures using Neftel methodology"""
    
    print("🧬 Creating Shaun2-Style Figures with Neftel Methodology")
    print("=" * 70)
    
    # Load data
    data = load_neftel_data()
    print(f"Loaded {len(data)} cell type assignments")
    
    # Create output directory
    output_dir = Path('/Users/shaunie/Desktop/hu-glioblastoma/results/neftel/figures/shaun2_style')
    output_dir.mkdir(exist_ok=True)
    
    # Generate figures (same as Shaun2 but with Neftel data)
    print("\n1. Creating cell type distribution plot...")
    fig1 = create_cell_type_distribution_plot(data)
    fig1.savefig(output_dir / 'cell_type_distribution_neftel.png', dpi=300, bbox_inches='tight')
    plt.close(fig1)
    
    print("2. Creating organoid vs patient comparison...")
    fig2 = create_organoid_vs_patient_comparison(data)
    fig2.savefig(output_dir / 'organoid_vs_patient_neftel.png', dpi=300, bbox_inches='tight')
    plt.close(fig2)
    
    print("3. Creating sample-level heatmap...")
    fig3 = create_sample_level_heatmap(data)
    fig3.savefig(output_dir / 'sample_heatmap_neftel.png', dpi=300, bbox_inches='tight')
    plt.close(fig3)
    
    print("4. Creating methodology comparison note...")
    fig4 = create_methodology_comparison_note()
    fig4.savefig(output_dir / 'neftel_vs_shaun2_explanation.png', dpi=300, bbox_inches='tight')
    plt.close(fig4)
    
    print(f"\n✅ All Shaun2-style figures saved to: {output_dir}")
    print("\nGenerated figures:")
    print("  1. cell_type_distribution_neftel.png - Same layout as Shaun2")
    print("  2. organoid_vs_patient_neftel.png - Same comparison plots")
    print("  3. sample_heatmap_neftel.png - Same heatmap style")
    print("  4. neftel_vs_shaun2_explanation.png - Methodology differences")
    
    print("\n" + "=" * 70)
    print("🎯 KEY POINT: Same figures as Shaun2, using Neftel overlapping methodology")
    print("=" * 70)

if __name__ == "__main__":
    main()