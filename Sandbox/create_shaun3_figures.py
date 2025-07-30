#!/usr/bin/env python3
"""
Create figures for Shaun3 exclusive cell type categories analysis
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

# Set style
plt.style.use('default')
sns.set_palette("husl")

def load_shaun3_data():
    """Load Shaun3 exclusive category data"""
    per_sample_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun3/Sox2-Cell-Populations-PerSample-Exclusive.csv'
    agg_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun3/Sox2-Cell-Populations-Exclusive.csv'
    
    per_sample = pd.read_csv(per_sample_path)
    aggregated = pd.read_csv(agg_path)
    
    return per_sample, aggregated

def create_stacked_bar_chart(aggregated_data):
    """Create stacked bar chart showing cell type distributions"""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Prepare data for plotting
    organoid_data = aggregated_data[aggregated_data['Dataset'] == 'Organoid'].copy()
    patient_data = aggregated_data[aggregated_data['Dataset'] == 'Patient'].copy()
    
    # Sort by total cell count for better visualization
    organoid_data = organoid_data.sort_values('Total_Population_Cells', ascending=True)
    patient_data = patient_data.sort_values('Total_Population_Cells', ascending=True)
    
    # Define colors for cell types
    cell_type_colors = {
        'Endothelial': '#FF6B6B',
        'Cycling_Oligodendrocytic': '#4ECDC4', 
        'Cycling_Astrocytic': '#45B7D1',
        'Cycling_Mesenchymal': '#96CEB4',
        'Cycling_Neural_Progenitor': '#FFEAA7',
        'Oligodendrocytic': '#DDA0DD',
        'Astrocytic': '#87CEEB',
        'Neural_Progenitor': '#F0E68C',
        'Mesenchymal': '#98D8C8'
    }
    
    # Plot 1: Organoid data
    y_pos = range(len(organoid_data))
    cumulative = np.zeros(len(organoid_data))
    
    for cell_type in organoid_data['Population']:
        values = organoid_data[organoid_data['Population'] == cell_type]['Percent_of_Dataset'].values
        color = cell_type_colors.get(cell_type, '#CCCCCC')
        ax1.barh(y_pos, values, left=cumulative, 
                label=cell_type, color=color, alpha=0.8)
        cumulative += values
    
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(organoid_data['Population'])
    ax1.set_xlabel('Percentage of Dataset')
    ax1.set_title('Organoid Cell Type Distribution\n(Mutually Exclusive Categories)', fontsize=14, fontweight='bold')
    ax1.set_xlim(0, 100)
    
    # Plot 2: Patient data  
    y_pos_patient = range(len(patient_data))
    cumulative_patient = np.zeros(len(patient_data))
    
    for cell_type in patient_data['Population']:
        values = patient_data[patient_data['Population'] == cell_type]['Percent_of_Dataset'].values
        color = cell_type_colors.get(cell_type, '#CCCCCC')
        ax2.barh(y_pos_patient, values, left=cumulative_patient,
                label=cell_type, color=color, alpha=0.8)
        cumulative_patient += values
    
    ax2.set_yticks(y_pos_patient)
    ax2.set_yticklabels(patient_data['Population'])
    ax2.set_xlabel('Percentage of Dataset')
    ax2.set_title('Patient Cell Type Distribution\n(Mutually Exclusive Categories)', fontsize=14, fontweight='bold')
    ax2.set_xlim(0, 100)
    
    # Add legend
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc='center right', bbox_to_anchor=(1.15, 0.5))
    
    plt.tight_layout()
    return fig

def create_sox2_heatmap(aggregated_data):
    """Create heatmap showing SOX2-negative percentages"""
    
    # Prepare data for heatmap
    heatmap_data = aggregated_data.pivot(index='Population', columns='Dataset', 
                                        values='Percent_SOX2_Negative_in_Population')
    
    # Reorder for better visualization (group cycling types together)
    order = ['Endothelial', 'Oligodendrocytic', 'Cycling_Oligodendrocytic',
             'Mesenchymal', 'Cycling_Mesenchymal', 'Astrocytic', 'Cycling_Astrocytic', 
             'Neural_Progenitor', 'Cycling_Neural_Progenitor']
    
    heatmap_data = heatmap_data.reindex(order)
    
    fig, ax = plt.subplots(figsize=(8, 10))
    
    # Create heatmap
    sns.heatmap(heatmap_data, annot=True, fmt='.1f', cmap='RdYlBu_r',
                center=50, ax=ax, cbar_kws={'label': '% SOX2-negative'})
    
    ax.set_title('SOX2-Negative Percentages by Cell Type\n(Exclusive Categories)', 
                fontsize=14, fontweight='bold')
    ax.set_xlabel('Dataset')
    ax.set_ylabel('Cell Type')
    
    plt.tight_layout()
    return fig

def create_cycling_comparison(aggregated_data):
    """Compare cycling vs non-cycling cells"""
    
    # Separate cycling and non-cycling data
    cycling_data = aggregated_data[aggregated_data['Population'].str.contains('Cycling_')]
    non_cycling_data = aggregated_data[~aggregated_data['Population'].str.contains('Cycling_') & 
                                     (aggregated_data['Population'] != 'Endothelial')]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Cell counts comparison
    cycling_summary = cycling_data.groupby('Dataset')['Total_Population_Cells'].sum()
    non_cycling_summary = non_cycling_data.groupby('Dataset')['Total_Population_Cells'].sum()
    
    x = np.arange(len(cycling_summary.index))
    width = 0.35
    
    bars1 = ax1.bar(x - width/2, cycling_summary.values, width, label='Cycling', color='#FF6B6B', alpha=0.7)
    bars2 = ax1.bar(x + width/2, non_cycling_summary.values, width, label='Non-Cycling', color='#4ECDC4', alpha=0.7)
    
    ax1.set_xlabel('Dataset')
    ax1.set_ylabel('Total Cell Count')
    ax1.set_title('Cycling vs Non-Cycling Cell Counts', fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(cycling_summary.index)
    ax1.legend()
    
    # Add value labels on bars
    for bar in bars1:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 50,
                f'{int(height):,}', ha='center', va='bottom')
    for bar in bars2:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 50,
                f'{int(height):,}', ha='center', va='bottom')
    
    # Plot 2: SOX2-negative percentages comparison
    cycling_sox2 = cycling_data.groupby('Dataset').apply(
        lambda x: (x['SOX2_Negative_in_Population'].sum() / x['Total_Population_Cells'].sum()) * 100)
    non_cycling_sox2 = non_cycling_data.groupby('Dataset').apply(
        lambda x: (x['SOX2_Negative_in_Population'].sum() / x['Total_Population_Cells'].sum()) * 100)
    
    bars3 = ax2.bar(x - width/2, cycling_sox2.values, width, label='Cycling', color='#FF6B6B', alpha=0.7)
    bars4 = ax2.bar(x + width/2, non_cycling_sox2.values, width, label='Non-Cycling', color='#4ECDC4', alpha=0.7)
    
    ax2.set_xlabel('Dataset')
    ax2.set_ylabel('% SOX2-negative')
    ax2.set_title('SOX2-negative: Cycling vs Non-Cycling', fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(cycling_sox2.index)
    ax2.legend()
    
    # Add value labels
    for bar in bars3:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 1,
                f'{height:.1f}%', ha='center', va='bottom')
    for bar in bars4:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 1,
                f'{height:.1f}%', ha='center', va='bottom')
    
    plt.tight_layout()
    return fig

def create_sample_level_plot(per_sample_data):
    """Create sample-level visualization"""
    
    # Calculate total cells per sample for normalization
    sample_totals = per_sample_data.groupby(['dataset', 'sample'])['total_cells'].sum()
    
    # Add percentage column
    per_sample_data = per_sample_data.copy()
    per_sample_data['sample_percentage'] = 0.0
    
    for (dataset, sample), total in sample_totals.items():
        mask = (per_sample_data['dataset'] == dataset) & (per_sample_data['sample'] == sample)
        per_sample_data.loc[mask, 'sample_percentage'] = (per_sample_data.loc[mask, 'total_cells'] / total) * 100
    
    # Focus on samples with reasonable cell counts (>100 cells)
    sample_sizes = per_sample_data.groupby(['dataset', 'sample'])['total_cells'].sum()
    large_samples = sample_sizes[sample_sizes > 100].index
    
    filtered_data = per_sample_data.set_index(['dataset', 'sample']).loc[large_samples].reset_index()
    
    # Create pivot table for heatmap
    pivot_data = filtered_data.pivot_table(index=['dataset', 'sample'], 
                                          columns='cell_type', 
                                          values='sample_percentage', 
                                          fill_value=0)
    
    fig, ax = plt.subplots(figsize=(12, max(8, len(pivot_data) * 0.3)))
    
    # Create heatmap
    sns.heatmap(pivot_data, annot=True, fmt='.1f', cmap='viridis', ax=ax,
                cbar_kws={'label': '% of Sample'})
    
    ax.set_title('Cell Type Distribution Across Samples\n(Mutually Exclusive Categories)', 
                fontsize=14, fontweight='bold')
    ax.set_xlabel('Cell Type')
    ax.set_ylabel('Sample')
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    
    plt.tight_layout()
    return fig

def main():
    """Generate all Shaun3 figures"""
    
    print("Creating Shaun3 figures for exclusive cell type categories...")
    
    # Load data
    per_sample_data, aggregated_data = load_shaun3_data()
    print(f"Loaded data: {len(per_sample_data)} per-sample rows, {len(aggregated_data)} aggregated rows")
    
    # Create output directory
    output_dir = Path('/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun3/figures')
    output_dir.mkdir(exist_ok=True)
    
    # Generate figures
    print("Creating stacked bar chart...")
    fig1 = create_stacked_bar_chart(aggregated_data)
    fig1.savefig(output_dir / 'cell_type_distribution_exclusive.png', dpi=300, bbox_inches='tight')
    plt.close(fig1)
    
    print("Creating SOX2 heatmap...")
    fig2 = create_sox2_heatmap(aggregated_data)
    fig2.savefig(output_dir / 'sox2_heatmap_exclusive.png', dpi=300, bbox_inches='tight')
    plt.close(fig2)
    
    print("Creating cycling comparison...")
    fig3 = create_cycling_comparison(aggregated_data)
    fig3.savefig(output_dir / 'cycling_vs_noncycling_comparison.png', dpi=300, bbox_inches='tight')
    plt.close(fig3)
    
    print("Creating sample-level plot...")
    fig4 = create_sample_level_plot(per_sample_data)
    fig4.savefig(output_dir / 'sample_level_distribution.png', dpi=300, bbox_inches='tight')
    plt.close(fig4)
    
    print(f"\n✅ All figures saved to: {output_dir}")
    print("Generated figures:")
    print("  1. cell_type_distribution_exclusive.png - Overall distribution comparison")
    print("  2. sox2_heatmap_exclusive.png - SOX2-negative percentages heatmap")
    print("  3. cycling_vs_noncycling_comparison.png - Cycling cell analysis")
    print("  4. sample_level_distribution.png - Per-sample breakdown")

if __name__ == "__main__":
    main()