#!/usr/bin/env python3
"""
Create visualization figures for per-sample data and organoid vs patient comparisons
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os

def setup_plotting():
    """Set up plotting parameters"""
    plt.style.use('default')
    sns.set_palette("Set2")
    plt.rcParams['figure.figsize'] = (12, 8)
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10

def create_cell_type_distribution_figure():
    """Create figure showing cell type distributions per sample"""
    
    print("Creating cell type distribution figure...")
    
    # Load the per-sample data
    data_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/OrganoidVsPatient-FourPop-PerSample-REAL.csv'
    df = pd.read_csv(data_path)
    
    # Create figure with subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Cell Type Distributions: Organoids vs Patient Tumors', fontsize=16, fontweight='bold')
    
    # Define colors for cell types
    colors = {
        'Mesenchymal': '#ff7f0e',
        'Astrocytic': '#2ca02c', 
        'Neural_Progenitor': '#1f77b4',
        'Oligodendrocytic': '#d62728'
    }
    
    # 1. Organoid samples stacked bar chart
    organoid_data = df[df['dataset'] == 'Organoid']
    organoid_pivot = organoid_data.pivot(index='sample', columns='cell_type', values='percentage')
    
    organoid_pivot.plot(kind='bar', stacked=True, ax=ax1, color=[colors[col] for col in organoid_pivot.columns])
    ax1.set_title('Organoid Samples - Cell Type Percentages', fontweight='bold')
    ax1.set_xlabel('Organoid Sample')
    ax1.set_ylabel('Percentage (%)')
    ax1.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.tick_params(axis='x', rotation=45)
    
    # 2. Patient tumor samples (top 10 by cell count)
    patient_data = df[df['dataset'] == 'GSE131928']
    top_tumors = patient_data.groupby('sample')['total_cells_in_sample'].first().nlargest(10)
    patient_subset = patient_data[patient_data['sample'].isin(top_tumors.index)]
    patient_pivot = patient_subset.pivot(index='sample', columns='cell_type', values='percentage')
    
    patient_pivot.plot(kind='bar', stacked=True, ax=ax2, color=[colors[col] for col in patient_pivot.columns])
    ax2.set_title('Patient Tumors (Top 10) - Cell Type Percentages', fontweight='bold')
    ax2.set_xlabel('Patient Tumor')
    ax2.set_ylabel('Percentage (%)')
    ax2.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax2.tick_params(axis='x', rotation=45)
    
    # 3. Box plot comparison between datasets
    plot_data = []
    for _, row in df.iterrows():
        plot_data.append({
            'Dataset': 'Organoid' if row['dataset'] == 'Organoid' else 'Patient',
            'Cell_Type': row['cell_type'],
            'Percentage': row['percentage']
        })
    
    plot_df = pd.DataFrame(plot_data)
    
    sns.boxplot(data=plot_df, x='Cell_Type', y='Percentage', hue='Dataset', ax=ax3)
    ax3.set_title('Cell Type Distribution Comparison', fontweight='bold')
    ax3.set_xlabel('Cell Type')
    ax3.set_ylabel('Percentage (%)')
    ax3.tick_params(axis='x', rotation=45)
    ax3.legend(title='Dataset')
    
    # 4. Sample size comparison
    sample_counts = df.groupby(['dataset', 'sample'])['total_cells_in_sample'].first().reset_index()
    sample_counts['Dataset'] = sample_counts['dataset'].map({'Organoid': 'Organoid', 'GSE131928': 'Patient'})
    
    sns.boxplot(data=sample_counts, x='Dataset', y='total_cells_in_sample', ax=ax4)
    ax4.set_title('Sample Size Distribution', fontweight='bold')
    ax4.set_xlabel('Dataset')
    ax4.set_ylabel('Total Cells per Sample')
    ax4.set_yscale('log')
    
    plt.tight_layout()
    
    # Save figure
    os.makedirs('/Users/shaunie/Desktop/hu-glioblastoma/figures', exist_ok=True)
    output_path = '/Users/shaunie/Desktop/hu-glioblastoma/figures/cell_type_distributions_per_sample.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved figure: {output_path}")
    
    plt.show()
    return fig

def create_sox2_analysis_figure():
    """Create figure showing Sox2 analysis per sample"""
    
    print("Creating Sox2 analysis figure...")
    
    # Load Sox2 data
    sox2_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/Sox2-Cell-Populations-PerSample-REAL.csv'
    df = pd.read_csv(sox2_path)
    
    # Create figure
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Sox2 Expression Analysis: Organoids vs Patient Tumors', fontsize=16, fontweight='bold')
    
    # Define colors
    colors = {
        'Mesenchymal': '#ff7f0e',
        'Astrocytic': '#2ca02c', 
        'Neural_Progenitor': '#1f77b4',
        'Oligodendrocytic': '#d62728'
    }
    
    # 1. Organoid Sox2 percentages
    organoid_sox2 = df[df['dataset'] == 'Organoid']
    organoid_pivot = organoid_sox2.pivot(index='sample', columns='cell_type', values='percent_SOX2_negative')
    
    organoid_pivot.plot(kind='bar', ax=ax1, color=[colors[col] for col in organoid_pivot.columns])
    ax1.set_title('Organoids - % Sox2-Negative by Cell Type', fontweight='bold')
    ax1.set_xlabel('Organoid Sample')
    ax1.set_ylabel('% Sox2-Negative')
    ax1.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.tick_params(axis='x', rotation=45)
    ax1.set_ylim(0, 100)
    
    # 2. Patient tumor Sox2 percentages (top 8 tumors)
    patient_sox2 = df[df['dataset'] == 'GSE131928']
    top_tumors_sox2 = patient_sox2.groupby('sample')['total_cells'].sum().nlargest(8)
    patient_subset_sox2 = patient_sox2[patient_sox2['sample'].isin(top_tumors_sox2.index)]
    patient_pivot_sox2 = patient_subset_sox2.pivot(index='sample', columns='cell_type', values='percent_SOX2_negative')
    
    patient_pivot_sox2.plot(kind='bar', ax=ax2, color=[colors[col] for col in patient_pivot_sox2.columns])
    ax2.set_title('Patient Tumors (Top 8) - % Sox2-Negative by Cell Type', fontweight='bold')
    ax2.set_xlabel('Patient Tumor')
    ax2.set_ylabel('% Sox2-Negative')
    ax2.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax2.tick_params(axis='x', rotation=45)
    ax2.set_ylim(0, 100)
    
    # 3. Box plot comparison of Sox2-negative percentages
    plot_data_sox2 = []
    for _, row in df.iterrows():
        plot_data_sox2.append({
            'Dataset': 'Organoid' if row['dataset'] == 'Organoid' else 'Patient',
            'Cell_Type': row['cell_type'],
            'Sox2_Negative_Percent': row['percent_SOX2_negative']
        })
    
    plot_df_sox2 = pd.DataFrame(plot_data_sox2)
    
    sns.boxplot(data=plot_df_sox2, x='Cell_Type', y='Sox2_Negative_Percent', hue='Dataset', ax=ax3)
    ax3.set_title('Sox2-Negative % Comparison by Cell Type', fontweight='bold')
    ax3.set_xlabel('Cell Type')
    ax3.set_ylabel('% Sox2-Negative')
    ax3.tick_params(axis='x', rotation=45)
    ax3.legend(title='Dataset')
    ax3.set_ylim(0, 100)
    
    # 4. Scatter plot: Total cells vs Sox2-negative percentage
    scatter_data = df.copy()
    scatter_data['Dataset'] = scatter_data['dataset'].map({'Organoid': 'Organoid', 'GSE131928': 'Patient'})
    
    for cell_type in scatter_data['cell_type'].unique():
        subset = scatter_data[scatter_data['cell_type'] == cell_type]
        ax4.scatter(subset['total_cells'], subset['percent_SOX2_negative'], 
                   label=cell_type, color=colors[cell_type], alpha=0.7, s=50)
    
    ax4.set_title('Sample Size vs Sox2-Negative %', fontweight='bold')
    ax4.set_xlabel('Total Cells in Sample')
    ax4.set_ylabel('% Sox2-Negative')  
    ax4.set_xscale('log')
    ax4.legend(title='Cell Type')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    output_path = '/Users/shaunie/Desktop/hu-glioblastoma/figures/sox2_analysis_per_sample.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved figure: {output_path}")
    
    plt.show()
    return fig

def create_summary_comparison_figure():
    """Create summary comparison figure between organoids and patients"""
    
    print("Creating summary comparison figure...")
    
    # Load both datasets
    cell_type_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/OrganoidVsPatient-FourPop-PerSample-REAL.csv'
    agg_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/OrganoidVsPatient-FourPop.csv'
    
    df_detailed = pd.read_csv(cell_type_path)
    df_agg = pd.read_csv(agg_path)
    
    # Create figure
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Summary: Organoids vs Patient Tumors Comparison', fontsize=16, fontweight='bold')
    
    # 1. Aggregated comparison bar chart
    cell_types = df_agg['cell_type_simplified'].values
    organoid_pcts = df_agg['Organoids_Percentage'].values
    patient_pcts = df_agg['GSE131928_Percentage'].values
    
    x = np.arange(len(cell_types))
    width = 0.35
    
    ax1.bar(x - width/2, organoid_pcts, width, label='Organoids', alpha=0.8, color='#2ca02c')
    ax1.bar(x + width/2, patient_pcts, width, label='Patients', alpha=0.8, color='#ff7f0e')
    
    ax1.set_title('Aggregated Cell Type Percentages', fontweight='bold')
    ax1.set_xlabel('Cell Type')
    ax1.set_ylabel('Percentage (%)')
    ax1.set_xticks(x)
    ax1.set_xticklabels(cell_types, rotation=45)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Sample count comparison
    sample_counts = df_detailed.groupby('dataset')['sample'].nunique()
    
    bars = ax2.bar(['Organoids', 'Patient Tumors'], 
                   [sample_counts['Organoid'], sample_counts['GSE131928']], 
                   color=['#2ca02c', '#ff7f0e'], alpha=0.8)
    
    ax2.set_title('Number of Samples', fontweight='bold')
    ax2.set_ylabel('Number of Samples')
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{int(height)}', ha='center', va='bottom', fontweight='bold')
    
    # 3. Total cell count comparison
    total_cells = df_detailed.groupby('dataset')['count'].sum()
    
    bars = ax3.bar(['Organoids', 'Patient Tumors'], 
                   [total_cells['Organoid'], total_cells['GSE131928']], 
                   color=['#2ca02c', '#ff7f0e'], alpha=0.8)
    
    ax3.set_title('Total Cell Counts', fontweight='bold')
    ax3.set_ylabel('Total Cells')
    ax3.set_yscale('log')
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                f'{int(height):,}', ha='center', va='bottom', fontweight='bold')
    
    # 4. Variability comparison (coefficient of variation)
    variability_data = []
    for dataset in ['Organoid', 'GSE131928']:
        dataset_data = df_detailed[df_detailed['dataset'] == dataset]
        for cell_type in dataset_data['cell_type'].unique():
            ct_data = dataset_data[dataset_data['cell_type'] == cell_type]
            cv = ct_data['percentage'].std() / ct_data['percentage'].mean() * 100
            variability_data.append({
                'Dataset': 'Organoid' if dataset == 'Organoid' else 'Patient',
                'Cell_Type': cell_type,
                'CV': cv
            })
    
    var_df = pd.DataFrame(variability_data)
    var_pivot = var_df.pivot(index='Cell_Type', columns='Dataset', values='CV')
    
    var_pivot.plot(kind='bar', ax=ax4, color=['#2ca02c', '#ff7f0e'], alpha=0.8)
    ax4.set_title('Cell Type Variability (Coefficient of Variation)', fontweight='bold')
    ax4.set_xlabel('Cell Type')
    ax4.set_ylabel('Coefficient of Variation (%)')
    ax4.tick_params(axis='x', rotation=45)
    ax4.legend(title='Dataset')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    output_path = '/Users/shaunie/Desktop/hu-glioblastoma/figures/summary_organoids_vs_patients.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved figure: {output_path}")
    
    plt.show()
    return fig

def create_heatmap_figure():
    """Create heatmap showing cell type percentages across all samples"""
    
    print("Creating heatmap figure...")
    
    # Load data
    data_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/OrganoidVsPatient-FourPop-PerSample-REAL.csv'
    df = pd.read_csv(data_path)
    
    # Create pivot table for heatmap
    heatmap_data = df.pivot(index='sample', columns='cell_type', values='percentage')
    
    # Sort samples: organoids first, then patients by cell count
    organoid_samples = df[df['dataset'] == 'Organoid']['sample'].unique()
    patient_samples = df[df['dataset'] == 'GSE131928'].groupby('sample')['total_cells_in_sample'].first().sort_values(ascending=False).index
    
    ordered_samples = list(organoid_samples) + list(patient_samples)
    heatmap_data = heatmap_data.reindex(ordered_samples)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 20))
    
    # Create heatmap
    sns.heatmap(heatmap_data, annot=True, fmt='.1f', cmap='YlOrRd', 
                cbar_kws={'label': 'Percentage (%)'}, ax=ax)
    
    ax.set_title('Cell Type Percentages Across All Samples', fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('Cell Type', fontsize=12)
    ax.set_ylabel('Sample', fontsize=12)
    
    # Add dividing line between organoids and patients
    ax.axhline(y=len(organoid_samples), color='black', linewidth=2)
    
    # Add dataset labels
    ax.text(-0.5, len(organoid_samples)/2, 'Organoids', rotation=90, 
           verticalalignment='center', fontsize=12, fontweight='bold')
    ax.text(-0.5, len(organoid_samples) + len(patient_samples)/2, 'Patient Tumors', 
           rotation=90, verticalalignment='center', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    
    # Save figure
    output_path = '/Users/shaunie/Desktop/hu-glioblastoma/figures/cell_type_heatmap_all_samples.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved figure: {output_path}")
    
    plt.show()
    return fig

def main():
    """Create all visualization figures"""
    
    print("Creating visualization figures for per-sample data...")
    print("="*60)
    
    # Set up plotting
    setup_plotting()
    
    # Create all figures
    fig1 = create_cell_type_distribution_figure()
    print()
    
    fig2 = create_sox2_analysis_figure() 
    print()
    
    fig3 = create_summary_comparison_figure()
    print()
    
    fig4 = create_heatmap_figure()
    print()
    
    print("="*60)
    print("SUCCESS: Created all visualization figures!")
    print("="*60)
    print("\nGenerated figures:")
    print("1. cell_type_distributions_per_sample.png - Cell type distributions and comparisons")
    print("2. sox2_analysis_per_sample.png - Sox2 expression analysis")
    print("3. summary_organoids_vs_patients.png - Overall comparison summary")  
    print("4. cell_type_heatmap_all_samples.png - Heatmap of all samples")
    print(f"\nFigures saved to: /Users/shaunie/Desktop/hu-glioblastoma/figures/")

if __name__ == "__main__":
    main()