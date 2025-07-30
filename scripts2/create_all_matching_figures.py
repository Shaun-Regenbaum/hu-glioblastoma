#!/usr/bin/env python3
"""
Create corresponding figures for each data file in results/Shaun2/
Matching the original Shaun structure with individual single-plot figures
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
    plt.rcParams['figure.figsize'] = (10, 8)
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10

def create_organoid_vs_patient_fourpop_figure():
    """Create OrganoidVsPatient-FourPop.png"""
    
    print("Creating OrganoidVsPatient-FourPop.png...")
    
    # Load aggregated data
    agg_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/OrganoidVsPatient-FourPop.csv'
    df_agg = pd.read_csv(agg_path)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Prepare data for grouped bar chart
    cell_types = df_agg['cell_type_simplified'].values
    organoid_pcts = df_agg['Organoids_Percentage'].values
    patient_pcts = df_agg['GSE131928_Percentage'].values
    
    x = np.arange(len(cell_types))
    width = 0.35
    
    # Create grouped bar chart
    bars1 = ax.bar(x - width/2, organoid_pcts, width, label='Organoids', 
                   alpha=0.8, color='#2ca02c')
    bars2 = ax.bar(x + width/2, patient_pcts, width, label='GSE131928 Patients', 
                   alpha=0.8, color='#ff7f0e')
    
    # Customize plot
    ax.set_title('Cell Type Distribution: Organoids vs Patient Tumors', 
                fontweight='bold', pad=20)
    ax.set_xlabel('Cell Type')
    ax.set_ylabel('Percentage (%)')
    ax.set_xticks(x)
    ax.set_xticklabels(cell_types, rotation=45, ha='right')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                   f'{height:.1f}%', ha='center', va='bottom', fontsize=10)
    
    plt.tight_layout()
    
    output_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/OrganoidVsPatient-FourPop.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

def create_sox2_cell_populations_figure():
    """Create Sox2-Cell Populations.png"""
    
    print("Creating Sox2-Cell Populations.png...")
    
    # Load Sox2 data
    sox2_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/Sox2-Cell Populations.csv'
    df_sox2 = pd.read_csv(sox2_path)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Prepare data for grouped bar chart
    cell_types = df_sox2[df_sox2['Dataset'] == 'Organoid']['Population'].values
    organoid_pcts = df_sox2[df_sox2['Dataset'] == 'Organoid']['Percent_SOX2_Negative_in_Population'].values
    patient_pcts = df_sox2[df_sox2['Dataset'] == 'Patient']['Percent_SOX2_Negative_in_Population'].values
    
    x = np.arange(len(cell_types))
    width = 0.35
    
    # Create grouped bar chart
    bars1 = ax.bar(x - width/2, organoid_pcts, width, label='Organoids', 
                   alpha=0.8, color='#2ca02c')
    bars2 = ax.bar(x + width/2, patient_pcts, width, label='GSE131928 Patients', 
                   alpha=0.8, color='#ff7f0e')
    
    # Customize plot
    ax.set_title('Sox2-Negative Cell Percentages by Population', 
                fontweight='bold', pad=20)
    ax.set_xlabel('Cell Type')
    ax.set_ylabel('% Sox2-Negative')
    ax.set_xticks(x)
    ax.set_xticklabels(cell_types, rotation=45, ha='right')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 1,
                   f'{height:.1f}%', ha='center', va='bottom', fontsize=10)
    
    plt.tight_layout()
    
    output_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/Sox2-Cell Populations.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

def create_foursample_fourpop_figure():
    """Create FourSample-FourPop.png using cell counts"""
    
    print("Creating FourSample-FourPop.png...")
    
    # Load cell counts data
    counts_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/FourSample-FourPop-CellCounts.csv'
    df_counts = pd.read_csv(counts_path)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Convert to percentages for stacked bar chart
    df_pct = df_counts.set_index('sample')
    df_pct = df_pct.div(df_pct.sum(axis=1), axis=0) * 100
    
    # Define colors for cell types
    colors = {
        'Astrocytic': '#2ca02c',
        'Mesenchymal': '#ff7f0e', 
        'Neural_Progenitor': '#1f77b4',
        'Oligodendrocytic': '#d62728'
    }
    
    # Create stacked bar chart
    df_pct.plot(kind='bar', stacked=True, ax=ax, 
                color=[colors[col] for col in df_pct.columns])
    
    ax.set_title('Cell Type Distribution Across Four Organoid Samples', 
                fontweight='bold', pad=20)
    ax.set_xlabel('Organoid Sample')
    ax.set_ylabel('Percentage (%)')
    ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.tick_params(axis='x', rotation=45)
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    output_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/FourSample-FourPop.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

def create_foursample_percentages_figure():
    """Create FourSample-FourPop-Percentages.png showing variability statistics"""
    
    print("Creating FourSample-FourPop-Percentages.png...")
    
    # Load percentages data
    pct_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/FourSample-FourPop-Percentages.csv'
    df_pct = pd.read_csv(pct_path)
    
    # Separate sample data from statistics
    sample_data = df_pct[~df_pct['sample'].isin(['range', 'rMAD'])].copy()
    stats_data = df_pct[df_pct['sample'].isin(['range', 'rMAD'])].copy()
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Left plot: Individual sample percentages
    sample_data_melted = sample_data.melt(id_vars=['sample'], 
                                         var_name='Cell_Type', 
                                         value_name='Percentage')
    
    sns.barplot(data=sample_data_melted, x='Cell_Type', y='Percentage', 
                hue='sample', ax=ax1)
    ax1.set_title('Cell Type Percentages by Organoid Sample', fontweight='bold')
    ax1.set_xlabel('Cell Type')
    ax1.set_ylabel('Percentage (%)')
    ax1.tick_params(axis='x', rotation=45)
    ax1.legend(title='Sample')
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Right plot: Variability statistics
    stats_melted = stats_data.melt(id_vars=['sample'], 
                                  var_name='Cell_Type', 
                                  value_name='Value')
    
    sns.barplot(data=stats_melted, x='Cell_Type', y='Value', 
                hue='sample', ax=ax2)
    ax2.set_title('Cell Type Variability Statistics', fontweight='bold')
    ax2.set_xlabel('Cell Type')
    ax2.set_ylabel('Statistical Value')
    ax2.tick_params(axis='x', rotation=45)
    ax2.legend(title='Statistic')
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    output_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/FourSample-FourPop-Percentages.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

def create_organoid_vs_patient_per_sample_figure():
    """Create OrganoidVsPatient-FourPop-PerSample.png"""
    
    print("Creating OrganoidVsPatient-FourPop-PerSample.png...")
    
    # Load per-sample data
    data_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/OrganoidVsPatient-FourPop-PerSample-REAL.csv'
    df = pd.read_csv(data_path)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create box plot showing variation by dataset and cell type
    plot_data = []
    for _, row in df.iterrows():
        plot_data.append({
            'Dataset': 'Organoid' if row['dataset'] == 'Organoid' else 'Patient',
            'Cell_Type': row['cell_type'],
            'Percentage': row['percentage']
        })
    
    plot_df = pd.DataFrame(plot_data)
    
    # Create box plot
    sns.boxplot(data=plot_df, x='Cell_Type', y='Percentage', hue='Dataset', ax=ax)
    
    ax.set_title('Cell Type Distribution Variation: Per-Sample Analysis', 
                fontweight='bold', pad=20)
    ax.set_xlabel('Cell Type')
    ax.set_ylabel('Percentage (%)')
    ax.tick_params(axis='x', rotation=45)
    ax.legend(title='Dataset')
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    output_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/OrganoidVsPatient-FourPop-PerSample.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

def create_sox2_per_sample_figure():
    """Create Sox2-Cell-Populations-PerSample.png"""
    
    print("Creating Sox2-Cell-Populations-PerSample.png...")
    
    # Load Sox2 per-sample data
    sox2_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/Sox2-Cell-Populations-PerSample-REAL.csv'
    df_sox2 = pd.read_csv(sox2_path)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create box plot for Sox2-negative percentages
    plot_data = []
    for _, row in df_sox2.iterrows():
        plot_data.append({
            'Dataset': 'Organoid' if row['dataset'] == 'Organoid' else 'Patient',
            'Cell_Type': row['cell_type'],
            'Sox2_Negative_Percent': row['percent_SOX2_negative']
        })
    
    plot_df = pd.DataFrame(plot_data)
    
    sns.boxplot(data=plot_df, x='Cell_Type', y='Sox2_Negative_Percent', 
                hue='Dataset', ax=ax)
    
    ax.set_title('Sox2-Negative Percentages: Per-Sample Analysis', 
                fontweight='bold', pad=20)
    ax.set_xlabel('Cell Type')
    ax.set_ylabel('% Sox2-Negative')
    ax.tick_params(axis='x', rotation=45)
    ax.legend(title='Dataset')
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, 100)
    
    plt.tight_layout()
    
    output_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/Sox2-Cell-Populations-PerSample.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

def create_gse131928_per_tumor_figure():
    """Create GSE131928_PerTumor_CellTypes.png"""
    
    print("Creating GSE131928_PerTumor_CellTypes.png...")
    
    # Load per-tumor data
    tumor_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/GSE131928_PerTumor_CellTypes-REAL.csv'
    df_tumor = pd.read_csv(tumor_path)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Left plot: Cell type percentages by age group
    sns.boxplot(data=df_tumor, x='cell_type', y='percentage', 
                hue='age_group', ax=ax1)
    ax1.set_title('Cell Type Percentages by Age Group (GSE131928)', fontweight='bold')
    ax1.set_xlabel('Cell Type')
    ax1.set_ylabel('Percentage (%)')
    ax1.tick_params(axis='x', rotation=45)
    ax1.legend(title='Age Group')
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Right plot: Sample size distribution
    sample_sizes = df_tumor.groupby(['sample', 'age_group'])['total_cells_in_sample'].first().reset_index()
    
    sns.boxplot(data=sample_sizes, x='age_group', y='total_cells_in_sample', ax=ax2)
    ax2.set_title('Sample Size Distribution by Age Group', fontweight='bold')
    ax2.set_xlabel('Age Group')
    ax2.set_ylabel('Total Cells per Sample')
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    output_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/GSE131928_PerTumor_CellTypes.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

def main():
    """Create all matching figures for data files in Shaun2"""
    
    print("Creating corresponding figures for each data file in results/Shaun2/...")
    print("="*80)
    
    # Set up plotting
    setup_plotting()
    
    # Create all figures
    create_organoid_vs_patient_fourpop_figure()
    print()
    
    create_sox2_cell_populations_figure()
    print()
    
    create_foursample_fourpop_figure()
    print()
    
    create_foursample_percentages_figure()
    print()
    
    create_organoid_vs_patient_per_sample_figure()
    print()
    
    create_sox2_per_sample_figure()
    print()
    
    create_gse131928_per_tumor_figure()
    print()
    
    print("="*80)
    print("✅ SUCCESS: Created corresponding figures for all data files!")
    print("="*80)
    print("\nGenerated figure-data pairs in results/Shaun2/:")
    print("1. OrganoidVsPatient-FourPop.csv → OrganoidVsPatient-FourPop.png")
    print("2. Sox2-Cell Populations.csv → Sox2-Cell Populations.png")
    print("3. FourSample-FourPop-CellCounts.csv → FourSample-FourPop.png")
    print("4. FourSample-FourPop-Percentages.csv → FourSample-FourPop-Percentages.png")
    print("5. OrganoidVsPatient-FourPop-PerSample-REAL.csv → OrganoidVsPatient-FourPop-PerSample.png")
    print("6. Sox2-Cell-Populations-PerSample-REAL.csv → Sox2-Cell-Populations-PerSample.png")
    print("7. GSE131928_PerTumor_CellTypes-REAL.csv → GSE131928_PerTumor_CellTypes.png")
    print(f"\n🔬 All figures use data from 31 GSE131928 tumors + 4 organoid samples")

if __name__ == "__main__":
    main()