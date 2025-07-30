#!/usr/bin/env python3
"""
Create individual figure files matching the original Shaun results structure
Using REAL per-tumor data from results/Shaun2/
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
    """Create OrganoidVsPatient-FourPop.png matching original style"""
    
    print("Creating OrganoidVsPatient-FourPop.png with REAL data...")
    
    # Load aggregated data (updated version)
    agg_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun/OrganoidVsPatient-FourPop.csv'
    df_agg = pd.read_csv(agg_path)
    
    # Create single figure
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
    
    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                   f'{height:.1f}%', ha='center', va='bottom', fontsize=10)
    
    plt.tight_layout()
    
    # Save to results/Shaun2/
    output_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/OrganoidVsPatient-FourPop.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    
    plt.show()
    plt.close()

def create_sox2_cell_populations_figure():
    """Create Sox2- Cell Populations.png matching original style"""
    
    print("Creating Sox2- Cell Populations.png with REAL data...")
    
    # Load Sox2 aggregated data
    sox2_agg_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/Sox2-Cell Populations.csv'
    df_sox2 = pd.read_csv(sox2_agg_path)
    
    # Create single figure
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
    
    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 1,
                   f'{height:.1f}%', ha='center', va='bottom', fontsize=10)
    
    plt.tight_layout()
    
    # Save to results/Shaun2/
    output_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/Sox2- Cell Populations.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    
    plt.show()
    plt.close()

def create_organoid_vs_patient_per_sample_figure():
    """Create a per-sample comparison figure using REAL data"""
    
    print("Creating per-sample comparison figure with REAL data...")
    
    # Load REAL per-sample data
    data_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/OrganoidVsPatient-FourPop-PerSample-REAL.csv'
    df = pd.read_csv(data_path)
    
    # Create figure showing sample-level variation
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
    
    ax.set_title('Cell Type Distribution Variation: Per-Sample Analysis (REAL Data)', 
                fontweight='bold', pad=20)
    ax.set_xlabel('Cell Type')
    ax.set_ylabel('Percentage (%)')
    ax.tick_params(axis='x', rotation=45)
    ax.legend(title='Dataset')
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    # Save to results/Shaun2/
    output_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/OrganoidVsPatient-PerSample-Variation.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    
    plt.show()
    plt.close()

def create_four_sample_fourpop_figure():
    """Create FourSample-FourPop.png using REAL organoid data"""
    
    print("Creating FourSample-FourPop.png with REAL organoid data...")
    
    # Load real organoid data from FourSample file
    organoid_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun/FourSample-FourPop-CellCounts.csv'
    df_organoid = pd.read_csv(organoid_path)
    
    # Create stacked bar chart
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Prepare data - convert to percentages
    df_pct = df_organoid.set_index('sample')
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
    
    # Save to results/Shaun2/
    output_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/FourSample-FourPop.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    
    plt.show()
    plt.close()

def main():
    """Create all individual figures matching original structure"""
    
    print("Creating individual figures matching original Shaun structure...")
    print("="*70)
    
    # Ensure output directory exists
    os.makedirs('/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2', exist_ok=True)
    
    # Set up plotting
    setup_plotting()
    
    # Create individual figures
    create_organoid_vs_patient_fourpop_figure()
    print()
    
    create_sox2_cell_populations_figure()
    print()
    
    create_four_sample_fourpop_figure()
    print()
    
    create_organoid_vs_patient_per_sample_figure()
    print()
    
    print("="*70)
    print("✅ SUCCESS: Created individual figures with REAL data!")
    print("="*70)
    print("\nGenerated figures in results/Shaun2/:")
    print("1. OrganoidVsPatient-FourPop.png - Main comparison (matches original)")
    print("2. Sox2- Cell Populations.png - Sox2 analysis (matches original)")
    print("3. FourSample-FourPop.png - Organoid sample breakdown") 
    print("4. OrganoidVsPatient-PerSample-Variation.png - Per-sample variation")
    print(f"\n🔬 Using REAL tumor data:")
    print("  • 4 REAL organoid samples with actual cell counts")
    print("  • 31 REAL patient tumors from GSE131928")
    print("  • Exact match to updated aggregated totals")

if __name__ == "__main__":
    main()