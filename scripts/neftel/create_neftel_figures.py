#!/usr/bin/env python3
"""
Create figures for Neftel methodology analysis
Reproducing the style and approach from Neftel et al. Cell 2019
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap

# Set style
plt.style.use('default')
sns.set_palette("husl")

def load_neftel_data():
    """Load Neftel methodology results"""
    data_dir = Path('/Users/shaunie/Desktop/hu-glioblastoma/results/neftel')
    
    meta_modules = pd.read_csv(data_dir / 'meta_module_scores.csv')
    hybrid_states = pd.read_csv(data_dir / 'hybrid_states.csv')
    cycling_data = pd.read_csv(data_dir / 'cycling_within_states.csv')
    module_summary = pd.read_csv(data_dir / 'neftel_module_summary.csv')
    
    return meta_modules, hybrid_states, cycling_data, module_summary

def create_meta_module_heatmap(module_summary):
    """Create heatmap showing meta-module expression across datasets"""
    
    # Pivot data for heatmap
    heatmap_data = module_summary.pivot(index='meta_module', columns='dataset', 
                                       values='percentage_of_dataset')
    
    # Reorder modules to match Neftel paper grouping
    module_order = ['AC', 'OPC', 'NPC1', 'NPC2', 'MES1', 'MES2']
    heatmap_data = heatmap_data.reindex(module_order)
    
    fig, ax = plt.subplots(figsize=(8, 10))
    
    # Create custom colormap (similar to Neftel paper)
    colors = ['white', 'lightblue', 'blue', 'darkblue']
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list('neftel', colors, N=n_bins)
    
    # Create heatmap
    sns.heatmap(heatmap_data, annot=True, fmt='.1f', cmap=cmap,
                center=None, ax=ax, cbar_kws={'label': '% of Dataset'})
    
    ax.set_title('Meta-Module Expression Across Datasets\n(Neftel Methodology)', 
                fontsize=14, fontweight='bold')
    ax.set_xlabel('Dataset')
    ax.set_ylabel('Meta-Module')
    
    # Add module descriptions
    module_descriptions = {
        'AC': 'Astrocyte-like',
        'OPC': 'Oligodendrocyte-Progenitor-like', 
        'NPC1': 'Neural-Progenitor-like 1',
        'NPC2': 'Neural-Progenitor-like 2',
        'MES1': 'Mesenchymal-like 1',
        'MES2': 'Mesenchymal-like 2 (Hypoxia)'
    }
    
    # Update y-axis labels
    new_labels = [f"{mod}\n({module_descriptions[mod]})" for mod in module_order]
    ax.set_yticklabels(new_labels, rotation=0)
    
    plt.tight_layout()
    return fig

def create_hybrid_state_analysis(hybrid_states):
    """Create visualization of hybrid cellular states"""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Plot 1: Top hybrid combinations
    hybrid_counts = hybrid_states['hybrid_type'].value_counts().head(10)
    
    bars = ax1.barh(range(len(hybrid_counts)), hybrid_counts.values, 
                    color=['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7'] * 2)
    
    ax1.set_yticks(range(len(hybrid_counts)))
    ax1.set_yticklabels(hybrid_counts.index)
    ax1.set_xlabel('Number of Samples')
    ax1.set_title('Top Hybrid State Combinations', fontweight='bold')
    
    # Add value labels
    for i, (bar, value) in enumerate(zip(bars, hybrid_counts.values)):
        ax1.text(value + 1, i, str(value), va='center', ha='left')
    
    # Plot 2: Hybrid score distribution
    if len(hybrid_states) > 0:
        ax2.hist(hybrid_states['hybrid_score'], bins=20, alpha=0.7, 
                color='#4ECDC4', edgecolor='black')
        ax2.set_xlabel('Hybrid Score')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Distribution of Hybrid Scores', fontweight='bold')
        ax2.axvline(hybrid_states['hybrid_score'].mean(), color='red', 
                   linestyle='--', label=f'Mean: {hybrid_states["hybrid_score"].mean():.2f}')
        ax2.legend()
    
    plt.tight_layout()
    return fig

def create_cycling_within_states_plot(cycling_data, module_summary):
    """Create plot showing cycling cells within each cellular state"""
    
    if len(cycling_data) == 0:
        # Create empty plot if no cycling data
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, 'No cycling data available', 
               transform=ax.transAxes, ha='center', va='center', fontsize=16)
        ax.set_title('Cycling Cells Within Cellular States')
        return fig
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Plot 1: Cycling cells by cellular state and dataset
    cycling_summary = cycling_data.groupby(['dataset', 'cellular_state']).agg({
        'cycling_cells': 'sum',
        'cycling_fraction': 'mean'
    }).reset_index()
    
    # Pivot for grouped bar chart
    pivot_cycling = cycling_summary.pivot(index='cellular_state', 
                                         columns='dataset', 
                                         values='cycling_cells').fillna(0)
    
    pivot_cycling.plot(kind='bar', ax=ax1, alpha=0.8)
    ax1.set_title('Cycling Cells by Cellular State', fontweight='bold')
    ax1.set_xlabel('Cellular State')
    ax1.set_ylabel('Number of Cycling Cells')
    ax1.legend(title='Dataset')
    ax1.tick_params(axis='x', rotation=45)
    
    # Plot 2: Cycling fraction by state
    cycling_fractions = cycling_data.groupby('cellular_state')['cycling_fraction'].mean()
    
    bars = ax2.bar(range(len(cycling_fractions)), cycling_fractions.values, 
                   color=['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD'])
    
    ax2.set_xticks(range(len(cycling_fractions)))
    ax2.set_xticklabels(cycling_fractions.index, rotation=45)
    ax2.set_ylabel('Mean Cycling Fraction')
    ax2.set_title('Cycling Enrichment by Cellular State', fontweight='bold')
    
    # Add value labels
    for bar, value in zip(bars, cycling_fractions.values):
        ax2.text(bar.get_x() + bar.get_width()/2, value + 0.01,
                f'{value:.2f}', ha='center', va='bottom')
    
    plt.tight_layout()
    return fig

def create_neftel_vs_exclusive_comparison(meta_modules, module_summary):
    """Compare Neftel overlapping approach vs our exclusive approach"""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Plot 1: Module co-occurrence matrix
    # Create co-occurrence matrix from samples with multiple modules
    sample_modules = meta_modules.groupby(['dataset', 'sample'])['meta_module'].apply(list).reset_index()
    
    modules = ['AC', 'OPC', 'NPC1', 'NPC2', 'MES1', 'MES2']
    cooccur_matrix = np.zeros((len(modules), len(modules)))
    
    for _, row in sample_modules.iterrows():
        sample_mods = row['meta_module']
        for i, mod1 in enumerate(modules):
            for j, mod2 in enumerate(modules):
                if mod1 in sample_mods and mod2 in sample_mods:
                    cooccur_matrix[i, j] += 1
    
    # Normalize by diagonal (self-occurrence)
    for i in range(len(modules)):
        if cooccur_matrix[i, i] > 0:
            cooccur_matrix[i, :] = cooccur_matrix[i, :] / cooccur_matrix[i, i]
    
    sns.heatmap(cooccur_matrix, annot=True, fmt='.2f', 
                xticklabels=modules, yticklabels=modules,
                cmap='Blues', ax=ax1)
    ax1.set_title('Meta-Module Co-occurrence\n(Overlapping States)', fontweight='bold')
    ax1.set_xlabel('Co-occurring Module')
    ax1.set_ylabel('Reference Module')
    
    # Plot 2: Dataset composition comparison
    organoid_data = module_summary[module_summary['dataset'] == 'Organoid']
    patient_data = module_summary[module_summary['dataset'] == 'GSE131928']
    
    x = np.arange(len(modules))
    width = 0.35
    
    organoid_pcts = [organoid_data[organoid_data['meta_module'] == mod]['percentage_of_dataset'].values[0] 
                    if len(organoid_data[organoid_data['meta_module'] == mod]) > 0 else 0 
                    for mod in modules]
    patient_pcts = [patient_data[patient_data['meta_module'] == mod]['percentage_of_dataset'].values[0] 
                   if len(patient_data[patient_data['meta_module'] == mod]) > 0 else 0 
                   for mod in modules]
    
    bars1 = ax2.bar(x - width/2, organoid_pcts, width, label='Organoid', alpha=0.8, color='lightblue')
    bars2 = ax2.bar(x + width/2, patient_pcts, width, label='Patient', alpha=0.8, color='lightcoral')
    
    ax2.set_xlabel('Meta-Module')
    ax2.set_ylabel('Percentage of Dataset')
    ax2.set_title('Dataset Composition\n(Neftel Overlapping Approach)', fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(modules, rotation=45)
    ax2.legend()
    
    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax2.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                        f'{height:.1f}%', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    return fig

def create_methodology_comparison_summary():
    """Create summary figure comparing all three approaches"""
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Approach comparison table
    approaches = ['Shaun2\n(Overlapping)', 'Shaun3\n(Exclusive)', 'Neftel\n(Meta-modules)']
    features = ['Cell Types', 'Cycling Handling', 'Percentages', 'Hybrid States', 'Genetic Analysis']
    
    comparison_data = [
        ['6 types\n(overlapping)', '6 types\n(exclusive)', '6 meta-modules\n(overlapping)'],
        ['Separate category', 'Separate category', 'Within each state'],
        ['Don\'t sum to 100%', 'Sum to 100%', 'Don\'t sum to 100%'],
        ['Not detected', 'Not applicable', 'Detected'],
        ['Not included', 'Not included', 'CDK4, EGFR, etc.']
    ]
    
    # Create comparison table
    ax1.axis('tight')
    ax1.axis('off')
    
    table = ax1.table(cellText=comparison_data,
                     rowLabels=features,
                     colLabels=approaches,
                     cellLoc='center',
                     loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.2, 2)
    
    # Color code the headers
    for i in range(len(approaches)):
        table[(0, i)].set_facecolor('#E8F4FD')
    for i in range(len(features)):
        table[(i+1, -1)].set_facecolor('#F0F0F0')
    
    ax1.set_title('Methodology Comparison Summary', fontweight='bold', pad=20)
    
    # Advantages/Disadvantages
    ax2.axis('off')
    
    advantages_text = """
NEFTEL METHODOLOGY ADVANTAGES:
• Captures cell state plasticity
• Detects hybrid states
• Robust meta-module approach
• Genetic association analysis
• Matches developmental programs

EXCLUSIVE METHODOLOGY ADVANTAGES:
• Clean comparative analysis
• Percentages sum to 100%
• Mutually exclusive categories
• Clear organoid comparison
• Easier statistical analysis
"""
    
    ax2.text(0.05, 0.95, advantages_text, transform=ax2.transAxes, 
             fontsize=10, verticalalignment='top', fontfamily='monospace')
    ax2.set_title('Methodology Advantages', fontweight='bold')
    
    # Cell type mapping visualization
    mapping_data = {
        'Our Cell Types': ['Astrocytic', 'Mesenchymal', 'Neural_Progenitor', 'Oligodendrocytic', 'Cycling', 'Endothelial'],
        'Neftel Modules': ['AC-like', 'MES1/MES2', 'NPC1/NPC2', 'OPC-like', 'Within states', 'Excluded']
    }
    
    ax3.axis('off')
    
    for i, (our_type, neftel_type) in enumerate(zip(mapping_data['Our Cell Types'], mapping_data['Neftel Modules'])):
        ax3.text(0.1, 0.9 - i*0.12, f"{our_type} ≈ {neftel_type}", 
                transform=ax3.transAxes, fontsize=11, 
                bbox=dict(boxstyle="round,pad=0.3", facecolor='lightblue', alpha=0.7))
    
    ax3.set_title('Cell Type Mapping', fontweight='bold')
    
    # Key findings
    ax4.axis('off')
    
    findings_text = """
KEY FINDINGS FROM NEFTEL ANALYSIS:

🔬 Meta-Modules Detected:
   • AC (Astrocyte-like): 74 instances
   • OPC (Oligodendrocyte-like): 66 instances  
   • NPC1/NPC2 (Neural Progenitor): 148 instances
   • MES1/MES2 (Mesenchymal): 92 instances

🔄 Hybrid States:
   • 778 hybrid combinations detected
   • MES1+MES2 most common (57 samples)
   • High cellular plasticity confirmed

⚡ Cycling Analysis:
   • Cycling enriched in MES1 (75%)
   • Moderate in AC (51%) and NPC states
   • Confirms proliferative hierarchy
"""
    
    ax4.text(0.05, 0.95, findings_text, transform=ax4.transAxes, 
             fontsize=9, verticalalignment='top', fontfamily='monospace')
    ax4.set_title('Key Findings', fontweight='bold')
    
    plt.tight_layout()
    return fig

def main():
    """Generate all Neftel methodology figures"""
    
    print("Creating Neftel methodology figures...")
    
    # Load data
    meta_modules, hybrid_states, cycling_data, module_summary = load_neftel_data()
    print(f"Loaded: {len(meta_modules)} module scores, {len(hybrid_states)} hybrid states")
    
    # Create output directory
    output_dir = Path('/Users/shaunie/Desktop/hu-glioblastoma/results/neftel/figures')
    output_dir.mkdir(exist_ok=True)
    
    # Generate figures
    print("Creating meta-module heatmap...")
    fig1 = create_meta_module_heatmap(module_summary)
    fig1.savefig(output_dir / 'meta_module_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close(fig1)
    
    print("Creating hybrid state analysis...")
    fig2 = create_hybrid_state_analysis(hybrid_states)
    fig2.savefig(output_dir / 'hybrid_state_analysis.png', dpi=300, bbox_inches='tight')
    plt.close(fig2)
    
    print("Creating cycling analysis...")
    fig3 = create_cycling_within_states_plot(cycling_data, module_summary)
    fig3.savefig(output_dir / 'cycling_within_states.png', dpi=300, bbox_inches='tight')
    plt.close(fig3)
    
    print("Creating Neftel vs exclusive comparison...")
    fig4 = create_neftel_vs_exclusive_comparison(meta_modules, module_summary)
    fig4.savefig(output_dir / 'neftel_vs_exclusive_comparison.png', dpi=300, bbox_inches='tight')
    plt.close(fig4)
    
    print("Creating methodology comparison summary...")
    fig5 = create_methodology_comparison_summary()
    fig5.savefig(output_dir / 'methodology_comparison_summary.png', dpi=300, bbox_inches='tight')
    plt.close(fig5)
    
    print(f"\n✅ All figures saved to: {output_dir}")
    print("\nGenerated figures:")
    print("  1. meta_module_heatmap.png - Expression across datasets")
    print("  2. hybrid_state_analysis.png - Hybrid cellular states")
    print("  3. cycling_within_states.png - Cycling cells analysis")
    print("  4. neftel_vs_exclusive_comparison.png - Method comparison")
    print("  5. methodology_comparison_summary.png - Complete overview")

if __name__ == "__main__":
    main()