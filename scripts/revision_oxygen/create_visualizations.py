#!/usr/bin/env python3
"""
Create visualizations for drug response data
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

class DrugResponseVisualizer:
    def __init__(self, base_path):
        self.base_path = Path(base_path)
        self.results_dir = self.base_path / "results"
        self.figures_dir = self.base_path / "figures"
        self.figures_dir.mkdir(exist_ok=True)
        
    def load_time_series_data(self):
        """Load the time series data"""
        file_path = self.results_dir / "all_drugs_time_series.csv"
        return pd.read_csv(file_path)
        
    def plot_dose_response_curves(self):
        """Create dose response curves for each drug"""
        data = self.load_time_series_data()
        
        # Create 2x2 subplot
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        axes = axes.flatten()
        
        drugs = ['DACTI', 'PLICA', 'GEMCI', 'VINC']
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
        
        for idx, drug in enumerate(drugs):
            ax = axes[idx]
            
            # Find all columns for this drug
            drug_cols = [col for col in data.columns if col.startswith(f"{drug}_")]
            
            if not drug_cols:
                continue
                
            # Group by concentration
            conc_groups = defaultdict(list)
            for col in drug_cols:
                parts = col.split('_')
                if len(parts) >= 3:
                    conc = float(parts[1])
                    conc_groups[conc].append(col)
            
            # Plot each concentration
            for i, (conc, cols) in enumerate(sorted(conc_groups.items())):
                # Calculate mean and std across replicates
                replicate_data = data[cols].values
                
                # Remove NaN values
                time_points = []
                means = []
                stds = []
                
                for t_idx in range(len(data)):
                    values = [val for val in replicate_data[t_idx] if not pd.isna(val)]
                    if values:
                        time_points.append(data.iloc[t_idx]['Time_h'])
                        means.append(np.mean(values))
                        stds.append(np.std(values))
                
                if time_points:
                    means = np.array(means)
                    stds = np.array(stds)
                    
                    # Use color gradients for concentrations
                    color_intensity = i / max(1, len(conc_groups) - 1)
                    color = plt.cm.viridis(color_intensity)
                    
                    ax.plot(time_points, means, color=color, linewidth=2, 
                           label=f'{conc:.3g} μM', alpha=0.8)
                    ax.fill_between(time_points, means - stds, means + stds, 
                                  color=color, alpha=0.2)
            
            ax.set_xlabel('Time (hours)', fontsize=12)
            ax.set_ylabel('O₂ (% air saturation)', fontsize=12)
            ax.set_title(f'{drug} Dose Response', fontsize=14, fontweight='bold')
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
            ax.grid(True, alpha=0.3)
            ax.set_xlim(0, 96)
            
        plt.suptitle('Drug Dose Response Curves', fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(self.figures_dir / 'dose_response_curves.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Created dose response curves")
        
    def plot_drug_comparison(self):
        """Create comparison plot with representative concentrations"""
        data = self.load_time_series_data()
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Representative concentrations for each drug (mid-range)
        drug_concs = {
            'DACTI': 0.02,     # Mid concentration
            'PLICA': 0.03,     # Mid concentration
            'GEMCI': 5.0,      # Mid concentration
            'VINC': 1.5        # Mid concentration
        }
        
        colors = {
            'DACTI': '#1f77b4',
            'PLICA': '#ff7f0e', 
            'GEMCI': '#2ca02c',
            'VINC': '#d62728',
            'DMSO': '#9467bd'
        }
        
        # Plot each drug at representative concentration
        for drug, target_conc in drug_concs.items():
            # Find columns matching this drug and concentration
            drug_cols = []
            for col in data.columns:
                if col.startswith(f"{drug}_"):
                    parts = col.split('_')
                    if len(parts) >= 3:
                        conc = float(parts[1])
                        if abs(conc - target_conc) < 0.001:  # Close match
                            drug_cols.append(col)
            
            if drug_cols:
                # Calculate mean across replicates
                replicate_data = data[drug_cols].values
                
                time_points = []
                means = []
                stds = []
                
                for t_idx in range(len(data)):
                    values = [val for val in replicate_data[t_idx] if not pd.isna(val)]
                    if values:
                        time_points.append(data.iloc[t_idx]['Time_h'])
                        means.append(np.mean(values))
                        stds.append(np.std(values))
                
                if time_points:
                    means = np.array(means)
                    stds = np.array(stds)
                    
                    ax.plot(time_points, means, color=colors[drug], linewidth=2.5,
                           label=f'{drug} ({target_conc} μM)')
                    ax.fill_between(time_points, means - stds, means + stds,
                                  color=colors[drug], alpha=0.2)
        
        # Add DMSO control (average of all DMSO)
        dmso_cols = [col for col in data.columns if col.startswith('DMSO_')]
        if dmso_cols:
            replicate_data = data[dmso_cols].values
            
            time_points = []
            means = []
            
            for t_idx in range(len(data)):
                values = [val for val in replicate_data[t_idx] if not pd.isna(val)]
                if values:
                    time_points.append(data.iloc[t_idx]['Time_h'])
                    means.append(np.mean(values))
            
            if time_points:
                ax.plot(time_points, means, '--', color=colors['DMSO'],
                       linewidth=2, label='DMSO Control')
        
        ax.set_xlabel('Time (hours)', fontsize=14)
        ax.set_ylabel('O₂ (% air saturation)', fontsize=14)
        ax.set_title('Drug Response Comparison', fontsize=16, fontweight='bold')
        ax.legend(fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 96)
        
        plt.tight_layout()
        plt.savefig(self.figures_dir / 'drug_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Created drug comparison plot")
        
    def plot_replicate_variability(self):
        """Plot showing variability across replicates"""
        data = self.load_time_series_data()
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        axes = axes.flatten()
        
        drugs = ['DACTI', 'PLICA', 'GEMCI', 'VINC']
        
        for idx, drug in enumerate(drugs):
            ax = axes[idx]
            
            # Find all columns for this drug
            drug_cols = [col for col in data.columns if col.startswith(f"{drug}_")]
            
            # Group by concentration and plot individual replicates
            conc_groups = defaultdict(list)
            for col in drug_cols:
                parts = col.split('_')
                if len(parts) >= 3:
                    conc = float(parts[1])
                    conc_groups[conc].append(col)
            
            # Pick one concentration to show replicate variability
            if conc_groups:
                # Use the concentration with most replicates
                best_conc = max(conc_groups.keys(), key=lambda x: len(conc_groups[x]))
                cols = conc_groups[best_conc]
                
                for i, col in enumerate(cols):
                    values = data[col].dropna()
                    times = data.loc[values.index, 'Time_h']
                    
                    ax.plot(times, values, alpha=0.7, linewidth=1.5,
                           label=f'Replicate {i+1}')
                
                # Add mean line
                replicate_data = data[cols].values
                time_points = []
                means = []
                
                for t_idx in range(len(data)):
                    values = [val for val in replicate_data[t_idx] if not pd.isna(val)]
                    if values:
                        time_points.append(data.iloc[t_idx]['Time_h'])
                        means.append(np.mean(values))
                
                if time_points:
                    ax.plot(time_points, means, 'k-', linewidth=3, alpha=0.8,
                           label='Mean')
                
                ax.set_xlabel('Time (hours)', fontsize=12)
                ax.set_ylabel('O₂ (% air saturation)', fontsize=12)
                ax.set_title(f'{drug} at {best_conc:.3g} μM - Replicates', fontsize=14)
                ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
                ax.grid(True, alpha=0.3)
                ax.set_xlim(0, 96)
        
        plt.suptitle('Replicate Variability', fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(self.figures_dir / 'replicate_variability.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Created replicate variability plot")
        
    def plot_endpoint_heatmap(self):
        """Create heatmap of endpoint values (8-hour timepoint)"""
        data = self.load_time_series_data()
        
        # Find the 8-hour timepoint
        target_time = 8.0
        time_diff = (data['Time_h'] - target_time).abs()
        endpoint_time_idx = time_diff.idxmin()
        actual_time = data.iloc[endpoint_time_idx]['Time_h']
        print(f"Using time point {actual_time}h for endpoint heatmap (target: {target_time}h)")
        
        # Organize data by drug and concentration
        drug_data = {}
        drugs = ['DACTI', 'PLICA', 'GEMCI', 'VINC']
        
        for drug in drugs:
            drug_cols = [col for col in data.columns if col.startswith(f"{drug}_")]
            
            conc_groups = defaultdict(list)
            for col in drug_cols:
                parts = col.split('_')
                if len(parts) >= 3:
                    conc = float(parts[1])
                    
                    # Get the value at the 8-hour timepoint
                    endpoint_value = data.iloc[endpoint_time_idx][col]
                    if not pd.isna(endpoint_value):
                        conc_groups[conc].append(endpoint_value)
            
            # Calculate mean endpoint for each concentration
            drug_data[drug] = {}
            for conc, values in conc_groups.items():
                if values:
                    drug_data[drug][conc] = np.mean(values)
        
        # Create heatmap
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Prepare data matrix
        all_concs = set()
        for drug_concs in drug_data.values():
            all_concs.update(drug_concs.keys())
        all_concs = sorted(all_concs)
        
        matrix = []
        y_labels = []
        
        for drug in drugs:
            if drug in drug_data:
                row = []
                for conc in all_concs:
                    if conc in drug_data[drug]:
                        row.append(drug_data[drug][conc])
                    else:
                        row.append(np.nan)
                matrix.append(row)
                y_labels.append(drug)
        
        matrix = np.array(matrix)
        
        # Create heatmap
        im = ax.imshow(matrix, aspect='auto', cmap='viridis_r', interpolation='nearest')
        
        # Set labels
        ax.set_xticks(range(len(all_concs)))
        ax.set_xticklabels([f'{c:.3g}' for c in all_concs], rotation=45, ha='right')
        ax.set_xlabel('Concentration (μM)', fontsize=12)
        
        ax.set_yticks(range(len(y_labels)))
        ax.set_yticklabels(y_labels)
        ax.set_ylabel('Drug', fontsize=12)
        
        ax.set_title('Endpoint O₂ Values (8-Hour Timepoint)', fontsize=16, fontweight='bold')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('O₂ (% air saturation)', fontsize=12)
        
        # Add text annotations
        for i in range(len(y_labels)):
            for j in range(len(all_concs)):
                if not np.isnan(matrix[i, j]):
                    text = ax.text(j, i, f'{matrix[i, j]:.1f}',
                                 ha="center", va="center", color="white", fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(self.figures_dir / 'endpoint_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Created endpoint heatmap")
        
    def plot_drug_vs_dmso(self, drug_name, drug_display_name):
        """Create individual plot comparing one drug's concentrations with DMSO"""
        data = self.load_time_series_data()
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Find all columns for this drug
        drug_cols = [col for col in data.columns if col.startswith(f"{drug_name.upper()}_")]
        
        if not drug_cols:
            print(f"No data found for drug: {drug_name}")
            return
            
        # Group by concentration
        conc_groups = defaultdict(list)
        for col in drug_cols:
            parts = col.split('_')
            if len(parts) >= 3:
                conc = float(parts[1])
                conc_groups[conc].append(col)
        
        # Plot each concentration
        colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(conc_groups)))
        
        for i, (conc, cols) in enumerate(sorted(conc_groups.items())):
            # Calculate mean and std across replicates
            replicate_data = data[cols].values
            
            # Remove NaN values
            time_points = []
            means = []
            stds = []
            
            for t_idx in range(len(data)):
                values = [val for val in replicate_data[t_idx] if not pd.isna(val)]
                if values:
                    time_points.append(data.iloc[t_idx]['Time_h'])
                    means.append(np.mean(values))
                    stds.append(np.std(values))
            
            if time_points:
                means = np.array(means)
                stds = np.array(stds)
                
                ax.plot(time_points, means, color=colors[i], linewidth=2.5, 
                       label=f'{conc:.3g} μM', alpha=0.8)
                ax.fill_between(time_points, means - stds, means + stds, 
                              color=colors[i], alpha=0.2)
        
        # Add DMSO control (average of all DMSO)
        dmso_cols = [col for col in data.columns if col.startswith('DMSO_')]
        if dmso_cols:
            replicate_data = data[dmso_cols].values
            
            time_points = []
            means = []
            stds = []
            
            for t_idx in range(len(data)):
                values = [val for val in replicate_data[t_idx] if not pd.isna(val)]
                if values:
                    time_points.append(data.iloc[t_idx]['Time_h'])
                    means.append(np.mean(values))
                    stds.append(np.std(values))
            
            if time_points:
                means = np.array(means)
                stds = np.array(stds)
                
                ax.plot(time_points, means, '--', color='black',
                       linewidth=3, label='DMSO Control', alpha=0.7)
                ax.fill_between(time_points, means - stds, means + stds,
                              color='black', alpha=0.1)
        
        ax.set_xlabel('Time (hours)', fontsize=14)
        ax.set_ylabel('O₂ (% air saturation)', fontsize=14)
        ax.set_title(f'{drug_display_name} vs DMSO Control', fontsize=16, fontweight='bold')
        ax.legend(fontsize=12, bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 96)
        
        plt.tight_layout()
        plt.savefig(self.figures_dir / f'{drug_name}_vs_dmso.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Created {drug_display_name} vs DMSO plot")
    
    def export_drug_vs_dmso_csv(self, drug_name, drug_display_name):
        """Export CSV data for individual drug vs DMSO comparison"""
        data = self.load_time_series_data()
        
        # Find all columns for this drug
        drug_cols = [col for col in data.columns if col.startswith(f"{drug_name.upper()}_")]
        dmso_cols = [col for col in data.columns if col.startswith('DMSO_')]
        
        if not drug_cols:
            print(f"No data found for drug: {drug_name}")
            return
            
        # Prepare export data
        export_data = {'Time_h': data['Time_h'].values}
        
        # Add drug concentration data
        conc_groups = defaultdict(list)
        for col in drug_cols:
            parts = col.split('_')
            if len(parts) >= 3:
                conc = float(parts[1])
                conc_groups[conc].append(col)
        
        # Add individual replicate columns and means for each concentration
        for conc in sorted(conc_groups.keys()):
            cols = conc_groups[conc]
            
            # Add individual replicate columns
            for i, col in enumerate(cols):
                export_data[f'{drug_display_name}_{conc}uM_rep{i+1}'] = data[col].values
            
            # Calculate and add mean
            replicate_values = data[cols].values
            means = []
            for t_idx in range(len(data)):
                values = [val for val in replicate_values[t_idx] if not pd.isna(val)]
                if values:
                    means.append(np.mean(values))
                else:
                    means.append(np.nan)
            export_data[f'{drug_display_name}_{conc}uM_mean'] = means
        
        # Add DMSO data
        if dmso_cols:
            # Add individual DMSO replicates
            for i, col in enumerate(dmso_cols):
                export_data[f'DMSO_rep{i+1}'] = data[col].values
            
            # Add DMSO mean
            replicate_values = data[dmso_cols].values
            means = []
            for t_idx in range(len(data)):
                values = [val for val in replicate_values[t_idx] if not pd.isna(val)]
                if values:
                    means.append(np.mean(values))
                else:
                    means.append(np.nan)
            export_data['DMSO_mean'] = means
        
        # Create DataFrame and save
        export_df = pd.DataFrame(export_data)
        csv_path = self.results_dir / f'{drug_name}_vs_dmso_data.csv'
        export_df.to_csv(csv_path, index=False)
        
        print(f"Exported {drug_display_name} vs DMSO data to: {csv_path}")
    
    def export_all_drug_vs_dmso_csvs(self):
        """Export CSV data for all individual drug vs DMSO comparisons"""
        drugs = [
            ('dacti', 'DACT'),
            ('plica', 'PLIC'), 
            ('gemci', 'GEMCI'),
            ('vinc', 'VINC')
        ]
        
        for drug_name, drug_display_name in drugs:
            try:
                self.export_drug_vs_dmso_csv(drug_name, drug_display_name)
            except Exception as e:
                print(f"Error exporting {drug_display_name} data: {e}")
    
    def plot_individual_drug_comparisons(self):
        """Create individual plots for all drugs vs DMSO"""
        self.plot_drug_vs_dmso('dacti', 'DACT')
        self.plot_drug_vs_dmso('plica', 'PLIC')
        self.plot_drug_vs_dmso('gemci', 'GEMCI')
        self.plot_drug_vs_dmso('vinc', 'VINC')
        
    def create_all_plots(self):
        """Create all visualization plots"""
        print("Creating visualizations...")
        
        try:
            self.plot_dose_response_curves()
        except Exception as e:
            print(f"Error creating dose response curves: {e}")
            
        try:
            self.plot_drug_comparison()
        except Exception as e:
            print(f"Error creating drug comparison: {e}")
            
        try:
            self.plot_replicate_variability()
        except Exception as e:
            print(f"Error creating replicate variability: {e}")
            
        try:
            self.plot_endpoint_heatmap()
        except Exception as e:
            print(f"Error creating endpoint heatmap: {e}")
            
        try:
            self.plot_individual_drug_comparisons()
        except Exception as e:
            print(f"Error creating individual drug comparisons: {e}")
            
        try:
            self.export_all_drug_vs_dmso_csvs()
        except Exception as e:
            print(f"Error exporting drug vs DMSO CSV data: {e}")
        
        print(f"All plots saved to: {self.figures_dir}")
        print(f"CSV data saved to: {self.results_dir}")


def main():
    visualizer = DrugResponseVisualizer("/Users/shaunie/Desktop/hu-glioblastoma/data/revision")
    visualizer.create_all_plots()


if __name__ == "__main__":
    main()