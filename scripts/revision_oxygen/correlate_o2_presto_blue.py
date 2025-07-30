#!/usr/bin/env python3
"""
Correlate O2 values with Presto Blue results
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

class O2PrestoBlueCorrelation:
    def __init__(self, base_path):
        self.base_path = Path(base_path)
        self.results_dir = self.base_path / "results"
        self.figures_dir = self.base_path / "figures"
        
    def load_o2_data(self, use_averaged=True):
        """Load O2 data from the final plate layout CSV"""
        if use_averaged:
            filename = "final_plate_layout_with_values_avg_last_10h.csv"
        else:
            filename = "final_plate_layout_with_values_latest_94.0h.csv"
        
        file_path = self.results_dir / filename
        return pd.read_csv(file_path)
    
    def load_presto_blue_data(self):
        """Load Presto Blue results"""
        file_path = self.base_path / "presto blue results with map.csv"
        df = pd.read_csv(file_path, index_col=0)
        
        # Convert to same format as O2 data
        presto_data = []
        rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        
        for row_idx, row_name in enumerate(rows):
            if row_name in df.index:
                for col_idx in range(24):
                    well_name = f"{row_name}{col_idx+1:02d}"
                    col_name = str(col_idx + 1)  # Column names are 1-24
                    
                    if col_name in df.columns:
                        value = df.loc[row_name, col_name]
                        
                        # Convert to numeric, handling 'NA' strings
                        if pd.notna(value) and str(value).strip().upper() != 'NA':
                            try:
                                numeric_value = float(value)
                                presto_data.append({
                                    'Well': well_name,
                                    'Row': row_name,
                                    'Column': col_idx + 1,
                                    'Presto_Blue_Value': numeric_value
                                })
                            except ValueError:
                                pass
        
        return pd.DataFrame(presto_data)
    
    def merge_datasets(self, use_averaged=True):
        """Merge O2 and Presto Blue data"""
        o2_data = self.load_o2_data(use_averaged)
        presto_data = self.load_presto_blue_data()
        
        print(f"O2 data: {len(o2_data)} wells")
        print(f"Presto Blue data: {len(presto_data)} wells")
        
        # Merge on Well
        merged = pd.merge(o2_data, presto_data, on='Well', how='inner', suffixes=('_o2', '_pb'))
        
        # Filter out wells without both values
        merged = merged.dropna(subset=['O2_Value', 'Presto_Blue_Value'])
        
        print(f"Merged data with both values: {len(merged)} wells")
        
        return merged
    
    def normalize_to_dmso(self, merged_data):
        """Normalize O2 and Presto Blue values to DMSO controls"""
        # Calculate DMSO means
        dmso_data = merged_data[merged_data['Content'].str.contains('DMSO', na=False)]
        
        if len(dmso_data) == 0:
            print("Warning: No DMSO controls found!")
            return merged_data
        
        dmso_o2_mean = dmso_data['O2_Value'].mean()
        dmso_pb_mean = dmso_data['Presto_Blue_Value'].mean()
        
        print(f"\nDMSO control values:")
        print(f"  O2 mean: {dmso_o2_mean:.2f}")
        print(f"  Presto Blue mean: {dmso_pb_mean:.2f}")
        print(f"  n = {len(dmso_data)} DMSO wells")
        
        # Normalize values (as percentage of DMSO)
        merged_data['O2_Normalized'] = (merged_data['O2_Value'] / dmso_o2_mean) * 100
        merged_data['Presto_Blue_Normalized'] = (merged_data['Presto_Blue_Value'] / dmso_pb_mean) * 100
        
        # Add a flag for DMSO wells
        merged_data['is_DMSO'] = merged_data['Content'].str.contains('DMSO', na=False)
        
        return merged_data
    
    def calculate_correlations(self, merged_data, use_normalized=False, exclude_dmso=False):
        """Calculate correlation statistics using numpy"""
        # Filter out DMSO if requested
        if exclude_dmso and 'is_DMSO' in merged_data.columns:
            data_for_corr = merged_data[~merged_data['is_DMSO']]
        else:
            data_for_corr = merged_data
            
        if use_normalized:
            o2_values = data_for_corr['O2_Normalized'].values
            pb_values = data_for_corr['Presto_Blue_Normalized'].values
        else:
            o2_values = data_for_corr['O2_Value'].values
            pb_values = data_for_corr['Presto_Blue_Value'].values
        
        # Pearson correlation using numpy
        correlation_matrix = np.corrcoef(o2_values, pb_values)
        pearson_r = correlation_matrix[0, 1]
        
        # Linear regression using numpy
        slope, intercept = np.polyfit(o2_values, pb_values, 1)
        
        # Calculate R-squared
        y_pred = slope * o2_values + intercept
        ss_res = np.sum((pb_values - y_pred) ** 2)
        ss_tot = np.sum((pb_values - np.mean(pb_values)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        
        stats_dict = {
            'pearson_r': pearson_r,
            'linear_slope': slope,
            'linear_intercept': intercept,
            'linear_r_squared': r_squared,
            'n_wells': len(data_for_corr)
        }
        
        return stats_dict
    
    def create_correlation_plot(self, merged_data, stats_dict, use_averaged=True, use_normalized=False):
        """Create scatter plot with correlation statistics"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
        
        if use_normalized:
            o2_values = merged_data['O2_Normalized'].values
            pb_values = merged_data['Presto_Blue_Normalized'].values
        else:
            o2_values = merged_data['O2_Value'].values
            pb_values = merged_data['Presto_Blue_Value'].values
        
        # Main scatter plot
        colors = []
        for _, row in merged_data.iterrows():
            if 'DMSO' in str(row['Content']):
                colors.append('red')
            elif 'DACTI' in str(row['Content']):
                colors.append('blue')
            elif 'PLICA' in str(row['Content']):
                colors.append('green')
            elif 'GEMCI' in str(row['Content']):
                colors.append('orange')
            elif 'VINC' in str(row['Content']):
                colors.append('purple')
            else:
                colors.append('gray')
        
        scatter = ax1.scatter(o2_values, pb_values, c=colors, alpha=0.7, s=50)
        
        # Add regression line
        x_line = np.linspace(o2_values.min(), o2_values.max(), 100)
        y_line = stats_dict['linear_slope'] * x_line + stats_dict['linear_intercept']
        ax1.plot(x_line, y_line, 'black', linestyle='--', alpha=0.8, linewidth=2)
        
        if use_normalized:
            ax1.set_xlabel('O₂ (% of DMSO control)', fontsize=12)
            ax1.set_ylabel('Presto Blue (% of DMSO control)', fontsize=12)
        else:
            ax1.set_xlabel('O₂ (% air saturation)', fontsize=12)
            ax1.set_ylabel('Presto Blue Value', fontsize=12)
        
        timepoint_desc = "avg last 10h" if use_averaged else "latest timepoint"
        norm_desc = " (DMSO-normalized)" if use_normalized else ""
        ax1.set_title(f'Cell Viability vs Oxygen Consumption{norm_desc}\\n({timepoint_desc})', fontsize=14)
        
        # Add statistics text
        stats_text = f"""Pearson r = {stats_dict['pearson_r']:.3f}
R² = {stats_dict['linear_r_squared']:.3f}
n = {stats_dict['n_wells']} wells

Negative correlation confirms:
Higher viability → More O₂ consumption"""
        
        ax1.text(0.05, 0.95, stats_text, transform=ax1.transAxes, 
                fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Add legend
        legend_elements = [
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=8, label='DMSO'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=8, label='DACTI'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=8, label='PLICA'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', markersize=8, label='GEMCI'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='purple', markersize=8, label='VINC')
        ]
        ax1.legend(handles=legend_elements, loc='lower right')
        
        # Drug-specific correlation analysis
        drug_stats = {}
        drug_colors = {'DMSO': 'red', 'DACTI': 'blue', 'PLICA': 'green', 'GEMCI': 'orange', 'VINC': 'purple'}
        
        y_pos = 0.95
        title_text = 'Drug-specific correlations:' if not use_normalized else 'Drug-specific correlations (vs DMSO baseline):'
        ax2.text(0.05, y_pos, title_text, transform=ax2.transAxes, 
                fontsize=12, fontweight='bold')
        y_pos -= 0.08
        
        # For normalized data, optionally skip DMSO or label it differently
        drugs_to_analyze = drug_colors.items()
        
        for drug, color in drugs_to_analyze:
            drug_data = merged_data[merged_data['Content'].str.contains(drug, na=False)]
            if len(drug_data) >= 3:  # Need at least 3 points for meaningful correlation
                if use_normalized:
                    drug_o2 = drug_data['O2_Normalized'].values
                    drug_pb = drug_data['Presto_Blue_Normalized'].values
                else:
                    drug_o2 = drug_data['O2_Value'].values
                    drug_pb = drug_data['Presto_Blue_Value'].values
                
                # Calculate correlation using numpy
                drug_corr_matrix = np.corrcoef(drug_o2, drug_pb)
                drug_r = drug_corr_matrix[0, 1] if not np.isnan(drug_corr_matrix[0, 1]) else 0
                drug_stats[drug] = {'r': drug_r, 'n': len(drug_data)}
                
                # Plot drug-specific data
                ax2.scatter(drug_o2, drug_pb, c=color, alpha=0.7, s=50, label=f'{drug} (n={len(drug_data)})')
                
                # Add text
                if use_normalized and drug == 'DMSO':
                    text = f'{drug}: r = {drug_r:.3f} (n={len(drug_data)}) *control variance'
                else:
                    text = f'{drug}: r = {drug_r:.3f} (n={len(drug_data)})'
                ax2.text(0.05, y_pos, text, 
                        transform=ax2.transAxes, fontsize=10, color=color)
                y_pos -= 0.06
        
        if use_normalized:
            ax2.set_xlabel('O₂ (% of DMSO control)', fontsize=12)
            ax2.set_ylabel('Presto Blue (% of DMSO control)', fontsize=12)
        else:
            ax2.set_xlabel('O₂ (% air saturation)', fontsize=12)
            ax2.set_ylabel('Presto Blue Value', fontsize=12)
        ax2.set_title('Drug-specific Correlations', fontsize=14)
        
        plt.tight_layout()
        
        # Save plot
        suffix = "_avg_last_10h" if use_averaged else "_latest"
        norm_suffix = "_normalized" if use_normalized else ""
        output_file = self.figures_dir / f'o2_presto_blue_correlation{suffix}{norm_suffix}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Created correlation plot: {output_file}")
        
        return drug_stats
    
    def export_correlation_results(self, merged_data, overall_stats, drug_stats, use_averaged=True, use_normalized=False):
        """Export correlation results to CSV"""
        # Overall results
        overall_df = pd.DataFrame([overall_stats])
        suffix = "_avg_last_10h" if use_averaged else "_latest"
        norm_suffix = "_normalized" if use_normalized else ""
        
        overall_file = self.results_dir / f'o2_presto_blue_overall_correlation{suffix}{norm_suffix}.csv'
        overall_df.to_csv(overall_file, index=False)
        
        # Drug-specific results
        drug_df = pd.DataFrame.from_dict(drug_stats, orient='index')
        drug_df.index.name = 'Drug'
        drug_file = self.results_dir / f'o2_presto_blue_drug_correlations{suffix}{norm_suffix}.csv'
        drug_df.to_csv(drug_file)
        
        # Full merged dataset
        merged_file = self.results_dir / f'o2_presto_blue_merged_data{suffix}{norm_suffix}.csv'
        merged_data.to_csv(merged_file, index=False)
        
        print(f"Exported correlation results:")
        print(f"  Overall: {overall_file}")
        print(f"  By drug: {drug_file}")
        print(f"  Full data: {merged_file}")
    
    def create_summary_plot(self, merged_data, use_averaged=True):
        """Create a summary plot showing assay concordance"""
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Calculate percentage changes relative to DMSO
        dmso_o2 = merged_data[merged_data['Content'].str.contains('DMSO', na=False)]['O2_Value'].mean()
        dmso_pb = merged_data[merged_data['Content'].str.contains('DMSO', na=False)]['Presto_Blue_Value'].mean()
        
        # Group by drug and calculate means
        drug_means = []
        drugs = ['DMSO', 'DACTI', 'PLICA', 'GEMCI', 'VINC']
        colors = ['red', 'blue', 'green', 'orange', 'purple']
        
        for drug, color in zip(drugs, colors):
            drug_data = merged_data[merged_data['Content'].str.contains(drug, na=False)]
            if len(drug_data) > 0:
                o2_mean = drug_data['O2_Value'].mean()
                pb_mean = drug_data['Presto_Blue_Value'].mean()
                
                # Calculate percent change from DMSO
                o2_change = ((o2_mean - dmso_o2) / dmso_o2) * 100
                pb_change = ((pb_mean - dmso_pb) / dmso_pb) * 100
                
                drug_means.append({
                    'drug': drug,
                    'o2_change': o2_change,
                    'pb_change': pb_change,
                    'color': color,
                    'n': len(drug_data)
                })
        
        # Plot
        for dm in drug_means:
            ax.scatter(dm['o2_change'], dm['pb_change'], 
                      s=200, c=dm['color'], alpha=0.7, edgecolors='black', linewidth=2)
            ax.annotate(f"{dm['drug']}\n(n={dm['n']})", 
                       (dm['o2_change'], dm['pb_change']), 
                       xytext=(5, 5), textcoords='offset points', fontsize=10)
        
        # Add diagonal line for perfect concordance
        ax_min = min(ax.get_xlim()[0], ax.get_ylim()[0])
        ax_max = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([ax_min, ax_max], [ax_min, ax_max], 'k--', alpha=0.3, label='Perfect concordance')
        
        # Add quadrant lines
        ax.axhline(y=0, color='gray', linestyle='-', alpha=0.3)
        ax.axvline(x=0, color='gray', linestyle='-', alpha=0.3)
        
        ax.set_xlabel('O₂ Change from DMSO (%)', fontsize=12)
        ax.set_ylabel('Presto Blue Change from DMSO (%)', fontsize=12)
        ax.set_title('Assay Concordance: O₂ Consumption vs Cell Viability', fontsize=14)
        
        # Add interpretation text
        ax.text(0.02, 0.98, 'Quadrants:\nBoth decrease = Drug effect\nDiscordant = Metabolic disruption', 
                transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        
        # Save
        timepoint = "avg_10h" if use_averaged else "latest"
        output_file = self.figures_dir / f'assay_concordance_summary_{timepoint}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Created assay concordance plot: {output_file}")
    
    def run_analysis(self, use_averaged=True, normalize=True):
        """Run complete correlation analysis"""
        print(f"=== O2 vs Presto Blue Correlation Analysis ===")
        timepoint_desc = "averaged last 10h" if use_averaged else "latest timepoint"
        print(f"Using O2 data: {timepoint_desc}")
        
        # Load and merge data
        merged_data = self.merge_datasets(use_averaged)
        
        if len(merged_data) == 0:
            print("No overlapping data found!")
            return
        
        # Normalize to DMSO if requested
        if normalize:
            merged_data = self.normalize_to_dmso(merged_data)
            print("\nData normalized to DMSO controls")
        
        # Calculate correlations
        overall_stats = self.calculate_correlations(merged_data, use_normalized=normalize)
        
        # Create plots
        drug_stats = self.create_correlation_plot(merged_data, overall_stats, use_averaged, use_normalized=normalize)
        
        # Create summary concordance plot (only for non-normalized data)
        if not normalize:
            self.create_summary_plot(merged_data, use_averaged)
        
        # Export results
        self.export_correlation_results(merged_data, overall_stats, drug_stats, use_averaged, use_normalized=normalize)
        
        # Print summary
        norm_desc = " (DMSO-normalized)" if normalize else ""
        print(f"\n=== CORRELATION SUMMARY{norm_desc} ===")
        print(f"Overall correlation (n={overall_stats['n_wells']} wells):")
        print(f"  Pearson r = {overall_stats['pearson_r']:.3f}")
        print(f"  R² = {overall_stats['linear_r_squared']:.3f}")
        
        # Interpret the correlation
        if overall_stats['pearson_r'] < -0.7:
            print("  → Strong negative correlation: Viability and O₂ consumption are tightly linked")
        elif overall_stats['pearson_r'] < -0.4:
            print("  → Moderate negative correlation: Viability and O₂ consumption are related")
        
        print(f"\nDrug-specific correlations:")
        for drug, stats in drug_stats.items():
            correlation_strength = ""
            if stats['r'] < -0.7:
                correlation_strength = " (maintains strong viability-O₂ relationship)"
            elif stats['r'] > -0.3 and stats['r'] < 0.3:
                correlation_strength = " (disrupts viability-O₂ relationship)"
            print(f"  {drug}: r = {stats['r']:.3f} (n={stats['n']}){correlation_strength}")

def main():
    analyzer = O2PrestoBlueCorrelation("/Users/shaunie/Desktop/hu-glioblastoma/data/revision")
    
    # Run analysis with DMSO normalization
    print("=== DMSO-NORMALIZED ANALYSIS ===")
    print("\n--- AVERAGED DATA ---")
    analyzer.run_analysis(use_averaged=True, normalize=True)
    
    print("\n--- LATEST TIMEPOINT ---")
    analyzer.run_analysis(use_averaged=False, normalize=True)
    
    print("\n" + "="*80 + "\n")
    
    # Run analysis without normalization for comparison
    print("=== RAW DATA ANALYSIS (NO NORMALIZATION) ===")
    print("\n--- AVERAGED DATA ---")
    analyzer.run_analysis(use_averaged=True, normalize=False)
    
    print("\n--- LATEST TIMEPOINT ---")
    analyzer.run_analysis(use_averaged=False, normalize=False)

if __name__ == "__main__":
    main()