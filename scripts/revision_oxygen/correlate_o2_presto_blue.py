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
    
    def calculate_correlations(self, merged_data):
        """Calculate correlation statistics using numpy"""
        o2_values = merged_data['O2_Value'].values
        pb_values = merged_data['Presto_Blue_Value'].values
        
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
            'n_wells': len(merged_data)
        }
        
        return stats_dict
    
    def create_correlation_plot(self, merged_data, stats_dict, use_averaged=True):
        """Create scatter plot with correlation statistics"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
        
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
        
        ax1.set_xlabel('O₂ (% air saturation)', fontsize=12)
        ax1.set_ylabel('Presto Blue Value', fontsize=12)
        
        timepoint_desc = "avg last 10h" if use_averaged else "latest timepoint"
        ax1.set_title(f'O₂ vs Presto Blue Correlation\\n({timepoint_desc})', fontsize=14)
        
        # Add statistics text
        stats_text = f"""Pearson r = {stats_dict['pearson_r']:.3f}
R² = {stats_dict['linear_r_squared']:.3f}
n = {stats_dict['n_wells']} wells"""
        
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
        ax2.text(0.05, y_pos, 'Drug-specific correlations:', transform=ax2.transAxes, 
                fontsize=12, fontweight='bold')
        y_pos -= 0.08
        
        for drug, color in drug_colors.items():
            drug_data = merged_data[merged_data['Content'].str.contains(drug, na=False)]
            if len(drug_data) >= 3:  # Need at least 3 points for meaningful correlation
                drug_o2 = drug_data['O2_Value'].values
                drug_pb = drug_data['Presto_Blue_Value'].values
                
                # Calculate correlation using numpy
                drug_corr_matrix = np.corrcoef(drug_o2, drug_pb)
                drug_r = drug_corr_matrix[0, 1] if not np.isnan(drug_corr_matrix[0, 1]) else 0
                drug_stats[drug] = {'r': drug_r, 'n': len(drug_data)}
                
                # Plot drug-specific data
                ax2.scatter(drug_o2, drug_pb, c=color, alpha=0.7, s=50, label=f'{drug} (n={len(drug_data)})')
                
                # Add text
                ax2.text(0.05, y_pos, f'{drug}: r = {drug_r:.3f} (n={len(drug_data)})', 
                        transform=ax2.transAxes, fontsize=10, color=color)
                y_pos -= 0.06
        
        ax2.set_xlabel('O₂ (% air saturation)', fontsize=12)
        ax2.set_ylabel('Presto Blue Value', fontsize=12)
        ax2.set_title('Drug-specific Correlations', fontsize=14)
        
        plt.tight_layout()
        
        # Save plot
        suffix = "_avg_last_10h" if use_averaged else "_latest"
        output_file = self.figures_dir / f'o2_presto_blue_correlation{suffix}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Created correlation plot: {output_file}")
        
        return drug_stats
    
    def export_correlation_results(self, merged_data, overall_stats, drug_stats, use_averaged=True):
        """Export correlation results to CSV"""
        # Overall results
        overall_df = pd.DataFrame([overall_stats])
        suffix = "_avg_last_10h" if use_averaged else "_latest"
        
        overall_file = self.results_dir / f'o2_presto_blue_overall_correlation{suffix}.csv'
        overall_df.to_csv(overall_file, index=False)
        
        # Drug-specific results
        drug_df = pd.DataFrame.from_dict(drug_stats, orient='index')
        drug_df.index.name = 'Drug'
        drug_file = self.results_dir / f'o2_presto_blue_drug_correlations{suffix}.csv'
        drug_df.to_csv(drug_file)
        
        # Full merged dataset
        merged_file = self.results_dir / f'o2_presto_blue_merged_data{suffix}.csv'
        merged_data.to_csv(merged_file, index=False)
        
        print(f"Exported correlation results:")
        print(f"  Overall: {overall_file}")
        print(f"  By drug: {drug_file}")
        print(f"  Full data: {merged_file}")
    
    def run_analysis(self, use_averaged=True):
        """Run complete correlation analysis"""
        print(f"=== O2 vs Presto Blue Correlation Analysis ===")
        timepoint_desc = "averaged last 10h" if use_averaged else "latest timepoint"
        print(f"Using O2 data: {timepoint_desc}")
        
        # Load and merge data
        merged_data = self.merge_datasets(use_averaged)
        
        if len(merged_data) == 0:
            print("No overlapping data found!")
            return
        
        # Calculate correlations
        overall_stats = self.calculate_correlations(merged_data)
        
        # Create plots
        drug_stats = self.create_correlation_plot(merged_data, overall_stats, use_averaged)
        
        # Export results
        self.export_correlation_results(merged_data, overall_stats, drug_stats, use_averaged)
        
        # Print summary
        print(f"\n=== CORRELATION SUMMARY ===")
        print(f"Overall correlation (n={overall_stats['n_wells']} wells):")
        print(f"  Pearson r = {overall_stats['pearson_r']:.3f}")
        print(f"  R² = {overall_stats['linear_r_squared']:.3f}")
        
        print(f"\nDrug-specific correlations:")
        for drug, stats in drug_stats.items():
            print(f"  {drug}: r = {stats['r']:.3f} (n={stats['n']})")

def main():
    analyzer = O2PrestoBlueCorrelation("/Users/shaunie/Desktop/hu-glioblastoma/data/revision")
    
    # Run analysis for both averaged and latest timepoint data
    print("=== AVERAGED DATA ANALYSIS ===")
    analyzer.run_analysis(use_averaged=True)
    
    print("\n" + "="*60 + "\n")
    
    print("=== LATEST TIMEPOINT ANALYSIS ===")
    analyzer.run_analysis(use_averaged=False)

if __name__ == "__main__":
    main()