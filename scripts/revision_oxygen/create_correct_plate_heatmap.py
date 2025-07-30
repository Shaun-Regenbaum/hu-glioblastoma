#!/usr/bin/env python3
"""
Create a heatmap that matches the actual 96-well plate layout
Fixed version that properly maps all data columns to wells
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

class CorrectPlateHeatmapGenerator:
    def __init__(self, base_path):
        self.base_path = Path(base_path)
        self.results_dir = self.base_path / "results"
        self.figures_dir = self.base_path / "figures"
        
    def load_expanded_layout(self):
        """Load the expanded plate layout"""
        layout_file = self.base_path / "full_expanded_plate_layout.csv"
        return pd.read_csv(layout_file, index_col=0)
    
    def load_time_series_data(self):
        """Load the time series data"""
        file_path = self.results_dir / "all_drugs_time_series.csv"
        return pd.read_csv(file_path)
    
    def create_well_to_column_mapping(self, data):
        """Create direct mapping from wells to data columns based on the layout order"""
        layout = self.load_expanded_layout()
        
        # Get all data columns (excluding Time_h)
        data_columns = [col for col in data.columns if col != 'Time_h']
        
        # Create ordered list of wells with content from the layout
        wells_with_content = []
        
        for row_name, row_data in layout.iterrows():
            for col_idx, cell_content in enumerate(row_data):
                if pd.notna(cell_content) and cell_content.strip():
                    well_name = f"{row_name}{col_idx+1:02d}"
                    wells_with_content.append((well_name, cell_content))
        
        print(f"Found {len(wells_with_content)} wells with content in layout")
        print(f"Found {len(data_columns)} data columns")
        
        # Map wells to columns in order
        well_to_column = {}
        
        # Debug: Show first few mappings
        for i, (well_name, content) in enumerate(wells_with_content):
            if i < len(data_columns):
                column_name = data_columns[i]
                well_to_column[well_name] = column_name
                if i < 10:  # Show first 10 for debugging
                    print(f"  {well_name} ({content}) -> {column_name}")
            
        print(f"Successfully mapped {len(well_to_column)} wells to columns")
        return well_to_column
    
    def create_plate_heatmap(self, use_last_hours=10, average_mode=True, specific_timepoint=None):
        """Create heatmap showing individual wells"""
        data = self.load_time_series_data()
        
        if specific_timepoint is not None:
            # Use specific timepoint
            time_diff = (data['Time_h'] - specific_timepoint).abs()
            timepoint_idx = time_diff.idxmin()
            actual_time = data.iloc[timepoint_idx]['Time_h']
            print(f"Using specific time point {actual_time}h (target: {specific_timepoint}h)")
            
            def get_well_value(column_name):
                return data.iloc[timepoint_idx][column_name]
            
            time_description = f"{actual_time}h"
            filename_suffix = f"_timepoint_{specific_timepoint}h"
            
        elif average_mode:
            # Average across last N hours
            max_time = data['Time_h'].max()
            min_time_for_avg = max_time - use_last_hours
            
            last_hours_mask = data['Time_h'] >= min_time_for_avg
            last_hours_data = data[last_hours_mask]
            
            print(f"Averaging across last {use_last_hours} hours ({min_time_for_avg:.1f}h to {max_time:.1f}h)")
            print(f"Using {len(last_hours_data)} timepoints for averaging")
            
            def get_well_value(column_name):
                values = last_hours_data[column_name].dropna()
                return values.mean() if len(values) > 0 else np.nan
            
            time_description = f"avg of last {use_last_hours}h"
            filename_suffix = f"_avg_last_{use_last_hours}h"
            
        else:
            # Use the latest timepoint within last N hours
            max_time = data['Time_h'].max() 
            min_time_for_selection = max_time - use_last_hours
            
            last_hours_mask = data['Time_h'] >= min_time_for_selection
            last_hours_data = data[last_hours_mask]
            latest_idx = last_hours_data['Time_h'].idxmax()
            
            actual_time = data.iloc[latest_idx]['Time_h']
            print(f"Using latest timepoint in last {use_last_hours} hours: {actual_time}h")
            
            def get_well_value(column_name):
                return data.iloc[latest_idx][column_name]
            
            time_description = f"{actual_time}h"
            filename_suffix = f"_latest_in_{use_last_hours}h"
        
        # Create the mapping
        well_to_column = self.create_well_to_column_mapping(data)
        
        # Create the heatmap data matrix (8 rows x 24 columns)
        plate_data = np.full((8, 24), np.nan)
        well_labels = np.full((8, 24), '', dtype=object)
        
        rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        
        # Fill in the data
        filled_wells = 0
        for row_idx, row_name in enumerate(rows):
            for col_idx in range(24):
                well_name = f"{row_name}{col_idx+1:02d}"
                
                if well_name in well_to_column:
                    column_name = well_to_column[well_name]
                    value = get_well_value(column_name)
                    
                    if not pd.isna(value):
                        plate_data[row_idx, col_idx] = value
                        filled_wells += 1
                        
                        # Create simplified label
                        if column_name.startswith('DMSO'):
                            well_labels[row_idx, col_idx] = 'DMSO'
                        else:
                            parts = column_name.split('_')
                            drug = parts[0]
                            conc = float(parts[1])
                            well_labels[row_idx, col_idx] = f"{drug}\n{conc:.3g}"
        
        print(f"Filled {filled_wells} wells with data")
        
        # Create the heatmap
        fig, ax = plt.subplots(figsize=(24, 8))
        
        # Create heatmap with proper color handling for NaN
        mask = np.isnan(plate_data)
        
        # Use viridis_r colormap (dark = low O2, bright = high O2)
        im = ax.imshow(plate_data, cmap='viridis_r', aspect='equal', interpolation='nearest')
        
        # Set up the axes
        ax.set_xticks(range(24))
        ax.set_xticklabels([f'{i+1:02d}' for i in range(24)])
        ax.set_xlabel('Column', fontsize=14)
        
        ax.set_yticks(range(8))
        ax.set_yticklabels(rows)
        ax.set_ylabel('Row', fontsize=14)
        
        # Add text annotations for each cell
        for i in range(8):
            for j in range(24):
                if not mask[i, j]:
                    # Add the drug/concentration label
                    ax.text(j, i-0.15, well_labels[i, j], ha='center', va='center', 
                           fontsize=7, color='white', weight='bold')
                    
                    # Add the oxygen value
                    ax.text(j, i+0.15, f'{plate_data[i, j]:.1f}', ha='center', va='center',
                           fontsize=8, color='white', weight='bold')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.6)
        cbar.set_label('O₂ (% air saturation)', fontsize=12)
        
        # Set title
        ax.set_title(f'96-Well Plate Layout - O₂ Values ({time_description})', 
                    fontsize=16, fontweight='bold', pad=20)
        
        # Add grid
        ax.set_xticks(np.arange(-0.5, 24, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, 8, 1), minor=True)
        ax.grid(which='minor', color='white', linestyle='-', linewidth=0.5)
        
        plt.tight_layout()
        
        # Save the plot
        output_file = self.figures_dir / f'corrected_plate_layout_heatmap{filename_suffix}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Created corrected plate layout heatmap: {output_file}")
        
        # Export CSV with the mapping
        self.export_plate_layout_csv(plate_data, well_labels, well_to_column, time_description, filename_suffix)
    
    def export_plate_layout_csv(self, plate_data, well_labels, well_to_column, timepoint, filename_suffix=""):
        """Export the plate layout data as CSV with column mapping info"""
        rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        
        plate_df_data = []
        for row_idx, row_name in enumerate(rows):
            for col_idx in range(24):
                well_name = f"{row_name}{col_idx+1:02d}"
                value = plate_data[row_idx, col_idx]
                label = well_labels[row_idx, col_idx]
                
                # Get the data column name
                data_column = well_to_column.get(well_name, '')
                
                plate_df_data.append({
                    'Well': well_name,
                    'Row': row_name,
                    'Column': col_idx + 1,
                    'Content': label if label else 'Empty',
                    'Data_Column': data_column,
                    'O2_Value': value if not np.isnan(value) else None,
                    'Timepoint': timepoint
                })
        
        plate_df = pd.DataFrame(plate_df_data)
        output_file = self.results_dir / f'corrected_plate_layout_with_values{filename_suffix}.csv'
        plate_df.to_csv(output_file, index=False)
        
        print(f"Exported corrected plate layout CSV: {output_file}")

def main():
    generator = CorrectPlateHeatmapGenerator("/Users/shaunie/Desktop/hu-glioblastoma/data/revision")
    
    # Create both versions for maximum data
    print("=== Creating averaged version (maximum data utilization) ===")
    generator.create_plate_heatmap(use_last_hours=10, average_mode=True)
    
    print("\n=== Creating latest timepoint version ===")
    generator.create_plate_heatmap(use_last_hours=10, average_mode=False)

if __name__ == "__main__":
    main()