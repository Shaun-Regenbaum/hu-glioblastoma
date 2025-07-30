#!/usr/bin/env python3
"""
Create a heatmap that matches the actual 96-well plate layout
Each cell in the heatmap represents one well (B1, B2, C1, etc.)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

class PlateHeatmapGenerator:
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
        """Create mapping from well positions to data columns"""
        well_to_column = {}
        
        # Process drug columns
        for col in data.columns:
            if col != 'Time_h' and not col.startswith('DMSO'):
                parts = col.split('_')
                if len(parts) >= 3:
                    drug = parts[0]
                    conc = float(parts[1])
                    replicate = int(parts[2])
                    
                    # Find the corresponding well based on our layout logic
                    # This maps back from the column naming to actual well positions
                    well_to_column[col] = {'drug': drug, 'concentration': conc, 'replicate': replicate}
        
        # Process DMSO columns
        for col in data.columns:
            if col.startswith('DMSO'):
                parts = col.split('_')
                if len(parts) >= 3:
                    conc = float(parts[1])
                    replicate = int(parts[2])
                    well_to_column[col] = {'drug': 'DMSO', 'concentration': conc, 'replicate': replicate}
        
        return well_to_column
    
    def map_columns_to_wells(self, data):
        """Map data columns to actual well positions using the expanded layout"""
        layout = self.load_expanded_layout()
        
        # Create a mapping from layout content to data columns
        content_to_columns = {}
        
        # Process drug columns to group by drug and concentration
        for col in data.columns:
            if col != 'Time_h':
                if col.startswith('DMSO'):
                    parts = col.split('_')
                    conc = float(parts[1])
                    content_key = f"DMSO_{conc}"
                else:
                    parts = col.split('_')
                    if len(parts) >= 3:
                        drug = parts[0]
                        conc = float(parts[1])
                        content_key = f"{drug}_{conc}"
                
                if content_key not in content_to_columns:
                    content_to_columns[content_key] = []
                content_to_columns[content_key].append(col)
        
        # Now map wells to specific columns based on the layout
        well_to_column = {}
        
        for row_name, row_data in layout.iterrows():
            for col_idx, cell_content in enumerate(row_data):
                if pd.notna(cell_content) and cell_content:
                    well_name = f"{row_name}{col_idx+1:02d}"
                    
                    # Find matching columns for this cell content
                    if cell_content in content_to_columns:
                        available_columns = content_to_columns[cell_content]
                        
                        # Assign columns in order (first available column for first well, etc.)
                        if available_columns:
                            assigned_column = available_columns.pop(0)
                            well_to_column[well_name] = assigned_column
        
        return well_to_column
    
    def create_plate_heatmap(self, use_last_hours=10, average_mode=True, specific_timepoint=None):
        """Create heatmap showing individual wells
        
        Args:
            use_last_hours: Number of hours from the end to consider (default 10)
            average_mode: If True, average across the last N hours. If False, use latest timepoint
            specific_timepoint: If provided, use this specific timepoint instead
        """
        data = self.load_time_series_data()
        layout = self.load_expanded_layout()
        
        if specific_timepoint is not None:
            # Use specific timepoint
            time_diff = (data['Time_h'] - specific_timepoint).abs()
            timepoint_idx = time_diff.idxmin()
            actual_time = data.iloc[timepoint_idx]['Time_h']
            print(f"Using specific time point {actual_time}h (target: {specific_timepoint}h)")
            
            # Get values at this timepoint
            def get_well_value(column_name):
                return data.iloc[timepoint_idx][column_name]
            
            time_description = f"{actual_time}h"
            
        elif average_mode:
            # Average across last N hours
            max_time = data['Time_h'].max()
            min_time_for_avg = max_time - use_last_hours
            
            # Filter to last N hours
            last_hours_mask = data['Time_h'] >= min_time_for_avg
            last_hours_data = data[last_hours_mask]
            
            print(f"Averaging across last {use_last_hours} hours ({min_time_for_avg:.1f}h to {max_time:.1f}h)")
            print(f"Using {len(last_hours_data)} timepoints for averaging")
            
            # Get average values
            def get_well_value(column_name):
                values = last_hours_data[column_name].dropna()
                return values.mean() if len(values) > 0 else np.nan
            
            time_description = f"avg of last {use_last_hours}h"
            
        else:
            # Use the latest timepoint within last N hours
            max_time = data['Time_h'].max() 
            min_time_for_selection = max_time - use_last_hours
            
            # Filter to last N hours and get the latest
            last_hours_mask = data['Time_h'] >= min_time_for_selection
            last_hours_data = data[last_hours_mask]
            latest_idx = last_hours_data['Time_h'].idxmax()
            
            actual_time = data.iloc[latest_idx]['Time_h']
            print(f"Using latest timepoint in last {use_last_hours} hours: {actual_time}h")
            
            # Get values at this timepoint
            def get_well_value(column_name):
                return data.iloc[latest_idx][column_name]
            
            time_description = f"{actual_time}h"
        
        # Map wells to data columns
        well_to_column = self.map_columns_to_wells(data)
        
        # Create the heatmap data matrix (8 rows x 24 columns for 96-well plate)
        plate_data = np.full((8, 24), np.nan)
        well_labels = np.full((8, 24), '', dtype=object)
        
        rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        
        for row_idx, row_name in enumerate(rows):
            for col_idx in range(24):
                well_name = f"{row_name}{col_idx+1:02d}"
                
                if well_name in well_to_column:
                    column_name = well_to_column[well_name]
                    value = get_well_value(column_name)
                    
                    if not pd.isna(value):
                        plate_data[row_idx, col_idx] = value
                        
                        # Create label for the well
                        if column_name.startswith('DMSO'):
                            well_labels[row_idx, col_idx] = 'DMSO'
                        else:
                            parts = column_name.split('_')
                            drug = parts[0]
                            conc = float(parts[1])
                            well_labels[row_idx, col_idx] = f"{drug}\n{conc}"
        
        # Create the heatmap
        fig, ax = plt.subplots(figsize=(20, 8))
        
        # Use a colormap that handles NaN values
        mask = np.isnan(plate_data)
        
        # Create heatmap
        im = ax.imshow(plate_data, cmap='viridis_r', aspect='equal', interpolation='nearest')
        
        # Set up the axes
        ax.set_xticks(range(24))
        ax.set_xticklabels([f'{i+1:02d}' for i in range(24)])
        ax.set_xlabel('Column', fontsize=14)
        
        ax.set_yticks(range(8))
        ax.set_yticklabels(rows)
        ax.set_ylabel('Row', fontsize=14)
        
        # Add well labels
        for i in range(8):
            for j in range(24):
                if not mask[i, j]:
                    # Add the drug/concentration label
                    ax.text(j, i, well_labels[i, j], ha='center', va='center', 
                           fontsize=8, color='white', weight='bold')
                    
                    # Add the oxygen value
                    ax.text(j, i+0.3, f'{plate_data[i, j]:.1f}', ha='center', va='center',
                           fontsize=7, color='white')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('O₂ (% air saturation)', fontsize=12)
        
        # Set title
        ax.set_title(f'96-Well Plate Layout - O₂ Values ({time_description})', 
                    fontsize=16, fontweight='bold', pad=20)
        
        # Add grid
        ax.set_xticks(np.arange(-0.5, 24, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, 8, 1), minor=True)
        ax.grid(which='minor', color='white', linestyle='-', linewidth=0.5)
        
        plt.tight_layout()
        
        # Save the plot with descriptive filename
        if specific_timepoint is not None:
            filename_suffix = f"_timepoint_{specific_timepoint}h"
        elif average_mode:
            filename_suffix = f"_avg_last_{use_last_hours}h"
        else:
            filename_suffix = f"_latest_in_{use_last_hours}h"
            
        output_file = self.figures_dir / f'plate_layout_heatmap{filename_suffix}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Created plate layout heatmap: {output_file}")
        
        # Also create a CSV of the plate layout with values
        self.export_plate_layout_csv(plate_data, well_labels, time_description, filename_suffix)
    
    def export_plate_layout_csv(self, plate_data, well_labels, timepoint, filename_suffix=""):
        """Export the plate layout data as CSV"""
        rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        
        # Create DataFrame with well names and values
        plate_df_data = []
        for row_idx, row_name in enumerate(rows):
            for col_idx in range(24):
                well_name = f"{row_name}{col_idx+1:02d}"
                value = plate_data[row_idx, col_idx]
                label = well_labels[row_idx, col_idx]
                
                plate_df_data.append({
                    'Well': well_name,
                    'Row': row_name,
                    'Column': col_idx + 1,
                    'Content': label if label else 'Empty',
                    'O2_Value': value if not np.isnan(value) else None,
                    'Timepoint_h': timepoint
                })
        
        plate_df = pd.DataFrame(plate_df_data)
        output_file = self.results_dir / f'plate_layout_with_values{filename_suffix}.csv'
        plate_df.to_csv(output_file, index=False)
        
        print(f"Exported plate layout CSV: {output_file}")

def main():
    generator = PlateHeatmapGenerator("/Users/shaunie/Desktop/hu-glioblastoma/data/revision")
    
    # Use maximum data - average across last 10 hours for best signal
    generator.create_plate_heatmap(use_last_hours=10, average_mode=True)
    
    # Also create one with just the latest timepoint for comparison
    generator.create_plate_heatmap(use_last_hours=10, average_mode=False)

if __name__ == "__main__":
    main()