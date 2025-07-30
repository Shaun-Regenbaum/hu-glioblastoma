#!/usr/bin/env python3
"""
Create a heatmap that matches the actual 96-well plate layout
Final version that properly maps data columns to wells based on layout content
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict

class FinalPlateHeatmapGenerator:
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
    
    def create_smart_well_to_column_mapping(self, data):
        """Create mapping using content-based matching"""
        layout = self.load_expanded_layout()
        
        # Group data columns by drug and concentration
        drug_columns = defaultdict(lambda: defaultdict(list))
        dmso_columns = defaultdict(list)
        
        for col in data.columns:
            if col == 'Time_h':
                continue
            
            if col.startswith('DMSO'):
                parts = col.split('_')
                conc = float(parts[1])
                dmso_columns[conc].append(col)
            else:
                parts = col.split('_')
                if len(parts) >= 3:
                    drug = parts[0]
                    conc = float(parts[1])
                    drug_columns[drug][conc].append(col)
        
        print(f"Found drug columns: {dict(drug_columns)}")
        print(f"Found DMSO columns: {dict(dmso_columns)}")
        
        # Create well-to-column mapping
        well_to_column = {}
        
        # Process each cell in the layout
        for row_name, row_data in layout.iterrows():
            for col_idx, cell_content in enumerate(row_data):
                if pd.notna(cell_content) and cell_content.strip():
                    well_name = f"{row_name}{col_idx+1:02d}"
                    
                    # Parse the cell content
                    if cell_content.startswith('DMSO'):
                        # DMSO cell: extract concentration
                        conc_str = cell_content.split('_')[1]
                        conc = float(conc_str)
                        
                        # Find available DMSO column for this concentration
                        if conc in dmso_columns and dmso_columns[conc]:
                            column = dmso_columns[conc].pop(0)
                            well_to_column[well_name] = column
                            
                    else:
                        # Drug cell: extract drug and concentration
                        parts = cell_content.split('_')
                        drug = parts[0]
                        conc = float(parts[1])
                        
                        # Find available column for this drug and concentration
                        if drug in drug_columns and conc in drug_columns[drug]:
                            if drug_columns[drug][conc]:
                                column = drug_columns[drug][conc].pop(0)
                                well_to_column[well_name] = column
        
        print(f"Successfully mapped {len(well_to_column)} wells to columns")
        
        # Show some example mappings
        rows = ['B', 'C', 'D', 'E', 'F', 'G', 'H']
        for row in rows[:2]:  # Show first 2 rows as examples
            for col in [1, 21, 22, 23, 24]:  # Show first column and DMSO columns
                well_name = f"{row}{col:02d}"
                if well_name in well_to_column:
                    print(f"  {well_name} -> {well_to_column[well_name]}")
        
        return well_to_column
    
    def create_plate_heatmap(self, use_last_hours=10, average_mode=True):
        """Create heatmap showing individual wells with correct mapping"""
        data = self.load_time_series_data()
        
        if average_mode:
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
            # Use the latest timepoint
            latest_idx = data['Time_h'].idxmax()
            actual_time = data.iloc[latest_idx]['Time_h']
            print(f"Using latest timepoint: {actual_time}h")
            
            def get_well_value(column_name):
                return data.iloc[latest_idx][column_name]
            
            time_description = f"{actual_time}h"
            filename_suffix = f"_latest_{actual_time}h"
        
        # Create the smart mapping
        well_to_column = self.create_smart_well_to_column_mapping(data)
        
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
        
        # Highlight DMSO columns (21-24)
        for dmso_col in [20, 21, 22, 23]:  # 0-indexed, so 21-24 becomes 20-23
            ax.axvline(x=dmso_col+0.5, color='red', linestyle='--', alpha=0.7, linewidth=2)
        
        plt.tight_layout()
        
        # Save the plot
        output_file = self.figures_dir / f'final_plate_layout_heatmap{filename_suffix}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Created final plate layout heatmap: {output_file}")
        
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
        output_file = self.results_dir / f'final_plate_layout_with_values{filename_suffix}.csv'
        plate_df.to_csv(output_file, index=False)
        
        print(f"Exported final plate layout CSV: {output_file}")

def main():
    generator = FinalPlateHeatmapGenerator("/Users/shaunie/Desktop/hu-glioblastoma/data/revision")
    
    # Create both versions
    print("=== Creating averaged version (maximum data utilization) ===")
    generator.create_plate_heatmap(use_last_hours=10, average_mode=True)
    
    print("\n=== Creating latest timepoint version ===")
    generator.create_plate_heatmap(use_last_hours=10, average_mode=False)

if __name__ == "__main__":
    main()