#!/usr/bin/env python3
"""
Create endpoint heatmap matching the actual plate layout
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

def load_time_series_data():
    """Load the time series data"""
    file_path = "data/revision/results/all_drugs_time_series.csv"
    return pd.read_csv(file_path)

def load_plate_layout():
    """Load the dose response brain orgs map"""
    file_path = "data/revision/dose reponse brain orgs map.csv"
    return pd.read_csv(file_path)

def create_plate_layout_heatmap():
    """Create heatmap matching actual plate layout"""
    data = load_time_series_data()
    layout = load_plate_layout()
    
    # Use final timepoint
    final_time_idx = len(data) - 1
    actual_time = data.iloc[final_time_idx]['Time_h']
    print(f"Using final time point {actual_time}h")
    
    # Create mapping from well position to oxygen value
    # Based on the layout: rows are a,b,c,d and columns are 1-12
    plate_data = np.full((4, 12), np.nan)  # 4 rows, 12 columns
    
    # Extract drug info from layout
    rows = ['a', 'b', 'c', 'd']
    drugs = ['dacti', 'plica', 'gemci', 'vinc']
    
    # Get concentrations from layout (columns 1-12)
    concentrations = []
    for col_idx in range(1, 13):  # columns 1-12
        if col_idx < len(layout.columns):
            conc_val = layout.iloc[0, col_idx]  # First row has concentrations
            if pd.notna(conc_val) and str(conc_val).strip():
                concentrations.append(float(conc_val))
            else:
                concentrations.append(np.nan)
        else:
            concentrations.append(np.nan)
    
    print(f"Concentrations: {concentrations}")
    
    # Fill the plate data
    for row_idx, (row_letter, drug) in enumerate(zip(rows, drugs)):
        for col_idx in range(12):
            if col_idx < len(concentrations) and not np.isnan(concentrations[col_idx]):
                conc = concentrations[col_idx]
                
                # Find matching wells in time series data
                # Look for wells like DACTI_0.002222222_1, DACTI_0.002222222_2, etc.
                drug_upper = drug.upper()
                matching_wells = []
                for well_name in data.columns:
                    if well_name.startswith(f"{drug_upper}_"):
                        parts = well_name.split('_')
                        if len(parts) >= 3:
                            well_conc = float(parts[1])
                            if abs(well_conc - conc) < 0.0001:  # Close match for floating point
                                matching_wells.append(well_name)
                
                # Calculate mean oxygen value for this position
                if matching_wells:
                    values = []
                    for well in matching_wells:
                        val = data.iloc[final_time_idx][well]
                        if not pd.isna(val):
                            values.append(val)
                    
                    if values:
                        plate_data[row_idx, col_idx] = np.mean(values)
                        print(f"Position {row_letter}{col_idx+1}: {drug_upper} {conc}μM = {np.mean(values):.1f}% (from {len(values)} wells)")
    
    # Create the heatmap
    fig, ax = plt.subplots(figsize=(14, 6))
    
    # Create heatmap
    im = ax.imshow(plate_data, aspect='auto', cmap='viridis_r', interpolation='nearest')
    
    # Set labels
    ax.set_xticks(range(12))
    ax.set_xticklabels([f'{i+1}' for i in range(12)])
    ax.set_xlabel('Column', fontsize=12)
    
    ax.set_yticks(range(4))
    ax.set_yticklabels(['A (DACTI)', 'B (PLICA)', 'C (GEMCI)', 'D (VINC)'])
    ax.set_ylabel('Row (Drug)', fontsize=12)
    
    ax.set_title(f'Endpoint O₂ Values - Plate Layout ({actual_time:.0f}h)', fontsize=16, fontweight='bold')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('O₂ (% air saturation)', fontsize=12)
    
    # Add text annotations
    for i in range(4):
        for j in range(12):
            if not np.isnan(plate_data[i, j]):
                text = ax.text(j, i, f'{plate_data[i, j]:.1f}',
                             ha="center", va="center", color="white", fontweight='bold')
    
    plt.tight_layout()
    
    # Save
    output_path = Path("data/revision/figures/endpoint_heatmap_plate_layout.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Created plate layout heatmap: {output_path}")
    
    # Print stats
    filled_cells = np.sum(~np.isnan(plate_data))
    total_cells = plate_data.size
    print(f"Filled {filled_cells}/{total_cells} cells ({filled_cells/total_cells*100:.1f}%)")

def main():
    print("Creating plate layout endpoint heatmap...")
    create_plate_layout_heatmap()

if __name__ == "__main__":
    main()