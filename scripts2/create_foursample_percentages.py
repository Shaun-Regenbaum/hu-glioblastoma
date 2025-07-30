#!/usr/bin/env python3
"""
Create FourSample-FourPop-Percentages.csv with statistics to match original Shaun structure
"""

import pandas as pd
import numpy as np

def create_foursample_percentages():
    """Create FourSample-FourPop-Percentages.csv with percentages and statistics"""
    
    print("Creating FourSample-FourPop-Percentages.csv with REAL data...")
    
    # Load the cell counts data
    counts_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/FourSample-FourPop-CellCounts.csv'
    df_counts = pd.read_csv(counts_path)
    
    # Convert to percentages
    df_pct = df_counts.set_index('sample')
    df_pct = df_pct.div(df_pct.sum(axis=1), axis=0) * 100
    
    # Calculate statistics for each cell type
    stats_data = []
    
    # Add the percentage rows
    for sample in df_pct.index:
        row_data = {'sample': sample}
        for cell_type in df_pct.columns:
            row_data[cell_type] = df_pct.loc[sample, cell_type]
        stats_data.append(row_data)
    
    # Calculate range (max - min) for each cell type
    range_row = {'sample': 'range'}
    for cell_type in df_pct.columns:
        cell_range = df_pct[cell_type].max() - df_pct[cell_type].min()
        range_row[cell_type] = cell_range
    stats_data.append(range_row)
    
    # Calculate rMAD (relative median absolute deviation) for each cell type
    rmad_row = {'sample': 'rMAD'}
    for cell_type in df_pct.columns:
        values = df_pct[cell_type].values
        median_val = np.median(values)
        mad = np.median(np.abs(values - median_val))
        rmad = (mad / median_val) * 100 if median_val != 0 else 0
        rmad_row[cell_type] = rmad
    stats_data.append(rmad_row)
    
    # Create final DataFrame
    df_final = pd.DataFrame(stats_data)
    
    # Save to results/Shaun2/
    output_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/FourSample-FourPop-Percentages.csv'
    df_final.to_csv(output_path, index=False)
    
    print(f"Saved: {output_path}")
    print("\nContents:")
    print(df_final.round(3))
    print("\nStatistics included:")
    print("- Percentages for each of the 4 organoid samples")
    print("- Range: Maximum - Minimum percentage for each cell type")
    print("- rMAD: Relative Median Absolute Deviation (variability measure)")

if __name__ == "__main__":
    create_foursample_percentages()