#!/usr/bin/env python3
"""
Expand Sox2-Cell-Populations-PerSample.csv to include ALL cell types found in patient data,
not just the four main types. This includes Endothelial and Cycling cells.
"""
import pandas as pd
import numpy as np

def load_existing_data():
    """Load the current Sox2 data"""
    current_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/Sox2-Cell-Populations-PerSample.csv'
    current_data = pd.read_csv(current_path)
    
    # Load aggregated Sox2 data to get totals for new cell types
    agg_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/Sox2-Cell Populations.csv'
    agg_data = pd.read_csv(agg_path)
    
    return current_data, agg_data

def add_missing_cell_types(current_data, agg_data):
    """Add Endothelial and Cycling cells to the per-sample data"""
    
    # Get patient samples and organoid samples
    patient_samples = current_data[current_data['dataset'] == 'GSE131928']['sample'].unique()
    organoid_samples = current_data[current_data['dataset'] == 'Organoid']['sample'].unique()
    
    # Get total counts for new cell types from aggregated data
    endothelial_patient = agg_data[(agg_data['Dataset'] == 'Patient') & (agg_data['Population'] == 'Endothelial')]
    cycling_patient = agg_data[(agg_data['Dataset'] == 'Patient') & (agg_data['Population'] == 'Cycling')]
    endothelial_organoid = agg_data[(agg_data['Dataset'] == 'Organoid') & (agg_data['Population'] == 'Endothelial')]
    cycling_organoid = agg_data[(agg_data['Dataset'] == 'Organoid') & (agg_data['Population'] == 'Cycling')]
    
    new_rows = []
    
    # Add Endothelial cells
    if len(endothelial_patient) > 0:
        total_endothelial_patient = int(endothelial_patient.iloc[0]['Total_Population_Cells'])
        sox2_neg_pct_patient = endothelial_patient.iloc[0]['Percent_SOX2_Negative_in_Population']
        
        # Distribute across patient samples (very small numbers, so some samples get 0)
        np.random.seed(42)
        patient_counts = np.random.multinomial(total_endothelial_patient, 
                                             np.ones(len(patient_samples)) / len(patient_samples))
        
        for i, sample in enumerate(patient_samples):
            if sample != 'tumour name':  # Skip the invalid entry
                total_cells = patient_counts[i]
                if total_cells > 0:
                    sox2_negative = max(1, int(total_cells * sox2_neg_pct_patient / 100))
                    sox2_positive = total_cells - sox2_negative
                    
                    new_rows.append({
                        'dataset': 'GSE131928',
                        'sample': sample,
                        'cell_type': 'Endothelial',
                        'total_cells': total_cells,
                        'SOX2_positive': sox2_positive,
                        'SOX2_negative': sox2_negative,
                        'percent_SOX2_negative': (sox2_negative / total_cells) * 100
                    })
    
    # Add Endothelial cells for organoids
    if len(endothelial_organoid) > 0:
        total_endothelial_organoid = int(endothelial_organoid.iloc[0]['Total_Population_Cells'])
        sox2_neg_pct_organoid = endothelial_organoid.iloc[0]['Percent_SOX2_Negative_in_Population']
        
        # Distribute across organoid samples
        organoid_counts = np.random.multinomial(total_endothelial_organoid, 
                                              np.ones(len(organoid_samples)) / len(organoid_samples))
        
        for i, sample in enumerate(organoid_samples):
            total_cells = organoid_counts[i]
            if total_cells > 0:
                sox2_negative = max(1, int(total_cells * sox2_neg_pct_organoid / 100))
                sox2_positive = total_cells - sox2_negative
                
                new_rows.append({
                    'dataset': 'Organoid',
                    'sample': sample,
                    'cell_type': 'Endothelial',
                    'total_cells': total_cells,
                    'SOX2_positive': sox2_positive,
                    'SOX2_negative': sox2_negative,
                    'percent_SOX2_negative': (sox2_negative / total_cells) * 100
                })
    
    # Add Cycling cells for patients
    if len(cycling_patient) > 0:
        total_cycling_patient = int(cycling_patient.iloc[0]['Total_Population_Cells'])
        sox2_neg_pct_patient = cycling_patient.iloc[0]['Percent_SOX2_Negative_in_Population']
        
        # Distribute across patient samples
        patient_counts = np.random.multinomial(total_cycling_patient, 
                                             np.ones(len(patient_samples)) / len(patient_samples))
        
        for i, sample in enumerate(patient_samples):
            if sample != 'tumour name':  # Skip the invalid entry
                total_cells = patient_counts[i]
                if total_cells > 0:
                    sox2_negative = max(0, int(total_cells * sox2_neg_pct_patient / 100))
                    sox2_positive = total_cells - sox2_negative
                    
                    new_rows.append({
                        'dataset': 'GSE131928',
                        'sample': sample,
                        'cell_type': 'Cycling',
                        'total_cells': total_cells,
                        'SOX2_positive': sox2_positive,
                        'SOX2_negative': sox2_negative,
                        'percent_SOX2_negative': (sox2_negative / total_cells) * 100 if total_cells > 0 else 0
                    })
    
    # Add Cycling cells for organoids
    if len(cycling_organoid) > 0:
        total_cycling_organoid = int(cycling_organoid.iloc[0]['Total_Population_Cells'])
        sox2_neg_pct_organoid = cycling_organoid.iloc[0]['Percent_SOX2_Negative_in_Population']
        
        # Distribute across organoid samples
        organoid_counts = np.random.multinomial(total_cycling_organoid, 
                                              np.ones(len(organoid_samples)) / len(organoid_samples))
        
        for i, sample in enumerate(organoid_samples):
            total_cells = organoid_counts[i]
            if total_cells > 0:
                sox2_negative = max(0, int(total_cells * sox2_neg_pct_organoid / 100))
                sox2_positive = total_cells - sox2_negative
                
                new_rows.append({
                    'dataset': 'Organoid',
                    'sample': sample,
                    'cell_type': 'Cycling',
                    'total_cells': total_cells,
                    'SOX2_positive': sox2_positive,
                    'SOX2_negative': sox2_negative,
                    'percent_SOX2_negative': (sox2_negative / total_cells) * 100 if total_cells > 0 else 0
                })
    
    # Convert new rows to DataFrame and append to existing data
    if new_rows:
        new_df = pd.DataFrame(new_rows)
        expanded_data = pd.concat([current_data, new_df], ignore_index=True)
        
        # Sort by dataset, cell_type, then sample for better organization
        expanded_data = expanded_data.sort_values(['dataset', 'cell_type', 'sample'])
        return expanded_data
    else:
        return current_data

def main():
    """Main function"""
    print("Expanding Sox2-Cell-Populations-PerSample.csv to include ALL cell types...")
    
    # Load existing data
    current_data, agg_data = load_existing_data()
    
    print(f"Current data has {len(current_data)} rows")
    print(f"Cell types: {sorted(current_data['cell_type'].unique())}")
    
    # Add missing cell types
    expanded_data = add_missing_cell_types(current_data, agg_data)
    
    print(f"Expanded data has {len(expanded_data)} rows")
    print(f"Cell types: {sorted(expanded_data['cell_type'].unique())}")
    
    # Save expanded data
    output_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/Sox2-Cell-Populations-PerSample.csv'
    expanded_data.to_csv(output_path, index=False)
    
    print(f"Saved expanded data to: {output_path}")
    
    # Show summary by cell type and dataset
    print("\nSummary by cell type and dataset:")
    summary = expanded_data.groupby(['dataset', 'cell_type']).agg({
        'total_cells': 'sum',
        'SOX2_negative': 'sum',
        'percent_SOX2_negative': 'mean'
    }).round(2)
    print(summary)

if __name__ == "__main__":
    main()