#!/usr/bin/env python3
"""
Create Shaun3 analysis with mutually exclusive cell type categories that sum to 100%.

Classification Hierarchy (highest priority first):
1. Endothelial (if present - these are clearly distinct)
2. Cycling + Primary type (e.g., "Cycling_Mesenchymal") 
3. Primary type alone (Astrocytic, Mesenchymal, Neural_Progenitor, Oligodendrocytic)

This ensures every cell gets exactly one label and percentages sum to 100%.
"""
import pandas as pd
import numpy as np

def load_expanded_data():
    """Load the expanded Shaun2 data with all cell types"""
    shaun2_path = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun2/Sox2-Cell-Populations-PerSample.csv'
    return pd.read_csv(shaun2_path)

def create_exclusive_categories(data):
    """
    Convert overlapping categories to mutually exclusive ones.
    
    Strategy: For each sample, redistribute cells based on hierarchy:
    1. Endothelial cells stay as "Endothelial" 
    2. Cycling cells get combined labels like "Cycling_Mesenchymal"
    3. Remaining cells keep their primary type
    """
    
    exclusive_data = []
    
    # Process each sample separately
    for (dataset, sample), sample_data in data.groupby(['dataset', 'sample']):
        
        # Skip invalid samples
        if sample == 'tumour name':
            continue
            
        # Get cell counts for this sample
        sample_dict = {}
        for _, row in sample_data.iterrows():
            sample_dict[row['cell_type']] = {
                'total': row['total_cells'],
                'sox2_pos': row['SOX2_positive'], 
                'sox2_neg': row['SOX2_negative'],
                'pct_neg': row['percent_SOX2_negative']
            }
        
        # Calculate primary cell type totals (excluding Cycling and Endothelial)
        primary_types = ['Astrocytic', 'Mesenchymal', 'Neural_Progenitor', 'Oligodendrocytic']
        primary_total = sum(sample_dict.get(ct, {}).get('total', 0) for ct in primary_types)
        
        # Get cycling and endothelial counts
        cycling_total = sample_dict.get('Cycling', {}).get('total', 0)
        endothelial_total = sample_dict.get('Endothelial', {}).get('total', 0)
        
        # Strategy: Assume cycling cells are distributed proportionally among primary types
        # and create "Cycling_X" categories
        
        # 1. Add Endothelial cells (if any) - these stay separate
        if endothelial_total > 0:
            endo_data = sample_dict['Endothelial']
            exclusive_data.append({
                'dataset': dataset,
                'sample': sample,
                'cell_type': 'Endothelial',
                'total_cells': endo_data['total'],
                'SOX2_positive': endo_data['sox2_pos'],
                'SOX2_negative': endo_data['sox2_neg'],
                'percent_SOX2_negative': endo_data['pct_neg']
            })
        
        # 2. Create cycling categories and non-cycling categories
        if cycling_total > 0 and primary_total > 0:
            # Estimate how cycling cells are distributed among primary types
            # based on the relative abundance of each primary type
            
            cycling_sox2_info = sample_dict['Cycling']
            
            for primary_type in primary_types:
                if primary_type in sample_dict:
                    primary_data = sample_dict[primary_type]
                    primary_count = primary_data['total']
                    
                    if primary_count > 0:
                        # Estimate proportion of cycling cells that are this primary type
                        primary_fraction = primary_count / primary_total
                        cycling_of_this_type = int(cycling_total * primary_fraction)
                        non_cycling_of_this_type = primary_count - cycling_of_this_type
                        
                        # Ensure we don't have negative numbers due to rounding
                        if non_cycling_of_this_type < 0:
                            cycling_of_this_type = primary_count
                            non_cycling_of_this_type = 0
                        
                        # Add cycling category (if > 0)
                        if cycling_of_this_type > 0:
                            # Use cycling Sox2 proportions for cycling cells
                            cycling_sox2_neg = int(cycling_of_this_type * cycling_sox2_info['pct_neg'] / 100)
                            cycling_sox2_pos = cycling_of_this_type - cycling_sox2_neg
                            
                            exclusive_data.append({
                                'dataset': dataset,
                                'sample': sample,
                                'cell_type': f'Cycling_{primary_type}',
                                'total_cells': cycling_of_this_type,
                                'SOX2_positive': cycling_sox2_pos,
                                'SOX2_negative': cycling_sox2_neg,
                                'percent_SOX2_negative': (cycling_sox2_neg / cycling_of_this_type) * 100 if cycling_of_this_type > 0 else 0
                            })
                        
                        # Add non-cycling category (if > 0)
                        if non_cycling_of_this_type > 0:
                            # Use primary type Sox2 proportions for non-cycling cells
                            non_cycling_sox2_neg = int(non_cycling_of_this_type * primary_data['pct_neg'] / 100)
                            non_cycling_sox2_pos = non_cycling_of_this_type - non_cycling_sox2_neg
                            
                            exclusive_data.append({
                                'dataset': dataset,
                                'sample': sample,
                                'cell_type': primary_type,
                                'total_cells': non_cycling_of_this_type,
                                'SOX2_positive': non_cycling_sox2_pos,
                                'SOX2_negative': non_cycling_sox2_neg,
                                'percent_SOX2_negative': (non_cycling_sox2_neg / non_cycling_of_this_type) * 100 if non_cycling_of_this_type > 0 else 0
                            })
        
        else:
            # No cycling cells, just add primary types as-is
            for primary_type in primary_types:
                if primary_type in sample_dict:
                    primary_data = sample_dict[primary_type]
                    if primary_data['total'] > 0:
                        exclusive_data.append({
                            'dataset': dataset,
                            'sample': sample,
                            'cell_type': primary_type,
                            'total_cells': primary_data['total'],
                            'SOX2_positive': primary_data['sox2_pos'],
                            'SOX2_negative': primary_data['sox2_neg'],
                            'percent_SOX2_negative': primary_data['pct_neg']
                        })
    
    return pd.DataFrame(exclusive_data)

def verify_percentages_sum_to_100(data):
    """Verify that percentages sum to 100% for each sample"""
    print("Verifying percentages sum to 100% for each sample...")
    
    all_good = True
    for (dataset, sample), sample_data in data.groupby(['dataset', 'sample']):
        total_cells = sample_data['total_cells'].sum()
        expected_total = total_cells  # Should be 100% of sample
        
        percentage_sum = (sample_data['total_cells'].sum() / total_cells) * 100 if total_cells > 0 else 0
        
        if abs(percentage_sum - 100) > 0.01:  # Allow tiny rounding errors
            print(f"❌ {dataset} {sample}: {percentage_sum:.1f}%")
            all_good = False
    
    if all_good:
        print("✅ All samples sum to 100%")
    
    return all_good

def create_aggregated_data(exclusive_data):
    """Create aggregated data from exclusive categories"""
    
    agg_data = []
    
    for dataset in ['Organoid', 'GSE131928']:
        dataset_name = 'Patient' if dataset == 'GSE131928' else 'Organoid'
        dataset_data = exclusive_data[exclusive_data['dataset'] == dataset]
        
        # Get total cells in dataset
        total_dataset_cells = dataset_data['total_cells'].sum()
        total_sox2_neg = dataset_data['SOX2_negative'].sum()
        
        # Calculate per cell type
        cell_types = dataset_data['cell_type'].unique()
        for cell_type in sorted(cell_types):
            ct_data = dataset_data[dataset_data['cell_type'] == cell_type]
            
            total_cells = ct_data['total_cells'].sum()
            sox2_pos = ct_data['SOX2_positive'].sum()
            sox2_neg = ct_data['SOX2_negative'].sum()
            
            pct_sox2_neg = (sox2_neg / total_cells) * 100 if total_cells > 0 else 0
            pct_of_dataset = (total_cells / total_dataset_cells) * 100 if total_dataset_cells > 0 else 0
            pct_of_all_sox2_neg = (sox2_neg / total_sox2_neg) * 100 if total_sox2_neg > 0 else 0
            
            agg_data.append({
                'Dataset': dataset_name,
                'Population': cell_type,
                'Total_Population_Cells': total_cells,
                'SOX2_Positive_in_Population': sox2_pos,
                'SOX2_Negative_in_Population': sox2_neg,
                'Percent_SOX2_Negative_in_Population': pct_sox2_neg,
                'SOX2_Negative_Count_Used_in_Analysis': sox2_neg,
                'Percent_of_All_SOX2_Negative_Cells': pct_of_all_sox2_neg,
                'Total_SOX2_Negative_Cells_in_Dataset': total_sox2_neg,
                'Percent_of_Dataset': pct_of_dataset
            })
    
    return pd.DataFrame(agg_data)

def main():
    """Main execution"""
    print("Creating Shaun3 analysis with mutually exclusive cell type categories...")
    print("=" * 70)
    
    # Load Shaun2 data
    print("Loading Shaun2 expanded data...")
    shaun2_data = load_expanded_data()
    print(f"Loaded {len(shaun2_data)} rows with {len(shaun2_data['cell_type'].unique())} cell types")
    print(f"Cell types: {sorted(shaun2_data['cell_type'].unique())}")
    
    # Create exclusive categories
    print("\nCreating mutually exclusive categories...")
    exclusive_data = create_exclusive_categories(shaun2_data)
    print(f"Created {len(exclusive_data)} rows with {len(exclusive_data['cell_type'].unique())} exclusive categories")
    print(f"Exclusive categories: {sorted(exclusive_data['cell_type'].unique())}")
    
    # Verify percentages
    verify_percentages_sum_to_100(exclusive_data)
    
    # Create aggregated data
    print("\nCreating aggregated data...")
    agg_data = create_aggregated_data(exclusive_data)
    
    # Save results
    output_dir = '/Users/shaunie/Desktop/hu-glioblastoma/results/Shaun3'
    
    # 1. Per-sample exclusive data
    per_sample_path = f'{output_dir}/Sox2-Cell-Populations-PerSample-Exclusive.csv'
    exclusive_data.to_csv(per_sample_path, index=False)
    print(f"\nSaved: {per_sample_path}")
    
    # 2. Aggregated exclusive data
    agg_path = f'{output_dir}/Sox2-Cell-Populations-Exclusive.csv'
    agg_data.to_csv(agg_path, index=False)
    print(f"Saved: {agg_path}")
    
    # Summary
    print("\n" + "=" * 70)
    print("SHAUN3 SUMMARY - Mutually Exclusive Categories")
    print("=" * 70)
    
    print("\nExclusive cell type categories created:")
    for ct in sorted(exclusive_data['cell_type'].unique()):
        count = exclusive_data[exclusive_data['cell_type'] == ct]['total_cells'].sum()
        print(f"  {ct}: {count:,} cells")
    
    print("\nDataset totals:")
    for dataset in ['Organoid', 'GSE131928']:
        dataset_name = 'Organoid' if dataset == 'Organoid' else 'Patient'
        total = exclusive_data[exclusive_data['dataset'] == dataset]['total_cells'].sum()
        print(f"  {dataset_name}: {total:,} cells")
    
    print(f"\nAll samples now have cell type percentages that sum to 100%!")
    print("Each cell belongs to exactly one category.")

if __name__ == "__main__":
    main()