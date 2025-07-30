#!/usr/bin/env python3
"""
Create additional CSV summaries and statistics
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

class CSVSummaryGenerator:
    def __init__(self, base_path):
        self.base_path = Path(base_path)
        self.results_dir = self.base_path / "results"
        
    def load_time_series_data(self):
        """Load the time series data"""
        file_path = self.results_dir / "all_drugs_time_series.csv"
        return pd.read_csv(file_path)
    
    def create_endpoint_summary(self):
        """Create summary of endpoint values (8h)"""
        data = self.load_time_series_data()
        
        # Get the 8-hour time point instead of the last one
        target_time = 8.0
        time_diff = (data['Time_h'] - target_time).abs()
        endpoint_time_idx = time_diff.idxmin()
        actual_time = data.iloc[endpoint_time_idx]['Time_h']
        print(f"Using time point {actual_time}h for endpoint analysis (target: {target_time}h)")
        
        results = []
        
        # Process drug columns
        for col in data.columns:
            if col != 'Time_h' and not col.startswith('DMSO'):
                parts = col.split('_')
                if len(parts) >= 3:
                    drug = parts[0]
                    conc = float(parts[1])
                    replicate = int(parts[2])
                    
                    endpoint_value = data.iloc[endpoint_time_idx][col]
                    
                    if not pd.isna(endpoint_value):
                        results.append({
                            'Drug': drug,
                            'Concentration_uM': conc,
                            'Replicate': replicate,
                            'Well': col,  # Will need to map back to actual well
                            'Endpoint_O2_Percent': endpoint_value
                        })
        
        # Process DMSO columns
        for col in data.columns:
            if col.startswith('DMSO'):
                parts = col.split('_')
                if len(parts) >= 3:
                    conc = float(parts[1])
                    replicate = int(parts[2])
                    
                    endpoint_value = data.iloc[endpoint_time_idx][col]
                    
                    if not pd.isna(endpoint_value):
                        results.append({
                            'Drug': 'DMSO',
                            'Concentration_uM': conc,  # Actually percentage for DMSO
                            'Replicate': replicate,
                            'Well': col,
                            'Endpoint_O2_Percent': endpoint_value
                        })
        
        df = pd.DataFrame(results)
        output_file = self.results_dir / "endpoint_values.csv"
        df.to_csv(output_file, index=False)
        
        print(f"Created endpoint values CSV: {output_file}")
        return df
    
    def create_drug_statistics(self):
        """Create statistics summary for each drug and concentration"""
        data = self.load_time_series_data()
        
        results = []
        
        # Process each drug
        drugs = ['DACTI', 'PLICA', 'GEMCI', 'VINC']
        
        for drug in drugs:
            drug_cols = [col for col in data.columns if col.startswith(f"{drug}_")]
            
            # Group by concentration
            conc_groups = defaultdict(list)
            for col in drug_cols:
                parts = col.split('_')
                if len(parts) >= 3:
                    conc = float(parts[1])
                    conc_groups[conc].append(col)
            
            # Calculate statistics for each concentration
            for conc, cols in conc_groups.items():
                # Get all values for this concentration across all time points
                all_values = []
                for col in cols:
                    values = data[col].dropna().values
                    all_values.extend(values)
                
                if all_values:
                    # Calculate endpoint statistics (8-hour time point)
                    target_time = 8.0
                    time_diff = (data['Time_h'] - target_time).abs()
                    endpoint_time_idx = time_diff.idxmin()
                    endpoint_values = []
                    for col in cols:
                        val = data.iloc[endpoint_time_idx][col]
                        if not pd.isna(val):
                            endpoint_values.append(val)
                    
                    results.append({
                        'Drug': drug,
                        'Concentration_uM': conc,
                        'Num_Replicates': len(cols),
                        'Num_Valid_Endpoints': len(endpoint_values),
                        'Mean_All_Values': np.mean(all_values),
                        'Std_All_Values': np.std(all_values),
                        'Min_All_Values': np.min(all_values),
                        'Max_All_Values': np.max(all_values),
                        'Mean_Endpoint': np.mean(endpoint_values) if endpoint_values else np.nan,
                        'Std_Endpoint': np.std(endpoint_values) if len(endpoint_values) > 1 else np.nan,
                        'Min_Endpoint': np.min(endpoint_values) if endpoint_values else np.nan,
                        'Max_Endpoint': np.max(endpoint_values) if endpoint_values else np.nan
                    })
        
        # Process DMSO
        dmso_cols = [col for col in data.columns if col.startswith('DMSO_')]
        conc_groups = defaultdict(list)
        for col in dmso_cols:
            parts = col.split('_')
            if len(parts) >= 3:
                conc = float(parts[1])
                conc_groups[conc].append(col)
        
        for conc, cols in conc_groups.items():
            all_values = []
            for col in cols:
                values = data[col].dropna().values
                all_values.extend(values)
            
            if all_values:
                target_time = 8.0
                time_diff = (data['Time_h'] - target_time).abs()
                endpoint_time_idx = time_diff.idxmin()
                endpoint_values = []
                for col in cols:
                    val = data.iloc[endpoint_time_idx][col]
                    if not pd.isna(val):
                        endpoint_values.append(val)
                
                results.append({
                    'Drug': 'DMSO',
                    'Concentration_uM': conc,  # Actually percentage
                    'Num_Replicates': len(cols),
                    'Num_Valid_Endpoints': len(endpoint_values),
                    'Mean_All_Values': np.mean(all_values),
                    'Std_All_Values': np.std(all_values),
                    'Min_All_Values': np.min(all_values),
                    'Max_All_Values': np.max(all_values),
                    'Mean_Endpoint': np.mean(endpoint_values) if endpoint_values else np.nan,
                    'Std_Endpoint': np.std(endpoint_values) if len(endpoint_values) > 1 else np.nan,
                    'Min_Endpoint': np.min(endpoint_values) if endpoint_values else np.nan,
                    'Max_Endpoint': np.max(endpoint_values) if endpoint_values else np.nan
                })
        
        df = pd.DataFrame(results)
        df = df.sort_values(['Drug', 'Concentration_uM'])
        
        output_file = self.results_dir / "drug_statistics.csv"
        df.to_csv(output_file, index=False)
        
        print(f"Created drug statistics CSV: {output_file}")
        return df
    
    def create_time_course_summary(self):
        """Create summary of key time points"""
        data = self.load_time_series_data()
        
        # Key time points to analyze
        key_times = [0, 6, 12, 24, 48, 72, 96]
        
        results = []
        
        for target_time in key_times:
            # Find closest time point in data
            closest_idx = (data['Time_h'] - target_time).abs().idxmin()
            actual_time = data.iloc[closest_idx]['Time_h']
            
            # Process each column
            for col in data.columns:
                if col != 'Time_h':
                    value = data.iloc[closest_idx][col]
                    
                    if not pd.isna(value):
                        if col.startswith('DMSO'):
                            parts = col.split('_')
                            drug = 'DMSO'
                            conc = float(parts[1])
                            replicate = int(parts[2])
                        else:
                            parts = col.split('_')
                            drug = parts[0]
                            conc = float(parts[1])
                            replicate = int(parts[2])
                        
                        results.append({
                            'Target_Time_h': target_time,
                            'Actual_Time_h': actual_time,
                            'Drug': drug,
                            'Concentration_uM': conc,
                            'Replicate': replicate,
                            'Column': col,
                            'O2_Percent': value
                        })
        
        df = pd.DataFrame(results)
        output_file = self.results_dir / "time_course_summary.csv"
        df.to_csv(output_file, index=False)
        
        print(f"Created time course summary CSV: {output_file}")
        return df
    
    def generate_all_summaries(self):
        """Generate all CSV summaries"""
        print("Generating CSV summaries...")
        
        try:
            self.create_endpoint_summary()
        except Exception as e:
            print(f"Error creating endpoint summary: {e}")
        
        try:
            self.create_drug_statistics()
        except Exception as e:
            print(f"Error creating drug statistics: {e}")
        
        try:
            self.create_time_course_summary()
        except Exception as e:
            print(f"Error creating time course summary: {e}")
        
        print("All CSV summaries created!")


def main():
    generator = CSVSummaryGenerator("/Users/shaunie/Desktop/hu-glioblastoma/data/revision")
    generator.generate_all_summaries()


if __name__ == "__main__":
    main()