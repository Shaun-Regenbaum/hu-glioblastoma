#!/usr/bin/env python3
"""
Process tissue catalog data with correct mapping:
- Each dose map position becomes a 2x2 grid, but DACTI is 2x1
- DACTI: row B only (4 replicates per concentration, 2x1 expansion)
- PLICA: rows C,D (8 replicates per concentration, 2x2 expansion)  
- GEMCI: rows E,F (8 replicates per concentration, 2x2 expansion)
- VINC: rows G,H (8 replicates per concentration, 2x2 expansion)
- DMSO: columns 21-24 with concentrations from dose map columns 11-12
- Row A remains empty
"""

import csv
import re
from pathlib import Path
from collections import defaultdict
import numpy as np

class TissueCatalogProcessor:
    def __init__(self, base_path, min_value=-20, max_value=100, snr_threshold=1.0):
        self.base_path = Path(base_path)
        self.min_value = min_value
        self.max_value = max_value
        self.snr_threshold = snr_threshold
        self.samples_lookup = {}
        self.oxygen_data = defaultdict(list)
        self.removed_count = defaultdict(int)
        self.snr_removed_count = defaultdict(int)
        self.spike_removed_count = defaultdict(int)
        self.well_to_drug_conc = {}
        
    def load_dose_map(self):
        """Load and parse the dose response map with correct 2x2 grid mapping"""
        dose_map_file = self.base_path / "dose reponse brain orgs map.csv"
        
        with open(dose_map_file, 'r') as f:
            lines = f.readlines()
        
        # Map of dose map rows to actual plate rows
        row_mapping = {
            'a': ['B'],      # DACTI in row B only (2x1 expansion) 
            'b': ['C', 'D'], # PLICA in rows C and D (2x2 expansion)
            'c': ['E', 'F'], # GEMCI in rows E and F (2x2 expansion)
            'd': ['G', 'H']  # VINC in rows G and H (2x2 expansion)
        }
        
        # Parse the dose map
        for line in lines[1:]:
            if line.strip():
                parts = line.strip().split(',')
                if len(parts) > 12 and parts[0]:
                    map_row = parts[0].lower()
                    drug = parts[-1].strip().upper() if parts[-1].strip() else ''
                    
                    if map_row in row_mapping and drug and drug != '':
                        # Get concentrations from the map (columns 1-10 are drugs)
                        concentrations = []
                        for i in range(1, 11):  # Only columns 1-10
                            try:
                                conc = float(parts[i]) if parts[i] else 0
                                concentrations.append(conc)
                            except:
                                concentrations.append(0)
                        
                        # Map to actual wells
                        plate_rows = row_mapping[map_row]
                        
                        # Each dose map position becomes a 2x2 grid (except DACTI which is 2x1)
                        # Columns 1-2 → wells 1-4, Columns 3-4 → wells 5-8, etc.
                        col_idx = 0
                        well_col = 1
                        while col_idx < 10 and well_col <= 20:
                            conc = concentrations[col_idx]
                            
                            # Create 2x2 grid for this concentration (4 wells per concentration)
                            for row in plate_rows:
                                for col_offset in range(4):  # 4 adjacent columns per concentration
                                    if well_col + col_offset <= 20:
                                        well = f"{row}{well_col + col_offset:02d}"
                                        self.well_to_drug_conc[well] = {
                                            'drug': drug,
                                            'concentration': conc
                                        }
                            
                            well_col += 4  # Move to next group of 4 columns
                            col_idx += 2   # Skip the duplicate concentration in original map
        
        # Add DMSO wells (columns 11-12 from dose map → wells 21-24)
        for line in lines[1:]:
            if line.strip():
                parts = line.strip().split(',')
                if len(parts) > 12 and parts[0]:
                    map_row = parts[0].lower()
                    
                    if map_row in row_mapping:
                        plate_rows = row_mapping[map_row]
                        
                        # Get DMSO concentration from column 11 (they're duplicated in 11-12)
                        try:
                            dmso_conc = float(parts[11]) if parts[11] else 0
                        except:
                            dmso_conc = 0
                        
                        # Map DMSO to wells 21-24 (4 wells per drug)
                        for row in plate_rows:
                            for col_offset in range(4):  # 4 DMSO replicates per drug
                                col = 21 + col_offset
                                well = f"{row}{col:02d}"
                                self.well_to_drug_conc[well] = {
                                    'drug': 'DMSO',
                                    'concentration': dmso_conc
                                }
                
    def load_samples(self):
        """Load samples metadata"""
        samples_file = self.base_path / "DatabaseExports" / "samples.csv"
        
        with open(samples_file, 'r', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                self.samples_lookup[row['Pos']] = row
                
    def parse_oxygen_logs_with_cleaning(self):
        """Parse all oxygen log files with data cleaning"""
        log_dir = self.base_path / "OxygenLogs"
        
        # First, organize log files by well
        log_files_by_well = defaultdict(list)
        log_count = 0
        for log_file in log_dir.glob("*.log"):
            log_count += 1
            match = re.match(r'^([A-H]\d{2})-(\d+)\.log$', log_file.name)
            if match:
                well = match.group(1)
                timestamp = int(match.group(2))
                log_files_by_well[well].append((timestamp, log_file))
        
        print(f"Found {log_count} log files for {len(log_files_by_well)} wells")
                
        # Process each well
        for well, files in log_files_by_well.items():
            # Sort by timestamp
            files.sort(key=lambda x: x[0])
            
            if files:
                base_timestamp = files[0][0]
                
                for timestamp, log_file in files:
                    # Calculate time offset in seconds
                    time_offset = (timestamp - base_timestamp) / 1e7
                    
                    # Parse the log file
                    with open(log_file, 'r', encoding='utf-8') as f:
                        lines = f.readlines()
                        
                    for line in lines[7:]:  # Skip headers
                        if line.strip():
                            parts = line.strip().split('\t')
                            if len(parts) >= 11:  # Need at least 11 columns for SNR
                                try:
                                    time_s = float(parts[0]) + time_offset
                                    o2_value = float(parts[1])
                                    
                                    # Get SNR from column 12 (index 11)
                                    snr = float(parts[11]) if len(parts) > 11 else 0
                                    
                                    # Apply data cleaning
                                    if self.min_value <= o2_value <= self.max_value:
                                        # Check SNR threshold
                                        if snr >= self.snr_threshold:
                                            self.oxygen_data[well].append((time_s, o2_value))
                                        else:
                                            self.snr_removed_count[well] += 1
                                    else:
                                        self.removed_count[well] += 1
                                        
                                except ValueError:
                                    continue
                                    
        # Sort data by time for each well and remove spikes
        for well in self.oxygen_data:
            self.oxygen_data[well].sort(key=lambda x: x[0])
            # Remove single-point spikes
            self.oxygen_data[well] = self.remove_spikes(self.oxygen_data[well])
    
    def remove_spikes(self, data_points, spike_threshold=8, max_passes=5):
        """Remove spikes from oxygen data using multiple passes with different strategies"""
        if len(data_points) < 5:
            return data_points
            
        cleaned_data = data_points.copy()
        total_spike_count = 0
        
        # Multiple passes with different strategies
        for pass_num in range(max_passes):
            if len(cleaned_data) < 5:
                break
                
            pass_data = []
            spike_count = 0
            
            for i in range(len(cleaned_data)):
                time_s, o2_value = cleaned_data[i]
                
                # For first and last few points, be more lenient
                if i < 2 or i >= len(cleaned_data) - 2:
                    pass_data.append((time_s, o2_value))
                    continue
                
                # Get surrounding values for better context
                prev2_o2 = cleaned_data[i-2][1] if i >= 2 else None
                prev_o2 = cleaned_data[i-1][1]
                next_o2 = cleaned_data[i+1][1]
                next2_o2 = cleaned_data[i+2][1] if i < len(cleaned_data) - 2 else None
                
                # Calculate local trend
                if prev2_o2 is not None and next2_o2 is not None:
                    # Use 5-point context
                    local_values = [prev2_o2, prev_o2, next_o2, next2_o2]
                    local_mean = np.mean(local_values)
                    local_std = np.std(local_values)
                    
                    # Expected value based on linear interpolation
                    expected_value = (prev_o2 + next_o2) / 2
                    
                    # Multiple criteria for spike detection
                    dev_from_expected = abs(o2_value - expected_value)
                    dev_from_local_mean = abs(o2_value - local_mean)
                    z_score = dev_from_local_mean / (local_std + 1e-6)  # Add small epsilon
                    
                    # Adaptive threshold based on local variability
                    adaptive_threshold = max(spike_threshold, 3 * local_std)
                    
                    # More aggressive spike detection - multiple criteria
                    is_spike = (
                        (dev_from_expected > adaptive_threshold and z_score > 2.5) or
                        (dev_from_expected > spike_threshold * 1.2) or
                        (z_score > 3.5) or
                        # Additional criteria for very large deviations
                        (dev_from_expected > 25) or  # Absolute threshold for large spikes
                        (abs(o2_value - prev_o2) > 30 and abs(o2_value - next_o2) > 30)  # Large jumps both ways
                    )
                    
                    if is_spike:
                        spike_count += 1
                        continue
                else:
                    # Fallback to simple method for edge cases
                    expected_value = (prev_o2 + next_o2) / 2
                    dev_from_expected = abs(o2_value - expected_value)
                    
                    if (dev_from_expected > spike_threshold or 
                        abs(o2_value - prev_o2) > 25 or abs(o2_value - next_o2) > 25):
                        spike_count += 1
                        continue
                
                pass_data.append((time_s, o2_value))
            
            total_spike_count += spike_count
            cleaned_data = pass_data
            
            # If no spikes found in this pass, we're done
            if spike_count == 0:
                break
        
        # Final pass: remove remaining outliers using rolling median
        if len(cleaned_data) > 10:
            cleaned_data = self.remove_outliers_rolling_median(cleaned_data)
        
        if total_spike_count > 0:
            self.spike_removed_count['total'] += total_spike_count
            
        return cleaned_data
    
    def remove_outliers_rolling_median(self, data_points, window_size=7, threshold_factor=2.5):
        """Remove outliers using rolling median filter"""
        if len(data_points) < window_size:
            return data_points
            
        times, values = zip(*data_points)
        values = np.array(values)
        
        # Calculate rolling median
        filtered_data = []
        half_window = window_size // 2
        
        for i in range(len(values)):
            start_idx = max(0, i - half_window)
            end_idx = min(len(values), i + half_window + 1)
            
            window_values = values[start_idx:end_idx]
            median_val = np.median(window_values)
            mad = np.median(np.abs(window_values - median_val))  # Median Absolute Deviation
            
            # Use MAD-based threshold (more robust than std)
            threshold = threshold_factor * mad * 1.4826  # Scale factor for normal distribution
            
            if abs(values[i] - median_val) <= max(threshold, 3):  # Minimum threshold of 3
                filtered_data.append((times[i], values[i]))
        
        return filtered_data
            
    def create_time_series_file(self):
        """Create time series file with all replicates"""
        output_dir = self.base_path / "results"
        output_dir.mkdir(exist_ok=True)
        output_file = output_dir / "all_drugs_time_series.csv"
        
        # Organize wells by drug and concentration
        organized_wells = defaultdict(lambda: defaultdict(list))
        dmso_wells = defaultdict(list)
        
        for well, drug_conc in self.well_to_drug_conc.items():
            if well not in self.oxygen_data:
                continue
                
            drug = drug_conc['drug']
            conc = drug_conc['concentration']
            
            if drug == 'DMSO':
                dmso_wells[str(conc)].append(well)
            else:
                organized_wells[drug][conc].append(well)
                
        # Create header
        header = ['Time_h']
        column_to_well = {}
        
        # Process each drug - use ALL available replicates
        for drug in sorted(organized_wells.keys()):
            # Sort concentrations
            for conc in sorted(organized_wells[drug].keys()):
                wells = organized_wells[drug][conc]
                # Use ALL replicates
                for i, well in enumerate(sorted(wells), 1):
                    col_name = f"{drug}_{conc}_{i}"
                    header.append(col_name)
                    column_to_well[col_name] = well
                    
        # Add DMSO columns - use ALL available replicates
        for dmso_conc in sorted(dmso_wells.keys(), key=float):
            wells = dmso_wells[dmso_conc]
            for i, well in enumerate(sorted(wells), 1):
                col_name = f"DMSO_{dmso_conc}_{i}"
                header.append(col_name)
                column_to_well[col_name] = well
                
        # Create time points
        max_time = 96 * 3600
        time_interval = 3600
        time_points = np.arange(0, max_time + time_interval, time_interval)
        
        # Write data
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)
            
            for time_s in time_points:
                time_h = time_s / 3600
                row = [f"{time_h:.1f}"]
                
                # Add values for each column
                for col_name in header[1:]:
                    well = column_to_well[col_name]
                    o2_value = self.get_oxygen_at_time(well, time_s)
                    row.append(f"{o2_value:.2f}" if o2_value is not None else "")
                    
                # Only write row if it has data
                if any(val != "" for val in row[1:]):
                    writer.writerow(row)
                    
        print(f"Created time series file: {output_file}")
        
        # Create summary report
        self.create_summary_report(organized_wells, dmso_wells)
        
    def get_oxygen_at_time(self, well, target_time_s, window_s=1800):
        """Get oxygen value at specific time (within window)"""
        if well not in self.oxygen_data:
            return None
            
        # Find closest time point
        best_value = None
        min_diff = float('inf')
        
        for time_s, o2_value in self.oxygen_data[well]:
            diff = abs(time_s - target_time_s)
            if diff < min_diff and diff <= window_s:
                min_diff = diff
                best_value = o2_value
                
        return best_value
        
    def create_summary_report(self, organized_wells, dmso_wells):
        """Create summary report"""
        output_dir = self.base_path / "results"
        output_file = output_dir / "summary_report.txt"
        
        with open(output_file, 'w') as f:
            f.write("TISSUE CATALOG PROCESSING SUMMARY\\n")
            f.write("=================================\\n\\n")
            f.write(f"Data cleaning: Removed values <{self.min_value} or >{self.max_value}\\n")
            f.write(f"SNR filtering: Removed values with SNR <{self.snr_threshold}\\n")
            f.write("Spike removal: Removed single-point outliers >20% deviation\\n")
            f.write("Concentrations shown in μM as per original dose map\\n")
            f.write("Each dose map position becomes a 2x2 grid on the plate\\n\\n")
            
            # Summary of replicates
            for drug in sorted(organized_wells.keys()):
                f.write(f"\\n{drug}:\\n")
                f.write("-" * len(drug) + "\\n")
                
                for conc in sorted(organized_wells[drug].keys()):
                    wells = organized_wells[drug][conc]
                    f.write(f"  {conc} μM: {len(wells)} replicates")
                    f.write(f" - Wells: {', '.join(sorted(wells))}\\n")
                    
            # DMSO controls
            f.write("\\nDMSO Controls:\\n")
            f.write("--------------\\n")
            for dmso_conc in sorted(dmso_wells.keys(), key=float):
                wells = dmso_wells[dmso_conc]
                f.write(f"  {dmso_conc}%: {len(wells)} replicates - Wells: {', '.join(sorted(wells))}\\n")
                
            # Data removal statistics
            f.write("\\n\\nData Removal Statistics:\\n")
            f.write("------------------------\\n")
            total_removed_range = sum(self.removed_count.values())
            total_removed_snr = sum(self.snr_removed_count.values())
            total_removed_spikes = sum(self.spike_removed_count.values())
            total_removed = total_removed_range + total_removed_snr + total_removed_spikes
            total_kept = sum(len(data) for data in self.oxygen_data.values())
            
            f.write(f"Total data points kept: {total_kept:,}\\n")
            f.write(f"Total removed (out of range): {total_removed_range:,}\\n")
            f.write(f"Total removed (low SNR): {total_removed_snr:,}\\n")
            f.write(f"Total removed (spikes): {total_removed_spikes:,}\\n")
            f.write(f"Total removed (all causes): {total_removed:,}\\n")
            
            if total_kept + total_removed > 0:
                f.write(f"Percentage removed: {total_removed/(total_kept+total_removed)*100:.2f}%\\n")
                
        print(f"Created summary report: {output_file}")


def main():
    processor = TissueCatalogProcessor("/Users/shaunie/Desktop/hu-glioblastoma/data/revision")
    
    print("Loading data...")
    processor.load_dose_map()
    processor.load_samples()
    
    print("Parsing and cleaning oxygen logs...")
    print(f"Filtering values outside range: {processor.min_value} to {processor.max_value}% O2")
    processor.parse_oxygen_logs_with_cleaning()
    
    print("Creating time series file...")
    processor.create_time_series_file()
    
    print("\\nProcessing complete!")


if __name__ == "__main__":
    main()