#!/usr/bin/env python3
"""
Simple exploration of GSE131928 files and create a manual research guide
"""

import os
import json
from pathlib import Path
import re

def explore_local_files():
    """Explore what GSE131928 files we have locally"""
    print("🔍 Exploring local GSE131928 files...")
    
    data_dir = Path('/Users/shaunie/Desktop/hu-glioblastoma/data')
    gse_files = []
    
    # Look for GSE131928 related files
    for file_path in data_dir.rglob('*'):
        if 'GSE131928' in file_path.name or 'gse131928' in file_path.name.lower():
            if file_path.is_file():
                try:
                    size = file_path.stat().st_size
                    gse_files.append({
                        'path': str(file_path),
                        'name': file_path.name,
                        'size': size,
                        'type': file_path.suffix
                    })
                except:
                    pass
    
    return gse_files

def read_metadata_file():
    """Read the GSE131928 metadata CSV file"""
    print("📄 Reading GSE131928 metadata file...")
    
    metadata_path = '/Users/shaunie/Desktop/hu-glioblastome/data/GSE131928_metadata.csv'
    
    try:
        with open(metadata_path, 'r') as f:
            content = f.read()
        
        # Extract key information
        lines = content.split('\n')
        
        info = {
            'title': '',
            'summary': '',
            'organism': '',
            'tissue': '',
            'cell_line': '',
            'platform': '',
            'samples': []
        }
        
        for line in lines:
            if line.startswith('title,'):
                info['title'] = line.split(',', 1)[1] if ',' in line else ''
            elif line.startswith('summary,'):
                info['summary'] = line.split(',', 1)[1] if ',' in line else ''
            elif 'organism' in line.lower():
                info['organism'] = line
            elif 'tissue' in line.lower():
                info['tissue'] = line
        
        return info, content
        
    except Exception as e:
        print(f"  Error reading metadata: {e}")
        return None, None

def create_research_guide():
    """Create a manual research guide for finding the original paper"""
    
    guide = {
        "research_strategy": {
            "step_1": "Search GSE131928 on NCBI GEO",
            "step_2": "Look for 'Related Articles' or 'Citation' section",
            "step_3": "Search key terms from our metadata",
            "step_4": "Check for supplementary data files",
            "step_5": "Look for author names and affiliations"
        },
        
        "search_terms": [
            "GSE131928",
            "single cell RNA-seq glioblastoma IDH-wildtype",
            "24131 single cells glioblastoma", 
            "28 patients glioblastoma single cell",
            "adult pediatric glioblastoma single cell",
            "Smartseq2 10X glioblastoma"
        ],
        
        "databases_to_check": [
            "NCBI GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131928)",
            "PubMed (https://pubmed.ncbi.nlm.nih.gov/)",
            "Google Scholar",
            "bioRxiv preprint server",
            "Nature, Cell, Science journals (common for this type of work)"
        ],
        
        "what_to_look_for": [
            "Original cell type classifications",
            "Clustering methodology used",
            "Marker genes for each cell type",
            "Number of clusters identified", 
            "Annotation method (manual vs automated)",
            "Supplementary tables with cell type assignments"
        ],
        
        "our_analysis_summary": {
            "our_cell_types": [
                "Astrocytic", "Mesenchymal", "Neural_Progenitor", 
                "Oligodendrocytic", "Cycling", "Endothelial"
            ],
            "our_methods": [
                "Scanpy pipeline", "Louvain clustering", "Leiden clustering",
                "Manual annotation based on marker genes"
            ],
            "note": "We re-analyzed the raw GSE131928 data with our own pipeline"
        }
    }
    
    return guide

def analyze_file_structure():
    """Analyze what types of files we have"""
    print("📁 Analyzing file structure...")
    
    data_dir = Path('/Users/shaunie/Desktop/hu-glioblastoma/data')
    file_analysis = {
        'h5_files': [],
        'csv_files': [],
        'xlsx_files': [],
        'tsv_files': [],
        'other_files': []
    }
    
    for file_path in data_dir.rglob('*'):
        if file_path.is_file():
            suffix = file_path.suffix.lower()
            name = file_path.name
            
            if suffix == '.h5':
                file_analysis['h5_files'].append(name)
            elif suffix == '.csv':
                file_analysis['csv_files'].append(name)
            elif suffix in ['.xlsx', '.xls']:
                file_analysis['xlsx_files'].append(name)
            elif suffix == '.tsv':
                file_analysis['tsv_files'].append(name)
            elif suffix in ['.txt', '.log', '.md', '.py']:
                continue  # Skip common files
            else:
                file_analysis['other_files'].append(name)
    
    return file_analysis

def main():
    """Main exploration function"""
    print("🧬 GSE131928 Simple Exploration")
    print("=" * 50)
    
    results = {}
    
    # 1. Explore local files
    print("\n1. LOCAL FILES")
    print("-" * 20)
    results['gse_files'] = explore_local_files()
    
    for file_info in results['gse_files']:
        print(f"  📄 {file_info['name']} ({file_info['size']} bytes)")
    
    # 2. Read metadata if available
    print("\n2. METADATA ANALYSIS")  
    print("-" * 20)
    metadata_info, raw_metadata = read_metadata_file()
    results['metadata'] = metadata_info
    
    if metadata_info:
        print(f"  Title: {metadata_info.get('title', 'Not found')}")
        summary = metadata_info.get('summary', 'Not found')
        if len(summary) > 150:
            summary = summary[:150] + "..."
        print(f"  Summary: {summary}")
    else:
        print("  ⚠️  Could not read metadata file")
    
    # 3. Analyze file structure
    print("\n3. FILE STRUCTURE")
    print("-" * 20)
    file_analysis = analyze_file_structure()
    results['file_structure'] = file_analysis
    
    for file_type, files in file_analysis.items():
        if files:
            print(f"  {file_type}: {len(files)} files")
            for file_name in files[:3]:  # Show first 3
                print(f"    - {file_name}")
            if len(files) > 3:
                print(f"    ... and {len(files) - 3} more")
    
    # 4. Create research guide
    print("\n4. RESEARCH GUIDE")
    print("-" * 20)
    guide = create_research_guide()
    results['research_guide'] = guide
    
    print("  Manual research steps:")
    for step, description in guide['research_strategy'].items():
        print(f"    {step}: {description}")
    
    print(f"\n  Key search terms:")
    for term in guide['search_terms'][:4]:
        print(f"    - \"{term}\"")
    
    # 5. Save results
    output_file = '/Users/shaunie/Desktop/hu-glioblastoma/Sandbox/GSE131928_simple_exploration.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print(f"\n💾 Results saved to: {output_file}")
    
    # 6. Create action items
    print("\n" + "=" * 50)
    print("🎯 NEXT STEPS")
    print("=" * 50)
    
    print("\n1. MANUAL RESEARCH (Recommended):")
    print("   → Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131928")
    print("   → Look for 'Related Articles' or 'Citation' links")
    print("   → Search PubMed with key terms from metadata")
    
    print("\n2. WHAT WE KNOW:")
    print("   → GSE131928 has 24,131 single cells from 28 GBM patients")
    print("   → Mixed Smartseq2 (7,930 cells) and 10X (16,201 cells)")  
    print("   → Adult and pediatric IDH-wildtype glioblastomas")
    
    print("\n3. OUR CONTRIBUTION:")
    print("   → We re-analyzed this raw data with our own methods")
    print("   → Created 6 cell type categories with biological relevance")
    print("   → Compared patient data to organoid model")
    
    print("\n✅ Exploration complete!")

if __name__ == "__main__":
    main()