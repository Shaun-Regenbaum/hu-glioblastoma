#!/usr/bin/env python3
"""
Explore the original GSE131928 paper and dataset to understand:
1. What cell types were originally identified
2. How they classified cells
3. What methodology they used
4. How our re-analysis compares to their original findings
"""

import requests
import pandas as pd
import numpy as np
from pathlib import Path
import re
import json
from urllib.parse import quote
import time

def search_ncbi_gse131928():
    """Search NCBI GEO for GSE131928 information"""
    print("🔍 Searching NCBI GEO for GSE131928...")
    
    # GEO API endpoint
    base_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    params = {
        'acc': 'GSE131928',
        'targ': 'self',
        'form': 'text',
        'view': 'full'
    }
    
    try:
        response = requests.get(base_url, params=params, timeout=10)
        if response.status_code == 200:
            print("✅ Retrieved GSE131928 metadata from NCBI")
            return response.text
        else:
            print(f"❌ Failed to retrieve GSE131928 data: {response.status_code}")
            return None
    except Exception as e:
        print(f"❌ Error accessing NCBI: {e}")
        return None

def search_pubmed_for_paper():
    """Search PubMed for the original GSE131928 paper"""
    print("🔍 Searching PubMed for GSE131928 associated paper...")
    
    # PubMed search API
    search_terms = [
        "GSE131928",
        "single cell RNA-seq glioblastoma IDH-wildtype",
        "24131 single cells glioblastoma",
        "28 patients glioblastoma single cell"
    ]
    
    results = {}
    
    for term in search_terms:
        try:
            # Use NCBI E-utilities API
            search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            params = {
                'db': 'pubmed',
                'term': quote(term),
                'retmode': 'json',
                'retmax': 5
            }
            
            response = requests.get(search_url, params=params, timeout=10)
            if response.status_code == 200:
                data = response.json()
                if 'esearchresult' in data and 'idlist' in data['esearchresult']:
                    results[term] = data['esearchresult']['idlist']
                    print(f"  Found {len(results[term])} papers for: {term}")
            
            time.sleep(0.5)  # Be nice to NCBI servers
            
        except Exception as e:
            print(f"  Error searching for '{term}': {e}")
    
    return results

def get_paper_details(pubmed_ids):
    """Get detailed information about papers from PubMed IDs"""
    print("📄 Retrieving paper details...")
    
    papers = []
    
    for pmid in pubmed_ids[:3]:  # Limit to first 3 papers
        try:
            # Use E-fetch to get paper details
            fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            params = {
                'db': 'pubmed',
                'id': pmid,
                'retmode': 'xml'
            }
            
            response = requests.get(fetch_url, params=params, timeout=10)
            if response.status_code == 200:
                # Extract basic info (simplified XML parsing)
                title_match = re.search(r'<ArticleTitle>(.*?)</ArticleTitle>', response.text, re.DOTALL)
                abstract_match = re.search(r'<AbstractText>(.*?)</AbstractText>', response.text, re.DOTALL)
                
                paper_info = {
                    'pmid': pmid,
                    'title': title_match.group(1) if title_match else 'Title not found',
                    'abstract': abstract_match.group(1) if abstract_match else 'Abstract not found'
                }
                papers.append(paper_info)
                print(f"  Retrieved: PMID {pmid}")
            
            time.sleep(0.5)
            
        except Exception as e:
            print(f"  Error fetching PMID {pmid}: {e}")
    
    return papers

def analyze_geo_metadata(geo_text):
    """Analyze the GEO metadata text for cell type information"""
    print("📊 Analyzing GEO metadata for cell type information...")
    
    if not geo_text:
        return {}
    
    analysis = {
        'title': '',
        'summary': '',
        'cell_types_mentioned': [],
        'methodology': [],
        'sample_info': []
    }
    
    lines = geo_text.split('\n')
    
    for line in lines:
        line = line.strip()
        
        # Extract title
        if line.startswith('!Series_title'):
            analysis['title'] = line.split('\t')[1] if '\t' in line else line
        
        # Extract summary
        elif line.startswith('!Series_summary'):
            analysis['summary'] = line.split('\t')[1] if '\t' in line else line
        
        # Look for cell type mentions
        cell_type_keywords = [
            'astrocyte', 'astrocytic', 'mesenchymal', 'oligodendrocyte', 'oligodendrocytic',
            'neural progenitor', 'progenitor', 'stem cell', 'endothelial', 'immune',
            'microglia', 'macrophage', 'cycling', 'proliferating', 'cell type', 'cluster'
        ]
        
        for keyword in cell_type_keywords:
            if keyword.lower() in line.lower():
                analysis['cell_types_mentioned'].append({
                    'keyword': keyword,
                    'context': line
                })
        
        # Look for methodology mentions
        method_keywords = [
            'clustering', 'louvain', 'leiden', 'seurat', 'scanpy', 'monocle',
            'UMAP', 'tSNE', 'PCA', 'marker genes'
        ]
        
        for keyword in method_keywords:
            if keyword.lower() in line.lower():
                analysis['methodology'].append({
                    'method': keyword,
                    'context': line
                })
    
    return analysis

def search_supplementary_data():
    """Look for supplementary data files that might contain cell type annotations"""
    print("📁 Looking for supplementary data files...")
    
    # Check if we have any supplementary files locally
    data_dir = Path('/Users/shaunie/Desktop/hu-glioblastoma/data')
    
    supplementary_files = []
    
    # Look for files that might contain cell type annotations
    patterns = [
        '*cell*type*',
        '*annotation*',
        '*cluster*',
        '*metadata*',
        '*GSE131928*'
    ]
    
    for pattern in patterns:
        for file_path in data_dir.rglob(pattern):
            if file_path.is_file():
                supplementary_files.append({
                    'path': str(file_path),
                    'name': file_path.name,
                    'size': file_path.stat().st_size
                })
    
    return supplementary_files

def create_comparison_summary():
    """Create a comparison between our analysis and what we find about the original"""
    print("📋 Creating comparison summary...")
    
    our_cell_types = [
        'Astrocytic',
        'Mesenchymal', 
        'Neural_Progenitor',
        'Oligodendrocytic',
        'Cycling',
        'Endothelial'
    ]
    
    our_methods = [
        'Louvain clustering',
        'Leiden clustering', 
        'Manual annotation based on marker genes',
        'UMAP visualization',
        'scanpy pipeline'
    ]
    
    summary = {
        'our_analysis': {
            'cell_types': our_cell_types,
            'methods': our_methods,
            'total_cells': 'Organoids: ~22,807 cells + Patients: ~12,079 cells',
            'samples': 'Organoid: 4 samples, Patient: ~45 tumor samples'
        },
        'findings': {
            'note': 'To be filled with search results'
        }
    }
    
    return summary

def main():
    """Main exploration function"""
    print("🧬 GSE131928 Original Paper Exploration")
    print("=" * 50)
    
    results = {
        'geo_metadata': None,
        'pubmed_search': {},
        'paper_details': [],
        'metadata_analysis': {},
        'supplementary_files': [],
        'comparison_summary': {}
    }
    
    # 1. Get GEO metadata
    print("\n1. RETRIEVING GEO METADATA")
    print("-" * 30)
    results['geo_metadata'] = search_ncbi_gse131928()
    
    # 2. Search PubMed for associated papers
    print("\n2. SEARCHING PUBMED")
    print("-" * 30)
    results['pubmed_search'] = search_pubmed_for_paper()
    
    # 3. Get details for top papers
    if results['pubmed_search']:
        print("\n3. GETTING PAPER DETAILS")
        print("-" * 30)
        # Get all unique PMIDs
        all_pmids = set()
        for term, pmids in results['pubmed_search'].items():
            all_pmids.update(pmids)
        
        if all_pmids:
            results['paper_details'] = get_paper_details(list(all_pmids))
    
    # 4. Analyze GEO metadata
    print("\n4. ANALYZING METADATA")
    print("-" * 30)
    if results['geo_metadata']:
        results['metadata_analysis'] = analyze_geo_metadata(results['geo_metadata'])
    
    # 5. Look for supplementary files
    print("\n5. SEARCHING SUPPLEMENTARY FILES")
    print("-" * 30)
    results['supplementary_files'] = search_supplementary_data()
    
    # 6. Create comparison summary
    print("\n6. CREATING COMPARISON SUMMARY")
    print("-" * 30)
    results['comparison_summary'] = create_comparison_summary()
    
    # Save results
    output_file = '/Users/shaunie/Desktop/hu-glioblastoma/Sandbox/GSE131928_exploration_results.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print(f"\n💾 Results saved to: {output_file}")
    
    # Print summary
    print("\n" + "=" * 50)
    print("🔍 EXPLORATION SUMMARY")
    print("=" * 50)
    
    if results['metadata_analysis']:
        print(f"\n📊 Dataset Title:")
        print(f"  {results['metadata_analysis'].get('title', 'Not found')}")
        
        print(f"\n📝 Summary:")
        summary = results['metadata_analysis'].get('summary', 'Not found')
        if len(summary) > 200:
            summary = summary[:200] + "..."
        print(f"  {summary}")
        
        print(f"\n🔬 Cell Types Mentioned: {len(results['metadata_analysis'].get('cell_types_mentioned', []))}")
        for mention in results['metadata_analysis'].get('cell_types_mentioned', [])[:5]:
            print(f"  - {mention['keyword']}")
        
        print(f"\n⚙️ Methods Mentioned: {len(results['metadata_analysis'].get('methodology', []))}")
        for method in results['metadata_analysis'].get('methodology', [])[:5]:
            print(f"  - {method['method']}")
    
    if results['paper_details']:
        print(f"\n📄 Found {len(results['paper_details'])} Associated Papers:")
        for paper in results['paper_details'][:3]:
            title = paper['title']
            if len(title) > 100:
                title = title[:100] + "..."
            print(f"  - PMID {paper['pmid']}: {title}")
    
    if results['supplementary_files']:
        print(f"\n📁 Found {len(results['supplementary_files'])} Supplementary Files:")
        for file_info in results['supplementary_files'][:10]:
            print(f"  - {file_info['name']} ({file_info['size']} bytes)")
    
    print(f"\n✅ Exploration complete! Check {output_file} for full results.")

if __name__ == "__main__":
    main()