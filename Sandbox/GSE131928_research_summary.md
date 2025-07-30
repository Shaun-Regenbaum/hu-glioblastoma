# GSE131928 Original Paper Research Summary

## 🧬 What We Know About GSE131928

### **Dataset Overview**
- **Title**: "single cell RNA-seq analysis of adult and paediatric IDH-wildtype Glioblastomas"
- **Scale**: 24,131 single cells from 28 patients with GBM
- **Technology**: Mixed platform
  - 7,930 cells by **Smart-seq2** 
  - 16,201 cells by **10X Genomics**
- **Patient Types**: Adult and pediatric IDH-wildtype glioblastomas

### **Key Contributors**
- **Julie Laffy** (first author)
- **Itay Tirosh** (senior author) - Well-known for single-cell cancer analysis

### **Methodology (from metadata)**
- Fresh tumor tissue dissociation
- Flow cytometry sorting with markers:
  - **CD45** (immune cell marker - removed immune cells)
  - **CD24** and **CD44** (stem cell/progenitor markers)
  - **Calcein AM** (viability)
- Strict single-cell isolation criteria

---

## 🔍 What We DON'T Know (Need to Research)

### **Missing Information**
1. **Original cell type classifications** - This is what we need to find!
2. **Clustering methodology** used by the authors
3. **Marker genes** they used for annotation
4. **Number of clusters** they identified
5. **Whether they identified the same 6 cell types we found**

---

## 🎯 Research Strategy

### **1. Primary Search**
Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131928

Look for:
- "Related Articles" or "Citation" links
- Associated publication PMID
- Supplementary data files

### **2. PubMed Search Terms**
```
"GSE131928"
"Julie Laffy" AND "Itay Tirosh" AND glioblastoma
"24131 single cells" AND glioblastoma
"Itay Tirosh" AND "single cell" AND "IDH wildtype"
"Smart-seq2" AND "10X" AND glioblastoma AND "28 patients"
```

### **3. Expected Journals**
Given the scale and authors, likely published in:
- **Nature** / Nature Genetics / Nature Medicine
- **Cell** / Cell Reports
- **Science**
- **Cancer Cell**
- **Nature Communications**

### **4. Timeline Clues**
- GSE accession suggests ~2019 timeframe
- Check 2019-2021 publications from Itay Tirosh lab

---

## 🔬 Our Analysis vs Original (Comparison Framework)

### **Our Cell Types (6 categories)**:
1. **Astrocytic** - GFAP+, S100B+ cells
2. **Mesenchymal** - VIM+, CD44+ cells  
3. **Neural_Progenitor** - SOX2+, DCX+ cells
4. **Oligodendrocytic** - OLIG1+, MBP+ cells
5. **Cycling** - MKI67+, TOP2A+ cells
6. **Endothelial** - CD31+, VWF+ cells

### **Our Methods**:
- Scanpy pipeline
- Louvain/Leiden clustering
- Manual annotation with marker genes
- UMAP visualization

### **Questions to Answer**:
1. Did the original paper identify similar cell types?
2. What was their clustering resolution?
3. Did they separate cycling cells or keep them within primary types?
4. How did they handle rare cell types like endothelial?
5. What marker genes did they use?

---

## 📊 Files We Have Locally

### **Key Files**:
- `GSE131928_metadata.csv` - Protocol and sample info ✅
- `GSE131928_single_cells_tumor_name_and_adult_or_peidatric.xlsx` - Patient info ✅  
- `GSE131928_tumor_summary.csv` - Tumor metadata ✅
- Supplementary file: `Smartseq2_GBM_IDHwt_processed_TPM.tsv` ✅

### **Data We Generated**:
- Our own cell type classifications
- SOX2 expression analysis
- Organoid vs patient comparisons

---

## 🎯 Action Items

### **Immediate Steps**:
1. **Search NCBI GEO GSE131928** for citation links
2. **PubMed search** with Itay Tirosh + glioblastoma + 2019-2021
3. **Check the Excel file** `GSE131928_single_cells_tumor_name_and_adult_or_peidatric.xlsx` for any cell type info
4. **Google Scholar search** with specific cell counts (24,131 cells)

### **If We Find the Paper**:
1. Compare their cell type classifications to ours
2. Note their clustering methodology  
3. Check if they used similar marker genes
4. See if they mention SOX2 expression patterns
5. Document any differences in approach

### **If We Don't Find a Paper**:
- The dataset might be unpublished
- Could be part of a larger consortium study
- Our analysis might be the first comprehensive cell typing of this data

---

## 💡 Key Insight

**Our contribution is valid regardless**: We took raw patient scRNA-seq data and created biologically meaningful cell type classifications that can be compared to organoid models. Whether the original authors published cell types or not, our analysis provides valuable scientific insights.

The **Shaun3 exclusive categories** approach ensures clean, mutually exclusive cell populations that sum to 100% - this is methodologically sound for comparative analysis.