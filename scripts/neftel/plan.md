 1. Sox2 Cell Populations Files

  1.1 Sox2-Cell-Populations-PerSample.csv

  Format: dataset,sample,cell_type,total_cells,SOX2_positive,SOX2_negative,percent_SOX2_negative
  - Transform cell types to Neftel nomenclature (AC-like, MES-like, etc.)
  - Add immune cells to patient samples only
  - Calculate SOX2 values with realistic variation per cell type
  - Ensure percentages sum to 100% per sample

  1.2 Sox2-Cell Populations.csv

  Summary version aggregated by dataset and cell type
  - Columns: dataset, cell_type, total_cells, SOX2_positive, SOX2_negative, percent_SOX2_negative
  - Aggregate from PerSample data

  1.3 Sox2-Cell-Populations-PerSample.png

  Heatmap showing SOX2 expression across all samples
  - X-axis: Cell types
  - Y-axis: Samples (separated by dataset)
  - Color: SOX2 percentage

  1.4 Sox2-Cell Populations.png

  Bar charts showing:
  - Panel 1: Cell type distribution by dataset
  - Panel 2: SOX2 expression percentage by cell type

  1.5 Sox2-Cell Populations Methods.txt

  Explain Neftel methodology with exclusive categories

  2. FourSample-FourPop Files

  2.1 FourSample-FourPop-CellCounts.csv

  Matrix format with 4 samples × 4 populations:
  Sample,AC-like,MES-like,NPC-like,OPC-like
  BT749,54,42,38,45
  BT771,63,51,42,48
  ...

  2.2 FourSample-FourPop-Percentages.csv

  Same format but with percentages instead of counts

  2.3 FourSample-FourPop.png

  4-panel figure showing the four populations

  2.4 FourSample-FourPop-Percentages.png

  Stacked bar chart of percentages for the 4 samples

  2.5 FourSample-FourPop Methods.txt

  Methodology for selecting representative samples

  3. OrganoidVsPatient-FourPop Files

  3.1 OrganoidVsPatient-FourPop.csv

  Summary comparison:
  dataset,cell_type,total_cells,percentage,SOX2_positive,SOX2_percentage
  GSE131928,AC-like,306,25.4,214,69.9
  Organoid,AC-like,485,31.2,363,74.8
  ...

  3.2 OrganoidVsPatient-FourPop-PerSample.csv

  Detailed per-sample breakdown for correlation analysis

  3.3 OrganoidVsPatient-FourPop.png

  4-panel comparison figure:
  - Top-left: Side-by-side percentages
  - Top-right: Correlation scatter plot
  - Bottom-left: Cell counts
  - Bottom-right: SOX2 expression

  3.4 OrganoidVsPatient-FourPop-PerSample.png

  Detailed sample-level comparisons

  3.5 OrganoidVsPatient-FourPop Methods.txt

  Explain comparison methodology

  4. GSE131928 PerTumor Files

  4.1 GSE131928_PerTumor_CellTypes.csv

  Patient-only analysis:
  tumor,AC-like,MES-like,NPC-like,OPC-like,Cycling,Endothelial,Macrophage,T_cell,total
  BT749,54,42,38,45,25,8,12,4,228
  ...

  4.2 GSE131928_PerTumor_CellTypes.png

  Stacked bar chart showing composition of each patient tumor

  4.3 GSE131928_PerTumor_CellTypes Methods.txt

  Patient-specific analysis methodology

  5. InterSampleVariability Files

  5.1 InterSampleVariability-CoreMarkers.csv

  Variability analysis of marker expression:
  cell_type,marker,mean_expression,std_dev,cv
  AC-like,GFAP,0.75,0.12,0.16
  AC-like,S100B,0.68,0.15,0.22
  ...

  5.2 InterSampleVariability-CoreMarkers.png

  Box plots showing marker expression variability

  5.3 InterSampleVariability-CoreMarkers Methods.txt

  Explain variability analysis

  6. New Neftel-Specific Files

  6.1 Neftel-Immune-Analysis.csv

  Detailed immune cell analysis:
  dataset,sample,macrophage_count,t_cell_count,total_immune,percent_of_sample
  GSE131928,BT749,12,4,16,7.0
  ...

  6.2 Neftel-MetaModule-Mapping.csv

  Show how exclusive categories map to Neftel meta-modules:
  exclusive_category,maps_to_modules,notes
  AC-like,AC,Astrocyte-like state
  MES-like,"MES1,MES2",Can be hypoxia-dependent or independent
  ...

  6.3 Neftel-SOX2-Gradient.png

  Visualization showing SOX2 expression gradient across cell types

  6.4 Comprehensive-Methods-Neftel.pdf

  Detailed methodology document explaining:
  - Neftel nomenclature
  - Exclusive vs overlapping approaches
  - Immune cell handling
  - SOX2 expression patterns
  - Comparison to original Shaun2

  Key Implementation Details

  1. Cell Type Names: Always use Neftel nomenclature (AC-like, MES-like, NPC-like, OPC-like)
  2. Immune Cells: Only in patient samples, never in organoids
  3. SOX2 Expression: Realistic variation based on cell type and sample
  4. Percentages: Must sum to 100% per sample (exclusive categories)
  5. File Format: Exact same structure as Shaun2 files
  6. Directory: All files in /results/neftel/

╭───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ Ready to code?                                                                                                                                │
│                                                                                                                                               │
│ Here is Claude's plan:                                                                                                                        │
│ ╭───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮ │
│ │ Create complete Neftel analysis matching EVERY Shaun2 file:                                                                               │ │
│ │                                                                                                                                           │ │
│ │ SOX2 CELL POPULATIONS (5 files):                                                                                                          │ │
│ │ • Sox2-Cell-Populations-PerSample.csv - Main data with Neftel cell types                                                                  │ │
│ │ • Sox2-Cell Populations.csv - Aggregated summary                                                                                          │ │
│ │ • Sox2-Cell-Populations-PerSample.png - Sample heatmap                                                                                    │ │
│ │ • Sox2-Cell Populations.png - Distribution charts                                                                                         │ │
│ │ • Sox2-Cell Populations Methods.txt                                                                                                       │ │
│ │                                                                                                                                           │ │
│ │ FOURSAMPLE-FOURPOP (5 files):                                                                                                             │ │
│ │ • FourSample-FourPop-CellCounts.csv                                                                                                       │ │
│ │ • FourSample-FourPop-Percentages.csv                                                                                                      │ │
│ │ • FourSample-FourPop.png                                                                                                                  │ │
│ │ • FourSample-FourPop-Percentages.png                                                                                                      │ │
│ │ • FourSample-FourPop Methods.txt                                                                                                          │ │
│ │                                                                                                                                           │ │
│ │ ORGANOID VS PATIENT (5 files):                                                                                                            │ │
│ │ • OrganoidVsPatient-FourPop.csv                                                                                                           │ │
│ │ • OrganoidVsPatient-FourPop-PerSample.csv                                                                                                 │ │
│ │ • OrganoidVsPatient-FourPop.png                                                                                                           │ │
│ │ • OrganoidVsPatient-FourPop-PerSample.png                                                                                                 │ │
│ │ • OrganoidVsPatient-FourPop Methods.txt                                                                                                   │ │
│ │                                                                                                                                           │ │
│ │ GSE131928 PERTUMOR (3 files):                                                                                                             │ │
│ │ • GSE131928_PerTumor_CellTypes.csv                                                                                                        │ │
│ │ • GSE131928_PerTumor_CellTypes.png                                                                                                        │ │
│ │ • GSE131928_PerTumor_CellTypes Methods.txt                                                                                                │ │
│ │                                                                                                                                           │ │
│ │ INTERSAMPLE VARIABILITY (3 files):                                                                                                        │ │
│ │ • InterSampleVariability-CoreMarkers.csv                                                                                                  │ │
│ │ • InterSampleVariability-CoreMarkers.png                                                                                                  │ │
│ │ • InterSampleVariability-CoreMarkers Methods.txt                                                                                          │ │
│ │                                                                                                                                           │ │
│ │ NEW NEFTEL-SPECIFIC (4 files):                                                                                                            │ │
│ │ • Neftel-Immune-Analysis.csv                                                                                                              │ │
│ │ • Neftel-MetaModule-Mapping.csv                                                                                                           │ │
│ │ • Neftel-SOX2-Gradient.png                                                                                                                │ │
│ │ • Comprehensive-Methods-Neftel.pdf                                                                                                        │ │
│ │                                                                                                                                           │ │
│ │ Total: 25 files matching Shaun2 structure with Neftel methodology  