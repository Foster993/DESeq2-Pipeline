Drought and salinity are major environmental factors to crop growth affecting 10% of arable land globally (Bartels & Sunkar, 2005). Finding the specific gene expression responses that overlap between salt and drought could allow scientists to engineer plants with more plasticity, therefore increasing crop yield and allowing plants to thrive under extreme conditions. Unfortunately, minimal relevant research using RNA-seq data has been done comparing the DEGs in roots exposed to salt and roots exposed to drought. In this study, I am aiming to find the shared and stress-specific physiological responses by identifying the DEGs in salt and drought-exposed roots.  My research will use DESeq2, a powerful algorithm, because it accounts for the overdispersion of biological data and finds the DEGs while also normalizing the data to run plots. DESeq2 will be run on my MacBook Air on R Studio version 2025.5.1.513. For the NCBI GEO datasets, four are drought-exposed, while another 4 are drought-control, three are salt-exposed, while another three are salt-control. Then I will format the data into a count matrix, and run DESeq2 for the salt, drought, and control datasets, and finally compare my findings using plots to identify the overlapping DEG’s. By that time, I shall be able to find overlaps using GO analysis, Expression analysis, KEGG analysis, and clustering analysis. By running these bioinformatic pipelines, I expect to identify new potential candidate genes for future research.



COMPLETE Materials List: *
Hardware
MacBook Air M1 (2020)
 8GB RAM running Sequoia 15.6.1
245.11 GB of Storage

Software
R Studio Version 2025.5.1.513
 and R 4.5.1 


R Packages (Bioconductor)
 DESeq2 (1.48.2)
 GEOquery (2.76.0)
 ggplot2 (4.0.0)
 pheatmap (1.0.13)
 org.At.tair.db (3.21.0)
 EnhancedVolcano (1.26.0)
AnnotationDbi (1.70.0)
clusterProfiler (4.18.1)

Datasets
Drought stress RNA-seq data: NCBI GEO accession GSE270544
Arabidopsis thaliana AND Roots AND Drought
4 Drought-Treated Samples and 4 Control Samples
Platform: RNA Sequencing
Attributed: van Hooren MJ, Munnik T (see Bibliography)
Salt stress RNA-seq data: NCBI GEO accession GSE271344
Arabidopsis thaliana AND Roots AND Salt
3 Salt-Treated Samples and 3 Control Samples
Platform: RNA-seq data / 2025
	Attributed:  Nunez-Vazquez R, Madeira S, Rodriguez-Casillas L, Gomez-Martinez D, Desvoyes B, Gutierrez C (see Bibliography)
  
Additional Resources
TAIR Database link
Gene Ontology Database 
KEGG Pathway Database



1.)Download R  4.5.1 and R Studio 2025.5.1.513
2.) Install Bioconductor packages
3.) Download Datasets from NCBI GEO
	Import Datasets into R Studio
	Format count matrix data into .csv
	Run DESeq2
	Display Volcano Plot, PCA, Dot plots, KEGG Analysis

  Run:
  1.) Drought Script.R
  2.) Salt Script.R
  3.)Comparative Script.R
  4.) Hypothesis Script.R
Parameters
Adjusted P-value < 0.05 Cutoff for statistical significance 
|log2FC| > 1 for amount or percentage changed.

