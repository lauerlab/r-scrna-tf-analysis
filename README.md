# TF Analysis Package

Simple R functions for transcription factor analysis using single-cell RNA-seq data.

## Installation

1. Clone this repository
2. Install required R packages:
```r
install.packages(c("Seurat", "dplyr", "igraph", "ggraph", "ggplot2", "stringr"))
```

## Files

- `tf_analysis.R` - Main functions
- `example_usage.R` - Example usage script
- `README.md` - This file

## Quick Start

```r
# Load functions
source("tf_analysis.R")

# Load your data
seurat_obj <- readRDS("your_seurat_object.rds")
transcription_factors <- c("SATB1", "BACH2", "TBX21", "EOMES")

# Load databases
db1 <- load_tf_database("hTFtarget.csv", "TF", "target", filter_column = "tissue", filter_value = "blood")
db2 <- load_tf_database("TFTG_Curate.txt", "TF", "gene", sep = "\t")

# Run analysis
results <- run_tf_analysis(seurat_obj, "group1", "group2", transcription_factors, list(db1, db2))
```

## TF Databases

Download these databases and place them in your working directory:

### 1. hTFtarget Database
- **File**: `hTFtarget.csv`
- **Download**: [Link to abstract](https://pubmed.ncbi.nlm.nih.gov/32858223/)

### 2. TFTG Curated Database
- **Download**: [Link to abstract](https://pubmed.ncbi.nlm.nih.gov/38707542/)

### 3. TFLink Database
- **Download**: [Link to abstract](https://pubmed.ncbi.nlm.nih.gov/36124642/)

## Main Functions

### `load_tf_database()`
Load any TF database by specifying column names:
```r
db <- load_tf_database(file_path, tf_column, target_column, sep = "auto")
```

### `perform_de_analysis()`
Differential expression analysis between two groups:
```r
de_results <- perform_de_analysis(seurat_obj, "group1", "group2")
```

### `create_network()`
Create TF-target network:
```r
network <- create_network(tf_database, tf_genes, target_genes)
```

### `run_tf_analysis()`
Complete analysis workflow:
```r
results <- run_tf_analysis(seurat_obj, ident1, ident2, transcription_factors, databases)
```

## Output

The analysis creates:
- `differential_expression.csv` - DE results
- `tf_genes_up.csv` - Upregulated TF genes
- `target_genes_up.csv` - Upregulated target genes
- `network_plot.png` - Network visualization
- `network_edges.csv` - Network edge data
- `network_vertices.csv` - Network node data
- `top_hubs.csv` - Top hub genes

## Example Usage

See `example_usage.R` for complete examples showing how to:
- Load different database formats
- Set up your Seurat object
- Run the analysis
- Access results

## Requirements

- R
- Seurat object with cell identities set
- TF-target databases (see above)
- List of transcription factor genes