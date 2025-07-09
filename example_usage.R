# =============================================================================
# Example Usage of TF Analysis Functions
# =============================================================================

# Load the functions
source("tf_analysis.R")

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# 1. Load your Seurat object (replace with your actual file)
# seurat_obj <- readRDS("your_seurat_object.rds")

# 2. Define transcription factors (replace with your actual TF list, better include ALL transcription factors that intersect with the full gene list in the object)
transcription_factors <- c("SATB1", "BACH2", "TBX21", "EOMES", "RUNX3", "GATA3")

# 3. Load TF databases using the universal function
# Replace file paths and column names with your actual database files

# Example 1: hTFtarget database
# Download from: [REPLACE WITH ACTUAL LINK]
db1 <- load_tf_database(
  file_path = "data/hTFtarget.csv",
  tf_column = "TF",
  target_column = "target",
  filter_column = "tissue",
  filter_value = "blood"
)

# Example 2: TFTG database
# Download from: [REPLACE WITH ACTUAL LINK]  
db2 <- load_tf_database(
  file_path = "data/TFTG_Curate.txt",
  tf_column = "TF",
  target_column = "gene",
  sep = "\t"
)

# Example 3: TFLink database
# Download from: [REPLACE WITH ACTUAL LINK]
db3 <- load_tf_database(
  file_path = "data/TFLink_SmallScale.tsv",
  tf_column = "Name.TF",
  target_column = "Name.Target",
  sep = "\t"
)

# Example 4: Custom database
# If you have your own database format:
# db4 <- load_tf_database(
#   file_path = "data/my_custom_db.csv",
#   tf_column = "transcription_factor",
#   target_column = "target_gene",
#   filter_column = "species",
#   filter_value = "human"
# )

# 4. Combine databases
databases <- list(db1, db2, db3)

# 5. Set up your Seurat object
# Make sure your Seurat object has the right identities set
# Example: Idents(seurat_obj) <- seurat_obj$infection_status

# 6. Run the complete analysis
results <- run_tf_analysis(
  seurat_obj = seurat_obj,
  ident1 = "Group_1",                    # Replace with your group 1
  ident2 = "Group_2",                 # Replace with your group 2
  transcription_factors = transcription_factors,
  databases = databases,
  output_dir = "tf_analysis_results",
  filter_mode = "both"              # Options: "both", "tf_only", "any"
)

# 7. View results
print("Network Statistics:")
print(results$stats)

print("Top Hub Genes:")
print(results$hubs)

print("TF Genes Upregulated:")
print(results$tf_genes)

print("Target Genes Upregulated:")
print(results$target_genes)

# 8. Display plots
print(results$plots$network)

# Access individual results
head(results$de_results)
results$plots$summary$degree_dist
results$plots$summary$top_hubs

# =============================================================================
# STEP-BY-STEP USAGE (if you want more control)
# =============================================================================

# Alternative: Run analysis step by step
# 1. Differential expression
de_results <- perform_de_analysis(seurat_obj, "Group_1", "Group_2")
de_genes <- extract_de_genes(de_results)

# 2. Get TF and target genes
tf_up <- intersect(de_genes$up, transcription_factors)
target_up <- setdiff(de_genes$up, transcription_factors)

# 3. Combine databases
combined_db <- combine_databases(databases)

# 4. Create network
network <- create_network(combined_db, tf_up, target_up, filter_mode = "both")

# 5. Analyze network
stats <- network_stats(network)
hubs <- get_top_hubs(network, n = 10)

# 6. Visualize
network_plot <- plot_network(network, title = "My TF Network")
summary_plots <- plot_network_summary(network)

print(network_plot)
print(summary_plots$degree_dist)
print(summary_plots$top_hubs)