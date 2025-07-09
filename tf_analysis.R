# =============================================================================
# TF Analysis Functions
# =============================================================================
# Simple R functions for transcription factor analysis with single-cell data
# 
# Required libraries: Seurat, dplyr, igraph, ggraph, ggplot2, stringr
# =============================================================================

# Load required libraries
library(Seurat)
library(dplyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(stringr)

# =============================================================================
# 1. DIFFERENTIAL EXPRESSION
# =============================================================================

#' Perform differential expression analysis
#' @param seurat_obj Seurat object with identities already set
#' @param ident1 First group
#' @param ident2 Second group
#' @param logfc_threshold Log fold change threshold (default: 0.25)
#' @return Data frame with results
perform_de_analysis <- function(seurat_obj, ident1, ident2, logfc_threshold = 0.25) {
  results <- FindMarkers(seurat_obj, ident.1 = ident1, ident.2 = ident2, 
                        test.use = "wilcox", logfc.threshold = logfc_threshold)
  results <- subset(results, p_val < 0.05)
  results$gene <- rownames(results)
  return(results)
}

#' Extract up/down regulated genes
#' @param de_results Results from perform_de_analysis
#' @param logfc_threshold Threshold for significance
#' @return List with up and down regulated genes
extract_de_genes <- function(de_results, logfc_threshold = 0.25) {
  up_genes <- de_results$gene[de_results$avg_log2FC > logfc_threshold]
  down_genes <- de_results$gene[de_results$avg_log2FC < -logfc_threshold]
  return(list(up = up_genes, down = down_genes))
}

# =============================================================================
# 2. DATABASE LOADER
# =============================================================================

#' Load TF database (universal function for any database format)
#' @param file_path Path to database file
#' @param tf_column Name of column containing TF names
#' @param target_column Name of column containing target gene names
#' @param sep File separator (default: auto-detect)
#' @param header Whether file has header (default: TRUE)
#' @param filter_column Optional column name to filter by
#' @param filter_value Optional value to filter by
#' @return Data frame with TF and target columns
load_tf_database <- function(file_path, tf_column, target_column, 
                             sep = "auto", header = TRUE, 
                             filter_column = NULL, filter_value = NULL) {
  
  # Auto-detect separator if not specified
  if (sep == "auto") {
    if (grepl("\\.csv$", file_path)) {
      sep <- ","
    } else if (grepl("\\.tsv$", file_path)) {
      sep <- "\t"
    } else if (grepl("\\.txt$", file_path)) {
      sep <- "\t"
    } else {
      sep <- ","
    }
  }
  
  # Load file
  db <- read.delim(file_path, header = header, sep = sep, stringsAsFactors = FALSE)
  
  # Check if columns exist
  if (!tf_column %in% colnames(db)) {
    stop(paste("Column", tf_column, "not found in database"))
  }
  if (!target_column %in% colnames(db)) {
    stop(paste("Column", target_column, "not found in database"))
  }
  
  # Apply filter if specified
  if (!is.null(filter_column) && !is.null(filter_value)) {
    if (filter_column %in% colnames(db)) {
      db <- db[db[[filter_column]] == filter_value, ]
      cat("Filtered by", filter_column, "==", filter_value, "\n")
    } else {
      warning(paste("Filter column", filter_column, "not found"))
    }
  }
  
  # Create standardized output
  result <- data.frame(
    TF = db[[tf_column]],
    target = db[[target_column]],
    stringsAsFactors = FALSE
  )
  
  # Clean target column (remove commas if present)
  result$target <- str_replace_all(result$target, ",", "")
  
  # Remove empty rows
  result <- result[!is.na(result$TF) & !is.na(result$target), ]
  result <- result[result$TF != "" & result$target != "", ]
  
  cat("Loaded", nrow(result), "TF-target interactions from", basename(file_path), "\n")
  
  return(result)
}

#' Combine multiple databases
#' @param db_list List of databases (each with TF and target columns)
#' @return Combined database
combine_databases <- function(db_list) {
  combined <- do.call(rbind, db_list)
  combined <- unique(combined)
  cat("Combined databases:", nrow(combined), "unique TF-target interactions\n")
  return(combined)
}

# =============================================================================
# 3. NETWORK FUNCTIONS
# =============================================================================

#' Create TF network
#' @param tf_database Database with TF and target columns
#' @param tf_genes TF genes of interest
#' @param target_genes Target genes of interest
#' @param filter_mode How to filter: "both" (both TF and target in gene sets), 
#'                    "tf_only" (only TF in tf_genes), "any" (TF or target in gene sets)
#' @return igraph network object
create_network <- function(tf_database, tf_genes, target_genes, filter_mode = "both") {
  
  all_genes <- union(tf_genes, target_genes)
  
  # Filter database based on mode
  if (filter_mode == "both") {
    edges <- tf_database[tf_database$TF %in% all_genes & tf_database$target %in% all_genes, ]
  } else if (filter_mode == "tf_only") {
    edges <- tf_database[tf_database$TF %in% tf_genes, ]
  } else if (filter_mode == "any") {
    edges <- tf_database[tf_database$TF %in% all_genes | tf_database$target %in% all_genes, ]
  } else {
    stop("filter_mode must be 'both', 'tf_only', or 'any'")
  }
  
  edges <- unique(edges)
  
  if (nrow(edges) == 0) {
    warning("No edges found with current filtering")
    return(NULL)
  }
  
  # Create nodes
  nodes <- data.frame(name = unique(c(edges$TF, edges$target)), stringsAsFactors = FALSE)
  nodes$type <- "other"
  nodes$type[nodes$name %in% tf_genes] <- "TF"
  nodes$type[nodes$name %in% target_genes] <- "target"
  nodes$type[nodes$name %in% tf_genes & nodes$name %in% target_genes] <- "both"
  
  # Create network
  network <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
  
  cat("Network created:", vcount(network), "nodes,", ecount(network), "edges\n")
  
  return(network)
}

#' Get network statistics
#' @param network igraph network
#' @return List of statistics
network_stats <- function(network) {
  if (is.null(network)) return(NULL)
  
  stats <- list(
    nodes = vcount(network),
    edges = ecount(network),
    density = edge_density(network),
    components = components(network)$no,
    diameter = ifelse(is_connected(network), diameter(network), NA),
    clustering = transitivity(network, type = "global")
  )
  
  # Node type counts
  node_types <- vertex_attr(network, "type")
  stats$node_types <- table(node_types)
  
  # Degree statistics
  degrees <- degree(network)
  stats$degree_mean <- mean(degrees)
  stats$degree_max <- max(degrees)
  
  return(stats)
}

#' Get top hub genes
#' @param network igraph network
#' @param n Number of top hubs to return
#' @return Data frame with top hub genes
get_top_hubs <- function(network, n = 10) {
  if (is.null(network)) return(NULL)
  
  degrees <- degree(network)
  top_degrees <- sort(degrees, decreasing = TRUE)[1:min(n, length(degrees))]
  
  hub_df <- data.frame(
    gene = names(top_degrees),
    degree = top_degrees,
    type = vertex_attr(network, "type")[names(top_degrees)],
    stringsAsFactors = FALSE
  )
  
  return(hub_df)
}

#' Plot network
#' @param network igraph network
#' @param title Plot title
#' @param layout Layout algorithm (default: "fr")
#' @param node_size Node size (default: 4)
#' @param show_labels Whether to show labels (default: TRUE)
#' @return ggplot object
plot_network <- function(network, title = "TF Network", layout = "fr", 
                        node_size = 4, show_labels = TRUE) {
  if (is.null(network)) {
    return(ggplot() + ggtitle("No network to display"))
  }
  
  colors <- c("TF" = "#1f78b4", "target" = "#33a02c", "both" = "#e31a1c", "other" = "#b2df8a")
  
  p <- ggraph(network, layout = layout) +
    geom_edge_link(alpha = 0.6, arrow = arrow(length = unit(0.15, "cm"))) +
    geom_node_point(aes(color = type), size = node_size) +
    scale_color_manual(values = colors) +
    theme_graph() +
    ggtitle(title) +
    theme(legend.position = "bottom")
  
  if (show_labels) {
    p <- p + geom_node_text(aes(label = name), repel = TRUE, size = 3)
  }
  
  return(p)
}

#' Plot network summary
#' @param network igraph network
#' @return List of ggplot objects
plot_network_summary <- function(network) {
  if (is.null(network)) return(list())
  
  colors <- c("TF" = "#1f78b4", "target" = "#33a02c", "both" = "#e31a1c", "other" = "#b2df8a")
  
  # Degree distribution
  degrees <- degree(network)
  node_types <- vertex_attr(network, "type")
  
  degree_df <- data.frame(degree = degrees, type = node_types)
  
  p1 <- ggplot(degree_df, aes(x = degree, fill = type)) +
    geom_histogram(bins = 20, alpha = 0.7) +
    facet_wrap(~type, scales = "free_y") +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    ggtitle("Degree Distribution by Node Type") +
    theme(legend.position = "none")
  
  # Top hubs
  top_hubs <- get_top_hubs(network, 15)
  
  p2 <- ggplot(top_hubs, aes(x = reorder(gene, degree), y = degree, fill = type)) +
    geom_col(alpha = 0.8) +
    coord_flip() +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    ggtitle("Top 15 Hub Genes") +
    xlab("Gene") + ylab("Degree") +
    theme(legend.position = "bottom")
  
  return(list(degree_dist = p1, top_hubs = p2))
}

# =============================================================================
# 4. MAIN WORKFLOW
# =============================================================================

#' Run complete TF analysis
#' @param seurat_obj Seurat object (pre-filtered)
#' @param ident1 First group
#' @param ident2 Second group
#' @param transcription_factors Vector of TF names
#' @param databases List of loaded databases
#' @param output_dir Output directory
#' @param filter_mode Network filtering mode
#' @return List of results
run_tf_analysis <- function(seurat_obj, ident1, ident2, transcription_factors, 
                           databases, output_dir = "results", filter_mode = "both") {
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("Starting TF analysis...\n")
  
  # 1. Differential expression
  cat("1. Performing differential expression analysis...\n")
  de_results <- perform_de_analysis(seurat_obj, ident1, ident2)
  de_genes <- extract_de_genes(de_results)
  
  # 2. Separate TF and non-TF genes
  cat("2. Separating TF and target genes...\n")
  tf_up <- intersect(de_genes$up, transcription_factors)
  target_up <- setdiff(de_genes$up, transcription_factors)
  
  cat("   TF genes upregulated:", length(tf_up), "\n")
  cat("   Target genes upregulated:", length(target_up), "\n")
  
  # 3. Combine databases
  cat("3. Combining databases...\n")
  combined_db <- combine_databases(databases)
  
  # 4. Create network
  cat("4. Creating network...\n")
  network <- create_network(combined_db, tf_up, target_up, filter_mode)
  
  # 5. Analyze network
  cat("5. Analyzing network...\n")
  stats <- network_stats(network)
  hubs <- get_top_hubs(network)
  
  # 6. Create plots
  cat("6. Creating plots...\n")
  network_plot <- plot_network(network)
  summary_plots <- plot_network_summary(network)
  
  # 7. Save results
  cat("7. Saving results...\n")
  write.csv(de_results, file.path(output_dir, "differential_expression.csv"), row.names = FALSE)
  write.csv(data.frame(gene = tf_up), file.path(output_dir, "tf_genes_up.csv"), row.names = FALSE)
  write.csv(data.frame(gene = target_up), file.path(output_dir, "target_genes_up.csv"), row.names = FALSE)
  
  if (!is.null(network)) {
    ggsave(file.path(output_dir, "network_plot.png"), network_plot, width = 12, height = 10)
    if (length(summary_plots) > 0) {
      ggsave(file.path(output_dir, "degree_distribution.png"), summary_plots$degree_dist, width = 10, height = 8)
      ggsave(file.path(output_dir, "top_hubs.png"), summary_plots$top_hubs, width = 10, height = 6)
    }
    
    # Save network data
    write.csv(as_data_frame(network, "edges"), file.path(output_dir, "network_edges.csv"), row.names = FALSE)
    write.csv(as_data_frame(network, "vertices"), file.path(output_dir, "network_vertices.csv"), row.names = FALSE)
  }
  
  if (!is.null(hubs)) {
    write.csv(hubs, file.path(output_dir, "top_hubs.csv"), row.names = FALSE)
  }
  
  cat("Analysis complete! Results saved to:", output_dir, "\n")
  
  # Return results
  return(list(
    de_results = de_results,
    de_genes = de_genes,
    tf_genes = tf_up,
    target_genes = target_up,
    network = network,
    stats = stats,
    hubs = hubs,
    plots = list(network = network_plot, summary = summary_plots)
  ))
}