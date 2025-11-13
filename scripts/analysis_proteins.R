library(tidyverse)
library(ggplot2)
library(limma)
library(UniprotR)
library(circlize)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(pheatmap)
library(phyloseq)
library(VennDiagram)

# **Load Data**
supernatant <- read.csv("supernatant_data.csv")
intracellular <- read.csv("intracellular_data.csv")
supernatant_abundance <- read.csv("Supernatant_Mutualproteins.csv")
intracellular_abundance <- read.csv("Intrasample_MutualProteins.csv")

# Create protein_info table for mapping
protein_info <- bind_rows(supernatant, intracellular) %>%
  distinct(From, Protein.names1, .keep_all = FALSE) %>%
  rename(Accession = From, ProteinName = Protein.names1)

write.csv(protein_info, "Combined_prot_list.csv")

# **Get the Ids**
supernatant_ids <- supernatant$From
intracellular_ids <- intracellular$From

intersection <- intersect(supernatant_ids, intracellular_ids)

output_dir <- "limma_outputs_vsn"

# make sure the folder exists
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

venn.plot <- venn.diagram(
  x = list(
    Supernatant   = supernatant_ids,
    Intracellular = intracellular_ids
  ),
  filename    = NULL,
  fill        = c("#298c8c", "#f1a226"),
  alpha       = 0.6,
  scaled      = FALSE,
  cex         = 1.2,
  cat.cex     = 1.0,
  cat.fontface = "bold",
  fontfamily   = "sans",
  cat.fontfamily = "sans",
  cat.pos     = c(30, -30),
  cat.dist    = c(0.06, 0.06),
  cat.just    = list(c(0.5, 1), c(0.5, 1)),
  margin      = 0.05
)

png(file.path(output_dir, "venn.png"), width = 2000, height = 2000, res = 300)
grid.newpage()
grid.draw(venn.plot)
dev.off()

supernatant_intersection <- supernatant |>
  filter(From %in% intersection)

intracellular_intersection <- intracellular |> 
  filter(From %in% intersection)

# **Create the abundance dataframe for the intersected Ids**
supernatant_abundance_1 <- supernatant_abundance[,-c(5:8)]
dim(supernatant_abundance_1)
intracellular_abundance_1 <- intracellular_abundance[,-c(5:8)]

intersection_abundance_super <- supernatant_abundance_1 |> 
  filter(Accession %in% intersection)
dim(intersection_abundance_super)
intersection_abundance_intra <- intracellular_abundance_1 |> 
  filter(Accession %in% intersection)

intersection_abundance <- merge(intersection_abundance_super, intersection_abundance_intra, 
                                by.x="Accession", by.y="Accession",
                                all.x = TRUE, all.y = TRUE)
dim(intersection_abundance)

intersection_abundance[is.na(intersection_abundance)] <- 0


# **Create the metadata**
metadata <- data.frame(
  run = c("Supernatant1", "Supernatant2",
  "Supernatant3", "Intracellular1",
  "Intracellular2","Intracellular3"),
  group = c("Supernatant", "Supernatant", "Supernatant", 
            "Intracellular", "Intracellular", "Intracellular")
)

metadata$group <- as.factor(metadata$group)
metadata <- as.data.frame(metadata)
rownames(intersection_abundance) <- intersection_abundance$Accession
intersection_abundance <- intersection_abundance[,-1]
colnames(intersection_abundance) <- c("Supernatant1", "Supernatant2",
                                      "Supernatant3", "Intracellular1",
                                      "Intracellular2","Intracellular3")

# Function to detect metadata columns dynamically
detect_metadata_column <- function(metadata, condition_names, sample_names) {
  # Convert to standard dataframe if metadata is a tibble
  metadata <- as.data.frame(metadata)  # Ensures proper renaming
  
  # Print available columns
  message("Checking metadata for valid columns...")
  message("Available columns in metadata: ", paste(colnames(metadata), collapse = ", "))
  
  # Detect condition column
  detected_condition <- condition_names[condition_names %in% colnames(metadata)][1]
  if (is.na(detected_condition)) {
    stop(paste("ERROR: No valid condition column found! Available columns:", paste(colnames(metadata), collapse = ", ")))
  }
  
  # Detect sample column (including "Run")
  detected_sample <- sample_names[sample_names %in% colnames(metadata)][1]
  
  if (is.na(detected_sample)) {
    message("ERROR: Could not find a valid sample ID column!")
    message("Expected sample ID column names: ", paste(sample_names, collapse = ", "))
    stop("Sample ID column not found! Please check metadata formatting.")
  }
  
  # Rename the detected sample column to "run"
  colnames(metadata)[colnames(metadata) == detected_sample] <- "run"
  
  # Debugging output
  message(paste("Detected condition column:", detected_condition))
  message(paste("Detected sample ID column:", detected_sample, " (Renamed to 'run')"))
  
  return(list(metadata = metadata, condition_col = detected_condition))
}

# Function to ensure condition levels are valid
validate_condition_levels <- function(condition) {
  condition <- droplevels(as.factor(condition))  # Convert to factor and drop unused levels
  
  # Print condition levels for debugging
  message("Final condition levels: ", paste(levels(condition), collapse=", "))
  
  # Check for NA values
  if (any(is.na(condition))) {
    stop("ERROR: Condition column contains NA values. Please check metadata!")
  }
  
  # Ensure exactly two unique condition levels
  if (length(levels(condition)) < 2) {
    stop("ERROR: At least two unique conditions are required!")
  }
  if (length(levels(condition)) > 2) {
    stop(paste("ERROR: More than two condition levels detected:", paste(levels(condition), collapse=", "), 
               "\nUse run_corncob_multiclass() instead."))
  }
  
  return(condition)
}

# Function to log messages
log_message <- function(log_file, message) {
  con <- file(log_file, open = "a")
  writeLines(paste(Sys.time(), "-", message), con)
  close(con)
}

generate_volcano_plot <- function(results_df, output_path, comparison_label, method_name) {
  # Check input validity
  if (!is.data.frame(results_df)) {
    stop("results_df must be a data frame.")
  }
  
  required_columns <- c("log2FoldChange", "pvalue", "padj")
  if (!all(required_columns %in% colnames(results_df))) {
    stop("results_df must contain the following columns: log2FoldChange, pvalue, padj.")
  }
  
  if (is.null(comparison_label) || comparison_label == "" || is.null(method_name) || method_name == "") {
    stop("comparison_label and method_name must be non-empty strings.")
  }
  
  # Filter out rows with NA values
  results_df <- results_df %>%
    filter(!is.na(log2FoldChange) & !is.na(pvalue) & !is.na(padj))
  
  if (nrow(results_df) == 0) {
    stop("Filtered results_df has no valid rows to plot.")
  }
  
  # Ensure output directory exists
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  
  # Define the filename
  output_file <- file.path(output_path, paste0("volcano_", method_name, "_", comparison_label, ".png"))
  
  # Define significance
  results_df$Significance <- ifelse(results_df$padj < 0.05, "Significant", "Not Significant")
  
  # Create volcano plot
  p <- ggplot(results_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Significance)) +
    geom_point(alpha = 0.7, na.rm = TRUE) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    theme_minimal() +
    labs(
      title = paste("Volcano Plot -", comparison_label, "(", method_name, ")"),
      x = "Log2 Fold Change",
      y = "-Log10(P-value)",
      color = "Significance"
    ) +
    theme(plot.title = element_text(hjust = 0.5, size = 14))
  
  # Save the plot
  ggsave(filename = output_file, plot = p, width = 8, height = 6, dpi = 300)
}

generate_heatmap <- function(count_matrix, sig_genes, metadata, condition_col, output_path, comparison_label, method_name) {
  # Ensure the directory exists
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  
  # Check inputs
  if (!is.matrix(count_matrix) && !is.data.frame(count_matrix)) {
    stop("count_matrix must be a matrix or data frame.")
  }
  
  if (!is.data.frame(sig_genes) || !("gene_id" %in% colnames(sig_genes))) {
    stop("sig_genes must be a data frame containing a column named 'gene_id'.")
  }
  
  if (!is.data.frame(metadata) || !(condition_col %in% colnames(metadata))) {
    stop("metadata must be a data frame containing the condition column specified.")
  }
  
  if (nrow(sig_genes) < 2) {
    message("Not enough significant genes for heatmap. Skipping...")
    return(NULL)
  }
  
  # **Filter count matrix for significant genes**
  # Match original IDs to protein names
  protein_names_map <- sig_genes %>%
    distinct(gene_id, label) %>%
    tibble::deframe()  # gene_id will become names, label will become values
  
  
  # Filter the matrix
  sig_gene_matrix <- count_matrix[rownames(count_matrix) %in% names(protein_names_map), , drop = FALSE]
  
  # Replace rownames with protein names
  rownames(sig_gene_matrix) <- protein_names_map[rownames(sig_gene_matrix)]
  
  
  if (nrow(sig_gene_matrix) < 2 || ncol(sig_gene_matrix) < 2) {
    message("Too few genes or samples after filtering. Skipping heatmap...")
    return(NULL)
  }
  
  # **Ensure sample names match between metadata and count matrix**
  valid_samples <- intersect(colnames(sig_gene_matrix), metadata$run)
  
  if (length(valid_samples) < 2) {
    message("ERROR: Sample name mismatch detected!")
    message("Samples in count_matrix: ", paste(colnames(sig_gene_matrix), collapse = ", "))
    message("Samples in metadata: ", paste(metadata$run, collapse = ", "))
    message("Valid samples after intersection: ", paste(valid_samples, collapse = ", "))
    return(NULL)
  }
  
  # **Filter metadata and count matrix**
  sig_gene_matrix <- sig_gene_matrix[, valid_samples, drop = FALSE]
  metadata_filtered <- metadata[metadata$run %in% valid_samples, , drop = FALSE]
  
  # **Ensure metadata has conditions**
  metadata_filtered[[condition_col]] <- as.character(metadata_filtered[[condition_col]])
  unique_conditions <- unique(metadata_filtered[[condition_col]])
  
  if (length(unique_conditions) < 2) {
    message("Not enough conditions to generate heatmap. Skipping...")
    return(NULL)
  }
  
  # **Define condition colors**
  num_conditions <- length(unique_conditions)
  condition_colors <- if (num_conditions < 3) {
    setNames(c("#298c8c", "#f1a226")[1:num_conditions], unique_conditions)
  } else {
    setNames(brewer.pal(min(9, num_conditions), "Set3"), unique_conditions)
  }
  
  # **Replace column names with host_subject_id if available**
  if ("host_subject_id" %in% colnames(metadata_filtered)) {
    sample_to_subject <- setNames(metadata_filtered$host_subject_id, metadata_filtered$run)
    colnames(sig_gene_matrix) <- sample_to_subject[colnames(sig_gene_matrix)]
  }
  
  # **Normalize data (log transform and scale)**
  sig_gene_matrix <- log2(sig_gene_matrix + 1)
  sig_gene_matrix <- t(scale(t(sig_gene_matrix), center = TRUE, scale = TRUE))
  
  # **Define heatmap color scale**
  heatmap_colors <- colorRamp2(c(-2, 0, 2), c("#2c7bb6", "#ffffbf", "#d7191c"))
  
  # **Create Condition Annotation (Show Names, Hide Legend)**
  condition_annotation <- HeatmapAnnotation(
    Condition = factor(metadata_filtered[[condition_col]], levels = unique_conditions),
    col = list(Condition = condition_colors),
    annotation_name_gp = gpar(fontsize = 11, fontface = "bold")  
  )
  
  # **Correct file path**
  output_file <- file.path(output_path, paste0("heatmap_", method_name, "_", gsub(" ", "_", comparison_label), ".png"))
  
  # **Generate heatmap**
  png(output_file, width = 2400, height = 1400, res = 300)
  draw(Heatmap(
    sig_gene_matrix,
    name = "Expression",
    col = heatmap_colors,
    cluster_rows = TRUE,  # Cluster genes
    cluster_columns = FALSE,  # Cluster samples within each condition
    show_row_names = TRUE,
    show_column_names = FALSE,
    column_names_rot = 90,
    row_names_gp = gpar(fontsize = 8.5),
    column_names_gp = gpar(fontsize = 8.5),
    column_split = factor(metadata_filtered[[condition_col]], levels = unique_conditions),
    top_annotation = condition_annotation
  ))
  dev.off()
}

run_limma_voom_binary <- function(count_matrix, metadata, condition_col, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # **Detect & Handle Non-Sample First Column**
  possible_non_sample_cols <- c("Taxonomic_Lineage", "OTU", "Tax", "Gene", "Genes", "Gene_Family", "Accession")
  
  first_col_name <- colnames(count_matrix)[1]
  if (first_col_name %in% possible_non_sample_cols) {
    message(paste("Detected non-sample column:", first_col_name, "- Moving to row names."))
    rownames(count_matrix) <- count_matrix[[first_col_name]]
    count_matrix <- count_matrix[, -1]
  } else {
    message("First column is already sample data. No changes needed.")
  }
  
  # Convert to matrix
  count_matrix <- as.matrix(count_matrix)
  
  # **Detect metadata columns**
  metadata_info <- detect_metadata_column(metadata, 
                                          condition_names = c("Condition", "group", "Class"),
                                          sample_names = c("Run", "run", "sampleid", "SampleID", "sample", "sample_id", "ID"))
  
  metadata <- metadata_info$metadata
  condition_col <- metadata_info$condition_col
  
  condition <- validate_condition_levels(metadata[[condition_col]])
  
  # **Ensure Sample Names Match**
  valid_samples <- intersect(metadata$run, colnames(count_matrix))
  if (length(valid_samples) < 2) {
    stop("ERROR: No matching sample names between metadata and count_matrix! Check formatting.")
  }
  
  # **Filter count_matrix and metadata**
  count_matrix <- count_matrix[, valid_samples, drop = FALSE]
  metadata_filtered <- metadata[metadata$run %in% valid_samples, , drop = FALSE]
  
  # **Prepare Limma-Voom Data**
  group <- factor(metadata_filtered[[condition_col]])
  design <- model.matrix(~ group)
  
  # # **Apply Voom Transformation**
  # v <- voom(count_matrix, design)
  
  # **Fit Model & Perform Differential Expression**
  fit <- lmFit(count_matrix, design)
  fit <- eBayes(fit)
  
  # **Extract Results**
  results_limma <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
  results_limma$gene_id <- rownames(results_limma)
  print(results_limma)
  # Join protein names to results using Accession
  results_limma <- results_limma %>%
    left_join(protein_info, by = c("gene_id" = "Accession"))
  
  # Set a label column for heatmaps and plots
  results_limma$label <- ifelse(!is.na(results_limma$ProteinName), results_limma$ProteinName, results_limma$gene_id)
  
  # Update rownames for plotting
  rownames(results_limma) <- results_limma$label
  
  # **Rename columns for compatibility**
  results_limma <- results_limma  |> 
    dplyr::rename(log2FoldChange = logFC, pvalue = P.Value, padj = adj.P.Val)
  
  # **Fix: Use condition names in the filename**
  comparison_label <- paste(levels(group)[2], "vs", levels(group)[1])
  method_name <- "limma"
  
  # **Save results**
  write.csv(results_limma, file = file.path(output_dir, paste0("limma_voom_binary_results_", comparison_label, ".csv")), row.names = FALSE)
  
  # **Create Directories for Plots**
  heatmap_output_dir <- file.path(output_dir, "heatmaps")
  volcano_output_dir <- file.path(output_dir, "volcano_plots")
  dir.create(heatmap_output_dir, showWarnings = FALSE)
  dir.create(volcano_output_dir, showWarnings = FALSE)
  
  # **Generate Volcano & Heatmap (With Handling Skipping)**
  sig_genes <- results_limma %>% filter(padj < 0.05)
  
  if (nrow(results_limma) > 0) {
    generate_volcano_plot(results_limma, volcano_output_dir, comparison_label, method_name)
  } else {
    message("Skipping volcano plot - No significant results found.")
  }
  
  if (nrow(sig_genes) > 1) {
    generate_heatmap(count_matrix, sig_genes, metadata_filtered, condition_col, heatmap_output_dir, comparison_label, method_name)
  } else {
    message("Skipping heatmap - Not enough significant genes.")
  }
}

run_limma_voom_binary(count_matrix = intersection_abundance, metadata = metadata, condition_col = "group", output_dir = "limma_outputs_vsn")



