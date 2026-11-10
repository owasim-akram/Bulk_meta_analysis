# Differential gene expression analysis
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135668
# Load packages
library(tidyverse)
library(DESeq2)
library(ggrepel)

# Load raw counts
counts_data <- readxl::read_excel("data/raw_counts/GSE135668_raw_counts.xlsx")
counts_data <- counts_data |>
    column_to_rownames("Geneid")

# Load metadata
metadata <- read.csv("data/metadata/GSE135668_metadata.csv")
glimpse(metadata)

# Prepare Sample information
col_data <- data.frame(condition = as.factor(metadata$Conditions),
                       row.names = colnames(counts_data))
glimpse(col_data)

# Match metadata with count data
metadata |>
    filter(Samples %in% colnames(counts_data)) |>
    arrange(match(Samples, colnames(counts_data)))

# Preapre count matrix
counts_data <- counts_data |>
    mutate_all(as.numeric) |>
    as.matrix()

# Create DESeq2 data object
dds <- DESeqDataSetFromMatrix(
    countData = counts_data,
    colData = col_data,
    design = ~condition
)

# Filter any counts less than 10
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)

# Analyze the dataset and compile results
res <- results(dds)
res$Gene_ID <- sub("gene:", "", rownames(res))

# Export results to a CSV file
write.csv(res, "outputs/DESeq2/GSE135668.csv", row.names = FALSE)

# Get Normalized counts
normalized_counts <- counts(dds, normalized=TRUE) |>
    as.data.frame()

normalized_counts$Gene_ID <- sub("gene:", "", rownames(normalized_counts))
rownames(normalized_counts) <- NULL

# Perform variance-stabilizing transformation (VST) for PCA
vsd <- vst(dds, blind = FALSE)

# Perform PCA on the transformed data
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

# Calculate the percentage of variance explained by each PC
percent_var <- round(100 * attr(pca_data, "percentVar"))

# Define custom colors for each conditions
custom_colors <- c("control" = "#1c9099",
                   "heat stress" = "#de2d26")

# Plot PCA with custom styling
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_text_repel(aes(label = name),
                    box.padding = 0.5,
                    max.overlaps = Inf,
                    size = 3) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    ggtitle("GSE135668") +
    scale_color_manual(values = custom_colors) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  # Add border
        legend.position = "right",
        aspect.ratio = 1
    ) +
    guides(color = guide_legend(title = "condition"))
# Export the plot
ggsave(filename = "figures/PCA/GSE135668.png",
       plot = pca_plot,
       width = 5,
       height = 5,
       units = "in",
       bg = "white",
       dpi = 600)
