library(Seurat)
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(HGNChelper)
library(cowplot)

input_single_cell_path <- "/mount/mdtaylor2/harmanci/Medullo_Multiome_RObject/Group3_seuratObj_updated_annotated_cc_060225.rda"
#input_single_cell_path <- "/mount/mdtaylor2/harmanci/Medullo_Multiome_RObject/Group4_seuratObj_updated_annotated_cc_060225.rda"                     # % cutoff for clonal vs subclonal

# Example of dynamic file naming by including the subtype_name in the final output figures
output_path <- '/mount/mdtaylor2/aliceyang/4Shubham'
subtype_name <- "Group 3" # "Group 3" or "Group 4"
category_name <- "Metabolism"  # "Catabolism" or "Metabolism" or "Immunotherapy"


load(input_single_cell_path)        # loads Seurat object "seuratObj" 
objects()

# ===================== Remove control =============================================
#seuratObj <- subset(seuratObj, subset = orig.ident != "control")
# if we remove control all the '_normal' in the seuratObj$final_cell_type2 are gone


# ============== Impute metadata data =============================================
meta <- seuratObj@meta.data

# --- Step 1: Impute the Tetraploid/Diploid value in the metadata dataframe ---
print("Ploidy status for MB0575 BEFORE imputation:")
table(meta$tetraploidy[meta$orig.ident == "MB0575"], useNA = "ifany")
# We find all rows where the sample is 'MB0575' and set 'tetraploidy' to 'Diploid'
meta$tetraploidy[meta$orig.ident == "MB0575"] <- "Diploid"


# --- Step2: Impute the i17q value in the metadata dataframe ---
table(meta$iso_chr17[meta$orig.ident == "MB0423"], useNA = "ifany")
# We find all rows where the sample is 'MB0423' and set 'iso_chr17' to 'yes'
meta$iso_chr17[meta$orig.ident == "MB0423"] <- "yes"


# --- Step3: Impute the gender value in the metadata dataframe ---
# Aggregate XIST expression and return summed counts ("pseudobulk") for each sample
#agg_expr_all <- AggregateExpression(
#object = loaded_data,
#features = "XIST",
#group.by = "sample_id",  # adjust if needed
#slot = "data",           # log-normalized expression
#verbose = FALSE
#)
table(meta$gender[meta$orig.ident == "HSMB0178"], useNA = "ifany")
meta$gender[meta$orig.ident == "HSMB0178"] <- "M"

table(meta$gender[meta$orig.ident == "MB2032"], useNA = "ifany")
meta$gender[meta$orig.ident == "MB2032"] <- "M"


table(meta$gender[meta$orig.ident == "MB2840"], useNA = "ifany")
meta$gender[meta$orig.ident == "MB2840"] <- "M"



# review sample_info
sample_info <- meta %>%
  select(orig.ident, gender, tetraploidy, iso_chr17) %>%
  distinct()  # One row per sample
print(sample_info)

# --- CRITICAL STEP ---
# Add the modified metadata back into the Seurat object
seuratObj@meta.data <- meta

# ===================== Sanity Check: Single Cell RNA Data Basic Stats ====================

# count number of subjects (including control)
length(unique(seuratObj$orig.ident)) # 35 group 3 (1 control) 84 (1 control)
table(seuratObj$orig.ident) # cells per samples

# Set default assay
DefaultAssay(seuratObj) <- "RNA"


Idents(seuratObj)
# Levels: G1 G2M S
# seuratObj = RenameIdents(seuratObj, 'G2M' = 'G2/M')


# ===================== Get Cell Type Annotations: Focus on Normal =============================================
# Ensure grouping column is a factor
#seuratObj$final_cell_type2 <- factor(seuratObj$final_cell_type2)

unique(seuratObj$final_cell_type2)
# Check the result
#table(seuratObj$status)


# ===================== Visualization =============================================

# We'll focus on the most relevant cell types for a clear plot
# You can change these to any cell types of interest
#unique(seuratObj$final_cell_type2)
#[1] "GC_GN"              "RL-VZ_cluster2"     "RL-VZ_cluster1"    
#[4] "OPC"                "RL-SVZ_photoR_cls1" "RL-SVZ_cycling"    
#[7] "UBC"                "Astro"              "immune"            
#[10] "RL-SVZ_non_cycling" "RL-SVZ_photoR_cls2" "Endothelial"       
#[13] "GCP"                "UBC_normal"         "GCP_normal"        
#[16] "GC_normal"          "RL-VZ_normal"       "eDCN_normal"       
#[19] "RL-SVZ_normal" 
cell_types_of_interest_photoR <- c("RL-SVZ_photoR_cls1", "RL-SVZ_photoR_cls2")
cell_types_of_interest_test <- c("RL-VZ_cluster1", "RL-VZ_cluster2", "RL-SVZ_cycling", "RL-SVZ_non_cycling", "UBC", "GCP", "GC_GN")
cell_types_of_interest_control <- c("RL-VZ_normal", "RL-SVZ_normal", "UBC_normal", "GCP_normal", "GC_normal")


seuratObj$merged_cell_type <- seuratObj$final_cell_type2
seuratObj$merged_cell_type[seuratObj$merged_cell_type %in% c("RL-SVZ_cycling", "RL-SVZ_non_cycling")] <- "RL-SVZ"




# ===================== Standardize Genes in HGNC Database =====================================
#genes_of_interest <- c("CD276", "GPC2", "ERBB2", "PVR", "IL13RA2")
genes_of_interest <- c("TMEM68", "DGAT1", "DGAT2", "PNPLA2", "PNPLA3", "DDHD2", "PLIN1","PLIN2", "PLIN3","PLIN4", "PLIN5", "BSCL2", "LDAF1")
#genes_of_interest <- c(
#"CERS1", "CERS2", "CERS3", "CERS4", "CERS5", "CERS6",
#"KDSR", "DEGS1", "DEGS2", "SGMS1", "SGMS2", "UGCG", 
#"ST3GAL5", "B4GALT5", "B4GALT6", "ST8SIA1", "ST8SIA2", "ST8SIA3", "ST8SIA5",
#"B4GALNT1", "B3GALT4", "ST3GAL2", "ST3GAL3",
#"NEU1", "NEU3", "NEU4", "GLB1", "GBA1", "ASAH1", "PSAP", "HEXA", "HEXB", "GALC", "GM2A"
#)

# Get the actual repeated values
unique(genes_of_interest[duplicated(genes_of_interest)])

# Check and standardize symbols
genes_of_interest_standard_status <- checkGeneSymbols(genes_of_interest)
genes_checked <- checkGeneSymbols(genes_of_interest)

# Report results
cat("=== Gene Symbol Validation Report ===\n\n")

# Valid symbols
valid_genes <- genes_checked[genes_checked$Approved, "x"]
cat("✅ Approved gene symbols:\n")
cat(paste(valid_genes, collapse = ", "), "\n\n")

# Invalid symbols
invalid_genes <- genes_checked[!genes_checked$Approved, ]
if (nrow(invalid_genes) > 0) {
  cat("❌ Invalid or outdated symbols:\n")
  cat(paste(invalid_genes$x, collapse = ", "), "\n\n")
  
  cat("🔁 Suggested corrections:\n")
  for (i in 1:nrow(invalid_genes)) {
    cat(sprintf("- %s → %s\n", invalid_genes$x[i], invalid_genes$Suggested.Symbol[i]))
  }
} else {
  cat("🎉 All gene symbols are approved!\n")
}


# Create the violin plot, splitting by status
#VlnPlot(seuratObj,
        #features = "DDHD2",
        #subset = final_cell_type2 %in% cell_types_of_interest, # Focus on specific cells
        #group.by = "final_cell_type2", # Group by the detailed cell types
        #split.by = "status",           # Split the violins by Tumor/Normal
        #pt.size = 0)                   # Hides individual dots for a cleaner look



# Visualize DDHD2 expression across your main functional categories
#VlnPlot(seuratObj,
        #features = "DDHD2",
        #group.by = "cell_category", # Use the broad categories
        #pt.size = 0.1,              # Show some dots
        #log = TRUE)                 # Use a log scale for expression if values are skewed

# ===================== 3. Generate Split Dot Plot with Continuous Color Scale =============================================

# Define cell type order (top to bottom)
y_axis_order <- c("RL-VZ", "RL-SVZ", "UBC", "GCP", "GN")

# Subset data for Normal and Tumor
normal_subset <- subset(seurat_subset, status == "Normal")
tumor_subset <- subset(seurat_subset, status == "Tumor")

# Set factor levels for both
normal_subset$plot_cell_type <- factor(normal_subset$plot_cell_type, levels = y_axis_order)
tumor_subset$plot_cell_type <- factor(tumor_subset$plot_cell_type, levels = y_axis_order)

# Create Normal plot with FIXED bottom margin
plot_normal <- DotPlot(
  normal_subset,
  features = genes_of_interest,
  group.by = "plot_cell_type",
  dot.scale = 8
) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Average\nExpression",
    limits = c(-2, 2)
  ) +
  scale_x_discrete(expand = c(0.05, 0.05)) +  # Force same x-axis expansion
  labs(y = "Cell Type", x = NULL, title = "Normal") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(10, 2, 50, 10)  # FIXED bottom margin (50)
  )

# Create Tumor plot with SAME FIXED bottom margin
plot_tumor <- DotPlot(
  tumor_subset,
  features = genes_of_interest,
  group.by = "plot_cell_type",
  dot.scale = 8
) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Average\n Scaled Expression",
    limits = c(-2, 2)
  ) +
  scale_x_discrete(expand = c(0.05, 0.05)) +  # Force same x-axis expansion
  labs(y = NULL, x = NULL, title = "Tumor") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(10, 10, 50, 2)  # SAME FIXED bottom margin (50)
  ) +
  guides(
    size = guide_legend(
      title = "% Expressing",
      override.aes = list(color = "grey50"),
      order = 1
    ),
    colour = guide_colorbar(
      title = "Average\nExpression",
      barwidth = 1.2,
      barheight = 10,
      order = 2
    )
  )

# Remove legend from tumor plot temporarily
plot_tumor_temp <- plot_tumor + theme(legend.position = "none")

# Align the two plots perfectly
aligned <- align_plots(plot_normal, plot_tumor_temp, align = "hv", axis = "tblr")

# Extract legend from original tumor plot
legend <- get_legend(plot_tumor)

# Combine aligned plots
plots_row <- plot_grid(
  aligned[[1]],  # Normal plot
  aligned[[2]],  # Tumor plot (no legend)
  ncol = 2,
  rel_widths = c(1, 1)
)

# Add legend
combined_plot <- plot_grid(plots_row, legend, ncol = 2, rel_widths = c(10, 1))

print(combined_plot)
# ===================== 4. Save the Plot =============================================
output_file <- file.path(
  output_path,
  paste0(
    "dotplot_",
    gsub(" ", "_", subtype_name),
    "_",
    category_name,
    "_Normal_Tumor_continuous_color8.png"
  )
)

ggsave(
  filename = output_file,
  plot = combined_plot,
  width = 16,
  height = 7,
  dpi = 300,
  bg = "white"
)

cat("✅ Split DotPlot with continuous color scale saved to:", output_file, "\n")