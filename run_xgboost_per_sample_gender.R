library(Seurat)
library(caret)
library(R6)
library(dplyr)
library(devtools)
library(ggplot2)
library(xgboost)

set.seed(123)  # for reproducibility

# genes located in the X chromosome that have been reported to escape X-inactivation
# http://bioinf.wehi.edu.au/software/GenderGenes/index.html
Xgenes<- c("ARHGAP4","STS","ARSD", "ARSL", "AVPR2", "BRS3", "S100G",
           "CHM", "CLCN4", "DDX3X","EIF1AX","EIF2S3", "GPM6B",
           "GRPR", "HCFC1", "L1CAM", "MAOA", "MYCLP1", "NAP1L3",
           "GPR143", "CDK16", "PLXNB3", "PRKX", "RBBP7", "RENBP",
           "RPS4X", "TRAPPC2", "SH3BGRL", "TBL1X","UBA1", "KDM6A",
           "XG", "XIST", "ZFX", "PUDP", "PNPLA4", "USP9X", "KDM5C",
           "SMC1A", "NAA10", "OFD1", "IKBKG", "PIR", "INE2", "INE1",
           "AP1S2", "GYG2", "MED14", "RAB9A", "ITM2A", "MORF4L2",
           "CA5B", "SRPX2", "GEMIN8", "CTPS2", "CLTRN", "NLGN4X",
           "DUSP21", "ALG13","SYAP1", "SYTL4", "FUNDC1", "GAB3",
           "RIBC1", "FAM9C","CA5BP1")

# genes belonging to the male-specific region of chromosome Y (unique genes)
# http://bioinf.wehi.edu.au/software/GenderGenes/index.html
Ygenes<-c("AMELY", "DAZ1", "PRKY", "RBMY1A1", "RBMY1HP", "RPS4Y1", "SRY",
          "TSPY1", "UTY", "ZFY","KDM5D", "USP9Y", "DDX3Y", "PRY", "XKRY",
          "BPY2", "VCY", "CDY1", "EIF1AY", "TMSB4Y","CDY2A", "NLGN4Y",
          "PCDH11Y", "HSFY1", "TGIF2LY", "TBL1Y", "RPS4Y2", "HSFY2",
          "CDY2B", "TXLNGY","CDY1B", "DAZ3", "DAZ2", "DAZ4")

### 1. DataLoader ###
DataLoader <- R6Class("DataLoader",
                      public = list(
                        filename = NULL,
                        base_filename = NULL,
                        meta_data = NULL,
                        
                        initialize = function(filename) {
                          self$filename <- filename
                          self$base_filename <- gsub("\\.rda$", "", filename)
                          self$load_data()
                        },
                        
                        load_data = function() {
                          if (!file.exists(self$filename)) stop(paste("❌ Missing file:", self$filename))
                          loaded <- load(self$filename)
                          self$meta_data <- get(loaded[1])
                        }
                      )
)


### 1. DataLoader ###
DataLoader <- R6Class("DataLoader",
                      public = list(
                        filename = NULL,
                        base_filename = NULL,
                        meta_data = NULL,
                        
                        initialize = function(filename) {
                          self$filename <- filename
                          self$base_filename <- gsub("\\.rda$", "", filename)
                          self$load_data()
                        },
                        
                        load_data = function() {
                          if (!file.exists(self$filename)) stop(paste("❌ Missing file:", self$filename))
                          loaded <- load(self$filename)
                          self$meta_data <- get(loaded[1])
                        }
                      )
)
# Load the dataset
filename <- "Primary_Group3_seuratObj_ann.rda"
#filename <- "Primary_Group4_seuratObj_ann.rda"
# Set group name dynamically
subgroup_name <- "group 3 MB"
gene_name <- "XIST"

# Create a DataLoader object
my_loader <- DataLoader$new(filename)

# Access the loaded data
loaded_data <- my_loader$meta_data

# Visualize QC metrics as a violin plot
#VlnPlot(loaded_data, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(loaded_data, features = "nFeature_RNA", ncol = 1, pt.size = 0.2) +
  theme(plot.margin = unit(c(1, 1, 1, 2), "cm"))
VlnPlot(loaded_data, features = "nCount_RNA", ncol = 1, pt.size = 0.2) +
  theme(plot.margin = unit(c(1, 1, 1, 2), "cm"))
VlnPlot(loaded_data, features = "percent.mt", ncol = 1, pt.size = 0.2) +
  theme(plot.margin = unit(c(1, 1, 1, 2), "cm"))



# Check gene names in the raw counts matrix
raw_genes <- rownames(loaded_data[["RNA"]]@counts) 


#Search for EIF2S3Y and XIST (case-insensitive)
grep(gene_name, raw_genes, value = TRUE, ignore.case = TRUE)
#grep("EIF2S3", raw_genes, value = TRUE, ignore.case = TRUE)

# FeaturePlot QC step: If XIST is expressed in female samples only, you should see specific clusters light up
plot_title <- paste(gene_name, "gene expression for", subgroup_name)
FeaturePlot(loaded_data, features = gene_name, reduction = "umap", pt.size = 0.5)  +
  ggtitle(plot_title)

# Aggregate XIST expression and return summed counts ("pseudobulk") for each sample
agg_expr_all <- AggregateExpression(
  object = loaded_data,
  features = "XIST",
  group.by = "sample_id",  # adjust if needed
  slot = "data",           # log-normalized expression
  verbose = FALSE
)
# sample_id

length(unique(loaded_data$sample_id))



# Result is a list of matrices (1 per assay); we usually use the first one
expr_matrix_all <- agg_expr_all[[1]]  # genes as rows, samples as columns

# Transpose: make rows = samples, columns = genes
expr_df_all <- as.data.frame(t(expr_matrix_all))
expr_df_all$orig.ident <- rownames(expr_df_all)  # add sample ID for merging

# input is 27 samples including control
# nrow(expr_df_all)
#anti_join(expr_df_all, metadata, by = "orig.ident")  # to debug missing samples

# Rename column V1 to XIST 
colnames(expr_df_all)[1] <- "XIST"
expr_df_all$sample_id <- sub("^MB", "MDT-AP-", expr_df_all$orig.ident)

# Load the metadata TSV file
metadata <- read.delim("20241218_multiome_metadata.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# View first rows
head(metadata)

# Filter for primary Group 3 or Group 4 and the repeated records are for different histology
primary_G3_metadata <- subset(metadata,
                              primary_trio_info == "Primary" & (
                                rna_dx == "Group 3" |
                                  methyl_dx == "MB, G3"
                              ))
length(primary_G3_metadata$unique_id) # 187
length(unique(primary_G3_metadata$unique_id)) # 117

primary_G4_metadata <- subset(metadata,
                              primary_trio_info == "Primary" & (
                                rna_dx == "Group 4" |
                                  methyl_dx == "MB, G4"
                              ))
length(primary_G4_metadata$unique_id) # 377
length(unique(primary_G4_metadata$unique_id)) # 267


# Match: metadata$unique_id == expr_df_all$sample_id??
# Rename metadata column for clean merge
colnames(metadata)[colnames(metadata) == "unique_id"] <- "orig.ident"

# Merge expression with metadata
merged_df_1 <- merge(expr_df_all, metadata, by = "sample_id")
#merged_df_1$"sample_id"
#merged_df_1$"sex"
#merged_df_1$"XIST"
merged_df_2 <- merge(expr_df_all, metadata, by = "orig.ident")
#merged_df_2$orig.ident
#merged_df_2$sex
#merged_df_2$XIST
#MDT merge 

# Select the three columns
new_df1 <- merged_df_1[, c("sample_id", "sex", "XIST")]
new_df1$sample_id <- sub("^MDT-AP-", "MB", new_df1$sample_id)

new_df2 <- merged_df_2[, c("orig.ident", "sex", "XIST")]
# Rename 'orig.ident' to 'sample_id'
colnames(new_df2)[colnames(new_df2) == "orig.ident"] <- "sample_id"

merged_df_all = rbind(new_df1, new_df2)

# Inspect merged data 12 samples?
head(merged_df_all)
dim(merged_df_all)
#merged_df_all$primary_trio_info # all trio_info are primary in the merged data

unlabeled_merged_df <- merged_df_all[!complete.cases(merged_df_all), ]
merged_df <- merged_df_all[complete.cases(merged_df_all), ]

####################Gender inference using xgboost method###########################################################

##########################CROSS VALIDATION TRAINING################################################################################
# Clean and prepare
df <- merged_df[!is.na(merged_df$sex), ]
df$sex <- toupper(df$sex)
df$label <- ifelse(df$sex == "F", 1, 0)  # XGBoost needs numeric labels

# 10-fold stratified CV
folds <- createFolds(df$label, k = 10)


results <- data.frame(Fold = integer(), Accuracy = numeric(), Sensitivity = numeric(), Specificity = numeric())
all_preds <- data.frame()

for (i in 1:10) {
  test_idx <- folds[[i]]
  train_df <- df[-test_idx, ]
  test_df <- df[test_idx, ]
  
  dtrain <- xgb.DMatrix(data = as.matrix(train_df$XIST), label = train_df$label)
  dtest  <- xgb.DMatrix(data = as.matrix(test_df$XIST))
  
  # Train model
  #model <- xgboost(data = dtrain, nrounds = 50, objective = "binary:logistic", verbose = 0)
  pos_weight <- sum(train_df$label == 0) / sum(train_df$label == 1)
  model <- xgb.train(
    data = dtrain,
    nrounds = 100,
    objective = "binary:logistic",
    params = list(scale_pos_weight = pos_weight),
    verbose = 0
  )
  # Predict probabilities
  probs <- predict(model, dtest)
  preds <- ifelse(probs > 0.2, 1, 0)
  
  # Add predictions
  test_df$predicted <- preds
  test_df$prob <- probs
  test_df$Fold <- i
  all_preds <- rbind(all_preds, test_df)
  
  # Evaluate
  TP <- sum(test_df$predicted == 0 & test_df$label == 0)
  TN <- sum(test_df$predicted == 1 & test_df$label == 1)
  FP <- sum(test_df$predicted == 0 & test_df$label == 1)
  FN <- sum(test_df$predicted == 1 & test_df$label == 0)
  
  eps <- 1e-10
  acc <- (TP + TN) / (TP + TN + FP + FN + eps)
  sens <- if ((TP + FN) == 0) 1 else TP / (TP + FN + eps) 
  spec <- if ((TN + FP) == 0) 1 else TN / (TN + FP + eps)
  
  results <- rbind(results, data.frame(Fold = i, Accuracy = acc, Sensitivity = sens, Specificity = spec))
}

#roc_obj <- roc(response = test_df$label, predictor = probs)

# Summary
cat("📊 5-Fold Cross-Validation Summary:\n")
print(round(results, 3))
cat("\nAverage Accuracy:", round(mean(results$Accuracy), 3), "\n")
cat("Average Sensitivity (Male):", round(mean(results$Sensitivity), 3), "\n")
cat("Average Specificity (Female):", round(mean(results$Specificity), 3), "\n")
#Average Accuracy: 0.967
#Average Sensitivity (Male): 1 
#Average Specificity (Female): 0.9 

