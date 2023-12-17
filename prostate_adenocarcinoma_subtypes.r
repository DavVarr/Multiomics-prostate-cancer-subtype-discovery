if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("curatedTCGAData")
BiocManager::install("TCGAutils")
BiocManager::install("TCGAbiolinks")
BiocManager::install("graph")

install.packages("SNFtool")
install.packages("caret")
install.packages("cluster")
install.packages("mclustcomp")
install.packages("ggplot2")
install.packages("mclust")
install.packages("NetPreProc")

library("curatedTCGAData")
library("TCGAutils")
library("TCGAbiolinks")
library("SNFtool")
library("caret")
library("cluster")
library("mclustcomp")
library("NetPreProc")
library("ggplot2")
library("mclust")

# The code follows the steps proposed by the project guidelines
# for this reason (since I chose to do the project individually and not in a group) steps 7, 8.d, 9 and 10 aren't mentioned

# Step 1) Download PRAD dataset from TCGA


assays <- c("miRNASeqGene", "RNASeq2Gene", "RPPAArray")
mo <- curatedTCGAData(diseaseCode = "PRAD",
                      assays = assays,
                      version = "2.0.1", dry.run = FALSE)


mo


# Step 2) Pre-processing steps

# Consider only primary solid tumors
# They are identified by the code "01" in the sample part of the barcode
primary <- TCGAutils::TCGAsampleSelect(colnames(mo), sampleCodes = "01")
mo <- mo[, primary, ]

# Check for replicates
check_rep <- anyReplicated(mo)
print(check_rep)

# Remove FFPE samples
# the information is store in the phenotype dataframe, which is accessible using colData()
no_ffpe <- which(as.data.frame(colData(mo))$patient.samples.sample.is_ffpe == "no")
mo <- mo[, no_ffpe, ]

# Obtain samples having all the considered omics
complete <- intersectColumns(mo)

# Extract assays in list
complete <- assays(complete)

# Transpose the matrices to obtain samples on rows and features on columns
complete <- lapply(complete, FUN = t)

# Check for NA values
for (i in 1:length(complete)){
  message(paste("omics",i,"missing:",sum(is.na(complete[[i]])),"elements"))
}
# Only the proteomics data contains NA values, remove features that contains them
complete[[3]] <- complete[[3]][, colSums(is.na(complete[[3]])) == 0]

# Remove features with near zero variance and retain top 100 features having higher variance:
nf <- 100
for (i in 1:length(complete)) {
  idx <- caret::nearZeroVar(complete[[i]])
  message(paste("Removed", length(idx), "features from", names(complete)[i]))
  if (length(idx) != 0) {
    complete[[i]] <- complete[[i]][, -idx]
  }
  
  if (ncol(complete[[i]]) <= nf) next
  
  vars <- apply(complete[[i]], 2, var)
  idx <- sort(vars, index.return = TRUE, decreasing = TRUE)$ix
  
  complete[[i]] <- complete[[i]][, idx[1:nf]]
}

# Perform features standardization using z-score:
zscore <- function(data){
  
  zscore_vec <- function(x) { return ((x - mean(x)) / sd(x))}
  data <- apply(data, 2, zscore_vec)
  
  
  return(data)
}

complete <- lapply(complete, zscore) 

# Clean barcodes retaining only "Project-TSS-Participant":
for (v in 1:length(complete)) {
  rownames(complete[[v]]) <- substr(rownames(complete[[v]]), 1, 12)
}


# Step 3) Download disease subtypes from TCGAbiolinks
subtypes <- as.data.frame(TCGAbiolinks::PanCancerAtlas_subtypes())
subtypes <- subtypes[subtypes$cancer.type == "PRAD", ]
# Retain only primary solid tumor
subtypes <- subtypes[TCGAutils::TCGAsampleSelect(subtypes$pan.samplesID, "01"), ]
# Select samples in common with omics data
sub_select <- substr(subtypes$pan.samplesID,1,12) %in% rownames(complete[[1]])
subtypes <- subtypes[sub_select, ]

#check for NA values in the iCluster subtype columnn
sum(is.na(subtypes$Subtype_Integrative))

# Set subtypes row index equal to the first 12 chars of the barcode
rownames(subtypes) <- substr(subtypes$pan.samplesID, 1, 12);

# The subtypes samples are less than the omics ones
print(length(rownames(complete[[1]])) > length(rownames(subtypes)))

# So, consider only the samples that are both in subtypes and in complete
for (v in 1:length(complete)) {
  complete[[v]] <- complete[[v]][rownames(subtypes),] 
}


# Step 4) Check if the patients in the multi-omics dataset and the subtypes dataset are in the same order
same_order <- TRUE
for (v in 1:length(complete)) {
  same_order <- same_order & identical(rownames(complete[[v]]), rownames(subtypes))
}

# Print the result
if (same_order) {
  message("The patients in the multi-omics dataset and the subtypes dataset are in the same order.")
} else {
  warning("The patients in the multi-omics dataset and the subtypes dataset are not in the same order.")
}

# Step 5) Similarity Network Fusion with the scaled exponential euclidean distance
W_list <- list();
# Compute the scaled euclidean distance for every omics as a similarity measure
for(i in 1:length(complete)){
  Dist <- (dist2(as.matrix(complete[[i]]), as.matrix(complete[[i]])))^(1/2)
  # Neigbour = 20  and  parameter alpha = 0.5 by default
  W_list[[i]] <- affinityMatrix(Dist)
  
}

# Integration using SNF 
# k is the number of neighbours considered to compute the local similarity matrix
# t is the number of iterations
W_snf <- SNF(W_list, K=20, t=20)


# Step 6) Integrate the similarity matrices from each data source using an average of the matrices
W_avg <- Reduce("+", W_list) / length(W_list)



# Step 8) Perform disease subtype discovery using PAM algorithm, which requires
# to convert the similarity matrices to normalized distance matrices.
k <- length(unique(subtypes$Subtype_Integrative)) #number of clusters = number of iCluster disease subtypes (3)

# Step 8.a) Run PAM algorithm on each omic data, individually
pam_singles_res <- list()
for (i in 1:length(W_list)){
  dist <- 1 - NetPreProc::Max.Min.norm(W_list[[i]])
  D <- as.dist(dist)
  pam_singles_res[[i]] <- pam(D,k=k)
}

# Step 8.b) Run PAM on the matrix obtained by averaging the 3 omics data sources 

W_avg_dist <- 1 - NetPreProc::Max.Min.norm(W_avg)
pam_avg_res <- pam(as.dist(W_avg_dist),k=k)

# Step 8.c) Run PAM on the matrix obtained with SNF
W_snf_dist <- 1 - NetPreProc::Max.Min.norm(W_snf)
pam_snf_res <- pam(as.dist(W_snf_dist),k=k)

# Step 11) Compare the clusterings obtained by each considered approach and show the results
# by using rand index, adjusted rand index and normalized mutual information
labels <- as.numeric(factor(subtypes$Subtype_Integrative, levels=unique(subtypes$Subtype_Integrative)))
types <- c("rand", "adjrand", "nmi1")

metrics_pam_avg <- mclustcomp(pam_avg_res$clustering, labels, types=types)
metrics_pam_snf <- mclustcomp(pam_snf_res$clustering, labels, types=types)
metrics_pam_miRNA <- mclustcomp(pam_singles_res[[1]]$clustering, labels, types=types)
metrics_pam_mRNA <- mclustcomp(pam_singles_res[[2]]$clustering, labels, types=types)
metrics_pam_protein <- mclustcomp(pam_singles_res[[3]]$clustering, labels, types=types)
# Create a data frame to store the evaluation metrics for each approach
metrics_df <- data.frame(
  Approach = c("SNF integration", "Avg integration", "miRNA", "mRNA", "proteins"),
  AdjRand = c(
    metrics_pam_snf$scores[which(metrics_pam_snf$types == "adjrand")],
    metrics_pam_avg$scores[which(metrics_pam_avg$types == "adjrand")],
    metrics_pam_miRNA$scores[which(metrics_pam_miRNA$types == "adjrand")],
    metrics_pam_mRNA$scores[which(metrics_pam_mRNA$types == "adjrand")],
    metrics_pam_protein$scores[which(metrics_pam_protein$types == "adjrand")]
  ),
  NMI = c(
    metrics_pam_snf$scores[which(metrics_pam_snf$types == "nmi1")],
    metrics_pam_avg$scores[which(metrics_pam_avg$types == "nmi1")],
    metrics_pam_miRNA$scores[which(metrics_pam_miRNA$types == "nmi1")],
    metrics_pam_mRNA$scores[which(metrics_pam_mRNA$types == "nmi1")],
    metrics_pam_protein$scores[which(metrics_pam_protein$types == "nmi1")]
  ),
  Rand = c(
    metrics_pam_snf$scores[which(metrics_pam_snf$types == "rand")],
    metrics_pam_avg$scores[which(metrics_pam_avg$types == "rand")],
    metrics_pam_miRNA$scores[which(metrics_pam_miRNA$types == "rand")],
    metrics_pam_mRNA$scores[which(metrics_pam_mRNA$types == "rand")],
    metrics_pam_protein$scores[which(metrics_pam_protein$types == "rand")]
  )
)

# Print the table
print(metrics_df)

# Convert the data frame to long format
metrics_df_long <- tidyr::gather(metrics_df, "Metric", "Score", -Approach)

# Create the bar plot
ggplot(metrics_df_long, aes(x = Approach, y = Score, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Approach", y = "Score", fill = "Metric") +
  ggtitle("Comparison of evaluation metrics for each approach")
