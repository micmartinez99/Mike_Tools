# The purpose of this script is to hold some helpful functions for future use

#----- Quick Regular Expressions

# When gene symbols are formatted by the convention I like (i.e., ENSROG0000128 - XYZ)
# This code will extract just the symbol portion
gsub("^[^-]+-(.*)$", "\\1", <whatever vector where the genes are stored>)

# To extract just the ENSEMBL portion
gsub("^(.*?)\\s-.*", "\\1", <whatever vector where the genes are stored>)

################################################################################

#----- Function to generate TPM values
# ! Gene lengths must be in terms of Kbp and NOT bp.
countToTPM <- function(count_df, length_df) {
  
  # Scale each count by its respective gene length
  rate <- count_df/length_df
  
  # Initialize an empty vector
  tpm <- c()
  
  # Iterate through the counts and scale by library size 
  for (i in 1:nrow(counts)) {
    tpm[i] <- rate[i]/sum(rate) * 1e6
  }
  
  return(tpm)
}

#----- Function to generate log2(TPM) values
log2TPMs <- function(TPMs) {
  
  # Add a small constant
  TPMs <- TPMs + 0.01
  
  # Log transform the tpm values
  logTPMs <- as.data.frame(t(apply(TPMs,1, log2)))
  
  return(logTPMs)
}

#----- Function to do exploratory PCA analysis
#! Must have a DESeq2 (dds) object that has been transformed with
# rlog() or vst()

# Function to plot variances to determine a suitable number of features
determineVarFeatures <- function(vsd) {
  
  # calculate gene expression level variance between samples
  var <- rev(rowVars(assay(vsd))[order(rowVars(assay(vsd)))])
  
  # reset plotting window to 1 row vs 1 columns
  par(mfrow=c(1,1))
  
  # plot variance for genes accross samples
  varPlot <- plot(var,
                  las = 1,
                  main="Sample gene expression variance",
                  xlab = "Gene", 
                  ylab = "Variance")
  
  # add vertical lines at specific gene number indexes
  abline(v=1000, col="red")
  abline(v=500, col="green")
  abline(v=250, col="blue")
  
  return(varPlot)
}

# Function to run PCA based on that number of variable features determined in
# determineVarFeatures() function

explorePCA <- function(vsd, nfeatures) {
  
  # Set a variable for the number of genes (features) to be used for PCA and clustering
  var_features_n <- nfeatures
  
  # Calculate the row variance
  rv <- rowVars(assay(vsd))
  
  # Order variance and select the rows (genes) with the most variance
  select <- order(rv, decreasing = TRUE)[1:var_features_n]
  
  # Subset vsd values for genes by top variance ranks
  vsd_sub <- assay(vsd)[select,]
  
  # Transpose the matrix
  vsd_sub <- t(vsd_sub)
  
  # Run principal component analysis
  pca <- prcomp(vsd_sub)
  
  # extract the variance explained by each PC
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  # subset for first 5 elemets
  percentVar <- percentVar[1:5]
  
  # add names to the percentVar vector
  names(percentVar) <- c("PC1", "PC2", "PC3", "PC4", "PC5")
  
  # Construct a data frame with PC loadings and add sample labels
  pca_df <- as.data.frame(pca$x)
  
  return(pca_df)
}











