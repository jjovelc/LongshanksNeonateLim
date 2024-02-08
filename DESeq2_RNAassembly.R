# Differential expression analysis
# Project: Long shank neonate lim
#
# By Juan Jovel (use it at your own risk)
# (juan.jovel@ucalgary.ca)
# 
# Last revision: Feb. 06, 2023

library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(gplots)
library(pheatmap)
library(biomaRt)
library(EnhancedVolcano)

rm(list = ls())

setwd('/Users/juanjovel/jj/data_analysis/coltonUnger/DE_analysis/mm10_plus_RNAspades_assembly')
data     <- read.table("all_samples_counts.tsv", sep = '\t', header = T, row.names = 1, stringsAsFactors = T)
metadata <- read.table("metadata.tsv", sep = '\t', header = T, row.names = 1)
metadata$LineTissue <- paste(metadata$Line, metadata$Tissue, sep = '_')

data_round  <- round(data, digits = 0)
data_10up   <- subset(data_round, rowMeans(data_round) >= 10)

# Because in an initial inspection of the PCA plot including all samples
# it was determined that replicate 4 in the CTL line was clustering with
# the LS1 line, it will be removed.

# Logical vector to retain rows not matching the criteria
keep_rows         <- !(metadata$Line == "CTL" & metadata$Replicate == "r4")
filtered_metadata <- metadata[keep_rows, ]
filtered_data     <- data_10up[,keep_rows]

data_matrix <- as.matrix(filtered_data) 


# Import data matrix into a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = data_matrix, colData = filtered_metadata,
                              design =~ LineTissue + Weight)

# Calculate size factors
dds <- estimateSizeFactors(dds)

# Apply a regularized logarithmic transformation
rld <- rlogTransformation(dds)

# Hierarchical clustering
makeHCheatmap <- function(rld, metadata, prefix){
  dist = dist(t(assay(rld)))
  dist_matrix <- as.matrix(dist)
  row.names(dist_matrix) <- metadata$Line
  # Retrieve the name of the dataframe
  df_name <- deparse(substitute(metadata))
  
  # Use gsub to create the new name
  file_name <- paste(prefix, 'HC_plot.png', sep = '_')
  png(file_name)
  pheatmap(dist_matrix, 
           color=colorRampPalette(brewer.pal(n=9,name = "BrBG"))(255), 
           clustering_distance_cols = dist, clustering_distance_rows = dist
  )
  dev.off()
}

makeHCheatmap(rld, filtered_metadata, "longShank")


# PCA plots
makePCA <- function(rld, metadata, prefix, color1, color2, color3) {
  df_name <- deparse(substitute(metadata))
  
  # Use gsub to create the new name
  file_name <- paste(prefix, "PCA_plot.png", sep = '_')
  
  data <- plotPCA(rld, intgroup = c("Line"), returnData = TRUE)
  percentVar <- round(100 * attr(data, "percentVar"))
  
  # Define custom color palette for groups
  custom_colors <- c(color1, color2, color3)
  
  PCA_EucDist <- ggplot(data, aes(x = PC1, y = PC2, color = group)) +
    xlab(paste("PC1 :", percentVar[1], "% variance")) +
    ylab(paste("PC2 :", percentVar[2], "% variance")) +
    ggtitle(paste(prefix, "PCA", sep = ' ')) +
    geom_point(aes(fill = Line), shape = 21, size = 15) + # Add filled points with black border
    scale_color_manual(values = custom_colors) + # Define point colors
    scale_fill_manual(values = custom_colors) +  # Define fill colors
    theme_bw() +
    theme(legend.position = "right")  # Adjust legend position
  
  # Add encircling ellipses to each group
  PCA_EucDist <- PCA_EucDist +
    stat_ellipse(aes(fill = group), geom = "polygon", level = 0.95, alpha = 0.2, 
                 show.legend = FALSE) +
  
  # Save the plot as a PNG file
  png(file_name, width = 1200, height = 1200)  # Specify width and height as needed
  print((PCA_EucDist) +
          geom_text(aes(label=rownames(metadata)), cex = 2, hjust=0.5, vjust=2, color="black"))
  dev.off()
}

makePCA(rld, filtered_metadata, "longShank", "darkgray", "dodgerblue", "goldenrod1")  



##### ANNOTATION OF DE FEATURES #####
library(biomaRt)

# Create remote connection
my_mart <- useEnsembl('ensembl', dataset = "mmusculus_gene_ensembl")


# Make list of attributes to retrieve
attributes <- c("ensembl_transcript_id",
                "ensembl_gene_id",
                "entrezgene_id",
                "external_gene_name",
                "wikigene_description",
                "name_1006",
                "definition_1006",
                "namespace_1003")

# Pull attributes table from Ensembl
# Function to extract records
pullRecords <- function(attributes, mart, filter_values){
  records <- getBM(attributes = attributes, filters = "ensembl_transcript_id", 
                   values = filter_values, mart = mart)
  
  return(records)
}


# Extract first hit only
getFirstMatch <- function(records, transcripts) {
  first_annotation_df <- data.frame()
  for (transcript in transcripts) {
    hits <- which(records$ensembl_transcript_id == transcript)
    if (length(hits) > 0) {
      first <- hits[1]
      first_annotation <- records[first,]
      first_annotation_df <- rbind(first_annotation_df, first_annotation)
    } else {
      first_annotation <- c(transcript, "_", "-", "-", "-", "-", "-", "-")
      first_annotation_df <- rbind(first_annotation_df, first_annotation)
    }
  }
  return(first_annotation_df)
}

renameColumns <- function(df){
  colnames(df)[6] <- "GO_group"
  colnames(df)[7] <- "GO_definition"
  colnames(df)[8] <- "Ontology"
  return(df)
}

makeVolcanoPlot <- function(df, vp_file){
  keyvals <- ifelse(
    df$log2FoldChange < -1, 'forestgreen',
    ifelse(df$log2FoldChange > 1, 'firebrick1',
           'dodgerblue1'))
  
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'firebrick1'] <- 'High'
  names(keyvals)[keyvals == 'dodgerblue1'] <- 'small FC'
  names(keyvals)[keyvals == 'forestgreen'] <- 'Low'
  
  evp <- EnhancedVolcano(df,
                         lab = df$external_gene_name,
                         x = 'log2FoldChange',
                         y = 'pvalue',
                         pCutoff = 0.01,
                         FCcutoff = 1,
                         colCustom = keyvals,
                         title = NULL,
                         subtitle = NULL,
                         colAlpha = 0.85,
                         shape = 20,
                         pointSize = 2,
                         labSize = 3)
  
  png(vp_file)
  print(evp)
  dev.off()
}



# Differential expression analysis
dds <- DESeq(dds, fitType = 'local')
# Running DESeq with the LRT:
# dds <- DESeq(dds, test="LRT", reduced=~Replicate)
# dds <- DESeq(dds, test="LRT", reduced=~1)
resultsNames(dds)
#res <- results(dds)

# Include normalized data in results table
allRes <- " LongShankNeonateLim_allRes.tsv"
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Transcript"
resdata <- cbind(geneName=rownames(resdata), resdata)
write.table(resdata, file=allRes, sep="\t", quote = F, row.names = F)


# Let's extract relevant contrasts from the DE analysis
contrasts <- list(
  control_ls1_fem = c(group = "LineTissue", "LS1_Fem", "CTL_Fem"),
  control_ls1_tib = c(group = "LineTissue", "LS1_Tib", "CTL_Tib"),
  control_ls1_rad = c(group = "LineTissue", "LS1_Rad", "CTL_Rad"),
  control_ls2_fem = c(group = "LineTissue", "LS2_Fem", "CTL_Fem"),
  control_ls2_tib = c(group = "LineTissue", "LS2_Tib", "CTL_Tib"),
  control_ls2_rad = c(group = "LineTissue", "LS2_Rad", "CTL_Rad")
)

for (contrast_name in names(contrasts)) {
  print(contrast_name)
  contrast_value <- contrasts[[contrast_name]]  # Use double square brackets to access the inner list
  result <- results(dds, contrast = contrast_value)
  sign_results <- subset(result, padj < 0.01)
  print("Significant results extracted...")
  sign_results <- cbind(transcript=rownames(sign_results), sign_results)
  if (nrow(sign_results)){
    transcripts       <- sign_results$transcript
    transcripts_clean <- gsub("\\.\\d+$", "", transcripts) # remove transc version
    print("Extract and clean transcripts names...")
    bt_transc_annotations <- pullRecords(attributes, my_mart, transcripts_clean)
    print("Annotation pulled down from Emsembl...")
    annotation_df <- getFirstMatch(bt_transc_annotations, transcripts_clean)
    print("First annotation record extracted for each hit...")
    annotation_df <- renameColumns(annotation_df)
    print("Rename GO columns...")
    tail(sign_results)
    tail(annotation_df)
    sign_results_annotated <- cbind(sign_results, annotation_df)
    print("Pasted significant results with annotation df...")
    file_name <- paste(contrast_name, "q0.01_annotated.tsv", sep = '_')
    write.table(sign_results_annotated, file_name, quote = F, sep = '\t', row.names = F)
    print("Results saved to a file...")
    vp_file <- paste(contrast_name, "volcanoPlot.png", sep = '_')
    makeVolcanoPlot(sign_results_annotated, vp_file)
    print("Volcano plot created...")
    print(paste("Number of significant transcripts deregulated", 
              nrow(sign_results_annotated) - 1))
    print("______________________")
  } else {
    print("No significant features were found...")
  }
}



