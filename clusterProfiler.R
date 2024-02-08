library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(gplots)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)

rm(list = ls())

setwd('/Users/juanjovel/jj/data_analysis/coltonUnger/DE_analysis/model_group')
data_file = 'all_samples_counts.tsv'
metadata_file = 'metadata.tsv'

data = read.table(data_file, sep = '\t', header = T, row.names = 1)
metadata = read.table(metadata_file, sep = '\t', header = T, row.names = 1)

data_round  <- round(data, digits = 0)
data_10up   <- subset(data_round, rowMeans(data_round) >= 10)

# Logical vector to retain rows not matching the criteria
keep_rows         <- !(metadata$Line == "CTL" & metadata$Replicate == "r4")
filtered_metadata <- metadata[keep_rows, ]
filtered_data     <- data_10up[,keep_rows]

data_matrix <- as.matrix(filtered_data) 

filtered_metadata$LineTissue <- paste0(filtered_metadata$Line, "_", filtered_metadata$Tissue)

# Import data matrix into a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = data_matrix, colData = filtered_metadata,
                              design =~ LineTissue + Weight)

# Calculate size factor (library size) for normalization
dds = estimateSizeFactors(dds)

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


# Extract results table and print it to a file
getAndPrintResults <- function(go_object, file_name){
  my_res <- go_object@result
  write.table(my_res, file_name,
              quote = F, row.names = F,
              sep = '\t')
}

##### DIFFERENTIAL EXPRESSION ANALYSIS ####
dds <- DESeq(dds, fitType = 'local')
  #####

# Let's extract relevant contrasts from the DE analysis
contrasts <- list(
  control_ls1_fem = c(group = "LineTissue", "LS1_Fem", "CTL_Fem"),
  control_ls1_tib = c(group = "LineTissue", "LS1_Tib", "CTL_Tib"),
  control_ls1_rad = c(group = "LineTissue", "LS1_Rad", "CTL_Rad"),
  control_ls2_fem = c(group = "LineTissue", "LS2_Fem", "CTL_Fem"),
  control_ls2_tib = c(group = "LineTissue", "LS2_Tib", "CTL_Tib"),
  control_ls2_rad = c(group = "LineTissue", "LS2_Rad", "CTL_Rad")
)

for (contrast_name in names(contrasts)){
  print(contrast_name)
  contrast_value <- contrasts[[contrast_name]]  # Use double square brackets to access the inner list
  result <- results(dds, contrast = contrast_value)
  sign_results <- subset(result, padj < 0.01)
  print("Significant results extracted...")
  sign_results <- cbind(transcript=rownames(sign_results), sign_results)
  transcripts       <- sign_results$transcript
  transcripts_clean <- gsub("\\.\\d+$", "", transcripts) # remove transc version
  print("Extract and clean transcripts names...")
  bt_transc_annotations <- pullRecords(attributes, my_mart, transcripts_clean)
  print("Annotation pulled down from Emsembl...")
  annotation_df <- getFirstMatch(bt_transc_annotations, transcripts_clean)
  print("First annotation record extracted for each hit...")
  annotation_df <- renameColumns(annotation_df)
  print("Rename GO columns...")
  sign_results_annotated <- cbind(sign_results, annotation_df)
  print("Pasted significant results with annotation df...")
  file_name <- paste(contrast_name, "q0.01_annotated.tsv", sep = '_')
  write.table(sign_results_annotated, file_name, quote = F, sep = '\t', row.names = F)
  print("Results saved to a file...")
  
  ego_file_name <- paste(contrast_name, "upRegulated_GO_simp_results.tsv", sep = '_')
  
  ### GO enrich ###
  sign_data_ann_up      <- subset(sign_results_annotated, log2FoldChange > 0)
  entz_list_up          <- sign_data_ann_up$entrezgene_id
  entz_list_up_na       <- entz_list_up[!is.na(entz_list_up)]
  entz_list_up_na_dedup <- entz_list_up_na[!duplicated(entz_list_up_na)]
  
  # ************************* #
  ego <- clusterProfiler::enrichGO(
    entz_list_up_na_dedup,
    org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "ALL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.1,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE,
    pool = FALSE
  )
  
  ego2 <- simplify(ego, cutoff=0.7, by="p.adjust",
                   select_fun=min)
  
  getAndPrintResults(ego2, ego_file_name)
  ego_dotplot_file_name <- paste(contrast_name, "upRegulated_GO_simp_dotplot.png", sep = '_')
  png(ego_dotplot_file_name)
  p <- dotplot(ego2, showCategory = 20, font.size = 10)
  print(p)
  dev.off()
  ego_network_file_name <- paste(contrast_name, "upRegulated_GO_simp_neywork.png", sep = '_')
  png(ego_network_file_name)
  # Network plot
  p <- cnetplot(ego2, categorySize="p.adjust", foldChange=0.5)
  print(p)
  dev.off()
  
  print(paste("Number of significant transcripts deregulated", 
              nrow(sign_results_annotated) - 1))
  print("______________________")
}

