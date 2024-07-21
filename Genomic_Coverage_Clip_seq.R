library(RCAS)
library(dplyr)
library(ggplot2)
library(knitr)
library(ggseqlogo)

# Function to run the analysis using mouse genome (This can be modified as per the user's input to other genomes)
run_analysis <- function(queryFilePath, gffFilePath, genomeVersion = 'mm10', sampleN = 10000, flankSize = 1000, quiet = FALSE, motifAnalysis = TRUE, goAnalysis = TRUE, printProcessedTables = FALSE) {
  
  # Importing the query regions in Bed format
  queryRegions1 <- importBed(queryFilePath, sampleN = sampleN)
  
  # Importing GFF file
  gff1 <- importGtf(filePath = gffFilePath)
  
  # Performing overlap analysis
  overlaps <- as.data.table(queryGff(queryRegions = queryRegions1, gffData = gff1))
  biotype_col <- grep('gene_type', colnames(overlaps), value = TRUE)
  df <- overlaps[, length(unique(queryIndex)), by = biotype_col]
  colnames(df) <- c("feature", "count")
  df$percent <- round(df$count / length(queryRegions1) * 100, 1)
  df <- df[order(count, decreasing = TRUE)]
  
  # Plotting percentage overlap
  p1 <- ggplot(df, aes(x = reorder(feature, -percent), y = percent)) + 
    geom_bar(stat = 'identity', aes(fill = feature)) + 
    geom_label(aes(y = percent + 0.5), label = df$count) + 
    labs(x = 'transcript feature', y = paste0('percent overlap (n = ', length(queryRegions1), ')')) + 
    theme_bw(base_size = 14) + 
    theme(axis.text.x = element_text(angle = 90))
  print(p1)
  
  # Extending annotation feature space
  txdbFeatures <- getTxdbFeaturesFromGRanges(gff1)
  summary <- summarizeQueryRegions(queryRegions = queryRegions1, txdbFeatures = txdbFeatures)
  df <- data.frame(summary)
  df$percent <- round((df$count / length(queryRegions1)), 3) * 100
  df$feature <- rownames(df)
  
  # Plotting annotation feature space
  p2 <- ggplot(df, aes(x = reorder(feature, -percent), y = percent)) + 
    geom_bar(stat = 'identity', aes(fill = feature)) + 
    geom_label(aes(y = percent + 3), label = df$count) + 
    labs(x = 'transcript feature', y = paste0('percent overlap (n = ', length(queryRegions1), ')')) + 
    theme_bw(base_size = 14) + 
    theme(axis.text.x = element_text(angle = 90))
  print(p2)
  
  setwd(dirname(queryFilePath))
  ggsave("percent.pdf", width = 7, height = 5)
  
  # Overlap counts between query regions and genes
  dt <- getTargetedGenesTable(queryRegions = queryRegions1, txdbFeatures = txdbFeatures)
  dt <- dt[order(transcripts, decreasing = TRUE)]
  if (printProcessedTables) {
    knitr::kable(dt[1:10,])
  }
  
  # Coverage profiles at feature boundaries
  cvgF <- getFeatureBoundaryCoverage(queryRegions = queryRegions1, featureCoords = txdbFeatures$transcripts, flankSize = flankSize, boundaryType = 'fiveprime', sampleN = sampleN)
  cvgT <- getFeatureBoundaryCoverage(queryRegions = queryRegions1, featureCoords = txdbFeatures$transcripts, flankSize = flankSize, boundaryType = 'threeprime', sampleN = sampleN)
  
  cvgF$boundary <- 'fiveprime'
  cvgT$boundary <- 'threeprime'
  df <- rbind(cvgF, cvgT)
  
  p3 <- ggplot(df, aes(x = bases, y = meanCoverage)) + 
    geom_line(aes(ymin = meanCoverage - standardError * 1.96, ymax = meanCoverage + standardError * 1.96)) + 
    geom_line(color = 'black') + 
    facet_grid(~ boundary) + 
    theme_bw(base_size = 14)
  print(p3)
  
  # Coverage profile for all transcript features
  cvgList <- calculateCoverageProfileList(queryRegions = queryRegions1, targetRegionsList = txdbFeatures, sampleN = sampleN)
  p4 <- ggplot(cvgList, aes(x = bins, y = meanCoverage)) + 
    geom_ribbon(aes(ymin = meanCoverage - standardError * 1.96, ymax = meanCoverage + standardError * 1.96), fill = 'lightgreen') + 
    geom_line(color = 'black') + 
    theme_bw(base_size = 14) +
    facet_wrap(~ feature, ncol = 3)
  print(p4)
  
  cvg_sublist <- subset(cvgList, feature %in% c('cds', 'threeUTRs', 'fiveUTRs'))
  cvg_sublist$feature[cvg_sublist$feature == "threeUTRs"] <- '3 UTR'
  cvg_sublist$feature[cvg_sublist$feature == "fiveUTRs"] <- '5 UTR'
  cvg_sublist$feature[cvg_sublist$feature == "cds"] <- 'CDS'
  p5 <- ggplot(cvg_sublist, aes(x = bins, y = meanCoverage)) + 
    geom_line(aes(ymin = 0, ymax = meanCoverage + standardError * 1.96)) + 
    geom_line(color = 'black') + 
    theme_bw(base_size = 14) +
    facet_wrap(~ feature, ncol = 3)
  print(p5)
  
  ggsave("CDS_5and3UTRs.pdf", width = 9, height = 5)
  write.table(cvg_sublist, "cvg_features.txt", sep = '\t')
  write.table(cvgList, "cvg_all_transcripts.txt", sep = '\t')
  
  # Motif analysis
  if (motifAnalysis) {
    motifResults1 <- runMotifDiscovery(queryRegions = queryRegions1, resizeN = 15, sampleN = sampleN, genomeVersion = genomeVersion, motifWidth = 6, motifN = 2, nCores = 1)
    write.table(motifResults1, "motif.txt")
    ggseqlogo::ggseqlogo(motifResults1$matches_query)
    ggsave("motif.pdf", width = 7, height = 5)
    summary <- getMotifSummaryTable(motifResults1)
    if (printProcessedTables) {
      knitr::kable(summary)
    }
  }
  
  # Functional enrichment analysis
  if (goAnalysis) {
    targetedGenes <- unique(overlaps$gene_id)
    res <- findEnrichedFunctions(targetGenes = targetedGenes, species = 'mouse')
    res <- res[order(res$p_value),]
    resGO <- res[grep('GO:BP', res$source),]
    if (printProcessedTables) {
      knitr::kable(subset(resGO[1:10,], select = c('p_value', 'term_name', 'source')))
    }
  }
}

# Example usage
run_analysis(queryFilePath = 'input.BED', gffFilePath = 'annotation.gtf', genomeVersion = 'mm10', quiet = TRUE, motifAnalysis = TRUE, goAnalysis = TRUE, printProcessedTables = TRUE)
