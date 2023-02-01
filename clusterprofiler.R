library(rio)
library(clusterProfiler)

setwd("")

# pathway databases

files <- list.files()
pathways <- lapply(files, function(file) {
  
  data <- read.gmt(file)
  name <- gsub(".txt", "", file)
  data$term <- paste(name, data$term, sep = "$")
  data
  
})
pathways <- do.call(rbind, pathways)

# import data

data <- import("TableS5.tsv")
gene <- data$Gene[which(data$adj.P.value < 0.05)] # significant genes
universe <- data$Gene # all genes in your data

# enrichment

res <- as.data.frame(enricher(gene, universe = universe, TERM2GENE = pathways, minGSSize = 1, maxGSSize = 1000000))

# number of upregulated and downregulated genes per pathway

updown <- lapply(res$geneID, function(genes) {
  
  genes <- unlist(strsplit(genes, split = "/"))
  up <- genes[which(genes %in% sapply(data$Gene[which(data$adj.P.value < 0.05 & data$Estimate > 0)], function(x) strsplit(x, split = ",")[1]))]
  down <- genes[which(genes %in% sapply(data$Gene[which(data$adj.P.value < 0.05 & data$Estimate < 0)], function(x) strsplit(x, split = ",")[1]))]
  data.frame(up = length(up), down = length(down))
  
})
updown <- do.call(rbind, updown)
res$up <- updown$up
res$down <- updown$down

# calculate jaccard index

jaccard <- lapply(res$ID, function(term) {
  
  genes <- unlist(strsplit(res$geneID[which(res$ID == term)], split = "/"))
  
  lapply(res$ID, function(term2) {
    
    genes2 <- unlist(strsplit(res$geneID[which(res$ID == term2)], split = "/"))
    jaccard <- length(intersect(genes, genes2)) / length(union(genes, genes2))
    data.frame(term, term2, jaccard)
    
  })
  
})
jaccard <- do.call(rbind, do.call(rbind, jaccard))

jaccard <- jaccard[which(jaccard$jaccard > 0.7), ] # threshold

# find groups

groups <- graph_from_data_frame(jaccard[, 1:2])
groups <- cluster_edge_betweenness(groups)
groups <- data.frame(term = groups$names, group = groups$membership)

# clean results

res$Group <- groups$group[match(res$ID, groups$term)]
split <- strsplit(res$ID, split = "$", fixed = T)
split <- do.call(rbind, split)
res$Term <- split[, 2]
res$Database <- split[, 1]
split <- strsplit(res$BgRatio, split = "/")
split <- do.call(rbind, split)
res$Overlap <- paste0(res$Count, "/", split[, 1])

res <- res[, c("Term", "Database", "Overlap", "pvalue", "p.adjust", "geneID", "up", "down", "Group")]
colnames(res) <- c("Term", "Database", "Overlap", "P.value", "Adjusted.P.value", "Genes", "Upregulated", "Downregulated", "Group")