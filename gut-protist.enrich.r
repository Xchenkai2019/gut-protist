##The primary resource code refers to https://github.com/YongxinLiu/MicrobiomeStatPlot/
#step1 Loaded package
library(edgeR)
library(tidyverse)

#step2 Enrichment analyses
setwd("F:/³ÌÐò/github/gut-protist/")
otu_relative <- apply(otu_table, 2, function(x){x/sum(x)})
  if (missing(threshold))
    threshold = 0.0005
  idx <- rowSums(otu_relative > threshold) >= 1
  otu <- as.data.frame(otu_table[idx, ])
  otu_relative <- as.data.frame(otu_relative[idx, ])

  ## construct a DGEList
  dge_list <- DGEList(counts = otu, group = design$group)

  ## Remove the lower abundance(In this case, it's usually useless)
  keep <- rowSums(dge_list$counts) >= 0
  dge_keep <- dge_list[keep, ,keep.lib.sizes = FALSE]

  ## scale the raw library sizes dgelist
  dge <- calcNormFactors(dge_keep)

  ## set the design_mat
  design_mat <- model.matrix(~0 + dge$samples$group)
  colnames(design_mat) <- gsub("([dge$samples$group])", 
                               "", colnames(design_mat))
  ## fit the GLM
  GLMC = estimateGLMCommonDisp(dge, design_mat)
  GLMT = estimateGLMTagwiseDisp(GLMC, design_mat)
  fit = glmFit(GLMT, design_mat)

  ## conducts likelihood ratio tests
  contrast_mat <- makeContrasts(
    contrasts = contrasts, 
    levels = colnames(design_mat))

  ## Fit a negative binomial generalized log-linear model to the read counts
  lrt = glmLRT(fit, contrast = contrast_mat)

  ## Multiple Testing
  de_lrt <- decideTestsDGE(lrt, adjust.method = "fdr", p.value = 0.05)

  ## extract values
  data <- lrt$table
  data$sign_level <- de_lrt@.Data

  ## enrichment status
  data$enrichment <- as.factor(ifelse(data$sign_level == 1, "enriched", ifelse(data$sign_level == -1, "depleted","nosig")))
  ## get the out's names
  data$otus <- rownames(data)

  ## negative logarithmes transformation
  data$neglog_p = -log(data$PValue)

  ## remove low foldchange
  idx <- data$logFC < 0
  data$neglog_p[idx] <- 0

  ## reorder OTUs according to taxonomy(Phylum, Class, Order)
  taxonomy <- taxonomy_table[order(taxonomy_table[, 3], 
                                   taxonomy_table[, 4], 
                                   taxonomy_table[, 5]), ]

  idx <- taxonomy[, 1] %in% data$otu
  taxonomy <- taxonomy[idx, ]

  idx <- match(taxonomy[, 1], data$otu)
  data <- data[idx, ]

  data$classification  <- taxonomy[, 3]
  data$otu <- factor(data$otu, levels = data$otu)

  ## calculating the relative abundance 
  ra <- apply(otu_relative, 1, mean)
  data$ra <- ra[match(data$otu, names(ra))]

  return(data)
} 

#step3 Set threshold line
threshold_line <- function(data, Pvalue){
  if (missing(Pvalue)){
    FDR <- min(data$neglog_p[data$enrichment == "enriched"])}
  else {
    FDR <- -log(Pvalue)}
 return(FDR) 
}

#step4 Data next processing
otu <- read.delim("otutab.txt", header = TRUE, 
                  sep = "\t", row.names = 1,
                  stringsAsFactors = FALSE)

taxonomy <- read.delim("taxonomy.txt", header = TRUE, 
                       sep = "\t", 
                       stringsAsFactors = FALSE)

design <- data.frame(
  samples = colnames(otu),
  group = rep(c("KO", "WT", "OE" ),each = 6))

## Running enrichment_analyses

data <- enrichment_analyses(otu_table = otu, design = design, taxonomy_table = taxonomy, contrasts = "KO-WT")

FDR <- threshold_line(data = data)

#step5 Figure
ggplot(data, aes(x = otu, y = neglog_p, color = classification, size = ra, shape = enrichment)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = FDR, linetype = 2, color = "red") +
  scale_shape_manual(values=c("enriched" = 17, 
                                "depleted" = 25, 
                                "nosig" = 20))+
  labs(x="OTUs", y = expression(~-log[10](P))) +
  guides(colour = "none") +
  basic_theme
