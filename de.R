# Functions to obtain differentially expressed genes from limma
# This version supports inputting already normalized data to the function
# Daisy Dan
# last updated 03/16/2023

# Plotting top limma results on count. 
plot_norm_count_group <- function (gene, res, count, metadata, x_axis_label){
  if(all(tolower(colnames(count)) == rownames(metadata))){
    d <- metadata
    d$count <- count[rownames(count) == gene,]
    ggplot(d, aes(x = Group, y = count, color = Group)) +
      geom_boxplot() + 
      geom_point(position=position_jitter(w = 0.1,h = 0)) +
      # geom_text_repel(aes(label = rownames(d))) +
      theme_bw() +
      ggtitle(paste0(res[res$gene_id == gene,]$gene_name, " Expression by Group")) +
      theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.position="none") + ylab("Normalized Count") + xlab("") + 
      scale_x_discrete(labels=x_axis_label) + scale_colour_manual(values=c(rgb(84, 136, 227, maxColorValue = 255),
                                                                           rgb(211,111,104, maxColorValue = 255)))
    # + scale_color_grey()
  }else{
    message("check if rownames and colnames matches")
  }
}

plot_combined <- function (gene, res, count, metadata, x_axis_label){
  if(all(tolower(colnames(count)) == rownames(metadata))){
    d <- metadata
    d$count <- count[rownames(count) == gene,]
    ggplot(d, aes(x = group.condition, y = count, color = group.condition)) +
      geom_boxplot() + 
      geom_point() +
      # position=position_jitter(w = 0.1,h = 0)
      # geom_text_repel(aes(label = rownames(d))) +
      theme_bw() +
      ggtitle(paste0(res[res$gene_id == gene,]$gene_name, " Expression by Group")) +
      theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.position="none") + ylab("Normalized Count") + xlab("") + 
      scale_x_discrete(labels=x_axis_label) + 
      scale_colour_manual(values=c(rgb(84, 136, 227, maxColorValue = 255),
                                   rgb(0, 41, 227, maxColorValue = 255),
                                   rgb(211,111,104, maxColorValue = 255),
                                   rgb(237, 57, 28, maxColorValue = 255))) +
    # scale_color_grey() +
    # ggtitle(res[res$gene_id == gene,]$gene_name) +
    # theme(plot.title = element_text(hjust = 0.5))+
    geom_line(aes(group = individual),size=1, color='gray', alpha=0.6)
  }else{ 
    message("check if rownames and colnames matches")
  }
}

subset_de_analysis <- function(metadata_subset, design, v, coef, block = NULL, correlation = NULL,
                               contrasts = NULL, padj.cutoff = 0.05, ens, plot = FALSE, x_axis_label){
  design_subset <- design[rownames(design) %in% rownames(metadata_subset),]
  v_subset <- v[,match(rownames(metadata_subset), colnames(v))]
  res_subset <- limma_de_subset(v_subset, design_subset, coef = coef, ens = ens, block = block, correlation = correlation,contrasts = contrasts, adj ="bonferroni")
  sig_res_subset   <- res_subset[res_subset$baconPadj < padj.cutoff,]
  print(dim(sig_res_subset))
  if(plot){
    if(nrow(sig_res_subset) > 0) {
      print(lapply(sig_res_subset$gene_id, plot_norm_count_group, res_subset, v_subset$E, metadata_subset, x_axis_label))
    }else{
      print(lapply(res_subset$gene_id[1:10], plot_norm_count_group, res_subset, v_subset$E, metadata_subset, x_axis_label))
    }
  }
  return(list(res_subset, sig_res_subset))
}

filter_count <- function(d0, design, metadata_subset){
  keep <- filterByExpr(d0, design)
  d_keep <- d0[keep, keep.lib.sizes = FALSE]
  v <- voom(d_keep, design, plot = FALSE)
  v_subset <- v[,match(rownames(metadata_subset), colnames(v))]
  return(v_subset)
}

make_phist <- function(p, title){
  hist(p, main = paste0("Histogram of ", title), col = "violet", 
       xlab = "p-value")
}

make_qqplot <- function(p, title){
  qqman::qq(p, 
            main = paste0("QQ plot of ", title),
            sub = paste0("lambda: ", 
                         round(P_lambda(p), 3)))
}


limma_de_subset <- function(v, design, coef,
                            adj = "fdr", plot = TRUE, contrasts = NULL, block = NULL,
                            correlation = NULL, ens, weights = FALSE){
  require(limma)
  require(bacon)
  require(tidyverse)
  require(dplyr)
  require(ensembldb)
  
  # d0 <- DGEList(count)
  # keep <- filterByExpr(d0, design)
  # d_keep <- d0[keep, keep.lib.sizes = FALSE]
  # d_keep <- calcNormFactors(d_keep, method = "TMM")
  # 
  # if(weights == TRUE){
  #   v <- voomWithQualityWeights(d_keep, design, plot = plot)
  # } else {
  #   v <- voom(d_keep, design, plot = plot)
  # }
  # 
  fit <- lmFit(v, design, block = block, correlation = correlation)
  
  if(!is.null(contrasts)){
    fit <- contrasts.fit(fit, contrasts = contrasts)
  }
  
  fit <- eBayes(fit)
  
  # coef <- paste0("group", coef)
  
  # Bacon adjustment
  set.seed(707)
  bc <- bacon((fit$coef/fit$stdev.unscaled/fit$sigma)[, coef])
  print(bc)
  
  if(plot){
    plotSA(fit)
    plotMD(fit, coef = coef, main = coef)
  }
  print(summary(decideTests(fit)))
  
  results <- limma::topTable(fit, coef = coef, adjust.method=adj, number = Inf)
  
  results$df <- fit[rownames(results),]$df.residual
  results$baconT <- bacon::tstat(bc)[rownames(results), 1]
  results$baconP <- bacon::pval(bc)[rownames(results), 1]
  results$baconPadj <- p.adjust(results$baconP, method = adj)
  
  make_qqplot(results$P.Value, "limma p-values")
  make_qqplot(results$baconP, "bacon p-values")
  
  make_phist(results$P.Value, "limma p-values")
  make_phist(results$baconP, "bacon p-values")
  
  
  # Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- results %>%
    data.frame() %>%
    rownames_to_column(var="gene_id") %>% 
    as_tibble()
  
  # Create a gene-level dataframe
  annotations_ahb <- ensembldb::genes(ens, return.type = "data.frame")  %>%
    dplyr::select(gene_name, gene_id, description) %>%
    dplyr::filter(gene_id %in% res_tbl$gene_id)
  
  res_tbl <- merge(res_tbl, annotations_ahb, by = "gene_id", all = TRUE)
  
  # annot_gene <- ens %>%
  #   dplyr::filter(external_gene_name %in% res_tbl$gene_id)
  # annot_gene$gene_id <- annot_gene$external_gene_name
  # res_tbl <- merge(res_tbl, annot_gene, by = "gene_id", all = TRUE)
  
  res_tbl <- res_tbl[order(res_tbl$baconPadj),]
  return(res_tbl)
}


GO_analysis <- function(res_gene, sig_res_gene, sig_res_log2FoldChange, OrgDb){
  require(enrichplot)
  ego <- enrichGO(gene = sig_res_gene,
                  universe = res_gene,
                  keyType = "ENSEMBL",
                  OrgDb = OrgDb,
                  ont = "ALL",
                  pAdjustMethod = "fdr",
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  # Output results from GO analysis to a table
  cluster_summary <- data.frame(ego)
  
  if(nrow(cluster_summary) > 0){
    a <- dotplot(ego, orderBy = "x", showCategory = 10,  font.size = 15)
    pw_b <- enrichplot::pairwise_termsim(ego)
    b <- emapplot(pw_b, showCategory = 10)
    
    # To color genes by log2 fold changes 
    c <- cnetplot(ego,
                  categorySize="p.adjust",
                  showCategory = 10,
                  foldChange= sig_res_log2FoldChange,
                  vertex.label.font=6)
    lapply(list(a,b,c),plot)
  }
  return(cluster_summary)
}
