setwd("/Users/giorgiomontesi/Desktop/Universita_di_Siena/A_PhD_Project/Biomarker_Prediction/ensembleBP/Codes")


a <- readRDS(file = "../Results/list_genes_1.rds")
gene_list <- a[[2]]


library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(tidyverse)
library(circlize)
library(reshape2)
library(RColorBrewer)
library(ggplot2)



#' @description
#' function to perform Reactome enrichment
#' @param gene_list list of genes in SYMBOL format on which to perform enrichment analysis
#' @param keyType Organism default "hgnc" referring to HomoSapiens
#' @param pval pvalue cutoff default 0.05
#' @param signif Logical TRUE: returns only significant enrichment results,
#'                      FALSE: returns all enrichment results
#' @returns DataFrame contsining enrichment results (all or only significant depending on signif)
enrich.React <- function(gene_list, keyType = "hgnc", pval = 0.05, signif = TRUE){
  
  reactome_results <- enrichPC(gene = gene_list, 
                               source = "reactome",  
                               keyType = keyType,    
                               pvalueCutoff = pval)
  
  if (signif == TRUE){
    padj <- pval
    return(reactome_results@result[reactome_results@result$p.adjust < padj,])
  } else {
      return(reactome_results@result)
    }
}

reactome_results_signif <- enrich.React(gene_list = gene_list)
reactome_results <- enrich.React(gene_list = gene_list, signif = F)


#' @description
#' function to perform GO Pathway enrichment
#' @param gene_list list of genes in on which to perform enrichment analysis
#' @param keyType format of gene_list, default "SYMBOL"
#' @param pval pvalue cutoff default 0.05
#' @param ont GO ontology on which to perform enrichment c("BP","MF","")
#' @param groupGO Logical TRUE: groups similar GO pathways together,
#'                       FALSE: returns all enriched GO pathways
#' @param signif Logical  TRUE: returns only significant enrichment results,
#'                       FALSE: returns all enrichment results
#' @param level level to group GO terms
#' @returns DataFrame containing enrichment results (all or only significant depending on signif)
enrich.GO <- function(gene_list, keyType = "SYMBOL", pval = 0.05, 
                      ont = "BP", groupGO = TRUE, signif = TRUE, level = 4){
  
  go_results <- enrichGO(gene = gene_list,
                         OrgDb = org.Hs.eg.db,
                         keyType = keyType,
                         ont = ont,
                         pvalueCutoff = pval, 
                         pool = 1)
  
  if (groupGO == TRUE & signif == TRUE){
    padj <- pval
    group_go <- groupGO(gene_list, OrgDb='org.Hs.eg.db', 
                        keyType = keyType, ont = ont, 
                        level = level)
    go_level <- group_go@result$ID
    go_signif <- go_results@result[go_results@result$p.adjust < padj, ]
    go_signif <- go_signif[go_signif$ID %in% go_level, ]
    return(go_signif)
    
  } else if (groupGO == TRUE & signif == FALSE){
    padj <- pval
    group_go <- groupGO(gene_list, OrgDb='org.Hs.eg.db', 
                        keyType = keyType, ont = ont, 
                        level = level)
    go_level <- group_go@result$ID
    go_signif <- go_results@result
    go_signif <- go_signif[go_signif$ID %in% go_level, ]
    return(go_signif)
    
  } else if (groupGO == FALSE & signif == TRUE){
    padj <- pval
    go_signif <- go_results@result[go_results@result$p.adjust < padj, ]
    return(go_signif)
    
  } else{
    return(go_results@result)
  }
}


go_grouped_signif <- enrich.GO(gene_list = gene_list, pval = 0.01)
go_grouped <- enrich.GO(gene_list = gene_list, pval = 0.01, signif = F)
go_results <- enrich.GO(gene_list = gene_list, pval = 0.01, signif = F, groupGO = F)
go_signif <- enrich.GO(gene_list = gene_list, pval = 0.01, signif = T, groupGO = F)




#' @description
#' function to perform circos plot
#' @param enrichment_results output of the enrichment. It is a df having in the rows
#'                           the enriched modules and in the columns the statistics of enrichemnt
#' @param palette_gene palette of colors to plot the genes default "Set2"
#' @param palette_modules palette of colors to plot the genes default "Paired"
#' @param transparency intensity of the color related to the chordDiagram, default =0.5
#' @param facing determines the orientation of the text relative to the circle circos.trackPlotRegion
#'               function,default ="clockwise"
#' @param cex determines the sixe of the text relative to the circle circos.trackPlotRegion function
#'            default = 0.7
#' @param legend Logical TRUE: shows the legend
#'                       FALSE: doesn't show the legend
#' @param legend_title assigns the title for the legend
#' @returns a circos plot displaying, for each significantly enriched modules,the genes
#'          from the list that are involved.
circos.plot <- function(enrichment_results,
                        palette_genes = "Set2", palette_modules = "Paired",
                        transparency = 0.5, facing = "clockwise",
                        cex = 0.7, legend = TRUE, legend_title ="circos_plot_legend"){
  
  
  geni_sel <- unique(unlist(str_split(enrichment_results$geneID, "/")))
  
  big_dat <- matrix(nrow = length(geni_sel), ncol = nrow(enrichment_results))
  rownames(big_dat) <- geni_sel
  colnames(big_dat) <- c(enrichment_results$ID)
  
  for (i in 1:ncol(big_dat)){
    big_dat[,i] <- as.numeric(rownames(big_dat) %in% 
                                unlist(str_split(enrichment_results$geneID[i], "/")))
  }
  
  
  big_dat <- as.data.frame(big_dat)
  df <- big_dat
  df$ID <- rownames(df)
  df <- melt(df, by = df$ID)
  
  colori_geni <- colorRampPalette(brewer.pal(brewer.pal.info[palette_genes, 1], 
                                             name = palette_genes))(length(unique(df$ID)))
  colori_moduli <- colorRampPalette(brewer.pal(brewer.pal.info[palette_modules, 1], 
                                               name = palette_modules))(length(unique(df$variable)))
  
  circos_enriched = function(){
    
    if (legend == T) {
      par(mar = c(0,0,0,15))
    } else {
      par(mar = c(0,0,0,0))
    }
    chordDiagram(df, 
                 transparency = 0.5, 
                 annotationTrack = c("grid", "axis"), 
                 preAllocateTracks = 1,
                 big.gap = 10,
                 grid.col = c(colori_geni, colori_moduli),
                 col = c(colori_geni, colori_moduli))
    
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
       xlim = get.cell.meta.data("xlim")
       ylim = get.cell.meta.data("ylim")
       sector.name = get.cell.meta.data("sector.index")
       circos.text(mean(xlim),
                   ylim[1] + .1,
                   sector.name,
                   facing = "clockwise",
                   niceFacing = TRUE,
                   adj = c(0, 0.5),
                   cex = 0.7)
      circos.axis(h = "top",
                  labels.cex = 0.2,
                  sector.index = sector.name,
                  track.index = 2)
    }, bg.border = NA)
    # 
    if (legend == T){ 
      # Aggiungi la legenda
      legend_df <- enrichment_results[, c("ID", "Description")]
      # Imposta le coordinate per la legenda
      par(xpd=TRUE)
      lgd.x <- par("usr")[2]*0.9
      lgd.y <- par("usr")[3]*0.1
      legend(x = lgd.x, y = lgd.y, 
           title = legend_title,
           title.cex = 0.7,
           legend = str_trunc(legend_df$Description, 40), 
           fill = colori_moduli, 
           col = "black", cex = 0.6, bty = "n", inset=c(0, 0.1))
    }
    
  }
  
  circos_enriched()
  
  
}


circos.plot(enrichment_results = reactome_results_signif,
            legend = T, 
            palette_genes = "Set1",palette_modules = "GnBu")

circos.plot(enrichment_results = go_grouped_signif,
            legend = T, 
            palette_genes = "Set1",palette_modules = "GnBu")




#'
#'
#'
#'
#'
#'
#'
bar.plot <- function(enrichment_results, pcutoff = 0.05, low.col = "indianred",
                     high.col = "lightblue"){
  
  library(viridis)
  
  enrichment_results <- enrichment_results[order(enrichment_results$p.adjust), ]
  
  ggplot(data = enrichment_results, aes(x = p.adjust, y = reorder(Description, -p.adjust), fill = p.adjust)) +
    geom_bar(stat = "identity") +
    #scale_fill_viridis(option = viridis.pal) +
    scale_fill_gradient(low = low.col, high = high.col) +
    geom_vline(xintercept = pcutoff, linetype="dashed", 
               color = "indianred") +
    theme_minimal()
  
  
}


bar.plot(enrichment_results = reactome_results[reactome_results$p.adjust< 0.5, ])








## TODO: vedi se riesci a inserire anche il KEGG enrichment, 
##       perchè non riconosce i SYMBOL ma solo gli ENTREZID
##       quindi potrebbe dare problemi nel grafico chordDiagram
##       visto che avresti i geni chiamati diversamente

# # Carica il pacchetto org.Hs.eg.db
# orgdb <- org.Hs.eg.db
# # Estrai gli Entrez IDs dai gene symbol
# entrez_ids <- mapIds(orgdb, keys = gene_list, column = "ENTREZID", keytype = "SYMBOL")
# # Rimuovi i geni che non hanno un Entrez ID associato
# entrez_ids <- entrez_ids[!is.na(entrez_ids)]
# # Utilizza gli Entrez ID per arricchire i pathway KEGG
# kegg_results <- enrichKEGG(gene = entrez_ids,
#                            organism = 'hsa',
#                            keyType = 'ncbi-geneid',
#                            pvalueCutoff = 0.05)



























