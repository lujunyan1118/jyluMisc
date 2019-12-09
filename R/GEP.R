
#' Function to run enrichment analysis using Camera from LIMMA
#'
#' A funtion to perform gene enrichment analysis using rotational test (fry or camera) implemented in limma package.
#' @param exprMat Gene expression matrix. Genes as rows and samples as columns.
#' @param design Design matrix as used by limma.
#' @param gmtFile A path to the gene signature file.
#' @param id A character vector of identifiers for genes used in the gmt file, ff the row names are not the same as the identfiers in the gmt file.
#' @param contrast The index of the colum in design matrix that should be used as contrast for enrichment analysis. Default is the last column.
#' @param method Wheter using "camera" (competitive) or "fry" (self-contained) method
#' @param pCut Cut-off for significance of enrichment.
#' @param ifDFR Whether raw p value or adjusted p value (FDR) should be used for significance cut-off
#' @param removePrefix Remove a certain string in the set names to make them shorter. For example, "HALLMARK_".
#' @param plotTitle The tilte used for enrichment bar plot.
#' @param insideLegned Whether to put the legend inside the bar plot to save space.
#' @export
#' @import limma
#'
runCamera <- function(exprMat, design, gmtFile, id = NULL,
                      contrast = ncol(design),  method = "camera", pCut = 0.05,
                      ifFDR = FALSE, removePrefix = NULL, plotTitle = "", insideLegend = FALSE) {
  
  #prepare indices
  if (is.null(id)) id <- rownames(exprMat)
  
  if (is.character(gmtFile)) {
    idx <- limma::ids2indices(piano::loadGSC(gmtFile)$gsc, id)
  } else {
    idx <- limma::ids2indices(gmtFile,id)
  }
  
  #run camera for fry
  if (method == "camera") {
    res <- limma::camera(exprMat, idx, design, contrast)
  } else if (method == "fry") {
    res <- limma::fry(exprMat, idx, design, contrast)
  }
  
  #plot enrichment results as bar plot
  
  plotTab <- res %>% rownames_to_column("Name")
  
  if (!is.null(removePrefix)) plotTab <- mutate(plotTab, Name = str_remove(Name, removePrefix))
  
  plotTab <- plotTab %>%
    mutate(Direction= factor(Direction, levels =c("Down","Up"))) %>%
    arrange(desc(Direction),desc(PValue)) %>%
    mutate(Name = factor(Name, levels = Name))
  
  if (ifFDR) {
    plotTab <- dplyr::filter(plotTab, FDR <= pCut)
  } else {
    plotTab <- dplyr::filter(plotTab, PValue <= pCut)
  }
  
  if (nrow(plotTab) == 0) {
    print("No sets passed the criteria")
    return(list(enrichTab = res, enrichPlot = NULL))
  } else {
    p <- ggplot(data = plotTab, aes(x = Name, y = -log10(PValue), fill = Direction)) +
      geom_bar(position = "dodge", stat = "identity", width = 0.5) +
      scale_fill_manual(values = c(Up = "blue", Down = "red")) +
      coord_flip() + xlab("") +
      ylab(expression(-log[10]*'('*italic(P)~value*')')) + ggtitle(plotTitle) +
      theme_bw() +
      theme(axis.text = element_text(size= 12))
    if (insideLegend) {
      p <- p + theme(legend.position = c(0.8,0.1))
    } else {
      p <- p + theme(legend.position = "right")
    }
    return(list(enrichTab = res, enrichPlot = p))
  }
}


#' Function to run enrichment analysis in R
#'
#' A funtion to perform GSEA or PAGE analysis
#' @param inputTab Ranked gene list for enrichment analysis, gene symbols as rownames and column stat contains the statistics for the ranking.
#' @param gmtFile A path to the gene signature file.
#' @param GSAmethod Method for enrichment analysis, currently supports GSEA and PAGE
#' @param nPerm Number of permutations for GSEA analysis
#' @export
#' @import piano

runGSEA <- function(inputTab,gmtFile,GSAmethod="gsea",nPerm=1000){
  require(piano)
  inGMT <- loadGSC(gmtFile,type="gmt")
  rankTab <- inputTab[order(inputTab[,1],decreasing = TRUE),,drop=FALSE] #re-rank by score
  
  if (GSAmethod == "gsea"){
    #readin geneset database
    #GSEA analysis
    res <- runGSA(geneLevelStats = rankTab,geneSetStat = GSAmethod,adjMethod = "fdr",
                  gsc=inGMT, signifMethod = 'geneSampling', nPerm = nPerm, verbose = FALSE)
    GSAsummaryTable(res)
  } else if (GSAmethod == "page"){
    res <- runGSA(geneLevelStats = rankTab,geneSetStat = GSAmethod,adjMethod = "fdr",
                  gsc=inGMT, signifMethod = 'nullDist', verbose = FALSE)
    GSAsummaryTable(res)
  }
}


#' Barplot for enrichment analysis result
#'
#' A function to plot enrichment analysis result from runGSEA function
#' @param resTab A list object of result tables from runGSEA function
#' @param pCut The cutoff for p-values of signatures to be plot
#' @param ifFDR Whether to use FDR cutoff or raw p value cutoff (default)
#' @param setName y axis label
#' @export
#' @import gridExtra ggplot2
plotEnrichmentBar <- function(resTab, pCut = 0.05, ifFDR = FALSE, setName = "Signatures") {
  pList <- list()
  rowNum <- c()
  for (i in names(resTab)) {
    plotTab <- resTab[[i]]
    if (ifFDR) {
      plotTab <- dplyr::filter(plotTab, `p adj (dist.dir.up)` <= pCut | `p adj (dist.dir.dn)` <= pCut)
    } else {
      plotTab <- dplyr::filter(plotTab, `p (dist.dir.up)` <= pCut | `p (dist.dir.dn)` <= pCut)
    }
    if (nrow(plotTab) == 0) {
      print("No sets passed the criteria")
      next
    } else {
      #firstly, process the result table
      plotTab <- apply(plotTab, 1, function(x) {
        statSign <- as.numeric(x[3])
        data.frame(Name = x[1], p = as.numeric(ifelse(statSign >= 0, x[4], x[6])), geneNum = ifelse(statSign >= 0, x[8], x[9]),
                   Direction = ifelse(statSign > 0, "Up", "Down"), stringsAsFactors = FALSE)
      }) %>% do.call(rbind,.)
      
      plotTab$Name <- sprintf("%s (%s)",plotTab$Name,plotTab$geneNum)
      plotTab <- plotTab[with(plotTab,order(Direction, p, decreasing=TRUE)),]
      plotTab$Direction <- factor(plotTab$Direction, levels = c("Down","Up"))
      plotTab$Name <- factor(plotTab$Name, levels = plotTab$Name)
      #plot the barplot
      pList[[i]] <- ggplot(data=plotTab, aes(x=Name, y= -log10(p), fill=Direction)) +
        geom_bar(position="dodge",stat="identity", width = 0.5) +
        scale_fill_manual(values=c(Up = "blue", Down = "red")) +
        coord_fixed(ratio = 0.5) + coord_flip() + xlab(setName) +
        ylab(expression(-log[10]*'('*p*')')) +
        ggtitle(i) + theme_bw() + theme(plot.title = element_text(face = "bold", hjust =0.5),
                                        axis.title = element_text(size=15))
      rowNum <-c(rowNum,nrow(plotTab))
    }
  }
  
  if (length(pList) == 0) {
    print("Nothing to plot")
  } else {
    rowNum <- rowNum
    grobList <- lapply(pList, ggplotGrob)
    grobList <- do.call(gridExtra::gtable_rbind,c(grobList,size="max"))
    panels <- grobList$layout$t[grep("panel", grobList$layout$name)]
    grobList$heights[panels] <- unit(rowNum, "null")
  }
  return(grobList)
}


