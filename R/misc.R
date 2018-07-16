# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


#' Load multiple data function
#'
#' This is a function to load a list of data in a specified folder
#' @param dataNames A list of the names of the RData files
#' @param folderDir The folder that contain the RData files
#' @param verbose Whether to show additional information
#' @export

loadRData <- function(dataNames, folderDir, verbose = TRUE) {
  #A function to load a list of data in a certain folder

  if (substr(folderDir,length(folderDir), length(folderDir)) != "/") {
    #if the last character is not "/", add one
    folderDir <- paste(folderDir,"/",sep = "")
  }

  for (eachName in dataNames) {
    if (eachName %in% list.files(folderDir)) {
      load(paste(folderDir,eachName,sep = ""),envir = .GlobalEnv)
      if (verbose) { print(sprintf("Loaded %s",eachName)) }
    } else print(sprintf("Can not find %s in %s",eachName, folderDir))
  }
}

#' Load the newst dataset
#'
#' Load the newest data object to R based on the file name
#' The file name should follow the format as "CPS1000_180126.RData"
#' @param prefix The prefix in the file name
#' @param dir The directory of the file
#' @param sep The separator the divides prefix and date
#' @export

loadData <- function(prefix, dir, sep = "_") {
  fileList <- grep(paste0("^",prefix,sep), list.files(dir), value = TRUE)
  dateStr <- lapply(fileList, function(x) {
    sp <- strsplit(x, split = sprintf("[.%s]",sep))[[1]]
    sp <- sp[length(sp)-1]
    if (nchar(sp) == 6) {
      sp <- as.Date(sp, format = "%y%m%d")
    } else if (nchar(sp) == 8) {
      sp <- as.Date(sp, format = "%Y%m%d")
    } else sp <- NA
    data.frame(fileName = x, date = sp, stringsAsFactors = FALSE)
  }) %>% bind_rows() %>% arrange(desc(date))
  if (nrow(filter(dateStr, !is.na(date))) == 0) stop("No file found")

  fileName <- dateStr$fileName[1]
  filePath <- file.path(dir,fileName)
  load(filePath,.GlobalEnv)
  print(sprintf("Loading data %s", filePath))
}

#' Generalized log transformation
#'
#' A function to perform generalized log transformation on a matrix, useful for dataset contians negative values
#' @param x The input matrix or datafram
#' @param q The quantile used for transformation
#' @export

glog <- function(x, q=0.03){
  c <- quantile(x,q,na.rm = TRUE)
  log((x+sqrt(x^2 +c^2))/2)
}

#' Calculate modified row z-scores
#'
#' Row centered by median and scaled by median absolute deviation
#' @param x The input matrix or datafram
#' @param center Centered by median
#' @param scale Scaled by 1.4826*MAD
#' @param censor Whether to censor the scaled value. If a positive integer value, the upper limit will be this value and the lower limit will be the negative value.
#' @param useMad Whether use row z score or modified z score (useMad = TRUE)
#' @export

mscale <- function(x, center = TRUE, scale = TRUE, censor = NULL, useMad = FALSE){
  if (useMad){
    x.scaled <- apply(x, 1, function(y) (y-median(y))/(1.4826*mad(y)))
  } else {
    x.scaled <- apply(x, 1, function(y) (y-mean(y))/(sd(y)))
  }
  if (!is.null(censor)) {
    x.scaled[x.scaled > censor] <- censor
    x.scaled[x.scaled < -censor] <- -censor
  }
  return(t(as.matrix(x.scaled)))
}


#' Find the intersect of multiple vectors
#'
#' A function to obtain the intersect of multiple vectors
#' @param ... Input vectors
#' @export

overlap <- function(...){
  vecList <- list(...)
  res <- vecList[[1]]
  for (i in seq(2,length(vecList))) {
    res <- intersect(res, vecList[[i]])
  }
  return(res)
}

#' Identify highly correlated variables in a dataset
#'
#' A function to identify highly correlated variables.
#' Return a Boolean vector of highly correlated variables that can be removed from the dataset
#' @param X The input dataframe for matrix. Columns as samples and rows as features.
#' @param cut The cutoff of correlation coefficient. The default value is 0.75
#' @param method The method used for the calculation of correlation coefficient, either "person" (default) or "spearman".
#' @export

highCor <- function(X, cut = 0.75, method = "pearson") {
  corMat <- cor(t(X), method = method, use = "pairwise.complete.obs")
  corMat[!lower.tri(corMat)] <- 0
  res <- apply(corMat, 1, function(x) any(abs(x) > cut))
  return(res)
}

#' Plot a list of ggplot object
#'
#' A function to plot a list of ggplot object in a pdf file
#' @param x A list of ggplot object
#' @param name The name of the output pdf file
#' @param ncol The number of columns in a page
#' @param nrow The number of rows in a page
#' @param figNum The number of total plots in a page
#' @param width The width of the page
#' @param height The height of the page
#' @export

makepdf <- function(x, name, ncol = 3, nrow = 2, figNum =NULL, width =20, height = 12) {
  require(gridExtra)
  if (is.null(figNum)) figNum = ncol*nrow
  if (length(x) == 0) return(NULL)
  pdf(name, width = width, height = height)
  for (i in seq(1, length(x), by = figNum)) {
    #print(names(x)[i])
    j  <- min(i + figNum - 1, length(x))
    do.call(grid.arrange, c(x[i:j], ncol = ncol, nrow = nrow))
  }
  dev.off() %>% invisible

}


#' Function to run enrichment analysis in R
#'
#' A funtion to perform GSEA or PAGE analysis
#' @param inputTab Ranked gene list for enrichment analysis, gene symbols as rownames and column stat contains the statistics for the ranking.
#' @param gmtFile A path to the gene signature file.
#' @param GSAmethod Method for enrichment analysis, currently supports GSEA and PAGE
#' @param nPerm Number of permutations for GSEA analysis
#' @export

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
#'
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




#' Model object for fitting the IC50 curve
#'
#' Use the drc package to perform IC50 fit. Can be directly used for geom_smooth() in ggplot2
#' @param formula Formula for the curve fitting.
#' @param data A data frame contain the raw concentration and the viability value. The viability should not be the percent viability value.
#' @param weigths Not used, mainly for geom_smooth() purpose
#' @param ... Parameters passed to logLogisticRegression()
#' @export
#' @import drc
fitIC50 <- function(formula, data = NULL, weights, ...) {
  if (! is.null(data) ) {
    modelFrame <- model.frame(formula, data)
  } else {
    modelFrame <- model.frame(formula)
  }
  parm_fit <- drm(modelFrame, fct = LL2.3u(), ...)
  newModel <- list(model = modelFrame, formula = formula, parm_fit = parm_fit)
  class(newModel) <- "fitIC50"
  return(newModel)
}


#' Predicted values based on IC50 fit
#'
#' Generic function for ic50 class generated from fitIC50 function.
#' @param object Object of class inheriting from "fitIC50"
#' @param newdata An optional data frame in which to look for variables with which to predict.If omitted, the fitted values are used.
#' @param se.fit Not used, mainly for geom_smooth purpose
#' @param level Not used, mainly for geom_smooth purpose
#' @param interval Not used, mainly for geom_smooth purpose
#' @export
#'
predict.fitIC50 <- function(object, newdata = NULL, se.fit = FALSE, level = 0.95 , interval = c("none", "confidence", "prediction")) {

  if (is.null(newdata))
    newdata <- object$model else
      newdata <- newdata

      parm_fit <- object$parm_fit

  res <- predict(parm_fit, newdata)

  return(res)
}
