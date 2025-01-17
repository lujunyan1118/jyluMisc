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

#' Return MAD or meanAD based on whether MAD equals zero
#'
#' Calculate mean absolute deviation or median absolute deviation for robust z-score calcualte
#'
#' @param x Input numeric vector
#' @export

meanAD <- function(x) {
  md <- mad(x, na.rm = TRUE)
  if (md == 0)
    {md <- (sum(abs(x-mean(x,na.rm=T)),na.rm = T)/length(x))*1.253314}
  return(md)
}

#' Calculate modified row z-scores
#'
#' Row centered by median and scaled by median absolute deviation
#' @param x The input matrix or datafram
#' @param center Centered by median
#' @param scale Scaled by SD (or 1.4826*MAD)
#' @param censor Whether to censor the scaled value. If a positive integer value, the upper limit will be this value and the lower limit will be the negative value.
#' @param useMad Whether use row z score or modified z score (useMad = TRUE)
#' @export

mscale <- function(x, center = TRUE, scale = TRUE, censor = NULL, useMad = FALSE){
  if (scale & center) {
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) (y-median(y,na.rm = T))/meanAD(y))
    } else {
      x.scaled <- apply(x, 1, function(y) (y-mean(y,na.rm=T))/sd(y,na.rm = T))
    }
  } else if (center & !scale) {
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) (y-median(y,na.rm=T)))
    } else {
      x.scaled <- apply(x, 1, function(y) (y-mean(y,na.rm=T)))
    }
  } else if (!center & scale) {
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) y/meanAD(y))
    } else {
      x.scaled <- apply(x, 1, function(y) y/sd(y,na.rm = T))
    }
  } else {
    x.scaled <- t(x)
  }

  if (!is.null(censor)) {
    if (length(censor) == 1) {
      x.scaled[x.scaled > censor] <- censor
      x.scaled[x.scaled < -censor] <- -censor
    } else {
      x.scaled[x.scaled > censor[2]] <- censor[2] #higher limit
      x.scaled[x.scaled < censor[1]] <- censor[1] #lower limit
    }
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
#' @import gridExtra

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


#' Convert tidy table to SummarizedExperiment object
#'
#' Function for converting a tidy table to SummarizedExperiment object.
#' @param tidyTable Input tidy table
#' @param rowID A character variable specifying the column in the tidy table that should be used as row identifiers in the SE object.
#' @param colID A character variable specifying the column in the tidy table that should be used as column identifiers in the SE object.
#' @param values A character or a vector of characters, specifying the columns that should be treated as measurements (assays).
#' @param annoRow A character or a vector of characters, specifying the columns that should be used a rowData
#' @param annoCol A character or a vector of characters, specifying the columns that should be used a colData.
#' @export
#' @import SummarizedExperiment
#'
tidyToSum <- function(tidyTable, rowID, colID, values, annoRow, annoCol) {
  #prepare value matrix
  matList <- lapply(values, function(n) {
    dplyr::select(tidyTable, !!rowID, !!colID, !!n) %>%
      spread(key = !!colID, value = !!n) %>%
      column_to_rownames(rowID) %>%
      as.matrix()
  })
  names(matList) <- values

  #prepare row annoation
  rowAnno <- tidyTable[,unique(c(rowID, annoRow))]
  rowAnno <- rowAnno[!duplicated(rowAnno[[rowID]]),] %>% data.frame(stringsAsFactors = FALSE) %>%
    column_to_rownames(rowID)
  #if rowID is also used as row annotation, then copy it back to the data table
  if (rowID %in% annoRow) rowAnno[[rowID]] <- rownames(rowAnno)
  rowAnno <- rowAnno[rownames(matList[[1]]),,drop = FALSE ]

  #prepare column annotation
  colAnno <- tidyTable[,unique(c(colID, annoCol))]
  colAnno <- colAnno[!duplicated(colAnno[[colID]]),] %>% data.frame(stringsAsFactors = FALSE) %>%
    column_to_rownames(colID)
  #if colID is also used as column annotation, then copy it back to the data table
  if (colID %in% annoCol) colAnno[[colID]] <- rownames(colAnno)
  colAnno <- colAnno[colnames(matList[[1]]), ,drop=FALSE]

  # assemble
  m <- SummarizedExperiment(assays = matList, colData= colAnno)
  rowData(m) <- rowAnno
  return(m)
}

#' Convert SummarizedExperiment object to a tidy table
#'
#' Function for converting a SummarizedExperiment object to a tidy table.
#' @param seObject Input SummarizedExperiment object.
#' @param rowID A character variable specifying the name of row identifier.
#' @param colID A character variable specifying the name of column identifier.
#' @export
#' @import SummarizedExperiment
#'
#'
sumToTidy <- function(seObject, rowID = "rowID", colID = "colID") {
  #check if assay names are present and correct
  if (any(c("",NA) %in% assayNames(seObject))) {
    stop("Check assay names!")
  }

  tidyTable <- lapply(assayNames(seObject),function(n) {
    valTab <- assays(seObject)[[n]] %>%
      as_tibble(rownames = rowID) %>%
      gather(key = !!colID, value = "val", -!!rowID) %>%
      mutate(assay = n)
  }) %>% bind_rows() %>%
    spread(key = assay, value = val)

  #append row annotations
  if (rowID %in% colnames(rowData(seObject))) rowData(seObject)[[rowID]] <- NULL #if the specified rowID is also present
  rowAnno <- rowData(seObject) %>%
    as_tibble(rownames = rowID)

  tidyTable <- left_join(tidyTable, rowAnno, by = rowID)

  #append column annotations
  if (colID %in% colnames(colData(seObject))) colData(seObject)[[colID]] <- NULL #if the specified colID is also present
  colAnno <- colData(seObject) %>%
    as_tibble(rownames = colID)
  tidyTable <- left_join(tidyTable, colAnno, by = colID)


  return(as_tibble(tidyTable))
}

#' Generalized log2 transformation
#'
#' Function for generalized log2 transformation, can be applied to 0 and negative values
#' @param x Input number
#' @export
#'
glog2 <- function(x) {
  (asinh(x)-log(2))/log(2)
}

#' Perform associations test between columns from two tables
#'
#' Function to perform associations test between the columns from two tables.
#' The testing method will be selected automatically based on the value types.
#'
#' @param tabX the first table that contains columns to be tested
#' @param tabY the second table that contains columns to be tested
#' @param joinID the column that can be used to join the two tables
#' @param correlation_method method used for correlation test
#' @param plot whether to return a P-value heatmap
#' @param pCut cut-off for p-value or adjusted p-value
#' @param ifFdr whether to show nominal P-value or adjusted P-value
#' @param pMax censoring the -log10(P-value) or -log10(adjusted P-value) above this threshold. Default value is 12
#' @param onlySignificant only show associations passed certain significance level
#' @return a table with test results
#' @export
#'
#'
testAssociation <- function(tabX, tabY, joinID, correlation_method = "pearson", plot = FALSE,
                            pCut = 0.01, ifFdr = FALSE, pMax = 12, onlySignificant = FALSE) {
  fullTab <- left_join(tabX, tabY, by = joinID)

  #get column names
  colNamesX <- colnames(tabX)
  colNamesX <- colNamesX[colNamesX != joinID]
  colNamesY <- colnames(tabY)
  colNamesY <- colNamesY[colNamesY != joinID]

  resTab <- lapply(colNamesX, function(colX) {
    lapply(colNamesY, function(colY) {
      var1 <- fullTab[[colX]]
      var2 <- fullTab[[colY]]

      p <- tryCatch({
        if (is.numeric(var1) & is.numeric(var2)) {
          cor.test(var1, var2, method = correlation_method, use = "pairwise.complete.obs")$p.value
        } else if (! (is.numeric(var1) | is.numeric(var2))) {
          chisq.test(factor(var1), factor(var2))$p.value
        } else if (is.numeric(var1) & !is.numeric(var2)) {
          aovRes <- summary(aov(var1 ~ factor(var2)))
          aovRes[[1]][5][1,1]
        } else if (is.numeric(var2) & !is.numeric(var1)) {
          aovRes <- summary(aov(var2 ~ factor(var1)))
          aovRes[[1]][5][1,1]
        }}, error = function(err) NA)

      data.frame(var1 = colX, var2 = colY,
                 p = p,
                 stringsAsFactors = FALSE)

    }) %>% bind_rows()
  }) %>% bind_rows() %>%
    arrange(p) %>%
    mutate(p.adj = p.adjust(p, method ="BH"))

  if (!plot) {
    return(resTab)

  } else {
    #plot pvalue heatmap

    plotTab <- resTab %>%
      mutate(pPlot = ifelse(rep(ifFdr, nrow(resTab)), -log10(p.adj), -log10(p))) %>%
      mutate(pPlot = ifelse(pPlot <= -log10(pCut),0,pPlot)) %>%
      mutate(pPlot = ifelse(pPlot >= pMax, pMax, pPlot))

    if (onlySignificant) {
      sigTab <- filter(plotTab, pPlot !=0)
      plotTab <- filter(plotTab, var1 %in% sigTab$var1 , var2 %in% sigTab$var2)
    }

    #if "PC", "Factor" are detected in var1, then re-level var1 based on values without suffix
    if (any(str_detect(plotTab$var1,"Factor|PC"))) {
      plotTab <- mutate(plotTab, indexNoSuffix = str_remove(var1, "Factor|PC")) %>%
        mutate(indexNoSuffix = as.numeric(indexNoSuffix)) %>%
        arrange(desc(indexNoSuffix)) %>%
        mutate(var1= factor(var1, levels = unique(var1)))
    }

    p <- ggplot(plotTab, aes(x=var2, y=var1, fill = pPlot)) +
      geom_tile(color = "black") +
      scale_fill_gradient(low = "grey", high= "red", name = ifelse(ifFdr, "-log10(adj. Pval)","-log10(Pval)")) +
      theme_void() +
      scale_y_discrete(expand = c(0,0)) +
      scale_x_discrete(expand = c(0,0)) +
      ggtitle("Association P-value heatmap") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(hjust = 1),
            plot.title = element_text(face = "bold", hjust = 0.5))

    return(list(resTab=resTab, plot = p))
  }
}
