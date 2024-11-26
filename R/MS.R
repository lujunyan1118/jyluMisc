# Functions for the analysis of mass-spectrometry data

#' Plot drop-out curve of an se object or matrix
#'
#' Plot the mean and missing value frequency (drop-out probability)
#' @param x a matrix or summarisedExepriment object
#' @param assayName the assay name or index if the input is a summarisedExperiment object
#' @export
#'

plotDropout <- function(x, useHex = FALSE, smooth = TRUE, assayName = 1) {
  if ("SummarizedExperiment" %in% class(x)) {
    exprMat <- assays(x)[[assayName]]
  } else {
    exprMat <- x
  }

  plotTab <- tibble(rowMean = rowMeans(exprMat, na.rm=TRUE),
                    perNA = rowMeans(is.na(exprMat)))
  p <- ggplot(plotTab, aes(x=rowMean, y=perNA)) +
    xlab("Mean expression") +
    ylab("Percentage of missing")

  if (useHex) {
    p <- p + geom_hex()
  } else {
    p <- p + geom_point()
  }

  if (smooth) {
    p <- p + geom_smooth()
  }

  return(p)
}
