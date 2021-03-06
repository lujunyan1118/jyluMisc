% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug.R
\name{fitIC50}
\alias{fitIC50}
\title{Model object for fitting the IC50 curve}
\usage{
fitIC50(formula, data = NULL, weights, ...)
}
\arguments{
\item{formula}{Formula for the curve fitting.}

\item{data}{A data frame contain the raw concentration and the viability value. The viability should not be the percent viability value.}

\item{...}{Parameters passed to logLogisticRegression()}

\item{weigths}{Not used, mainly for geom_smooth() purpose}
}
\description{
Use the drc package to perform IC50 fit. Can be directly used for geom_smooth() in ggplot2
}
