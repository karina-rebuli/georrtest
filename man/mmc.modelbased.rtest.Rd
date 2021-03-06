% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mmc.modelbased.rtest.R
\name{mmc.modelbased.rtest}
\alias{mmc.modelbased.rtest}
\title{Model-based MMC for testing significance of Pearson correlation coeficient on spatial
autodependent data}
\usage{
mmc.modelbased.rtest(geodata, X, Y, coords, shuffle = 1, lambda = 1,
  likfit.X = NULL, likfit.Y = NULL, cov.model = "exp", kappa = NULL,
  nugget = 0, cov.pars.X = NULL, cov.pars.Y = NULL, nSimMMC = 1000, ...)
}
\arguments{
\item{geodata}{An object of geodata class. Mandatory if X or Y are not setted.}

\item{X}{X data vector. Mandatory.}

\item{Y}{Y data vector. Mandatory.}

\item{coords}{Data grid (locations). Mandatory.}

\item{shuffle}{If it should be permuted one or both datasets. Allowed values: "single", "s" or 1 and "double", "d" or 2. Default \code{1}.}

\item{lambda}{Box-Cox transformation parameter. A numeric with shuffle lenght. Default \code{1}.}

\item{likfit.X}{likfit{geoR} X data returned object. Default \code{NULL}.}

\item{likfit.Y}{likfit{geoR} Y data returned object. Default \code{NULL}.}

\item{cov.model}{A character with shuffle lenght. It's a string indicating the type of the
correlation function. Allowed values: see cov.spatial{geoR} manual. Default \code{"exp"}.}

\item{kappa}{A numeric with shuffle lenght, it's the smoothness parameter of the correlation function.
Only required by the following cov.model: "matern", "powered.exponential", "cauchy", "gencauchy"
and "gneiting.matern". Default \code{NULL}.}

\item{nugget}{tau^2 parameter for cov.model. A numeric with shuffle lenght. Default \code{0}.}

\item{cov.pars.X}{numeric(2) vector with sigma^2 and phi from X model. Default \code{NULL}.}

\item{cov.pars.Y}{numeric(2) vector with sigma^2 and phi from Y model. Default \code{NULL}.}

\item{nSimMMC}{Number of MMC simulations to be done. Allowed integers from 1e2 to 1e4.Default \code{1e3}.}
}
\value{
r.hat Sample Pearson correlation coefficient.

p.value The empirical p-value obtained.

r.perm A vector with estimated Pearson correlation coefficient in each MMC permutation.

differences.phis A vector with differences between phi OLS estimated with original data and with choosen permuted and
transformed data used as MMC permutation.
}
\description{
Model-based Monte Carlo method for testing Pearson correlation coeficient (r) significance under null
hypothesis (\eqn{H_0: r = 0}) between spatial autodependent data. It uses grf\{geoR\} to simulate GRFs
with specified parameters and evaluates r between X and Y to each simulation.
}
\note{
Warning: The function can be time consuming!

Warning: If likfit._ is not provided, this function will try to estimate spatial model parameters with default settings
with likfit\{geoR\} function. It's not a recommended use! Since the autodependence structure affects the correlation coefficient
distribution, a bad parameters estimation can lead to misleading results.
}
\author{
Karina Brotto rebuli and Paulo Justiniano Ribeiro Jr.
}
\keyword{coeficient}
\keyword{correlation}
\keyword{data,}
\keyword{empirical}
\keyword{geostatistical}
\keyword{mmc,}
\keyword{significance,}
\keyword{variograms}
