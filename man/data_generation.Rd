% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{data_generation}
\alias{data_generation}
\title{Simulated data generation}
\usage{
data_generation(
  p = 150,
  n = 200,
  s = 5,
  d = 2,
  An = 5,
  M = 10,
  eta = 20,
  cov_type = "auto",
  res_type = "normal",
  seed = 111
)
}
\arguments{
\item{p}{Input integer. Default is 150. The dimensionality for the covariate.}

\item{n}{Input integer. Default is 200. The sample size for each datasets, including all the target and sources.}

\item{s}{Input integer. Default is 5. The sparsity level in the simulation.}

\item{d}{Input scalar. Default is 2. Control the scale of the Laplace distribution, larger d leads to smaller similarity between the target and informative sets. Details are in the section for simulation.}

\item{An}{Input integer. Default is 5. The number of informative set.}

\item{M}{Input integer. Default is 10. The total number of sources.}

\item{eta}{Input scalar. Default is 20. The signal-to-noise ratio for residual.}

\item{cov_type}{Input string. Default is 'auto'. Control the covariance matrix of datasets. The value 'auto' leads to a homogeneous design with auto-covariance; the value 'toeplitz' leads to heterogeneous designs with Toeplitz covariance matrix. No other value is permitted.}

\item{res_type}{Input string. Default is 'normal'. Control the type of distribution for residual. Only one from 'normal' and 'cauchy' is allowed.}

\item{seed}{Input integer. The default is 111.  The random seed.}
}
\value{
NULL. But will generate corresponding pairs of (X,y) in global environment.
}
\description{
Data generation function to implement the simulation settings in our paper.
}
