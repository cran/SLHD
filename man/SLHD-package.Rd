\name{SLHD-package}
\alias{SLHD-package}
\alias{SLHD}
\docType{package}
\title{
Sliced Latin hypercube design
}
\description{
Generate the optimal Latin hypercube designs and the optimal sliced Latin hypercube designs for computer experiments. 
}
\details{
\tabular{ll}{
Package: \tab SLHD\cr
Type: \tab Package\cr
Version: \tab 2.1-1\cr
Date: \tab 2015-01-26\cr
License: \tab LGPL-2.1\cr
}
This package contains functions for generating the optimal Latin hypercube designs (LHDs) when \eqn{t=1} and the optimal sliced Latin hypercube designs (SLHDs) when \eqn{t>1}. The maximin distance criterion is adopted as the optimality criterion. 

When \eqn{t=1}, the maximin-distance LHD is popularly used for designing computer experiments with quantitative factors.  

When \eqn{t>1}, the maximin-distance SLHD is a special class of LHD which can be partitioned into several slices (blocks), each of which is also a LHD of smaller size. The optimal SLHD structure guarantees the uniformity (space-filling property) in each slice as well as in the whole design. The SLHD is very important in designing computer experiments with quantitative and qualitative factors, where each slice is used as a design for quantitative factors under one of the t different level combinations of qualitative factors.

Important function in this package is \code{\link{maximinSLHD}}. 
}
\author{
Shan Ba

Maintainer: Shan Ba <shanbatr@gmail.com>
}
\references{
Ba, S., Brenneman, W. A. and Myers, W. R. (2015), "Optimal Sliced Latin Hypercube Designs," \emph{Technometrics}. 
}

\keyword{ package }

