\name{blocksdesign-package}
\alias{blocksdesign-package}
\alias{blocksdesign}
\docType{package}
\title{
Block designs package
}
\description{
The blocksdesign package provides for the construction of block designs for unstructured treatment
sets with arbitrary levels of replication not necessarily all equal and arbitrary depth of nesting. 
}
\details{
\tabular{ll}{
Package: \tab blocksdesign\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-09-16\cr
License: \tab GPL(>= 2)\cr
}
The main function in the \pkg{blocksdesign} package is \code{\link{blocks}} which provides functionality for building general nested block designs for unstructured treatment sets where treatments can have arbitrary replication and blocks can be nested to any arbitrary depth of nesting. Block sizes in the same stratum are as equal as possible in size and never differ by more than a single plot. Each nested blocks stratum is optimized subject to the treatments in the preceding higher strata remaining unchanged. In general, block designs are  optimized algorithmically by a swapping algorithm that maximizes the determinant of the information matrix (D-optimality) but certain special lattice designs (see vignette) are constructed algebraically.
}
\author{
R. N. Edmondson

Maintainer: R. N. Edmondson <rodney.edmondson@gmail.com>
}
\references{
EDMONDSON, R. N. (1998). Trojan square and incomplete Trojan square designs for crop research. Journal of Agricultural Science, Cambridge, 131, 135-142\cr
EDMONDSON, R. N. (2002) Generalised incomplete Trojan square designs. Biometrika, 89, 877-891\cr
JOHN, J. A. & WILLIAMS, E. R. (1995). Cyclic and Computer Generated Designs. Monographs on Statistics and Applied Probability 38. Chapman and Hall.\cr
}
\keyword{ design }

