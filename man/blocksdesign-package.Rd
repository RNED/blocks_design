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
The main function in the \pkg{blocksdesign} package is \code{\link{blocks}} which provides functionality for building general nested block designs for unstructured treatment sets where treatments can have arbitrary levels of replication and blocks can be nested to any arbitrary depth of nesting. Block sizes in the same stratum are always as equal as possible in size and never differ by more than a single plot. Each nested blocks stratum is optimized subject to the treatments in the preceding higher blocks strata remaining unchanged. General block designs are optimized algorithmically by a swapping algorithm that maximizes the determinant of the information matrix (D-optimality) but certain special lattice designs (see vignette) are constructed algebraically using the \pkg{crossdes} package. 

Design outputs include tables showing the treatments allocation to blocks, the treatments-by-blocks incidence matrices for each blocks stratum and the A-efficiency of the optimized design together with upper efficiency bounds where available (see John and Williams 1995, Chapter 2).
}
\author{
R. N. Edmondson

Maintainer: R. N. Edmondson <rodney.edmondson@gmail.com>
}
\references{

JOHN, J. A. & WILLIAMS, E. R. (1995). Cyclic and Computer Generated Designs.Second Edition. Monographs on Statistics and Applied Probability 38. Chapman and Hall.\cr
}
\keyword{ design }

