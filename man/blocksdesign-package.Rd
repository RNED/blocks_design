\name{blocksdesign-package}
\alias{blocksdesign-package}
\alias{blocksdesign}
\docType{package}
\title{
Blocks design package
}
\description{
The \code{blocksdesign} package provides functionality for the construction of block designs for unstructured treatment sets with arbitrary levels of 
replication and arbitrary depth of nesting. 
}

\details{
\tabular{ll}{
Package: \tab blocksdesign\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-09-16\cr
License: \tab GPL(>= 2)\cr
}

\code{\link[blocksdesign]{blocks}} builds general nested block designs for unstructured treatment sets where treatments can have any number of levels of replication and blocks can be nested to any depth of nesting. Blocks in the same stratum are always  equal in size if the number of blocks exactly divides the number of plots otherwise they never differ by more than a single unit. Nested blocks are optimized hierarchically whereby each succesive set of nested blocks is optimized conditionally within the blocks of the preceding set. General block designs are optimized algorithmically by a swapping algorithm that finds a local D-optimal design but for certain special design sizes, block designs are constructed algebraically by using MOLS from the \pkg{crossdes} package function \code{\link[crossdes]{MOLS}}. See the vignette for further discussion of block designs.

Design outputs include a dataframe containing the treatment design, a table containing the treatment allocations to blocks, tables showing the blocks-by-treatments stratum incidence matrices and a table showing the A-efficiencies for each stratum together with efficiency upper bounds, where available (see John and Williams 1995, Chapter 2).

\code{\link[blocksdesign]{efficiencies}} is a secondary function which takes the design dataframe output from \code{blocks} and constructs tables of efficiency factors for each pairwise treatment difference in each stratum, if required.

\code{\link[blocksdesign]{upper_bounds}} is a subsidiary function that provides A-efficiency upper bounds for regular block designs with equally replicated treatments and equal block sizes. 
}

\author{
R. N. Edmondson

Maintainer: R. N. Edmondson <rodney.edmondson@gmail.com>
}
\references{

John, J. A. & Williams, E. R. (1995). Cyclic and Computer Generated Designs.Second Edition. Monographs on Statistics and Applied Probability 38. Chapman and Hall.\cr

Patterson H.D., Hunter E.A., (1983). The efficiency of incomplete block designs in National List and Recommended List cereal variety trials. J. Agric. Sci. Camb., 101, pp. 427-433.\cr
}
\keyword{ design }

