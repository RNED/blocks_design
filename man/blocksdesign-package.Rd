\name{blocksdesign-package}
\alias{blocksdesign-package}
\alias{blocksdesign}
\docType{package}
\title{
Blocks design package
}
\description{
The blocksdesign package provides functionality for construction of block designs for unstructured treatment sets with arbitrary levels of 
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
\code{\link[blocksdesign]{blocks}} builds general nested block designs for unstructured treatment sets where treatments can have any number of levels of replication and blocks can be nested to any depth of nesting. Blocks in the same stratum are always  equal in size where possible and never differ by more than a single unit. Nested blocks are optimized hierarchically with each nested set of blocks optimized conditionally within the blocks of the preceding stratum. General block designs are optimized algorithmically by a swapping algorithm that finds a local D-optimal design but certain special lattice designs are constructed algebraically by using MOLS from the \pkg{crossdes} package function \code{\link[crossdes]{MOLS}}. Design outputs include a dataframe containing the block and treatment design, a table containing the treatments allocation to blocks, tables showing the blocks-by-treatments incidence matrices per stratum and a table showing the A-efficiencies together with efficiency upper bounds for each stratum, where available (see John and Williams 1995, Chapter 2).

\code{\link[blocksdesign]{pairwise_efficiencies}} is a secondary function which uses the design dataframe to provide tables of efficiency factors for pairwise treatment differences, if required.

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

