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

The purpose of block designs is to group experimental units into homogeneous blocks to provide maximum precision for treatment comparisons (see Edmondson 2005). The \code{blocksdesign} package supports the design of block experiments for unstructured treatment sets and provides a general purpose tool for the construction of practical nested block designs.

The function \code{\link[blocksdesign]{blocks}} builds general nested block designs for unstructured treatment sets where treatments can have any number of levels of replication and blocks can be nested to any feasible depth of nesting. If the number of blocks exactly divides the number of plots the blocks are always  equal in size, otherwise they never differ by more than a single unit. Where designs have one or more levels of nesting, blocks are optimized hierarchically with each successive set of nested blocks optimized within the blocks of the preceding set. General block designs are optimized algorithmically by a swapping algorithm that finds a locally D-optimal design. Special designs where the number of treatments t = v*v is equal to the square of the block size, the number of replicates is k+1 and there are k mutually orthogonal latin squares (MOLS) of size v x v, are constructed algebraically from the rows, columns and letters of the latin squares. If v is a prime-power, the MOLS are constructed using the \pkg{crossdes} package function \code{\link[crossdes]{MOLS}}. 

Design efficiency is assessed by an A-optimality criterion which can be compared with an upper-bound when designs have equal block sizes and equal replication (see John and Williams 1995, Chapter 2).  

Design outputs comprise a dataframe showing the treatment design, a table showing the treatment allocations to blocks, tables showing the blocks-by-treatments incidence matrices for each stratum in the design and a table showing the A-efficiencies for each stratum together with efficiency upper bounds, where available.

The \code{\link[blocksdesign]{efficiencies}} function is a secondary function which takes the design dataframe output from \code{blocks} and constructs tables of efficiency factors for each pairwise treatment difference in each stratum, if required.

The \code{\link[blocksdesign]{upper_bounds}} function is a subsidiary function that provides A-efficiency upper bounds for regular block designs with equally replicated treatments and equal block sizes. 

The vignette associated with the  package provides additional discussion of block designs and .
}

\author{
R. N. Edmondson

Maintainer: R. N. Edmondson <rodney.edmondson@gmail.com>
}
\references{

Edmondson, R. N. (2005) Past developments and future opportunities in the design and analysis of crop experiments. The Journal of Agricultural Science, Cambridge, 143, pp 27-33\cr

John, J. A. and Williams, E. R. (1995). Cyclic and Computer Generated Designs. Chapman and Hall, London
}
\keyword{ design }

