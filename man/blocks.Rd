% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blocks.r
\name{blocks}
\alias{blocks}
\title{Block designs for unstructured treatment sets}
\usage{
blocks(treatments, replicates = 1, rows = NULL, columns = NULL,
  searches = NULL, seed = sample(10000, 1), jumps = 1)
}
\arguments{
\item{treatments}{a partition of the total required number of treatments into equally replicated treatment sets.}

\item{replicates}{a set of treatment replication numbers with one replication number for each partitioned treatment set.}

\item{rows}{the number of rows nested in each preceding block for each level of nesting from the top-level block downwards. The top-level block is a
single super-block which need not be defined explicitly.}

\item{columns}{the number of columns nested in each preceding block for each level of nesting from the top-level block downwards. The \code{rows} and 
\code{columns} parameters must be of equal length unless the \code{columns} parameter is null, in which case the design has a
single column block for each level of nesting and the design becomes a simple nested row blocks design.}

\item{searches}{the maximum number of local optima searched for a design optimization. The default number decreases
as the design size increases.}

\item{seed}{an integer initializing the random number generator. The default is a random seed.}

\item{jumps}{the number of pairwise random treatment swaps used to escape a local maxima. The default is a single swap.}
}
\value{
\item{Treatments}{A table showing the replication number of each treatment in the design.}
\item{Design}{Data frame giving the optimized block and treatment design in plot order.}
\item{Plan}{Data frame showing a plan view of the treatment design in the bottom level of the design.}
\item{BlocksEfficiency}{The D-efficiencies and the A-efficiencies of the blocks in each nested level of the design together with A-efficiency upper-bounds, where available.}
\item{seed}{Numerical seed used for random number generator.}
\item{searches}{Maximum number of searches used for each level.}
\item{jumps}{Number of random treatment swaps used to escape a local maxima.}
}
\description{
Constructs randomized nested block designs for unstructured treatment sets with any
feasible depth of nesting and up to two crossed block structures for each level of nesting.
}
\details{
Constructs randomized nested block designs with arbitrary depth of nesting for arbitrary unstructured treatment sets.

The \code{treatments} parameter is a set of numbers that partitions the total number of treatments into equally 
replicated treatment sets while the \code{replicates} parameter is a matching set of numbers that defines the replication of each equally 
replicated treatment set.

The \code{rows} parameter, if any, defines the number of nested row blocks for each level of nesting from the highest to the lowest. The
first number, if any, is the number of nested row blocks in the first-level of nesting, the second number, if any, is the number of nested row blocks in
the second-level of nesting and so on down to any required feasible depth of nesting.

The \code{columns} parameter, if any, defines the numbers of nested columns for each level of nesting,
where the first number is the 
number of column blocks crossed with the first set of nested row blocks, the second is the number of column blocks crossed with the second
 set of nested row blocks and so on for each level of the \code{rows} parameter.
 
If the \code{rows} and \code{columns} parameters are both defined they must be of equal length. If the number of columns for any
particular level of nesting is one, then that particular level of nesting will have a simple set of nested row blocks.
If both the \code{rows} parameter and the \code{columns} parameter are null, the default block design will be a set of orthogonal
main blocks equal in number to the highest common factor of the replication numbers. If the \code{rows} parameter is defined but the \code{columns}
parameter is null, the design will be a simple nested blocks design with nested block levels defined by the levels of the  \code{rows} parameter.

Block sizes are as nearly equal as possible and will never differ by more than a single plot for any particular block classification. 
Row blocks and column blocks always contain at least two plots per block and this restriction will constrain the permitted numbers of 
rows and columns for the various nested levels of a block design.

Unreplicated treatments are allowed and any simple nested block design can be augmented by any number of single unreplicated treatments 
to give augmented blocks that never differ in size by more than a single plot. General crossed block designs are more complex and currently 
the algorithm will only accommodate single unreplicated treatments in a crossed block design if the block sizes of the replicated part of 
the design are all equal in each nested level of the design.

For any particular level of nesting, the algorithm first optimizes the row blocks conditional on the next higher-level of blocks
and then optimizes the columns blocks, if any, conditional on the rows blocks.

Special designs:

Trojan designs are row-and-column designs for p replicates of v*p treatments arranged in p-rows and p-columns where v < p and 
where every row x column intersection contains v plots. Trojan designs have orthogonal rows and columns and optimal rows x columns
blocks and exist whenever p is prime or prime-power. Trojan designs are constructed algebraically from mutually
orthogonal Latin squares (MOLS).

Square lattice designs are resolvable incomplete block designs for r replicates of p*p treatments arranged in blocks of size p where
r < p+2 for prime or prime power p or r < 4 for general p. Square lattice designs are constructed algebraically from Latin squares or MOLS.

Lattice designs and Trojan designs based on prime-power MOLS require the \code{\link[crossdes]{MOLS}} package.

All other designs are constructed algorithmically using the D-optimality criterion.

Comment:

Row-and-column designs may contain useful treatment information in the individual row-by-column intersection blocks but the \code{blocks} function 
does not currently
optimize the efficiency of these blocks except for the special case of Trojan designs.

Row-and-column design with 2 complete treatment replicates, 2 complete rows and 2 complete columns will always confound one treatment contrast in the
rows-by-columns interaction. For these designs, it is impossible to nest a non-singular block design in the rows-by-columns intersections and instead
we suggest a randomized nested blocks design with four incomplete main blocks.

Outputs:

The principle design outputs comprise:parameters must be of the same length unless the \code{columns} parameter is null,
\itemize{
 \item  A data frame showing the allocation of treatments to blocks with successive nested strata arranged in standard block order. \cr
 \item  A table showing the replication number of each treatment in the design. \cr
 \item  A table showing the block levels and the achieved D-efficiency and A-efficiency factor for each nested level together
  with A-efficiency upper bounds, where available. \cr
 \item  A plan showing the allocation of treatments to blocks or to rows and to columns in the bottom level of the design.\cr
}
}
\examples{

## The number of searches in the following examples have been limited for fast execution.  
## In practice, the number of searches may need to be increased for optimum results.
## Designs should be rebuilt several times to check that a near-optimum design has been found.  


## Unstructured treatments partitioned into equally replicated treatment sets

# 3 treatments x 2 replicates + 2 treatments x 4 replicates in two complete randomized blocks
blocks(treatments=c(3,2),replicates=c(2,4),searches=10)

# 50 treatments x 4 replicates with 4 main blocks and 5 nested sub-blocks in each main block
blocks(treatments=50,replicates=4,rows=c(4,5))

# as above but with 20 additional single replicate treatments, one single treatment per sub-block
blocks(treatments=c(50,20),replicates=c(4,1),rows=c(4,5))

# 6 replicates of 6 treatments in 4 blocks of size 9 (non-binary block design)
blocks(treatments=6,replicates=6,rows=4)

# 4 replicates of 13 treatments arranged in a 13 x 4 Youden rectangle
blocks(treatments=13,replicates=4,rows=13,columns=4)

# 64 treatments x 2 replicates with nested 8 x 8 row-and-column designs in two main blocks
blocks(treatments=64,replicates=2,rows=c(2,8),columns=c(1,8),searches=10)

# 128 treatments x 2 replicates with two main blocks and 3 levels of nesting
\dontrun{blocks(128,2,c(2,2,2,2))}

# Durban et al example
 blocks(treatments=272,replicates=2,rows=c(16),columns=c(34))


}
\references{
Sailer, M. O. (2013). crossdes: Construction of Crossover Designs. R package version 1.1-1. https://CRAN.R-project.org/package=crossdes

Edmondson R. N. (1998). Trojan square and incomplete Trojan square designs for crop research. Journal of Agricultural Science, Cambridge, 131, pp.135-142

Cochran, W.G., and G.M. Cox. 1957. Experimental Designs, 2nd ed., Wiley, New York.
}
