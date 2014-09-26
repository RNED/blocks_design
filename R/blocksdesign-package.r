#' @name blocksdesign-package
#' @title Blocks design package
#' @aliases blocksdesign
#' @docType package
#' 
#' @description The \code{blocksdesign} package provides functionality for the construction of block designs for unstructured 
#' treatment sets with arbitrary levels of replication and arbitrary depth of nesting.
#'  
#' @details
#' 
#' Block designs group experimental units into homogeneous blocks to provide maximum precision for treatment comparisons. 
#' The most basic type of block design are complete randomised blocks where each block contains one or more complete sets of treatments. 
#' Complete randomized blocks are excellent for small designs but for larger designs, the variability within blocks may become too large 
#' for reliable treatment comparison and then it becomes desirable to subdivide the complete blocks into smaller incomplete blocks.
#'
#' Traditionally, nested block designs in large experiments have used a single set of nested blocks contained within a set of complete main blocks.
#' The complexity of design and analysis of designs with more than a single level of nesting have made hierarchical nesting infeasible for 
#' practical experiments. However, modern software such as the \code{lme4} mixed model analysis package and the availability of
#' modern design algorithms have removed these restrictions and the use of hierarchical nesting for large block designs is now
#' entirely feasible . 
#' 
#' The advantage of hierarchical nesting is that random variability can be captured across
#' a range of block sizes which allows for more realistic modelling of block effects compared with a single level of nesting. 
#' The \code{blocksdesign} package is a general purpose tool for the construction of general block designs where
#' treatments can have any number of levels of replication and blocks can be nested to any feasible depth of nesting. 
#' Where designs have one or more levels of nesting, blocks are optimized hierarchically with each successive set of nested blocks
#' optimized within the blocks of the preceding set. 
#'
#' The main package function is \code{\link[blocksdesign]{blocks}} which is used to generate the actual required design. The output from 
#' \code{blocks} includes a data frame showing the block and treatment factors for each plot, 
#' a data frame showing the allocation of treatments to plots for each block in the design,
#' a blocks-by-treatments incidence matrices for each stratum in the design and an A-efficiency factor for each stratum in the design
#'  together with an upper bound, where available.
#'  
#' The secondary package function \code{\link[blocksdesign]{efficiencies}} takes  the  the design data frame from the \code{blocks} 
#' function and constructs tables of efficiency factors for each pairwise treatment difference in each stratum, if required.
#' 
#' The \code{\link[blocksdesign]{upper_bounds}} function is a subsidiary function that estimates
#' A-efficiency upper bounds for regular block designs with equally replicated treatments and equal block sizes. 
#'
#' Further discussion of hierarchical nesting is in the package vignette at: vignette("blocksdesign")
#' 
NULL


