#' Fission yeast genotypes
#'
#' A matrix containing the genotypes of 
#' 46 different fission yeast strains.
#'
#' @format A matrix with 46 rows and 4481 columns: Each 
#' row corresponds to a strain, and each column to 
#' a genomic locus. Alleles are encoded as 0 or 1.
#' 
#' @source unpublished data.
"RFQTLgenotype"

#' chromosomal positions of the markers in RFQTLgenotype
#'
#' @format A dataframe with three columns, each
#' row corresponding to one locus.
#' \describe{
#'   \item{chr}{chromosome number}
#'   \item{start}{starting position of locus in bp}
#'   \item{end}{end position of locus in bp}
#' }
#' 
#' @source unpublished data.
"markerPositions"

#' chromosome Vector
#'
#' @format  A vector containing the chromosomes of the loci in RFQTLgenotype
#' 
#' @source unpublished data.
"chrVec"

#' Fission yeast phenotypes
#'
#' A dataframe containing the phenotype of 
#' 46 different fission yeast strains for two traits
#'
#' @format A matrix with 46 rows and 2 columns: Each 
#' row corresponds to a strain, and each column to 
#' a phenotype. 
#' 
#' @source unpublished data.
"RFQTLphenotype"

#' pre-computed permuted RF scores
#'
#' A matrix containing null distribution of importance scores
#' generated through permutations.
#'
#' @source unpublished data.
"perm1"


#' pre-computed permuted RF scores
#'
#' A matrix containing null distribution of importance scores
#' generated through permutations.
#'
#' @source unpublished data.
"perm2"


