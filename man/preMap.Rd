% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RFQTL.R
\name{preMap}
\alias{preMap}
\title{Prepare mapping data}
\usage{
preMap(
  genotype,
  maxNAs = floor(nrow(genotype) * 0.5),
  propVar = 0.75,
  phenotype,
  sampleInfo,
  scale = T
)
}
\arguments{
\item{genotype}{A binary matrix with strain (rows) X marker (columns)
entries. Do not include markers with MAF=1.}

\item{maxNAs}{Maximum number of missing values that are allowed for
any marker. If this limit is exceeded, the marker is excluded 
from the genotype matrix.}

\item{propVar}{Proportion of genotypic variance that is included
in the population-structure covariates.}

\item{phenotype}{A vector or a matrix with numerical values. In 
case of a matrix, each row should correspond to a trait, and each
column to a strain. The matrix may contain NAs.}

\item{sampleInfo}{A vector with integers indicating to which strain
(row in the genotype) the measurement(s) in phenotype correspond.}

\item{scale}{Logical; Should centering and scaling of the phenotype be performed?}
}
\value{
List object containing the genotype, the groupings of the 
 markers, the mapping-covariates and the phenotype.
}
\description{
Genotype and phenotype data are pre-processed and reformatted
 for QTL-mapping with random forest.
}
\details{
The genotype matrix is checked and markers with too many missing
values are removed. Positions of remaining missing values and 
minor allele frequencies (MAF) are computed for later imputation
of missing genotypes. Highly correlated markers are collapsed to 
marker groups. The groupings are also stored and returned.
Population structure covariates, representing \code{propVar}%
of genetic variance are generated using emma. Phenotype values
are optionally scaled and centered. Genotype and phenotype matrices
are reordered according to the same sample order.
}
