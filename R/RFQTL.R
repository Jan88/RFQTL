#' @title Prepare mapping data
#' 
#' @description
#' Genotype and phenotype data are pre-processed and reformatted
#'  for QTL-mapping with random forest. 
#'
#' @details
#' The genotype matrix is checked and markers with too many missing
#' values are removed. Positions of remaining missing values and 
#' minor allele frequencies (MAF) are computed for later imputation
#' of missing genotypes. Highly correlated markers are collapsed to 
#' marker groups. The groupings are also stored and returned.
#' Population structure covariates, representing \code{propVar}%
#' of genetic variance are generated using emma. Phenotype values
#' are optionally scaled and centered. Genotype and phenotype matrices
#' are reordered according to the same sample order.
#' 
#' @param genotype A binary matrix with strain (rows) X marker (columns)
#'  entries. Do not include markers with MAF=0.
#' @param maxNAs Maximum number of missing values that are allowed for
#'  any marker. If this limit is exceeded, the marker is excluded 
#'  from the genotype matrix.
#' @param propVar Proportion of genotypic variance that is included
#'  in the population-structure covariates.
#' @param phenotype A vector or a matrix with numerical values. In 
#'  case of a matrix, each row should correspond to a trait, and each
#'  column to a strain. The matrix may contain NAs. 
#' @param sampleInfo A vector with integers indicating to which strain
#'  (row in the genotype) the measurement(s) in phenotype correspond.
#' @param scale Logical; Should centering and scaling of the phenotype be performed?
#'
#' @return List object containing the genotype, the groupings of the 
#'  markers, the mapping-covariates and the phenotype.
#' @export
#'
#' @examples
preMap <- function(genotype, maxNAs=floor(nrow(genotype)*0.5), propVar=0.75, phenotype, sampleInfo, scale=T){
   
   ###check input
   
   if(length(unique(sampleInfo))<nrow(genotype)){
      stop("Don't include genotype-data for strains for which no trait data is supplied. Details can be found in the manual.")
   }
   if(!is.integer(sampleInfo)){
      stop("sampleInfo has to be an integer vector.")
   }
   if(any(apply(genotype,2,FUN=function(m){length(unique(m))==1}))){
      stop("Don't include markers with MAF=1.")
   }
   if(any(!is.numeric(maxNAs),maxNAs>length(sampleInfo),maxNAs<0)){
      stop("maxNAs has to be a single integer>=0 and <= to the number of samples")
   }
   if(!is.numeric(propVar)){
      stop("propVar has to be numeric.")
   }
   if(propVar>=1|propVar<0){
      stop("propVar has to be at least zero and has to be smaller than one.")
   }
   if(any(!is.numeric(phenotype),(!is.vector(phenotype))&(!is.matrix(phenotype)))){
      stop("phenotype has to be numeric and either a vector or a matrix")
   }
   if(length(sampleInfo)!=ifelse(is.matrix(phenotype),ncol(phenotype),length(phenotype))){
      stop("sampleInfo has to be as long as there are samples")
   }
   if(ifelse(is.vector(phenotype),any(is.na(phenotype)),any(apply(is.na(phenotype),2,all)))){
      stop("Don't include samples without any trait data.")
   }
   
   
   ###reduce genotype while retaining information
   
   #find markers with too many NAs
   nNA <- apply(genotype,2,FUN=function(x){sum(is.na(x))})
   
   #group redundant markers into genotype groups, leave out markers with too many NAs
   genotype2group <- rep(NA,ncol(genotype))
   reversedGenotype <- abs(genotype-1)
   
   redGeno <- NULL
   nGroups <- 0
   
   for(i in 1:ncol(genotype)){
      if(nNA[i]>maxNAs){
         next
      }
      if(nGroups>0){
         id2G <- apply(redGeno,2,FUN=function(x){
            identical(x,genotype[,i])|identical(x,reversedGenotype[,i])
         })
         if(any(id2G)){
            genotype2group[i] <- which(id2G)
         }else{
            redGeno <- cbind(redGeno,genotype[,i])
            nGroups <- nGroups+1
            genotype2group[i] <- nGroups
         }
      }else{
         redGeno <- cbind(redGeno,genotype[,i])
         genotype2group <- 1
         nGroups <- 1
      }
   }
   group2genotype <- lapply(1:nGroups,FUN=function(x){which(genotype2group==x)})
   
   ###estimate the population-structure using the reduced genotype
   K<-emma.kinship(t(redGeno),use="pairwise.complete.obs")
   KPCA <- eigen(K)
   useEigen <- 1:which(cumsum(KPCA[[1]])>=sum(KPCA[[1]])*propVar)[1]
   populationStructure <- KPCA[[2]][,useEigen,drop=F]
   colnames(populationStructure) <- paste0("popStr",1:ncol(populationStructure))
   
   ###parse the genotype and the population-structure covariates according to sampleInfo
   mappingGenotype <- redGeno[sampleInfo,,drop=F]
   populationStructure <- populationStructure[sampleInfo,,drop=F]
   
   ###prepare a list with missing genotypes for imputation
   NApos <- extractNAs(mappingGenotype)
   
   ###scale the phenotype if needed
   if(scale){
      if(is.vector(phenotype)){
         phenotype <- scale(phenotype,center=T,scale=T)[,1]
      }else{
         phenotype <- t(apply(phenotype,1,FUN=function(x){
            out <- scale(x,center=T,scale=T)[,1]
            return(out)
         }))
      }
   }
   
   ###return output
   return(list(genotype=mappingGenotype,
               genotype2group=genotype2group,
               group2genotype=group2genotype,
               mappingCovariates=populationStructure,
               NAlist=NApos,
               phenotype=phenotype))
}

#' @title Allele frequencies and missing genotypes.
#' 
#' @description
#' Function that calculates the allele frequencies (1 or 0) for 
#' each marker (column) and identifies missing values.
#'
#' @param genotype Genotype matrix with one row per individual, one
#'  column per marker. Entries should be binary marker information (1,0,NA).
#'
#' @return List with one vector for each marker. The first element 
#'  is the frequency of allele 1. The following numbers indicate the
#'  indices of any missing values for the respective samples.
#' @export
#'
#' @examples
extractNAs <- function(genotype){
   lapply(1:ncol(genotype),FUN=function(m){
     out <- which(is.na(genotype[,m]))
     return(out)
   })
}

#' @title Impute missing genotypes
#' 
#' @description
#' Function that replaces missing genotypes according to the
#' local allele frequencies.
#'
#' @param genotype Genotype as a matrix with 1,0,NA.
#' @param NAlist List object with one vector per predictor, containing
#'  the frequency of the 1-allele, followed by the indices of 
#'  individuals with missing values. As created by preMap.
#'
#' @return Matrix with the same dimensions as the input-genotype.
#'  NAs are replaced with alleles (0,1) according to the allele
#'  frequencies in the input (NAlist).
#' @export
#'
#' @examples
replaceGenoNAs <- function(genotype,NAlist){
   sapply(1:ncol(genotype),FUN=function(x){ #per predictor
      lNA <- length(NAlist[[x]])
      if(lNA==0){
         return(genotype[,x]) #if no values are missing, return the original vector
      }
      genotype[NAlist[[x]],x] <- sample(genotype[-NAlist[[x]],x],size=lNA,replace=T) #sample NAs according to NAlist and existing alleles
      return(genotype[,x]) #return the extended genotype
   })
}

#' @title Combined RF score
#' 
#' @description
#' Extracts a combined importance score from a randomForest object.
#'
#' @param rf A randomForest object that contains both permutation 
#'  importance (PI) and the increase in node purity (RSS). Random Forests
#'  that were generated with importance=TRUE are suited for this.
#' @param predictorsOfInterest Predictors for which the combined 
#'  score should be computed. As a default all predictors are selected.
#' @param normalize Should the scores be normalized with the trait variance?
#'  This makes the importance scores comparable across different traits.
#'
#' @return The output consists of a numeric vector containing a score 
#'  for each predictor of interest. The score is the product of RSS and PI.
#'  Negative RSS- or PI-scores are set to zero first. These scores are
#'  influenced by correlation between predictors and therefore contain
#'  a marker-specific bias.
#' @export
#'
#' @examples
cScore <- function(rf=NULL,predictorsOfInterest=NULL,normalize=TRUE){
   if(class(rf)!="randomForest"){
      stop("A randomForest object is required as input!")
   }
   if(is.null(predictorsOfInterest)){
      predictorsOfInterest <- 1:nrow(rf$importance)
   }else{
      if((!is.numeric(predictorsOfInterest))|(!is.vector(predictorsOfInterest))){
         stop("Predictors of interest have to be entered as a numeric vector!")
      }
   }
   if(ncol(rf$importance)<2){
      stop("A randomForest object with both permutation importance and increase in node purity (importance=TRUE) is required as input!")
   }
   if(normalize){ #divide scores by the trait variance. default
      if(is.null(rf$y)|rf$type!="regression"){
         stop("Normalization can only be performed if the randomForest object was created by supervised regression!")
      }
      tVar <- var(rf$y)
      PI <- rf$importance[predictorsOfInterest,1]/tVar
      PI[which(PI<0)] <- 0
      RSS <- rf$importance[predictorsOfInterest,2]/tVar
      RSS[which(RSS<0)] <- 0
   }else{ #don't normalize
      PI <- rf$importance[predictorsOfInterest,1]
      PI[which(PI<0)] <- 0
      RSS <- rf$importance[predictorsOfInterest,2]
      RSS[which(RSS<0)] <- 0
   }
   combinedScore <- PI*RSS
   names(combinedScore) <- rownames(rf$importance)[predictorsOfInterest]
   return(combinedScore)
}

#' @title Permutation scheme
#' 
#' @description
#' Function that generates a permutation scheme.
#'
#' @param permutationGroups List with an entry for each sample indicating
#'  within which group it may be permuted. Entries may contain more than
#'  one integer. If this case one is chosen randomly.
#'
#' @return Permuted order of individuals.
#' @export
#'
#' @examples
pscheme <- function(permutationGroups){
   
   permutationGroups <- sapply(permutationGroups,FUN=function(pG){
      if(length(pG)==1){
         return(pG)
      }else{
         return(sample(pG,size=1))
      }
   })
   out <- rep(NA,length(permutationGroups))
   for(i in unique(permutationGroups)){
      pos <- which(permutationGroups==i)
      out[pos] <- sample(pos,size=length(pos))
   }
   return(out)
}

#' @title Main function for RF QTL mapping
#' 
#' @description
#' Wrapper function for qtl mapping with random forest.
#'
#' @param mappingData A list containing all necessary information to map
#'  a trait as returned by preMap.
#' @param nTrait If the phenotype in mapping data is a matrix, nTrait 
#'  specifies which of the multiple traits provided should be mapped.
#'  Should be a single integer.
#' @param permutationGroups Optional; list with an entry for each sample
#'  indicating within which group it may be permuted. Entries may 
#'  contain more than one integer. If this is the case one is chosen
#'  randomly. Samples with the same integer may be permuted. If this
#'  list is not provided, all samples are assumed to be in the 
#'  same permutation-group.
#' @param permute Logical indicating whether permutations should be
#'  performed (if FALSE the real trait-values are used).
#' @param nPermutations Number of permutations to be performed. 
#'  If permute is false, nPermutations is ignored.
#' @param nforest Number of randomized imputations for missing genotypes
#' @param ntree Number of trees to grow per subforest
#' @param file Optional; if the results are to be saved as a file,
#'  this argument specifies the file paths
#' @param nCl Integer; number of cores to use
#' @param clType Character; type of cluster to create. 
#'  Has to be a cluster-type implemented in snow.
#' @param pMat Optional; matrix containing predefined permutation 
#'  schemes. The number of rows has to equal the number of samples and 
#'  the number of columns has to equal the number of permutations.
#'
#' @return For unpermuted traits: a numerical vector with importance 
#'  scores for all predictors aside from those in exclude. For permuted
#'  traits: a matrix of scores with one column per permutation
#' @export
#'
#' @examples
rfMapper <- function(mappingData, nTrait=NULL, permutationGroups=NULL, permute, nPermutations=1000, nforest, ntree, file=NULL, nCl=1, clType="SOCK", pMat=NULL){
   if(!is.logical(permute)){
      stop("permute has to be TRUE or FALSE.")
   }
   if(any(!c("genotype","genotype2group","group2genotype","mappingCovariates","NAlist","phenotype")%in%names(mappingData))){
      stop("mappingData has to generated with preMap.")
   }
   if(ntree<1|nforest<1){
      stop("Choose parameters so at least one tree is built.")
   }
   if(!is.null(file)){
      dirString <- strsplit(x = file,split = "")[[1]]
      if(grepl("/",file)){
         dir <- paste(dirString[1:tail(which(dirString=="/"),n=1)],collapse = "")
         if(!file.exists(dir)){
            stop("Directory doesn't exist.")
         }
      }
      if(!grepl(".RData",file,ignore.case=T)){
         stop("Specify a RData file.")
      }
   }
   if(!permute|nCl<=1){
      out <- rfMapperNonPar(mappingData=mappingData,
                            nTrait=nTrait,
                            permutationGroups=permutationGroups,
                            permute=permute,
                            nPermutations=nPermutations,
                            nforest=nforest,
                            ntree=ntree,
                            file=file,
                            pMat=pMat)
      return(out)
   }else{
      out <- rfMapperPar(mappingData=mappingData,
                         nTrait=nTrait,
                         permutationGroups=permutationGroups,
                         nPermutations=nPermutations,
                         nforest=nforest,
                         ntree=ntree,
                         file=file,
                         nCl=nCl,
                         clType=clType,
                         pMat=pMat)
      return(out)
   }
}

#' @title Parallel RF mapping
#' 
#' @description
#' Parallelized implementation of trait mapping with Random Forest.
#'
#' @param mappingData A list containing all necessary information to 
#'  map a trait as returned by preMap.
#' @param nTrait If the phenotype in mapping data is a matrix, nTrait
#'  specifies which of the multiple traits provided, should be mapped.
#'  Should be a single integer.
#' @param permutationGroups Optional; list with an entry for each 
#'  sample indicating within which group it may .
#' @param file Optional; if the results are to be saved as a file,
#'  this argument specifies the exact path.
#' @param permute Logical, should permutations be performed (if false
#'  the real trait-values are used)?
#' @param nPermutations Number of permutations to be performed. If
#'  permute is false, nPermutations is irrelevant.
#' @param nforest Number of randomizations for missing genotypes.
#' @param ntree Number of trees to grow per subforest.
#' @param pMat Optional; matrix containing predefined permutation 
#'  schemes. The number of rows has to equal the number of samples 
#'  and the number of columns has to equal the number of permutations.
#' @param nCl Integer; number of cores to use.
#' @param clType Character; type of cluster to create. 
#'  Has to be a cluster-type implemented in snow.
#'
#' @return A matrix of scores with one column per permutation.
#' @export
#'
#' @examples
rfMapperPar <- function(mappingData, nTrait=NULL, permutationGroups=NULL, 
                        file=NULL, permute=T, nPermutations=1000, nforest,
                        ntree, pMat=NULL, nCl, clType="SOCK"){
   
   genotype <- mappingData$genotype
   phenotype <- mappingData$phenotype
   
   if(!is.vector(phenotype)){
      if(!is.null(nTrait)&is.numeric(nTrait)&length(nTrait)==1){
         if(nTrait%%1!=0|nTrait<1|nTrait>nrow(phenotype)){
            stop("Specify nTrait as a natural number between 1 and the number traits supplied.")
         }
         phenotype <- phenotype[nTrait,]
      }else{
         stop("Specify which trait to map.")
      }
   }
   
   NAlist <- mappingData$NAlist
   mappingCovariates <- mappingData$mappingCovariates
   exclude <- (ncol(genotype)+1):(ncol(genotype)+ncol(mappingCovariates))
   genotype <- cbind(genotype,mappingCovariates)
   if(is.null(permutationGroups)){
      permutationGroups <- as.list(rep(1,nrow(genotype)))
   }
   
   if(any(is.na(phenotype))&is.null(pMat)){ #remove individuals with missing trait-values if pMat
      if(all(is.na(phenotype))){
         stop("No data available for the specified trait.")
      }
      missing <- which(is.na(phenotype))
      phenotype <- phenotype[-missing]
      genotype <- genotype[-missing,]
      NAlist <- lapply(NAlist,FUN=function(markerVec){
         aFreq <- markerVec[1]
         markerVec <- markerVec[-1]
         if(length(markerVec)==0){
            return(aFreq)
         }
         markerVec <- sapply(markerVec,FUN=function(sampleInd){
            sampleInd-sum(sampleInd>missing)
         })
         return(c(aFreq,markerVec))
      })
      permutationGroups <- permutationGroups[-missing]
   }
   
   #make the cluster
   require(snow)
   cl <- makeCluster(nCl,type=clType)
   clusterExport(cl,
                 list = c("cScore",
                          "pscheme",
                          "genotype",
                          "phenotype",
                          "permutationGroups",
                          "file",
                          "nforest",
                          "ntree",
                          "pMat",
                          "NAlist",
                          "replaceGenoNAs",
                          "exclude"),
                 environment())  
   clusterEvalQ(cl,library(randomForest))
   
   
   
   out <- parSapply(cl=cl,X=1:nPermutations,FUN=function(nPerm){
      if(!is.null(pMat)){
         pScheme <- pMat[,nPerm] #use a permutation scheme from pMat
      }else{
         pScheme <- pscheme(permutationGroups) #generate a permutation scheme
      }
      permPhenotype <- phenotype[pScheme] #permute the trait vector
      permGenotype <- genotype 
      permGenotype[,exclude] <- permGenotype[pScheme,exclude] #preserve the association between population structure and trait values
      if(any(is.na(phenotype))&!is.null(pMat)){ #remove individuals with missing trait-values
         missing <- which(is.na(phenotype))
         phenotype <- phenotype[-missing]
         permGenotype <- permGenotype[-missing,]
         NAlist <- lapply(NAlist,FUN=function(markerVec){
            aFreq <- markerVec[1]
            markerVec <- markerVec[-1]
            if(length(markerVec)==0){
               return(aFreq)
            }
            markerVec <- sapply(markerVec,FUN=function(sampleInd){
               sampleInd-sum(sampleInd>missing)
            })
            return(c(aFreq,markerVec))
         })
      }
      rfs <- lapply(1:nforest,FUN=function(nF){
         permGenotype[,-exclude] <- replaceGenoNAs(genotype=permGenotype[,-exclude],NAlist=NAlist)
         rf <- randomForest(x=permGenotype, y=permPhenotype, ntree=ntree, importance=TRUE, keep.forest=F)
         return(rf)
      })
      
      bigRf <- do.call("combine",rfs) #combine the subforests that were generated with genotypes that differ for the NA positions
      if(!is.null(exclude)){
         scores <- cScore(bigRf)[-exclude] #don't keep the importance scores of the population-structure, since it is not permuted
      }else{
         scores <- cScore(bigRf)
      }
      return(scores)
   })
   stopCluster(cl)
   if(!is.null(file)){
      save(out,file=file)
   }
   return(out) #return the scores as a vector for unpermuted traits and as a matrix for permuted traits
}

#' @title Non-parallel RF mapping
#' 
#' @description
#' Non-parallelized implementation of trait mapping with Random Forest.
#'
#' @param mappingData A list containing all necessary information to
#'  map a trait as returned by preMap.
#' @param nTrait If the phenotype in mapping data is a matrix, nTrait
#'  specifies which of the multiple traits provided, should be mapped.
#'  Should be a single integer.
#' @param permutationGroups Optional; list with an entry for each sample
#'  indicating within which group it may be permuted. Entries may contain
#'  more than one integer. If this is the case one is chosen randomly.
#' @param permute Logical, should permutations be performed (if false 
#'  the real trait-values are used)?
#' @param nPermutations Number of permutations to be performed. If 
#'  permute is false, nPermutations is irrelevant.
#' @param nforest Number of randomizations for missing genotypes
#' @param ntree Number of trees to grow per subforest.
#' @param file Optional; if the results are to be saved as a file, 
#'  this argument specifies the exact path.
#' @param pMat Optional; matrix containing predefined permutation 
#'  schemes. The number of rows has to equal the number of samples 
#'  and the number of columns has to equal the number of permutations.
#'
#' @return A numerical vector with importance scores for all predictors,
#'  not including covariates.
#' @export
#'
#' @examples
rfMapperNonPar <- function(mappingData, nTrait=NULL, permutationGroups=NULL,
                           permute, nPermutations=1000, nforest, ntree, 
                           file=NULL, pMat=NULL){
   
   genotype <- mappingData$genotype
   phenotype <- mappingData$phenotype
   if(!is.vector(phenotype)){
      if(!is.null(nTrait)&is.numeric(nTrait)&length(nTrait)==1){
         if(nTrait%%1!=0|nTrait<1|nTrait>nrow(phenotype)){
            stop("Specify nTrait as a natural number between 1 and the number traits supplied.")
         }
         phenotype <- phenotype[nTrait,]
      }else{
         stop("Specify which trait to map.")
      }
   }
   NAlist <- mappingData$NAlist
   mappingCovariates <- mappingData$mappingCovariates
   exclude <- (ncol(genotype)+1):(ncol(genotype)+ncol(mappingCovariates))
   genotype <- cbind(genotype,mappingCovariates)
   if(is.null(permutationGroups)){
      permutationGroups <- as.list(rep(1,nrow(genotype)))
   }
   
   if(any(is.na(phenotype))&is.null(pMat)){ #remove individuals with missing trait-values if pMat is not provided
      if(all(is.na(phenotype))){
         stop("No data available for the specified trait.")
      }
      missing <- which(is.na(phenotype))
      phenotype <- phenotype[-missing]
      genotype <- genotype[-missing,]
      NAlist <- lapply(NAlist,FUN=function(markerVec){
         aFreq <- markerVec[1]
         markerVec <- markerVec[-1]
         if(length(markerVec)==0){
            return(aFreq)
         }
         markerVec <- sapply(markerVec,FUN=function(sampleInd){
            sampleInd-sum(sampleInd>missing)
         })
         return(c(aFreq,markerVec))
      })
      permutationGroups <- permutationGroups[-missing]
   }
   if(permute){ #do permutations
      out <- sapply(1:nPermutations,FUN=function(nPerm){
         if(!is.null(pMat)){
            pScheme <- pMat[,nPerm] #use a permutation scheme from pMat
         }else{
            pScheme <- pscheme(permutationGroups) #generate a permutation scheme
         }
         permPhenotype <- phenotype[pScheme] #permute the trait vector
         permGenotype <- genotype
         permGenotype[,exclude] <- permGenotype[pScheme,exclude] #preserve the association between population structure and trait values
         if(any(is.na(phenotype))&!is.null(pMat)){ #remove individuals with missing trait-values
            missing <- which(is.na(phenotype))
            phenotype <- phenotype[-missing]
            permGenotype <- permGenotype[-missing,]
            NAlist <- lapply(NAlist,FUN=function(markerVec){
               aFreq <- markerVec[1]
               markerVec <- markerVec[-1]
               if(length(markerVec)==0){
                  return(aFreq)
               }
               markerVec <- sapply(markerVec,FUN=function(sampleInd){
                  sampleInd-sum(sampleInd>missing)
               })
               return(c(aFreq,markerVec))
            })
         }
         rfs <- lapply(1:nforest,FUN=function(nF){
            permGenotype[,-exclude] <- replaceGenoNAs(genotype=permGenotype[,-exclude],NAlist=NAlist)
            rf <- randomForest(x=permGenotype, y=permPhenotype, ntree=ntree, importance=TRUE)
            return(rf)
         })
         bigRf <- do.call("combine",rfs) #combine the subforests that were generated with genotypes that differ for the NA positions
         scores <- cScore(bigRf)[-exclude] #don't keep the importance scores of the population-structure, since it is not permuted
      })
   }else{ #map real traits
      rfs <- lapply(1:nforest,FUN=function(nF){
         genotype[,-exclude] <- replaceGenoNAs(genotype=genotype[,-exclude],NAlist=NAlist)
         rf <- randomForest(x=genotype, y=phenotype, ntree=ntree, importance=TRUE)
         return(rf)
      })
      bigRf <- do.call("combine",rfs)
      out <- cScore(bigRf)[-exclude]
   }
   
   if(!is.null(file)){
      save(out,file=file)
   }
   return(out) #return the scores as a vector for unpermuted traits and as a matrix for permuted traits
}


#' @title Group Markers
#' 
#' @description
#' Function that joins consecutive markers that were significantly
#' linked to the same trait to a QTL.
#'
#' @param QTLmat Matrix consisting of two columns (target, predictor)
#'  and one row per linked predictor.
#' @param chrVec Character-vector containing the chromosome on which 
#'  a given marker is located.
#'
#' @return List with one entry for each regulated trait. Each entry 
#'  contains the target and a matrix with the start and end of the QTL-regions.
#' @export
#'
#' @examples
joinConsecutive <- function(QTLmat, chrVec){
   
   #identify traits with at least one qtl
   targets <- sort(unique(QTLmat[,1]))
   
   #join consecutive significant predictors per trait
   QTLperTarget <- lapply(targets,FUN=function(target){
      targetMat <- QTLmat[which(QTLmat[,"target"]==target),] #subset for the trait
      if(is.vector(targetMat)){ #if there is only a single significant predictor, there is no joining to do
         outMat <-matrix(c(targetMat[2],targetMat[2]),nrow=1)
         colnames(outMat) <- c("start", "end")
         return(list(target=target,predictors=outMat))
      }
      predictors <- sort(targetMat[,2])
      col <- predictors[1]
      outMat <- NULL
      for(i in 2:length(predictors)){
         #join loci if they are on the same chromosome and are direct neighbors
         if(chrVec[col[1]]==chrVec[predictors[i]]&(predictors[i]-col[length(col)])==1){
            col <- c(col,predictors[i])
         }else{
            outMat <- rbind(outMat,c(min(col),max(col)))
            col <- predictors[i] #start with a new region
         }
      }
      outMat <- rbind(outMat,c(min(col),max(col)))
      colnames(outMat) <- c("start","end")
      #return one matrix per trait containing start and end points of all qtl
      return(list(target=target,predictors=outMat)) 
   })
   #return a list containg these matrices for all regulated traits
   return(QTLperTarget)
}

#' @title Join close QTL
#' 
#' @description
#' Function that joins different QTL to a single one if they 
#' contain highly correlated markers.
#'
#' @param QTLlist Output from joinConsecutive.
#' @param corMat Squared matrix of predictor-correlation coefficients.
#' @param corThreshold Threshold indicating how strong QTL-regions have
#'  to be correlated to be joined.
#' @param distThreshold Threshold indicating how close QTL-regions 
#'  have to be to be joined.
#' @param chrVec Character vector containing the information on which
#'  chromosome a marker is located.
#'
#' @return List with an entry for each regulated trait. Each entry 
#'  contains the target and a matrix with the start and end of the 
#'  QTL regions after combining close and correlated QTL.
#' @export
#'
#' @examples
joinNear <- function(QTLlist, corMat, corThreshold, distThreshold, chrVec){
   
   out <- lapply(QTLlist,FUN=function(targetList){
      qtl <- targetList[[2]] #matrix with QTL-regions for this trait
      if(nrow(qtl)==1){ #if there is only one region, no regions can be joined
         return(list(target=targetList[[1]],predictors=qtl))
      }
      prog <- T #progress-variable
      while(prog){ #as long as there is progress
         nextProg <- F #loop stops if there is no change
         #each qtl is tested to its right neighbor
         for(i in 1:(nrow(qtl)-1)){
            #are the QTL correlated and close enough?
            if(max(corMat[qtl[i,1]:qtl[i,2],qtl[i+1,1]:qtl[i+1,2]])>=corThreshold&(qtl[i+1,1]-qtl[i,2])<=distThreshold){
               qtl[i,] <- c(qtl[i,1],qtl[i+1,2]) #end of the combined QTL is the end of the second QTL
               qtl <- qtl[-(i+1),] #the second QTL is deleted as it is part of the combined QTL now
               qtl <- matrix(as.vector(qtl),ncol=2,byrow=F) #if there is one QTL remaining it still has the same format
               nextProg <- T
               break #break for loop because indices and matrix size are changed
            }
         }
         prog <- nextProg #prepare for the next iteration
         if(nrow(qtl)==1){ #if only one QTL remains, no further joining can be done
            prog <- F 
         }
      }
      colnames(qtl) <- c("start","end")
      #return the target and the QTL-matrix with start and end points
      return(list(target=targetList[[1]],predictors=qtl)) 
   })
   return(out)
}


#' @title Join correlated QTL
#' 
#' @description
#' Loin QTL-regions that are distant from each other but contain
#' highly correlated markers to QTL-groups using hierarchical clustering.
#'
#' @param QTLlist Output of joinNear.
#' @param corMat Squared matrix of predictor-correlation coefficients.
#' @param corThreshold Threshold indicating how strong QTL-regions have to
#'  be correlated to be joined.
#'
#' @return List with an entry for each QTL. Each entry contains the 
#'  target and a matrix with the start and end of the QTL-regions after
#'  grouping highly correlated distant QTL.
#' @export
#'
#' @examples
joinCorrelated <- function(QTLlist, corMat, corThreshold){
   
   out <- lapply(QTLlist,FUN=function(targetList){
      qtl <- targetList[[2]] #matrix with QTL-regions for this trait
      if(nrow(qtl)==1){ #if there is only one region, no regions can be joined
         return(list(list(target=targetList[[1]],predictors=qtl)))
      }
      mat <- sapply(1:nrow(qtl), function(i) { #correlation matrix
         sapply(1:nrow(qtl), function(j) {
            return(max(corMat[qtl[i,1]:qtl[i,2],qtl[j,1]:qtl[j,1]], na.rm=T))
         })
      })
      colnames(mat) <- rownames(mat) <- 1:nrow(qtl)
      diag(mat) <- 0 # remove the diagonal
      mat[is.na(mat)] <- 0
      #here we convert the correlation-matrix into distances for clustering 
      dcorm <- as.dist((1 - mat)/2) 
      dcor_cut <- (1-corThreshold)/2 #convert the threshold accordingly
      dcorm_hclust <- hclust(dcorm)
      group <- cutree(hclust(dcorm), h=dcor_cut) #apply the threshold to the tree
      group <- split(as.numeric(names(group)), group)
      out <- lapply(group,FUN=function(gr){ #join regions to QTL-groups
         outMat <- matrix(qtl[gr,],nrow=length(gr),byrow = F)
         colnames(outMat) <- c("start","end")
         list(target=targetList[[1]],predictors=outMat)
      })
      return(out)
   })
   out <- do.call(c,out) #parse list so that one QTL-group is one entry regardless of the target
   names(out) <- NULL
   return(out)
}

#' @title Estimate empirical p-value
#' 
#' @description
#' Function that uses finished permutations and real scores 
#' to estimate an empirical p-value.
#'
#' @param path File path where the .RData files with the matrices of 
#'  empirical values can be found. Thse have to conform to predictors 
#'  (rows) X permutations (columns).
#' @param markersPerIteration Total number of null distributions that 
#'  are loaded at the same time. The possible number is dependent on 
#'  available working memory and size of the null distributions.
#' @param scores Real predictor-scores to be compared to the empirical 
#'  null distributions. If scores is a matrix, it has to conform to 
#'  traits (rows) X predictors (columns). Scores for the same predictor
#'  in matrix, it has to conform to traits (rows) X predictors (columns).
#'  Scores for the same predictor in different traits (rows) are 
#'  compared to the same null distribution.
#' @param printProg Logical; if TRUE the last finished predictor is printed.
#' @param pCorrection Which type, if any, of multiple testing correction.
#'  should be performed? Has to be either "none", "fdr" or "bonferroni".
#'
#' @return A matrix with empirical p-values is returned.
#' @export
#'
#' @examples
pEst <- function(path,markersPerIteration,scores,printProg=T,pCorrection="none"){
   
   #if scores is a vector (i.e. only one trait), it is transformed into a matrix
   if(is.vector(scores)){
      scores <- matrix(scores,nrow=1)
   }
   
   #identify data objects
   files <- list.files(path)
   files <- files[grep(files,pattern = ".RData")] 
   
   #separate the predictors into packages according to markersPerIteration
   blocks <- seq(1,ncol(scores),markersPerIteration)
   
   #get null distributions per block and compare it to the actual scores
   pMatBlock <- lapply(1:length(blocks),FUN=function(blockN){
      if(blockN<length(blocks)){
         preds <- blocks[blockN]:(blocks[blockN+1]-1)
      }else{
         preds <- blocks[blockN]:ncol(scores)
      }
      #load the null distributions and bind them into a single matrix
      nullMat <- lapply(files,FUN=function(file){
         if("out"%in%ls()){
            rm(out)
         }
         load(paste0(path,file))
         if(!"out"%in%ls()){
            return(NULL)
         }
         out <- out[preds,]
         return(out)
      })
      nullMat <- do.call(cbind, nullMat)
      
      #if no files were available or did not contain distributions stop here
      if(is.null(dim(nullMat))){
         stop("No null distributions can be found in the chosen directory.")
      }
      
      #compute the size of the null distribution to calculate the smallest possible p-value
      nullSize <- ncol(nullMat)
      
      #compare scores to the empirical null distributions
      pMat <- lapply(1:length(preds),FUN=function(predN){
         null <- nullMat[predN,]
         lnull <- length(null)
         pPred <- sapply(scores[,preds[predN]],FUN=function(scr){
            max(sum(null>=scr)/lnull,1/lnull)
         })
         return(pPred)
      })
      pMat <- do.call(cbind,pMat)
      print(preds[length(preds)])
      return(pMat)
   })
   #bind all matrices to a single one
   pMatBlock <- do.call(cbind,pMatBlock)
   
   if(pCorrection!="none"){
      pMatBlock[1:length(pMatBlock)] <- p.adjust(pMatBlock[1:length(pMatBlock)],method = pCorrection)
   }
   
   #return results
   if(nrow(pMatBlock)==1){
      return(pMatBlock[1,])
   }else{
      return(pMatBlock)
   }
}


#' @title Group QTL
#' 
#' @description
#' Function that identifies significant marker-trait associations
#' and joins close markers to a single QTL if they are consecutive
#' and/or highly correlated and joins distant QTL to QTL-groups if
#' they contain highly correlated markers.
#'
#' @param pmat Matrix with p-values, traits (rows) X loci (columns).
#' @param sigThreshold Significance threshold to applied on the p-values.
#' @param corThreshold Threshold for absolute correlation values above
#'  which signifcant loci are treated as linked to each other. If they 
#'  are close they are joined to a single region, otherwise they are 
#'  grouped together without including inbetween loci.
#' @param distThreshold Distance threshold up until which loci containing
#'  highly correlated markers including all markers in between are linked.
#' @param genotype Allele-information for all samples (rows) and 
#'  markers (columns).
#' @param chrVec Character-vector containing the chromosome on which
#'  a given marker is located.
#'
#' @return List with QTLgroups. Each entry is a list of the target
#'  and a matrix with the QTL-regions involved.
#' @export
#'
#' @examples
QTLgrouper <- function(pmat, sigThreshold, corThreshold, 
                       distThreshold, genotype, chrVec){
   
   if(is.vector(pmat)){
      pmat <- matrix(pmat,nrow=1)
   }
   if(ncol(pmat)<ncol(genotype)){
      stop("Number of markers in pmat and genotype have to be identical.")
   }
   QTLmat <- which(pmat<=sigThreshold,arr.ind=T)
   
   #check if any significant QTL are present
   if(nrow(QTLmat)==0){
      return("no significant loci")
   }
   colnames(QTLmat) <- c("target","predictor")
   
   #square threshold and correlation matrix to account for extreme negative correlation
   corThreshold <- corThreshold^2
   corMat <- cor(genotype,use="pair")^2
   
   #join consecutive markers regulating the same trait to a single QTL
   QTLlist <- joinConsecutive(QTLmat, chrVec)
   
   #join markers to QTL if they are close and correlated
   QTLlist <- joinNear(QTLlist,corMat,corThreshold,distThreshold,chrVec)
   
   #link distant QTL that are correlated
   QTLlist <- joinCorrelated(QTLlist,corMat,corThreshold)
   
   #find the predictor with the smallest p-value in each qtl
   QTLlist <- minPV(QTLlist,pmat)
   return(QTLlist)
}

#' @title Lead Marker
#' 
#' @description
#' Identify the predictor with the smallest p-value in each QTL.
#'
#' @param QTLlist List-object with an entry for each QTL as returned
#'  by joinCorrelated.
#' @param pmat Matrix containing the p-values for each predictor (columns) 
#'  and trait (rows).
#'
#' @return Object similar to QTLlist with added information regarding the most 
#'  significant predictor in each QTL.
#' @export
#'
#' @examples
minPV <- function(QTLlist,pmat){
   out <- lapply(QTLlist,FUN=function(qtl){
      target <- qtl$target
      predictors <- apply(qtl$predictors,1,FUN=function(preds){return(preds[1]:preds[2])})
      if(is.list(predictors)){
         predictors <- do.call("c",predictors)
      }
      sigPred <- predictors[which.min(pmat[target,predictors])]
      qtl <- c(qtl,sigPred)
      names(qtl)[length(qtl)] <- "mostSignificantPredictor"
      qtl <- c(qtl,pmat[target,sigPred])
      names(qtl)[length(qtl)] <- "minP"
      return(qtl)
   })
   return(out)
}

#' @title Write .qtl file
#' 
#' @description
#' Writes a QTL-list to a .qtl file as a table.
#'
#' @param QTLlist QTL stored as a list-object as returned by QTLgrouper.
#' @param traitNames Optional; a vector with names for the traits 
#'  which correspond to the rows in the p-value matrix used for QTLgrouper.
#' @param path Complete filepath ending in .qtl, specifying where 
#'  the results should be saved.
#' @param markerPositions Optional; Matrix with three columns 
#'  specifying the chromosome, start and end for each marker.
#' @param digits Optional; p-values are formatted to scientific 
#'  convention and can be rounded according to digits.
#'
#' @return The output is written in a file, specified by path.
#' @export
#'
#' @examples
writeQTL <- function(QTLlist,traitNames=NULL,path,
                     markerPositions=NULL,digits=NULL){
   if(identical(QTLlist,"no significant loci")){
      stop("no significant loci")
   }
   metaInfo <- lapply(QTLlist,names)
   if(length(QTLlist)>1){
      allId <- sapply(1:(length(metaInfo)-1),FUN=function(x){
         identical(metaInfo[x],metaInfo[x+1])
      })
      if(!all(allId)){
         stop("corrupted QTL-list")
      }
   }
   if(!is.null(markerPositions)){
      if(!is.matrix(markerPositions)&!is.data.frame(markerPositions)){
         stop("markerPositions has to be a matrix or data frame.")
      }
      if(ncol(markerPositions)!=3){
         stop("markerPositions has to contain three columns, specifying the chromosome, start and end for each marker.")
      }
   }
   if(!is.character(path)){
      stop("Specify an exact path.")
   }
   if(!grepl(".qtl",path)){
      stop("path has to end in .qtl")
   }
   dirString <- strsplit(x = path,split = "")[[1]]
   if(grepl("/",path)){
      dir <- paste(dirString[1:tail(which(dirString=="/"),n=1)],collapse = "")
      if(!file.exists(dir)){
         stop("Directory doesn't exist.")
      }
   }
   
   col_names <- c("QTL","trait number","trait name","first marker","last marker","chromosome","first position","last position","smallest p-value","most significant predictor")
   out <- lapply(QTLlist,FUN=function(x){
      qtlMat <- matrix(NA,ncol=length(col_names),nrow=nrow(x$predictors))
      qtlMat[,1] <- 0
      qtlMat[1,1] <- 1
      qtlMat[,2] <- x$target
      if(!is.null(traitNames)){
         qtlMat[,3] <- traitNames[x$target]
      }
      qtlMat[,4:5] <- x$predictors
      if(!is.null(markerPositions)){
         qtlMat[,6] <- markerPositions[x$predictors[,1],1]
         qtlMat[,7] <- markerPositions[x$predictors[,1],2]
         qtlMat[,8] <- markerPositions[x$predictors[,2],3]
      }
      qtlMat[,9] <- format(x$minP,scientific=T,digits=digits)
      qtlMat[,10] <- x$mostSignificantPredictor
      return(qtlMat)
   })
   out <- do.call("rbind",out)
   colnames(out) <- col_names
   out[,1] <- cumsum(as.numeric(out[,1]))
   write.table(x = out,file = path,append = F,quote = F,sep = "\t",col.names = T,row.names = F)
}

#' @title Load .qtl file
#' 
#' @description
#' Reconstructs a QTLlist from a .qtl file.
#'
#' @param path The location of the .qtl file that should be loaded.
#'
#' @return A list-object as it is returned by QTLgrouper.
#' @export
#'
#' @examples
readQTL <- function(path){
   if(!is.character(path)|!grepl(".qtl",path)){
      stop("Specify a .qtl file.")
   }
   if(!file.exists(path)){
      stop("Specify an existing .qtl file.")
   }
   qtlMat <- read.table(path,header=T,sep="\t",as.is = T)
   out <- lapply(1:max(as.numeric(qtlMat[,1])),FUN=function(i){
      submat <- qtlMat[as.numeric(qtlMat[,1])==i,,drop=F]
      target <- submat[1,2]
      predictors <- submat[,4:5,drop=F]
      predictors <- as.matrix(predictors)
      rownames(predictors) <- c()
      colnames(predictors) <- c("start","end")
      mostSignificantPredictor <- submat[1,10]
      minP <- as.numeric(submat[1,9])
      return(list(target=target,predictors=predictors,mostSignificantPredictor=mostSignificantPredictor,minP=minP))
   })
   return(out)
}

#' @title emma kinship
#' 
#' @description
#' Gratefully taken from emma 1.1.2 (url:http://mouse.cs.ucla.edu/emma/news.html,license: LGPL)
#' Used to generate thepopulation structure covariates for RF QTL mapping.
emma.kinship <- function(snps, method="additive", use="all") {
   n0 <- sum(snps==0,na.rm=TRUE)
   nh <- sum(snps==0.5,na.rm=TRUE)
   n1 <- sum(snps==1,na.rm=TRUE)
   nNA <- sum(is.na(snps))
   
   stopifnot(n0+nh+n1+nNA == length(snps))
   
   if ( method == "dominant" ) {
      flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
      snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
   }
   else if ( method == "recessive" ) {
      flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
      snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
   }
   else if ( ( method == "additive" ) && ( nh > 0 ) ) {
      dsnps <- snps
      rsnps <- snps
      flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
      dsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
      flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
      rsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
      snps <- rbind(dsnps,rsnps)
   }
   
   if ( use == "all" ) {
      mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
      snps[is.na(snps)] <- mafs[is.na(snps)]
   }
   else if ( use == "complete.obs" ) {
      snps <- snps[rowSums(is.na(snps))==0,]
   }
   
   n <- ncol(snps)
   K <- matrix(nrow=n,ncol=n)
   diag(K) <- 1
   
   for(i in 2:n) {
      for(j in 1:(i-1)) {
         x <- snps[,i]*snps[,j] + (1-snps[,i])*(1-snps[,j])
         K[i,j] <- sum(x,na.rm=TRUE)/sum(!is.na(x))
         K[j,i] <- K[i,j]
      }
   }
   return(K)
}

#' @title RF epistasis detection
#' 
#' @description
#' function that detects epistasis between markers using random forest. 
#' For details see https://www.biorxiv.org/content/early/2018/06/21/353193
#' Note that all submitted markers in mappingData will be used to build the RF. 
#' The parameters markerInds1 and markerInds2, or pairMat can be used to
#' restrict interaction testing to a subset of marker pairs. If neither
#' of them are specified, all markers are tested.
#'
#' @param mappingData Dataset that was prepared with the function preMap.
#' @param markerInds1,markerInds2 The indices of the markers that should 
#'  be tested for interaction. All markers in markerInds1 will be tested 
#'  against all markers in markerInds2. Will be ignored if pairMat is submitted.
#' @param pairMat An integer matrix with two columns indicating the
#'  marker pairs that should be tested for interaction. Each row should
#'  be the indices of the two markers that should be tested. 
#' @param mtry Number of predictors (markers) that are randomly sampled
#'  as candidates at each split. Will be passed on to the Random Forest
#'  algorithm. Default is a third of the number of markers.
#' @param ntree Total number of trees. Defaults to 30000. More trees
#'  lead to higher power.
#' @param nodesize Minimum size of terminal nodes. Setting 
#'  this number larger causes smaller trees to be grown 
#'  (and thus take less time). Default is 5.
#' @param minTest Minimal number of slopes that need to be 
#'  available for a marker pair in order for the splitA test 
#'  to be applied. Default is 5 and should not be lowered for
#'  statistical reasons.
#' @param nthreads Number of threads that should be used 
#'  for building the forest.
#' @param fullTable logical; Should the output table include
#'  mean slopes, selection frequencies, and ungrouped 
#'  p-values? Normally not needed.
#' @param npermut Number of different imputations that should
#'  be created in the case of missing genotypes. ntree should
#'   be a multiple of npermut. Defaults to 100.
#' @param Rcutoff Absolute pearson correlation 
#'  coefficient R (0 > Rcutoff <= 1) that should be used as a 
#'  threshold for marker pairs to be excluded. 
#' 
#' @return Table with the marker pair indices and their 
#'  respective p-values from the pairedSF, the splitA, the selA, 
#'  and their ensemble. If there were multiple phenotypes, the
#'  results for each phenotype are returned as a list. If fullTable
#'  is TRUE, intermediate values, including the pair of mean
#'  slopes for each marker pair, the conditinal selection 
#'  frequencies, and ungrouped p-values.
#' @export
#'
#' @examples
RFepistasis <- function(mappingData, 
                        markerInds1 = NULL, markerInds2 = NULL, pairMat = NULL, 
                        mtry = NULL, ntree = 30000, npermut = 100,
                        nodesize = 5, minTest = 5, nthreads = 1,fullTable = F,
                        Rcutoff = 0.9){
   RFE = require(RandomForestExtended)
   if(!RFE){
      stop("RFepistasis depends on the package RandomForestExtended. 
           It can be accessed at http://cellnet-sb.cecad.uni-koeln.de/resources/RandomForestExtended.")
   }
   missingG = any(is.na(mappingData$genotype))
   # check input
   if(any(!c("genotype","genotype2group","group2genotype","mappingCovariates",
             "NAlist","phenotype")%in%names(mappingData))){
      stop("mappingData has to be generated with preMap.")
   }
   if(!is.null(mtry) & !is.numeric(mtry)){
      stop("mtry has to be either NULL or a positive integer.")
   }
   if(is.numeric(mtry)){
      if(mtry < 0){
         stop("mtry has to be either NULL or a positive integer.")
      }
   } 
   if(ntree<1 | !is.numeric(ntree)){
      stop("ntree has to be a positive integer.")
   }
   if(ntree<100 & is.numeric(ntree)){
      warning("number of trees is likely too small (<100)")
   }
   if(minTest<1 | !is.numeric(minTest)){
      stop("minTest has to be a positive integer.")
   }
   if(nodesize<1 | !is.numeric(nodesize)){
      stop("nodesize has to be a positive integer.")
   }
   if(nthreads<1 | !is.numeric(nthreads)){
      stop("nthreads has to be a positive integer.")
   }
   if(!is.logical(fullTable)){
      stop("fullTable has to be either TRUE or FALSE")
   }
   if(missingG){
      if(npermut < 1 | !is.numeric(npermut))
         stop("there are missing genotypes, but npermut is not a positive integer.")
      if(ntree %% npermut != 0)
         warning("there are missing genotypes, but ntree is not a multiple of npermut:\n final RF size may differ slightly from ntree")
   }
   if(Rcutoff<0 | Rcutoff>1 | !is.numeric(Rcutoff)){
      stop("Rcutoff has to be a number between 0 and 1.")
   }
   if(!is.null(markerInds1) & !is.integer(markerInds1)){
      stop("markerInds1 has to be either NULL or an integer vector")
   }
   if(!is.null(markerInds2) & !is.integer(markerInds2)){
      stop("markerInds2 has to be either NULL or an integer vector")
   }
   if(!is.null(pairMat) & !is.matrix(pairMat)){
      stop("pairMat has to be either NULL or an integer matrix with two columns")
   }
   if(is.null(markerInds1) & is.null(markerInds2) & is.null(pairMat)){
      message("no markerInds or pairMat submitted --> using all markers")
   }
   if(is.matrix(pairMat)){
      if(!is.integer(pairMat) & !ncol(pairMat) == 2){
         stop("the pairMat matrix has to be an integer matrix with two columns")
      }
   } 
   
   # initialize some variables
   genotype = mappingData$genotype
   phenotypes = mappingData$phenotype
   if(is.vector(phenotypes))
      phenotypes = matrix(phenotypes, nrow = 1)
   popStr = mappingData$mappingCovariates
   NAlist = mappingData$NAlist
   
   if(is.null(pairMat)){
      if(is.null(markerInds1)){
         markerInds1 = 1:ncol(genotype)
      } 
      if(is.null(markerInds2)){
         markerInds2 = 1:ncol(genotype)
      } 
      #make a matrix with unique marker combinations
      pairMat <- expand.grid(markerInds1, markerInds2)
      pairMat <- t(apply(pairMat,1,function(x){ c(min(x), max(x)) }))
      pairMat <- pairMat[(pairMat[,1] != pairMat[,2]),]
      pairMat <- unique(pairMat)
   }else{
      markerInds1 = unique(pairMat[,1])
      markerInds2 = unique(pairMat[,2])
   }
   interestingMarkers <- union(markerInds1,markerInds2)
   interestingMarkers <- sort(interestingMarkers)
   cors <- cor(genotype[,interestingMarkers], use = "pairwise.complete.obs") 
   ngenotype <- ncol(genotype)
   if(is.null(mtry)){
      mtry = floor(ncol(genotype)/3)
   }
   mHash <- rep(NA,ngenotype) # create a Hash in order to be able to re-map the original indices of the alleles
   ind <- 1
   for(i in 1:ngenotype){
      if(i %in% interestingMarkers){
         mHash[i] <- ind
         ind <- ind+1
      }
   }
   # rowHash <- matrix(NA, ncol = ngenotype, nrow = ngenotype)
   # for(i in 1:nrow(pairMat)){
   #    rowHash[pairMat[i,1],pairMat[i,2]] <- i
   #    rowHash[pairMat[i,2],pairMat[i,1]] <- i
   # }
   
   # iterate through the traits and call interactions
   phenoRes = lapply(seq_len(nrow(phenotypes)),function(p){
      phenotype = phenotypes[p,]
      # remove samples with missing phenotypes
      missingP <- !is.na(phenotype)
      phenotype <- phenotype[missingP]
      
      # grow RF
      if(missingG){
         #impute missing genotypes and grow the forest 
         rf = lapply(1:npermut, FUN = function(x){
            imputed = replaceGenoNAs(genotype, NAlist)
            imputed = imputed[missingP,]
            rf <- RandomForestExtended::randomForest(x=imputed,y=phenotype,
                                                     ntree=ceiling(ntree/npermut),mtry=mtry,
                                                     nodesize=nodesize,keep.forest=T,
                                                     importance=F, nthreads=nthreads)
         })
         rf = combineForestsExtended(rf)
      } else{
         mappingGenotype = genotype[missingP,]
         rf <- RandomForestExtended::randomForest(x=mappingGenotype,y=phenotype,
                                                  ntree=ntree,mtry=mtry,
                                                  nodesize=nodesize,keep.forest=T,
                                                  importance=F, nthreads=nthreads)
      }
      ntree = rf$ntree 
      
      #organize the forest-matrices into an array
      forestArray <- array(0,dim = c(dim(rf$forest$nodestatus),5))
      forestArray[,,1] <- rf$forest$leftDaughter # rows: node, cols: tree
      forestArray[,,2] <- rf$forest$rightDaughter # rows: node, cols: tree
      forestArray[,,3] <- rf$forest$nodestatus # whether it is a terminal node or not (-3 means not terminal, -1 terminal, and 0 not used)
      forestArray[,,4] <- rf$forest$bestvar # index of the alleles that were used for the split at each node (and in each tree)
      forestArray[,,5] <- rf$forest$oobpred
      rm(rf)
      
      #create an object to store the distributions of side specific differences in means
      #maybe use sparse matrix representation here? for example library('Matrix') or library('slam')
      mdim <- length(interestingMarkers)
      left <- array(dim=c(mdim,mdim,3)) # dimensions: sum, sumofsquares, count
      left[,,3] <- 0
      right <- array(dim=c(mdim,mdim,3))
      right[,,3] <- 0
      counts <- array(0,dim=c(mdim,mdim,2))
      
      # walk through the forest and collect the distribution of side-specific means
      for (i in 1:ntree){
         tree <- forestArray[,i,]
         dim(tree) = c(length(tree)/5, 5)
         tS <- treeStructure(tree)
         ## for the paired frequency test:
         markers <- tree[,4]
         markers <- unique(markers[markers!=0]) # get all the markers that were used in this tree
         markers <- markers[markers %in% interestingMarkers]
         markers <- mHash[markers]
         counts[markers,markers,1] = counts[markers,markers,1] + 1 # increase the counter for the marker pairs that occured together
         counts[-markers,markers,2] = counts[-markers,markers,2] + 1 # increase the counter for the cases where only one of the markers was used (not used in row)
         
         ## left side:
         #build a matrix: col 1&2: nodes, col3&4: respective markers
         pairs <- cbind(tS$nl, matrix(mHash[tree[tS$nl,4]], ncol = 2,byrow = F))
         pairs = pairs[apply(pairs,1,function(x) !any(is.na(x))),,drop = F] # exclude pairs where one of the markers is not in interestingMarkers  
         # for asymmetry tests
         if(nrow(pairs) > 0){
            for(j in 1:nrow(pairs)){
               stored <- left[pairs[j,3],pairs[j,4],3]
               childMeans <- tree[tree[pairs[j,2],1:2],5] 
               if(any(childMeans==0)){ # if there are no oob samples for this node, the oob prediction is zero. Here these cases are excluded.
                  next
               }
               if(stored >= 1){
                  left[pairs[j,3],pairs[j,4],1] <- left[pairs[j,3],pairs[j,4],1] + 
                     (childMeans[2] - childMeans[1])
                  left[pairs[j,3],pairs[j,4],2] <- left[pairs[j,3],pairs[j,4],2] + 
                     (childMeans[2] - childMeans[1]) * (childMeans[2] - childMeans[1])
               } else {
                  left[pairs[j,3],pairs[j,4],1] <- childMeans[2] - childMeans[1]
                  left[pairs[j,3],pairs[j,4],2] <- (childMeans[2] - childMeans[1]) * (childMeans[2] - childMeans[1])
               }
               left[pairs[j,3],pairs[j,4],3] <- stored + 1
            }
         }
         
         ## right side: (same stuff)
         pairs <- cbind(tS$nr, matrix(mHash[tree[tS$nr,4]], ncol = 2, byrow = F))
         pairs = pairs[apply(pairs,1,function(x) !any(is.na(x))),,drop = F] # exclude pairs where one of the markers is not in interestingMarkers
         if(nrow(pairs) > 0){
            for(j in 1:nrow(pairs)){
               stored <- right[pairs[j,3],pairs[j,4],3]
               childMeans <- tree[tree[pairs[j,2],1:2],5] 
               if(any(childMeans == 0)){
                  next
               }
               if(stored >= 1){
                  right[pairs[j,3],pairs[j,4],1] <- right[pairs[j,3],pairs[j,4],1] + 
                     (childMeans[2] - childMeans[1])
                  right[pairs[j,3],pairs[j,4],2] <- right[pairs[j,3],pairs[j,4],2] + 
                     (childMeans[2] - childMeans[1]) * (childMeans[2] - childMeans[1])
               } else{
                  right[pairs[j,3],pairs[j,4],1] <- (childMeans[2] - childMeans[1])
                  right[pairs[j,3],pairs[j,4],2] <- (childMeans[2] - childMeans[1]) * (childMeans[2] - childMeans[1])
               }
               right[pairs[j,3],pairs[j,4],3] <- stored+1
            }
         }
         
      }
      # get indices of marker pairs which have occured more than minTest times in the tree 
      toTest <- which((left[,,3] >= minTest) | (right[,,3] >= minTest), arr.ind = T)
      temp = matrix(interestingMarkers[toTest],ncol = 2)
      # toTest = toTest[!is.na(rowHash[temp]),]
      
      nTrialsLeft <- rowSums(left[,,3]) # this is needed for the binom test
      nTrialsRight <- rowSums(right[,,3]) # for each marker, how many splits were there on the left/right side in total
      
      if(nrow(toTest) < 1){
         warning("no tests possible for this phenotype")
         return(NULL)
      } else {  # perform ttest, frequency test, and binomial test
         toTest <- cbind(toTest, NA, NA, NA, NA, NA, NA)
         toTest[,3:8] <- t(apply(toTest[,1:2], 1, FUN = function(row){
            if (abs(cors[row[1], row[2]]) >= Rcutoff) { # exclude markers in LD
               return(c(1, 0, 1, 0, 0, 1))
            }
            
            # paired frequency test:
            freqmat <- c(counts[row[1],row[2],1],
                         counts[row[1],row[2],2], 
                         counts[row[2],row[1],2],
                         0)
            freqmat[4] <- ntree - sum(freqmat)
            dim(freqmat) <- c(2, 2)
            freq_res <- fisher.test(freqmat, alternative = "greater")$p.value
            
            # prepare the counts for sel. asymmetry:
            nleft <- left[row[1],row[2],3]
            nright <- right[row[1],row[2],3]
            
            # binomial test:
            ntrialsL <- nTrialsLeft[row[1]] 
            ntrialsR <- nTrialsRight[row[1]]
            
            if(ntrialsL == 0 | ntrialsR == 0){
               return(c(NA, NA, NA, NA, NA, freq_res))
            }
            
            # t-test:
            if((nright >= minTest) & (nleft >= minTest)){ # we only want to compute the ttest if we have at least 5 slopes on the right and left side.
               n1 = left[row[1],row[2],3]
               m1 = left[row[1],row[2],1] / n1
               v1 = (left[row[1],row[2],2] / n1) - (m1 * m1)
               n2 = right[row[1],row[2],3]
               m2 = right[row[1],row[2],1] / n2
               v2 = (right[row[1],row[2],2] / n2) - (m2 * m2)
               se = sqrt((((n1 - 1) * v1) + ((n2 - 1) * v2)) / (n1 + n2 - 2))
               t  = sqrt(n1 * n2 / (n1 + n2)) * ((m1 - m2) / se)
               df <- n1 + n2 - 2
               tt_res = c(2 * pt(-abs(t), df = df), m1 - m2)
            } else 
               tt_res <- c(NA,NA)
            
            chi_res <- tryCatch(unlist(prop.test(x = c(nleft,nright), 
                                                 n = c(ntrialsL,ntrialsR), 
                                                 correct = F, 
                                                 alternative = "two.sided")[3:4], 
                                       use.names = FALSE),
                                warning=function(w) c(NA, NA, NA))
            
            return(c(tt_res, chi_res, freq_res))
         }))
      }
      
      # put all results in one table:
      toTest <- toTest[as.logical(rowSums(is.finite(toTest[,3:8]))),] # remove rows where no tests were performed
      toTest[,1] = interestingMarkers[toTest[,1]]
      toTest[,2] = interestingMarkers[toTest[,2]]
      
      pairMat <- cbind(pairMat, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      colnames(pairMat) <- c("markerA", "markerB",
                             "p_splitA_AbeforeB", "slope_diff_AbeforeB", 
                             "p_splitA_BbeforeA", "slope_diff_BbeforeA", 
                             "p_selA_AbeforeB", "prop_BleftOfA", "prop_BrightOfA", 
                             "p_selA_BbeforeA", "prop_AleftOfB", "prop_ArightOfB", 
                             "p_pairedSF_AbeforeB","p_pairedSF_BbeforeA")
      for(i in seq_len(nrow(toTest))){
         mA = toTest[i,1]
         mB = toTest[i,2]
         pairMatInd = which(min(mA,mB) == pairMat[,1] & max(mA,mB) == pairMat[,2])
         if(mA < mB){
            pairMat[pairMatInd,c(3,4,7,8,9,13)] <- toTest[i,3:8]
            # pairMat[rowHash[toTest[i,1],toTest[i,2]],3] <- toTest[i,3]
            # pairMat[rowHash[toTest[i,1],toTest[i,2]],4] <- toTest[i,4]
            # pairMat[rowHash[toTest[i,1],toTest[i,2]],7] <- toTest[i,5]
            # pairMat[rowHash[toTest[i,1],toTest[i,2]],8] <- toTest[i,6]
            # pairMat[rowHash[toTest[i,1],toTest[i,2]],9] <- toTest[i,7]
            # pairMat[rowHash[toTest[i,1],toTest[i,2]],13] <- toTest[i,8]
         }else{
            pairMat[pairMatInd,c(5,6,10,11,12,14)] <- toTest[i,3:8]
            # pairMat[rowHash[toTest[i,2],toTest[i,1]],5] <- toTest[i,3]
            # pairMat[rowHash[toTest[i,2],toTest[i,1]],6] <- toTest[i,4]
            # pairMat[rowHash[toTest[i,2],toTest[i,1]],10] <- toTest[i,5]
            # pairMat[rowHash[toTest[i,2],toTest[i,1]],11] <- toTest[i,6]
            # pairMat[rowHash[toTest[i,2],toTest[i,1]],12] <- toTest[i,7]
            # pairMat[rowHash[toTest[i,2],toTest[i,1]],14] <- toTest[i,8]
         }
      }
      res <- cbind(pairMat,
                   "splitA" = apply(pairMat[,c(3,5)],1,combinePvalues) ,
                   "selA" = apply(pairMat[,c(7,10)],1,combinePvalues),
                   "pairedSF" = apply(X = pairMat[,13:14],MARGIN = 1,FUN = function(x){
                      bool <- is.na(x)
                      ifelse(any(bool), yes = x[!bool], no = x[1])
                   }))
      res <- cbind(res,"ensemble" = apply(res[,15:17], 1, combinePvalues))
      if(fullTable)
         return(res)
      else
         return(res[,c(1:2,15:18)])
   })
   return(phenoRes)
}

#' @title Reformat RF tree structure
#' 
#' @description
#' Function that reformats the structure of a RF tree for use in RFepistasis.
#'
#' @param treeMat Matrix representing the structure of a tree of a forest.
#'  Each row represents one node. The first two columns are the indices of 
#'  the rows of the left and right child nodes, respectively. 
#'  The third column indicates whether it is a terminal node 
#'  or not (-3 means not terminal, -1 terminal, and 0 not used). 
#'  The fourth column is the index in the genotype of the marker 
#'  that was used for the split at this node. The last column is 
#'  the mean out-of-bag (oob) trait value for this node.
#' 
#' @return List of length 2. $nl: matrix with two columns containing
#'  indices of predictors. Predictors in the second column were used
#'  somewhere in the tree on the left side of the respective predictor
#'  in clumn a. $nr: same as $nl, but for pairs where predictors in
#'  column two were used somewhere on the right side of the predictors in column 1.
#' @export
#'
#' @examples
treeStructure <- function(treeMat){
   rownames(treeMat)=1:nrow(treeMat)
   treeMat <- treeMat[treeMat[,3] == -3,, drop = F] # get only nodes which have children
   L <- vector("list",nrow(treeMat))
   names(L) <- rownames(treeMat)
   R <- L
   nodes <- rev(rownames(treeMat))
   for(node in nodes){ # go through nodes and store all non-terminal left/right children in L/R
      dirChildren <- as.character.default(treeMat[node,1:2])
      logic <- (dirChildren %in% nodes)
      if(logic[1]){
         L[[node]] <- c(as.integer(dirChildren[1]), L[[dirChildren[1]]], R[[dirChildren[1]]])
      }
      if(logic[2]){
         R[[node]] <- c(as.integer(dirChildren[2]), L[[dirChildren[2]]], R[[dirChildren[2]]])
      } 
   }
   
   ## transform above lists into matrices with valid marker combinations, first column mother, second column dauther node
   nl <- c(unlist(sapply(names(L), FUN = function(x){
      rep(as.numeric(x), length(L[[x]]))
   },USE.NAMES = F)), unlist(L, use.names = FALSE))
   if(is.null(nl))
      nl = matrix(data = NA, nrow = 0, ncol = 2)
   else
      dim(nl) <- c(length(nl)/2, 2)
   nr <- c(unlist(sapply(names(R), FUN = function(x){
      rep(as.numeric(x), length(R[[x]]))
   },USE.NAMES = F)), unlist(R, use.names = FALSE))
   if(is.null(nr))
      nr = matrix(data = NA, nrow = 0, ncol = 2)
   else
      dim(nr) <- c(length(nr)/2, 2)
   
   return(list(nl = nl, nr = nr))
}

#' @title Fisher combine p-values
#' 
#' @description
#' Ffunction that combines p-values using the fisher method. 
#' Ignores NA in the input. Returns NA is no valid p-values were
#' in the input (e.g., only NAs). If only one valid p-value is
#' submitted, this p-value is returned unchanged. Used in RFepistasis.
#'
#' @param p Numeric vector of p-values that should be combined
#' 
#' @return Combined p-value.
#' @export
#'
#' @examples
combinePvalues <- function(p){
   if(!is.numeric(p) | length(p)<1){
      stop("p has to be a numeric vector.")
   }
   pvals = p[!is.na(p)]
   validpvals = length(pvals)
   if(validpvals==0){ # there are no pvalues --> return NA
      return(NA)
   }
   if(validpvals==1){ # there's only one pvalue --> return it
      return(pvals)
   }else{ # there are several pvalues --> combine using fisher's method
      pcomb <- pchisq( -2*sum(log(pvals)) , df = 2*validpvals , lower.tail=FALSE)
      return(pcomb)
   }
}


#' @title Combine RFs
#' 
#' @description
#' Function that combines several randomForest objects 
#' to one big Random Forest. Workaround for a bug in the
#' function  \code{combine} in \code{RandomForestExtended}.
#'
#' @param ... Two or more objects of class randomForest, to be combined into one.
#' 
#' @return An object of class randomForest.
#' @export
#'
#' @examples
combineForestsExtended = function (...) {
   pad0 <- function(x, len) c(x, rep(0, len - length(x)))
   padm0 <- function(x, len) rbind(x, matrix(0, nrow = len - 
                                                nrow(x), ncol = ncol(x)))
   rflist <- list(...)
   areForest <- sapply(rflist, function(x) inherits(x, "randomForest"))
   if (any(!areForest)) 
      rflist <- rflist[[1]]
   areForest <- sapply(rflist, function(x) inherits(x, "randomForest"))
   if (any(!areForest)) 
      stop("Argument must be a list of randomForest objects")
   rf <- rflist[[1]]
   classRF <- rf$type == "classification"
   trees <- sapply(rflist, function(x) x$ntree)
   ntree <- sum(trees)
   rf$ntree <- ntree
   nforest <- length(rflist)
   haveTest <- !any(sapply(rflist, function(x) is.null(x$test)))
   vlist <- lapply(rflist, function(x) rownames(RandomForestExtended::importance(x)))
   numvars <- sapply(vlist, length)
   if (!all(numvars[1] == numvars[-1])) 
      stop("Unequal number of predictor variables in the randomForest objects.")
   for (i in seq_along(vlist)) {
      if (!all(vlist[[i]] == vlist[[1]])) 
         stop("Predictor variables are different in the randomForest objects.")
   }
   haveForest <- sapply(rflist, function(x) !is.null(x$forest))
   if (all(haveForest)) {
      nrnodes <- max(sapply(rflist, function(x) x$forest$nrnodes))
      rf$forest$nrnodes <- nrnodes
      rf$forest$ndbigtree <- unlist(sapply(rflist, function(x) x$forest$ndbigtree))
      rf$forest$nodestatus <- do.call("cbind", lapply(rflist, 
                                                      function(x) padm0(x$forest$nodestatus, nrnodes)))
      rf$forest$bestvar <- do.call("cbind", lapply(rflist, 
                                                   function(x) padm0(x$forest$bestvar, nrnodes)))
      rf$forest$xbestsplit <- do.call("cbind", lapply(rflist, 
                                                      function(x) padm0(x$forest$xbestsplit, nrnodes)))
      rf$forest$nodepred <- do.call("cbind", lapply(rflist, 
                                                    function(x) padm0(x$forest$nodepred, nrnodes)))
      tree.dim <- dim(rf$forest$treemap)
      if (classRF) {
         rf$forest$treemap <- array(unlist(lapply(rflist, 
                                                  function(x) apply(x$forest$treemap, 2:3, pad0, 
                                                                    nrnodes))), c(nrnodes, 2, ntree))
      }
      else {
         rf$forest$leftDaughter <- do.call("cbind", lapply(rflist, 
                                                           function(x) padm0(x$forest$leftDaughter, nrnodes)))
         rf$forest$rightDaughter <- do.call("cbind", lapply(rflist, 
                                                            function(x) padm0(x$forest$rightDaughter, nrnodes)))
         rf$forest$nodepops <- do.call("cbind", lapply(rflist, 
                                                       function(x) padm0(x$forest$nodepops, nrnodes)))
         rf$forest$oobpred <- do.call("cbind", lapply(rflist, 
                                                      function(x) padm0(x$forest$oobpred, nrnodes)))
         rf$forest$oobmse <- do.call("cbind", lapply(rflist, function(x) padm0(x$forest$oobmse, 
                                                                               nrnodes)))
      }
      rf$forest$ntree <- ntree
      if (classRF) 
         rf$forest$cutoff <- rflist[[1]]$forest$cutoff
   }
   else {
      rf$forest <- NULL
   }
   if (classRF) {
      rf$votes <- 0
      rf$oob.times <- 0
      areVotes <- all(sapply(rflist, function(x) any(x$votes > 
                                                        1, na.rf = TRUE)))
      if (areVotes) {
         for (i in 1:nforest) {
            rf$oob.times <- rf$oob.times + rflist[[i]]$oob.times
            rf$votes <- rf$votes + ifelse(is.na(rflist[[i]]$votes), 
                                          0, rflist[[i]]$votes)
         }
      }
      else {
         for (i in 1:nforest) {
            rf$oob.times <- rf$oob.times + rflist[[i]]$oob.times
            rf$votes <- rf$votes + ifelse(is.na(rflist[[i]]$votes), 
                                          0, rflist[[i]]$votes) * rflist[[i]]$oob.times
         }
         rf$votes <- rf$votes/rf$oob.times
      }
      rf$predicted <- factor(colnames(rf$votes)[max.col(rf$votes)], 
                             levels = levels(rf$predicted))
      if (haveTest) {
         rf$test$votes <- 0
         if (any(rf$test$votes > 1)) {
            for (i in 1:nforest) rf$test$votes <- rf$test$votes + 
                  rflist[[i]]$test$votes
         }
         else {
            for (i in 1:nforest) rf$test$votes <- rf$test$votes + 
                  rflist[[i]]$test$votes * rflist[[i]]$ntree
         }
         rf$test$predicted <- factor(colnames(rf$test$votes)[max.col(rf$test$votes)], 
                                     levels = levels(rf$test$predicted))
      }
   }
   else {
      rf$predicted <- 0
      for (i in 1:nforest) rf$predicted <- rf$predicted + rflist[[i]]$predicted * 
            rflist[[i]]$ntree
      rf$predicted <- rf$predicted/ntree
      if (haveTest) {
         rf$test$predicted <- 0
         for (i in 1:nforest) rf$test$predicted <- rf$test$predicted + 
               rflist[[i]]$test$predicted * rflist[[i]]$ntree
         rf$test$predicted <- rf$test$predicted/ntree
      }
   }
   have.imp <- !any(sapply(rflist, function(x) is.null(x$importance)))
   if (have.imp) {
      rf$importance <- rf$importanceSD <- 0
      for (i in 1:nforest) {
         rf$importance <- rf$importance + rflist[[i]]$importance * 
            rflist[[i]]$ntree
         rf$importanceSD <- rf$importanceSD + rflist[[i]]$importanceSD^2 * 
            rflist[[i]]$ntree
      }
      rf$importance <- rf$importance/ntree
      rf$importanceSD <- sqrt(rf$importanceSD/ntree)
      haveCaseImp <- !any(sapply(rflist, function(x) is.null(x$localImportance)))
      if (haveCaseImp) {
         rf$localImportance <- 0
         for (i in 1:nforest) {
            rf$localImportance <- rf$localImportance + rflist[[i]]$localImportance * 
               rflist[[i]]$ntree
         }
         rf$localImportance <- rf$localImportance/ntree
      }
   }
   have.prox <- !any(sapply(rflist, function(x) is.null(x$proximity)))
   if (have.prox) {
      rf$proximity <- 0
      for (i in 1:nforest) rf$proximity <- rf$proximity + rflist[[i]]$proximity * 
            rflist[[i]]$ntree
      rf$proximity <- rf$proximity/ntree
   }
   if (classRF) {
      rf$confusion <- NULL
      rf$err.rate <- NULL
      if (haveTest) {
         rf$test$confusion <- NULL
         rf$err.rate <- NULL
      }
   }
   else {
      rf$mse <- rf$rsq <- NULL
      if (haveTest) 
         rf$test$mse <- rf$test$rsq <- NULL
   }
   keep.inbag <- !is.null(rf$inbag)
   if (keep.inbag) {
      for (i in 2:length(rflist)) {
         rf$inbag <- cbind(rf$inbag, rflist[[i]]$inbag)
      }
   }
   else {
      rf$inbag <- NULL
   }
   rf
}
