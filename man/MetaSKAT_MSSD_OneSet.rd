 \name{MetaSKAT_MSSD_OneSet}
 \alias{MetaSKAT_MSSD_OneSet}
 \alias{MetaSKAT_MSSD_ALL}
 \title{Meta analysis SKAT with summary data from each study cohort.}
 \description{
    Meta analysis SKAT with Meta SSD (MSSD) and Info files. MetaSKAT_MSSD_OneSet computes a p-value of a given set,
    and MetaSKAT_MSSD_ALL computes p-values of all sets.
 }
 \usage{

MetaSKAT_MSSD_OneSet(Cohort.Info, SetID, combined.weight=TRUE, weights.beta=c(1,25),
method="davies", r.corr=0, is.separate = FALSE, Group_Idx=NULL, MAF.cutoff=1)

MetaSKAT_MSSD_ALL(Cohort.Info, ...)

 }
\arguments{
    \item{Cohort.Info}{output object from Open_MSSD_File_2Read function. }
	\item{SetID}{a character value of set id to test.} 
    \item{combined.weight}{a logical value (default=TRUE) for the type of weighting. 
	If it is TRUE, a weight for each SNP is computed using MAFs that are common across studies. 
	If it is FALSE, cohort specific weights will be used based on study specific MAFs.} 
     \item{weights.beta}{a numeric vector of parameters of beta weights (default=c(1,25))}
     \item{method}{a method to compute a p-value (default= "davies"). 
      "davies" represents an exact method that  computes the p-value by inverting the characteristic function of the mixture chisq for SKAT statistic,
       and "optimal" represents the optimal test (SKAT-O) that is based on an optimal linear combination  of burden and SKAT statistics.}
      \item{r.corr}{the \eqn{\rho} parameter of new class of kernels with compound symmetric correlation structure 
      for genotype effects (default= 0). If r.corr=0, it does the SKAT test.  If r.corr=1, it does the burden test. 
      If r.corr=a vector of grid values between 0 and 1,   it does SKAT-O. See the manual of SKAT.}
      \item{is.separate}{a logical value (default=FALSE) for homogeneous(=FALSE) or heterogeneous(=TRUE) genetic effects
      of a SNP set across studies.
      When FALSE, it is assumed that all studies share the same causal variants with the same effect size. 
      When TRUE, it is assumed that studies/groups may have different causal variants. }  
      \item{Group_Idx}{a vector of group indicator (default=NULL). 
      If a vector of integers are specified, it assumes causal variants are the same for studies with the same group index, and different for studies  with different group indexes.
      When NULL, studies are assumed to be in different groups with different group indexes.
      When is.separate=FALSE, it will be ignored.}
      \item{MAF.cutoff}{a cutoff of the MAFs of SNPs (default=1). Any SNPs with MAFs > MAF.cutoff will be excluded from the analysis.}
      \item{...}{the same parameters of MetaSKAT_MSSD_OneSet after SetID.}
}
\value{
	MetaSKAT_MSSD_OneSet and MetaSKAT_wZ return the same object. See MetaSKAT_wZ for details.
	MetaSKAT_MSSD_ALL returns a dataframe with SetIDs (first column) and p-values (second column). 
}

\details{
	Please see MetaSKAT_wZ for details.

}


\author{Seunggeun Lee}
