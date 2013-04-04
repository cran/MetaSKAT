 \name{Meta_Null_Model}
 \alias{Meta_Null_Model}
 \title{Get parameters and residuals from the H0 model}
 \description{
     Compute model parameters and residuals for MetaSKAT. 
     It only can be used when individual level data are available.
 }
 \usage{

Meta_Null_Model(y.list, x.list, n.cohort, out_type="C",  n.Resampling=0)
 }
\arguments{
      \item{y.list}{ a list object for phenotypes. Each element should be a vector of phenotypes. If you have 3 cohorts, it should have 3 elements.}
      \item{x.list}{a list object for covariates. Each element should be a vector or a  matrix of covariates. If you have 3 cohorts, it should have 3 elements. If you don't have any covariates in the particular cohort,  please set the element as "intercept". See the examples.}
      \item{n.cohort}{a numeric value of the number of cohort.}
      \item{out_type}{an indicator of the outcome type. "C" for the continuous outcome and "D" for the dichotomous outcome.} 
      \item{n.Resampling}{internal use only.}     
}
\value{
	It returns an object that has model parameters and residuals of the NULL model of no association between genetic variables and outcome phenotypes. 
}


\author{Seunggeun Lee}

\examples{


data(Example)
attach(Example)

#############################################################
#	Compute a p-value of the first gene

obj<-Meta_Null_Model(y.list, x.list, n.cohort=3, out_type="D")
MetaSKAT_wZ(Z.list[[1]], obj)$p.value

#############################################################
#	If you want to use the intercept-only model for the 2nd cohort
x.list[[2]]<-"intercept"
obj<-Meta_Null_Model(y.list, x.list, n.cohort=3, out_type="D")
MetaSKAT_wZ(Z.list[[1]], obj)$p.value


}


