 \name{Generate_Meta_Files}
 \alias{Generate_Meta_Files}
 \title{Generate summary statistic files}
 \description{
    Generate Meta SSD (MSSD) and Meta Info (MInfo) files. 
    Both files are needed to run MetaSKAT with summary statistics.
 }
 \usage{

	Generate_Meta_Files(obj, File.Bed, File.Bim, File.SetID, File.MSSD
	, File.MInfo, N.Sample, data=NULL)

 }
 \arguments{
    \item{obj}{returned object from SKAT_Null_Model. }
    \item{File.Bed}{name of the binary ped file (BED).}
    \item{File.Bim}{name of the binary bim file (BIM).}
    \item{File.SetID}{name of the SNP set ID file. The first column must be Set ID, and the second column must be SNP ID. There should be no header!!}
    \item{File.MSSD}{name of MSSD file that will be generated.}
    \item{File.MInfo}{name of MInfo file that will be generated.}
    \item{N.Sample}{number of samples.}
    \item{data}{an optional data frame containing the variables in the model (default=NULL).
If it is NULL, the variables are taken from environment(formula)}
}
\details{
	It generates summary statistic files (MSSD and MInfo files) from plink formated data files. 
	To run meta analysis, each study should provide both MSSD and MInfo files. 
	The MSSD is a binary file with between-SNP information matrices, 
	and MInfo is a text formated file with information on cohorts and SNPsets. 
	
}


\author{Seunggeun Lee}
