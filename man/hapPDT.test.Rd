\name{hapPDT.test}
\alias{hapPDT.test}
\title{Haplotype-based pedigree disequilibrium test for rare variants in pedigrees.}
\description{
To detecting rare variants for binary traits using general pedigrees, the pedigree disequilibrium tests are proposed by collapsing rare haplotypes with/without weights.}
\usage{ hapPDT.test(preHapList,aff=2,unaff=1, trace=FALSE)
}
    
   

\arguments{   
\item{preHapList}{output from whap.prehap function}  
\item{aff}{indicates the values that represents affected status in preHapList data; default is 2.}
\item{unaff}{indicates the values that represents unaffected status in preHapList data; default is 1.} 
\item{trace}{indicates whether or not the intermediate outcomes should be printed; default is FALSE.}  
}

  
       
\value{
\item{TDT}{ Transmission/disequilibrium matrix for each pedigrees.}   
\item{Sib}{ Discordant sib pairs matrix for each pedigrees.}   
\item{PDT}{ Pedigree disequilibrium matrix for each pedigrees, which is the sum of TDT and Sib.}    
\item{W}{ Weights used in Weighted haplotype PDT test.}              
\item{test.v1}{ Weighted haplotype PDT test statistic with weights W.}    
\item{test.v0}{ Unweighted haplotype PDT test statistic with weights=1.}    
\item{pvalue.v1}{ P value of Weighted haplotype PDT test statistic (test.v1).}    
\item{pvalue.v0}{ P value of Unweighted haplotype PDT test statistic (test.v0).}    

 

}

 
\references{    

Guo W , Shugart YY, Does Haplotype-based Collapsing Tests Gain More Power than Variant-based Collapsing Tests for Detecting Rare Variants in Pedigrees  
(manuscript). 
 
} 

 
 