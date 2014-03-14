\name{rvPDT.test.permu}
\alias{rvPDT.test.permu}
\title{Variants-based pedigree disequilibrium test for rare variants in pedigrees.}
\description{
To detecting rare variants for binary traits using general pedigrees, the pedigree disequilibrium tests are extended by collapsing rare variants with/without weights.}
\usage{rvPDT.test.permu(ped, aff=2,unaff=1, snpCol, hfreq=NULL,
useFamWeight=TRUE, nperm=1000,trace=FALSE)    
}  
   

\arguments{   
\item{ped}{input data, has same format with PLINK but having column names. The PED file is a white-space (space or tab) delimited file: the first six columns are mandatory: 
     FID: Family ID;
     IID: Individual ID;
     FA: Paternal ID;
     MO: Maternal ID;
     SEX: Sex (1=male; 2=female; other=unknown);
     PHENO: Phenotype; 
Genotypes (column 7 onwards) should also be white-space delimited; they are coded as 0, 1 and 2, indicating the number of coding allele, and NA is for missing genotype. } 
\item{aff}{indicates the values that represents affected status in ped data; default is 2.}
\item{unaff}{indicates the values that represents unaffected status in ped data; default is 1.} 
\item{snpCol}{indicates the columns of variants in ped data.}
\item{hfreq}{indicates the frequencies of variants that used in calculating weights; when it is NULL, the frequencies are estimated by ped data.} 
\item{useFamWeight}{indicates whether the family weights need to be used in the test.}    
\item{nperm}{indicates the times of permutation; default is 1000.}    
\item{trace}{indicates wether or not the intermediate outcomes should be printed; default is FALSE.}  
}

  
       
\value{
\item{TDT}{ Transmission/disequilibrium matrix for each pedigrees.}   
\item{Sib}{ Discordant sib pairs matrix for each pedigrees.}   
\item{PDT}{ Pedigree disequilibrium matrix for each pedigrees, which is the sum of TDT and Sib.}    
\item{W}{ Weights used in Weighted rvPDT test.}              
\item{test.v1}{ Weighted rvPDT test statistic with weights W.}    
\item{test.v0}{ Unweighted rvPDT test statistic with weights=1.}    
\item{pvalue.v1}{ P value of weighted rvPDT test statistic (test.v1).}    
\item{pvalue.v0}{ P value of unweighted rvPDT test statistic (test.v0).}  
}

 
\references{    

Guo W , Shugart YY, Does Haplotype-based Collapsing Tests Gain More Power than Variant-based Collapsing Tests for Detecting Rare Variants in Pedigrees  
(manuscript). 
 
} 

 
 