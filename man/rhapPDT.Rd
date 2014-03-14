\name{rhapPDT}
\alias{rhapPDT}
\title{Calling haplotype-based and variant-based pedigree disequilibrium test for rare variants in pedigrees.}
\description{
To detecting rare variants for binary traits using general pedigrees, the pedigree disequilibrium tests are proposed by collapsing rare haplotypes/variants with/without weights.}
\usage{ rhapPDT(ped, map, aff=2, unaff=1, mu=1.04,
merlinFN.prefix="merlin", nperm=1000, trace=TRUE) 
}
    
   

 

\arguments{  
\item{ped}{input data, has same format with PLINK but having column names. The PED file is a white-space (space or tab) delimited file, and the first six columns are mandatory: 
     FID: Family ID;
     IID: Individual ID;
     FA: Paternal ID;
     MO: Maternal ID;
     SEX: Sex (1=male; 2=female; other=unknown);
     PHENO: Phenotype; 
Genotypes (column 7 onwards) should also be white-space delimited; they are coded as 0, 1 and 2, indicating the number of coding allele, and NA is for missing genotype. } 
\item{map}{input data, has same format with MAP file required by MERLIN. The MAP file is a white-space (space or tab) delimited file with 3 columns as follows,  
CHROMOSOME: chromosome (1-22, X, Y or 0 if unplaced)  
MARKER: marker name in PED file that is usually rs# or snp identifier         
POSITION: Genetic distance (morgans)
The data file and map file can include different sets of markers, but markers that are absent from the map file will be ignored by MERLIN. 
} 
\item{aff}{indicates the values that represents affected status in "PHENO" column of PED data; default is 2.}
\item{unaff}{indicates the values that represents unaffected status in "PHENO" column of PED data; default is 1.}  
\item{mu}{indicates mu value that defines causal in the training data; default is 1.04.}  
\item{merlinFN.prefix}{ Requests that output file of MERLIN names should be derived from outFN.prefix. For example, when it is set to be "merlin" as default, estimated haplotypes should be stored in a file called merlin.chr.}   
\item{nperm}{indicates the times of permutation; default is 1000.}     
\item{trace}{ Indicates whether or not the intermediate outcomes should be printed; default is FALSE.}  
}

  
       
\value{ 
\item{hPDT_v0}{ P value of unweighted haplotype PDT test statistic.}    
\item{hPDT_v1}{ P value of weighted haplotype PDT test statistic.} 
\item{rvPDT_v0}{ P value of unweighted rvPDT test statistic.}   
\item{rvPDT_v1}{ P value of weighted rvPDT test statistic.}  
}

 
\references{    

Guo W , Shugart YY, Does Haplotype-based Collapsing Tests Gain More Power than Variant-based Collapsing Tests for Detecting Rare Variants in Pedigrees  
(manuscript). 
 
} 

 
\examples{
#ped<-read.table("MLIP.ped",head=1,stringsAsFactors=FALSE)
#map<-read.table("MLIP.map",head=1,stringsAsFactors=FALSE)  
#test<-rhapPDT(ped, map, trace=TRUE) 
#test 
#$hPDT_v0
#[1] 0.4231359

#$hPDT_v1
#[1] 0.1481145

#$rvPDT_v0
#[1] 0.03237073

#$rvPDT_v1
#[1] 0.162997
}

 