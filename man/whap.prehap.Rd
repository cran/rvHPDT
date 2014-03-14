\name{whap.prehap}
\alias{whap.prehap}
\title{Prepare haplotype pairs for hPDT tests in pedigree data.}
\description{
Before running hPDT test, haplotype pairs are inferred by calling MERLIN in linux for all pedigree members, and then perpare some internal statistics. 
Require the R package of "gregmisc" and MERLIN software.
}  
\usage{ 
whap.prehap(ped,map, merlinDir="", outFN.prefix="merlin",aff=2,trace=FALSE)  
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
\item{merlinDir}{indicates the directory of Merlin, for example, merlinDir="./Merlin/"; use the default="" when Merlin is in current directory or your bin directory.}  
\item{outFN.prefix}{Requests that output file of MERLIN names should be derived from outFN.prefix. For example, when it is set to be "merlin" as default, estimated haplotypes should be stored in a file called merlin.chr.}  
\item{aff}{indicates the values that represents affected status in ped data; default is 2.} 
\item{trace}{indicates whether or not the intermediate outcomes should be printed; default is FALSE.}  
} 
 
       
\value{
\item{SNPname}{ SNP names of testing.}   
\item{hapData}{ Haplotype data for each individuals.}   
\item{freq}{ Estimated frequencies of haplotypes.}    
\item{trans}{ Transmission matrix of haplotypes.}              
\item{hapScore}{ Score matrix of haplotypes.}    
}

 
\references{    

Guo W , Shugart YY, Does Haplotype-based Collapsing Tests Gain More Power than Variant-based Collapsing Tests for Detecting Rare Variants in Pedigrees  
(manuscript). 
 
} 

 
 