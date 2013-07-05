 
##################################################################################################### 
#   1.  hapPDT.test
##################################################################################################### 

 

'hapPDT.test'<-function(preHapList,aff=2,unaff=1, trace=FALSE){#   score<-sqrt(Nhap*freq*(1-freq)) 

  mydaoshu<-function(x){y<-NULL; for (i in 1:length(x)) {
              if (is.na(x[i])) y[i]<-NA else {if (x[i]==0) y[i]<-NA else y[i]<-1/x[i]}}
              return(y)
     }  
    
   hapData<- convert.factors.to.strings.in.dataframe(preHapList$hapData)
   trans<-convert.factors.to.strings.in.dataframe(preHapList$trans) 
   freq<-convert.factors.to.strings.in.dataframe(preHapList$freq)  
   hapScoreMatrix<-preHapList$hapScore  
   hapScore<-hapScoreMatrix[,"score"]  
   names(hapScore)<-hapScoreMatrix[,"hapName"] 

   haplist<-unique(as.vector(trans[,c("Trans","Untrans")]))   #remove the most common one when sum up
   del.hap<-which(haplist==names(which.max(freq))) 
   for (h in 1:length(haplist)) if ("?" %in% unlist(strsplit(haplist[h],split=""))) del.hap<-c(del.hap,h)
   if (length(del.hap)>=1) haplist<-haplist[-del.hap]

   weight<-mydaoshu(as.numeric(hapScore[haplist])) 



   famlist<-unique(as.vector(trans[,"FID"]))
   nfam<-length(famlist)
   nhap<-length(haplist)
   
   # (1) make work tdt matrix for each family
   tdtM<-array(NA,dim=c(nfam,nhap))
   colnames(tdtM)<- haplist  
   for (f in 1:nfam){# f=1;h=1
            setS<-which(trans[,"FID"]==famlist[f] & trans[,"PHENO"]==aff) 
            for (h in 1:nhap) {
                setS.trans<-which(trans[,"FID"]==famlist[f] & trans[,"PHENO"]==aff & trans[setS,"Trans"]==haplist[h]) 
                setS.untrans<-which(trans[,"FID"]==famlist[f] & trans[,"PHENO"]==aff & trans[setS,"Untrans"]==haplist[h])  
                tdtM[f,h]<-sum(as.numeric(trans[setS.trans,"weight"]),na.rm=TRUE)- sum(as.numeric(trans[setS.untrans,"weight"]),na.rm=TRUE)
            }
   }   
   # here please make it to multiple siblings later for S test in DSP.


   # (2) make work sib pairs matrix (one is aff, and the other is unaff) for each family 
   dspM<-array(NA,dim=c(nfam,nhap))
   colnames(dspM)<- haplist   
   for (f in 1:nfam){# f=1;h=1;p=1
            setS<-which(hapData[,"FID"]==famlist[f]  ) 
            parents<-unique(hapData[setS,c("FA","MO")])  
            parents<-parents[which(!( is.na(parents[,"FA"]) | is.na(parents[,"MO"]) | parents[,"FA"]==0 | parents[,"MO"]==0) ),]
  
            if (is.null(dim(parents))) npar<-1 else npar<-nrow(parents)  
            if (npar>=1){
                dspM.f<-array(0,dim=c(npar,nhap)) 
                for (p in 1:npar){
                   fid<-parents[p,"FA"]
                   mid<-parents[p,"MO"]
                  
                   for (h in 1:nhap) { 
                        setS.aff<-which(hapData[,"FA"]==fid & hapData[,"MO"]==mid & hapData[,"FID"]==famlist[f]  & hapData[,"PHENO"]==aff )
                        setS.unaff<-which(hapData[,"FA"]==fid & hapData[,"MO"]==mid & hapData[,"FID"]==famlist[f]  & hapData[,"PHENO"]==unaff ) 
                        setS.aff.h<-which(hapData[,"FA"]==fid & hapData[,"MO"]==mid & hapData[,"FID"]==famlist[f]  & hapData[,"PHENO"]==aff & hapData[,"hapName"]==haplist[h] )
                        setS.unaff.h<-which(hapData[,"FA"]==fid & hapData[,"MO"]==mid & hapData[,"FID"]==famlist[f]  & hapData[,"PHENO"]==unaff & hapData[,"hapName"]==haplist[h] ) 
                        
                        naff<-sum(as.numeric(hapData[setS.aff,"weight"]),na.rm=TRUE)/2
                        nunaff<-sum(as.numeric(hapData[setS.unaff,"weight"]),na.rm=TRUE)/2
                        if (naff>=1 & nunaff>=1)  
                          dspM.f[p,h]<-nunaff*sum(as.numeric(hapData[setS.aff.h,"weight"]),na.rm=TRUE)-
                                       naff*sum(as.numeric(hapData[setS.unaff.h,"weight"]),na.rm=TRUE)
                         # print(c(naff,nunaff,dspM.f[p,h],setS.aff,setS.unaff))
                    } #h 
               }# npar
               dspM[f,]<-colSums(dspM.f,na.rm=TRUE)
            }  else dspM[f,]<-rep(0,nhap)  #if npar
    }#f

   # (3) sum up two matrix, which have same famID order and same dimension.
    pdtM<-tdtM+dspM
    rownames(pdtM)<-famlist 
    weight[is.na(weight)]<-0
    pdt.test<-sum(pdtM %*% weight ,na.rm=TRUE)/sqrt(sum( (pdtM %*% weight)^2,na.rm=TRUE))
    pdt.pvalue<- 2* pnorm(q= abs(pdt.test), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)  # 2*pnorm(q= 1.96, lower.tail = FALSE )=0.05

    pdt.test.v0<-sum(pdtM,na.rm=TRUE)/sqrt(sum( (pdtM)^2,na.rm=TRUE))
    pdt.pvalue.v0<- 2* pnorm(q= abs(pdt.test.v0), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)  # 2*pnorm(q= 1.96, lower.tail = FALSE )=0.05


 
    names(weight)<-haplist
    ans<-list(TDT=tdtM,Sib=dspM,PDT=pdtM,W=weight,test.v1=pdt.test,pvalue.v1=pdt.pvalue,test.v0=pdt.test.v0,pvalue.v0=pdt.pvalue.v0) 
    return(ans)      
       
}

 

 
##################################################################################################### 
#   2.  rvPDT.test
##################################################################################################### 

 
'rvPDT.test'<-function(ped, aff=2,unaff=1, trace=FALSE){  

 
     mydaoshu<-function(x){y<-NULL; for (i in 1:length(x)) {
              if (is.na(x[i])) y[i]<-NA else {if (x[i]==0) y[i]<-NA else y[i]<-1/x[i]}}
              return(y)
     }  

   snpset=1:(ncol(ped)-6)
   snplist<-colnames(ped)[6+snpset]  
   famlist<-unique(as.vector(ped[,"FID"]))
   nfam<-length(famlist)
   nsnp<-length(snplist)
   ped<- convert.factors.to.strings.in.dataframe(ped[,c(1:6,6+snpset)])  
   all.freq<-NULL
   Nmiss<-NULL
   for (i in 1:length(snplist)){
       this.all<-ped[which((is.na(ped[,"FA"]) & is.na(ped[,"MO"])) |(ped[,"FA"]==0 & ped[,"MO"]==0)), 6+i]  
      # this.all<-ped[, 6+i]   
      Nmiss[i]<-length(!is.na(this.all))  
      all.freq[i]<-(sum(as.numeric(this.all),na.rm=TRUE))/(2*Nmiss[i])  # @@@@@@@ 3/12  @@@@@@@ 3/26

    }  
   names(all.freq)<-snplist 
   weight<-mydaoshu(as.numeric(sqrt(Nmiss*all.freq*(1-all.freq))) )   
   names(weight)<-snplist
 
   # (1) make work tdt, dsp  pairs matrix (one is aff, and the other is unaff) for each family 
   tdtM<-array(NA,dim=c(nfam,nsnp))
   colnames(tdtM)<- snplist  
   dspM<-array(NA,dim=c(nfam,nsnp))
   colnames(dspM)<- snplist   
   for (f in 1:nfam){# f=1;h=1;p=1
            setS<-which(ped[,"FID"]==famlist[f]  ) 
            parents<-unique(ped[setS,c("FA","MO")])
            parents<-parents[which(!( is.na(parents[,"FA"]) | is.na(parents[,"MO"]) | parents[,"FA"]==0 | parents[,"MO"]==0) ),]
            if (is.null(dim(parents))) {npar<-1; parents<-t(as.matrix(parents))} else npar<-nrow(parents)   
 
            if (npar>=1){
                tdtM.f<-array(0,dim=c(npar,nsnp)) 
                dspM.f<-array(0,dim=c(npar,nsnp)) 

                for (p in 1:npar){
                   fid<-parents[p,"FA"]
                   mid<-parents[p,"MO"]
                   work.data<-ped[which(ped[,"FA"]==fid & ped[,"MO"]==mid & ped[,"FID"]==famlist[f]),] 
                   naff<-length(which(as.numeric(work.data[,"PHENO"])==aff))
                   nunaff<-length(which(as.numeric(work.data[,"PHENO"])==unaff))
                   if (naff>=1 & nunaff>=1){  
                      for (h in 1:nsnp)  
                      dspM.f[p,h]<-nunaff*sum(as.numeric(ped[which(as.numeric(work.data[,"PHENO"])==aff),6+h]),na.rm=TRUE)-
                                     naff*sum(as.numeric(ped[which(as.numeric(work.data[,"PHENO"])==unaff),6+h]),na.rm=TRUE)
                   }
                   if (naff>=1){  
                      for (h in 1:nsnp)  
                      tdtM.f[p,h]<-2*sum(as.numeric(work.data[which(as.numeric(work.data[,"PHENO"])==aff),6+h]),na.rm=TRUE)-
                                     naff*sum(as.numeric(ped[which(as.character(ped[,"IID"])==fid &  ped[,"FID"]==famlist[f]),6+h]), 
                                              as.numeric(ped[which(as.character(ped[,"IID"])==mid &  ped[,"FID"]==famlist[f]),6+h]),na.rm=TRUE )  #NA leads error here
                   }  
               } #p th pair
               dspM[f,]<-colSums(dspM.f,na.rm=TRUE)
               tdtM[f,]<-colSums(tdtM.f,na.rm=TRUE)

            }  else {dspM[f,]<-rep(0,nsnp);tdtM[f,]<-rep(0,nsnp)}

    }

   # (2) sum up two matrix, which have same famID order and same dimension.
    pdtM<-tdtM+dspM
    rownames(pdtM)<-famlist 
    weight[is.na(weight)]<-0 
    pdt.test<-sum(pdtM %*% weight ,na.rm=TRUE)/sqrt(sum( (pdtM %*% weight)^2,na.rm=TRUE))
    pdt.pvalue<- 2* pnorm(q= abs(pdt.test), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)  # 2*pnorm(q= 1.96, lower.tail = FALSE )=0.05
 
    pdt.test.v0<-sum(pdtM,na.rm=TRUE)/sqrt(sum( (pdtM)^2,na.rm=TRUE))
    pdt.pvalue.v0<- 2* pnorm(q= abs(pdt.test.v0), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)  # 2*pnorm(q= 1.96, lower.tail = FALSE )=0.05 
   
    ans<-list(TDT=tdtM,Sib=dspM,PDT=pdtM,W=weight,test.v1=pdt.test,pvalue.v1=pdt.pvalue,test.v0=pdt.test.v0,pvalue.v0=pdt.pvalue.v0)  
    return(ans)      
 
}

 
  



 
##################################################################################################### 
#   3.  whap.prehap subroutine for hPDT.test
##################################################################################################### 

   
'whap.prehap'<-function(ped,map,outFN.prefix="merlin", aff=2,trace=FALSE){ 
 

   #  try(system("rm merlin*") )   
   #  install.packages("gregmisc") 

   library(gregmisc)
   snpset=1:(ncol(ped)-6) 
   weight.alpha=1.28 
   snp.colname<-colnames(ped)[6+snpset] 
   if (trace) print(paste("whap.prehap-0: start to work on ",getwd(),sep=""))
 

   #(1) format data for merlin
   ped<- convert.factors.to.strings.in.dataframe(ped[,c(1:6,6+snpset)])
   merlin.ped<-array(NA,dim=c(nrow(ped),5+(ncol(ped)-6)*2))
   for (i in 1:4) merlin.ped[,i]<-as.character(ped[,i]) 
   merlin.ped[,5]<-as.numeric(ped[,5]) 
   for (j in 1:(ncol(ped)-6)){
      for (i in 1:nrow(ped)){
         x<-as.numeric(ped[i,j+6])
         if (is.na(x)) geno2c<-c(NA,NA)  else {# 0 is missing alleles in merlin.ped
            if (x==0) geno2c<-c(2,2)
            if (x==1) geno2c<-c(1,2)
            if (x==2) geno2c<-c(1,1)  # 1 is the rare allele
            }
         merlin.ped[i,5+2*j-1]<-geno2c[1]
         merlin.ped[i,5+2*j]<-geno2c[2]
      }
    }
  

	  merlin.map<-map[which(map$MARKER %in% colnames(ped)),]
       merlin.dat<-cbind("M",colnames(ped)[-(1:6)]) 	 
       
       pedFN<-paste(outFN.prefix,".ped",sep="")
       datFN<-paste(outFN.prefix,".dat",sep="")
       mapFN<-paste(outFN.prefix,".map",sep="") 
       chrFN<-paste(outFN.prefix,".chr",sep="") 
     
     if (trace) print(paste("whap.prehap-1: Will write merlin files to",pedFN,sep=""))
 

    
       write.table(merlin.ped,file=pedFN,col.names=FALSE,row.names=FALSE,quote=FALSE)
       write.table(merlin.dat,file=datFN,col.names=FALSE,row.names=FALSE,quote=FALSE) 
       write.table(merlin.map,file=mapFN,col.names=TRUE,row.names=FALSE,quote=FALSE)

       #(2) call merlin 
     
       system(paste("merlin -d",datFN,"-p",pedFN,"-m",mapFN,"-x NA --best --horizontal --prefix ",outFN.prefix,sep=" ")) 
       if (trace) print("whap.prehap-2: get merlin data for haplotyping")

   

     #(3.1) read merlin haplotype --best output 
    
     myremove<-function(u){ t<-unlist(strsplit(u,NULL));  y<-t[which(!t %in% LETTERS)];return(as.numeric(y))}
     mystrsplit<-function(x){
        y1<-NULL;y2<-NULL
        for (i in 1:length(x)) {
          t<-unlist(strsplit(x[i],"\\,"))
          if (length(t)==1) {y1[i]=t; y2[i]=t} else {y1[i]=myremove(t[1]);y2[i]=myremove(t[2])}
        }
        ans<-list();ans[[1]]<-y1;ans[[2]]<-y2
        return(ans)
     } 
    mypaste<-function(x) paste(x,collapse="")  
  
      

     if (file.exists(chrFN)){ #make fam.hap= rewrite hap matrix according weights.
       if (trace) print(paste("whap.prehap-3.1: Read merlin output:",chrFN,sep="")) 

       hap<-read.table(file=chrFN,header=FALSE,colClasses = "character",stringsAsFactors=FALSE,fill=TRUE)   
       hap<- convert.factors.to.strings.in.dataframe(hap)


       fam.start.line<-which(hap[,1]=="FAMILY") 
       fam.hap<-NULL
       for (i in 1:length(fam.start.line)){#i=1;j=1;k=1
         last.line<-ifelse(i==length(fam.start.line),nrow(hap), fam.start.line[i+1]-1)
         FID<-hap[fam.start.line[i],2]
         fam.hap.v1<- cbind(FID,hap[(fam.start.line[i]+1):last.line,])        
         colnames(fam.hap.v1)<-c("FID","IID","hapSource",colnames(ped)[6+snpset])
  
         amb.loci<-rep(NA,ncol(fam.hap.v1)) # #### &&&### loci with A1,2 format in haplotype
         for (k in 4:ncol(fam.hap.v1)) amb.loci[k]<-max(unlist(lapply(fam.hap.v1[,k],nchar)))  
         which.amb<-which(amb.loci>=3)

  
         if (length(which.amb)>=1){ 
             amb.permu.set<- permutations(n=2,r=length(which.amb),repeats.allowed=TRUE)   #combinations for ambigous loci
          
             for (j in 1:nrow(amb.permu.set)){ 
                fam.hap.v2<-fam.hap.v1  
                for (k in 1:length(which.amb)){
                   fam.hap.v2[,which.amb[k]]<-mystrsplit(fam.hap.v1[,which.amb[k]])[[amb.permu.set[j,k]]]                  
                }#kth amb loci
             fam.hap.v2$weight<-1/length(amb.permu.set)
             fam.hap.v2$subFID<-j   
             fam.hap<-rbind(fam.hap,  fam.hap.v2) 
            }#j for permu
         } else     fam.hap<-rbind(fam.hap,  cbind(fam.hap.v1,weight=1,subFID=1)) 
 
       } #i family
   

     if (trace) print("whap.prehap-3.2:: get initial hap data with equal weights")


     
       fam.hap.paste=apply(fam.hap[,4:(ncol(fam.hap)-2)],1,mypaste)   
       all.hap<-unique(fam.hap.paste) 
       snp.col2<- 3:ncol(hap)
       hap.paste<-NULL; for (i in 1:nrow(hap)) hap.paste[i]<- paste(hap[i,snp.col2],collapse="") 
       hap.freq<-NULL
       for (i in 1:length(all.hap)){ 
          hap.freq[i]<-length(which(hap.paste==all.hap[i]))  +1    #use all family members to estimate hap freq.     
       }     
       hap.freq<- hap.freq/sum(hap.freq)                                          
       names(hap.freq)<-all.hap

    
    hap.freq.corr<-hap.freq 
    for (i in 1:length(all.hap))  if (hap.freq[i]==0) hap.freq.corr[i]<-0.00001
    hap.freq.corr<-hap.freq.corr/sum(hap.freq.corr)  
    if (trace) print("whap.prehap-4:: get hap.freq from founders")

   
    amb.family.ID<-unique(fam.hap[which(fam.hap$subFID>=2),"FID"])
    if (  length(amb.family.ID)>=1){ 
     
       this.snp.col<-4:(ncol(fam.hap)-2) 
       for (f in 1:length(amb.family.ID)){# f=1;a=1;b=1
         this.fam<-fam.hap[which(fam.hap$hapSource=="(FOUNDER)"  & fam.hap$FID==amb.family.ID[f]),]
         this.subFID<-unique(this.fam$subFID)
         this.prob<-NULL 
         for (b in 1:nrow(this.fam))  this.fam$hap[b]<- paste(this.fam[b,this.snp.col],collapse="")
         for (a in 1:length(this.subFID)){
             
              this.ind<-unique(this.fam[which(this.fam$subFID==this.subFID[a]),"IID"])
              this.prob[a]<-1 #calculate prob for individuals for sub.subFID a
              for (c in 1:length(this.ind)){
                this.ind.c.hap<- this.fam[which(this.fam$subFID==this.subFID[a]  & this.fam$IID==this.ind[c]  ),"hap"]
                if (length(this.ind.c.hap)!=2) print("warning: no 2 haplotypes for an individual!")
                het.coeff<-ifelse(this.ind.c.hap[1]==this.ind.c.hap[2],1,2)
                this.prob[a]<-this.prob[a]*het.coeff*prod(hap.freq.corr[this.ind.c.hap],na.rm=TRUE) 
                
              }
         }
         this.prob<-this.prob/sum(this.prob) 
         for (a in 1:length(this.subFID)){
              fam.hap[which(fam.hap$FID==amb.family.ID[f]  &  fam.hap$subFID==this.subFID[a]),"weight"]<-this.prob[a] 
          }#a
        }#f   
      }#if  
      fam.hap<- convert.factors.to.strings.in.dataframe(fam.hap) 
      fam.hap<-merge(fam.hap,ped[,c(1:6)],by=c("FID","IID"),all.x=TRUE,all.y=FALSE)
      snp.col3<- which( colnames(fam.hap) %in% snp.colname )  
      #for (i in 1:nrow(fam.hap)) fam.hap$hapName[i]<- paste(fam.hap[i,snp.col3],collapse="")   
      fam.hap$hapName=apply(fam.hap[,snp.col3],1,mypaste)

      if (trace) print(paste("whap.prehap-5:: get hap.freq with freq.weights on ",length(amb.family.ID)," ambigous families.",sep="") )

  

    ######################################################################
    #   Section 2: based on weighted fam.hap, make transmission matrix### 
    ######################################################################

    # (4) make transmission matrix 
    snp.col4<- which( colnames(fam.hap) %in% snp.colname )  
    child.list<-which(fam.hap$hapSource!="(FOUNDER)")
    trans<-NULL
    for (i in 1:length(child.list)){# i=1
       lc<-as.vector(fam.hap[child.list[i],])
       famID<-as.character(lc["FID"])
       if (as.character(lc["hapSource"])=="(MATERNAL)") parentID<-lc["MO"] else parentID<-lc["FA"] 
       parentID<-as.character(unlist(parentID))
 

       pair<-paste(parentID,lc["IID"],sep="->")
       subFID.i<- unique(fam.hap[which(fam.hap[,"FID"]==famID & fam.hap[,"IID"]==parentID),"subFID"]) 
       for (j in 1:length(subFID.i)){
         parent.line<- which(fam.hap[,"FID"]==famID & fam.hap[,"IID"]==parentID & fam.hap[,"subFID"]==subFID.i[j] )
         parent.weight<-unique(fam.hap[parent.line,"weight"])
         if (length(parent.weight)!=1) {print(fam.hap[ which(fam.hap[,"FID"]==famID),]);print(parent.line);
                                        print(paste("famID,subFID,child=",famID,subFID.i[j],child.list[i],sep=",")); stop("STOP : wrong weights for subFamily")}
        
         Trans.hap<-as.character(lc["hapName"])
         Untrans.hap<- ifelse(fam.hap[parent.line,"hapName"][1]==Trans.hap,  fam.hap[parent.line,"hapName"][2],fam.hap[parent.line,"hapName"][1])

 
         ans<-as.character(unlist(c(famID,pair,Trans.hap,Untrans.hap,parent.weight,lc["PHENO"],lc["subFID"]))) 
         trans<-rbind(trans,ans)  
       }
     if (i %% 500==0  & trace) print(paste("In tdt matrix, finish ", i, " child",collpase=""))
    }
    colnames(trans)<-c("FID","pair","Trans","Untrans","weight","PHENO","subFID")
    rownames(trans)<-NULL
    if (trace) print("whap.prehap-6: Get trans matrix")

   # (5)  get haplotype score for future weights.

   hap.list<-unique(as.vector(trans[,c("Trans","Untrans")]))
   hap.score<-NULL
   for (i in 1:length(hap.list)){
     tar.hap<-hap.list[i]
     T<-sum(as.numeric(trans[which(trans[,"Trans"]==tar.hap  &  as.numeric(trans[,"PHENO"])==aff),"weight"]),na.rm=TRUE)
     U<-sum(as.numeric(trans[which(trans[,"Untrans"]==tar.hap &  as.numeric(trans[,"PHENO"])==aff),"weight"]),na.rm=TRUE)
 
     tdt<-(T-U)/sqrt(T+U) 
     ans<-c(tar.hap,T,U,tdt,T+U) 
     hap.score<-rbind(hap.score,ans)
   }
   colnames(hap.score)<-c("hapName","T","U","tdt","N")
   Nhap<-sum(as.numeric(hap.score[,"N"]),na.rm=TRUE)
   freq<-hap.freq[as.character(hap.score[,"hapName"])]  
  
   score<- sqrt(Nhap*freq*(1-freq))  
   score.alpha<-score
   set1<-which(as.numeric(hap.score[,"tdt"]) >=weight.alpha) 
   set2<-which(as.numeric(hap.score[,"tdt"]) <=-1*weight.alpha)
   score.alpha[set2]<- -1*score.alpha[set2]
   score.alpha[-union(set1,set2)]<- 0 


   hap.score<-cbind(hap.score,freq,score,score.alpha) 


   miss.hap<-NULL  #set score for ??????? haplotype
   charAlikeB<-function(A,B,missing="?"){
           if (nchar(A)==nchar(B)){ ans<-TRUE
                   Alist<-unlist(strsplit(A,split=NULL))
                   Blist<-unlist(strsplit(B,split=NULL))
                   for (s in 1:nchar(A)) if (!(Alist[s]==Blist[s] | Alist[s]==missing | Blist[s]==missing)) ans<-FALSE
               } else    ans<-FALSE
     return(ans)
   }  
   miss.hap<-NULL; for (i in 1:length(hap.list))  if ("?" %in% unlist(strsplit(hap.list[i],split=NULL)) ) miss.hap<-c(miss.hap,hap.list[i])
   if (length(miss.hap)>=1){
      for (i in 1:length(miss.hap)){
         tar.hap<-miss.hap[i]
         tar.hap.inc<-NULL
         for (j in 1:length(hap.list))
            if ( !(hap.list[j] %in% miss.hap) &  charAlikeB(tar.hap, hap.list[j]) )  tar.hap.inc<-c(tar.hap.inc,hap.list[j]) 
         tar.hap.inc.i<-which(hap.score[,"hapName"] %in% tar.hap.inc)
         Escore<-sum(freq[tar.hap.inc.i]*score[tar.hap.inc.i],na.rm=TRUE)
         Escore.alpha<-sum(freq[tar.hap.inc.i]*score.alpha[tar.hap.inc.i],na.rm=TRUE)
         hap.score[which(hap.score[,"hapName"]==tar.hap),"score"]<-Escore
         hap.score[which(hap.score[,"hapName"]==tar.hap),"score.alpha"]<-Escore.alpha  
     } 
    }     
  
     
 
 } else {fam.hap<-NA   
         hap.freq<-NA  
         trans<-NA
         hap.score<-NA 
         print(getwd()) 
         print(paste("Did not find merlin.chr files for haplotyping output from merlin file:",chrFN,sep=""))
        }#when no .chr file found

  
   if (trace) print("whap.prehap-7: Get haplotype score based on freq and tdt") 
   if (trace) print("whap.prehap-8: Done prehap!") 
   return(list(SNPname=snp.colname,hapData=fam.hap,freq=hap.freq, trans=trans,hapScore=hap.score))


}




   'convert.factors.to.strings.in.dataframe' <- function(dataframe){
        class.data  <- sapply(dataframe, class)
        factor.vars <- class.data[class.data == "factor"]
        for (colname in names(factor.vars))
        { dataframe[,colname] <- as.character(dataframe[,colname])}
        return (dataframe)
       }

##################################################################################################### 
#   4: main-- call 2 tests
#####################################################################################################  
 
'rhapPDT'<-function(ped, map,  aff=2, unaff=1, merlinFN.prefix="merlin", trace=TRUE){  

 
       count.snp<-NULL; for (k in 7:ncol(ped)) count.snp[k]<-sum(as.numeric(ped[,k]),na.rm=TRUE)
       bad.snp<-which(count.snp==0)
       if (length(bad.snp)>=1) ped<-ped[,-bad.snp]
       print(paste("delete ", length(bad.snp), " snps", sep=""))
       Nsnp<-ncol(ped)-6 
      
 

       prehap<-try(whap.prehap(ped=ped,map=map,aff=aff,outFN.prefix=merlinFN.prefix,trace=trace))  
       if (trace) try(print(paste(c("Begin test on hapdata with dimension=",dim(prehap$hapData)),collapse="  ")))  
       t1<-try(hapPDT.test(preHapList=prehap, aff=aff,unaff=unaff, trace=trace))  
       t2<-try(rvPDT.test(ped=ped, aff=aff,unaff=unaff, trace=trace))   
       
       if (is.list(t1)) hap.p<-c(t1$pvalue.v1,t1$pvalue.v0) else hap.p<-c(NA,NA)
       if (is.list(t2)) rvpdt.p<-c(t2$pvalue.v1,t2$pvalue.v0) else rvpdt.p<-c(NA,NA)
       ans<-list(hPDT_v0=hap.p[2],hPDT_v1=hap.p[1],rvPDT_v0=rvpdt.p[2],rvPDT_v1=rvpdt.p[1])
       if (trace) try(print(ans))
       
      
       return(ans)
}  

 

 
 
############################################################################  
#        end
############################################################################ 
   

