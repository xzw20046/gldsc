library(xlsx)

real.result<-readRDS('/Users/zewei/Desktop/HKUgrad/G-LDSC/plots/summary.result2.Rdata')
pheno<-read.csv('/Users/zewei/Desktop/HKUgrad/G-LDSC/ukbb.sumstat.csv',header = T)
i=1
gLDSC_h2<-gLDSC_intercept<-c()
for (i in 1:15) {
  ldsc_ukb=real.result[[1]][[i]]
  gls_ukb=real.result[[2]][[i]]
  
  Annotation=ldsc_ukb$Category
  sLDSC_Enrichment=ldsc_ukb$Enrichment
  sLDSC_Enrichment_SD=ldsc_ukb$Enrichment_std_error
  sLDSC_P=ldsc_ukb$Enrichment_p
  
  gLDSC_Enrichment=gls_ukb[[3]]
  gLDSC_Enrichment_SD=sqrt(gls_ukb[[10]])
  gLDSC_P=gls_ukb[[7]]
  
  gLDSC_intercept<-c(gLDSC_intercept,gls_ukb[[4]])
  gLDSC_h2<-c(gLDSC_h2,gls_ukb[[2]][1])
  
  write_file=as.data.frame(cbind(Annotation,gLDSC_Enrichment,gLDSC_Enrichment_SD,gLDSC_P,
                                 sLDSC_Enrichment,sLDSC_Enrichment_SD,sLDSC_P))
  names(write_file)=c('Annotation','gLDSC Enrich','gLDSC Enrich_SD','gLDSC P',
                      'sLDSC Enrich','sLDSC Enrich_SD','sLDSC P')
  #write.xlsx(write_file,"/Users/zewei/Desktop/HKUgrad/G-LDSC/plots/UKBB_detail.xlsx",
  #           sheetName=pheno$Phenotype2[i],
  #           col.names=TRUE,append=TRUE,row.names = F)
}

h2.out=ldsc_ukb=real.result[[1]][['H2']]
h2.out<-cbind(h2.out,gLDSC_h2,gLDSC_intercept)
h2.out<-as.data.frame(h2.out)
h2.out<-cbind(pheno$Phenotype2,h2.out)
names(h2.out)<-c('Phenotype','s-LDSC h2','var','s-LDSC intercept','g-LDSC h2','g-LDSC intecept')
for (i in 2:6) {
  h2.out[,i]=as.numeric(h2.out[,i])
  h2.out[,i]=round(h2.out[,i],3)
}
h2.out
h2.out2<-h2.out[,c('Phenotype','g-LDSC h2','g-LDSC intecept','s-LDSC h2','s-LDSC intercept')]
write.xlsx(h2.out2,"/Users/zewei/Desktop/HKUgrad/G-LDSC/plots/UKBB_h2_intercept.xlsx",
                      col.names=TRUE,append=TRUE,row.names = F)
