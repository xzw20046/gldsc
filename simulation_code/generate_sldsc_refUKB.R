ref.pannel<-readRDS('/Users/zewei/local/str-HDL/data/chr10UKBAFR.pannel.Rdata')
attach(ref.pannel)
#make sldsc refpannel
outpath='/Users/zewei/local/ldsc_new_ref'
pre.anno<-read.delim('/Users/zewei/local/ldsc_dhs_ref/baseline.10.annot',header = T)

#
used.snp=unlist(snp.list)
anno2=pre.anno[match(used.snp,pre.anno$SNP),]

#ldscore
ldsc.file=NULL
for (i in 1:length(x.gls)) {
  ldsc.file=rbind(ldsc.file,x.gls[[i]])
}
names(ldsc.file)=names(anno2)[-1:-4]
ldsc.file1=anno2[,1:4]
ldsc.file1[,'MAF']=0.1
ldsc.file1=ldsc.file1[,c('CHR','SNP','BP','CM','MAF')]

ldsc.file2=cbind(ldsc.file1,ldsc.file)

w.file=apply(anno2[,-1:-4],2,sum)

#output
write.table(anno2,paste0(outpath,'/baselineUKBAFR.10.annot'),quote = F,row.names = F,sep = '\t')
write.table(ldsc.file2,paste0(outpath,'/baselineUKBAFR.10.l2.ldscore'),quote = F,row.names = F,sep = '\t')
write.table(t(w.file),paste0(outpath,'/baselineUKBAFR.10.l2.M_5_50'),quote = F,row.names = F,col.names = F,sep = '\t')


#mk model misspecifile ref for ldsc
ref.pannel<-load('/Users/zewei/local/str-HDL/data/chr10ref.EAout.pannel.Rdata')
#attach(ref.pannel)
#make sldsc refpannel
outpath='/Users/zewei/local/ldsc_new_ref'
pre.anno<-read.delim('/Users/zewei/local/ldsc_dhs_ref/baseline.10.annot',header = T)
#match(colnames(Atotal),colnames(pre.anno))

pre.anno<-pre.anno[,c(1:4,match(colnames(Atotal),colnames(pre.anno)))]
#
used.snp=unlist(snp.list)
anno2=pre.anno[match(used.snp,pre.anno$SNP),]

#ldscore
ldsc.file=NULL
for (i in 1:length(x.gls)) {
  ldsc.file=rbind(ldsc.file,x.gls[[i]])
}
names(ldsc.file)=names(anno2)[-1:-4]
ldsc.file1=anno2[,1:4]
ldsc.file1[,'MAF']=0.1
ldsc.file1=ldsc.file1[,c('CHR','SNP','BP','CM','MAF')]

ldsc.file2=cbind(ldsc.file1,ldsc.file)

w.file=apply(anno2[,-1:-4],2,sum)

#output
write.table(anno2,paste0(outpath,'/baselineEAout.10.annot'),quote = F,row.names = F,sep = '\t')
write.table(ldsc.file2,paste0(outpath,'/baselineEAout.10.l2.ldscore'),quote = F,row.names = F,sep = '\t')
write.table(t(w.file),paste0(outpath,'/baselineEAout.10.l2.M_5_50'),quote = F,row.names = F,col.names = F,sep = '\t')




