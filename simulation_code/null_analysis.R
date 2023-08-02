#
devtools::install_github("xzw20046/gldsc")

run_analy=F
if(fun_analy){
  load('/Users/zewei/local/str-HDL/data/chr10ref.pannel.Rdata')
  #load('/Users/zewei/local/str-HDL/data/chr10ref.DHS.pannel.Rdata')
  ref.gwas=read.delim(paste0(gwas.out,'/chr10simgwas_test',1,'.sumstat'),header = T,sep = ' ')
  source('/Users/zewei/Desktop/HKUgrad/G-LDSC/R/mlfun.R')
  simu<-50
  p.hole<-0.05/24
  
  registerDoParallel(5)
  
  x.sum<-Reduce(sum,lapply(x.gls, sum))
  
  
  
  gwas.path<-'/Users/zewei/local/GWAS_error_50times'
  
  
  gwas.out='/Users/zewei/local/str-HDL/sim_gwas/null'
  #gwas.out<-paste0('/Users/zewei/local/str-HDL/sim_gwas/',
  #                 'DHS_newN','c',0.05,'h',0.4,'e',1,'N',20000)
  #section 1
  #gls pre
  gls.res.summary<-gls.res.summary.first<-list()
  gls.ll<-c()
  true.ll<-c()
  #ldsc pre
  ldsc.ll<-c()
  ldsc.enrichment<-0
  ldsc.enrichment.emp<-NA
  ldsc.taus<-NA
  #simulation start
  for (simu.time in 1:simu) {
    #gls
    #gwas.df<-read.table(paste0(gwas.out,'/chr10simgwas_test',simu.time,'.sumstat'),header = T)
    gwas.error<-read.table(paste0(gwas.path,'/gwas_',simu.time,'.qassoc'),header = T)
    gwas.df=ref.gwas
    gwas.df$Z=gwas.error$T
    gwas.df$N=gwas.error$NMISS
    y.sum<-sum(gwas.df[match(unlist(snp.list), gwas.df$SNP),'Z']**2-1)
    raw.Ntau=y.sum/x.sum
    
    gwasN<-median(gwas.df$N,na.rm = T)
    gls.left<-gls.right<-0
    raw.result.all<-list()
    for (chr in 2) {
      chr.exact<-substr(chrs[chr],nchar('ldblk_1kg_')+1,nchar(chrs[chr])-nchar('.hdf5'))
      
      #file2<-gzfile(paste(anno.path,ldscs[chr], sep = '/'),'rt')
      #ldsc<-read.table(file2, header = T,sep = '')
      #close(file2)
      
      
      block.names.all<-names(LDSM)[grep(x = names(LDSM), pattern = paste0(chr.exact,'_'))]
      total.block<-length(block.names.all)
      
      if(total.block>0){
        raw.result<-foreach (i=1:total.block) %dopar% gls.left.right(i=i,x=x.gls,y=gwas.df,block.names.all = block.names.all,
                                                                     LDmatrix = LDmatrix,LDSM = LDSM,snp.list = snp.list,
                                                                     raw.Ntau = raw.Ntau,intercept = T
        )
        for (batch in 1:length(raw.result)) {
          gls.left<-gls.left+raw.result[[batch]][[1]]
          gls.right<-gls.right+raw.result[[batch]][[2]]
        }
        raw.result.all[[chr.exact]]<-raw.result
      }
    }
    print(simu.time)
    
    #first time result
    result.one<-gls.estimator(left = gls.left,right = gls.right,
                              A=Atotal,M.anno = M.anno,anno.total = anno.total,N = gwasN,intercept = T)
    
    gls.res.summary.first<-c(gls.res.summary.first,result.one)
    
  }
  
  #saveRDS(gls.res.summary,paste0(gwas.out,'/gls.result.Rdata'))
  
  raw.result<-gls.res.summary.first
  
  
  #saveRDS(raw.result,paste0(gwas.out,'/gls.result.2cat.Rdata'))
  saveRDS(raw.result,paste0(gwas.out,'/gls.result.Rdata'))
  
  
  simu=50
  for (simu.time in 1:simu) {
    #gls
    #gwas.df<-read.table(paste0(gwas.out,'/chr10simgwas_test',simu.time,'.sumstat'),header = T)
    gwas.error<-read.table(paste0(gwas.path,'/gwas_',simu.time,'.qassoc'),header = T)
    gwas.df=ref.gwas
    gwas.df$Z=gwas.error$T
    gwas.df$N=gwas.error$NMISS
    write.table(gwas.df,paste0(gwas.out,'/chr10simgwas_test',simu.time,'.sumstat'),
                row.names = F,quote = F)
    
  }
}

#pull result here
remake.ldsc.res=F
if(remake.ldsc.res){
  gwas.out='/Users/zewei/local/str-HDL/sim_gwas/null'
  ldsc.out<-'/Users/zewei/local/str-HDL/result/null'
  #start_time <- Sys.time()
  #ldsc pre
  ldsc.ll<-list()
  #simulation start
  for (i in 1:simu.time) {
    #ldsc
    ldsc.temp<-read.table(paste0(ldsc.out,'/ldsc_v2_',i,
                                 '.results'),header = T,stringsAsFactors = F)
    ldsc.ll[[i]]<-ldsc.temp
  }

  saveRDS(ldsc.ll,paste0(gwas.out,'/ldsc.result_v2.Rdata'))
}
  
gwas.out='/Users/zewei/local/str-HDL/sim_gwas/null'
result<-readRDS(paste0(gwas.out,'/gls.result.Rdata'))
ldsc.result<-readRDS(paste0(gwas.out,'/ldsc.result_v2.Rdata'))

ldsc.e=NULL
gls.e=NULL
total.stat.n<-7
for (i in 1:simu) {
  ldsc.temp=ldsc.result[[i]]
  test1=as.vector(ldsc.temp[,3]/ldsc.temp[1,3]/M.anno)
  ldsc.e=rbind(ldsc.e,test1)
  gls.e<-cbind(gls.e, result[[total.stat.n*(i-1)+3]])
}
  
ldsc.mean=round(apply(ldsc.e,2,mean),2)
ldsc.sd=round(apply(ldsc.e,2,sd),2)
gls.mean=round(apply(gls.e,1,mean),2)
gls.sd=round(apply(gls.e,1,sd),2)

null.anal.res=as.data.frame(cbind(gls.mean,gls.sd,ldsc.mean,ldsc.sd))

names(null.anal.res)=c('g-LDSC Mean','(SD)','s-LDSC Mean','(SD)')
null.anal.res[,2]=paste0('(',null.anal.res[,2],')')
null.anal.res[,4]=paste0('(',null.anal.res[,4],')')

#sup table for null analysis: null.anal.res
#library(xlsx)
#write.xlsx(null.anal.res,'/Users/zewei/Desktop/HKUgrad/G-LDSC/plots/null_table.xlsx',row.names = T,col.names = T)
