

load('/Users/zewei/local/str-HDL/data/chr10ref.pannel.Rdata')
#load('/Users/zewei/local/str-HDL/data/chr10ref.DHS.pannel.Rdata')

source('/Users/zewei/Desktop/HKUgrad/G-LDSC/R/mlfun.R')
simu<-50
p.hole<-0.05/24
simu.name='DHS_newN'

registerDoParallel(5)

x.sum<-Reduce(sum,lapply(x.gls, sum))

#change DHSenrichment to c(1,2,3)
#N in c(5000,50000,500000)
for (N in c(500,2500,5000,10000,20000,50000)) {
  for (prob.causal in c(0.05)) {
    for (heritability in c(0.04)) {
      for (DHSenrichment in c(1,2,3)) {
        gwas.out<-paste0('/Users/zewei/local/str-HDL/sim_gwas/',
                         simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N',N)
        ldsc.out<-paste0('/Users/zewei/local/str-HDL/result/',
                         simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N',N)
        start_time <- Sys.time()
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
          gwas.df<-read.table(paste0(gwas.out,'/chr10simgwas_test',simu.time,'.sumstat'),header = T)
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
        
        
        Sys.time()-start_time
        
        
        raw.result<-gls.res.summary.first
        
        
        #saveRDS(raw.result,paste0(gwas.out,'/gls.result.2cat.Rdata'))
        saveRDS(raw.result,paste0(gwas.out,'/gls.result.Rdata'))
      }
    }
  }
  
}



