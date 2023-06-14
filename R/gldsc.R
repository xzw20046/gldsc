#' Title
#'
#' @param panel Path of LDSM.Rdata file.
#' @param gwas Path of GWAS summary input file.
#' @param jackknife Whether use Jackknife to estimate standard error.
#' @param intercept Whether consider confounding bias in the analysis.
#' @param numCores Number of cores used in analysis.
#'
#' @return
#' @export
#'
#' @examples
gldsc<-function(panel,gwas,out,jackknife,intercept,numCores=4){
  start_time=Sys.time()
  message('Analysis begin')
  ref.pannel<-readRDS(panel)
  message('load reference panel...')
  snplist.pan<-unlist(ref.pannel[['snp.list']])
  message(paste0('read ',length(snplist.pan),' SNPs from reference panel'))
  registerDoParallel(numCores)
  x.sum<-Reduce(sum,lapply(ref.pannel[['x.gls']], sum))
  
  
  #gls
  gwas.df<-read.table(gwas,header = T,fill = TRUE)
  message('read GWAS...')
  gwas.df<-na.omit(gwas.df)
  message(paste0('read ',nrow(gwas.df),' SNPs from GWAS'))
  snp.used=intersect(snplist.pan,gwas.df$SNP)
  message(paste0('After merging with reference panel ',length(snp.used),' SNPs remain'))
  #
  merged.loc<-match(snp.used,snplist.pan)
  
  gwas.df<-gwas.df[match(snp.used, gwas.df$SNP),]
  y.sum<-sum(gwas.df[,'Z']**2-1,na.rm = T)
  raw.Ntau=y.sum/x.sum
  
  gwasN<-median(gwas.df$N,na.rm = T)
  gls.left<-gls.right<-0
  raw.result.all<-list()
  for (chr in 1:length(ref.pannel[['chrs']])) {
    chr.exact<-substr(ref.pannel[['chrs']][chr],nchar('ldblk_1kg_')+1,nchar(ref.pannel[['chrs']][chr])-nchar('.hdf5'))
    block.names.all<-names(ref.pannel[['LDSM']])[grep(x = names(ref.pannel[['LDSM']]), pattern = paste0(chr.exact,'_'))]
    total.block<-length(block.names.all)
    
    if(total.block>0){
      raw.result<-foreach (i=1:total.block) %dopar% gls.left.right(i=i,x=ref.pannel[['x.gls']],y=gwas.df,
                                                                   block.names.all = block.names.all,
                                                                   LDmatrix = ref.pannel[['LDmatrix']],LDSM = ref.pannel[['LDSM']],
                                                                   snp.list = ref.pannel[['snp.list']],
                                                                   raw.Ntau = raw.Ntau,intercept = intercept
      )
      for (batch in 1:length(raw.result)) {
        gls.left<-gls.left+raw.result[[batch]][[1]]
        gls.right<-gls.right+raw.result[[batch]][[2]]
      }
      raw.result.all[[chr.exact]]<-raw.result
    }
    message(paste0(chr.exact,' calculation complited'))
  }
  
  
  #first time result
  result.one<-gls.estimator(left = gls.left,right = gls.right,
                            A=ref.pannel[['Atotal']],M.anno = ref.pannel[['M.anno']],anno.total = ref.pannel[['anno.total']],N = gwasN,intercept = intercept)
  message('point estimate completed')
  #result
  if(jackknife==F){
    message(Sys.time()-start_time)
    return(result.one)
  }else if(jackknife==T){
    message('start Jackknife estimation')
    #jackknife starts here
    jack.tau<-jack.h2<-jack.e<-jack.intersecp<-jack.stat<-NA
    jack.time=0
    for (chr in 1:length(ref.pannel[['chrs']])) {
      chr.exact<-substr(ref.pannel[['chrs']][chr],nchar('ldblk_1kg_')+1,nchar(ref.pannel[['chrs']][chr])-nchar('.hdf5'))
      raw.result.temp=raw.result.all[[chr.exact]]
      total.block2<-length(raw.result.temp)
      jack.result<-foreach (i=1:total.block2) %dopar% oneout(i,raw.result=raw.result.temp,left.total = gls.left,right.total = gls.right,
                                                             A=ref.pannel[['Atotal']],M.anno = ref.pannel[['M.anno']],
                                                             anno.total = ref.pannel[['anno.total']],N = gwasN,intercept = intercept)
      for (i in 1:total.block2) {
        jack.tau<-rbind(jack.tau,as.vector(jack.result[[i]][[1]]))
        jack.h2<-rbind(jack.h2,as.vector(jack.result[[i]][[2]]))
        jack.e<-rbind(jack.e,as.vector(jack.result[[i]][[3]]))
        #jack.intersecp<-c(jack.intersecp,jack.result[[i]][[4]])
        #jack.stat<-rbind(jack.stat,as.vector(jack.result[[i]][[5]]))
        jack.time=jack.time+1
      }
    }
    jack.tau.var<-apply(jack.tau[-1,], 2, var)*(jack.time-1)^2/jack.time
    jack.h2.var<-apply(jack.h2[-1,], 2, var)*(jack.time-1)^2/jack.time
    jack.e.var<-apply(jack.e[-1,], 2, var)*(jack.time-1)^2/jack.time
    #jack.intersecp.var<-var(jack.intersecp[-1])*(jack.time-1)^2/jack.time
    #jack.stat.var<-apply(jack.stat[-1,], 2, var)*(jack.time-1)^2/jack.time
    
    
    result.one[['tau.var.jack']]=jack.tau.var
    result.one[['h2.var.jack']]=jack.h2.var
    result.one[['e.var.jack']]=jack.e.var
    #result.one[['inter.var.jack']]=jack.intersecp.var
    #result.one[['P']]=pnorm(result.one[['e.stat']]/sqrt(jack.stat.var),lower.tail=F)
    message(Sys.time()-start_time)
    return(result.one)
  }else{print('Jackknife should be TRUE of FALSE')}
}