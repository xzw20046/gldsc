#' function of cell type specific enricment
#'
#' @param panel Path of LDSM.Rdata file.
#' @param gwas Path of GWAS summary input file.
#' @param cts_path Path of cell type specific annotation file.
#' @param jackknife Whether use Jackknife to estimate standard error.
#' @param intercept Whether consider confounding bias in the analysis.
#' @param numCores Number of cores used in analysis.
#'
#' @return
#' @export
#'
#' @examples
gLDSC_cts<-function(panel,gwas,out,cts_path,jackknife,intercept,numCores=4){
  start_time=Sys.time()
  message('Analysis begin')
  ref.pannel<-readRDS(panel)
  message('load reference panel...')
  snplist.pan<-unlist(ref.pannel[['snp.list']])
  M=length(snplist.pan)
  message(paste0('read ',M,' SNPs from reference panel'))
  #load cts annotation matrix
  cts.files<- list.files(cts_path)
  if(any(grepl(x = cts.files, pattern = "*.annot.gz"))){
    ctss <- cts.files[grep(x = cts.files, pattern = "*.annot.gz")]
    #ldscs <- anno.files[grep(x = anno.files, pattern = "*.ldscore.gz")]
  }else{
    error.message <- "Ref pannel (Cell type ANNOfile) input error"
    stop(error.message)
  }  
  #use a function to modify LDSMatrix and LDScore
  message('modifying LD Score Matrix...')
  chr.length=length(ref.pannel[['chrs']])
  block.names.total=names(ref.pannel[['LDSM']])
  cts.total=c()
  for (chr in 1:chr.length){
    #read cts annotation by chrs
    file1<-gzfile(paste(cts_path,ctss[chr], sep = '/'),'rt')
    cts2<-read.table(file1, header = T,sep = '')
    close(file1)
    chr.exact=cts2[1,1]    
    cts.name=names(cts2)[5]
    
    #select the block names in this chr
    block.name.chr=block.names.total[grep(x=block.names.total, pattern=paste0('chr',chr.exact,'_'))]
    
    #adj ref by block
    for (block in block.name.chr){
      snp.list.temp=ref.pannel[['snp.list']][[block]]
      snp.loc.temp=match(snp.list.temp,cts2$SNP)
      if (sum(is.na(snp.loc.temp))>0){
        error.message <- "Cell type anno does not match with baseline anno "
        stop(error.message)
      }
      cts.temp=cts2[snp.loc.temp,]
      
      #message(cts.name)
      cts.vector=as.numeric(as.vector(cts.temp[,5]))
      #modefy score statrs
      LDSM=ref.pannel[['LDSM']][[block]]
      LDSC=ref.pannel[['x.gls']][[block]]
      LDm=ref.pannel[['LDmatrix']][[block]]
      middle=LDm*as.vector(sqrt(cts.vector))
      middle2=crossprod(middle)
      #message(paste0('LDSM:',dim(LDSM),' middle:', dim(middle),dim(middle2)))
      ref.pannel[['LDSM']][[block]]=LDSM+middle2
      ref.pannel[['x.gls']][[block]]=cbind(LDSC,(LDm**2)%*%cts.vector)
      cts.total=c(cts.total,cts.vector)
    }
    message(paste0('finish transforming reference matrixs of CHR ',chr.exact))
  }
  #cts.total=as.numeric(cts.total)
  ref.pannel[['Atotal']]=cbind(ref.pannel[['Atotal']],cts.total)
  colnames(ref.pannel[['Atotal']])[ncol(ref.pannel[['Atotal']])]=cts.name
  #names(ref.pannel[['Atotal']])[ncol(ref.pannel[['Atotal']])]=cts.name
  #ref.pannel[['M.Anno']][[cts.name]]=mean(cts.total)
  ref.pannel[['anno.total']][[cts.name]]=sum(cts.total)
  ref.pannel[['M.Anno']][[cts.name]]=ref.pannel[['anno.total']][[cts.name]]/M
  
  
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
    result.one[['Prob of SNP']]=ref.pannel[['M.Anno']]
    #result.one[['inter.var.jack']]=jack.intersecp.var
    #result.one[['P']]=pnorm(result.one[['e.stat']]/sqrt(jack.stat.var),lower.tail=F)
    message(Sys.time()-start_time)
    return(result.one)
  }else{print('Jackknife should be TRUE of FALSE')}
}
