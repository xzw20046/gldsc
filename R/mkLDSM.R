#' Title
#'
#' @param LD.path The LD matrix files to use in calculating LD score matrix. Under this folder all LD matrixs file should be in .hdf5 format.
#' @param anno.path The location of the annotation. The input format here remain the same as .annot in ldsc.
#' @param maf.path The location of the MAF of SNPs. The input format here remain the same as .M_5_50 in ldsc.
#' @param annotation Specify a vector of annotation name to be used in the analysis. 'all' means use all of them.
#' @param snp Specify a list of SNPs 
#' @param MAF MAF for filtering SNPs in the analysis. Defalt set as 0.05. 
#' @param numCores Number of cores used in the analysis.
#'
#' @return A Rdata file contains all LD score matrix information used in gldsc analysis.
#' @export
#'
#' @examples
mkLDSM<-function(LD.path,anno.path,maf.path,annotation='all',snp=NULL,MAF,numCores=4){
  #register cores
  registerDoParallel(numCores)
  #Check reference files#####
  LD.files <- list.files(LD.path)
  if(any(grepl(x = LD.files, pattern = "ldblk_*"))){
    chrs <- LD.files[grep(x = LD.files, pattern = "ldblk_*")]
  }else{
    error.message <- "Ref pannel (LDfiles) input error"
    stop(error.message)
  }
  
  anno.files<- list.files(anno.path)
  if(any(grepl(x = anno.files, pattern = "*.annot.gz"))){
    annos <- anno.files[grep(x = anno.files, pattern = "*.annot.gz")]
    ldscs <- anno.files[grep(x = anno.files, pattern = "*.ldscore.gz")]
  }else{
    error.message <- "Ref pannel (ANNOfile) input error"
    stop(error.message)
  }
  
  maf.files<- list.files(maf.path)
  if(any(grepl(x = maf.files, pattern = "1000G.mac5eur.*"))){
    mafs <- maf.files[grep(x = maf.files, pattern = "1000G.mac5eur.*")]
  }else{
    error.message <- "Ref pannel (MAFfiles) input error"
    stop(error.message)
  }
  
  if(length(chrs)==length(annos)&length(chrs)==length(mafs)){
    chr.length<-length(chrs)
  }else{
    error.message <- "Ref pannel (number of chr) input error"
    stop(error.message)
  }
  
  
  #Preparation#####
  if(is.null(snp)==F){
    gwas.df<-read.table(snp,header = T)
    gwas.df<-na.omit(gwas.df)
    gwas.snp<-as.vector(gwas.df[,'SNP'])
    print(paste('load',length(gwas.snp),'SNPs from custom snplist',sep = ' '))
  }
  
  #Starts#####
  snp.list<-LDmatrix<-LDSM<-x.gls<-list()
  Atotal<-NA
  start_time <- Sys.time()
  for (chr in 1:chr.length) {
    
    chr.exact<-substr(chrs[chr],nchar('ldblk_1kg_')+1,nchar(chrs[chr])-nchar('.hdf5'))
    #read ref by chr
    LDinfo<-h5ls(file = paste(LD.path,chrs[chr], sep = '/'),recursive=F)
    
    file1<-gzfile(paste(anno.path,annos[chr], sep = '/'),'rt')
    file3<-gzfile(paste(maf.path,mafs[chr], sep = '/'),'rt')
    
    anno2<-read.table(file1, header = T,sep = '')
    maf<-read.table(file3,header = T,sep = '')
    
    close(file1)
    #close(file2)
    close(file3)
    
    #select annotation
    if(annotation[1]=='all'){used.a<-names(anno2)[-4:-1]
    }else{used.a<-annotation}
    
    #select chosen annotation
    anno<-anno2[,c(names(anno2)[1:4],used.a)]
    maf.snp<-as.vector(maf[maf$FRQ<(1-MAF)&maf$FRQ>MAF,'SNP'])
    #ldsc.snp<-ldsc$SNP
    #
    
    #generate by block
    if(is.null(snp)){gwas.snp=maf.snp}
    raw.result<-foreach (batch=1:nrow(LDinfo)) %dopar% refpannel.par(batch,LDinfo=LDinfo,LD.path=LD.path,chrs=chrs[chr],
                                                                     anno=anno,maf.snp=maf.snp,gwas.snp=gwas.snp,
                                                                     used.a=used.a)
    #record result
    for (batch in 1:length(raw.result)) {
      if(raw.result[[batch]][[1]]==1){
        batch.name<-LDinfo[batch,'name']
        snp.list[[paste(chr.exact,batch.name,sep = '_')]]<-raw.result[[batch]][[2]]
        LDmatrix[[paste(chr.exact,batch.name,sep = '_')]]<-raw.result[[batch]][[3]]
        LDSM[[paste(chr.exact,batch.name,sep = '_')]]<-raw.result[[batch]][[4]]
        x.gls[[paste(chr.exact,batch.name,sep = '_')]]<-raw.result[[batch]][[5]]
        
        Anno.temp<-anno[match(raw.result[[batch]][[2]],anno[,'SNP']),used.a]
        Anno.temp<-as.matrix(Anno.temp)
        Atotal<-rbind(Atotal,Anno.temp)
      }
    }
    print(paste0('LD Score Matrix of CHR ',chr,' has been complited'))
  }
  
  Atotal<-Atotal[-1,]
  snplist<-unlist(snp.list)
  M<-length(snplist)
  anno.total<-apply(Atotal, 2, sum)
  M.anno<-anno.total/M
  print(paste(M,'SNPs left for analysis',sep = ' '))
  print(Sys.time()-start_time)
  
  #output
  result<-list()
  result[['LDinfo']]<-LDinfo
  result[['LDmatrix']]<-LDmatrix
  result[['LDSM']]<-LDSM
  result[['snp.list']]<-snp.list
  result[['x.gls']]<-x.gls
  result[['anno.total']]<-anno.total
  result[['M.anno']]<-M.anno
  result[['Atotal']]<-Atotal
  result[['chrs']]<-chrs
  return(result)
  #saveRDS(result,paste0(out,'.Rdata'))
}