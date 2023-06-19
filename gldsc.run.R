#input parameters
args <- commandArgs(trailingOnly = TRUE)

LD.path <- gsub(x = args[grep(x = args, pattern = "LDpath=")],
                     pattern = "LDpath=", replacement = "")
anno.path <- gsub(x = args[grep(x = args, pattern = "annopath=")],
                     pattern = "annopath=", replacement = "")
maf.path <- gsub(x = args[grep(x = args, pattern = "mafpath=")],
                     pattern = "mafpath=", replacement = "")
annotation <- gsub(x = args[grep(x = args, pattern = "annotation=")],
                     pattern = "annotation=", replacement = "")
snp <- gsub(x = args[grep(x = args, pattern = "snplist=")],
                     pattern = "snplist=", replacement = "")
out <- gsub(x = args[grep(x = args, pattern = "out=")],
                     pattern = "out=", replacement = "")
MAF <- gsub(x = args[grep(x = args, pattern = "maf=")],
                     pattern = "maf=", replacement = "")
numCores <- gsub(x = args[grep(x = args, pattern = "cores=")],
                     pattern = "cores=", replacement = "")
panel <- gsub(x = args[grep(x = args, pattern = "panel=")],
                     pattern = "panel=", replacement = "")
gwas <- gsub(x = args[grep(x = args, pattern = "gwas=")],
                     pattern = "gwas=", replacement = "")
jackknife <- gsub(x = args[grep(x = args, pattern = "jackknife=")],
                     pattern = "jackknife=", replacement = "")
intercept <- gsub(x = args[grep(x = args, pattern = "intercept=")],
                     pattern = "intercept=", replacement = "")
function.path <- gsub(x = args[grep(x = args, pattern = "function=")],
                     pattern = "function=", replacement = "")
cell.path <- gsub(x = args[grep(x = args, pattern = "CELLanno=")],
                     pattern = "CELLanno=", replacement = "")                     
#default setting
#snp <- ifelse(length(snp)==0, NULL, snp)
#bug1
snp <- NULL
MAF <- ifelse(length(MAF)==0, 0.05, as.numeric(MAF))
numCores <- ifelse(length(numCores)==0, 4, as.numeric(numCores))
intercept <- ifelse(length(intercept)==0, T, as.logical(intercept))
annotation <- ifelse(length(annotation)==0, 'all',
                     as.vector(unlist(strsplit(annotation,split = ','))))
jackknife <- ifelse(length(jackknife)==0, T, as.logical(jackknife))

source(function.path)
if(length(LD.path)>0 & length(panel)==0){
  message('Begin panel caluclation')
  estimated.panel<-mkLDSM(LD.path,anno.path,maf.path,annotation,snp,out,MAF,numCores)
  saveRDS(estimated.panel,paste0(out,'/LDSM.pannel.Rdata'))
}else if(length(panel)>0 & length(LD.path)==0){
  if (length(cell.path)==0){
  message('Begin gLDSC analysis')
  result<-gldsc(panel,gwas,out,jackknife,intercept,numCores)
  saveRDS(result,paste0(out,'.result.Rdata'))    
  }else {
    message('Begin gLDSC cell type analysis')
    result<-gLDSC_cts(panel,gwas,out,cell.path,jackknife,intercept,numCores)
    saveRDS(result,paste0(out,'.result.Rdata'))
  }

}else {message('No paremeter input')}



