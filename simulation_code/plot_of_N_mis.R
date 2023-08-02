library(ggplot2)
library("gridExtra")
library("cowplot")
library(ggpubr)
library(latex2exp)

load('/Users/zewei/Desktop/HKUgrad/G-LDSC/plots/plot.data.v2.Rdata')
p.hole=24/0.05
simu.name='DHS_newN'
simu.time=50
remake.ldsc.res=F
if (remake.ldsc.res){
  for (N in c(500,2500,5000,10000,20000,50000)){
    for (prob.causal in c(0.05)) {
      for (heritability in c(0.04,0.4)) {
        for (DHSenrichment in c(1,2,3)) {
          gwas.out<-paste0('/Users/zewei/local/str-HDL/sim_gwas/',
                           simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N',N)
          ldsc.out<-paste0('/Users/zewei/local/str-HDL/result/',
                           simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N',N)
          #ldsc pre
          ldsc.ll<-c()
          ldsc.enrichment<-0
          ldsc.enrichment.emp<-NA
          ldsc.taus<-NA
          #simulation start
          for (i in 1:simu.time) {
            #ldsc
            ldsc.temp<-read.table(paste0(ldsc.out,'/ldsc_AFR_',i,
                                         '.results'),header = T,stringsAsFactors = F)
            #ldsc.category<-ldsc.temp[,1]
            ldsc.temp[,'Enrich_test']<-ifelse(ldsc.temp$Enrichment_p<p.hole,1,0)
            #ldsc.temp<-as.matrix(ldsc.temp)
            ldsc.enrichment<-ldsc.enrichment+ldsc.temp[,-1]
            #enrichment column, may have some problem
            ldsc.enrichment.emp<-rbind(ldsc.enrichment.emp,as.vector(ldsc.temp[,3]/ldsc.temp[1,3]/M.anno))
            #ldsc.enrichment.emp<-rbind(ldsc.enrichment.emp,as.vector(ldsc.temp[,5]))
            tau.ldsc<-ldsc.temp$Coefficient
            ldsc.taus<-rbind(ldsc.taus,tau.ldsc)
          }
          #ldsc res pre
          ldsc.enrichment<-ldsc.enrichment/simu.time
          ldsc.enrichment.emp<-ldsc.enrichment.emp[-1,]
          ldsc.taus.emp<-ldsc.taus[-1,]
          
          
          
          ###analysis ldsc result
          ldsc.h2<-read.table(paste0(ldsc.out,'/ldsc_h2_AFR.txt'),header = F,sep = '(')
          
          for (i in 1:nrow(ldsc.h2)) {
            ldsc.h2[i,'var']<-substr(ldsc.h2[i,2],1,4)
            ldsc.h2[i,'h2']<-strsplit(ldsc.h2[i,1], ': ')[[1]][[2]]
          }
          ldsc.h2[,'var']<-as.numeric(ldsc.h2[,'var'])
          ldsc.h2[,'h2']<-as.numeric(ldsc.h2[,'h2'])
          
          ldsc.h2<-ldsc.h2[,c('h2','var')]
          
          raw.result<-list()
          raw.result[['ldsc.res']]=ldsc.enrichment.emp
          raw.result[['ldsc.h2']]=ldsc.h2
          raw.result[['ldsc.taus']]=ldsc.taus.emp
          saveRDS(raw.result,paste0(gwas.out,'/ldsc.result_AFR.Rdata'))
        }
      }
    }
  }
}


#grab data mis
for (N in c(500,2500,5000,10000,20000,50000)){
  plot.result.grabber<-NA
  for (prob.causal in c(0.05)) {
    for (heritability in c(0.04,0.4)) {
      for (DHSenrichment in c(1,2,3)) {
        gwas.out<-paste0('/Users/zewei/local/str-HDL/sim_gwas/',
                         simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N',N)
        result<-readRDS(paste0(gwas.out,'/gls.result_UKB.Rdata'))
        ldsc.result<-readRDS(paste0(gwas.out,'/ldsc.result_mis.Rdata'))
        gls.res<-result
        ldsc.enrichment.emp<-ldsc.result$ldsc.res
        ldsc.taus.emp<-ldsc.result$ldsc.taus
        
        M<-ncol(ldsc.enrichment.emp)
        #statistics of gls
        total.stat.n<-7
        
        gls.res.summary<-gls.res.summary.first<-list()
        gls.ll<-c()
        true.ll<-c()
        
        
        gls.taus.res<-NA
        gls.taus.var<-NA
        gls.h2.res<-NA
        gls.e.res<-NA
        gls.pvalue<-NA
        gls.h2.var.res<-NA
        gls.stat<-NA
        gls.intercept<-NA
        for (i in 1:simu.time) {
          gls.taus.res<-cbind(gls.taus.res, gls.res[[total.stat.n*(i-1)+1]])
          #gls.taus.var<-cbind(gls.taus.var, gls.res[[total.stat.n*(i-1)+2]])
          gls.h2.res<-cbind(gls.h2.res, gls.res[[total.stat.n*(i-1)+2]])
          gls.e.res<-cbind(gls.e.res, gls.res[[total.stat.n*(i-1)+3]])
          gls.pvalue<-cbind(gls.pvalue, gls.res[[total.stat.n*(i-1)+7]])
          #gls.h2.var.res<-cbind(gls.h2.var.res, gls.res[[total.stat.n*(i-1)+6]])
          gls.stat<-cbind(gls.stat, gls.res[[total.stat.n*(i-1)+6]])
          gls.intercept<-cbind(gls.intercept,gls.res[[total.stat.n*(i-1)+4]])
        }
        gls.taus.res<-gls.taus.res[,-1]
        #gls.taus.var<-gls.taus.var[,-1]
        gls.h2.res<-gls.h2.res[,-1]
        gls.e.res<-gls.e.res[,-1]
        gls.pvalue<-gls.pvalue[,-1]
        #gls.h2.var.res<-gls.h2.var.res[,-1]
        gls.stat<-gls.stat[,-1]
        gls.intercept<-gls.intercept[,-1]
        
        
        
        gls.h2.res.mean<-apply(gls.h2.res,1,mean)
        gls.e.res.mean<-apply(gls.e.res,1,mean)
        gls.stat.mean<-apply(gls.stat,1,mean)
        gls.pvalue.mean<-apply(gls.pvalue,1,mean)
        #gls.enrich.test.mean<-apply(gls.enrich.test,1,mean)
        #gls.h2.var.res.mean<-apply(gls.h2.var.res,1,mean)
        #gls.taus.res.mean<-apply(gls.taus.res, 1, mean)
        #gls.taus.sd.mean<-apply(sqrt(gls.taus.var), 1, mean)
        #gls.taus.sd.emp<-apply(gls.taus.res, 1, sd)
        
        
        #box plot of enrichment
        gls.box<-ldsc.box<-sample.box<-h2.box<-gls.tau.box<-ldsc.tau.box<-c()
        for (i in 1:M) {
          gls.box<-c(gls.box,as.vector(gls.e.res[i,]))
          ldsc.box<-c(ldsc.box,as.vector(ldsc.enrichment.emp[,i]))
          gls.tau.box<-c(gls.tau.box,as.vector(gls.taus.res[i,]))
          ldsc.tau.box<-c(ldsc.tau.box,as.vector(ldsc.taus.emp[,i]))
        }
        #e.box<-c(gls.box,ldsc.box,sample.box)
        e.box<-c(gls.box,ldsc.box)
        tau.box<-c(gls.tau.box,ldsc.tau.box)
        #group.box<-c(rep('gls',M*simu),rep('ldsc',M*simu),rep('sample',M*simu))
        #cate.box<-rep(1:M,times=3,each=simu)
        group.box<-c(rep('gls',M*simu.time),rep('ldsc',M*simu.time))
        cate.box<-rep(1:M,times=2,each=simu.time)
        
        
        box.data<-as.data.frame(cbind(e.box,tau.box,group.box,cate.box))
        box.data$e.box<-as.numeric(box.data$e.box)
        box.data$tau.box<-as.numeric(box.data$tau.box)
        number_of_cate=which(rownames(gls.e.res)=='DHS_Trynka')
        
        ##take info for plot 1 out#####
        ldsc.h2<-as.vector(ldsc.result$ldsc.h2[,1])
        gls.h2<-as.vector(gls.h2.res[1,])
        est.h2.temp<-c(gls.h2,ldsc.h2)
        
        plot.result.grabber.temp=box.data[box.data$cate.box==number_of_cate,]
        plot.result.grabber.temp[,'causal']=rep(prob.causal,2*simu.time)
        plot.result.grabber.temp[,'h2']=rep(heritability,2*simu.time)
        plot.result.grabber.temp[,'Enrichment']=rep(DHSenrichment,2*simu.time)
        plot.result.grabber.temp[,'Est.h2']=est.h2.temp
        #add tau
        plot.result.grabber.temp[,'Est.tau1']=box.data[box.data$cate.box==1,'tau.box']
        plot.result.grabber<-rbind(plot.result.grabber,plot.result.grabber.temp)
      }
    }
  }
  plot.result.grabber<-na.omit(plot.result.grabber)
  assign(paste0('mis_plot.N',N),plot.result.grabber)
}

#plot4 comparison in h2 with sample size
#plot.4.data<-rbind(plot.1.data,plot.1.lessN2.data,plot.1.lessN1.data)
mis_plot.4.data<-rbind(mis_plot.N50000,mis_plot.N20000,mis_plot.N10000,
                   mis_plot.N5000,mis_plot.N2500,mis_plot.N500)
mis_plot.4.data$group.box[mis_plot.4.data$group.box=='gls']='g-LDSC'
mis_plot.4.data$group.box[mis_plot.4.data$group.box=='ldsc']='s-LDSC'
mis_plot.4.data[,'N']<-rep(c('n = 50,000','n = 20,000','n = 10,000',
                         'n = 5,000','n = 2,500','n = 500'),times=1,each=600)
mis_plot.4.data[,'N']<-factor(mis_plot.4.data[,'N'],levels = c('n = 500','n = 2,500','n = 5,000',
                                                       'n = 10,000','n = 20,000','n = 50,000'))
#control of show how many N
mis_plot.4.data2<-mis_plot.4.data[which(mis_plot.4.data$N %in% c('n = 2,500','n = 5,000',
                                                     'n = 10,000','n = 50,000')),]
mis_plot.4.data1<-mis_plot.4.data[which(mis_plot.4.data$N %in% c('n = 500','n = 2,500','n = 5,000',
                                                     'n = 10,000')),]

mis_plot.4.data1<-mis_plot.4.data
mis_plot.4.data2<-mis_plot.4.data

mis_plot.data<-mis_plot.4.data2[which(mis_plot.4.data1$h2=='0.4'&
                                mis_plot.4.data1$causal=='0.05'),]




#grab data AFR
for (N in c(500,2500,5000,10000,20000,50000)){
  plot.result.grabber<-NA
  for (prob.causal in c(0.05)) {
    for (heritability in c(0.04,0.4)) {
      for (DHSenrichment in c(1,2,3)) {
        gwas.out<-paste0('/Users/zewei/local/str-HDL/sim_gwas/',
                         simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N',N)
        result<-readRDS(paste0(gwas.out,'/gls.result_UKBAFR.Rdata'))
        ldsc.result<-readRDS(paste0(gwas.out,'/ldsc.result_AFR.Rdata'))
        gls.res<-result
        ldsc.enrichment.emp<-ldsc.result$ldsc.res
        ldsc.taus.emp<-ldsc.result$ldsc.taus
        
        M<-ncol(ldsc.enrichment.emp)
        #statistics of gls
        total.stat.n<-7
        
        gls.res.summary<-gls.res.summary.first<-list()
        gls.ll<-c()
        true.ll<-c()
        
        
        gls.taus.res<-NA
        gls.taus.var<-NA
        gls.h2.res<-NA
        gls.e.res<-NA
        gls.pvalue<-NA
        gls.h2.var.res<-NA
        gls.stat<-NA
        gls.intercept<-NA
        for (i in 1:simu.time) {
          gls.taus.res<-cbind(gls.taus.res, gls.res[[total.stat.n*(i-1)+1]])
          #gls.taus.var<-cbind(gls.taus.var, gls.res[[total.stat.n*(i-1)+2]])
          gls.h2.res<-cbind(gls.h2.res, gls.res[[total.stat.n*(i-1)+2]])
          gls.e.res<-cbind(gls.e.res, gls.res[[total.stat.n*(i-1)+3]])
          gls.pvalue<-cbind(gls.pvalue, gls.res[[total.stat.n*(i-1)+7]])
          #gls.h2.var.res<-cbind(gls.h2.var.res, gls.res[[total.stat.n*(i-1)+6]])
          gls.stat<-cbind(gls.stat, gls.res[[total.stat.n*(i-1)+6]])
          gls.intercept<-cbind(gls.intercept,gls.res[[total.stat.n*(i-1)+4]])
        }
        gls.taus.res<-gls.taus.res[,-1]
        #gls.taus.var<-gls.taus.var[,-1]
        gls.h2.res<-gls.h2.res[,-1]
        gls.e.res<-gls.e.res[,-1]
        gls.pvalue<-gls.pvalue[,-1]
        #gls.h2.var.res<-gls.h2.var.res[,-1]
        gls.stat<-gls.stat[,-1]
        gls.intercept<-gls.intercept[,-1]
        
        
        
        gls.h2.res.mean<-apply(gls.h2.res,1,mean)
        gls.e.res.mean<-apply(gls.e.res,1,mean)
        gls.stat.mean<-apply(gls.stat,1,mean)
        gls.pvalue.mean<-apply(gls.pvalue,1,mean)
        #gls.enrich.test.mean<-apply(gls.enrich.test,1,mean)
        #gls.h2.var.res.mean<-apply(gls.h2.var.res,1,mean)
        #gls.taus.res.mean<-apply(gls.taus.res, 1, mean)
        #gls.taus.sd.mean<-apply(sqrt(gls.taus.var), 1, mean)
        #gls.taus.sd.emp<-apply(gls.taus.res, 1, sd)
        
        
        #box plot of enrichment
        gls.box<-ldsc.box<-sample.box<-h2.box<-gls.tau.box<-ldsc.tau.box<-c()
        for (i in 1:M) {
          gls.box<-c(gls.box,as.vector(gls.e.res[i,]))
          ldsc.box<-c(ldsc.box,as.vector(ldsc.enrichment.emp[,i]))
          gls.tau.box<-c(gls.tau.box,as.vector(gls.taus.res[i,]))
          ldsc.tau.box<-c(ldsc.tau.box,as.vector(ldsc.taus.emp[,i]))
        }
        #e.box<-c(gls.box,ldsc.box,sample.box)
        e.box<-c(gls.box,ldsc.box)
        tau.box<-c(gls.tau.box,ldsc.tau.box)
        #group.box<-c(rep('gls',M*simu),rep('ldsc',M*simu),rep('sample',M*simu))
        #cate.box<-rep(1:M,times=3,each=simu)
        group.box<-c(rep('gls',M*simu.time),rep('ldsc',M*simu.time))
        cate.box<-rep(1:M,times=2,each=simu.time)
        
        
        box.data<-as.data.frame(cbind(e.box,tau.box,group.box,cate.box))
        box.data$e.box<-as.numeric(box.data$e.box)
        box.data$tau.box<-as.numeric(box.data$tau.box)
        number_of_cate=which(rownames(gls.e.res)=='DHS_Trynka')
        
        ##take info for plot 1 out#####
        ldsc.h2<-as.vector(ldsc.result$ldsc.h2[,1])
        gls.h2<-as.vector(gls.h2.res[1,])
        est.h2.temp<-c(gls.h2,ldsc.h2)
        
        plot.result.grabber.temp=box.data[box.data$cate.box==number_of_cate,]
        plot.result.grabber.temp[,'causal']=rep(prob.causal,2*simu.time)
        plot.result.grabber.temp[,'h2']=rep(heritability,2*simu.time)
        plot.result.grabber.temp[,'Enrichment']=rep(DHSenrichment,2*simu.time)
        plot.result.grabber.temp[,'Est.h2']=est.h2.temp
        #add tau
        plot.result.grabber.temp[,'Est.tau1']=box.data[box.data$cate.box==1,'tau.box']
        plot.result.grabber<-rbind(plot.result.grabber,plot.result.grabber.temp)
      }
    }
  }
  plot.result.grabber<-na.omit(plot.result.grabber)
  assign(paste0('AFR_plot.N',N),plot.result.grabber)
}

#plot4 comparison in h2 with sample size
#plot.4.data<-rbind(plot.1.data,plot.1.lessN2.data,plot.1.lessN1.data)
AFR_plot.4.data<-rbind(AFR_plot.N50000,AFR_plot.N20000,AFR_plot.N10000,
                       AFR_plot.N5000,AFR_plot.N2500,AFR_plot.N500)
AFR_plot.4.data$group.box[AFR_plot.4.data$group.box=='gls']='g-LDSC'
AFR_plot.4.data$group.box[AFR_plot.4.data$group.box=='ldsc']='s-LDSC'
AFR_plot.4.data[,'N']<-rep(c('n = 50,000','n = 20,000','n = 10,000',
                             'n = 5,000','n = 2,500','n = 500'),times=1,each=600)
AFR_plot.4.data[,'N']<-factor(AFR_plot.4.data[,'N'],levels = c('n = 500','n = 2,500','n = 5,000',
                                                               'n = 10,000','n = 20,000','n = 50,000'))
#control of show how many N
AFR_plot.4.data2<-AFR_plot.4.data[which(AFR_plot.4.data$N %in% c('n = 2,500','n = 5,000',
                                                                 'n = 10,000','n = 50,000')),]
AFR_plot.4.data1<-AFR_plot.4.data[which(AFR_plot.4.data$N %in% c('n = 500','n = 2,500','n = 5,000',
                                                                 'n = 10,000')),]

AFR_plot.4.data2<-AFR_plot.4.data
AFR_plot.4.data1<-AFR_plot.4.data


AFR_plot.data<-AFR_plot.4.data2[which(AFR_plot.4.data1$h2=='0.4'&
                                        AFR_plot.4.data1$causal=='0.05'),]


#grab data
for (N in c(500,2500,5000,10000,20000,50000)){
  plot.result.grabber<-NA
  for (prob.causal in c(0.05)) {
    for (heritability in c(0.04,0.4)) {
      for (DHSenrichment in c(1,2,3)) {
        gwas.out<-paste0('/Users/zewei/local/str-HDL/sim_gwas/',
                         simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N',N)
        result<-readRDS(paste0(gwas.out,'/gls.result.Rdata'))
        ldsc.result<-readRDS(paste0(gwas.out,'/ldsc.result_v2.Rdata'))
        gls.res<-result
        ldsc.enrichment.emp<-ldsc.result$ldsc.res
        ldsc.taus.emp<-ldsc.result$ldsc.taus
        
        M<-ncol(ldsc.enrichment.emp)
        #statistics of gls
        total.stat.n<-7
        
        gls.res.summary<-gls.res.summary.first<-list()
        gls.ll<-c()
        true.ll<-c()
        
        
        gls.taus.res<-NA
        gls.taus.var<-NA
        gls.h2.res<-NA
        gls.e.res<-NA
        gls.pvalue<-NA
        gls.h2.var.res<-NA
        gls.stat<-NA
        gls.intercept<-NA
        for (i in 1:simu.time) {
          gls.taus.res<-cbind(gls.taus.res, gls.res[[total.stat.n*(i-1)+1]])
          #gls.taus.var<-cbind(gls.taus.var, gls.res[[total.stat.n*(i-1)+2]])
          gls.h2.res<-cbind(gls.h2.res, gls.res[[total.stat.n*(i-1)+2]])
          gls.e.res<-cbind(gls.e.res, gls.res[[total.stat.n*(i-1)+3]])
          gls.pvalue<-cbind(gls.pvalue, gls.res[[total.stat.n*(i-1)+7]])
          #gls.h2.var.res<-cbind(gls.h2.var.res, gls.res[[total.stat.n*(i-1)+6]])
          gls.stat<-cbind(gls.stat, gls.res[[total.stat.n*(i-1)+6]])
          gls.intercept<-cbind(gls.intercept,gls.res[[total.stat.n*(i-1)+4]])
        }
        gls.taus.res<-gls.taus.res[,-1]
        #gls.taus.var<-gls.taus.var[,-1]
        gls.h2.res<-gls.h2.res[,-1]
        gls.e.res<-gls.e.res[,-1]
        gls.pvalue<-gls.pvalue[,-1]
        #gls.h2.var.res<-gls.h2.var.res[,-1]
        gls.stat<-gls.stat[,-1]
        gls.intercept<-gls.intercept[,-1]
        
        
        
        gls.h2.res.mean<-apply(gls.h2.res,1,mean)
        gls.e.res.mean<-apply(gls.e.res,1,mean)
        gls.stat.mean<-apply(gls.stat,1,mean)
        gls.pvalue.mean<-apply(gls.pvalue,1,mean)
        #gls.enrich.test.mean<-apply(gls.enrich.test,1,mean)
        #gls.h2.var.res.mean<-apply(gls.h2.var.res,1,mean)
        #gls.taus.res.mean<-apply(gls.taus.res, 1, mean)
        #gls.taus.sd.mean<-apply(sqrt(gls.taus.var), 1, mean)
        #gls.taus.sd.emp<-apply(gls.taus.res, 1, sd)
        
        
        #box plot of enrichment
        gls.box<-ldsc.box<-sample.box<-h2.box<-gls.tau.box<-ldsc.tau.box<-c()
        for (i in 1:M) {
          gls.box<-c(gls.box,as.vector(gls.e.res[i,]))
          ldsc.box<-c(ldsc.box,as.vector(ldsc.enrichment.emp[,i]))
          gls.tau.box<-c(gls.tau.box,as.vector(gls.taus.res[i,]))
          ldsc.tau.box<-c(ldsc.tau.box,as.vector(ldsc.taus.emp[,i]))
        }
        #e.box<-c(gls.box,ldsc.box,sample.box)
        e.box<-c(gls.box,ldsc.box)
        tau.box<-c(gls.tau.box,ldsc.tau.box)
        #group.box<-c(rep('gls',M*simu),rep('ldsc',M*simu),rep('sample',M*simu))
        #cate.box<-rep(1:M,times=3,each=simu)
        group.box<-c(rep('gls',M*simu.time),rep('ldsc',M*simu.time))
        cate.box<-rep(1:M,times=2,each=simu.time)
        
        
        box.data<-as.data.frame(cbind(e.box,tau.box,group.box,cate.box))
        box.data$e.box<-as.numeric(box.data$e.box)
        box.data$tau.box<-as.numeric(box.data$tau.box)
        number_of_cate=which(rownames(gls.e.res)=='DHS_Trynka')
        
        ##take info for plot 1 out#####
        ldsc.h2<-as.vector(ldsc.result$ldsc.h2[,1])
        gls.h2<-as.vector(gls.h2.res[1,])
        est.h2.temp<-c(gls.h2,ldsc.h2)
        
        plot.result.grabber.temp=box.data[box.data$cate.box==number_of_cate,]
        plot.result.grabber.temp[,'causal']=rep(prob.causal,2*simu.time)
        plot.result.grabber.temp[,'h2']=rep(heritability,2*simu.time)
        plot.result.grabber.temp[,'Enrichment']=rep(DHSenrichment,2*simu.time)
        plot.result.grabber.temp[,'Est.h2']=est.h2.temp
        #add tau
        plot.result.grabber.temp[,'Est.tau1']=box.data[box.data$cate.box==1,'tau.box']
        plot.result.grabber<-rbind(plot.result.grabber,plot.result.grabber.temp)
      }
    }
  }
  plot.result.grabber<-na.omit(plot.result.grabber)
  assign(paste0('plot.N',N),plot.result.grabber)
}

#plot4 comparison in h2 with sample size
#plot.4.data<-rbind(plot.1.data,plot.1.lessN2.data,plot.1.lessN1.data)
plot.4.data<-rbind(plot.N50000,plot.N20000,plot.N10000,
                   plot.N5000,plot.N2500,plot.N500)
plot.4.data$group.box[plot.4.data$group.box=='gls']='g-LDSC'
plot.4.data$group.box[plot.4.data$group.box=='ldsc']='s-LDSC'
plot.4.data[,'N']<-rep(c('n = 50,000','n = 20,000','n = 10,000',
                         'n = 5,000','n = 2,500','n = 500'),times=1,each=600)
plot.4.data[,'N']<-factor(plot.4.data[,'N'],levels = c('n = 500','n = 2,500','n = 5,000',
                                                       'n = 10,000','n = 20,000','n = 50,000'))
#control of show how many N
plot.4.data2<-plot.4.data[which(plot.4.data$N %in% c('n = 2,500','n = 5,000',
                                                     'n = 10,000','n = 50,000')),]
plot.4.data1<-plot.4.data[which(plot.4.data$N %in% c('n = 500','n = 2,500','n = 5,000',
                                                     'n = 10,000')),]
plot.4.data2<-plot.4.data
plot.4.data1<-plot.4.data


plot.data<-plot.4.data2[which(plot.4.data1$h2=='0.4'&
                                plot.4.data1$causal=='0.05'),]


#make a boxplot here
plot.data[,'Panel']='1000G_EUR'
mis_plot.data[,'Panel']='UKBB_EUR'
AFR_plot.data[,'Panel']='UKBB_AFR'

boxplot.data=rbind(plot.data,mis_plot.data,AFR_plot.data)

boxplot.data[,'Class']=paste(boxplot.data$group.box,boxplot.data$Panel,sep = '_')
boxplot.data$Enrichment=as.factor(boxplot.data$Enrichment)
boxplot.data$Class=factor(boxplot.data$Class,levels=c('g-LDSC_1000G_EUR','g-LDSC_UKBB_EUR','g-LDSC_UKBB_AFR',
                                                      's-LDSC_1000G_EUR','s-LDSC_UKBB_EUR','s-LDSC_UKBB_AFR'))





boxplot.mis1=ggplot(data=boxplot.data,aes(x=Enrichment,y=e.box, fill=Class))+
  geom_boxplot()+facet_wrap(. ~ `N`,nrow = 2)+
  stat_summary(fun.y = mean, color = "darkred", position = position_dodge(0.75),
               geom = "point", shape = c(18), size = 3,
               show.legend = FALSE)+geom_abline(slope=0, intercept=c(1,2,3),linetype = "dashed")+
  xlab('True enrichment')+ylab('Estimated enrichment')+
  coord_cartesian(ylim=c(-2, 10))+
  theme_light()+
  scale_fill_manual(values = c('#F8766D','orange','purple','#619CFF','green','yellow'))+
  #theme(axis.title.x=element_blank(),
  #      axis.text.x=element_blank(),
  #      axis.ticks.x=element_blank())+
  theme(legend.position='top')+theme(legend.title=element_blank())

boxplot.mis1



#used plot 
boxplot.data.used<-boxplot.data[boxplot.data$N=='n = 50,000',]
boxplot.mis.used=ggplot(data=boxplot.data.used,aes(x=Enrichment,y=e.box, fill=Class))+
  geom_boxplot()+
  stat_summary(fun.y = mean, color = "darkred", position = position_dodge(0.75),
               geom = "point", shape = c(18), size = 3,
               show.legend = FALSE)+geom_abline(slope=0, intercept=c(1,2,3),linetype = "dashed")+
  xlab('True enrichment')+ylab('Estimated enrichment')+
  coord_cartesian(ylim=c(0, 5))+
  theme_light()+
  scale_fill_manual(values = c('#F8766D','orange','purple','#619CFF','green','yellow'))+
  #theme(axis.title.x=element_blank(),
  #      axis.text.x=element_blank(),
  #      axis.ticks.x=element_blank())+
  theme(legend.position='top')+theme(legend.title=element_blank())

boxplot.mis.used

#save
location='/Users/zewei/Desktop/HKUgrad/G-LDSC/plots/figure/'
ggsave(paste0(location,"figure3.pdf"),plot = boxplot.mis.used, width = 174, height = 185, units = "mm")
ggsave(paste0(location,"figure3.svg"),plot = boxplot.mis.used, width = 174, height = 185, units = "mm")


#
mis_plot.data2<-mis_plot.4.data2[which(mis_plot.4.data1$h2=='0.04'&
                                        mis_plot.4.data1$causal=='0.05'),]
AFR_plot.data2<-AFR_plot.4.data2[which(AFR_plot.4.data1$h2=='0.04'&
                                         AFR_plot.4.data1$causal=='0.05'),]
plot.data2<-plot.4.data2[which(plot.4.data1$h2=='0.04'&
                                plot.4.data1$causal=='0.05'),]

plot.data2[,'Panel']='1000G_EUR'
mis_plot.data2[,'Panel']='UKBB_EUR'
AFR_plot.data2[,'Panel']='UKBB_AFR'


boxplot.data2=rbind(plot.data2,mis_plot.data2,AFR_plot.data2)

boxplot.data2[,'Class']=paste(boxplot.data2$group.box,boxplot.data2$Panel,sep = '_')
boxplot.data2$Enrichment=as.factor(boxplot.data2$Enrichment)

boxplot.mis2=ggplot(data=boxplot.data2,aes(x=Enrichment,y=e.box, fill=Class))+
  geom_boxplot()+facet_wrap(. ~ `N`,nrow = 2)+
  stat_summary(fun.y = mean, color = "darkred", position = position_dodge(0.75),
               geom = "point", shape = c(18), size = 3,
               show.legend = FALSE)+geom_abline(slope=0, intercept=c(1,2,3),linetype = "dashed")+
  xlab('True enrichment')+ylab('Estimated enrichment')+
  coord_cartesian(ylim=c(-6, 20))+
  theme_light()+
  scale_fill_manual(values =  c('#F8766D','orange','purple','#619CFF','green','yellow'))+
  #theme(axis.title.x=element_blank(),
  #      axis.text.x=element_blank(), 
  #      axis.ticks.x=element_blank())+
  theme(legend.position='top')+theme(legend.title=element_blank())

boxplot.mis2
