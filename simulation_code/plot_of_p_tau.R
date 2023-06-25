library(ggplot2)
library("gridExtra")
library("cowplot")
library(ggpubr)
library(latex2exp)
library(scales)


p.hole<-0.05
simu.time<-50
#M.anno
#ldsc
#ldsc only
remake.ldsc.res=F
if(remake.ldsc.res==T){
  for (resultN in c(20000)) {
    simu.name<-paste0('new_evalup')
    for (prob.causal in c(0.05)) {
      for (heritability in c(0.4)) {
        for (DHSenrichment in seq(0.5,3,by=0.1)) {
          gwas.out<-paste0('/Users/zewei/local/str-HDL/sim_gwas/',
                           simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment)
          ldsc.out<-paste0('/Users/zewei/local/str-HDL/result/',
                           simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment)
          #start_time <- Sys.time()
          #ldsc pre
          ldsc.ll<-c()
          ldsc.enrichment<-0
          ldsc.enrichment.emp<-NA
          ldsc.taus<-NA
          ldsc.p<-NA
          #simulation start
          for (i in 51:200) {
            #ldsc
            ldsc.temp<-read.table(paste0(ldsc.out,'/ldsc_v2_',i,
                                         '.results'),header = T,stringsAsFactors = F)
            #ldsc.category<-ldsc.temp[,1]
            ldsc.temp[,'Enrich_test']<-ifelse(ldsc.temp$Enrichment_p<p.hole,1,0)
            #ldsc.temp<-as.matrix(ldsc.temp)
            ldsc.enrichment<-ldsc.enrichment+ldsc.temp[,-1]
            ldsc.enrichment.emp<-rbind(ldsc.enrichment.emp,as.vector(ldsc.temp[,5]))
            tau.ldsc<-ldsc.temp$Coefficient
            ldsc.taus<-rbind(ldsc.taus,tau.ldsc)
            p.ldsc<-ldsc.temp$Enrichment_p
            ldsc.p<-rbind(ldsc.p,p.ldsc)
          }
          #ldsc res pre
          ldsc.enrichment<-ldsc.enrichment/simu.time
          ldsc.enrichment.emp<-ldsc.enrichment.emp[-1,]
          ldsc.taus.emp<-ldsc.taus[-1,]
          ldsc.p.emp<-ldsc.p[-1,]
          
          
          
          
          ###analysis ldsc result
          ldsc.h2<-read.table(paste0(ldsc.out,'/ldsc_h2_v2.txt'),header = F,sep = '(')
          
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
          raw.result[['ldsc.p']]=ldsc.p.emp
          saveRDS(raw.result,paste0(gwas.out,'/ldsc.result_v2.Rdata'))
        }
      }
    }
    
  }
}

#analysis
p.analysis<-NA
for (resultN in c(20000)) {
  simu.name<-paste0('new_evalup')
  for (prob.causal in c(0.05)) {
    for (heritability in c(0.4)) {
      for (DHSenrichment in seq(0.5,3,by=0.1)) {
        
        gwas.out<-paste0('/Users/zewei/local/str-HDL/sim_gwas/',
                         simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment)
        result<-readRDS(paste0(gwas.out,'/gls.result.Rdata'))
        ldsc.result<-readRDS(paste0(gwas.out,'/ldsc.result_v2.Rdata'))
        #gls.res<-result$gls.res
        ldsc.enrichment.emp<-ldsc.result$ldsc.res
        ldsc.taus.emp<-ldsc.result$ldsc.taus
        
        M<-ncol(ldsc.enrichment.emp)
        
        gls.p.temp<-ldsc.p.temp<-NA
        total.stat.n<-7
        for (i in 1:150) {
          gls.p.temp<-cbind(gls.p.temp, result[[total.stat.n*(i-1)+7]])
        }
        gls.p.temp<-gls.p.temp[,-1]
        dhs.p.temp<-gls.p.temp['DHS_Trynka',]
        #dhs.p.temp<-pnorm(dhs.p.temp/sd(dhs.p.temp),lower.tail = F)
        
        ldsc.p.temp<-ldsc.result$ldsc.p
        colnames(ldsc.p.temp)<-rownames(gls.p.temp)
        ldsc.p.temp<-ldsc.p.temp[,'DHS_Trynka']
        
        res.temp<-cbind(dhs.p.temp,ldsc.p.temp,rep(DHSenrichment,150))
        res.temp<-as.data.frame(res.temp)
        names(res.temp)<-c('gldsc','sldsc','E')
        p.analysis<-rbind(p.analysis,res.temp)
      }
    }
  }
}

p.analysis<-p.analysis[-1,]
p.analysis[,'g-LDSC']=p.analysis$gldsc<0.05
p.analysis[,'s-LDSC']=p.analysis$sldsc>=0.05

p.analysis.1<-aggregate(p.analysis$`g-LDSC`,list(p.analysis$E),mean)
p.analysis.2<-aggregate(p.analysis$`s-LDSC`,list(p.analysis$E),mean)

p.a.data<-as.data.frame(rbind(p.analysis.1,p.analysis.2))
names(p.a.data)<-c('E','Prob')
p.a.data[,'Methods']<-rep(c('g-LDSC','s-LDSC'),each=nrow(p.analysis.1))
p.a.data[,'text']=ifelse(p.a.data$E==1,percent(round(p.a.data$Prob,3)),'')

#p.a.data[p.a.data$Methods=='g-LDSC','Prob']=1-p.a.data[p.a.data$Methods=='g-LDSC','Prob']
p.a.data=p.a.data[c(1:6,seq(8,20,by=2),27:32,seq(34,46,by=2)),]
plot.pr.of.reject<-ggplot(data=p.a.data, aes(x=E, y=Prob, group=Methods, col=Methods,label=text))+
  geom_line()+ xlim(c(0.4,2.5))+
  xlab('True fold enrichment')+ylab('Pr(rejected at P < 0.05)')+
  geom_point()+theme_light()+
  scale_colour_manual(values = c('#F8766D','#619CFF'))+
  theme(legend.position='none')+theme(legend.title=element_blank())+
  geom_text(size=3,colour='black',nudge_y = 0.025,nudge_x = -0.1)

plot.pr.of.reject





##new table1
simu.time<-50
p.hole<-0.05/24
#ldsc
#ldsc only
remake.ldsc.res=T
if(remake.ldsc.res==T){
  for (resultN in c(20000)) {
    simu.name<-paste0('DHS_newpoly')
    for (prob.causal in c(0.01,0.05,1)) {
      for (heritability in c(0.4)) {
        for (DHSenrichment in c(1,2,3)) {
          gwas.out<-paste0('/Users/zewei/local/str-HDL/sim_gwas/',
                           simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N20000')
          ldsc.out<-paste0('/Users/zewei/local/str-HDL/result/',
                           simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N20000')
          #start_time <- Sys.time()
          #ldsc pre
          ldsc.ll<-c()
          ldsc.enrichment<-0
          ldsc.enrichment.emp<-NA
          ldsc.taus<-NA
          #simulation start
          for (i in 1:simu.time) {
            #ldsc
            ldsc.temp<-read.table(paste0(ldsc.out,'/ldsc_v2_',i,
                                         '.results'),header = T,stringsAsFactors = F)
            #ldsc.category<-ldsc.temp[,1]
            ldsc.temp[,'Enrich_test']<-ifelse(ldsc.temp$Enrichment_p<p.hole,1,0)
            #ldsc.temp<-as.matrix(ldsc.temp)
            ldsc.enrichment<-ldsc.enrichment+ldsc.temp[,-1]
            #ldsc.enrichment.emp<-rbind(ldsc.enrichment.emp,as.vector(ldsc.temp[,5]))
            ldsc.enrichment.emp<-rbind(ldsc.enrichment.emp,as.vector(ldsc.temp[,3]/ldsc.temp[1,3]/M.anno))
            tau.ldsc<-ldsc.temp$Coefficient
            ldsc.taus<-rbind(ldsc.taus,tau.ldsc)
          }
          #ldsc res pre
          ldsc.enrichment<-ldsc.enrichment/simu
          ldsc.enrichment.emp<-ldsc.enrichment.emp[-1,]
          ldsc.taus.emp<-ldsc.taus[-1,]
          
          Sys.time()-start_time
          
          
          ###analysis ldsc result
          ldsc.h2<-read.table(paste0(ldsc.out,'/ldsc_h2_v2.txt'),header = F,sep = '(')
          
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
          saveRDS(raw.result,paste0(gwas.out,'/ldsc.result_v2.Rdata'))
        }
      }
    }
    
  }
}


for (resultN in c(20000)) {
  simu.name<-paste0('DHS_newpoly')
  plot.result.grabber<-NA
  for (prob.causal in c(0.01,0.05,1)) {
    for (heritability in c(0.4)) {
      for (DHSenrichment in c(1,2,3)) {
        ##start to read result from here#####
        gwas.out<-paste0('/Users/zewei/local/str-HDL/sim_gwas/',
                         simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N20000')
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
  assign(paste0('newpoly.N',resultN),plot.result.grabber)
  
}

newpoly.N20000[,'bias2']=(newpoly.N20000$e.box-newpoly.N20000$Enrichment)**2
temp<-aggregate(newpoly.N20000[,10],list(newpoly.N20000[,5],
                                         newpoly.N20000[,7],
                                         newpoly.N20000[,3]),mean)
temp[,'rmse']<-sqrt(temp$x)

temp2<-aggregate(newpoly.N20000[,1],list(newpoly.N20000[,5],
                                         newpoly.N20000[,7],
                                         newpoly.N20000[,3]),sd)

#temp
#temp2

