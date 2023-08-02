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
            ldsc.temp<-read.table(paste0(ldsc.out,'/ldsc_v2_',i,
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

plot.data<-plot.4.data1[which(plot.4.data1$h2=='0.4'&
                               plot.4.data1$causal=='0.05'),]



plot.h2.N<-ggplot(data=plot.data,aes(x=h2,y=Est.h2, fill=group.box))+
  geom_boxplot()+facet_wrap(. ~ `N`,nrow = 1)+
  stat_summary(fun.y = mean, color = "darkred", position = position_dodge(0.75),
               geom = "point", shape = c(18), size = 3,
               show.legend = FALSE)+geom_abline(slope=0, intercept=c(0.4),linetype = "dashed")+
  xlab('Causal rate')+ylab('Estimated heritability')+
  theme_light()+
  scale_fill_manual(values = c('#F8766D','#619CFF'))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position='top')+theme(legend.title=element_blank())

plot.data<-plot.4.data1[which(plot.4.data1$h2=='0.04'&
                               plot.4.data1$causal=='0.05'),]
plot.h2.N.004<-ggplot(data=plot.data,aes(x=h2,y=Est.h2, fill=group.box))+
  geom_boxplot()+facet_wrap(. ~ `N`,nrow = 1)+
  stat_summary(fun.y = mean, color = "darkred", position = position_dodge(0.75),
               geom = "point", shape = c(18), size = 3,
               show.legend = FALSE)+geom_abline(slope=0, intercept=c(0.04),linetype = "dashed")+
  xlab('Causal rate')+ylab('Estimated heritability')+
  theme_light()+
  scale_fill_manual(values = c('#F8766D','#619CFF'))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position='top')+theme(legend.title=element_blank())
plot.h2.N.004

#plot e vs N
plot.data2<-plot.4.data2[which(plot.4.data2$h2=='0.4'&
                                 plot.4.data2$causal=='0.05'),]
#plot.data[,'N']<-c(rep('n = 50,000',1800),rep('n = 25,000',1800),rep('n = 20,000',1800),rep('n = 15,000',1800),
#                   rep('n = 10,000',1800),rep('n = 5,000',1800))
#plot.data[,'N']<-factor(plot.data[,'N'],levels = c('n = 5,000','n = 10,000','n = 15,000','n = 20,000','n = 25,000','n = 50,000'))
#plot.data<-plot.data[plot.data$h2=='0.4'&
#                       plot.data$causal=='0.05'&(plot.data$N!='n = 25,000'&plot.data$N!='n = 50,000'),]

mean.data<-aggregate(plot.data2[,c(1)],list(plot.data2$group.box,plot.data2$N,
                                           plot.data2$Enrichment),mean)
sd.data<-aggregate(plot.data2[,c(1)],list(plot.data2$group.box,plot.data2$N,
                                         plot.data2$Enrichment),sd)
names(mean.data)=c('Method','Sample size','True enrichment','Est enrichment')
mean.data$Method[mean.data$Method=='gls']='g-LDSC'
mean.data$Method[mean.data$Method=='ldsc']='s-LDSC'


plot.e.N=ggplot(mean.data, aes(x=`True enrichment`, y=`Est enrichment`, color=Method))+
  geom_point(size = 1.5,position=position_dodge(0.5))+ facet_grid(. ~ `Sample size`)+xlim(0,4)+ylim(-1,6)+
  #geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.1))+
  geom_errorbar(aes(ymin=`Est enrichment`-1.96*sd.data$x, ymax=`Est enrichment`+1.96*sd.data$x), width=.2,
                position=position_dodge(0.5))+geom_abline(slope=1, intercept=0,linetype = "dashed")+
  panel_border()+xlab('True fold enrichment')+ylab('Estimated fold enrichment')+
  scale_fill_discrete(name = "Methods")+
  theme_light()+
  #theme_bw()+background_grid(major = 'none', minor = "none")+
  scale_colour_manual(values = c('#F8766D','#619CFF'))+
  theme(legend.position='top')+theme(legend.title=element_blank())

plot.e.N

plot.4.data2<-plot.4.data[which(plot.4.data$N %in% c('n = 2,000','n = 20,000',
                                                     'n = 10,000','n = 50,000')),]
plot.data2<-plot.4.data2[which(plot.4.data2$h2=='0.04'&
                                 plot.4.data2$causal=='0.05'),]

mean.data2<-aggregate(plot.data2[,c(1)],list(plot.data2$group.box,plot.data2$N,
                                            plot.data2$Enrichment),mean)
sd.data2<-aggregate(plot.data2[,c(1)],list(plot.data2$group.box,plot.data2$N,
                                          plot.data2$Enrichment),sd)
names(mean.data2)=c('Method','Sample size','True enrichment','Est enrichment')
mean.data2$Method[mean.data2$Method=='gls']='g-LDSC'
mean.data2$Method[mean.data2$Method=='ldsc']='s-LDSC'


plot.e.N.004=ggplot(mean.data2, aes(x=`True enrichment`, y=`Est enrichment`, color=Method))+
  geom_point(size = 1.5,position=position_dodge(0.5))+ facet_grid(. ~ `Sample size`)+xlim(0,4)+#ylim(-6,8)+
  #geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.1))+
  geom_errorbar(aes(ymin=`Est enrichment`-1.96*sd.data2$x, ymax=`Est enrichment`+1.96*sd.data2$x), width=.2,
                position=position_dodge(0.5))+geom_abline(slope=1, intercept=0,linetype = "dashed")+
  panel_border()+xlab('True fold enrichment')+ylab('Estimated fold enrichment')+
  scale_fill_discrete(name = "Methods")+
  theme_light()+
  #theme_bw()+background_grid(major = 'none', minor = "none")+
  scale_colour_manual(values = c('#F8766D','#619CFF'))+
  theme(legend.position='top')+theme(legend.title=element_blank())

plot.e.N.004
plot.e.N
plot.h2.N
plot.h2.N.004

