library(Metrics)

simu.name='Confounding'
confound.res=NULL
total.stat.n<-7
for (N in c(5000,10000,20000)) {
  for (prob.causal in c(0.05)) {
    for (heritability in c(0.4)) {
      for (DHSenrichment in c(1,3)) {
        for (confound in c(1,1.05,1.1,1.2)) {
          gwas.out<-paste0('/Users/zewei/local/str-HDL/sim_gwas/',
                           simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N',N,'I',confound)
          ldsc.out<-paste0('/Users/zewei/local/str-HDL/result/',
                           simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N',N,'I',confound)
          
          ldsc.a<-read.table(paste0(ldsc.out,'/ldsc_inter_v2.txt'),header = F,sep = '(')
          result<-readRDS(paste0(gwas.out,'/gls.result.Rdata'))
          ldsc.confound<-gls.confound<-c()
          for (i in 1:nrow(ldsc.a)) {
            ldsc.confound=c(ldsc.confound,as.numeric(strsplit(ldsc.a[i,1], ': ')[[1]][[2]]))
            gls.confound=c(gls.confound,result[[total.stat.n*(i-1)+4]])
          }
          bias=c(gls.confound,ldsc.confound)
          Methods=rep(c('g-LDSC','s-LDSC'),each=nrow(ldsc.a))
          data.temp=as.data.frame(cbind(bias,Methods))
          data.temp$N=N
          data.temp$Poly=prob.causal
          data.temp$h2=heritability
          data.temp$E=DHSenrichment
          data.temp$true_bias=confound
          
          confound.res=rbind(confound.res,data.temp)
        }
      }
    }
  }
}
confound.res$N=paste0('n = ',confound.res$N)
confound.res$N<-factor(confound.res$N,levels = c('n = 5000','n = 10000','n = 20000'))
confound.res$bias=as.numeric(confound.res$bias)

#plot based on different true bias
plot.data=confound.res[confound.res$true_bias==1,]
plot.bias.1<-ggplot(data=plot.data,aes(x=true_bias,y=bias, fill=Methods))+
  geom_boxplot()+facet_wrap(. ~ `N`,nrow = 1)+
  stat_summary(fun.y = mean, color = "darkred", position = position_dodge(0.75),
               geom = "point", shape = c(18), size = 3,
               show.legend = FALSE)+geom_abline(slope=0, intercept=c(1),linetype = "dashed")+
  xlab('Causal rate')+ylab('Estimated Inflation')+
  theme_light()+
  scale_fill_manual(values = c('#F8766D','#619CFF'))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position='top')+theme(legend.title=element_blank())

plot.bias.1
#plot of inflation 1.1
plot.data=confound.res[confound.res$true_bias==1.1,]
plot.bias.1.1<-ggplot(data=plot.data,aes(x=true_bias,y=bias, fill=Methods))+
  geom_boxplot()+facet_wrap(. ~ `N`,nrow = 1)+
  stat_summary(fun.y = mean, color = "darkred", position = position_dodge(0.75),
               geom = "point", shape = c(18), size = 3,
               show.legend = FALSE)+geom_abline(slope=0, intercept=c(1.1^2),linetype = "dashed")+
  xlab('Causal rate')+ylab('Estimated Inflation')+
  theme_light()+
  scale_fill_manual(values = c('#F8766D','#619CFF'))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position='top')+theme(legend.title=element_blank())

plot.bias.1.1
#inflation 1.2
plot.data=confound.res[confound.res$true_bias==1.2,]
plot.bias.1.2<-ggplot(data=plot.data,aes(x=true_bias,y=bias, fill=Methods))+
  geom_boxplot()+facet_wrap(. ~ `N`,nrow = 1)+
  stat_summary(fun.y = mean, color = "darkred", position = position_dodge(0.75),
               geom = "point", shape = c(18), size = 3,
               show.legend = FALSE)+geom_abline(slope=0, intercept=c(1.2),linetype = "dashed")+
  xlab('Causal rate')+ylab('Estimated Inflation')+
  theme_light()+
  scale_fill_manual(values = c('#F8766D','#619CFF'))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position='top')+theme(legend.title=element_blank())

plot.bias.1.2

#inflation 1.05
plot.data=confound.res[confound.res$true_bias==1.05,]
plot.bias.1.05<-ggplot(data=plot.data,aes(x=true_bias,y=bias, fill=Methods))+
  geom_boxplot()+facet_wrap(. ~ `N`,nrow = 1)+
  stat_summary(fun.y = mean, color = "darkred", position = position_dodge(0.75),
               geom = "point", shape = c(18), size = 3,
               show.legend = FALSE)+geom_abline(slope=0, intercept=c(1.05^2),linetype = "dashed")+
  xlab('Causal rate')+ylab('Estimated Inflation')+
  theme_light()+
  scale_fill_manual(values = c('#F8766D','#619CFF'))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position='top')+theme(legend.title=element_blank())

plot.bias.1.05



#E of inflation



confound.res=NULL
total.stat.n<-7
load('/Users/zewei/Desktop/HKUgrad/G-LDSC/plots/plot.data.v2.Rdata')
simu.name='Confounding'
p.hole=24/0.05
simu.time=50
remake.ldsc.res=F
if(remake.ldsc.res){
  for (N in c(500,5000,10000,20000)) {
    for (prob.causal in c(0.05)) {
      for (heritability in c(0.4)) {
        for (DHSenrichment in c(1,3)) {
          for (confound in c(1,1.05,1.1)) {
            gwas.out<-paste0('/Users/zewei/local/str-HDL/sim_gwas/',
                             simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N',N,'I',confound)
            ldsc.out<-paste0('/Users/zewei/local/str-HDL/result/',
                             simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N',N,'I',confound)
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
}



for (N in c(500,5000,10000,20000)) {
  plot.result.grabber<-NA
  for (prob.causal in c(0.05)) {
    for (heritability in c(0.4)) {
      for (DHSenrichment in c(1,3)) {
        for (confound in c(1,1.05,1.1)) {
          gwas.out<-paste0('/Users/zewei/local/str-HDL/sim_gwas/',
                           simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N',N,'I',confound)
          ldsc.out<-paste0('/Users/zewei/local/str-HDL/result/',
                           simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N',N,'I',confound)
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
          plot.result.grabber.temp[,'Inflation']=confound
          #add tau
          plot.result.grabber.temp[,'Est.tau1']=box.data[box.data$cate.box==1,'tau.box']
          plot.result.grabber<-rbind(plot.result.grabber,plot.result.grabber.temp)
        }
      }
    }
  }
  plot.result.grabber<-na.omit(plot.result.grabber)
  assign(paste0('plot.confound',N),plot.result.grabber)
}
plot.confound.data<-rbind(plot.confound500,plot.confound5000,plot.confound10000,plot.confound20000)
plot.confound.data$group.box[plot.confound.data$group.box=='gls']='g-LDSC'
plot.confound.data$group.box[plot.confound.data$group.box=='ldsc']='s-LDSC'
plot.confound.data[,'N']<-rep(c('n = 5,00','n = 5000','n = 10,000','n = 20,000'
                         ),times=1,each=800)
plot.confound.data[,'N']<-factor(plot.confound.data[,'N'],levels = c('n = 500','n = 5,000',
                                                       'n = 10,000','n = 20,000'))
#plot.confound.data$Inflation=round(plot.confound.data$Inflation^2,2)
plot.confound.data$Inflation=rep(c('Lambda = 1','Lambda = 1.1','Lambda = 1.21', 'Lambda = 1.44'),each=100,times=8)
plot.confound.data$Inflation=factor(plot.confound.data$Inflation,levels = c('Lambda = 1','Lambda = 1.1',
                                                                            'Lambda = 1.21', 'Lambda = 1.44'))


#control of show how many N
plot.data.confounding<-plot.confound.data[which(plot.confound.data$h2=='0.4'&
                                plot.confound.data$causal=='0.05'&plot.confound.data$N=='n = 20,000'),]
plot.h2.confound<-ggplot(data=plot.confound.data,aes(x=h2,y=Est.h2, fill=group.box))+
  geom_boxplot()+facet_wrap(. ~ `Inflation`,nrow = 1)+
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

plot.h2.confound        

#plto of e
mean.data<-aggregate(plot.data.confounding[,c(1)],list(plot.data.confounding$group.box,plot.data.confounding$Inflation,
                                            plot.data.confounding$Enrichment),mean)
sd.data<-aggregate(plot.data.confounding[,c(1)],list(plot.data.confounding$group.box,plot.data.confounding$Inflation,
                                          plot.data.confounding$Enrichment),sd)
names(mean.data)=c('Method','Inflation','True enrichment','Est enrichment')
mean.data$Method[mean.data$Method=='gls']='g-LDSC'
mean.data$Method[mean.data$Method=='ldsc']='s-LDSC'


plot.e.confound=ggplot(mean.data, aes(x=`True enrichment`, y=`Est enrichment`, color=Method))+
  geom_point(size = 1.5,position=position_dodge(0.5))+ facet_grid(. ~ `Inflation`)+xlim(0,4)+ylim(-1,6)+
  #geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.1))+
  geom_errorbar(aes(ymin=`Est enrichment`-1.96*sd.data$x, ymax=`Est enrichment`+1.96*sd.data$x), width=.2,
                position=position_dodge(0.5))+geom_abline(slope=1, intercept=0,linetype = "dashed")+
  panel_border()+xlab('True fold enrichment')+ylab('Estimated fold enrichment')+
  scale_fill_discrete(name = "Methods")+
  theme_light()+
  #theme_bw()+background_grid(major = 'none', minor = "none")+
  scale_colour_manual(values = c('#F8766D','#619CFF'))+
  theme(legend.position='top')+theme(legend.title=element_blank())

plot.e.confound

#lambda analysis

simu.name='Confounding'
#lambda.res=NULL
total.stat.n<-7
for (N in c(500,5000,10000,20000)) {
  plot.result.grabber<-NA
  for (prob.causal in c(0.05)) {
    for (heritability in c(0.4)) {
      for (DHSenrichment in c(1,3)) {
        for (confound in c(1,1.05,1.1,1.2)) {
          gwas.out<-paste0('/Users/zewei/local/str-HDL/sim_gwas/',
                           simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N',N,'I',confound)
          ldsc.out<-paste0('/Users/zewei/local/str-HDL/result/',
                           simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N',N,'I',confound)
          result<-readRDS(paste0(gwas.out,'/gls.result_v2.Rdata'))
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
          plot.result.grabber.temp[,'Inflation']=confound
          #add tau
          plot.result.grabber.temp[,'Est.tau1']=box.data[box.data$cate.box==1,'tau.box']
          plot.result.grabber<-rbind(plot.result.grabber,plot.result.grabber.temp)
        }
      }
    }
  }
  plot.result.grabber<-na.omit(plot.result.grabber)
  assign(paste0('plot.lambda',N),plot.result.grabber)
}
plot.lambda.data<-rbind(plot.lambda500,plot.lambda5000,plot.lambda10000,plot.lambda20000)
plot.lambda.data$group.box[plot.lambda.data$group.box=='gls']='g-LDSC'
plot.lambda.data$group.box[plot.lambda.data$group.box=='ldsc']='s-LDSC'
plot.lambda.data[,'N']<-rep(c('n = 20,000','n = 10,000',
                                'n = 5,000','n = 500'),times=1,each=800)
plot.lambda.data[,'N']<-factor(plot.lambda.data[,'N'],levels = c('n = 500','n = 5,000',
                                                                     'n = 10,000','n = 20,000'))

aggregate(plot.confound.data$e.box,list(plot.confound.data$group.box,plot.confound.data$Inflation,
                                        plot.confound.data$N),mean)
aggregate(plot.lambda.data$e.box,list(plot.lambda.data$group.box,plot.lambda.data$Inflation,
                                        plot.lambda.data$N),mean)
View(cbind(aggregate(plot.confound.data$e.box,list(plot.confound.data$group.box,plot.confound.data$Inflation,
                                                   plot.confound.data$N),mean),
           aggregate(plot.lambda.data$e.box,list(plot.lambda.data$group.box,plot.lambda.data$Inflation,
                                                 plot.lambda.data$N),mean)))


pambda.data1=plot.confound.data[plot.confound.data$group.box=='g-LDSC',]
pambda.data2=plot.lambda.data[plot.lambda.data$group.box=='g-LDSC',]

View(cbind(aggregate(pambda.data1$e.box,list(pambda.data1$Enrichment,pambda.data1$Inflation,
                                             pambda.data1$N),mean),
           aggregate(pambda.data2$e.box,list(pambda.data2$Enrichment,pambda.data2$Inflation,
                                             pambda.data2$N),mean)))


View(cbind(aggregate(pambda.data1$e.box,list(pambda.data1$Enrichment,pambda.data1$Inflation,
                                             pambda.data1$N),sd),
           aggregate(pambda.data2$e.box,list(pambda.data2$Enrichment,pambda.data2$Inflation,
                                             pambda.data2$N),sd)))

#detail of confound


simu.name='Confounding'
#lambda.res=NULL
total.stat.n<-7
plot.of.lambdaE=NULL
for (N in c(2000,10000,50000)){
  for (prob.causal in c(0.05)) {
    for (heritability in c(0.4)) {
      for (DHSenrichment in c(3)) {
        for (alpha in c(0,0.2/10000,1/10000,5/10000)) {
          #set.seed(seedno+N)
          confound=1+N*alpha
          gwas.out<-paste0('/Users/zewei/local/str-HDL/sim_gwas/',
                           simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N',N,'I',confound)
          ldsc.out<-paste0('/Users/zewei/local/str-HDL/result/',
                           simu.name,'c',prob.causal,'h',heritability,'e',DHSenrichment,'N',N,'I',confound)
          result<-readRDS(paste0(gwas.out,'/gls.result.Rdata'))
          result2<-readRDS(paste0(gwas.out,'/gls.result_v2.Rdata'))
          
          gls.e.res1<-NULL
          gls.e.res2<-NULL
          gls.lambda1<-NULL
          gls.lambda2<-NULL
          for (i in 1:simu.time) {
            gls.e.res1<-c(gls.e.res1, result[[total.stat.n*(i-1)+3]][11])
            gls.e.res2<-c(gls.e.res2, result2[[total.stat.n*(i-1)+3]][11])
            gls.lambda1<-c( gls.lambda1,result[[total.stat.n*(i-1)+4]])
            gls.lambda2<-c( gls.lambda2,result2[[total.stat.n*(i-1)+4]])
          }
          plot.of.lambda.temp=as.data.frame(cbind(gls.e.res1,gls.e.res2))
          plot.of.lambda.temp[,'N']=N
          plot.of.lambda.temp[,'h2']=heritability
          plot.of.lambda.temp[,'Poly']=prob.causal
          plot.of.lambda.temp[,'True_E']=DHSenrichment
          plot.of.lambda.temp[,'alpha']=alpha
          plot.of.lambda.temp[,'lambda1']=gls.lambda1
          plot.of.lambda.temp[,'lambda2']=gls.lambda2
          
          plot.of.lambdaE=rbind(plot.of.lambdaE,plot.of.lambda.temp)
        }
      }
    }
  }
}
names(plot.of.lambdaE)=c('Org','Iter','N','h2','poly','True_E','alpha','lambda1','lambda2')
plot(plot.of.lambdaE$Org~plot.of.lambdaE$Iter)

plot.of.lambdaE[,'Class']=paste0(plot.of.lambdaE$N,'_',plot.of.lambdaE$alpha)

ggplot(plot.of.lambdaE, aes(x=Org, y=Iter)) +
  geom_point(size=1.5)+ geom_abline(slope=1, intercept=c(0),linetype = "dashed")+facet_wrap(. ~ `Class`,nrow = 3)+
  geom_smooth(method = "lm", fill = NA)+coord_cartesian(ylim=c(-2, 10),xlim = c(-4,10))


aggregate(plot.of.lambdaE$Org,list(plot.of.lambdaE$N,plot.of.lambdaE$Confound,plot.of.lambdaE$True_E),mean)
aggregate(plot.of.lambdaE$Iter,list(plot.of.lambdaE$N,plot.of.lambdaE$Confound,plot.of.lambdaE$True_E),mean)


View(cbind(aggregate(plot.of.lambdaE$Org,list(plot.of.lambdaE$alpha,plot.of.lambdaE$True_E,plot.of.lambdaE$N),mean),
           aggregate(plot.of.lambdaE$Iter,list(plot.of.lambdaE$alpha,plot.of.lambdaE$True_E,plot.of.lambdaE$N),mean)[,4]))


View(cbind(aggregate(plot.of.lambdaE$Org,list(plot.of.lambdaE$alpha,plot.of.lambdaE$True_E,plot.of.lambdaE$N),sd),
           aggregate(plot.of.lambdaE$Iter,list(plot.of.lambdaE$alpha,plot.of.lambdaE$True_E,plot.of.lambdaE$N),sd)[,4]))



tables4_1=cbind(aggregate(plot.of.lambdaE$Org,list(plot.of.lambdaE$alpha,plot.of.lambdaE$True_E,plot.of.lambdaE$N),rmse,actual=3),
                aggregate(plot.of.lambdaE$Iter,list(plot.of.lambdaE$alpha,plot.of.lambdaE$True_E,plot.of.lambdaE$N),rmse,actual=3)[,4])


tables4_2=cbind(aggregate(plot.of.lambdaE$Org,list(plot.of.lambdaE$alpha,plot.of.lambdaE$True_E,plot.of.lambdaE$N),sd),
                aggregate(plot.of.lambdaE$Iter,list(plot.of.lambdaE$alpha,plot.of.lambdaE$True_E,plot.of.lambdaE$N),sd)[,4])

tables4=cbind(tables4_1,tables4_2[,4:5])
tables4=as.data.frame(tables4)
for (i in 4:7) {
  tables4[,i]=round(tables4[,i],3)
}

names(tables4)=c('Alpha','True_Enrichment','N','g-LDSC RMSE','g-LDSC_iter RMSE','g-LDSC SE','g-LDSC_iter SE')
tables4[,'Lambda']=1+tables4$Alpha*tables4$N

tables4=tables4[c('Alpha','N','Lambda','g-LDSC RMSE','g-LDSC SE','g-LDSC_iter RMSE','g-LDSC_iter SE')]

write.xlsx(tables4,'/Users/zewei/Desktop/HKUgrad/G-LDSC/plots/figure/TableS4.xlsx')







temp1=cbind(aggregate(plot.of.lambdaE$lambda1,list(plot.of.lambdaE$alpha,plot.of.lambdaE$True_E,plot.of.lambdaE$N),mean),
           aggregate(plot.of.lambdaE$lambda2,list(plot.of.lambdaE$alpha,plot.of.lambdaE$True_E,plot.of.lambdaE$N),mean)[,4])
temp1[,'true_lambda']=(temp1$Group.3*temp1$Group.1+1)
View(temp1)


summary(plot.of.lambdaE$Iter[plot.of.lambdaE$N==10000&plot.of.lambdaE$alpha==5e-4])
