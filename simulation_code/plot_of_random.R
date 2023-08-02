
simu.name='randomE'
remake.ldsc.res=T
if(remake.ldsc.res){
  for (N in c(5000,10000,20000)){
    for (prob.causal in c(0.05)) {
      for (heritability in c(0.4)) {
  gwas.out<-paste0('/Users/zewei/local/str-HDL/sim_gwas/',
                   simu.name,'c',prob.causal,'h',heritability,'N',N)
  ldsc.out<-paste0('/Users/zewei/local/str-HDL/result/',
                   simu.name,'c',prob.causal,'h',heritability,'N',N)
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
    }
  }
}

#pull resutl
result.random=NULL
for (N in c(5000,10000,20000)){
  for (prob.causal in c(0.05)) {
    for (heritability in c(0.4)) {
      gwas.out<-paste0('/Users/zewei/local/str-HDL/sim_gwas/',
                       simu.name,'c',prob.causal,'h',heritability,'N',N)
      ldsc.out<-paste0('/Users/zewei/local/str-HDL/result/',
                       simu.name,'c',prob.causal,'h',heritability,'N',N)
      result<-readRDS(paste0(gwas.out,'/gls.result.Rdata'))
      ldsc.result<-readRDS(paste0(gwas.out,'/ldsc.result_v2.Rdata'))
      
      total.stat.n<-7
      for (i in 1:simu) {
        ldsc.temp=ldsc.result[[i]]
        test1=as.vector(ldsc.temp[,3]/ldsc.temp[1,3]/M.anno)
        true.e.temp=as.vector(read.table(paste0(gwas.out,'/true_E_',simu.time,'.txt'),header = T)[,1])
        
        temp.res=cbind(c(as.vector(result[[total.stat.n*(i-1)+3]]),test1),
                       rep(true.e.temp,2),
                       rep(rownames(result[[total.stat.n*(i-1)+3]]),2),
                       rep(c('g-LDSC','s-LDSC'),each=length(test1))
        )
        temp.res=as.data.frame(temp.res)
        names(temp.res)=c('Est_E','True_E','Anno','Method')
        temp.res[,'N']=paste0('n = ',N)
        temp.res[,'H2']=heritability
        temp.res[,'Poly']=prob.causal
        
        result.random=rbind(result.random,temp.res)
      }

    }
  }
}
result.random$Est_E=as.numeric(result.random$Est_E)
result.random$True_E=as.numeric(result.random$True_E)
result.random$N<-factor(result.random$N,levels = c('n = 5000','n = 10000','n = 20000'))



random.e.plot<-ggplot(result.random, aes(x=True_E, y=Est_E, color=Method, shape=Method)) +
  geom_point(size=1.5)+ geom_abline(slope=1, intercept=c(0),linetype = "dashed")+facet_wrap(. ~ `N`,nrow = 1)+
  geom_smooth(method = "lm", fill = NA)+coord_cartesian(ylim=c(-2, 6))
  
random.e.plot

plot.data=result.random[result.random$N=='n = 10000',]
output1=aggregate(plot.data$Est_E,list(plot.data$Anno,plot.data$Method),mean)
output2=aggregate(plot.data$Est_E,list(plot.data$Anno,plot.data$Method),sd)

Anno=plot.data$Anno[1:53]
True_E=round(as.numeric(plot.data$True_E[1:53]),3)
g_E=round(output1[1:53,3],3)
s_E=round(output1[-1:-53,3],3)
g_SD=round(output2[1:53,3],3)
s_SD=round(output2[-1:-53,3],3)

tables2=as.data.frame(cbind(Anno,True_E,g_E,g_SD,s_E,s_SD))
names(tables2)=c('Annotation','True Enrichment','g-LDSC Enrich','g-LDSC SD','s-LDSC Enrich','s-LDSC SD')
write.xlsx(tables2,'/Users/zewei/Desktop/HKUgrad/G-LDSC/plots/random_table.xlsx',row.names = F,col.names = T)
