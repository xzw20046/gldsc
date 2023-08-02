##
library(Rcpp)
library(RcppEigen)
library(inline)
library(emulator)
library(doParallel)
library(rhdf5)


#gwas.error<-read.table('/Users/zewei/local/str-HDL/data/error_estimate_phase1.qassoc.ref.txt',header = T)
#gwas.error<-na.omit(gwas.error)




#gwas.error<-gwas.df#no need to run for the second time
#load('/Users/zewei/local/str-HDL/data/chr10ref.DSout.pannel.Rdata')
load('/Users/zewei/local/str-HDL/data/chr10ref.pannel.Rdata')
source('/Users/zewei/Desktop/software/str-HDL/mlfun.v2.r')

gwas.path<-'/Users/zewei/local/GWAS_error_50times'
simu.name='randomE'


##
heritability<-0.4
heritability
prob.causal<-0.05
nref=379
N<-20000
simu<-50
DHSenrichment=1

Nanno<-length(M.anno)

total.snp<-unlist(snp.list)
M<-length(total.snp)
gwas.df2<-gwas.df[match(total.snp,gwas.df[,'SNP']),]

Atotal2<-read.table('/Users/zewei/local/str-HDL/data/53Amatrix.txt')
P=ncol(Atotal2)
##randomly assign enrichment and causal SNP
seedno<-202301023
set.seed(seedno)
#set.seed(20211105)
#loopstart here
for (N in c(5000,10000,20000)){
  for (prob.causal in c(0.05)) {
    for (heritability in c(0.4)) {
        #set.seed(seedno+N)
        gwas.out<-paste0('/Users/zewei/local/str-HDL/sim_gwas/',
                         simu.name,'c',prob.causal,'h',heritability,'N',N)
        ldsc.out<-paste0('/Users/zewei/local/str-HDL/result/',
                         simu.name,'c',prob.causal,'h',heritability,'N',N)
        system(paste0('mkdir ',gwas.out))
        system(paste0('mkdir ',ldsc.out))
        #

        Atotal.temp2<-Atotal2
        


          causal.snp.loc<-sample(1:M, round(prob.causal*M))
          
          #
          Atotal.temp<-Atotal.temp2[causal.snp.loc,]
          Atotal.temp<-as.matrix(Atotal.temp) # causal A matrix
          Atotal.temp2<-as.matrix(Atotal.temp2)#full A matrix
          anno.total.temp<-apply(Atotal.temp2, 2, sum)
          #M.anno.temp<-anno.total.temp/round(prob.causal*M)
          M.anno.temp<-anno.total.temp/round(M)
          
          #generate true taus
          raw.tau=100^rbinom(P,1,0.1)*runif(P)
          raw.hv=crossprod(Atotal.temp)%*%raw.tau
          Vh=heritability*raw.hv/raw.hv[1]
          tau=heritability*raw.tau/raw.hv[1]
          
          enrichment=Vh/heritability/M.anno
          enrichment
          
          per.snp.her=Atotal.temp%*%tau
          eff.var<-per.snp.her
          if(sum(eff.var<0)){stop('bad')}
          
          #simulation joint effect
          true.effect<-rep(0,length(total.snp))
          true.effect[causal.snp.loc]<-eff.var
          
          Vh.true=t(Atotal)%*%true.effect
          h2.true=Vh.true[1]
          enrichment.true<-Vh.true/h2.true/M.anno
          
          
          #simulation starts here
          sample.h2.res<-NA
          sample.e.res<-NA
          sample.stat<-NA
          #set.seed(20211105)
          for (simu.time in 1:simu) {
          joint.effect.size<-NA
          for (i in 1:length(true.effect)) {
            joint.effect.size[i]<-rnorm(1,0,sqrt(true.effect[i]))
          }
          
          #simulation GWAS
          margin.eff.raw<-c()
          for (batch in 1:length(snp.list)) {
            snp.pos<-match(snp.list[[batch]],total.snp)
            eff.temp<-joint.effect.size[snp.pos]
            margin.ef.temp<-LDmatrix[[batch]]%*%eff.temp
            margin.eff.raw<-c(margin.eff.raw,margin.ef.temp)
            rm(snp.pos);rm(eff.temp);rm(margin.ef.temp)
          }
          #gwas.error.temp<-gwas.error %>% filter(SNP %in% total.snp)
          gwas.error<-read.table(paste0(gwas.path,'/gwas_',simu.time,'.qassoc'),header = T)
          
          gwas.error.temp<-gwas.error[match(total.snp,gwas.error[,'SNP']),]
          #old version
          random.eff<-gwas.error.temp[,'T']*sqrt(1/N)
          sim.gwas<-gwas.error.temp
          sim.gwas[,'BETA']<-margin.eff.raw+random.eff
          sim.gwas[,'NMISS']<-rep(N,nrow(sim.gwas))
          sim.gwas[,'SE']<-rep(sqrt(1/N),nrow(sim.gwas))
          sim.gwas<-sim.gwas[,-match(c('R2','T','P'),names(sim.gwas))]
          sim.gwas[,'Z']<-sim.gwas[,'BETA']/sim.gwas[,'SE']
          sim.gwas[,'A1']<-gwas.df2[,'A1']
          sim.gwas[,'A2']<-gwas.df2[,'A2']
          #gwas.df.old<-gwas.df
          gwas.df<-sim.gwas[,c('SNP','A1','A2','Z','NMISS')]
          names(gwas.df)<-c('SNP','A1','A2','Z','N')
          #var(gwas.df$Z)
          #summary(gwas.df$Z)
          #new version
          
          
          
          write.table(gwas.df,paste0(gwas.out,'/chr10simgwas_test',simu.time,'.sumstat'),
                      row.names = F,quote = F)
          write.table(enrichment,paste0(gwas.out,'/true_E_',simu.time,'.txt'),
                      row.names = F,quote = F)
        }
        #sample.h2.res<-sample.h2.res[,-1]
        #sample.e.res<-sample.e.res[,-1]
        #sample.stat<-sample.stat[,-1]
        print('Done 1')
    }
  }
}



