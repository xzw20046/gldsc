#' @export
crossprodCpp<-'using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Lower;
const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
const int m(A.rows()), n(A.cols());
MatrixXd AtA(MatrixXd(n, n).setZero().
selfadjointView<Lower>().rankUpdate(A.adjoint()));
//MatrixXd AAt(MatrixXd(m, m).setZero().
//selfadjointView<Lower>().rankUpdate(A));
return wrap(AtA);'

fcprd <- cxxfunction(signature(AA = "matrix"), crossprodCpp, "RcppEigen")

#
LDSCmatrix<-function(LDmatrix,Amatrix){
  cat<-ncol(Amatrix)
  res<-diag(LDmatrix)
  for (an in 1:cat){
    Rnew<-(LDmatrix**2)%*%as.vector(Amatrix[,an])
    res<-cbind(res,Rnew)
  }
  return(as.data.frame(res))
}

##LDSM function V2####
LDSMsum2<-function(LDmatrix,Amatrix){
  middle<-apply(Amatrix, 1, sum)
  Rnew<-LDmatrix*as.vector(sqrt(middle))
  res<-fcprd(Rnew)
  return(res)
}

#LDSM manual tau version
LDSMman<-function(LDmatrix,betavar){
  Rnew<-LDmatrix*(sqrt(betavar))
  res<-fcprd(Rnew)
  return(res)
}


##eigen cut function
eigen.cut<-function(vector,percentage){
  total<-sum(vector)
  k=1
  while (sum(vector[1:k])<percentage*total) {
    k=k+1
  }
  return(k)
}

gls.left.right<-function(i,x,y,block.names.all,LDmatrix,raw.Ntau,
                         LDSM,snp.list,cutpoint=0.99,intercept=T){
  block.name<-block.names.all[i]
  snp.list.temp<-snp.list[[block.name]]
  
  used.snp<-intersect(snp.list.temp,y$SNP)
  
  #message(paste0(i,' OK'))
  if(length(used.snp)<2){return(list(0,0))} else{
    used.loc<-match(used.snp, snp.list.temp)
    y.temp<-as.vector(y[match(used.snp,y$SNP),'Z']**2-1)
    #message(paste0(i,' Y OK'))
    x.temp<-as.matrix(x[[block.name]])[used.loc,]
    #message(paste0(i,' X OK'))
    #raw.Ntau<-mean(y.temp)/mean(apply(x.temp, 1, sum))
    
    if(intercept==T){x.temp<-cbind(rep(1,nrow(x.temp)),x.temp)}
    
    
    ldsm.temp<-LDSM[[block.name]][used.loc,used.loc]
    #message(paste0(i,'ldsm OK'))
    ldsm.temp2<-(raw.Ntau*ldsm.temp+LDmatrix[[block.name]][used.loc,used.loc])**2
    #message(paste0(i,'ldsm2 OK'))
    eigdi<-eigen(ldsm.temp2,symmetric = T)
    eigval<-eigdi$values
    #
    cut.point<-eigen.cut(eigval,cutpoint)
    eigval.inv<-(1/eigval[1:cut.point])
    
    omiga.inv<-crossprod(t(eigdi$vectors[,1:cut.point])*sqrt(eigval.inv))
    
    omiga.inv2<-t(x.temp)%*%omiga.inv
    
    gls.left<-omiga.inv2%*%x.temp
    
    gls.right<-as.vector(omiga.inv2%*%y.temp)
    
    return(list(gls.left,gls.right))
  }
}


##input left right output set of taus
gls.estimator<-function(left,right,A,M.anno,anno.total,N,intercept=T){
  #point estimator
  if(intercept==T){
    test.tau<-solve(left)%*%right
    
    tau.all<-as.vector(test.tau)/N
    tau<-tau.all[-1]
    constant=N*tau.all[1]+1
    
    est.h<-crossprod(A)%*%tau
    h.total<-est.h[1]
    e<-est.h/M.anno/h.total
    
    #hypothesis test
    M<-anno.total[1]
    h0.test<-(est.h/anno.total)-(h.total-est.h)/(M-anno.total)
    
    #variance in theory
    cov.tau.all<-solve(left)/N^2
    cov.tau<-cov.tau.all[-1,-1]
    constant.var=N^2*cov.tau.all[1,1]
    
    hv.var<-2*crossprod(A)%*%cov.tau%*%crossprod(A)
    h0.theory.var<-c()
    for (i in 2:length(anno.total)) {
      C<-anno.total[i]
      h0.theory.var[i]<-1/(C^2*(M-C)^2)*(C^2*hv.var[1,1]+M^2*hv.var[i,i]-2*M*C*hv.var[1,i])
    }
    P.theory<-pnorm(h0.test/sqrt(h0.theory.var), lower.tail = F)
    
    result<-list(tau, est.h, e, constant, sqrt(constant.var), h0.test, P.theory)
    names(result)<-c('Taus', 'Partition H2', 'Enrichment', 'intercept', 'intercept SD', 'e.stat','P')
    return(result)
  }else{
    test.tau<-solve(left)%*%right
    
    tau<-as.vector(test.tau)/N
    
    est.h<-crossprod(A)%*%tau
    h.total<-est.h[1]
    e<-est.h/M.anno/h.total
    
    #hypothesis test
    M<-anno.total[1]
    h0.test<-(est.h/anno.total)-(h.total-est.h)/(M-anno.total)
    
    #variance in theory
    cov.tau<-solve(left)/N^2
    
    hv.var<-2*crossprod(A)%*%cov.tau%*%crossprod(A)
    h0.theory.var<-c()
    for (i in 2:length(anno.total)) {
      C<-anno.total[i]
      h0.theory.var[i]<-1/(C^2*(M-C)^2)*(C^2*hv.var[1,1]+M^2*hv.var[i,i]-2*M*C*hv.var[1,i])
    }
    P.theory<-pnorm(h0.test/sqrt(h0.theory.var), lower.tail = F)
    
    result<-list(tau, est.h, e, 1, 0, h0.test, P.theory)
    names(result)<-c('Taus', 'Partition H2', 'Enrichment', 'intercept', 'intercept SD', 'e.stat','P')
    return(result)
  }
}

oneout<-function(i,raw.result,left.total,right.total,A,M.anno,anno.total,N,intercept){
  left=left.total-raw.result[[i]][[1]]
  right=right.total-raw.result[[i]][[2]]
  #the gls.generator
  if(intercept==T){
    test.tau<-solve(left)%*%right
    
    tau.all<-as.vector(test.tau)/N
    tau<-tau.all[-1]
    constant=tau.all[1]+1
    
    est.h<-crossprod(A)%*%tau
    h.total<-est.h[1]
    e<-est.h/M.anno/h.total
    
    #hypothesis test
    #M<-anno.total[1]
    #h0.test<-(est.h/anno.total)-(h.total-est.h)/(M-anno.total)
    
    #variance in theory
    #cov.tau.all<-solve(left)/N^2
    #cov.tau<-cov.tau.all[-1,-1]
    #constant.var=N^2*cov.tau.all[1,1]
    
    #hv.var<-2*crossprod(A)%*%cov.tau%*%crossprod(A)
    #h0.theory.var<-c()
    #for (i in 2:length(anno.total)) {
    #  C<-anno.total[i]
    #  h0.theory.var[i]<-1/(C^2*(M-C)^2)*(C^2*hv.var[1,1]+M^2*hv.var[i,i]-2*M*C*hv.var[1,i])
    #}
    #P.theory<-1-pnorm(h0.test/sqrt(h0.theory.var))
    
    result<-list(tau, est.h, e, constant)
    names(result)<-c('Taus', 'Partition H2', 'Enrichment')
    return(result)
  }else{
    test.tau<-solve(left)%*%right
    
    tau<-as.vector(test.tau)/N
    
    est.h<-crossprod(A)%*%tau
    h.total<-est.h[1]
    e<-est.h/M.anno/h.total
    
    #hypothesis test
    #M<-anno.total[1]
    #h0.test<-(est.h/anno.total)-(h.total-est.h)/(M-anno.total)
    
    #variance in theory
    #cov.tau<-solve(left)/N^2
    
    #hv.var<-2*crossprod(A)%*%cov.tau%*%crossprod(A)
    #h0.theory.var<-c()
    #for (i in 2:length(anno.total)) {
    #  C<-anno.total[i]
    #  h0.theory.var[i]<-1/(C^2*(M-C)^2)*(C^2*hv.var[1,1]+M^2*hv.var[i,i]-2*M*C*hv.var[1,i])
    #}
    #P.theory<-1-pnorm(h0.test/sqrt(h0.theory.var))
    
    result<-list(tau, est.h, e, 1)
    names(result)<-c('Taus', 'Partition H2', 'Enrichment', 'intercept')
    return(result)
  }
}



refpannel.par<-
  function(batch,LDinfo=LDinfo,LD.path=LD.path,chrs=chrs[chr],
           anno=anno,maf.snp=maf.snp,gwas.snp=gwas.snp,
           used.a=used.a){
    batch.name<-LDinfo[batch,'name']
    
    #merge all snps
    LDinfo.temp<-h5read(file = paste(LD.path,chrs,sep = '/'),name = LDinfo[batch,'name'])
    snp1<-LDinfo.temp$snplist
    
    #snp1.intersect<-Reduce(intersect, list(snp1,anno[,'SNP'],maf.snp,gwas.snp,ldsc.snp))
    snp1.intersect<-Reduce(intersect, list(snp1,anno[,'SNP'],maf.snp,gwas.snp))
    
    batch.used=0
    if(length(snp1.intersect)>0){
      Anno.temp<-anno[match(snp1.intersect,anno[,'SNP']),used.a]
      
      Anno.temp<-as.matrix(Anno.temp)
      
      LD.temp<-LDinfo.temp$ldblk
      LD.temp<-LD.temp[match(snp1.intersect,snp1),match(snp1.intersect,snp1)]
      
      LDSM.temp<-LDSMsum2(LD.temp,Anno.temp)
      x.ldsc.temp<-LDSCmatrix(LD.temp,Anno.temp)[,-1]
      batch.used=1
    }else{snp1.intersect<-LD.temp<-LDSM.temp<-x.ldsc.temp<-batch.used<-0}
    return(list(batch.used,snp1.intersect,LD.temp,LDSM.temp,
                x.ldsc.temp))
  }

refpannel.par.man<-
  function(batch,LDinfo,LD.path,chrs,
           anno,maf.snp,gwas.snp,
           used.a,betavar){
    batch.name<-LDinfo[batch,'name']
    
    #merge all snps
    LDinfo.temp<-h5read(file = paste(LD.path,chrs,sep = '/'),name = LDinfo[batch,'name'])
    
    snp1<-LDinfo.temp$snplist
    
    #snp1.intersect<-Reduce(intersect, list(snp1,anno[,'SNP'],maf.snp,gwas.snp,ldsc.snp))
    snp1.intersect<-Reduce(intersect, list(snp1,anno[,'SNP'],maf.snp,gwas.snp))
    
    batch.used=0
    if(length(snp1.intersect)>0){
      Anno.temp<-anno[match(snp1.intersect,anno[,'SNP']),used.a]
      
      Anno.temp<-as.matrix(Anno.temp)
      
      LD.temp<-LDinfo.temp$ldblk
      LD.temp<-LD.temp[match(snp1.intersect,snp1),match(snp1.intersect,snp1)]
      
      betavar<-as.vector(betavar)
      betavar2=betavar[match(snp1.intersect,snp1)]
      LDSM.temp<-LDSMman(LD.temp,betavar=betavar2)
      x.ldsc.temp<-LDSCmatrix(LD.temp,Anno.temp)[,-1]
      batch.used=1
    }else{snp1.intersect<-LD.temp<-LDSM.temp<-x.ldsc.temp<-batch.used<-0}
    return(list(batch.used,snp1.intersect,LD.temp,LDSM.temp,
                x.ldsc.temp))
  }
