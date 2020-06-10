
library(mvtnorm)
#Algorithm 1
k_split <- function(df, k) {
  folds <- split(sample(nrow(df), nrow(df), replace=F), as.factor(1:k))
  lapply(folds, function(idxs) df[idxs, ])
}
fun1=function(x,y) (x*y)
ff=function(x){
  aa=outer(x,x,"*")
  sum(aa[upper.tri(aa)])
}
fun2=function(idx,x) mean(x[idx])
fun3=function(x,KK){
  LL=de_mean(x,KK)
  term1=ff(x)
  temp2=x*LL$me2;term2=sum(temp2[upper.tri(temp2)])
  temp3=x*LL$me1;term3=sum(temp3[upper.tri(temp3)])
  temp4=LL$me1*LL$me2;term4=sum(temp4[upper.tri(temp4)])
  term1-term2-term3+term4
  
}

generate.df=function(model,type,p,n1=50,n2=50){
  mu1=rep(0,p)
  if (model==1){
    sig1=diag(p);sig2=2*diag(p)
    if (type==1){
      x1=rmvnorm(n1,mu1,sig1);x2=rmvnorm(n1,mu1,sig2)
      x3=rmvnorm(n2,mu1,sig2);x4=rmvnorm(n2,mu1,sig2)
      x=(rbind(x1,x2,x3,x4))
    }
    if (type==2){
      df=25
      x1=rmvt(n1,sig1,df);x2=rmvt(n1,sig2,df)
      x3=rmvt(n2,sig2,df);x4=rmvt(n2,sig2,df)
      x=(rbind(x1,x2,x3,x4))
    }
  }
  ###model 2
  if (model==2){
    sig2=diag(p)
    sig2[row(sig2)-1 == col(sig2)|row(sig2)+1 == col(sig2)] <- 0.5
    sig1=diag(p)
    sig1[row(sig1)-1 == col(sig1)|row(sig1)+1 == col(sig1)] <- -0.5
    if (type==1){
      x1=rmvnorm(n1,mu1,sig1);x2=rmvnorm(n1,mu1,sig2)
      x3=rmvnorm(n2,mu1,sig2);x4=rmvnorm(n2,mu1,sig2)
      x=(rbind(x1,x2,x3,x4))
    }
    if (type==2){
      df=30
      x1=rmvt(n1,sig1,df);x2=rmvt(n1,sig2,df)
      x3=rmvt(n2,sig2,df);x4=rmvt(n2,sig2,df)
      x=(rbind(x1,x2,x3,x4))
    }
  }
  
  ###model 3
  
  if (model==3){
    
    sig.m301=0.5^abs(outer(c(1:p),c(1:p),"-"))
    sig.m302=(-0.5)^abs(outer(c(1:p),c(1:p),"-"))
    sig1=sig.m301
    sig2=sig.m302
    if (type==1){
      x1=rmvnorm(n1,mu1,sig1);x2=rmvnorm(n1,mu1,sig1)
      x3=rmvnorm(n2,mu1,sig2);x4=rmvnorm(n2,mu1,sig2)
      x=(rbind(x1,x2,x3,x4))
    }
    if (type==2){
      df=25
      x1=rmvt(n1,sig1,df);x2=rmvt(n1,sig1,df)
      x3=rmvt(n2,sig2,df);x4=rmvt(n2,sig2,df)
      x=(rbind(x1,x2,x3,x4))
    }
  }
  ####model 4
  if (model==4){
    n1=n2=60;n3=80 
    sig.m301=0.8^abs(outer(c(1:p),c(1:p),"-"))
    sig.m302=(-0.8)^abs(outer(c(1:p),c(1:p),"-"))
    sig1=sig.m301
    sig2=sig.m302
    if (type==1){
      x1=rmvnorm(n1,mu1,sig1);x2=rmvnorm(n2,mu1,sig2)
      x3=rmvnorm(n3,mu1,diag(p))
      x=(rbind(x1,x2,x3))
    }
    if (type==2){
      df=25
      x1=rmvt(60,sig1,df);x2=rmvt(60,sig2,df)
      x3=rmvt(n3,diag(p),df)
      x=(rbind(x1,x2,x3))
    }
  }
  return(x)
}


Bchosen=function(n,m=as.integer(0.9*n)){
  ll=c(1:n)
  m1=as.integer(m/2);m2=m-m1
  ss1=array(0,dim=c(n,n,m1))
  ss2=array(0,dim=c(n,n,m2))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      set.seed(i+j)
      temp=sample(ll[-c(i,j)],m)
      temp1=sample(temp,m1)
      temp2=setdiff(temp,temp1)
      ss1[i,j,]=temp1;ss1[j,i,]=temp1
      ss2[i,j,]=temp2;ss2[j,i,]=temp2
      
    }
  }
  return(list(B1=ss1,B2=ss2))
}
##
##KK=Bchosen(n)
de_mean=function(x,KK){
  #x is 1dim
  n=length(x)
  me1=me2=matrix(0,n,n)
  B1=KK$B1;B2=KK$B2
  me1=apply(B1,c(1,2),fun2,x=x);diag(me1)=0
  me2=apply(B2,c(1,2),fun2,x=x);diag(me2)=0
  return(list(me1=me1,me2=me2))
}




AAA1=function(x,KK){
  #x is a n*p data matrix
  p=ncol(x);n=nrow(x)
  x=t(t(x)-matrix(apply(x,2,mean),p,n))
  ind=matrix(0,p,3)
  for (i in 1:(p)){
    newx=x*matrix(x[,i],n,p)
    us=apply(newx,2,fun3,KK=KK)
    ind[i,1]=which(us==us[rank(us)==p])
    ind[i,2]=which(us==us[rank(us)==(p-1)])
    ind[i,3]=which(us==us[rank(us)==(p-2)]);
    print(i)
  }
  ind
}




####conduct a test

p=50;n=200;type=1;model1=1;P=diag(n)-rep(1,n)%*%t(rep(1,n))/n
KK=Bchosen(n)
model2=2
###model 2
km1=gmm1=skm1=rfm1=prop1=NULL
ad=0;mu1=mu2=c(0,rep(0,p-1));M1=M2=NULL
pos1=c(rep(1,50),rep(2,50),rep(2,50),rep(2,50))
pos2=c(rep(2,50),rep(1,50),rep(1,50),rep(1,50)) 


  x=generate.df(model2,type,p)
  
  ind=AAA1(x,KK)
  cr1=cr2=rep(0,p)
  cr1[1]=cr2[1]=ifelse((2)%in%ind[1,],1,0)
  cr1[p]=cr2[p]=ifelse((p-1)%in%ind[p,],1,0)
  for (i in 2:(p-1)){
    cr1[i]=ifelse((i+1)%in%ind[i,1]|(i-1)%in%ind[i,1],1,0)
    cr2[i]=ifelse((i+1)%in%ind[i,1:2]|(i-1)%in%ind[i,1:2],1,0)
  }

mean(cr1)
 mean(cr2)

         ###clustering
         newx1=x*matrix(x[,ind[,1]],n,p)+x*matrix(x[,ind[,2]],n,p)
       #S=P%*%newx1%*%t(newx1)%*%P/n
       S=newx1%*%t(newx1)/n
       temp05=kmeans(svd(S)$u[,1],2)$cluster
       prop1_1=temp05-pos1
       prop1_2=temp05-pos2
      1-max(length(which(prop1_1==0)), length(which(prop1_2==0)))/n
       
  



