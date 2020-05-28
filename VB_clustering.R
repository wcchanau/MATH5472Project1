rm(list=ls())
gc()
while(!require('rstudioapi')) install.packages('rstudioapi')
while(!require('xtable')) install.packages('xtable')
while(!require('ggplot2')) install.packages('ggplot2')
while(!require('mvtnorm')) install.packages('mvtnorm')
while(!require('dplyr')) install.packages('dplyr')
setwd(dirname(getSourceEditorContext()$path))
set.seed(1)

n=5000
K=5
sigma2=4
mu=t(rmvnorm(K,sigma = diag(sigma2,2)))
c=sample(K,n,replace = TRUE,prob=rep(1/K,K))
y=apply(mu[,c],2,function(mean) rmvnorm(1,mean))

parameter=list(m=y[,sample.int(n,K)],s2=array(diag(0.5,2),c(2,2,K)),phi=matrix(1/K,n,K),L=c(0))
repeat{
  parameter$phi=exp(t(y)%*%parameter$m-matrix(apply(parameter$s2,3,sum)+diag(t(parameter$m)%*%parameter$m),n,K,byrow = TRUE)/2)%>%
    apply(1,function(x) x/sum(x))%>%t
  parameter$m=t(t(y%*%parameter$phi)/(1/sigma2+colSums(parameter$phi)))
  parameter$s2=sapply(1/(1/sigma2+colSums(parameter$phi)),diag,nrow=2,simplify = 'array')
  parameter$L=c(parameter$L,
                -sum(colSums(apply(parameter$s2,3,diag)-parameter$m^2)/(2*sigma2))+ # \sum_{k=1}^K E(log p(\mu_k);m_k,s_k)
                  -sum(parameter$phi)*log(K) + # \sum_{i=1}^n E(log p(c_i);\phi_i)
                  -sum(sapply(1:K,function(x) ((y-parameter$m[,x])^2)%*%parameter$phi[,x]))/2 - # \sum_{i=1}^n E(log p(x_i|c_i,\bm{\mu});\phi_i,\bm{m},\bm{s})
                  sum(parameter$phi*log(parameter$phi))- # \sum_{i=1}^n E(log q(c_i;\phi_i))
                  -sum(apply(parameter$s2,3,function(x) log(det(x))))/2 # \sum_{k=1}^K E(log q(\mu_k;m_k,s_k))
                )
  if(abs(diff(tail(parameter$L,2)))<(10^-8)) break
}

# reorder
new_order=apply((dist(t(cbind(mu,parameter$m)))%>%as.matrix)[1:K,-1:-K],1,which.min)
parameter$phi=parameter$phi[,new_order]
parameter$m=parameter$m[,new_order]
parameter$s2=parameter$s2[,,new_order]

# Testing
n_testing=15000
c_testing=sample(K,n_testing,replace = TRUE,prob=rep(1/K,K))
y_testing=apply(mu[,c_testing],2,function(mean) rmvnorm(1,mean))
phi_testing=exp(t(y_testing)%*%parameter$m-matrix(apply(parameter$s2,3,sum)+diag(t(parameter$m)%*%parameter$m),n_testing,K,byrow = TRUE)/2)%>%
  apply(1,function(x) x/sum(x))%>%t

# Error
1-sum(apply(parameter$phi,1,which.max)==c)/n
1-sum(apply(phi_testing,1,which.max)==c_testing)/n_testing

# Table
table(c,estimate=apply(parameter$phi,1,which.max))%>%xtable
table(c_testing,estimate=apply(phi_testing,1,which.max))%>%xtable

# Plot
rbind(
  data.frame(x=c,y=t(y),TE="Actual",training="Training"),
  data.frame(x=apply(parameter$phi,1,which.max),y=t(y),TE="Predict",training="Training"),
  data.frame(x=c_testing,y=t(y_testing),TE="Actual",training="Testing"),
  data.frame(x=apply(phi_testing,1,which.max),y=t(y_testing),TE="Predict",training="Testing"))%>%
  ggplot(aes(x=y.1,y=y.2,color=as.factor(x)))+geom_point()+
  facet_grid(TE~training)+
  geom_point(data = data.frame(y=t(mu)), colour="black", size = 4, shape = 13)+  
  labs(color="Cluster") +
  xlab("")+
  ylab("")


ggplot(data.frame(x=c,y=t(y)),aes(x=y.1,y=y.2,color=as.factor(x)))+
  geom_point()+
  labs(color="Cluster") +
  ggtitle("True clustering")+
  xlab("")+
  ylab("")+
  geom_point(data = data.frame(y=t(mu)), colour="black", size = 4, shape = 13)

ggplot(data.frame(x=apply(parameter$phi,1,which.max),y=t(y)),aes(x=y.1,y=y.2,color=as.factor(x)))+
  geom_point()+
  labs(color="Cluster") +
  ggtitle("Estimated clustering")+
  xlab("")+
  ylab("")+
  geom_point(data = data.frame(y=t(parameter$m)), colour="black", size = 4, shape = 13)

ggplot(data.frame(x=c_testing,y=t(y_testing)),aes(x=y.1,y=y.2,color=as.factor(x)))+
  geom_point()+
  labs(color="Cluster") +
  ggtitle("True clustering (testing)")+
  xlab("")+
  ylab("")+
  geom_point(data = data.frame(y=t(mu)), colour="black", size = 4, shape = 13)

ggplot(data.frame(x=apply(phi_testing,1,which.max),y=t(y_testing)),aes(x=y.1,y=y.2,color=as.factor(x)))+
  geom_point()+
  labs(color="Cluster") +
  ggtitle("Estimated clustering (testing)")+
  xlab("")+
  ylab("")+
  geom_point(data = data.frame(y=t(parameter$m)), colour="black", size = 4, shape = 13)

ggplot(data.frame(x=2:length(parameter$L)-1,y=parameter$L[-1]),aes(x=x,y=y))+
  geom_line() + 
  xlab("Iterate")+
  ylab("ELBO")
