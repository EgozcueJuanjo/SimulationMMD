#############################################
# Function simulating count-samples from a
# Multinomial-Multinomial Distribution (MMD)
# input parameters:
# nsam number of samples to be simulated
# M number of trials in the first sampling
# p true probabilities for D features D=length(p)
# N number of trials in the second sampling
# optran option: 0, N fixed (default)
#                1, N random Poisson
#                2, N random triangular
# pardist: parameters for distribution of N
#          (mean-mode, min, max)
# output:
#  matrix of MMD simulated counts
#       rows: features
#       columns: cases, individuals or patients
# uses function rmultinom in R-package "nnet"
#............................................
# Written by JJ Egozcue, April 2020
#############################################
rMMD <- function(nsam,M,N,p,optran=0,pardist=NULL){
  library(nnet)
  library(EnvStats)
  #  first multinomial sampling
  count1 <- rmultinom(n=nsam,size=M,prob=p)
  D <- length(p)
  # second multinomial sampling
  q <- count1/outer(rep(1,D),colSums(count1))
  # non random N
  if(optran==0){count2 <- apply(q,2, "rmultinom",n=1,size=N)}
  # random N (Poisson)
  if(optran!=0){
    # Poisson mean equal N or given in pardist[1]
    if(optran==1){
      lam <-ifelse(is.null(pardist),N,pardist[1])
      Nvec <- rpois(n=nsam,lambda=lam)}
    # triangular pardist=(mode,min,max)
    if(optran==2){Nvec <- floor(
      rtri(n=nsam,min=pardist[2],max=pardist[3],mode=pardist[1]))}
    # simulation with random N
    extq <- rbind(q,Nvec)
    count2 <- apply(extq,2,function(qN){qpr <- qN[1:length(qN)-1];
    siz<-qN[length(qN)];rmultinom(n=1,size=siz,prob=qpr)})
  }
  return(count2)
}
############################################

############################################
# function simulating the distribution MAN
# Multinomial-Asymptotic Normal
# input:
#  nsam number of samples to be computed
#  M  number of trials in the multinomial sampling
#  p  multinomial probabilities in the first
#     multinomial sampling. Its length define the
#     number of features considered D=length(p)
# output:
#  simp  (nsam,D)-matrix containing the nsam simulated
#        samples expressed in proportions, including
#        zeros produced in the first sampling.
# uses R-libraries "nnet" and "MAS"
#      auxiliary function "buildCM"
#####################################################
# written by JJ Egozcue, May 2020
#####################################################
rMAN <- function(nsam,M,p){
  library(nnet)
  library(MASS)
  #  first multinomial sampling
  # individuals by rows, features by columns
  D <- length(p)
  count1 <- t(rmultinom(n=nsam,size=M,prob=p))
  # second sampling: asymptotic normal distribution
  # probability parameters
  q <- count1/outer(rowSums(count1),rep(1,D))
  # loop for samples (as many as nsam)
  simp <- matrix(0,nrow=nsam,ncol=D)
  for(isim in 1:nsam){
    fz <- (q[isim,]!=0)   # position of non zero
    d <- sum(fz==TRUE)  # features in second sampling
    d1 <- d-1
    if(d1<=0){print(c("d-1 error",d))}
    # build contast matrix (depends on number zeros)
    V <- buildCM(d)
    # ilr mean and variance
    ilrqmean <- t(V) %*% as.matrix(log(q[isim,fz==TRUE]))
    Sig <- t(V) %*% diag(1/q[isim,fz==TRUE])%*% V
    # normal simulation
    ilrsim <- mvrnorm(n=1,mu=ilrqmean,Sigma=Sig)
    # back to proportions, ilr^(-1)
    simpnc <- exp( matrix(ilrsim,nrow=1) %*% t(V))
    clssimp <- sum(simpnc)
    simpd <- simpnc/clssimp
    simp[isim,fz==TRUE] <- simpd
  }
  return(simp)
}
############################################
# buildCM builds a (d,d-1)-contrast matrix
# for ilr coordinates
############################################
buildCM <- function(d){
  V <- matrix(0,nrow=d,ncol=(d-1))
  for(j in 1:(d-1)){
    vdiag = sqrt((d-j)/(d-j+1))
    vout = -sqrt((d-j)/(d-j+1)) /(d-j)
    V[j,j]=vdiag
    V[(j+1):d,j]=vout
  }
  return(V)
}
############################################



# Example: use functions rMMD, rMAN

pseqnc <- seq(from=1,to=1000,length=103); spseq<- sum(pseqnc)
pseq <- pseqnc/spseq
simMMD <- rMMD(nsam=180,M=100,N=5500,p=pseq,optran=2,
                 pardist=c(5500,800,20000) )  
simMAN <- rMAN(nsam=180,M=100,p=pseq)
# zero patterns (use R-package zCompositions)
library(zCompositions)
zPatterns(X=t(simMMD),label=0)
zPatterns(X=simMAN,label=0)
