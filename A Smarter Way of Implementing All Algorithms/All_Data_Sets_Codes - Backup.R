## Important Packages ##
library(ISLR2)
library(survival)
library(randomForestSRC)

## BrainCancer Dataset ##
datas <- na.omit(BrainCancer)
datas$sex <- as.factor(datas$sex)
datas$diagnosis <- as.factor(datas$diagnosis)
datas$loc <- as.factor(datas$loc)
datas$stereo <- as.factor(datas$stereo)
datas <- datas[,c(1:6,8,7)]

## Specifications alpha_trn = 0.6, alpha1 = 0.5

## GBSG Dataset ##
datas <- na.omit(gbsg)
names(datas)[10] <- "time"
datas$meno <- as.factor(datas$meno)
datas$grade <- as.factor(datas$grade)
datas$hormon <- as.factor(datas$hormon)

## Veteran Dataset ##
datas <- na.omit(veteran)
datas$trt <- as.factor(datas$trt)
datas$celltype <- as.factor(datas$celltype)
datas$prior <- as.factor(datas$prior)
datas <- datas[,c(1,2,5,6,7,8,3,4)]

## LUNG Dataset ##
datas <- na.omit(lung)
datas$sex <- as.factor(datas$sex)
datas$ph.ecog <- as.factor(datas$ph.ecog)
datas <- datas[,c(1,4:10,2,3)]
datas$status <- datas$status-1

## Rats Dataset ##
datas <- na.omit(rats)
datas$sex <- as.factor(datas$sex)
datas <- datas[,c(1,2,5,3,4)]

## Simulated Linear Dataset ##
n = 500;p = 5
x = rnorm(n*p,0,1)
x = matrix(x,nrow = n)
c0 = 5
beta = runif(p,0.5,1.5)
times = rexp(n,c0*(exp(x%*%beta)))  
time.censor = rexp(n,c0*(exp(x%*%beta))) 
summary(time.censor)
delta <- ifelse(times<time.censor, 1, 0)
time <- ifelse(times<time.censor, times, time.censor)
datas <- data.frame("x" = x,time = time,status = delta)

## Simulated Non-Linear Dataset 1 ##
n = 500;p = 5
x = rnorm(n*p,0,1)
x = matrix(x,nrow = n)
c0 = 5
beta = runif(p,0.5,1.5)
times = rexp(n,c0*((x%*%beta)^2))  
time.censor = rexp(n,c0*((x%*%beta)^2)) 
summary(time.censor)
delta <- ifelse(times<time.censor, 1, 0)
time <- ifelse(times<time.censor, times, time.censor)
datas <- data.frame("x" = x,time = time,status = delta)

## Simulated Non-Linear Dataset 2 ##
n = 500;p = 5
x = rnorm(n*p,0,1)
x = matrix(x,nrow = n)
c0 = 5
beta = runif(p,0.5,1.5)
times = rexp(n,c0*(abs(x%*%beta)))  
time.censor = rexp(n,c0*(abs(x%*%beta))) 
summary(time.censor)
delta <- ifelse(times<time.censor, 1, 0)
time <- ifelse(times<time.censor, times, time.censor)
datas <- data.frame("x" = x,time = time,status = delta)
