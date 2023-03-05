rm(list = ls())
## Important Packages required for datasets ##
library(ISLR2)
library(survival)
library(randomForestSRC)

## Creating a big list where all the datas are stored ##
data_list = list()

## BrainCancer Dataset ##
datas <- na.omit(BrainCancer)
datas$sex <- as.factor(datas$sex)
datas$diagnosis <- as.factor(datas$diagnosis)
datas$loc <- as.factor(datas$loc)
datas$stereo <- as.factor(datas$stereo)
datas <- datas[,c(1:6,8,7)]

# First component
data_list[[1]] <- datas

## GBSG Dataset ##
datas <- na.omit(gbsg)
names(datas)[10] <- "time"
datas$meno <- as.factor(datas$meno)
datas$grade <- as.factor(datas$grade)
datas$hormon <- as.factor(datas$hormon)

# Second component
data_list[[2]] <- datas

## Veteran Dataset ##
datas <- na.omit(veteran)
datas$trt <- as.factor(datas$trt)
datas$celltype <- as.factor(datas$celltype)
datas$prior <- as.factor(datas$prior)
datas <- datas[,c(1,2,5,6,7,8,3,4)]

# Third component
data_list[[3]] <- datas

## LUNG Dataset ##
datas <- na.omit(lung)
datas$sex <- as.factor(datas$sex)
datas$ph.ecog <- as.factor(datas$ph.ecog)
datas <- datas[,c(1,4:10,2,3)]
datas$status <- datas$status-1

# Fourth component
data_list[[4]] <- datas

## Rats Dataset ##
datas <- na.omit(rats)
datas$sex <- as.factor(datas$sex)
datas <- datas[,c(1,2,5,3,4)]

# Fifth component
data_list[[5]] <- datas

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

# Sixth component
data_list[[6]] <- datas

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

# Seventh component
data_list[[7]] <- datas


### Running all the datasets at once

# Required Libraries

library(ISLR2)
library(survival)
library(randomForestSRC)
library(survcomp)
library(reticulate)
library(tensorflow)
library(keras)
library(survivalmodels)
library(svMisc)

# Setting seed values for numpy and torch

np <- reticulate::import("numpy")
np$random$seed(1L)
torch <- reticulate::import("torch")
torch$manual_seed(1L)

## Running the bigggg loop
for(data_loop in 1:7)
{
  datas <- as.data.frame(data_list[[data_loop]])
  
  # Now running the big loop
  DeepSurv_only = NULL
  cox_only = NULL
  rsf_only = NULL
  cox_deepsurv = NULL
  cox_rsf = NULL
  cox_forest_boosted = NULL
  set.seed(108)
  # Components of ensemble form of 
  # hybrid cox+survival tree model
  ### Bagged Cox+Survival Tree Model
  ind_time_status <- NULL
  ind_time_status[1] <- which(names(datas) == "time")
  ind_time_status[2] <- which(names(datas) == "status")
  ind_cov <- (1:ncol(datas))[-ind_time_status]
  # No of models
  b = 50
  # Train test splitting ratio
  alpha_trn = 0.7
  # No of sampled covariates
  alpha1 = 0.5
  s_sam = round((ncol(datas)-2)*alpha1)
  # No of sampled observations
  alpha2 = (1-exp(-1))
  n_sam = round(alpha_trn*alpha2*nrow(datas))
  nodesz = round(sqrt(ncol(datas))) + 1
  # The LOOP
  n_repeat_models = 30
  for(sim in 1:n_repeat_models)
  {
    ###############
    #ONLY DEEPSURV#
    ###############
    
    # define the training and testing datasets
    train_ind <- sample(1:nrow(datas),size = alpha_trn*(nrow(datas)))
    datas_trn <- datas[train_ind,]
    datas_tst <- datas[-train_ind,]
    
    # fitting the DeepSurv model
    torch$manual_seed(1L)
    mod.deepsurv <- deepsurv(formula = Surv(time, status) ~ ., data = datas_trn,num_nodes = round(sqrt(nrow(datas_trn))),epochs = 10)
    ypred.deepsurv <- predict(object = mod.deepsurv, newdata = datas_tst,type = "risk")
    
    # And calculate the concordance index
    DeepSurv_only[sim] = 1-get.cindex(datas_tst$time,datas_tst$status,ypred.deepsurv)
    
    ###############
    ####OnlyCox####
    ###############
    
    # Now we fit the CoxPh model
    mod.cox <- coxph(formula = Surv(time, status) ~ ., data = datas_trn)
    risk_tst_cox <- predict(object = mod.cox,newdata = datas_tst)
    
    # And calculate the concordance index
    cox_only[sim] = 1-get.cindex(datas_tst$time,datas_tst$status,(risk_tst_cox))
    
    ###############
    ####ONLYRSF####
    ###############
    
    # training the rsf model
    ntree <- b
    mod.rsf <- rfsrc(formula = Surv(time, status) ~ ., data = datas_trn,
                     ntree = ntree,nsplit = 0,splitrule = "logrank",importance = TRUE,nodesize = nodesz)
    risk_tst_rsf <- predict(object = mod.rsf,newdata = datas_tst)
    
    # And calculate the concordance index
    rsf_only[sim] <- 1-get.cindex(datas_tst$time,datas_tst$status,(risk_tst_rsf$predicted))
    
    ###############
    ##CoxDeepSurv##
    ###############
    
    # Now we modify the train dataset for deepsurv model fitting
    risk_trn_cox <- predict(object = mod.cox,newdata = datas_trn)
    datas_trn1 <- data.frame(datas_trn,"risk"= (risk_trn_cox))
    torch$manual_seed(1L)
    mod.deepsurv1 <- deepsurv(formula = Surv(time, status) ~ ., data = datas_trn1,num_nodes = round(sqrt(nrow(datas_trn1))),epochs = 10)
    
    ### Now comes the test data manipulation part
    datas_tst1 <- data.frame(datas_tst,"risk"= (risk_tst_cox))
    
    # Now we make the final predictions using DeepSurv model
    ypred.deepsurv1 <- predict(object = mod.deepsurv1, newdata = datas_tst1,type = "risk")
    cox_deepsurv[sim] = 1-get.cindex(datas_tst$time,datas_tst$status,(ypred.deepsurv1))
    
    ###############
    ####COXphRSF###
    ###############
    
    # Now we make the final predictions using DeepSurv model
    mod.rsf1 <- rfsrc(formula = Surv(time, status) ~ ., data = datas_trn1,ntree = ntree,
                      nsplit = 0,splitrule = "logrank",importance = TRUE,nodesize = nodesz)
    ypred.rsf1 <- predict(object = mod.rsf1, newdata = datas_tst1,type = "risk")
    cox_rsf[sim] <- 1-get.cindex(datas_tst1$time,datas_tst1$status,(ypred.rsf1$predicted))
    
    ###############
    ##Boosted Cox##
    #####Tree######
    ###############
    
    y_pred_mat = matrix(nrow = b,ncol = nrow(datas_tst))
    for(i in 1:b)
    {
      ## Data sub-selection ##
      obs_sampled = sample(1:nrow(datas_trn),size = n_sam,replace = FALSE)
      cov_sampled = sample(1:(ncol(datas_trn)-2),size = s_sam,replace = FALSE)
      data_trn_sampled = datas_trn[obs_sampled,c(cov_sampled,ind_time_status)]
      data_tst_sampled = datas_tst[,c(cov_sampled,ind_time_status)]
      ## Defining Dummy Variables ##
      ## for factor variates ##
      data_trn_sampled = as.data.frame(model.matrix(~.,data = data_trn_sampled))
      data_tst_sampled = as.data.frame(model.matrix(~.,data = data_tst_sampled))
      ## Hybrid Model Fitting ##
      mod.cox <- coxph(formula = Surv(time, status) ~ ., data = data_trn_sampled)
      risk_trn_cox <- predict(object = mod.cox,newdata = data_trn_sampled)
      datas_trn1 <- data.frame(data_trn_sampled,"risk"= (risk_trn_cox))
      mod.rsf1 <- rfsrc(formula = Surv(time, status) ~ ., data = datas_trn1,ntree = 1,
                        nsplit = 0,splitrule = "logrank",importance = TRUE,mtry = ncol(datas_trn1)-2,
                        bootstrap = "none",samptype = "swor",nodesize = nodesz)
      ## nodesize obtained after tuning
      risk_tst_cox <- predict(object = mod.cox,newdata = data_tst_sampled)
      datas_tst1 <- data.frame(data_tst_sampled,"risk"= (risk_tst_cox))
      ## Putting the fitted value in a matrix
      ypred.rsf1 <- predict(object = mod.rsf1, newdata = datas_tst1,type = "risk")
      y_pred_mat[i,] <- ypred.rsf1$predicted
    }
    y_pred_mean_final <- apply(y_pred_mat, 2, mean)
    ## C-index
    cox_forest_boosted[sim] = 1-get.cindex(datas_tst$time,datas_tst$status,y_pred_mean_final)
    
    # To get idea of the computational load print the sim value
    progress(value = sim,max.value = n_repeat_models)
  }
  
  ## Writing down the results
  summary_results <- cbind(DeepSurv_only,rsf_only,cox_only,cox_rsf,cox_deepsurv,cox_forest_boosted)
  colnames(summary_results) <- c('DeepSurv_only','rsf_only','cox_only','cox_rsf','cox_deepsurv','cox_forest_boosted')
  summary_all_methods <- as.data.frame(apply(summary_results,2,summary))
  setwd(dir = "E:\\Dekstop\\Tanujit Sir Project\\Loop_Results_All_datasets")
  write.csv(x = summary_all_methods,file = paste("data",data_loop,".csv"))
  
  print(paste(data_loop,"th dataset complete"))
}
