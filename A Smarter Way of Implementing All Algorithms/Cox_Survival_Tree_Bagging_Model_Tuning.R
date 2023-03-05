rm(list = ls())

# Required Libraries

library(ISLR2)
library(survival)
library(randomForestSRC)

# Any Survival data with censoring index and survival times coded as "status" & "time" 
# respectively put your dataset here || ||  
#                                 || ||
#                                 \/ \/

## BrainCancer Dataset ##
datas <- na.omit(BrainCancer)
datas$sex <- as.factor(datas$sex)
datas$diagnosis <- as.factor(datas$diagnosis)
datas$loc <- as.factor(datas$loc)
datas$stereo <- as.factor(datas$stereo)
datas <- datas[,c(1:6,8,7)]


# Running the cox model as a beginning
cox_only1 = NULL
set.seed(108)

for(sim in 1:20)
{
  # define the training and testing datasets
  train_ind <- sample(1:nrow(datas),size = 0.7*(nrow(datas)))
  datas_trn <- datas[train_ind,]
  datas_tst <- datas[-train_ind,]
  
  ###############
  ####OnlyCox####
  ###############
  
  # Now we fit the CoxPh model
  mod.cox <- coxph(formula = Surv(time, status) ~ ., data = datas_trn)
  risk_tst_cox <- predict(object = mod.cox,newdata = datas_tst)
  
  # And calculate the concordance index
  cox_only1[sim] = 1-get.cindex(datas_tst$time,datas_tst$status,(risk_tst_cox))
}

summary(cox_only1)

# Also the RSF model as another model in consideration
rsf_only1 = NULL
set.seed(108)

for(sim in 1:20)
{
  # define the training and testing datasets
  train_ind <- sample(1:nrow(datas),size = 0.7*(nrow(datas)))
  datas_trn <- datas[train_ind,]
  datas_tst <- datas[-train_ind,]
  
  ###############
  ####ONLYRSF####
  ###############
  
  # training the rsf model
  ntree <- round(sqrt(nrow(datas_trn)))
  mod.rsf <- rfsrc(formula = Surv(time, status) ~ ., data = datas_trn,ntree = ntree,nsplit = 0,splitrule = "logrank",importance = TRUE,nodesize = 3)
  risk_tst_rsf <- predict(object = mod.rsf,newdata = datas_tst)
  
  # And calculate the concordance index
  rsf_only1[sim] <- 1-get.cindex(datas_tst$time,datas_tst$status,(risk_tst_rsf$predicted))
}

summary(rsf_only1)

## Now the bagged cox model
ind_time_status <- NULL
ind_time_status[1] <- which(names(datas) == "time")
ind_time_status[2] <- which(names(datas) == "status")
ind_cov <- (1:ncol(datas))[-ind_time_status]

# No of models
b = 50
# No of sampled covariates
alpha1 = 0.5
s_sam = round((ncol(datas)-2)*alpha1)
# No of sampled observations
alpha2 = 1/exp(1)
n_sam = round(nrow(datas_trn)*alpha2)
y_pred_mat = matrix(nrow = b,ncol = nrow(datas_tst))

cox_boosted = NULL
set.seed(108)

for(sim in 1:20)
{
  # define the training and testing datasets
  train_ind <- sample(1:nrow(datas),size = 0.7*(nrow(datas)))
  datas_trn <- datas[train_ind,]
  datas_tst <- datas[-train_ind,]
  
  for(i in 1:b)
  {
    ## Data sub-selection ##
    obs_sampled = sample(1:nrow(datas_trn),size = nrow(datas_trn),replace = TRUE)
    obs_sampled = unique(obs_sampled)
    cov_sampled = sample(1:(ncol(datas_trn)-2),size = s_sam,replace = FALSE)
    data_trn_sampled = datas_trn[obs_sampled,c(cov_sampled,ind_time_status)]
    data_tst_sampled = datas_tst[,c(cov_sampled,ind_time_status)]
    ## Cox Model Fitting ##
    mod.cox <- coxph(formula = Surv(time, status) ~ ., data = data_trn_sampled)
    risk_tst_cox <- predict(object = mod.cox,newdata = data_tst_sampled,type = "risk")
    y_pred_mat[i,] <- risk_tst_cox
  }
  y_pred_final <- apply(y_pred_mat, 2, mean)
  ## C-index
  cox_boosted[sim] = 1-get.cindex(datas_tst$time,datas_tst$status,y_pred_final)
  ## Printing progress
  print(noquote(paste(sim,"/",20,"itrations done with value obtained:",round(cox_boosted[sim],3))))
}

summary(cox_boosted)

### Bagged Cox+Survival Tree Model
ind_time_status <- NULL
ind_time_status[1] <- which(names(datas) == "time")
ind_time_status[2] <- which(names(datas) == "status")
ind_cov <- (1:ncol(datas))[-ind_time_status]

# No of models
b = 50
# No of sampled covariates
alpha1 = 0.5
s_sam = round((ncol(datas)-2)*alpha1)
# No of sampled observations
alpha2 = 1/exp(1)
n_sam = round(nrow(datas_trn)*alpha2)
y_pred_mat = matrix(nrow = b,ncol = nrow(datas_tst))

cox_forest_boosted = NULL
set.seed(108)

for(sim in 1:20)
{
  # define the training and testing datasets
  train_ind <- sample(1:nrow(datas),size = 0.7*(nrow(datas)))
  datas_trn <- datas[train_ind,]
  datas_tst <- datas[-train_ind,]
  
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
                      bootstrap = "none",samptype = "swor",nodesize = 12)
    risk_tst_cox <- predict(object = mod.cox,newdata = data_tst_sampled)
    datas_tst1 <- data.frame(data_tst_sampled,"risk"= (risk_tst_cox))
    ## Putting the fitted value in a matrix
    ypred.rsf1 <- predict(object = mod.rsf1, newdata = datas_tst1,type = "risk")
    y_pred_mat[i,] <- ypred.rsf1$predicted
  }
  y_pred_final <- apply(y_pred_mat, 2, mean)
  ## C-index
  cox_forest_boosted[sim] = 1-get.cindex(datas_tst$time,datas_tst$status,y_pred_final)
  ## Printing progress
  print(noquote(paste(sim,"/",20,"itrations done with value obtained:"
                      ,round(cox_forest_boosted[sim],2))))
}

summary(cox_forest_boosted)

###########################
## Hyperparameter Tuning ##
###########################

ind_time_status <- NULL
ind_time_status[1] <- which(names(datas) == "time")
ind_time_status[2] <- which(names(datas) == "status")
ind_cov <- (1:ncol(datas))[-ind_time_status]

# No of models
Performance_Mat = matrix(nrow = 10,ncol = 10)
b = 30

for(alpha1 in 1:9)
{
  for(j in 1:10)
  {
    alpha2 = 1/exp(1)
    # No of sampled covariates
    s_sam = alpha1
    #s_sam = round((ncol(datas)-2)*alpha1/10)
    # No of sampled observations
    n_sam = round(nrow(datas_trn)*alpha2)
    y_pred_mat = matrix(nrow = b,ncol = nrow(datas_tst))
    set.seed(108)
    
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
                        bootstrap = "none",samptype = "swor",nodesize = 15) # 15 is calculated as 
      # the opt node size..
      risk_tst_cox <- predict(object = mod.cox,newdata = data_tst_sampled)
      datas_tst1 <- data.frame(data_tst_sampled,"risk"= (risk_tst_cox))
      ## Putting the fitted value in a matrix
      ypred.rsf1 <- predict(object = mod.rsf1, newdata = datas_tst1,type = "risk")
      y_pred_mat[i,] <- ypred.rsf1$predicted
    }
    y_pred_final <- apply(y_pred_mat, 2, mean)
    ## C-index
    Performance_Mat[alpha1,j] = 1-get.cindex(datas_tst$time,datas_tst$status,y_pred_final)
  }
  print(noquote(paste(alpha1*j,"/",100,"iterations done.")))
}

Performance_Mat <- na.omit(Performance_Mat)
Performance_Med_Vec <- apply(Performance_Mat,MARGIN = 1,median)
which(Performance_Med_Vec == max(Performance_Med_Vec))

# Further tuning other parameters

ind_time_status <- NULL
ind_time_status[1] <- which(names(datas) == "time")
ind_time_status[2] <- which(names(datas) == "status")
ind_cov <- (1:ncol(datas))[-ind_time_status]

# No of models
Performance_Vec = NULL
b = 30
# After tuning we get
alpha1 = 1;alpha2 = 1/exp(1)

# No of sampled covariates
s_sam = round((ncol(datas)-2)*alpha1)
# No of sampled observations
n_sam = round(0.7*(nrow(datas))*alpha2)
y_pred_mat = matrix(nrow = b,ncol = nrow(datas)-round(0.7*(nrow(datas))))
set.seed(108)
ind_size = 20
Performance_Mat = matrix(nrow = 20,ncol = 1)

for(ind2 in 1:5)
{
  # define the training and testing datasets
  train_ind <- sample(1:nrow(datas),size = 0.7*(nrow(datas)))
  datas_trn <- datas[train_ind,]
  datas_tst <- datas[-train_ind,]
  
  for(ind1 in 1:ind_size)
  {
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
      nd_size <- round((n_sam/2)*ind1/20)
      mod.rsf1 <- rfsrc(formula = Surv(time, status) ~ ., data = datas_trn1,ntree = 1,
                        nsplit = 0,splitrule = "logrank",importance = TRUE,mtry = ncol(datas_trn1)-2,
                        bootstrap = "none",samptype = "swor",nodesize = nd_size)
      risk_tst_cox <- predict(object = mod.cox,newdata = data_tst_sampled)
      datas_tst1 <- data.frame(data_tst_sampled,"risk"= (risk_tst_cox))
      ## Putting the fitted value in a matrix
      ypred.rsf1 <- predict(object = mod.rsf1, newdata = datas_tst1,type = "risk")
      y_pred_mat[i,] <- ypred.rsf1$predicted
    }
    y_pred_final <- apply(y_pred_mat, 2, mean)
    ## C-index
    Performance_Vec[ind1] = 1-get.cindex(datas_tst$time,datas_tst$status,y_pred_final)
    if(ind1 %% 10 == 0){print(noquote(paste(ind1,"/",ind_size,"iterations done.")))}
  }
  Performance_Mat = cbind(Performance_Mat,Performance_Vec)
}

Performance_Med <- apply(Performance_Mat[,-1], 1, median)
nodesize_vec <- round((n_sam/2)*(1:ind_size)/20)
C_Val_Tuning <- data.frame(nodesize_vec,Performance_Mat[,-1],Performance_Med)

# Optimal choice of Nodesize
nodesize_vec[which(Performance_Med == sort(Performance_Med,decreasing = TRUE)[1])]
nodesize_vec[which(Performance_Med == sort(Performance_Med,decreasing = TRUE)[2])]
nodesize_vec[which(Performance_Med == sort(Performance_Med,decreasing = TRUE)[3])]

### After doing all the tuning now we implement them in the final model
### for this stage
### Tuned parameters : alpha1 = 1;nodesize = 32

### Bagged Cox+Survival Tree Model
ind_time_status <- NULL
ind_time_status[1] <- which(names(datas) == "time")
ind_time_status[2] <- which(names(datas) == "status")
ind_cov <- (1:ncol(datas))[-ind_time_status]

# No of models
b = 30
# No of sampled covariates
alpha1 = 1
s_sam = round((ncol(datas)-2)*alpha1)
# No of sampled observations
alpha2 = 1/exp(1)
n_sam = round(nrow(datas_trn)*alpha2)
y_pred_mat = matrix(nrow = b,ncol = nrow(datas_tst))

cox_forest_boosted = NULL
set.seed(108)

for(sim in 1:20)
{
  # define the training and testing datasets
  train_ind <- sample(1:nrow(datas),size = 0.7*(nrow(datas)))
  datas_trn <- datas[train_ind,]
  datas_tst <- datas[-train_ind,]
  
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
                      bootstrap = "none",samptype = "swor",nodesize = 32)
    ## nodesize obtained after tuning
    risk_tst_cox <- predict(object = mod.cox,newdata = data_tst_sampled)
    datas_tst1 <- data.frame(data_tst_sampled,"risk"= (risk_tst_cox))
    ## Putting the fitted value in a matrix
    ypred.rsf1 <- predict(object = mod.rsf1, newdata = datas_tst1,type = "risk")
    y_pred_mat[i,] <- ypred.rsf1$predicted
  }
  y_pred_final <- apply(y_pred_mat, 2, mean)
  ## C-index
  cox_forest_boosted[sim] = 1-get.cindex(datas_tst$time,datas_tst$status,y_pred_final)
  ## Printing progress
  print(noquote(paste(sim,"/",20,"itrations done with value obtained:"
                      ,round(cox_forest_boosted[sim],2))))
}

summary(cox_forest_boosted)


### Boosted Cox+Random Forest #
### with different weights ####

ind_time_status <- NULL
ind_time_status[1] <- which(names(datas) == "time")
ind_time_status[2] <- which(names(datas) == "status")
ind_cov <- (1:ncol(datas))[-ind_time_status]

# No of models
b = 25
# No of sampled covariates
alpha1 = 0.5
s_sam = round((ncol(datas)-2)*alpha1)
# No of sampled observations
alpha2 = 0.9
n_sam = round(nrow(datas_trn)*alpha2)
y_pred_mat = matrix(nrow = b,ncol = nrow(datas_tst))

cox_forest_boosted = NULL
set.seed(108)

for(sim in 1:20)
{
  # define the training and testing datasets
  train_ind <- sample(1:nrow(datas),size = 0.7*(nrow(datas)))
  datas_trn <- datas[train_ind,]
  datas_tst <- datas[-train_ind,]
  w <- NULL
  for(i in 1:b)
  {
    ## Data sub-selection ##
    obs_sampled = sample(1:nrow(datas_trn),size = n_sam,replace = FALSE)
    cov_sampled = sample(1:(ncol(datas_trn)-2),size = s_sam,replace = FALSE)
    data_trn_sampled = datas_trn[obs_sampled,c(cov_sampled,ind_time_status)]
    data_tst_sampled = datas_tst[,c(cov_sampled,ind_time_status)]
    ## Weight Vector
    w[i] <- sum(data_trn_sampled$status == 1)/n_sam
    ## Hybrid Model Fitting ##
    mod.cox <- coxph(formula = Surv(time, status) ~ ., data = data_trn_sampled)
    risk_trn_cox <- predict(object = mod.cox,newdata = data_trn_sampled)
    datas_trn1 <- data.frame(data_trn_sampled,"risk"= (risk_trn_cox))
    mod.rsf1 <- rfsrc(formula = Surv(time, status) ~ ., data = datas_trn1,ntree = 1,
                      nsplit = 0,splitrule = "logrank",importance = TRUE,mtry = ncol(datas_trn1)-2,
                      bootstrap = "none",samptype = "swor",nodesize = 15)
    risk_tst_cox <- predict(object = mod.cox,newdata = data_tst_sampled)
    datas_tst1 <- data.frame(data_tst_sampled,"risk"= (risk_tst_cox))
    ## Putting the fitted value in a matrix
    ypred.rsf1 <- predict(object = mod.rsf1, newdata = datas_tst1,type = "risk")
    y_pred_mat[i,] <- ypred.rsf1$predicted
  }
  y_pred_final_weighted <- t(y_pred_mat)%*%w/sum(w)
  ## C-index
  cox_forest_boosted[sim] = 1-get.cindex(datas_tst$time,datas_tst$status,y_pred_final_weighted)
  ## Printing progress
  print(noquote(paste(sim,"/",20,"itrations done with value obtained:"
                      ,round(cox_forest_boosted[sim],2))))
}

summary(cox_forest_boosted)

### Boosted Random Survival Tree + Cox Model with different weights
## Now the boosting model
ind_time_status <- NULL
ind_time_status[1] <- which(names(datas) == "time")
ind_time_status[2] <- which(names(datas) == "status")
ind_cov <- (1:ncol(datas))[-ind_time_status]

# No of models
b = 50
# No of sampled covariates
alpha1 = 0.5
s_sam = round((ncol(datas)-2)*alpha1)
# No of sampled observations
alpha2 = 1/exp(1)
n_sam = round(nrow(datas_trn)*alpha2)
y_pred_mat = matrix(nrow = b,ncol = nrow(datas_tst))

cox_forest_boosted = NULL
set.seed(108)

for(sim in 1:20)
{
  # define the training and testing datasets
  train_ind <- sample(1:nrow(datas),size = 0.7*(nrow(datas)))
  datas_trn <- datas[train_ind,]
  datas_tst <- datas[-train_ind,]
  w <- NULL
  for(i in 1:b)
  {
    ## Data sub-selection ##
    obs_sampled = sample(1:nrow(datas_trn),size = n_sam,replace = FALSE)
    cov_sampled = sample(1:(ncol(datas_trn)-2),size = s_sam,replace = FALSE)
    data_trn_sampled = datas_trn[obs_sampled,c(cov_sampled,ind_time_status)]
    data_tst_sampled = datas_tst[,c(cov_sampled,ind_time_status)]
    ## Weight Vector
    w[i] <- sum(data_trn_sampled$status == 1)/n_sam
    ## Hybrid Model Fitting ##
    mod.rsf <- rfsrc(formula = Surv(time, status) ~ ., data = data_trn_sampled,ntree = 1,
                     nsplit = 0,splitrule = "logrank",importance = TRUE,mtry = ncol(datas_trn1)-3,
                     bootstrap = "none",samptype = "swor",nodesize = 12)
    risk_trn_rsf <- predict(object = mod.rsf,newdata = data_trn_sampled)
    datas_trn1 <- data.frame(data_trn_sampled,"risk"= (risk_trn_rsf$predicted))
    mod.cox1 <- coxph(formula = Surv(time, status) ~ ., data = datas_trn1)
    risk_tst_rsf <- predict(object = mod.rsf,newdata = data_tst_sampled)
    datas_tst1 <- data.frame(data_tst_sampled,"risk"= (risk_tst_rsf$predicted))
    ## Putting the fitted value in a matrix
    ypred.cox1 <- predict(object = mod.cox1, newdata = datas_tst1,type = "risk")
    y_pred_mat[i,] <- ypred.cox1
  }
  y_pred_final_weighted <- t(y_pred_mat)%*%w/sum(w)
  ## C-index
  cox_forest_boosted[sim] = 1-get.cindex(datas_tst$time,datas_tst$status,y_pred_final_weighted)
  ## Printing progress
  print(noquote(paste(sim,"/",20,"itrations done with value obtained:"
                      ,round(cox_forest_boosted[sim],2))))
}

summary(cox_forest_boosted)

### A different method for bagging with stratified sampling
### of censored and observed samples
### Boosted Cox+Random Forest
## Now the boosting model

ind_time_status <- NULL
ind_time_status[1] <- which(names(datas) == "time")
ind_time_status[2] <- which(names(datas) == "status")
ind_cov <- (1:ncol(datas))[-ind_time_status]

# No of models
b = 50
# No of sampled covariates
alpha1 = 1
s_sam = round((ncol(datas)-2)*alpha1)
# No of sampled observations
alpha2 = 1/exp(1)
n_sam = round(nrow(datas_trn)*alpha2)
y_pred_mat = matrix(nrow = b,ncol = nrow(datas_tst))

cox_forest_boosted = NULL
set.seed(108)

for(sim in 1:20)
{
  # define the training and testing datasets
  train_ind <- sample(1:nrow(datas),size = 0.7*(nrow(datas)))
  datas_trn <- datas[train_ind,]
  datas_tst <- datas[-train_ind,]
  ind_cen <- which(datas_trn$status == 0)
  ind_obs <- which(datas_trn$status == 1)
  cen_ratio <- sum(datas_trn$status)/nrow(datas_trn)
  for(i in 1:b)
  {
    ## Data sub-selection ##
    obs_sampled1 = sample(ind_cen,size = n_sam*(1-cen_ratio),replace = FALSE)
    obs_sampled2 = sample(ind_obs,size = n_sam*cen_ratio,replace = FALSE)
    obs_sampled = c(obs_sampled1,obs_sampled2)
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
                      bootstrap = "none",samptype = "swor",nodesize = 23)
    ## nodesize obtained after tuning
    risk_tst_cox <- predict(object = mod.cox,newdata = data_tst_sampled)
    datas_tst1 <- data.frame(data_tst_sampled,"risk"= (risk_tst_cox))
    ## Putting the fitted value in a matrix
    ypred.rsf1 <- predict(object = mod.rsf1, newdata = datas_tst1,type = "risk")
    y_pred_mat[i,] <- ypred.rsf1$predicted
  }
  y_pred_final <- apply(y_pred_mat, 2, mean)
  ## C-index
  cox_forest_boosted[sim] = 1-get.cindex(datas_tst$time,datas_tst$status,y_pred_final)
  ## Printing progress
  print(noquote(paste(sim,"/",20,"itrations done with value obtained:"
                      ,round(cox_forest_boosted[sim],2))))
}

summary(cox_forest_boosted)


### Tuning for this different type of bagging model with 
### stratified sampling techinque
###########################
## Hyperparameter Tuning ##
###########################

ind_time_status <- NULL
ind_time_status[1] <- which(names(datas) == "time")
ind_time_status[2] <- which(names(datas) == "status")
ind_cov <- (1:ncol(datas))[-ind_time_status]

# No of models
Performance_Mat = matrix(nrow = 10,ncol = 20)
b = 30
alpha2 = 1/exp(1)
# No of sampled covariates
s_sam = round((ncol(datas)-2)*alpha1/10)
# No of sampled observations
n_sam = round(nrow(datas_trn)*alpha2)
y_pred_mat = matrix(nrow = b,ncol = nrow(datas_tst))
set.seed(108)

for(alpha1 in 1:10)
{
  # No of sampled covariates
  s_sam = round((ncol(datas)-2)*alpha1/10)
  
  for(j in 1:20)
  {
    train_ind <- sample(1:nrow(datas),size = 0.7*(nrow(datas)))
    datas_trn <- datas[train_ind,]
    datas_tst <- datas[-train_ind,]
    ind_cen <- which(datas_trn$status == 0)
    ind_obs <- which(datas_trn$status == 1)
    cen_ratio <- sum(datas_trn$status)/nrow(datas_trn)
    for(i in 1:b)
    {
      ## Data sub-selection ##
      obs_sampled1 = sample(ind_cen,size = n_sam*(1-cen_ratio),replace = FALSE)
      obs_sampled2 = sample(ind_obs,size = n_sam*cen_ratio,replace = FALSE)
      obs_sampled = c(obs_sampled1,obs_sampled2)
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
                        bootstrap = "none",samptype = "swor",nodesize = 19)
      risk_tst_cox <- predict(object = mod.cox,newdata = data_tst_sampled)
      datas_tst1 <- data.frame(data_tst_sampled,"risk"= (risk_tst_cox))
      ## Putting the fitted value in a matrix
      ypred.rsf1 <- predict(object = mod.rsf1, newdata = datas_tst1,type = "risk")
      y_pred_mat[i,] <- ypred.rsf1$predicted
    }
    y_pred_final <- apply(y_pred_mat, 2, mean)
    ## C-index
    Performance_Mat[alpha1,j] = 1-get.cindex(datas_tst$time,datas_tst$status,y_pred_final)
  }
  print(noquote(paste(alpha1*j,"/",200,"iterations done.")))
}

Performance_Mat <- na.omit(Performance_Mat)
Performance_Med_Vec <- apply(Performance_Mat,MARGIN = 1,median)
alpha1_grid <- round((ncol(datas)-2)*1:10/10)
alpha1_grid[which(Performance_Med_Vec == max(Performance_Med_Vec))]

# Further tuning other parameters

ind_time_status <- NULL
ind_time_status[1] <- which(names(datas) == "time")
ind_time_status[2] <- which(names(datas) == "status")
ind_cov <- (1:ncol(datas))[-ind_time_status]

# No of models
Performance_Vec = NULL
b = 30
# After tuning we get
alpha1 = 1;alpha2 = 1/exp(1)

# No of sampled covariates
s_sam = round((ncol(datas)-2)*alpha1)
# No of sampled observations
n_sam = round(0.7*(nrow(datas))*alpha2)
y_pred_mat = matrix(nrow = b,ncol = nrow(datas)-round(0.7*(nrow(datas))))
set.seed(108)
ind_size = 20
Performance_Mat = matrix(nrow = 20,ncol = 1)

for(ind2 in 1:5)
{
  # define the training and testing datasets
  train_ind <- sample(1:nrow(datas),size = 0.7*(nrow(datas)))
  datas_trn <- datas[train_ind,]
  datas_tst <- datas[-train_ind,]
  ind_cen <- which(datas_trn$status == 0)
  ind_obs <- which(datas_trn$status == 1)
  cen_ratio <- sum(datas_trn$status)/nrow(datas_trn)
  
  for(ind1 in 1:ind_size)
  {
    for(i in 1:b)
    {
      ## Data sub-selection ##
      obs_sampled1 = sample(ind_cen,size = n_sam*(1-cen_ratio),replace = FALSE)
      obs_sampled2 = sample(ind_obs,size = n_sam*cen_ratio,replace = FALSE)
      obs_sampled = c(obs_sampled1,obs_sampled2)
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
      nd_size <- round((n_sam/2)*ind1/20)
      mod.rsf1 <- rfsrc(formula = Surv(time, status) ~ ., data = datas_trn1,ntree = 1,
                        nsplit = 0,splitrule = "logrank",importance = TRUE,mtry = ncol(datas_trn1)-2,
                        bootstrap = "none",samptype = "swor",nodesize = nd_size)
      risk_tst_cox <- predict(object = mod.cox,newdata = data_tst_sampled)
      datas_tst1 <- data.frame(data_tst_sampled,"risk"= (risk_tst_cox))
      ## Putting the fitted value in a matrix
      ypred.rsf1 <- predict(object = mod.rsf1, newdata = datas_tst1,type = "risk")
      y_pred_mat[i,] <- ypred.rsf1$predicted
    }
    y_pred_final <- apply(y_pred_mat, 2, mean)
    ## C-index
    Performance_Vec[ind1] = 1-get.cindex(datas_tst$time,datas_tst$status,y_pred_final)
    if(ind1 %% 10 == 0){print(noquote(paste(ind1,"/",ind_size,"iterations done.")))}
  }
  Performance_Mat = cbind(Performance_Mat,Performance_Vec)
}

Performance_Med <- apply(Performance_Mat[,-1], 1, median)
nodesize_vec <- round((n_sam/2)*(1:ind_size)/20)
C_Val_Tuning <- data.frame(nodesize_vec,Performance_Mat[,-1],Performance_Med)

# Optimal choice of Nodesize
nodesize_vec[which(Performance_Med == sort(Performance_Med,decreasing = TRUE)[1])]
nodesize_vec[which(Performance_Med == sort(Performance_Med,decreasing = TRUE)[2])]
nodesize_vec[which(Performance_Med == sort(Performance_Med,decreasing = TRUE)[3])]

### Tuning for different type of  
### prediction ensembling using ranking 
### for stratified sampling techinque
###########################
## Hyperparameter Tuning ##
###########################

ind_time_status <- NULL
ind_time_status[1] <- which(names(datas) == "time")
ind_time_status[2] <- which(names(datas) == "status")
ind_cov <- (1:ncol(datas))[-ind_time_status]

# No of models
Performance_Mat = matrix(nrow = 10,ncol = 20)
b = 50
# No of sampled covariates
alpha2 = 1/exp(1)
# No of sampled covariates
s_sam = round((ncol(datas)-2)*alpha1)
# No of sampled observations
n_sam = round(0.7*(nrow(datas))*alpha2)
y_pred_mat = matrix(nrow = b,ncol = nrow(datas)-round(0.7*(nrow(datas))))
set.seed(108)

for(alpha1 in 1:10)
{
  # No of sampled covariates
  s_sam = round((ncol(datas)-2)*alpha1/10)
  
  for(j in 1:20)
  {
    train_ind <- sample(1:nrow(datas),size = round(0.7*(nrow(datas))))
    datas_trn <- datas[train_ind,]
    datas_tst <- datas[-train_ind,]
    ind_cen <- which(datas_trn$status == 0)
    ind_obs <- which(datas_trn$status == 1)
    cen_ratio <- sum(datas_trn$status)/nrow(datas_trn)
    for(i in 1:b)
    {
      ## Data sub-selection ##
      obs_sampled1 = sample(ind_cen,size = n_sam*(1-cen_ratio),replace = FALSE)
      obs_sampled2 = sample(ind_obs,size = n_sam*cen_ratio,replace = FALSE)
      obs_sampled = c(obs_sampled1,obs_sampled2)
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
                        bootstrap = "none",samptype = "swor",nodesize = 19)
      risk_tst_cox <- predict(object = mod.cox,newdata = data_tst_sampled)
      datas_tst1 <- data.frame(data_tst_sampled,"risk"= (risk_tst_cox))
      ## Putting the fitted value in a matrix
      ypred.rsf1 <- predict(object = mod.rsf1, newdata = datas_tst1,type = "risk")
      y_pred_mat[i,] <- ypred.rsf1$predicted
    }
    y_pred_rank <- apply(y_pred_mat, 1, rank)
    y_pred_rank <- t(y_pred_rank)
    y_pred_rank_final <- apply(y_pred_rank, 2, mean)
    ## C-index
    Performance_Mat[alpha1,j] = 1-get.cindex(datas_tst$time,datas_tst$status,y_pred_rank_final)
  }
  print(noquote(paste(alpha1*j,"/",200,"iterations done.")))
}

Performance_Mat <- na.omit(Performance_Mat)
Performance_Med_Vec <- apply(Performance_Mat,MARGIN = 1,median)
alpha1_grid <- round((ncol(datas)-2)*1:10/10)
alpha1_grid[which(Performance_Med_Vec == max(Performance_Med_Vec))]

# Further tuning other parameters

ind_time_status <- NULL
ind_time_status[1] <- which(names(datas) == "time")
ind_time_status[2] <- which(names(datas) == "status")
ind_cov <- (1:ncol(datas))[-ind_time_status]

# No of models
Performance_Vec = NULL
b = 50
# After tuning we get
alpha1 = 1;alpha2 = 1/exp(1)

# No of sampled covariates
s_sam = round((ncol(datas)-2)*alpha1)
# No of sampled observations
n_sam = round(0.7*(nrow(datas))*alpha2)
y_pred_mat = matrix(nrow = b,ncol = nrow(datas)-round(0.7*(nrow(datas))))
set.seed(108)
ind_size = 20
Performance_Mat = matrix(nrow = 20,ncol = 1)

for(ind2 in 1:5)
{
  # define the training and testing datasets
  train_ind <- sample(1:nrow(datas),size = round(0.7*(nrow(datas))))
  datas_trn <- datas[train_ind,]
  datas_tst <- datas[-train_ind,]
  ind_cen <- which(datas_trn$status == 0)
  ind_obs <- which(datas_trn$status == 1)
  cen_ratio <- sum(datas_trn$status)/nrow(datas_trn)
  
  for(ind1 in 1:ind_size)
  {
    for(i in 1:b)
    {
      ## Data sub-selection ##
      obs_sampled1 = sample(ind_cen,size = n_sam*(1-cen_ratio),replace = FALSE)
      obs_sampled2 = sample(ind_obs,size = n_sam*cen_ratio,replace = FALSE)
      obs_sampled = c(obs_sampled1,obs_sampled2)
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
      nd_size <- round((n_sam/2)*ind1/20)
      mod.rsf1 <- rfsrc(formula = Surv(time, status) ~ ., data = datas_trn1,ntree = 1,
                        nsplit = 0,splitrule = "logrank",importance = TRUE,mtry = ncol(datas_trn1)-2,
                        bootstrap = "none",samptype = "swor",nodesize = nd_size)
      risk_tst_cox <- predict(object = mod.cox,newdata = data_tst_sampled)
      datas_tst1 <- data.frame(data_tst_sampled,"risk"= (risk_tst_cox))
      ## Putting the fitted value in a matrix
      ypred.rsf1 <- predict(object = mod.rsf1, newdata = datas_tst1,type = "risk")
      y_pred_mat[i,] <- ypred.rsf1$predicted
    }
    y_pred_rank <- apply(y_pred_mat, 1, rank)
    y_pred_rank <- t(y_pred_rank)
    y_pred_rank_final <- apply(y_pred_rank, 2, mean)
    ## C-index
    Performance_Vec[ind1] = 1-get.cindex(datas_tst$time,datas_tst$status,y_pred_rank_final)
    if(ind1 %% 10 == 0){print(noquote(paste(ind1,"/",ind_size,"iterations done.")))}
  }
  Performance_Mat = cbind(Performance_Mat,Performance_Vec)
}

Performance_Med <- apply(Performance_Mat[,-1], 1, median)
nodesize_vec <- round((n_sam/2)*(1:ind_size)/20)
C_Val_Tuning <- data.frame(nodesize_vec,Performance_Mat[,-1],Performance_Med)

# Optimal choice of Nodesize
nodesize_vec[which(Performance_Med == sort(Performance_Med,decreasing = TRUE)[1])]
nodesize_vec[which(Performance_Med == sort(Performance_Med,decreasing = TRUE)[2])]
nodesize_vec[which(Performance_Med == sort(Performance_Med,decreasing = TRUE)[3])]
nodesize_vec[which(Performance_Med == sort(Performance_Med,decreasing = TRUE)[4])]
nodesize_vec[which(Performance_Med == sort(Performance_Med,decreasing = TRUE)[5])]
nodesize_vec[which(Performance_Med == sort(Performance_Med,decreasing = TRUE)[6])]

### After doing all the tuning now we implement them in the final model
### for this stage
### Tuned parameters : alpha1 = 1;nodesize = 32

### Bagged Cox+Survival Tree Model with different risk ensembling
ind_time_status <- NULL
ind_time_status[1] <- which(names(datas) == "time")
ind_time_status[2] <- which(names(datas) == "status")
ind_cov <- (1:ncol(datas))[-ind_time_status]

# No of models
b = 50
# No of sampled covariates
alpha1 = 0.7
s_sam = round((ncol(datas)-2)*alpha1)
# No of sampled observations
alpha2 = 1/exp(1)
n_sam = round(0.7*nrow(datas)*alpha2)
y_pred_mat = matrix(nrow = b,ncol = nrow(datas_tst))

cox_forest_boosted = NULL
set.seed(108)

for(sim in 1:20)
{
  # define the training and testing datasets
  train_ind <- sample(1:nrow(datas),size = round(0.7*(nrow(datas))))
  datas_trn <- datas[train_ind,]
  datas_tst <- datas[-train_ind,]
  
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
                      bootstrap = "none",samptype = "swor",nodesize = 12)
    ## nodesize obtained after tuning
    risk_tst_cox <- predict(object = mod.cox,newdata = data_tst_sampled)
    datas_tst1 <- data.frame(data_tst_sampled,"risk"= (risk_tst_cox))
    ## Putting the fitted value in a matrix
    ypred.rsf1 <- predict(object = mod.rsf1, newdata = datas_tst1,type = "risk")
    y_pred_mat[i,] <- ypred.rsf1$predicted
  }
  y_pred_rank <- apply(y_pred_mat, 1, rank)
  y_pred_rank <- t(y_pred_rank)
  y_pred_rank_final <- apply(y_pred_rank, 2, mean)
  ## C-index
  cox_forest_boosted[sim] = 1-get.cindex(datas_tst$time,datas_tst$status,y_pred_rank_final)
  ## Printing progress
  print(noquote(paste(sim,"/",20,"itrations done with value obtained:"
                      ,round(cox_forest_boosted[sim],2))))
}

summary(cox_forest_boosted)
