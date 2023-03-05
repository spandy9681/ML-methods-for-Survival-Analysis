rm(list = ls())
set.seed(1)
np <- reticulate::import("numpy")
np$random$seed(1L)
torch <- reticulate::import("torch")
torch$manual_seed(1L)

# Required Libraries

library(ISLR2)
library(survival)
library(randomForestSRC)
library(survcomp)
library(tensorflow)
library(keras)
library(survivalmodels)
library(svMisc)

# Any Survival data with censoring index and survival times coded as "status" & "time" 
# respectively put your dataset here || ||  
#                                 || || ||
#                                 \/ \/ \/

## BrainCancer Dataset ##
datas <- na.omit(gbsg)
names(datas)[10] <- "time"
datas$meno <- as.factor(datas$meno)
datas$grade <- as.factor(datas$grade)
datas$hormon <- as.factor(datas$hormon)

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
                   ntree = ntree,nsplit = 0,splitrule = "logrank",importance = TRUE,nodesize = 1)
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
                    nsplit = 0,splitrule = "logrank",importance = TRUE,nodesize = 1)
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
                      bootstrap = "none",samptype = "swor",nodesize = 1)
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

summary(DeepSurv_only)
summary(rsf_only)
summary(cox_only)
summary(cox_rsf)
summary(cox_deepsurv)
summary(cox_forest_boosted)
boxplot(DeepSurv_only,rsf_only,cox_only,cox_rsf,cox_deepsurv,cox_forest_boosted,names = c('DeepSurv_only','rsf_only','cox_only','cox_rsf','cox_deepsurv','cox_forest_boosted'))
abline(h = median(cox_forest_boosted),lty = 2)
summary_results <- cbind(DeepSurv_only,rsf_only,cox_only,cox_rsf,cox_deepsurv,cox_forest_boosted)
colnames(summary_results) <- c('DeepSurv_only','rsf_only','cox_only','cox_rsf','cox_deepsurv','cox_forest_boosted')
summary_all_methods <- as.data.frame(apply(summary_results,2,summary))
setwd(dir = "E:\\Dekstop\\Tanujit Sir Project\\Loop_Results_All_datasets")
write.csv(x = summary_all_methods,file = paste("data",data_loop,".csv"))

data_loop = 1
range_results <- apply(summary_results, 2, range)
range_results[2,]-range_results[1,]
plot(cox_only,type = "o")
lines(cox_rsf,type = "o",col = "red",lty = 2)
lines(cox_forest_boosted,type = "o",col = "blue",lty = 3)
sum(cox_forest_boosted > cox_rsf & cox_forest_boosted > cox_only)
sum(cox_only > cox_forest_boosted & cox_only > cox_rsf)
sum(cox_rsf > cox_only & cox_rsf > cox_forest_boosted)
