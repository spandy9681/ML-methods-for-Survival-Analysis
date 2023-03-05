rm(list = ls())
set.seed(1)

# Required Libraries

library(ISLR2)
library(survival)
library(randomForestSRC)
library(survcomp)
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
cox_only = NULL
rsf_only = NULL
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
alpha1 = 1
s_sam = round((ncol(datas)-2)*alpha1)
# No of sampled observations
alpha2 = (1-exp(-1))
n_sam = round(alpha_trn*alpha2*nrow(datas))

# The LOOP
n_repeat_models = 30
for(sim in 1:n_repeat_models)
{
  # define the training and testing datasets
  train_ind <- sample(1:nrow(datas),size = alpha_trn*(nrow(datas)))
  datas_trn <- datas[train_ind,]
  datas_tst <- datas[-train_ind,]
  
  ###############
  ####OnlyCox####
  ###############
  
  # Now we fit the CoxPh model
  mod.cox <- coxph(formula = Surv(time, status) ~ ., data = datas_trn)
  
  # And calculate the concordance index
  risk_trn_cox <- predict(object = mod.cox,newdata = datas_trn,type = "risk")
  risk_tst_cox <- predict(object = mod.cox,newdata = datas_tst,type = "risk")
  cox_only[sim] = 1-get.cindex(datas_tst$time,datas_tst$status,(risk_tst_cox))
  
  ### Now comes the test data manipulation part
  datas_trn1 <- data.frame(datas_trn,"risk"= (risk_trn_cox))
  datas_tst1 <- data.frame(datas_tst,"risk"= (risk_tst_cox))
  
  ###############
  ####ONLYRSF####
  ###############
  
  # training the rsf model
  ntree <- b
  mod.rsf <- rfsrc(formula = Surv(time, status) ~ ., data = datas_trn,
                   ntree = ntree,nsplit = 0,splitrule = "logrank",importance = TRUE,nodesize = 18)
  risk_tst_rsf <- predict(object = mod.rsf,newdata = datas_tst,type = "risk")
  
  # And calculate the concordance index
  rsf_only[sim] <- 1-get.cindex(datas_tst$time,datas_tst$status,(risk_tst_rsf$predicted))
  
  ###############
  ####COXphRSF###
  ###############
  
  # Now we make the final predictions using DeepSurv model
  mod.rsf1 <- rfsrc(formula = Surv(time, status) ~ ., data = datas_trn1,ntree = ntree,
                    nsplit = 0,splitrule = "logrank",importance = TRUE,nodesize = 18)
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
    data_trn_sampled = as.data.frame(model.matrix(~.-1,data = data_trn_sampled))
    data_tst_sampled = as.data.frame(model.matrix(~.-1,data = data_tst_sampled))
    ## Hybrid Model Fitting ##
    mod.cox <- coxph(formula = Surv(time, status) ~ ., data = data_trn_sampled)
    risk_trn_cox <- predict(object = mod.cox,newdata = data_trn_sampled,type = "risk")
    datas_trn1 <- data.frame(data_trn_sampled,"risk"= (risk_trn_cox))
    mod.rsf1 <- rfsrc(formula = Surv(time, status) ~ ., data = datas_trn1,ntree = 1,
                      nsplit = 0,splitrule = "logrank",importance = TRUE,mtry = ncol(datas_trn1)-2,
                      bootstrap = "none",samptype = "swor",nodesize = 10)
    ## nodesize obtained after tuning
    risk_tst_cox <- predict(object = mod.cox,newdata = data_tst_sampled,type = "risk")
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

summary(rsf_only)
summary(cox_only)
summary(cox_rsf)
summary(cox_forest_boosted)

boxplot(rsf_only,cox_only,cox_rsf,cox_forest_boosted,names = c('rsf_only','cox_only','cox_rsf','cox_forest_boosted'))
abline(h = median(cox_forest_boosted),lty = 2)
summary_results <- cbind(DeepSurv_only,rsf_only,cox_only,cox_rsf,cox_deepsurv,cox_forest_boosted)
colnames(summary_results) <- c('DeepSurv_only','rsf_only','cox_only','cox_rsf','cox_deepsurv','cox_forest_boosted')

# Just an experiment
mod.cox <- coxph(formula = Surv(time, status) ~ .-1, data = datas_trn)
surv.obj <- survfit(formula = mod.cox,newdata = datas_trn)
surv.obj$time
ind <- sort(unique(datas_trn$time))
