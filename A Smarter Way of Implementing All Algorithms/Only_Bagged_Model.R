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
library(mgcv)

## Simulated Non-Linear Dataset ##
n = 500;p = 6
x = rnorm(n*p,0,1)
x = matrix(x,nrow = n)
c0 = 5
beta = runif(p,0.5,1.5)
h = exp(-(x^2%*%beta))
times = rexp(n,c0*(exp(h)))  
time.censor = rexp(n,c0*(exp(h))) 
summary(time.censor)
delta <- ifelse(times<time.censor, 1, 0)
time <- ifelse(times<time.censor, times, time.censor)
datas <- data.frame("x" = x,time = time,status = delta)

# Now running the big loop
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
alpha_trn = 0.6
# No of sampled covariates
alpha1 = 0.5
s_sam = round((ncol(datas)-2)*alpha1)
# No of sampled observations
alpha2 = (1-exp(-1))
n_sam = round(alpha_trn*alpha2*nrow(datas))

# The LOOP
n_repeat_models = 1
for(sim in 1:n_repeat_models)
{
  # define the training and testing datasets
  train_ind <- sample(1:nrow(datas),size = alpha_trn*(nrow(datas)))
  datas_trn <- datas[train_ind,]
  datas_tst <- datas[-train_ind,]
  
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
    ind_time_status_sampled = NULL
    ind_time_status_sampled[1] = which(names(data_trn_sampled) == "time")
    ind_time_status_sampled[2] = which(names(data_trn_sampled) == "status")
    ## Hybrid Model Fitting ##
    formula = paste("s(",names(data_trn_sampled)[-ind_time_status_sampled],")")
    formula = paste(formula,collapse = "+")
    formula = paste("time ~",formula,collapse = "")
    formula = as.formula(formula)
    m4 <- mgcv::gam(formula = formula, 
                    data = data_trn_sampled, family = "cox.ph", weights = status)
    pred.gam.trn <- predict(object = m4,newdata = data_trn_sampled)
    datas_trn1 <- data.frame(data_trn_sampled,"pred"= (pred.gam.trn))
    mod.rsf1 <- rfsrc(formula = Surv(time, status) ~ ., data = datas_trn1,ntree = 100,
                      nsplit = 0,splitrule = "logrank",importance = TRUE,mtry = ncol(datas_trn1)-2,
                      bootstrap = "none",samptype = "swor")
    ## nodesize obtained after tuning
    pred.gam.tst <- predict(object = m4,newdata = data_tst_sampled)
    datas_tst1 <- data.frame(data_tst_sampled,"pred"= (pred.gam.tst))
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
summary(cox_forest_boosted)
boxplot(cox_forest_boosted,names = c('cox_forest_boosted'))
abline(h = median(cox_forest_boosted),lty = 2)

# Using GAM model for Cox regression
formula = paste("s(",names(datas)[-ind_time_status],")")
formula = paste(formula,collapse = "+")
formula = paste("time ~",formula,collapse = "")
formula = as.formula(formula)
m4 <- mgcv::gam(formula = formula, 
                data = datas_trn, family = "cox.ph", weights = status)
summary(m4)
pred.gam.trn <- predict(object = m4,newdata = datas_trn)
datas_trn1 <- data.frame(datas_trn,"pred"= (pred.gam.trn))
mod.rsf1 <- rfsrc(formula = Surv(time, status) ~ ., data = datas_trn1,ntree = 100,
                  nsplit = 0,splitrule = "logrank",importance = TRUE,mtry = ncol(datas_trn1)-2,
                  bootstrap = "none",samptype = "swor")
## nodesize obtained after tuning
pred.gam.tst <- predict(object = m4,newdata = datas_tst)
datas_tst1 <- data.frame(datas_tst,"pred"= (pred.gam.tst))
## Putting the fitted value in a matrix
ypred.rsf1 <- predict(object = mod.rsf1, newdata = datas_tst1,type = "risk")
1-get.cindex(datas_tst$time,datas_tst$status,ypred.rsf1$predicted)
1-get.cindex(datas_tst$time,datas_tst$status,pred.gam.tst)

## Using survival probs
surv.cox <- survfit(mod.cox,newdata = datas_trn)
dim(surv.cox$surv)
surv.cox$time
plot(surv.cox)
quants <- quantile(surv.cox$time,probs = c(0.25,0.5,0.75))
ind.times <- NULL
ind.times[1] <- which.min(abs(surv.cox$time - quants[1]))
ind.times[2] <- which.min(abs(surv.cox$time - quants[2]))
ind.times[3] <- which.min(abs(surv.cox$time - quants[3]))

length(surv.cox$time)
dim(surv.cox$surv)
dim(datas_trn)

dim(surv.cox$surv[ind.times,])
