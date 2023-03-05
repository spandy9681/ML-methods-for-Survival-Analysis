rm(list = ls())
set.seed(1)

# Packages
library(ISLR2)
library(survival)
library(randomForestSRC)
library(muhaz)

## Simulated Non-Linear Dataset ##
n = 500;p = 4
x = rnorm(n*p,0,1)
x = matrix(x,nrow = n)
c0 = 5
beta = runif(p,0.5,1.5)
h = (exp(-x^2%*%beta))
times = rexp(n,c0*(exp(-h)))  
time.censor = rexp(n,c0*(exp(h))) 
summary(time.censor)
delta <- ifelse(times<time.censor, 1, 0)
time <- ifelse(times<time.censor, times, time.censor)
datas <- data.frame("x" = x,time = time,status = delta)

## All the basic analysis one by one
Surv.Obj <- Surv(datas$time,datas$status)

# Kaplan-Meire Estimator
result.km <- survfit(Surv.Obj ~ 1, conf.type="log-log")
plot(result.km)

# Non-parametric hazard function estimation
result.simple <- muhaz(datas$time,datas$status,bw.method="global", b.cor="none")
plot(result.simple)
hist(datas$time,probability = TRUE)
pairs(datas)

# Cox Hazard Model
mod.cox <- coxph(Surv(time,status) ~ 1,data = datas)
summary(mod.cox)

# residuals
rr.0 <- residuals(mod.cox,type = "martingale")
par(mfrow=c(1,1))
plot(rr.0 ~ datas$x.1)
smoothSEcurve(rr.0, datas$x.1)    

# Comparison Between Cox/Cox-RSF/RSF
# define the training and testing datasets
alpha_trn = 0.7
set.seed(108)
train_ind <- sample(1:nrow(datas),size = alpha_trn*(nrow(datas)))
datas_trn <- datas[train_ind,]
datas_tst <- datas[-train_ind,]

###############
####OnlyCox####
###############

# Now we fit the CoxPh model
mod.cox <- coxph(formula = Surv(time, status) ~ ., data = datas_trn)
risk_tst_cox <- predict(object = mod.cox,newdata = datas_tst)

# And calculate the concordance index
1-get.cindex(datas_tst$time,datas_tst$status,(risk_tst_cox))

###############
####ONLYRSF####
###############

# training the rsf model
ntree <- 50
mod.rsf <- rfsrc(formula = Surv(time, status) ~ ., data = datas_trn,ntree = ntree,
                 nsplit = 0,splitrule = "logrank",importance = TRUE,nodesize = 1)
risk_tst_rsf <- predict(object = mod.rsf,newdata = datas_tst)

# And calculate the concordance index
1-get.cindex(datas_tst$time,datas_tst$status,(risk_tst_rsf$predicted))

###############
####COXphRSF###
###############

# Now we make the final predictions using DeepSurv model
surv.cox.trn <- survfit(mod.cox,newdata = datas_trn)
quants.trn <- quantile(surv.cox$time,probs = c(0.25,0.5,0.75))
ind.times.trn <- NULL
ind.times.trn[1] <- which.min(abs(surv.cox.trn$time - quants[1]))
ind.times.trn[2] <- which.min(abs(surv.cox.trn$time - quants[2]))
ind.times.trn[3] <- which.min(abs(surv.cox.trn$time - quants[3]))
risk_trn_cox_mat <- (surv.cox$surv[ind.times.trn,])

datas_trn1 <- data.frame(datas_trn,"risk"= t(risk_trn_cox_mat))

## Doing the same for test dataset
surv.cox.tst <- survfit(mod.cox,newdata = datas_tst)
quants.tst <- quantile(surv.cox.tst$time,probs = c(0.25,0.5,0.75))
ind.times.tst <- NULL
ind.times.tst[1] <- which.min(abs(surv.cox.tst$time - quants[1]))
ind.times.tst[2] <- which.min(abs(surv.cox.tst$time - quants[2]))
ind.times.tst[3] <- which.min(abs(surv.cox.tst$time - quants[3]))
risk_tst_cox_mat <- (surv.cox.tst$surv[ind.times.tst,])
datas_tst1 <- data.frame(datas_tst,"risk" = t(risk_tst_cox_mat))

## Combining the random forest model
mod.rsf1 <- rfsrc(formula = Surv(time, status) ~ ., data = datas_trn1,ntree = ntree,
                  nsplit = 0,splitrule = "logrank",importance = TRUE)
ypred.rsf1 <- predict(object = mod.rsf1, newdata = datas_tst1,type = "risk")
1-get.cindex(datas_tst1$time,datas_tst1$status,(ypred.rsf1$predicted))

###############
##Boosted Cox##
#####Tree######
###############

# No of models
ind_time_status <- NULL
ind_time_status[1] <- which(names(datas) == "time")
ind_time_status[2] <- which(names(datas) == "status")
ind_cov <- (1:ncol(datas))[-ind_time_status]

# Model parameters
b = 50
# Train test splitting ratio
# No of sampled covariates
alpha1 = 1
s_sam = round((ncol(datas)-2)*alpha1)
# No of sampled observations
alpha2 = (1-exp(-1))
n_sam = round(alpha_trn*alpha2*nrow(datas))

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
1-get.cindex(datas_tst$time,datas_tst$status,y_pred_mean_final)

## A different proposal
mod.cox <- coxph(formula = Surv(time, status) ~ ., data = datas_trn)
risk_tst_cox <- predict(object = mod.cox,newdata = datas_tst)
