rm(list = ls())

# Required Libraries

library(ggplot2)
library(ISLR2)
library(survival)
library(randomForestSRC)
library(viridis)

# Any Survival data with censoring index and survival times coded as "status" & "time" 
# respectively put your dataset here || ||  
#                                 || ||
#                                 \/ \/

## Simulated Non-Linear Dataset 1 ##
x1 <- seq(-1,1,0.1)
x2 <- seq(-1,1,0.1)
d <- expand.grid(x1,x2)
x = as.matrix(d)
c0 = 5
beta = c(4,4)
h = exp((-x^2%*%beta))
#h = x%*%beta
n <- nrow(x)
times = rexp(n,c0*exp(h))
time.censor = rexp(n,c0*exp(h))
summary(time.censor)
delta <- ifelse(times<time.censor, 1, 0)
time <- ifelse(times<time.censor, times, time.censor)
datas <- data.frame("x" = x,time = time,status = delta)

### Plotting the original risk surface
data_risk <- data.frame("x" = x,"true_risk" = h)
ggplot(data = data_risk, mapping = aes(x.Var1, x.Var2)) + 
  geom_point(aes(colour = true_risk), shape = 19) +
  scale_color_viridis(option = "D") +
  ggtitle("True Risk Surface")


## Only running the model once
ind_time_status <- NULL
ind_time_status[1] <- which(names(datas) == "time")
ind_time_status[2] <- which(names(datas) == "status")
ind_cov <- (1:ncol(datas))[-ind_time_status]

# No of models
b = 50
# Train test splitting ratio
alpha_trn = 0.6
# No of sampled covariates
alpha1 = 1
s_sam = round((ncol(datas)-2)*alpha1)
# No of sampled observations
alpha2 = (1-exp(-1))
n_sam = round(alpha_trn*alpha2*nrow(datas))
cox_forest_boosted = NULL
set.seed(108)
# define the training and testing datasets
train_ind <- sample(1:nrow(datas),size = alpha_trn*(nrow(datas)))
datas_trn <- datas[train_ind,]
datas_tst <- datas
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
                    bootstrap = "none",samptype = "swor",nodesize = 50)
  risk_tst_cox <- predict(object = mod.cox,newdata = data_tst_sampled)
  datas_tst1 <- data.frame(data_tst_sampled,"risk"= (risk_tst_cox))
  ## Putting the fitted value in a matrix
  ypred.rsf1 <- predict(object = mod.rsf1, newdata = datas_tst1,type = "risk")
  y_pred_mat[i,] <- ypred.rsf1$predicted
}
y_pred_final <- apply(y_pred_mat, 2, mean)
## C-index
1-get.cindex(datas_tst$time,datas_tst$status,y_pred_final)
get.cindex(h,datas_tst$status,y_pred_final)

## Now again making the plot
data_risk_pred <- data.frame("x" = x,"pred_risk" = y_pred_final)
ggplot(data = data_risk_pred, mapping = aes(x.Var1, x.Var2)) + 
  geom_point(aes(colour = pred_risk), shape = 19) +
  scale_color_viridis(option = "D") +
  ggtitle("Prediction Using Our 
  Proposed Model")


## Comparing with cox-ph and RSF
# training the rsf model
ntree <- b
mod.rsf <- rfsrc(formula = Surv(time, status) ~ ., data = datas_trn,ntree = ntree,nsplit = 0,splitrule = "logrank",importance = TRUE,nodesize = 50)
risk_tst_rsf <- predict(object = mod.rsf,newdata = datas_tst)

# And calculate the concordance index
1-get.cindex(datas_tst$time,datas_tst$status,(risk_tst_rsf$predicted))
get.cindex(h,datas_tst$status,(risk_tst_rsf$predicted))

data_risk_pred_rsf <- data.frame("x" = x,"pred_risk" = risk_tst_rsf$predicted)
ggplot(data = data_risk_pred_rsf, mapping = aes(x.Var1, x.Var2)) + 
  geom_point(aes(colour = pred_risk), shape = 19) +
  scale_color_viridis(option = "D") +
  ggtitle("Using Random Survival 
       Forest Model")

### CoxPh ###
mod.cox <- coxph(formula = Surv(time, status) ~ x.Var1^2+x.Var2^2, data = datas_trn)
risk_tst_cox <- predict(object = mod.cox,newdata = datas_tst)

1-get.cindex(datas_tst$time,datas_tst$status,(risk_tst_cox))
get.cindex(h,datas_tst$status,(risk_tst_cox))

data_risk_pred_cox <- data.frame("x" = x,"pred_risk" = risk_tst_cox)
ggplot(data = data_risk_pred_cox, mapping = aes(x.Var1, x.Var2)) + 
  geom_point(aes(colour = pred_risk), shape = 19) +
  scale_color_viridis(option = "D") +
  ggtitle("Using CoxPh Model")
