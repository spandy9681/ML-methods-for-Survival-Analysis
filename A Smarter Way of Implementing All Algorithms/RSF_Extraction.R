# Playing with RSFs
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
x1 <- seq(-1,1,0.05)
x2 <- seq(-1,1,0.05)
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

# RSF Model
ntree <- 1
mod.rsf <- rfsrc(formula = Surv(time, status) ~ ., data = datas,ntree = ntree,nsplit = 0,splitrule = "logrank",importance = TRUE,nodedepth = 2,bootstrap = "none",samptype = "swor")
risk_tst_rsf <- predict(object = mod.rsf,newdata = datas)

# Extracting the tree information
plot(get.tree(mod.rsf,1))
mod.rsf$forest$nativeArray
direc <- mod.rsf$forest$nativeArray$contPT
direc <- as.vector(na.omit(direc))

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
datas_trn <- datas
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
  data_trn_sampled = as.data.frame(model.matrix(~.-1,data = data_trn_sampled))
  data_tst_sampled = as.data.frame(model.matrix(~.-1,data = data_tst_sampled))
  ## Hybrid Model Fitting ##
  mod.rsf <- rfsrc(formula = Surv(time, status) ~ ., data = data_trn_sampled,ntree = 1,
                    nsplit = 0,splitrule = "logrank",importance = TRUE,mtry = ncol(data_trn_sampled)-2,
                    bootstrap = "none",samptype = "swor",nodedepth = 1)
  ## Extracting information from the tree
  direc <- mod.rsf$forest$nativeArray$contPT
  direc <- as.vector(na.omit(direc))
  var <- mod.rsf$forest$nativeArray$parmID
  ## Splitting the dataset into each nodes
  ind_trn <- which(data_trn_sampled[,var[1]] <= direc)
  ind_tst1 <- which(data_tst_sampled[,var[1]] <= direc)
  ind_tst2 <- which(data_tst_sampled[,var[1]] > direc)
  ## Fitting seperate cox models for each node
  mod.cox1 <- coxph(formula = Surv(time, status) ~ ., data = data_trn_sampled[ind_trn,])
  mod.cox2 <- coxph(formula = Surv(time, status) ~ ., data = data_trn_sampled[-ind_trn,])
  ## Then predicting the risk values for the test dataset in a piecewise manner
  risk_tst_cox1 <- predict(object = mod.cox1,newdata = data_tst_sampled[ind_tst1,])
  risk_tst_cox2 <- predict(object = mod.cox2,newdata = data_tst_sampled[ind_tst2,])
  ## Now filling the matrix with the predicted risks in the appropriate places
  y_pred_mat[i,ind_tst1] <- risk_tst_cox1
  y_pred_mat[i,ind_tst2] <- risk_tst_cox2
}

y_pred_final <- apply(y_pred_mat, 2, mean)
## C-index
1-get.cindex(datas_tst$time,datas_tst$status,y_pred_final)

## Now again making the plot
data_risk_pred <- data.frame("x" = x,"pred_risk" = y_pred_final)
ggplot(data = data_risk_pred, mapping = aes(x.Var1, x.Var2)) + 
  geom_point(aes(colour = pred_risk), shape = 19) +
  scale_color_viridis(option = "D") +
  ggtitle("Prediction Using Our 
  Proposed Model")


## Bagged Cox model
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
alpha2 = 0.4
n_sam = round(alpha_trn*alpha2*nrow(datas))
cox_forest_boosted = NULL
set.seed(108)
# define the training and testing datasets
train_ind <- sample(1:nrow(datas),size = alpha_trn*(nrow(datas)))
datas_trn <- datas
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
  data_trn_sampled = as.data.frame(model.matrix(~.-1,data = data_trn_sampled))
  data_tst_sampled = as.data.frame(model.matrix(~.-1,data = data_tst_sampled))
  ## Cox Model Fitting ##
  mod.cox <- coxph(formula = Surv(time, status) ~ ., data = data_trn_sampled)
  risk_tst_cox <- predict(object = mod.cox,newdata = data_tst_sampled,type = "risk")
  y_pred_mat[i,] <- risk_tst_cox
}
y_pred_final <- apply(y_pred_mat, 2, mean)
1-get.cindex(datas_tst$time,datas_tst$status,y_pred_final)

## Now again making the plot
data_risk_pred <- data.frame("x" = x,"pred_risk" = y_pred_final)
ggplot(data = data_risk_pred, mapping = aes(x.Var1, x.Var2)) + 
  geom_point(aes(colour = pred_risk), shape = 19) +
  scale_color_viridis(option = "D") +
  ggtitle("Prediction Using Our 
  Proposed Model")

## Now again making the plot
data_risk_pred <- data.frame("x" = x,"pred_risk" = y_pred_mat[4,])
ggplot(data = data_risk_pred, mapping = aes(x.Var1, x.Var2)) + 
  geom_point(aes(colour = pred_risk), shape = 19) +
  scale_color_viridis(option = "D") +
  ggtitle("Prediction Using Our 
  Proposed Model")
