set.seed(108)
# define the training and testing datasets
train_ind <- sample(1:nrow(datas),size = alpha_trn*(nrow(datas)))
datas_trn <- datas[train_ind,]
datas_trn <- cbind(datas_trn,"var.new"=rep(0,nrow(datas_trn)))
mod.cox <- coxph(formula = Surv(time, status) ~ ., data = datas_trn)
risk_trn_cox <- predict(object = mod.cox,newdata = datas_trn)
datas_trn1 <- data.frame(datas_trn,"risk"= (risk_trn_cox))
mod.rsf1 <- rfsrc(formula = Surv(time, status) ~ ., data = datas_trn1,ntree = 10,
                  nsplit = 0,splitrule = "logrank",importance = TRUE,mtry = ncol(datas_trn1)-2,
                  nodesize = 5)
ypred.rsf1 <- predict(object = mod.rsf1, newdata = datas_trn1,type = "risk")
## Then again repeat the process
datas_trn[,ncol(datas_trn)] <- ypred.rsf1$predicted
mod.cox <- coxph(formula = Surv(time, status) ~ ., data = datas_trn)
risk_trn_cox <- predict(object = mod.cox,newdata = datas_trn)
datas_trn1 <- data.frame(datas_trn,"risk"= (risk_trn_cox))
mod.rsf1 <- rfsrc(formula = Surv(time, status) ~ ., data = datas_trn1,ntree = 1,
                  nsplit = 0,splitrule = "logrank",importance = TRUE,mtry = ncol(datas_trn1)-2,
                  bootstrap = "none",samptype = "swor",nodesize = 50)
ypred.rsf1 <- predict(object = mod.rsf1, newdata = datas_trn1,type = "risk")
ypred.rsf1.tst <- predict(object = mod.rsf1, newdata = datas_tst,type = "risk")
