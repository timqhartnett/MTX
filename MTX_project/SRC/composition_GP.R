#Materials Informatics --2-27-21
#Gaussian Process Regression MnNiSi
#pvb5e@virginia.edu
#University of Virginia
#Charlottesville, VA 22903
##### CLEAN WORKSPACE ########
# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014") 
# Clean workspace
rm(list=ls())

##libraries
library(parallel)
library(GPfit)
library(ggplot2)

#### DIRECTORY DEFINITIONS ############
DATA_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/Data'
MODEL_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/Models'
PLOT_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/Plots'
ANALYSIS_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/analysis'

#DATA
library(readxl)
comp.dataset <- read_xlsx(paste(DATA_DIR,'composition_5_27_21.xlsx',sep = '/'))
reduced.dataset <- as.data.frame(na.omit(comp.dataset[,-c(1,12:14)]))
colnames(reduced.dataset)[10] <- 'response'
training_y <- reduced.dataset[, 10]
input <- reduced.dataset[,1:9]

#response <- elemental_features[,c(7)]
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
x1 = range01(input[,c(1)])
x2 = range01(input[,c(2)])
x3 = range01(input[,c(3)])
x4 = range01(input[,c(4)])
x5 = range01(input[,c(5)])
x6 = range01(input[,c(6)])
x7 = range01(input[,c(7)])
x8 = range01(input[,c(8)])
x9 = range01(input[,c(9)])

full_input<-as.data.frame(cbind(x1, x2, x3, x4, x5, x6, x7, x8, x9))
input_train <- full_input

input_all = cbind(input_train, training_y)

#test train split
set.seed(42)
sample = sample.int(n = nrow(input_all),size = floor(0.8*nrow(input_all)),replace = F)
train = input_all[sample,]
test = input_all[-sample,]

x_train <- train[,c(1:9)]
y_train <-train[,c(10)]

x_test <- test[,c(1:9)]
y_test <- test[,c(10)]

## GP TRAINING ##
d=3
nugget_optim <- function(param){
  nug = param[1]
  GPmodel_nug = GP_fit(x_train, y_train, control = c(200*d, 80*d, 2*d), nug_thres = nug)
  print(nug)
  GPmodel_resub = predict(GPmodel_nug, x_train)
  ll = dnorm(y_train, mean=GPmodel_resub$Y_hat, sd=GPmodel_resub$MSE, log=T)
  sumll = sum(ll)
  return(sumll)
}

nugget_value = function(x){return(nugget_optim(x))}

ptm <- proc.time()
slope = lapply(seq(1,20, by=1), nugget_value)
print(proc.time() - ptm)  #### SYSTEM TIME MONITORING ####

plot(seq(1,20, by=1), slope, type="o") ### Look for the nugget with the largest slope
slope2 = as.matrix(slope)
slope3 = as.data.frame(cbind(seq(1:20), slope))
slope3.df <- data.frame(matrix(unlist(slope3), nrow=20))
slope3.df.order <-slope3.df[order(slope3.df$X2, decreasing=TRUE),]
best_nug <- slope3.df.order[1,1]

## BEST GP MODEL ##
best_nug = 7
GPmodel_optim = GP_fit(x_train, y_train, control = c(200*d, 80*d, 2*d), nug_thres = best_nug)
save(GPmodel_optim, file = paste(MODEL_DIR,'5feature.gp.rda',sep = '/'))

#performance
load(paste(MODEL_DIR,'5feature.gp.rda',sep = '/'))
GPmodel_insample = predict(GPmodel_optim, x_train)

me = as.data.frame(GPmodel_insample$complete_data)
me = cbind(y_train, me)
p = qplot(me$y_train, me$Y_hat)+xlim(50,600)+ylim(50,600)+
  geom_errorbar(aes(x=me$y_train, ymin=me$Y_hat-me$MSE, ymax=me$Y_hat+me$MSE), width=0.25)
p+geom_abline()

gp.train.performance = cbind(y_train,me[,7:8])
colnames(gp.train.performance) <- c('Theating','prediction.mean','prediction.sd')
gp.train.performance$residuals = (gp.train.performance$Theating-gp.train.performance$prediction.mean)**2
gp.train.rmse = sqrt(mean(gp.train.performance$residuals))
train.rsqr = cor(gp.train.performance$Theating,gp.train.performance$prediction.mean)**2

#GPmodel_optim = GP_fit(x_train, y_train, control = c(200*d, 80*d, 2*d), nug_thres = best_nug)
GPmodel_out_of_sample = predict(GPmodel_optim, x_test)

me_test = as.data.frame(GPmodel_out_of_sample$complete_data)
me_test = cbind(y_test, me_test)
p2_test = qplot(me_test$y_test, me_test$Y_hat)+xlim(-300,1000)+ylim(-300,1000)+
  geom_errorbar(aes(x=me_test$y_test, ymin=me_test$Y_hat-me_test$MSE, ymax=me_test$Y_hat+me_test$MSE), width=0.25)
p2_test+geom_abline()

gp.test.performance = cbind(y_test,me_test[,7:8])
colnames(gp.test.performance) <- c('Theating','prediction.mean','prediction.sd')
gp.test.performance$residuals = (gp.test.performance$Theating-gp.test.performance$prediction.mean)**2
gp.test.rmse = sqrt(mean(gp.test.performance$residuals))
gp.test.rsqrd = cor(gp.test.performance$Theating,gp.test.performance$prediction.mean)**2
print(rsq(me_test$y_test, me_test$Y_hat))

gp.performance = rbind(gp.train.performance,gp.test.performance)
write.csv(gp.performance,file = paste(ANALYSIS_DIR,'gp_DFT.csv',sep = '/'))
