#tim DFT features models MTX compounds Multi-caloric Magnets
#check against multiple training sets
#3-4-21
#University of Virginia
#Charlottesville, VA 22903

##### CLEAN WORKSPACE ########
# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014") 
# Clean workspace
rm(list=ls())

#directroy definitions
DATA_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/Data'
MODEL_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/Models'
PLOT_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/Plots'
ANALYSIS_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/analysis'

#feature reduction by material science assumptions
dft.dataset <- read.csv(paste(DATA_DIR,'DFT_combined_4-7-21.csv',sep = '/'))
raw.features <- dft.dataset[,2:31]
reduced = na.omit(raw.features[,c(5:6,10:12,22,25,27,28)])
colnames(reduced)[9] = "response"

#### Support Vector Machines functions for bootstrap optimized model training ####
#libraries
library(matrixStats) # for rowSds
library(e1071) # the e1071 library contains implementations for a number of statistical learning methods.
library(parallel)

fn.cv.svr.rbf = function(data){
  tuneResult.svr.rbf = tune(svm, response ~ ., data = data,
                            kernel = "radial",
                            ranges = list(gamma = c(0.001,0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), 
                                          cost = c(0.001, 0.01, 0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100)),
                            tunecontrol = tune.control(cross =9))
  return(tuneResult.svr.rbf$best.parameters)
}

fn.train = function(data){
  parameter.svr.rbf = fn.cv.svr.rbf(data)
  svr.rbf.ecorr = svm(response ~ ., data = data, kernel = "radial", 
                      gamma =parameter.svr.rbf[1], cost = parameter.svr.rbf[2])
  return(svr.rbf.ecorr)
}

fn.boot.mean.oob = function(data,number_cores) {
  values <- c(30,50,100,150)
  oob.mean.err <- NULL
  oob.sd.err <- NULL
  for (i in 1:length(values)){
    print(values[i])
    oob <- NULL
    data.boot<-NULL
    for (j in 1:values[i]){
      n = dim(data)[1]
      bootsamples <- sample(n, replace = TRUE)
      data.boot[[j]] <- data[bootsamples,]
      oob[[j]] <- data[-bootsamples,]
    }
    oob <<- oob
    data.boot <<- data.boot
    oob.r2 <- NULL
    models <- mclapply(data.boot,fn.train,mc.cores = number_cores)
    models <<- models
    for (j in 1:length(models)){
      oob.r2 <- cor(oob[[j]]$response,predict(models[[j]],oob[[j]]))**2
    }
    oob.r2<<-oob.r2
    oob.mean.err[[i]] <- mean(oob.r2)
    oob.sd.err[[i]] <- sd(oob.r2)
    oob.mean <<- oob.mean.err
    if (i ==1 ){
      model <- models
    }else if (mean(oob.r2) == max(oob.mean.err)){
      model <- models
    }
  }
  oob.err.total <- as.data.frame(cbind(values,oob.mean.err,oob.sd.err))
  colnames(oob.err.total) <- c('number_bootstraps','mean R^2 error','SD R^2 error')
  oob.err.total <<- oob.err.total
  return(model)
} 

fn.svr.predict = function(model,data){#SVR bootstrapped prediction function is required, 
  predict.data = NULL
  i=1
  repeat{
    predict.data[[i]] = predict(model[[i]], data)
    i=i+1
    if(i>length(model)) break ()
  }
  if (nrow(data) == 1){ #previous rowmeans immplementation had a bug for single observation
    predict.data.mean = mean(predict.data)
  }
  else {
    predict.data.df = as.data.frame(predict.data)
    predict.data.mean = rowMeans(as.matrix(predict.data.df))
  }
  return(predict.data.mean)
}

fn.rsqrd = function(model,data){
  predictions <- fn.svr.predict(model,data)
  Rsqrd <- cor(predictions,data[,ncol(data)])**2
  return(Rsqrd)
}

fn.ttsplit = function(data,seed){
  #test train split
  set.seed(seed)
  print(seed)
  sample = sample.int(n = nrow(data),size = floor(0.8*nrow(data)),replace = F)
  train = data[sample,]
  test = data[-sample,]
  return(list(train,test))
}

fn.robust.ttsplit = function(data,number_of_seeds,number_of_cores = 1){
  seeds = sample.int(4242,number_of_seeds)
  rsqrd <- NULL
  test_sets <- NULL
  models <- NULL
  for (i in 1:length(seeds)){
    print(paste('evaluating seed',seeds[i],sep = ' '))
    split = fn.ttsplit(data = data, seed = seeds[[i]])
    train = split[[1]]
    test = split[[2]]
    model = fn.boot.mean.oob(data = train,number_cores = number_of_cores)
    train_rsqrd <- fn.rsqrd(model,train)
    print(paste('train R^2:',train_rsqrd,sep = ' '))
    rsqrd[[i]] <- fn.rsqrd(model,test)
    print(paste('test R^2:',rsqrd[[i]],sep = ' '))
    test_sets[[i]] <- test
    models[[i]] <- model
  }
  test_sets <<-test_sets
  models <<-models
  return(rsqrd)
}

fn.temp.predict = function(data, model){
  predict.data = NULL
  i=1
  repeat{
    predict.data[[i]] = predict(model[[i]], data)
    i=i+1
    if(i>length(model)) break ()
  }
  predict.data.df = as.data.frame(predict.data)
  predict.data.mean = rowMeans(as.matrix(predict.data.df))
  predict.data.sds = rowSds(as.matrix(predict.data.df))
  data.predictions <- as.data.frame(cbind(predict.data.mean, predict.data.sds))
  names(data.predictions) <- c("prediction.mean","prediction.sd")
  return(data.predictions)
}

sample_performance <- fn.robust.ttsplit(reduced,number_of_seeds = 20,number_of_cores = 10)

best_model <- models[[4]]
save(best_model,file = paste(MODEL_DIR,'best.svr.dft.6-10-21.rda',sep = '/'))

best_sample <- c("1","12","16","21","23","24","25")
theating.train <- reduced[!(rownames(reduced)) %in% best_sample,]
theating.test <- reduced[best_sample,]


#performance
load(file = paste(MODEL_DIR,'best.svr.dft.6-10-21.rda',sep = '/'))
svr.train.predictions <- fn.temp.predict(theating.train,best_model)
svr.train.performance <- cbind(theating.train$response,svr.train.predictions)
colnames(svr.train.performance)[1] <- 'T_heating'
svr.train.performance$residuals = (svr.train.performance$T_heating-svr.train.performance$prediction.mean)**2
svr.train.rmse = sqrt(mean(svr.train.performance$residuals))
train.rsqr = cor(svr.train.performance$T_heating, svr.train.performance$prediction.mean)**2

svr.test.predictions <- fn.temp.predict(theating.test,best_model)
svr.test.performance <- cbind(theating.test$response,svr.test.predictions)
colnames(svr.test.performance)[1] <- 'T_heating'
svr.test.performance$residuals = (svr.test.performance$T_heating-svr.test.performance$prediction.mean)**2
svr.test.rmse = sqrt(mean(svr.test.performance$residuals))
test.rsqr = cor(svr.test.performance$T_heating, svr.test.performance$prediction.mean)**2

svr.performance = rbind(svr.train.performance,svr.test.performance)

write.csv(svr.performance, file = paste(ANALYSIS_DIR,'bestDFTmodel_6-7-21svr.csv',sep = '/'))

#svr
library(ggplot2)
labels = c(rep('train',25),rep('test',7))
svr.performance['label'] <- labels

p = ggplot(svr.performance,aes(x=T_heating,y=prediction.mean,col=label))+geom_point()
p_error = p+geom_errorbar(data = svr.performance,aes(x=T_heating, ymin=prediction.mean-prediction.sd*1.9, ymax=prediction.mean+prediction.sd*1.9, width=0.25))
p_final = p_error+geom_abline()
p_final

#Hysteresis
reduced = na.omit(raw.features[,c(5:6,10:12,22,25,27,30)])
colnames(reduced)[9] <- "response"
sample_performance <- fn.robust.ttsplit(reduced,number_of_seeds = 20,number_of_cores = 10)

best_model <- models[[8]]
save(best_model,file = paste(MODEL_DIR,'best.hyst.dft.6-7-21.rda',sep = '/'))

print(which(sample_performance==max(sample_performance)))
best_sample <- c("2","4","9","11","14","30")
theating.train <- reduced[!(rownames(reduced)) %in% best_sample,]
theating.test <- reduced[best_sample,]

#performance
load(file = paste(MODEL_DIR,'best.hyst.dft.6-7-21.rda',sep = '/'))
svr.train.predictions <- fn.temp.predict(theating.train,best_model)
svr.train.performance <- cbind(theating.train$response,svr.train.predictions)
colnames(svr.train.performance)[1] <- 'T_heating'
svr.train.performance$residuals = (svr.train.performance$T_heating-svr.train.performance$prediction.mean)**2
svr.train.rmse = sqrt(mean(svr.train.performance$residuals))
train.rsqr = cor(svr.train.performance$T_heating, svr.train.performance$prediction.mean)**2

svr.test.predictions <- fn.temp.predict(theating.test,best_model)
svr.test.performance <- cbind(theating.test$response,svr.test.predictions)
colnames(svr.test.performance)[1] <- 'T_heating'
svr.test.performance$residuals = (svr.test.performance$T_heating-svr.test.performance$prediction.mean)**2
svr.test.rmse = sqrt(mean(svr.test.performance$residuals))
test.rsqr = cor(svr.test.performance$T_heating, svr.test.performance$prediction.mean)**2

svr.performance = rbind(svr.train.performance,svr.test.performance)

write.csv(svr.performance, file = paste(ANALYSIS_DIR,'bestDFThyst_6-10-21svr.csv',sep = '/'))
