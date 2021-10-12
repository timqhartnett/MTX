#tim DFT features models MTX compounds Multi-caloric Magnets
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

#### Support Vector Machines functions for bootstrap optimized model training ####
#libraries
library(matrixStats) # for rowSds
library(e1071) # the e1071 library contains implementations for a number of statistical learning methods.
library(parallel)
library(matrixStats)

### function definitions #### this is poor coding formating...
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

fn.svr.jacknife.regression = function(data,test.data,nboot,number_cores) {
  print('global feature importance determination by leave one out methodology')
  features = colnames(data[,-ncol(data)])
  print(features)
  oob.mean.err <- NULL
  oob.sd.err <- NULL
  test.mean.err <- NULL
  test.sd.err <- NULL
  for (i in 1:length(features)){
    print(paste('feature withheld', features[i],sep = ' = '))
    tmp.data = data[,-i]
    oob <- NULL
    data.boot<-NULL
    for (j in 1:nboot){
      n = dim(tmp.data)[1]
      bootsamples <- sample(n, replace = TRUE)
      data.boot[[j]] <- tmp.data[bootsamples,]
      oob[[j]] <- tmp.data[-bootsamples,]
    }
    oob <<- oob
    models <- mclapply(data.boot,fn.train,mc.cores = number_cores)
    oob.rms.err <- NULL
    test.rms.err <- NULL
    for (j in 1:length(oob)){
      oob.rms.err[[j]] <- sqrt(sum((oob[[j]][,ncol(oob[[j]])]-predict(models[[j]],oob[[j]]))**2)/nrow(oob[[j]]))
      test.rms.err[[j]] <- sqrt(sum(((test.data[,ncol(test.data)])-predict(models[[j]],test.data))**2)/nrow(test.data))
    }
    oob.mean.err[[i]] <- mean(oob.rms.err)
    oob.sd.err[[i]] <- sd(oob.rms.err)
    test.mean.err[[i]] <- mean(test.rms.err)
    test.sd.err[[i]] <- sd(test.rms.err)
  }
  oob.err.total <- as.data.frame(cbind(features,oob.mean.err,oob.sd.err,test.mean.err,test.sd.err))
  colnames(oob.err.total) <- c('number_bootstraps','OOB.mean.RMS.error','OOB.SD.RMS.error','Test.mean.RMS.error','Test.sd.RMS.error')
  return(oob.err.total)
} 


##### Tt ###### 
#feature reduction by material science assumptions
dft.dataset <- read.csv(paste(DATA_DIR,'DFT_combined_4-7-21.csv',sep = '/'))
raw.features <- dft.dataset[,2:31]
reduced = na.omit(raw.features[,c(5:6,10:12,22,25,27,28)])
colnames(reduced)[9] = "response"

best_sample <- c("1","12","16","21","23","24","25")
theating.train <- reduced[!(rownames(reduced)) %in% best_sample,]
theating.test <- reduced[best_sample,]

feature.importance.tt <- fn.svr.jacknife.regression(theating.train,theating.test,nboot=100,number_cores = 10)
write.csv(feature.importance.tt,paste(ANALYSIS_DIR,'svr_best_dft_importance_6-8-21.csv',sep = '/'))
#### hysteresis ####
reduced = na.omit(raw.features[,c(5:6,10:12,22,25,27,30)])
colnames(reduced)[9] <- "response"

best_sample <- c("2","4","9","11","14","30")
theating.train <- reduced[!(rownames(reduced)) %in% best_sample,]
theating.test <- reduced[best_sample,]
feature.importance.hyst <- fn.svr.jacknife.regression(theating.train,theating.test,nboot=150,number_cores = 10)

write.csv(feature.importance.hyst,paste(ANALYSIS_DIR,'hyst_dft_importance_6-8-21.csv',sep = '/'))
