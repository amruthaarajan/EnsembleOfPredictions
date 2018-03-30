#********************************************************************************
#********************************************************************************
#                    Acquiring the Predicted Probabilities
#  We run the code for T0 that we have (3 models) to get predicted probabilities
#             for each sample id, for each model through Xiao's CV
#********************************************************************************
#********************************************************************************





#*****************    Model 1: ViResPred     *****************
#predictions from https://www.synapse.org/#!Synapse:syn7208407
#*************************************************************

#install.packages("ada", repos = "http://cran.us.r-project.org")
library(data.table)
library(e1071) #install e1071
library(ada) #install ada
library(randomForest) #install randomForest
library(class)

#==========shuffle==========
## Input 
## data: dataframe that needs to be shuffled
## k: number of folds

## Return
## new indeces
set.seed(588)

shuffle <- function(data, k){
  ## Randomly shuffle the data
  indeces <- seq(1, length(data[,1]), by = 1)
  indeces <- sample(indeces)
  ## k fold
  cut(indeces, breaks = k, labels = FALSE)
}

#==========Feature selection=================
## Input:
## file.name : data frame

## Return:
## a vector of p-values

this.is.t.test <- function(file.name){
  data.feature <- file.name$SYMPTOMATIC_SC2
  #print(data.feature)
  pos.index <- which(data.feature == 1)
  neg.index <- which(data.feature == 0)
  
  #print(apply(file.name,2,function(x){length(unique(x))}))
  
  apply (file.name[,names(file.name) != "SYMPTOMATIC_SC2"], 2 , function (x){
    t.test(x[neg.index], x[pos.index])$p.value})
}

#==========Model=================
# it gets train and test data as input to train the model and return the prediction result

####KNN######
nb_classification<-function(train_fulldata,test_fulldata) {
  train_data<-train_fulldata[, !names(train_fulldata) %in% c("SYMPTOMATIC_SC2")]
  train_result<- train_fulldata$SYMPTOMATIC_SC2
  test_data<- test_fulldata[, !names(test_fulldata) %in% c("SYMPTOMATIC_SC2")]
  test_result<-test_fulldata$SYMPTOMATIC_SC2
  #test
  print(dim(train_data))
  print(dim(test_data))
  knn_model<-knn(train=as.matrix(train_data), test=as.matrix(test_data), cl=train_result, k=5, prob=TRUE)
  knn_model
  prob <- attr(knn_model, "prob")
  pred <- as.numeric(as.character(knn_model))
  ind <- which(pred == 0)
  prob[ind] <- 1 - prob[ind]
  #pred <- prob
  return(as.numeric(as.vector(prob)))
}

#==========kfoldCV==========
## Input 
## data: dataframe contains predictors and target variable
## round: number of rounds we want to repeat cross-validation

## Model:the function name of the model:
## logistic, svm.run, ada.run, forest, knn.run   e.g. logistic(train, test)
## K: number of folds for crossValidation
## X: number of sample from each of y=0 and y=1 

## Return
## average true_auroc & average true_auprc

kfoldCV <- function(data, round, Model, k, X) {
  ## Remove column TIMEHOURS
  #data <- data[, !names(data) %in% c("TIMEHOURS")]
  
  for(j in 1:round) {
    ## For each round, shuffle the data first
    fold.indeces <-  shuffle(data, k)
    
    ## Perform k fold cross validation
    predict <- NULL
    target <- NULL
    for (i in 1:k) {
      index <- which(fold.indeces == i, arr.ind = TRUE)
      train <- data[-index, ]
      test <- data[index, ]
      
      #### bagging
      ## UnderSampling of both class
      dataP = train[which(train$SYMPTOMATIC_SC2 == 1),]
      dataN = train[which(train$SYMPTOMATIC_SC2 == 0),]
      sampleDataP = dataP[sample(x=1:nrow(dataP),size=X,replace=TRUE),]
      sampleDataN = dataN[sample(x=1:nrow(dataN),size=X,replace=TRUE),]
      train <- rbind(sampleDataP,sampleDataN)
      train <- train[sample(nrow(train)),]
      
      
      print("feature selection")
      ####################### feature selection #######################
      ## Define the p-value threshold
      p.threshold <- 0.1
      response <- "SYMPTOMATIC_SC2"
      
      p.value <- this.is.t.test(train)
      
      column.id <- c(which(p.value < p.threshold), which(names(train) == response))
      
      train <- train[, column.id]
      test <- test[, column.id]
      print(dim(train))
      
      #################################################################
      
      ## Train classification model on data using selected features
      ## predicted values
      print("predict model begin")
      new_pred_value<-nb_classification(train,test)
      print(new_pred_value)
      print("predict finish")
      
      predict <- c(predict, new_pred_value)
      print(predict)
      ## actual values
      target <- c(target, test$SYMPTOMATIC_SC2)
      print(target)
      write.csv(predict, "predict.csv")
      write.csv(target, "target.csv")
    }
  }
}

#==========Main Function==========
# read the reult file and data file, and combine them to one table
result_file<-read.table("ViralChallenge_training_CLINICAL.tsv",sep="\t",header=TRUE)

data_file <- fread("ViralChallenge_training_EXPRESSION_RMA.tsv", data.table=FALSE)
rownames(data_file) <- data_file$FEATUREID
data_file <- data_file[, -1]
final_data_file<-as.data.frame(t(data_file))
final_data_file$SYMPTOMATIC_SC2<-result_file$SYMPTOMATIC_SC2

accuracy<- kfoldCV(data = final_data_file,round = 1,Model = "NB",k = 5,X = round(dim(final_data_file)[1]/2))




#*****************    Model 2: Agnita     *****************
#Bin Zhang: submission2-aganita_viral_code
#Prediction from https://www.synapse.org/#!Synapse:syn7203047
#*************************************************************

#install.packages("ada", repos = "http://cran.us.r-project.org")
library(data.table)
library(e1071) #install e1071
library(ada) #install ada
library(randomForest) #install randomForest
library(class)
library(glmnet)

#==========shuffle==========
## Input 
## data: dataframe that needs to be shuffled
## k: number of folds

## Return
## new indeces
set.seed(588)

shuffle <- function(data, k){
  ## Randomly shuffle the data
  indeces <- seq(1, length(data[,1]), by = 1)
  indeces <- sample(indeces)
  ## k fold
  cut(indeces, breaks = k, labels = FALSE)
}

#==========Feature selection=================
## Input:
## file.name : data frame

## Return:
## a vector of p-values

this.is.t.test <- function(file.name){
  data.feature <- file.name$SYMPTOMATIC_SC2
  #print(data.feature)
  pos.index <- which(data.feature == 1)
  neg.index <- which(data.feature == 0)
  
  #print(apply(file.name,2,function(x){length(unique(x))}))
  
  apply (file.name[,names(file.name) != "SYMPTOMATIC_SC2"], 2 , function (x){
    t.test(x[neg.index], x[pos.index])$p.value})
}

#==========Model=================
# it gets train and test data as input to train the model and return the prediction result

####Lasso Regression######
nb_classification<-function(train_fulldata,test_fulldata) {
  train_data<-train_fulldata[, !names(train_fulldata) %in% c("SYMPTOMATIC_SC2")]
  train_result<- train_fulldata$SYMPTOMATIC_SC2
  test_data<- test_fulldata[, !names(test_fulldata) %in% c("SYMPTOMATIC_SC2")]
  test_result<-test_fulldata$SYMPTOMATIC_SC2
  #test
  print(dim(train_data))
  print(dim(test_data))
  
  lr_model<-cv.glmnet(as.matrix(train_data), as.factor(train_result), alpha=1, family= "binomial")
  pred<-predict(lr_model, as.matrix(test_data), s=lr_model$lambda.min, type="response")
  return(as.numeric(as.vector(pred)))
}

#==========kfoldCV==========
## Input 
## data: dataframe contains predictors and target variable
## round: number of rounds we want to repeat cross-validation

## Model:the function name of the model:
## logistic, svm.run, ada.run, forest, knn.run   e.g. logistic(train, test)
## K: number of folds for crossValidation
## X: number of sample from each of y=0 and y=1 

## Return
## average true_auroc & average true_auprc

kfoldCV <- function(data, round, Model, k, X) {
  ## Remove column TIMEHOURS
  #data <- data[, !names(data) %in% c("TIMEHOURS")]
  
  for(j in 1:round) {
    ## For each round, shuffle the data first
    fold.indeces <-  shuffle(data, k)
    
    ## Perform k fold cross validation
    predict <- NULL
    target <- NULL
    for (i in 1:k) {
      index <- which(fold.indeces == i, arr.ind = TRUE)
      train <- data[-index, ]
      test <- data[index, ]
      
      #### bagging
      ## UnderSampling of both class
      dataP = train[which(train$SYMPTOMATIC_SC2 == 1),]
      dataN = train[which(train$SYMPTOMATIC_SC2 == 0),]
      sampleDataP = dataP[sample(x=1:nrow(dataP),size=X,replace=TRUE),]
      sampleDataN = dataN[sample(x=1:nrow(dataN),size=X,replace=TRUE),]
      train <- rbind(sampleDataP,sampleDataN)
      train <- train[sample(nrow(train)),]
      
      
      print("feature selection")
      ####################### feature selection #######################
      ## Define the p-value threshold
      p.threshold <- 0.1
      response <- "SYMPTOMATIC_SC2"
      
      p.value <- this.is.t.test(train)
      
      column.id <- c(which(p.value < p.threshold), which(names(train) == response))
      
      train <- train[, column.id]
      test <- test[, column.id]
      print(dim(train))
      
      #################################################################
      
      ## Train classification model on data using selected features
      ## predicted values
      print("predict model begin")
      new_pred_value<-nb_classification(train,test)
      print(new_pred_value)
      print("predict finish")
      
      predict <- c(predict, new_pred_value)
      print(predict)
      ## actual values
      target <- c(target, test$SYMPTOMATIC_SC2)
      print(target)
      write.csv(predict, "predict2.csv")
      write.csv(target, "target2.csv")
    }
  }
}

#==========Main Function==========
# read the reult file and data file, and combine them to one table
result_file<-read.table("ViralChallenge_training_CLINICAL.tsv",sep="\t",header=TRUE)

data_file <- fread("ViralChallenge_training_EXPRESSION_RMA.tsv", data.table=FALSE)
rownames(data_file) <- data_file$FEATUREID
data_file <- data_file[, -1]
final_data_file<-as.data.frame(t(data_file))
final_data_file$SYMPTOMATIC_SC2<-result_file$SYMPTOMATIC_SC2

accuracy<- kfoldCV(data = final_data_file,round = 1,Model = "NB",k = 5,X = round(dim(final_data_file)[1]/2))



#*****************    Model 3: Shosty     *****************
#dingyuyang_3654505_42777971_individual_1 (1) 
##Prediction from https://www.synapse.org/#!Synapse:syn7205111
#*************************************************************

#==========shuffle==========
## Input 
## data: dataframe that needs to be shuffled
## k: number of folds

## Return
## new indeces
set.seed(588)

shuffle <- function(data, k){
  ## Randomly shuffle the data
  indeces <- seq(1, length(data[,1]), by = 1)
  indeces <- sample(indeces)
  ## k fold
  cut(indeces, breaks = k, labels = FALSE)
}

#==========Feature selection=================
## Input:
## file.name : data frame

## Return:
## a vector of p-values

this.is.t.test <- function(file.name){
  data.feature <- file.name$SYMPTOMATIC_SC2
  #print(data.feature)
  pos.index <- which(data.feature == 1)
  neg.index <- which(data.feature == 0)
  
  #print(apply(file.name,2,function(x){length(unique(x))}))
  
  apply (file.name[,names(file.name) != "SYMPTOMATIC_SC2"], 2 , function (x){
    t.test(x[neg.index], x[pos.index])$p.value})
}

#==========Model=================
# it gets train and test data as input to train the model and return the prediction result

####SVM######
nb_classification<-function(train_fulldata,test_fulldata) {
  train_data<-train_fulldata[, !names(train_fulldata) %in% c("SYMPTOMATIC_SC2")]
  train_result<- train_fulldata$SYMPTOMATIC_SC2
  test_data<- test_fulldata[, !names(test_fulldata) %in% c("SYMPTOMATIC_SC2")]
  test_result<-test_fulldata$SYMPTOMATIC_SC2
  #test
  print(dim(train_data))
  print(dim(test_data))
  
  svm.model <- svm (x = train_data,y= train_result, scale=FALSE,kernel="linear",cost=2,class.weights=c("0"=2,"1"=1),cachesize=600,probability=TRUE)
  predict.value <- predict(svm.model,  test_data)
  return(as.numeric(as.vector(predict.value)))
  #pred = predict(fit, type="response")
}

#==========kfoldCV==========
## Input 
## data: dataframe contains predictors and target variable
## round: number of rounds we want to repeat cross-validation

## Model:the function name of the model:
## logistic, svm.run, ada.run, forest, knn.run   e.g. logistic(train, test)
## K: number of folds for crossValidation
## X: number of sample from each of y=0 and y=1 

## Return
## average true_auroc & average true_auprc

kfoldCV <- function(data, round, Model, k, X) {
  ## Remove column TIMEHOURS
  #data <- data[, !names(data) %in% c("TIMEHOURS")]
  
  for(j in 1:round) {
    ## For each round, shuffle the data first
    fold.indeces <-  shuffle(data, k)
    
    ## Perform k fold cross validation
    predict <- NULL
    target <- NULL
    for (i in 1:k) {
      index <- which(fold.indeces == i, arr.ind = TRUE)
      train <- data[-index, ]
      test <- data[index, ]
      
      #### bagging
      ## UnderSampling of both class
      dataP = train[which(train$SYMPTOMATIC_SC2 == 1),]
      dataN = train[which(train$SYMPTOMATIC_SC2 == 0),]
      sampleDataP = dataP[sample(x=1:nrow(dataP),size=X,replace=TRUE),]
      sampleDataN = dataN[sample(x=1:nrow(dataN),size=X,replace=TRUE),]
      train <- rbind(sampleDataP,sampleDataN)
      train <- train[sample(nrow(train)),]
      
      
      print("feature selection")
      ####################### feature selection #######################
      ## Define the p-value threshold
      p.threshold <- 0.1
      response <- "SYMPTOMATIC_SC2"
      
      p.value <- this.is.t.test(train)
      
      column.id <- c(which(p.value < p.threshold), which(names(train) == response))
      
      train <- train[, column.id]
      test <- test[, column.id]
      print(dim(train))
      
      #################################################################
      
      ## Train classification model on data using selected features
      ## predicted values
      print("predict model begin")
      new_pred_value<-nb_classification(train,test)
      print(new_pred_value)
      print("predict finish")
      
      predict <- c(predict, new_pred_value)
      print(predict)
      ## actual values
      target <- c(target, test$SYMPTOMATIC_SC2)
      print(target)
      write.csv(predict, "predict3.csv")
      write.csv(target, "target3.csv")
    }
  }
}

#==========Main Function==========
# read the reult file and data file, and combine them to one table
result_file<-read.table("ViralChallenge_training_CLINICAL.tsv",sep="\t",header=TRUE)

data_file <- fread("ViralChallenge_training_EXPRESSION_RMA.tsv", data.table=FALSE)
rownames(data_file) <- data_file$FEATUREID
data_file <- data_file[, -1]
final_data_file<-as.data.frame(t(data_file))
final_data_file$SYMPTOMATIC_SC2<-result_file$SYMPTOMATIC_SC2

accuracy<- kfoldCV(data = final_data_file,round = 1,Model = "NB",k = 5,X = round(dim(final_data_file)[1]/2))

#create the prediction probability table for the models
pred1 <- read.csv("predict1.csv")
pred1 <- pred1[,-1]
pred2 <- read.csv("predict2.csv")
pred2 <- pred2[,-1]
pred3 <- read.csv("predict3.csv")
pred3 <- pred3[,-1]
predProbs <- read.csv("target1.csv")
predProbs[1] <- predProbs[2]
predProbs[2] <- pred1
predProbs[3] <- pred2
predProbs[4] <- pred3
colnames(predProbs) <- c("Target", "KNN", "Lasso", "SVM")
write.csv(predProbs, "PredictionProb.csv")
