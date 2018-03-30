#********************************************************************************
#********************************************************************************
#                                Ensemble Model
#********************************************************************************
#********************************************************************************
#install.packages("devtools")
library("devtools")
#install_github("ecpolley/SuperLearner")
library(SuperLearner)
library(data.table)

predProbs <- read.csv("PredictionProb.csv", header = TRUE)
predProbs <- predProbs[,-1]

#Create Average Ensemble Model
allPredProbs <- predProbs[,-1]
avgsEnsemble <- apply(allPredProbs, 1, function(x) (mean(x)))


#Create Superlearner Ensemble Model
SL.library <- c("SL.knn",
                "SL.glm",
                "SL.svm")
method <- "method.NNLS"
family <- "binomial"

outcome <- predProbs[,1]
allPredProbs <- predProbs[,-1]

fit <- SuperLearner(Y = outcome, X = allPredProbs,
                    family = family,
                    SL.library = SL.library,
                    method = method)

#pred <- predict(fit, newdata = )



#Test the Average Ensemble Model
result_file<-read.table("ViralChallenge_test_Phase1_CLINICAL.tsv",sep="\t",header=TRUE)
data_file <- fread("ViralChallenge_test_Phase1_EXPRESSION_RMA.tsv", data.table=FALSE)
rownames(data_file) <- data_file$FEATUREID
data_file <- data_file[, -1]
final_data_file<-as.data.frame(t(data_file))
final_data_file$SYMPTOMATIC_SC2<-result_file$SYMPTOMATIC_SC2
