library(bnlearn)

library(gRain)
library(caret)
library(expss)
set.seed(12345)
rm(list = ls())


data <- data("asia")

data <- asia

###### Task 1 ####################
# Function to give info about two BN:s and how "similar" they are
checkInfoBns <- function(bn1, bn2){
  plot(bn1, main="BN 1")
  plot(bn2, main="BN 2")
  
  print(vstructs(bn1))
  print(vstructs(bn2))
  
  classbn1 <- cpdag(bn1)
  classbn2 <- cpdag(bn2)
  plot(classbn1, main ="cpdag bn1")
  plot(classbn2, main ="cpdag bn2")
  
  print(all.equal(bn1,bn2))
  print(all.equal(classbn1, classbn2))
}


# Trying hill climb with different iss.

bn1 <- hc(data,iss=5, score="bde")

bn2 <- hc(data, iss=2, score="bde")

par(mfrow = c(2,2))

checkInfoBns(bn1,bn2)
## Resulted in different BN:s. Different iss gives a different weight of the prior to the score.

# Trying hillclimb with different type of scores
bn1 <- hc(data, score = "bic")
bn2 <- hc(data, score = "bde")

checkInfoBns(bn1,bn2)

## Resulted in similar networks


# Checking with different random restarts


bn1 <- hc(data, restart = 0)
bn2 <- hc(data, restart = 100)

checkInfoBns(bn1,bn2)

## Resulted in Similar networks

# Trying with different start structure
bnstart1 <- model2network("[A|B:D][S|T:L:B][T][L][B|D][E][X][D]")
bnstart2 <- model2network("[A][S][T][L][B][E][X][D]")

bn1 <- hc(data, start = bnstart1, score="bde", iss =2)
bn2 <- hc(data, start = bnstart2, score="bde", iss =2)

checkInfoBns(bn1,bn2)

## Different start bns resulted in different end BN:s. Problably because the algorithm will start in i different state and thus end up in different
## local maximum/minimum.
### NOTE: Can use iss as argument if start is used. Why? Answer: iss can only be used with certain kind of scores eg. bde. Bic can't use 
### iss since it does not match with dirchlet distribution(?)

####### Task 2 #########################################


missclass=function(X,X1) {
  n=length(X)
  return(1-sum(diag(table(X,X1)))/n)
}

my_predict_bn <- function(grain_obj, data, predictors, to_predict){
  
  predictions <- rep(NA,nrow(data))
  # Compile gets the data ready. Prepares data => faster algorithm. 
  compiled_bn <- compile(grain_obj)
  
  for(row in 1:nrow(data)){
    z<-NULL
    for(p in predictors){
      if(data[row, p]=="no"){
        z<-c(z,"no")
      }
      else{
        z<-c(z,"yes")
      }
    }
    # Setting evidence that the states are actual observations.
    curr_evidence <- setEvidence(object = compiled_bn, nodes = predictors, states = z)
    
    # Gives conditional posterior distribution on the node to predict 
    conditional_dist <- querygrain(curr_evidence, nodes=to_predict)
    
    predictions[row] <-if(conditional_dist$S["yes"] >= 0.50) "yes" else "no"
    
  }
  return(predictions)
}

# Setting sample size
sample_size <- floor(0.8* nrow(data))

train_indices <- sample(seq_len(nrow(data)), sample_size)

# creating train and test
train <- data[train_indices,]
test <- data[-train_indices,]

# training network with tabu algorithm
bn_trained <- hc(train)
bn_true_model <- model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")

# Ploting the different bns
par(mfrow = c(1,2))
plot(bn_trained, main = "Trained Model")
plot(bn_true_model, main = "True Model")

#Fitting the network
bn_trained_fit <- bn.fit(bn_trained, data = train)
bn_true_dist_fit <- bn.fit(bn_true_model, data = train)

# Converting to "grain" object
bn_trained_grain <- as.grain(bn_trained_fit)
bn_true_grain <- as.grain(bn_true_dist_fit)



predictors <- c("A","D", "X", "E", "B", "L", "T")
to_predict <- c("S")

res_bn_trained <- my_predict_bn(bn_trained_grain,
                                data = test, 
                                predictors = predictors,
                                to_predict = to_predict)
res_bn_true <- my_predict_bn(bn_true_grain,
                             data = test,
                             predictors = predictors,
                             to_predict = to_predict)

trained_confm <- cro(res_bn_trained, test$S)
true_confm <- cro(res_bn_true, test$S)

print(trained_confm)
print(true_confm)

print(missclass(res_bn_trained, test$S))

########### Task 3 #################
# Predict only on parents and children of node S
mb_trained = mb(bn_trained_fit, to_predict)
mb_true = mb(bn_true_dist_fit, to_predict)

res_mb_bn_trained <- my_predict_bn(bn_trained_grain,
                                   data = test, 
                                   predictors = mb_trained,
                                   to_predict = to_predict)
res_mb_bn_true <- my_predict_bn(bn_true_grain,
                                data = test,
                                predictors = mb_true,
                                to_predict = to_predict)

trained_mb_confm <- cro(res_mb_bn_trained, test$S)
true_mb_confm <- cro(res_mb_bn_true, test$S)

print(trained_mb_confm)
print(true_mb_confm)

# Markov blanket is enough for the predictions. The predictor is independent of all Nodes that are not
# in MB given MB

########## Task 4 ###################
# Condition: All features are independet given the predictor. Tail to tail.
# https://towardsdatascience.com/introduction-to-naive-bayes-classification-4cffabb1ae54

bn_naive <- model2network("[S][A|S][T|S][L|S][B|S][D|S][E|S][X|S]")
plot(bn_naive)

bn_naive_fit <- bn.fit(bn_naive, data = train)

bn_naive_grain <- as.grain(bn_naive_fit)

bn_naive_res <- my_predict_bn(bn_naive_grain, data = test, predictors = predictors, to_predict = to_predict)

naive_confm <- cro(bn_naive_res, test$S)
print(naive_confm)
missclass_naive <- missclass(bn_naive_res, test$S)
print(true_confm)

