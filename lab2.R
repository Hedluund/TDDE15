
rm(list = ls())
set.seed(12345)

############## Implemented Functions ############# ====

calc_alpha <- function(state_space, emission_model, transition_model, observation) {
  
  alpha <- NULL
  
  for(state in state_space){
    
    # calculating the probability of the observation given
    # a certain state
    prob_of_state <- emission_model[state, observation]
    
    
  }
}


############## Task 1 ############################# ==== 
# Build a hidden markov model describing a circle of states
library(HMM)

states = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
symbols = c("Xa", "Xb", "Xc", "Xd", "Xe", "Xf", "Xg", "Xh", "Xi", "Xj")


transitionProbMatrix <- t(matrix(data = c(.5,.5,0,0,0,0,0,0,0,0,
                                          0,.5,.5,0,0,0,0,0,0,0,
                                          0,0,.5,.5,0,0,0,0,0,0,
                                          0,0,0,.5,.5,0,0,0,0,0,
                                          0,0,0,0,.5,.5,0,0,0,0,
                                          0,0,0,0,0,.5,.5,0,0,0,
                                          0,0,0,0,0,0,.5,.5,0,0,
                                          0,0,0,0,0,0,0,.5,.5,0,
                                          0,0,0,0,0,0,0,0,.5,.5,
                                          .5,0,0,0,0,0,0,0,0,.5),
                                 nrow = length(states),
                                 ncol = length(symbols)))




emissionProbMatrix <- matrix(data = c(.2,.2,.2,0,0,0,0,0,.2,.2,
                                      .2,.2,.2,.2,0,0,0,0,0,.2,
                                      .2,.2,.2,.2,.2,0,0,0,0,0,
                                      0,.2,.2,.2,.2,.2,0,0,0,0,
                                      0,0,.2,.2,.2,.2,.2,0,0,0,
                                      0,0,0,.2,.2,.2,.2,.2,0,0,
                                      0,0,0,0,.2,.2,.2,.2,.2,0,
                                      0,0,0,0,0,.2,.2,.2,.2,.2,
                                      .2,0,0,0,0,0,.2,.2,.2,.2,
                                      .2,.2,0,0,0,0,0,.2,.2,.2),
                             nrow = length(states),
                             ncol= length(symbols))


circle_hmm <- initHMM(States = states,
                      Symbols = symbols,
                      emissionProbs = emissionProbMatrix,
                      transProbs = transitionProbMatrix)

print(circle_hmm)


############## Task 2 ########################## ====
# simulate the hidden markov model 100 times


res_sim_circle_hmm <- simHMM(circle_hmm, 100)
print(res_sim_circle_hmm)



############## Task 3 and  4 ########  ====
classifyState <- function(probabilityvector){
  possible_states <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
  
  index_of_max <- which.max(probabilityvector)
  
  most_prob_state <- possible_states[index_of_max]
  return(most_prob_state)
  
}

calcAccuracyOfStates <- function(predicted_states, true_states){
  result_predictions <-c()
  
  amount_states_predicted <- length(predicted_states)
  for (i in 1:amount_states_predicted) {
    if(predicted_states[i] == true_states[i]){
      result_predictions <- cbind(result_predictions, 1)
    }else {
      result_predictions <- cbind(result_predictions, 0)
    }
  }
  correct_classifications <- sum(result_predictions)
  
  classificationsrate <- correct_classifications/amount_states_predicted
  
  return(classificationsrate)
}

## computing the smooth distribution with the posterior function
smooth_dist <- posterior(circle_hmm, res_sim_circle_hmm$observation)
# calculating most problable state
smooth_dist_most_prob_states <- apply(X = smooth_dist,MARGIN=2 , FUN= classifyState)
# calculating accuracy
smooth_dist_accuracy <- calcAccuracyOfStates(smooth_dist_most_prob_states, res_sim_circle_hmm$states)

## computing the filtered distribution with the forwards function.Forward gives log probabilities, thus i run exp on the vector.
filtered_dist <- exp(forward(circle_hmm, res_sim_circle_hmm$observation))
# normalizing by using prob.table
filtered_dist_norm <- apply(X = filtered_dist, MARGIN = 2, prop.table)
# calculating most problable state
filtered_dist_most_prob_state <- apply(X = filtered_dist_norm, MARGIN = 2, FUN= classifyState)
#calculating accuracy
filtered_dist_accuracy <- calcAccuracyOfStates(filtered_dist_most_prob_state, res_sim_circle_hmm$states)

## computing the most probable path using the viterbi algorithm
most_prob_path <- viterbi(circle_hmm, res_sim_circle_hmm$observation)
# calculating accuracy
most_prob_path_accuracy <- calcAccuracyOfStates(most_prob_path, res_sim_circle_hmm$states)


############## Task 5 ########  #######################====


res_sim_circle_hmm_2 <- simHMM(circle_hmm, 300)



## computing the smooth distribution with the posterior function
smooth_dist_2 <- posterior(circle_hmm, res_sim_circle_hmm_2$observation)
# calculating most problable state
smooth_dist_most_prob_states_2 <- apply(X = smooth_dist_2,MARGIN=2 , FUN= classifyState)
# calculating accuracy
smooth_dist_accuracy_2 <- calcAccuracyOfStates(smooth_dist_most_prob_states_2, res_sim_circle_hmm_2$states)

print(paste("Smoothing distribution accuracy 100 samples: ", smooth_dist_accuracy))
print(paste("Smoothing distribution accuracy 300 samples: ", smooth_dist_accuracy_2))

## computing the filtered distribution with the forwards function.Forward gives log probabilities, thus i run exp on the vector.
filtered_dist_2 <- exp(forward(circle_hmm, res_sim_circle_hmm_2$observation))
# normalizing by using prob.table
filtered_dist_norm_2 <- apply(X = filtered_dist_2, MARGIN = 2, prop.table)
# calculating most problable state
filtered_dist_most_prob_state_2 <- apply(X = filtered_dist_norm_2, MARGIN = 2, FUN= classifyState)
#calculating accuracy
filtered_dist_accuracy_2 <- calcAccuracyOfStates(filtered_dist_most_prob_state_2, res_sim_circle_hmm_2$states)

print(paste("Filtered distribution accuracy 100 samples: ", filtered_dist_accuracy))
print(paste("Filtered distribution accuracy 300 samples: ", filtered_dist_accuracy_2))

## computing the most probable path using the viterbi algorithm
most_prob_path_2 <- viterbi(circle_hmm, res_sim_circle_hmm_2$observation)
# calculating accuracy
most_prob_path_accuracy_2 <- calcAccuracyOfStates(most_prob_path_2, res_sim_circle_hmm_2$states)

print(paste("Most problable path accuracy 100 samples: ", most_prob_path_accuracy))
print(paste("Most problable path accuracy 300 samples: ", most_prob_path_accuracy_2))

## Same result as in task 4. Smooth distribution performs better than filtered and viterbi most problable path



## Smooth can use the whole dataset of the observations when predicting the state, while the filtering only use 
## the data up to the timepoint of the prediction. That is why it get better predictions. And the viterbi
## has the condition of having to follow a viable path which tthe smooth and the filtering does not have to. 
## That is why that most often performs the worst.




############## Task 6 ######## ====

library(entropy)
entropy_filtered <- apply(filtered_dist, 2,entropy.empirical)
plot(entropy_filtered, type = "o", main = "Filtered distribution entropy. 100 samples")
entropy_filtered_2 <- apply(filtered_dist_2, 2,entropy.empirical)
plot(entropy_filtered_2, type = "o", main = "Filtered distribution entropy. 100 samples")

entropy_smooth <- apply(smooth_dist, 2,entropy.empirical)
plot(entropy_smooth, type = "o", main = "Smoothed distribution entropy. 100 samples")
entropy_smooth_2 <- apply(smooth_dist_2, 2,entropy.empirical)
plot(entropy_smooth_2, type = "o", main = "Filtered distribution entropy. 300 samples")



############## Task 7 ######## ==== 
last_sensor_input <- filtered_dist_norm[,100]

prob_101 <- transitionProbMatrix %*% last_sensor_input

print("The probabilies for the hidden state at time step 101 ")
print( prob_101)