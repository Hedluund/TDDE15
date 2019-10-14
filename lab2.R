#######################################
# Jesper Hedlund                      #
# Lab2                                #
# TDDE15 - advanced Machine Learning  #
#######################################
rm(list=ls())
set.seed(12345)

############ Functions ###################

transitionModel <- function(z_state){
  # samples a var to choose from wich ditribution to sample from
  choose_model <- sample(c(1, 2, 3), size = 1)
  
  if(choose_model == 1){
    znext_state <-rnorm(1, mean = z_state, sd=1)
  }else if(choose_model == 2){
    znext_state <-rnorm(1, mean = z_state+1, sd=1)
  }else {
    znext_state <-rnorm(1, mean = z_state+2, sd=1)
  }
 
  return(znext_state)
}

emissionModel <- function(z_state, sd=1) {
  choose_model <- sample(c(1, 2, 3), size = 1)
  
  if(choose_model == 1){
    x_obs <-rnorm(1, mean = z_state, sd=sd)
  }else if(choose_model == 2){
    x_obs <-rnorm(1, mean = z_state-1, sd=sd)
  }else {
    x_obs <-rnorm(1, mean = z_state+1, sd=sd)
  }
  
  return(x_obs)
}

getWeight <- function(observation, particle, sd=1,static = FALSE){
  weight <- (dnorm(observation,particle,sd) +
               dnorm(observation,particle-1,sd) +
               dnorm(observation,particle+1,sd))/3
  #weight <- dnorm(x = particle, mean =(3*obs+1-1)/3 , sd=sqrt((1/3)^2*(1+1+1)))
  
  if(static == FALSE){
    return(weight)
  }else{
    return(rep(1, length(particle)))
  }
  
}


simulate_observations_and_states <- function(init_state, num_sim, sd=1){
  all_states <- rep(0,num_sim)
  observations <- rep(0,num_sim)
  
  curr_state <- init_state
  for( i in 1:num_sim){
    all_states[i] <- curr_state 
    next_state <- transitionModel(curr_state)

    observations[i] <- emissionModel(next_state, sd=sd)
    
    # setting up state for next iteration
    curr_state <- next_state
    
  }
  
  
  return(list("states"=all_states,"observations"=observations))
}


# particle filter function
myParticleFilter <- function(observations, sd=1, T = 100, M=100, static_weights = FALSE){
  
  particle_matrix <- matrix(nrow = T, ncol = T)
  weight_matrix <- matrix(nrow = T, ncol = T)
  
  # create init probabilities
  particles_uniform <- runif(M, 0 ,100)
  
  # calculate of the uniform distributed particles
  weights_init <- getWeight(observations[1], particles_uniform,sd=sd, static_weights)
 
  
  # sample init particles from the weights of the uniforma particles and use as init particles
  particles_init <- sample(particles_uniform, length(particles_uniform), replace = TRUE, prob = weights_init)
  
  particle_matrix[1,] <- particles_init
  weight_matrix[1,] <- weights_init
  
  # renaming init particles and weights to fit the loop
  particles <- particles_init
  weights <- weights_init
  for(t in 2:T){
    
    # setting observation for time t
    observation <- observations[t]
    for(m in 1:M){
      
      particles[m] <- transitionModel(particles[m])
      
      #calculating weights
      m_weigth <- getWeight(observation, particles[m], sd = sd, static= static_weights) 
      weights[m] <- m_weigth
    }
    #sampling new particles given weights from perv particles
    newparticles <- sample(particles, length(particles), replace = TRUE, prob = weights)
    
    #updating
    particles <- newparticles
    
    #saving particles and weights
    particle_matrix[t,] <- particles
    weight_matrix[t,] <- weights
  }
  
  return(list("p_matrix" = particle_matrix, "w_matrix" = weight_matrix))
}



plotObsTrueAndParticle <- function(particles, observation, state, index){
  points(x = rep(index, length(particles)), y = particles)
  points(x = index, y = observation, col="green", pch=19)
  points(x = index, y = state, col=rgb(1,0,0,0.8), pch = 19)
}

################ init variables ##############################


T <- 100
M <- 100
start_state <- runif(1,0, 100)


############## start #########################################

####### task 1 ##########################
states_and_obs <- simulate_observations_and_states(start_state, 100)
simulated_states <- states_and_obs$states
simulated_observations <- states_and_obs$observations

plot(simulated_states, type="o")
points(simulated_observations, type="o", col = rgb(1,0,0,0.3), pch=19)
legend("topleft", c("states", "observations"),col= c("black", "red"), pch= 19)

res_sd1 <- myParticleFilter(simulated_observations, sd = 1)
weights_sd1 <- res_sd1$w_matrix
particles_sd1 <- res_sd1$p_matrix

# plotting variables
particles_to_plot <- c(1, 30, 60, 100)
xGrid <- c(0, 101)
yGrid <- c(0, 200)

# plotting task 1
plot(x=NULL,
     y=NULL,
     xlim = xGrid,
     ylim = yGrid,
     pch=20,
     col=rgb(0.7,0.8,0),
     main=paste("time="),
     xlab = "time",
     ylab = "sd = 1")

for(i in particles_to_plot){
  plotObsTrueAndParticle(particles_sd1[i, ], observation = simulated_observations[i], state = simulated_states[i], index = i)
}


########### Task 2 #############################

states_and_obs_sd5 <- simulate_observations_and_states(start_state, 100, sd=5)
simulated_states_sd5 <- states_and_obs_sd5$states
simulated_observations_sd5 <- states_and_obs_sd5$observations

res_sd5 <- myParticleFilter(simulated_observations_sd5, sd = 5)
weights_sd5 <- res_sd5$w_matrix
particles_sd5 <- res_sd5$p_matrix

plot(x=NULL,
     y=NULL,
     xlim = xGrid,
     ylim = yGrid,
     pch=20,
     col=rgb(0.7,0.8,0),
     main=paste("sd= 5"),
     xlab = "time",
     ylab = "values")

for(i in particles_to_plot){
  plotObsTrueAndParticle(particles_sd5[i, ], observation = simulated_observations_sd5[i], state = simulated_states_sd5[i], index = i)
}

#### sd = 50
states_and_obs_sd50 <- simulate_observations_and_states(start_state, 100, sd=50)
simulated_states_sd50 <- states_and_obs_sd50$states
simulated_observations_sd50 <- states_and_obs_sd50$observations

res_sd50 <- myParticleFilter(simulated_observations_sd50, sd = 50)
weights_sd50 <- res_sd50$w_matrix
particles_sd50 <- res_sd50$p_matrix

plot(x=NULL,
     y=NULL,
     xlim = xGrid,
     ylim = yGrid,
     pch=20,
     col=rgb(0.7,0.8,0),
     main=paste("sd= 50"),
     xlab = "time",
     ylab = "values")

for(i in particles_to_plot){
  plotObsTrueAndParticle(particles_sd50[i, ], observation = simulated_observations_sd50[i], state = simulated_states_sd50[i], index = i)
}

# we can see that the particles gain more spread when we increase the standard deviation. This seems reasonable since we 
# can consider the standard deviations to be how certain we are about the output of the sensors. And with less certainty the spread
# of the particles will increas. Also we can see that they converge towards t=100 but they converge only up to a limit. 
# They do not converge towards zero when the standard deviation is increased.


############ Task 3 ####################


res_static <- myParticleFilter(simulated_observations, sd = 1,static_weights = TRUE)
weights_static <- res_static$w_matrix
particles_static <- res_static$p_matrix

plot(x=NULL,
     y=NULL,
     xlim = xGrid,
     ylim = yGrid,
     pch=20,
     col=rgb(0.7,0.8,0),
     main=paste("weights static"),
     xlab = "time",
     ylab = "values")

for(i in particles_to_plot){
  plotObsTrueAndParticle(particles_static[i, ], observation = simulated_observations[i], state = simulated_states[i], index = i)
}


# The particles converge to a number, but it can be any particle since we sample all particle s with equal probabilities. 
# So if a some particle get sampled several times(we use replacement) the number of that particle will probplably increase. 
# The particles wil end upp clustering some random particle(s)