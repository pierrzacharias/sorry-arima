##########################################################
# Gaussian mixture simmulation using Gibbs sampling
#########################################################
set.seed(150)

library(coda)
library(bquote)

  Z_given_mu <- function(X,Z,mu,pi_1){
    # Z variables latentes
    
    remove_i <- (c(1:length(X)) !=i) # on enleve la   coordonnee  i
    estimate_sigma_1 <- sd(X[(remove_i & Z == 1)]) # classe 1
    estimate_sigma_2 <- sd(X[(remove_i & Z == 2)]) # classe 2 
                  
    for (i in 1:length(Z)){
                  
        proba1 <- pi_1 * dnorm(X[i],mu[1],estimate_sigma_1) / (pi_1 * dnorm(X[i],mu[1],estimate_sigma_1) + (1-pi_1) * dnorm(X[i],mu[2],estimate_sigma_2))
        Z[i] = sample(1:2, size=1,prob=c(proba1, 1-proba1),replace=TRUE)
    }
    return(Z)
  }

  mu_given_Z = function(X, Z, mu_prior){
    # Z variables latentes  
    # mu_prior contient parametres loi a priori de mu
    mu = rep(0,2)
    sigma = rep(0,2)
      
    for(j in 1:2){
        
      sample_j_size = sum(Z==j)
      sample_j_mean = mean(X[Z==j])
      sigma[j] = sd(X[Z==j])
      precision_j = 1 / sigma[j]^2
        
      precision_post = sample_j_size * precision_j + mu_prior$precision
      mean_post = (sample_j_mean * sample_j_size * precision_j + mu_prior$mean * mu_prior$precision ) / precision_post
        
      mu[j] = rnorm(1,mean_post,sqrt(1/precision_post)) # on tire mu selon la loi normale a posteriori
    }
    return(list(mu = mu, sigma = sigma))
  }

  echantillonneur_gibbs <- function(X,N_simu){
    
    # initalisation
    Z <- sample(c(1,2),length(X),replace=TRUE) # matrice des assignations
    
    mu = rep(0,2)
    sigma = rep(0,2)
    for(j in 1:2){
      mu[j] <- mean(X[Z == j])
      sigma[j] <- sd(X[Z == j])
    }
    pi_1 <- sum(Z==1)/length(X)
    pi_1_vect <- pi_1
    
    mu_prior = list(mean = 0, precision = 0.1)
    
    # vecteur pour stocker les iterations
    Z_matrix <- Z
    mu_1_vect <- mu[1]
    sigma_1_vect <- sigma[1]
    mu_2_vect <- mu[2]
    sigma_2_vect <- sigma[2]
    
    
    for (k in 1:N_simu){ # iterations
      
      Z <- Z_given_mu(X,Z,mu,pi_1)
      
      
      Z_matrix <- rbind(Z_matrix,Z)
      
      pi_1 <- sum(Z==1)/length(X)
      pi_1_vect <- c(pi_1_vect, pi_1)
      
      param_post <- mu_given_Z(X, Z, mu_prior)
      
      mu = param_post$mu
      sigma = param_post$sigma
      
      mu_1_vect <- c(mu_1_vect,mu[1])
      mu_2_vect <- c(mu_2_vect,mu[2])
      sigma_1_vect <- c(sigma_1_vect,sigma[1])
      sigma_2_vect <- c(sigma_2_vect,sigma[2])
      r
    }
    return(list(
      Z = Z_matrix,
      mu_1 = mu_1_vect,
      sigma_1 = sigma_1_vect,
      mu_2 = mu_2_vect,
      sigma_2 = sigma_2_vect,
      pi_1 = pi_1_vect
    ))
  }

#############################################################################
#             Test
###########################################################################


rmix = function(n,pi,mu,s){
  # generate from mixture of normals
  # n number of samples
  # pi mixture proportions
  # mu mixture means
  # s mixture standard deviations
  z = sample(1:length(pi),prob=pi,size=n,replace=TRUE)
  x = rnorm(n,mu[z],s[z])
  return(x)
}
X = rmix(n=1000,pi=c(0.5,0.5),mu=c(-2,2),s=c(1,1))
hist(X,breaks=100,freq=FALSE)

# other data
#data(faithful)
#X <- faithful$eruptions
#hist(X,breaks=100,freq=FALSE)



res <- echantillonneur_gibbs(X,2000)

plot(res$mu_1,ylim=c(-4,4),type="l",col="red",xlab = "nombre de simulations", 
     ylab = expression(paste("", mu, 1)))
lines(res$mu_2,col="blue")
legend(res$mu_1,legend=c(expression(paste("Phase Angle ", mu_1)), expression(paste("Phase Angle ", mu_1))), col=c("red", "blue"))


burn = 100

Z <- as.mcmc(res$Z)
mu_1 <- as.mcmc(res$mu_1[-(1:burn)])
mu_2 <- as.mcmc(res$mu_2) 
sigma_1 <- as.mcmc(res$sigma_1)
sigma_2 <- as.mcmc(res$sigma_2)

pi_1_vect <- as.mcmc(res$pi_1)
plot(mu_1)
mean(mu_1)
mu_1est <- mean(as.mcmc(res$mu_1[-(1:burn)]))
mu_1est

mu_2est <- mean(mu_2[50:200])
sigma_1est <- mean(sigma_1[50:200])
sigma_2est <- mean(sigma_2[50:200])
pi_1est <- mean(pi_1_vect[50:200])
hist(X,breaks=100,freq=FALSE)
curve(pi_1est*dnorm(x,mu_1est,sigma_1est),add=TRUE,col='red')
curve((1-pi_1est)*dnorm(x,mu_2est,sigma_2est),add=TRUE,col='blue')

################### intervale de confiances
quantile(res$mu[-(1:10)],c(0.05,0.95))










