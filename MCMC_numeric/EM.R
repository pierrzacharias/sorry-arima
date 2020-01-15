library(coda)
data(faithful)
X <- faithful$eruptions
hist(X,breaks=100,freq=FALSE)

mu_sigma_given <- function(X,Z,i,pi_1){
  
  keep <- (c(1:length(Z)) != i)  # toutes les coordonnés sauf i
  sample1 <- X[keep & (Z == 1)]  # données associées à la classe 1 sans la variable i
  sample2 <- X[keep & (Z == 2)]  # données associées à la classe 1 sans la variable i
  
  mu_1 <- mean(sample1)
  sigma_1 <- sd(sample1)
  
  mu_2 <- mean(sample2)
  sigma_2 <- sd(sample2)
  
  proba1 <- pi_1 * dnorm(X[i],mu_1,sigma_1) / ( pi_1 * dnorm(X[i],mu_1,sigma_1) + (1-pi_1) * dnorm(X[i],mu_2,sigma_2))
  
  if ( runif(1) < proba1 ){
        Z[i]=1
    } else {
      Z[i]=2
      }
  return(Z)
}

mixturegibbs <- function(X,N_simu){
  
  # initalisation
  Z <- sample(c(1,2),length(X),replace=TRUE) # matrice assignations
  mu_1 <- mean(X[Z == 1])
  sigma_1 <- sd(X[Z == 1])
  mu_2 <- mean(X[Z == 2])
  sigma_2 <- sd(X[Z == 2])
  
  pi_1 <- sum(Z==1)/length(X)
  
  
  # vecteur pour stocker les itérations
  Z_matrix <- Z
  mu_1_vect <- mu_1
  sigma_1_vect <- sigma_1
  mu_2_vect <- mu_2
  sigma_2_vect <- sigma_2
  pi_1_vect <- pi_1
  
  for (k in 1:N_simu){ # itérations
    
    for (i in 1:length(Z)){
      Z <- singleupdate(X,Z,i,pi_1)
    }
    
    Z_matrix <- rbind(Z_matrix,Z)
    
    mu_1 <- mean(X[Z==1])
    mu_1_vect <- c(mu_1_vect,mu_1)
    
    sigma_1 <- sd(X[Z==1])
    sigma_1_vect <- c(sigma_1_vect,sigma_1)
    
    mu_2 <- mean(X[Z==2])
    mu_2_vect <- c(mu_2_vect,mu_2)
    
    sigma_2 <- sd(X[Z==2])
    sigma_2_vect <- c(sigma_2_vect,sigma_2)
   
    pi_1 <- sum(Z==1)/length(X)
    pi_1_vect <- c(pi_1_vect, pi_1)
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

res <- mixturegibbs(X,200)

Z <- as.mcmc(res$Z)
mu_1 <- as.mcmc(res$mu_1)
mu_2 <- as.mcmc(res$mu_2) 
sigma_1 <- as.mcmc(res$sigma_1)
sigma_2 <- as.mcmc(res$sigma_2)

pi_1_vect <- as.mcmc(res$pi_1)
plot(mu_1)
   
mu_1est <- mean(mu_1[50:200])
mu_2est <- mean(mu_2[50:200])
sigma_1est <- mean(sigma_1[50:200])
sigma_2est <- mean(sigma_2[50:200])
pi_1est <- mean(pi_1_vect[50:200])
hist(X,breaks=100,freq=FALSE)
curve(pi_1est*dnorm(x,mu_1est,sigma_1est),add=TRUE,col='red')
curve((1-pi_1est)*dnorm(x,mu_2est,sigma_2est),add=TRUE,col='blue')
   
   
   
   
   
   
   
   
   
   
   
   
   