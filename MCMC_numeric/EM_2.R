##########################################################
# Gaussian mixture simmulation using Gibbs sampling
#########################################################
set.seed(150)

library(coda)

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
    
    mu_prior = list(mean = 0, precision = 0.5)
    
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
X = rmix(n=1000,pi=c(0.5,0.5),mu=c(-2,3),s=c(1,sqrt(0.5)))

png(filename="/home/pierre/Documents/Courses/project_serie_temporelles/time_series_project/MCMC_numeric/simu_gaussian/data.png",
    width=600, height=350)
hist(X,breaks=100,freq=FALSE)
dev.off()
# other data
#data(faithful)
#X <- faithful$eruptions
#hist(X,breaks=100,freq=FALSE)

N_simu = 1000

res <- echantillonneur_gibbs(X,N_simu)

plot(res$mu_1,ylim=c(-4,4),type="l",col="red",xlab = "nombre de simulations", 
     ylab = "tirages")
lines(res$mu_2,col="blue")
legend("bottomleft", 
       legend = c(expression(paste("", mu, 1)), expression(paste("", mu, 2))), 
       col = c("red", 
               "blue"), 
       pch = c(0,0), 
       bty = "n", 
       pt.cex = 0, 
       cex = 2, 
       text.col =  c("red", 
                     "blue"), 
       horiz = F , 
       # inset = c(0.1, 0.1)
       )

legend(1:length(res$mu_1),legend=c(expression(paste("estimation de la densité de ", mu,1)), expression(paste("Phase Angle ", mu_1))), col=c("red", "blue"))


burn = 100

Z <- as.mcmc(res$Z)
mu_1 <- as.mcmc(res$mu_1) 
mu_2 <- as.mcmc(res$mu_2) 
sigma_1 <- as.mcmc(res$sigma_1)
sigma_2 <- as.mcmc(res$sigma_2)

plot(as.mcmc(res$mu_1[-(1:burn)]))


####### trace 
png(filename="/home/pierre/Documents/Courses/project_serie_temporelles/time_series_project/MCMC_numeric/simu_gaussian/mu_plot_dens.png",
      width=600, height=350)
densplot(as.mcmc(res$mu_1[-(1:burn)]),main = expression(paste("estimation de la densité de ", mu,1)))
dev.off()

png(filename="/home/pierre/Documents/Courses/project_serie_temporelles/time_series_project/MCMC_numeric/simu_gaussian/mu_plot_trac.png",
    width=600, height=350)
traceplot(as.mcmc(res$mu_1[-(1:burn)]),main = expression(paste("trace de ", mu,1)))
dev.off()

#### estimateur 
mu_1_mean <- mean(as.mcmc(res$mu_1[-(1:burn)])) 
mu_1_mean
mu_2_mean <- mean(as.mcmc(res$mu_2[-(1:burn)]))
mu_2_mean

HPDinterval(as.mcmc(res$mu_1[-(1:burn)]), 0.95)

HPDinterval(as.mcmc(res$mu_2[-(1:burn)]), 0.95)

gelman.plot(as.mcmc(res$mu_1[100:400]),as.mcmc(res$mu_1[500:800]))

pi_1_vect <- as.mcmc(res$pi_1)
plot(as.mcmc(res$mu_1[-(1:burn)]))
mu_1 <- as.mcmc(res$mu_1) 
mean(mu_1)


sigma_1_mean <- mean(sigma_1[-(1:burn)])
sigma_2_mean <- mean(sigma_2[-(1:burn)])
pi_1_mean <- mean(pi_1_vect[-(1:burn)])
hist(X,breaks=100,freq=FALSE)
curve(pi_1_mean*dnorm(x,mu_1_mean,sigma_1_mean),add=TRUE,col='red')
curve((1-pi_1_mean)*dnorm(x,mu_2_mean,sigma_2_mean),add=TRUE,col='blue')

################### intervale de confiances
plot(quantile(res$mu_1[-(1:burn)],c(0.05,0.95)))

1.96 * sd(res$mu_1[-(1:burn)]) / sqrt(N_simu)

boxplot(res$mu_1[-(1:burn)])










