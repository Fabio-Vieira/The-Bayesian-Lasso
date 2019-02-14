###################################################################################
###################################################################################
library(care)

#loading the data
data("efron2004")
dim(efron2004$x) # 442 10
colnames(efron2004$x)
length(efron2004$y) # 442
#head(efron2004)

##################################################################################
##################################################################################
#Full conditional posterior distributions
library(MASS)
library(statmod)

updateBeta <- function(Y,X,Sig2,Tau){
  D <- solve(diag(Tau))
  sig.post <- Sig2 * solve(t(X)%*%X + D)
  m.post <- sig.post %*% t(X) %*% Y
  return(mvrnorm(1,m.post,sig.post))
}

updateSig2 <- function(a,b,Y,X,Beta,Tau){
  N <- length(Y)
  p <- length(Beta)
  D <- solve(diag(Tau))
  a.post <- (N - 1 + p)*0.5 + a
  b.post <- 0.5*(t(Y - X%*%Beta)%*%(Y - X%*%Beta) + t(Beta)%*%D%*%Beta) + b
  return(1/rgamma(1,a.post,b.post))
}

updateTau <- function(Lambda,Sig2,Beta){
  dispersion <- Lambda^2
  mu <- sqrt(Lambda^2 * Sig2 / Beta^2)
  return(1/rinvgauss(length(mu),mean = mu, dispersion = dispersion))
}

updateLambda <- function(p,r,Delta,Tau){#p is length of Beta
  a.post <- p + r
  b.post <- sum(Tau^2)*0.5 + Delta
  return(sqrt(rgamma(1,a.post,b.post)))
}

#################################################################################
#Number of iterations
Niter <- 10000

#Variables
Y <- efron2004$y
X <- efron2004$x

#Creating matrices that will store MCMC samples
Beta.out <- array(NA, dim = c(Niter, dim(efron2004$x)[2]))
Tau.out <- array(NA, dim = c(Niter, dim(efron2004$x)[2]))
Sig.out <- array(NA, dim = Niter)
Lambda.out <- array(NA, dim = Niter)

#Initial Values
Beta.out[1,] <- rep(1,dim(efron2004$x)[2])
Tau.out[1,] <- rep(1,dim(efron2004$x)[2])
Sig.out[1] <- 1
Lambda.out[1] <- 1

#Gibbs sampler
for(i in 2:Niter){
  Beta.out[i,] <- updateBeta(Y,X,Sig.out[i-1],Tau.out[i-1,])
  Tau.out[i,] <- updateTau(Lambda.out[i-1],Sig.out[i-1],Beta.out[i,])
  Sig.out[i] <- updateSig2(0.01,0.01,Y,X,Beta.out[i,],Tau.out[i,])
  Lambda.out[i] <- updateLambda(length(Beta.out[1,]),1,1.78,Tau.out[i,])
  print(i)
}

###################################################################################
#Burn-in phase
Nburn <- 1000
Beta <- Beta.out[-(1:Nburn),]

plot(Beta[,1],type = 'l')
median(Beta[,1])
quantile(Beta[,1],probs=c(0.025,0.975))

plot(Beta[,2],type = 'l')
median(Beta[,2])
quantile(Beta.out[,2],probs=c(0.025,0.975)) #Doesn't include zero

plot(Beta[,3],type = 'l') #Doesn't include zero
median(Beta[,3])
quantile(Beta.out[,3],probs=c(0.025,0.975))

plot(Beta[,4],type = 'l') #Doesn't include zero
median(Beta[,4])
quantile(Beta.out[,4],probs=c(0.025,0.975))

plot(Beta[,5],type = 'l')
median(Beta[,5])
quantile(Beta.out[,5],probs=c(0.025,0.975))

plot(Beta[,6],type = 'l')
median(Beta[,6])
quantile(Beta.out[,6],probs=c(0.025,0.975))

plot(Beta[,7],type = 'l')
median(Beta[,7])
quantile(Beta.out[,7],probs=c(0.025,0.975))

plot(Beta[,8],type = 'l')
median(Beta[,8])
quantile(Beta.out[,8],probs=c(0.025,0.975))

plot(Beta[,9],type = 'l') #Doesn't include zero
median(Beta[,9])
quantile(Beta.out[,9],probs=c(0.025,0.975))

plot(Beta[,10],type = 'l')
median(Beta[,10])
quantile(Beta.out[,10],probs=c(0.025,0.975))

#################################################################################
#Plot similar to figure 2 of the paper, showing the same results were obtained
plot(1, type="n", xlab="", ylab="", xlim=c(-1, 1), ylim=c(11, 1),
     main = "Diabetes data intervals")
abline(v = 0, lty = 2)

for(i in 10:1){
  lines(x = c(quantile(Beta[,i], probs = c(0.025)), 
              quantile(Beta[,i], probs = c(0.975))), y = c(i,i))
  points(median(Beta[,i]), y = i, pch = 10)
}

