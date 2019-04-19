### 1. Gausian Mixture Model (1-dim, 2-class)

#1) stopping rule에 likelihood 이용

MyGMM1 <- function(X, threshold=10^(-10), max.iter=1000){
  
  pi <- 1/2; mu <- sort(sample(X, 2)); sigma <- rep(var(X), 2)
  err <- 1; iter <- 0; lik <- 0; n <- length(X)
  
  while(err > threshold & iter < max.iter ){
    
    Y <- pi*dnorm(X, mu[1], sigma[1]) / (pi*dnorm(X, mu[1], sigma[1])+(1-pi)*dnorm(X, mu[2], sigma[2]))
    pi <- sum(Y)/n
    mu <- c(sum(Y*X)/sum(Y), sum((1-Y)*X)/sum(1-Y))
    sigma <- sqrt(c(sum(Y*(X-mu[1])^2)/sum(Y), sum((1-Y)*(Xmu[2])^2)/sum(1-Y)))
    
    old.lik <- lik
    lik <- sum(Y*log(dnorm(X, mu[1], sigma[1])) + (1-Y)*log(dnorm(X, mu[2], sigma[2])))
    err <- abs(lik-old.lik)
    iter <- iter + 1
    
  }
  
  return(list(Y=ifelse(Y>0.5, 1, 0), theta=c(pi=pi, mu=mu, sigma=sigma), niter=iter))
  
}


#2) stopping rule에 pi값 이용

MyGMM2 <- function(X, threshold=10^(-10), max.iter=1000) {
  
  pi <- 1/2; mu <- sort(sample(X, 2)); sigma <- rep(var(X), 2)
  err <- 1; iter <- 0; n <- length(X)
  
  while(err > threshold & iter < max.iter ){
    
    Y <- pi*dnorm(X, mu[1], sigma[1])/(pi*dnorm(X, mu[1], sigma[1])+(1-pi)*dnorm(X, mu[2], sigma[2]))
    old.pi <- pi
    pi <- sum(Y)/n
    mu <- c(sum(Y*X)/sum(Y), sum((1-Y)*X)/sum(1-Y))
    sigma <- sqrt(c(sum(Y*(X-mu[1])^2)/sum(Y), sum((1-Y)*(Xmu[2])^2)/sum(1-Y)))
    err <- abs(pi-old.pi); iter <- iter + 1
    
  }
  
  return(list(Y=ifelse(Y>0.5, 1, 0), theta=c(pi=pi, mu=mu, sigma=sigma)))
  
}


X <- c(X1 = rnorm(20000, -5, 3), X2 = rnorm(10000, 5, 2))
Class <- as.factor(rep(c(1, 2), c(20000,10000)))

library(ggplot2)
ggplot() + geom_density(aes(X))
ggplot() + geom_density(aes(X, linetype=Class)) + scale_linetype_manual(values=c(1, 2))

GMM1 <- MyGMM1(X); GMM1[c(2,3)]
GMM2 <- MyGMM2(X); GMM2[c(2,3)]





### 2. Bivariate Normal with Missing Values

W1 <- c(8, 11, 16, 18, 6, 4, 20, 25, 9, 13)
W2 <- c(10, 14, 16, 15, 20, 4, 18, 22, NA, NA)

n <- length(W1)

# initial value 설정
W2[c(9,10)] <- mean(W2, na.rm=T)
T <- c(sum(W1), sum(W2)); mu <- T/n
sigma <- matrix(c(sum(W1^2)-T[1]^2/n, rep(sum(W1*W2)-T[1]*T[2]/n, 2), sum(W2^2)-T[2]^2/n)/n, nrow=2)
rho <- 0.5

err <- 1; iter <- 0

while(err > 10^(-10) & iter < 100) {
  
  W2[9:10] <- mu[2]+(sigma[1,2]/sigma[1,1])*(W1[9:10]-mu[1])
  W2sq <- c(W2[1:8]^2, W2[9:10]^2 + sigma[2,2]*(1-rho^2))
  
  old.mu <- mu; old.sigma <- sigma
  T <- c(sum(W1), sum(W2)); mu <- T/n
  sigma <- matrix(c(sum(W1^2)-T[1]^2/n, rep(sum(W1*W2)-T[1]*T[2]/n, 2), sum(W2sq)-T[2]^2/n)/n, nrow=2)
  rho <- sigma[1,2]/sqrt(sigma[1,1]*sigma[2,2])
  
  err <- sum(abs(old.mu-mu), abs(old.sigma-sigma))
  iter <- iter + 1

}

list(mu=mu, sigma=sigma)
iter





### 3. Missing Values in Response Variable

data <- read.csv("Movie.csv"); dim(data)
data <- data[complete.cases(data), -c(1, 2)]; dim(data)

original.Y <- data$Gross
n <- nrow(data)
p <- seq(0.01, 0.1, 0.01) # p : missing proportion
missing <- list(); for(i in 1:10) missing[[i]] <- sample(n, n*p[i])


## Fit Y Using Linear Regression

ggplot(data) + geom_histogram(aes(Gross), color="grey")
ggplot(data) + geom_histogram(aes(Gross^(1/3)), color="grey")
ggplot(data) + geom_qq(aes(sample=Gross^(1/3)))


result.lm <- data.frame(p=p, error=0, n.iter=0, MAE=0)

for(i in 1:10){
  
  data$Gross <- original.Y
  data$Gross[missing[[i]]] <- mean(original.Y[-missing[[i]]])
  
  err <- 1; iter <- 0
  
  while(err > 10^(-10) & iter < 1000){
    
    lm.step <- step(lm((Gross)^(1/3) ~ ., data))
    old.Y <- data$Gross[missing[[i]]]
    data$Gross[missing[[i]]] <- fitted(lm.step)[missing[[i]]]^3
    
    err <- sum(abs(old.Y-data$Gross[missing[[i]]]))
    iter <- iter + 1

  }
  
  result.lm[i, -1] <- c(err, iter, mean(abs(original.Y-data$Gross)))
  
}

result.lm


## Fit Y Using Random Forest

library(randomForest)
result.rf <- data.frame(p=p, error=0, n.iter=0, MAE=0)

for(i in 1:10){
  
  data$Gross <- original.Y
  data$Gross[missing[[i]]] <- mean(original.Y[-missing[[i]]])
  
  err <- 1; iter <- 0
  
  while(err > 10^(-10) & iter < 1000){
    
    rf <- randomForest(Gross ~ ., data)
    old.Y <- data$Gross[missing[[i]]]
    data$Gross[missing[[i]]] <- predict(rf)[missing[[i]]]
    
    err <- sum(abs(old.Y-data$Gross[missing[[i]]]))
    iter <- iter + 1
    
  }
  
  result.rf[i, -1] <- c(err, iter, mean(abs(original.Y-data$Gross)))
  
}

result.rf





### 4. Multinomial with Complex Cell Structure

nO <- 176; nA <- 182; nB <- 60; nAB <- 17; n <- 435
p <- (nA*(3/4)+nAB/2)/n; q <- (nB*(3/4)+nAB/2)/n; r <- 1-p-q

err <- 1; iter <- 0

while(err > 10^(-10) & iter < 100){
  
  nAA <- nA*p^2/(p^2+2*p*r); nAO <- nA*2*p*r/(p^2+2*p*r)
  nBB <- nB*q^2/(q^2+2*q*r); nBO <- nB*2*q*r/(q^2+2*q*r)
  nA2 <- nAA+nAO/2+nAB/2; nB2 <- nBB+nBO/2+nAB/2; nO2 <- nO+nAO/2+nBO/2
  
  old.p <- p; old.q <- q
  p <- nA2/n; q <- nB2/n; r <- 1-p-q
  
  err <- sum(abs(p-old.p)+abs(q-old.q))
  iter <- iter + 1
  
}

c(p=p, q=q)
iter
