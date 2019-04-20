##### 3-1. Bivariate Normal with Missing Values

W1 <- c(8, 11, 16, 18, 6, 4, 20, 25, 9, 13)
W2 <- c(10, 14, 16, 15, 20, 4, 18, 22)
W2[9:10] <- mean(W2)

mu1 <- mean(W1)
mu2 <- mean(W2)
S <- matrix(c(sum(W1^2)-sum(W1)^2/10, rep(sum(W1*W2)-sum(W1)*sum(W2)/10, 2), sum(W2^2)-sum(W2)^2/10)/10, nrow=2)


### Initialization

n <- 30; m <- 100

for(i in 1:n){
  
  w9 <- c(); w10 <- c()
  
  for(t in 1:m){
    
    w9[t] <- rnorm(1, mu2+(S[1,2]/S[1,1])*(W1[9]-mu1), sqrt(S[2,2]-S[1,2]^2/S[1,1]))
    W2[9] <- w9[t]; mu2 <- mean(W2)
    S <- matrix(c(sum(W1^2)-sum(W1)^2/10, rep(sum(W1*W2)-sum(W1)*sum(W2)/10, 2), sum(W2^2)-sum(W2)^2/10)/10, nrow=2)
    
    w10[t] <- rnorm(1, mu2+(S[1,2]/S[1,1])*(W1[10]-mu1), sqrt(S[2,2]-S[1,2]^2/S[1,1]) )
    W2[10] <- w10[t]; mu2 <- mean(W2)
    S <- matrix(c(sum(W1^2)-sum(W1)^2/10, rep(sum(W1*W2)-sum(W1)*sum(W2)/10, 2), sum(W2^2)-sum(W2)^2/10)/10, nrow=2)
    
  }
  
  W2[9:10] <- c(mean(w9), mean(w10)); mu2 <- mean(W2)
  S <- matrix(c(sum(W1^2)-sum(W1)^2/10, rep(sum(W1*W2)-sum(W1)*sum(W2)/10, 2), sum(W2^2)-sum(W2)^2/10)/10, nrow=2)

}

mu20 <- mu2; S0 <- S


### Iterations

m <- 10^5; w9 <- c(); w10 <- c()

for(t in 1:m){
  w9[t] <- rnorm(1, mu2+(S[1,2]/S[1,1])*(W1[9]-mu1), sqrt(S[2,2]-S[1,2]^2/S[1,1]))
  W2[9] <- w9[t]; mu2 <- mean(W2)
  S <- matrix(c(sum(W1^2)-sum(W1)^2/10, rep(sum(W1*W2)-sum(W1)*sum(W2)/10, 2), sum(W2^2)-sum(W2)^2/10)/10, nrow=2)

  w10[t] <- rnorm(1, mu2+(S[1,2]/S[1,1])*(W1[10]-mu1), sqrt(S[2,2]-S[1,2]^2/S[1,1]) )
  W2[10] <- w10[t]; mu2 <- mean(W2)
  S <- matrix(c(sum(W1^2)-sum(W1)^2/10, rep(sum(W1*W2)-sum(W1)*sum(W2)/10, 2), sum(W2^2)-sum(W2)^2/10)/10, nrow=2)

}

W2[9:10] <- c(mean(w9), mean(w10)); mu2 <- mean(W2)
S <- matrix(c(sum(W1^2)-sum(W1)^2/10, rep(sum(W1*W2)-sum(W1)*sum(W2)/10, 2), sum(W2^2)-sum(W2)^2/10)/10, nrow=2)

L <- function(x, mu2, S){
  prod(dnorm(x, mu2+(S[1,2]/S[1,1])*(W1[9:10]-mu1), sqrt(S[2,2]-S[1,2]^2/S[1,1])))
}

r <- 0; err <- 1; s12 <- S[1,2]; s22 <- S[2,2]

system.time(
  
  while(err > 10^(-10)){
    
    r <- r+1
    w <- apply(cbind(w9, w10), 1, function(x){ L(x, mu2[r], S) }) / apply(cbind(w9, w10), 1, function(x){ L(x, mu20, S0) })
    W2[9:10] <- c(sum(w9*w)/sum(w), sum(w10*w)/sum(w))
    mu2[r+1] <- mean(W2); old.S <- S
    S <- matrix(c(sum(W1^2)-sum(W1)^2/10, rep(sum(W1*W2)-sum(W1)*sum(W2)/10, 2), sum(W2^2)-sum(W2)^2/10)/10, nrow=2)
    s12[r+1] <- S[1,2]; s22[r+1] <- S[2,2]
    err <- sum(abs(mu2[r+1]-mu2[r]), abs(S-old.S))
    
  }
)

r; mu2; s12; s22
mu2.IS <- mu2; s12.IS <- s12; S22.IS <- s22


### No Important Sampling

W1 <- c(8, 11, 16, 18, 6, 4, 20, 25, 9, 13)
W2 <- c(10, 14, 16, 15, 20, 4, 18, 22)
W2[9:10] <- mean(W2)

mu1 <- mean(W1); mu2 <- mean(W2)
S <- matrix(c(sum(W1^2)-sum(W1)^2/10, rep(sum(W1*W2)-sum(W1)*sum(W2)/10, 2), sum(W2^2)-sum(W2)^2/10)/10, nrow=2)

## Initialization

n <- 16; m <- 100

for(i in 1:n){
  
  w9 <- c(); w10 <- c()
  
  for(t in 1:m){
    
    w9[t] <- rnorm(1, mu2+(S[1,2]/S[1,1])*(W1[9]-mu1), sqrt(S[2,2]-S[1,2]^2/S[1,1]))
    W2[9] <- w9[t]; mu2 <- mean(W2)
    S <- matrix(c(sum(W1^2)-sum(W1)^2/10, rep(sum(W1*W2)-sum(W1)*sum(W2)/10, 2), sum(W2^2)-sum(W2)^2/10)/10, nrow=2)

    w10[t] <- rnorm(1, mu2+(S[1,2]/S[1,1])*(W1[10]-mu1), sqrt(S[2,2]-S[1,2]^2/S[1,1]) )
    W2[10] <- w10[t]
    S <- matrix(c(sum(W1^2)-sum(W1)^2/10, rep(sum(W1*W2)-sum(W1)*sum(W2)/10, 2), sum(W2^2)-sum(W2)^2/10)/10, nrow=2)

  }
  
  W2[9:10] <- c(mean(w9), mean(w10)); mu2 <- mean(W2)
  S <- matrix(c(sum(W1^2)-sum(W1)^2/10, rep(sum(W1*W2)-sum(W1)*sum(W2)/10,b2), sum(W2^2)-sum(W2)^2/10)/10, nrow=2)
}

## Iterations

r <- 0; err <- 1; m <- 10^5; s12 <- S[1,2]; s22 <- S[2,2]

system.time(
  
  while(r < 30){
    
    r <- r+1
    
    for(t in 1:m){
      
      w9[t] <- rnorm(1, mu2[r]+(S[1,2]/S[1,1])*(W1[9]-mu1), sqrt(S[2,2]-S[1,2]^2/S[1,1]))
      W2[9] <- w9[t]; mu2[r] <- mean(W2)
      S <- matrix(c(sum(W1^2)-sum(W1)^2/10, rep(sum(W1*W2)-sum(W1)*sum(W2)/10, 2), sum(W2^2)-sum(W2)^2/10)/10, nrow=2)

      w10[t] <- rnorm(1, mu2[r]+(S[1,2]/S[1,1])*(W1[10]-mu1), sqrt(S[2,2]-S[1,2]^2/S[1,1]) )
      W2[10] <- w10[t]; mu2[r] <- mean(W2)
      S <- matrix(c(sum(W1^2)-sum(W1)^2/10, rep(sum(W1*W2)-sum(W1)*sum(W2)/10, 2), sum(W2^2)-sum(W2)^2/10)/10, nrow=2)

    }
    
    W2[9:10] <- c(mean(w9), mean(w10))
    mu2[r+1] <- mean(W2); old.S <- S
    S <- matrix(c(sum(W1^2)-sum(W1)^2/10, rep(sum(W1*W2)-sum(W1)*sum(W2)/10, 2), sum(W2^2)-sum(W2)^2/10)/10, nrow=2)
    s12[r+1] <- S[1,2]; s22[r+1] <- S[2,2]
    err <- sum(abs(mu2[r+1]-mu2[r]), abs(S-old.S))
    
  }
)

mu2; s12; s22
mu2.MC <- mu2; s12.MC <- s12; s22.MC <- s22


### Comparison

library(ggplot2)

df <- data.frame(mu2=c(mu2.IS, mu2.MC), method=rep(c("IS", "without IS"), each=32))
Iteration <- rep(0:31, 2)

ggplot(df) + geom_hline(yintercept=14.61) +
  geom_point(aes(Iteration, mu2, shape=method, color=method), size=2) +
  ylab(expression(hat(mu)[2])) + ylim(10, 18)

df2 <- data.frame(s12=c(s12.IS, s12.MC), s22=c(s22.IS, s22.MC),
                  method=rep(c("IS", "without IS"), each=22))
Iteration <- rep(0:21, 2)

ggplot(df2) + geom_hline(yintercept=20.885) +
  geom_point(aes(Iteration, s12, shape=method, color=method), size=2) +
  ylab(expression(hat(sigma)[12]))

ggplot(df2) + geom_hline(yintercept=23.573) +
  geom_point(aes(Iteration, s22, shape=method, color=method), size=2) +
  ylab(expression(hat(sigma)[22]))





##### 3-2. Probit-Normal Linear Mixed Model with Salamander Data

library(glmm)
library(truncnorm)
library(mvtnorm)

data(salamander)
n <- 120; data <- salamander[1:n,]
attach(data)

W <- Mate
X <- model.matrix(~Cross-1)
Z.f <- model.matrix(~factor(Female)-1)
Z.m <- model.matrix(~factor(Male)-1)

gibbs.sampling <- function(m){
  
  y <- matrix(nrow=n, ncol=m+1); y[,1] <- W
  
  for(j in 2:(m+1)){
    
    y[,j] <- y[,j-1]
    
    for(i in 1:n){
      mu.i <- X[i,]%*%beta+t(V[-i,i])%*%solve(V[-i,-i])%*%(y[-i,j]-X[-i,]%*%beta)
      sigma.i <- sqrt(V[i,i]-t(V[-i,i])%*%solve(V[-i,-i])%*%V[-i,i])
      y[i,j] <- ifelse(W[i]==1, rtruncnorm(1, 0, Inf, mu.i, sigma.i),
                       rtruncnorm(1, -Inf, 0, mu.i, sigma.i))
    }
  }
  
  invisible(y[,-1])
}


### Initialization

r <- 0; m <- 30; S <- 16
beta <- rep(0, 4); sigma.f <- sqrt(0.2); sigma.m <- sqrt(0.2)

for(s in 1:S){
  
  V <- diag(n)+sigma.f^2*Z.f%*%t(Z.f)+sigma.m^2*Z.m%*%t(Z.m)
  y <- gibbs.sampling(m)
  y.bar <- rowMeans(y); mu <- X%*%beta; q <- ncol(Z.f)
  beta <- rowMeans(solve(t(X)%*%solve(V)%*%X)%*%t(X)%*%solve(V)%*%y)
  cov <- diag(0, n); for(j in 1:m){ cov <- cov +(y[,j]-y.bar)%*%t(y[,j]-y.bar)/m }
  
  sigma.f <- sqrt(mean(diag(sigma.f^2*t(Z.f)%*%solve(V)%*%(y.bar%*%t(y.bar)+covy.bar%*%t(mu)-
                       mu%*%t(y.bar)+mu%*%t(mu))%*%(sigma.f^2*solve(V)%*%Z.f)+sigma.f^2*diag(q)-
                       sigma.f^4*t(Z.f)%*%solve(V)%*%Z.f)))
  sigma.m <- sqrt(mean(diag(sigma.m^2*t(Z.m)%*%solve(V)%*%(y.bar%*%t(y.bar)+covy.mu%*%t(y.bar)+
                       mu%*%t(mu))%*%(sigma.m^2*solve(V)%*%Z.m)+sigma.m^2*diag(q)-
                       sigma.m^4*t(Z.m)%*%solve(V)%*%Z.m)))
  
}

beta0 <- beta; sigma.f0 <- sigma.f; sigma.m0 <- sigma.m
c(beta0, sigma.f0, sigma.m0)


### Iterations

y <- gibbs.sampling(m=30)

L <- function(beta, sigma.f, sigma.m){
  V <- diag(n)+sigma.f^2*Z.f%*%t(Z.f)+sigma.m^2*Z.m%*%t(Z.m)
  apply(y, 2, function(x){ prod(dmvnorm(x, X%*%beta, V)) })
}

beta <- matrix(nrow=65+1, ncol=4); beta[1,] <- beta0
sigma.f <- sigma.f0; sigma.m <- sigma.m0


gibbs.sampling <- function(m){
  
  y <- matrix(nrow=n, ncol=m+1); y[,1] <- W
  
  for(j in 2:(m+1)){
    
    y[,j] <- y[,j-1]
    
    for(i in 1:n){
      mu.i <- X[i,]%*%beta[r,]+t(V[-i,i])%*%solve(V[-i,-i])%*%(y[-i,j]-X[-i,]%*%beta[r,])
      sigma.i <- sqrt(V[i,i]-t(V[-i,i])%*%solve(V[-i,-i])%*%V[-i,i])
      y[i,j] <- ifelse(W[i]==1, rtruncnorm(1, 0, Inf, mu.i, sigma.i),
                       rtruncnorm(1, -Inf, 0, mu.i, sigma.i))
    }
  }
  
  invisible(y[,-1])
}


for(r in 1:65){
  
  V <- diag(n)+sigma.f[r]^2*Z.f%*%t(Z.f)+sigma.m[r]^2*Z.m%*%t(Z.m)
  w <- L(beta[r,], sigma.f[r], sigma.m[r])/L(beta0, sigma.f0, sigma.m0)
  y.bar <- rowMeans(y); y.hat <- y%*%w/sum(w); mu <- X%*%beta[r,]
  beta[r+1,] <- solve(t(X)%*%solve(V)%*%X)%*%t(X)%*%solve(V)%*%y.hat
  cov <- diag(0, n); for(j in 1:m){ cov <- cov + (y[,j]-y.bar)%*%t(y[,j]-y.bar)*w[j]/sum(w) }
  
  sigma.f[r+1] <- sqrt(mean(diag(sigma.f[r]^2*t(Z.f)%*%solve(V)%*%(y.hat%*%t(y.hat)+covy.hat%*%t(mu)-
                            mu%*%t(y.hat)+mu%*%t(mu))%*%(sigma.f[r]^2*solve(V)%*%Z.f)+sigma.f[r]^2*diag(q)-
                            sigma.f[r]^4*t(Z.f)%*%solve(V)%*%Z.f)))
  sigma.m[r+1] <- sqrt(mean(diag(sigma.m[r]^2*t(Z.m)%*%solve(V)%*%(y.hat%*%t(y.hat)+covy.hat%*%t(mu)-
                            mu%*%t(y.hat)+mu%*%t(mu))%*%(sigma.m[r]^2*solve(V)%*%Z.m)+sigma.m[r]^2*diag(q)-
                            sigma.m[r]^4*t(Z.m)%*%solve(V)%*%Z.m)))
  
}