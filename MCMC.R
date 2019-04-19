#####----- Metropolis-Hastings Algorithm -----#####

### GMM (Gausian Mixture Model)


#1. Independence Chains

n <- 100; d <- 0.7; y <- c(rnorm(d*n, 7, 0.5), rnorm((1-d)*n, 10, 0.5))

library(ggplot2)
ggplot() + geom_histogram(aes(y), color="white")

n.iter <- 10000; delta <- matrix(nrow=n.iter, ncol=2); delta1 <- runif(1)
L <- function(delta){ prod((delta*dnorm(y, 7, 0.5)+(1-delta)*dnorm(y, 10, 0.5))) }

for(t in 1:n.iter){
  delta[t,] <- delta1; delta.star <- c(rbeta(1, 1, 1), rbeta(1, 2, 10))
  R <- sapply(delta.star, L)/sapply(delta[t,], L)
  delta1 <- ifelse(R>=1, delta.star, ifelse(rbinom(2, 1, R), delta.star, delta[t,]))
}

ggplot() + geom_line(aes(1:n.iter, delta[,1])) +
  labs(title="Independence Chain with Beta(1,1)", x="t", y=expression(delta^(t)))
ggplot() + geom_line(aes(1:n.iter, delta[,2])) +
  labs(title="Independence Chain with Beta(2,10)", x="t", y=expression(delta^(t)))

delta.new <- delta[-(1:500), ]

rbind(apply(delta.new, 2, summary), apply(delta.new, 2, sd))
ggplot() + geom_histogram(aes(delta.new[,1]), bins=20, color="white") +
  labs(title="Independence Chain with Beta(1,1)", x=expression(delta^(t)))
ggplot() + geom_histogram(aes(delta.new[,2]), bins=20, color="white") +
  labs(title="Independence Chain with Beta(2,10)", x=expression(delta^(t)))



#2. Random Walk Chains

u <- matrix(nrow=n.iter, ncol=2); u1 <- runif(1, -1, 1); b <- c(1, 0.01)

myf <- function(u){ L(exp(u)/(1+exp(u)))*(exp(u)/(1+exp(u))^2) }

for(t in 1:n.iter){
  u[t,] <- u1; u.star <- u[t,] + runif(2, -b, b)
  R <- sapply(u.star, myf)/sapply(u[t,], myf)
  u1 <- ifelse(R>=1, u.star, ifelse(rbinom(2, 1, R), u.star, u[t,]))
}

delta <- exp(u)/(1+exp(u))

ggplot() + geom_line(aes(1:n.iter, delta[,1])) +
  labs(title="Random Walk Chain with b=1", x="t", y=expression(delta^(t)))
ggplot() + geom_line(aes(1:n.iter, delta[,2])) +
  labs(title="Random Walk Chain with b=0.01", x="t", y=expression(delta^(t)))

delta.new <- delta[-(1:500), ]

rbind(apply(delta.new, 2, summary), apply(delta.new, 2, sd))
ggplot() + geom_histogram(aes(delta.new[,1]), bins=20, color="white") +
  labs(title="Histogram of Random Walk Chain with b=1", x=expression(delta^(t)))
ggplot() + geom_histogram(aes(delta.new[,2]), bins=20, color="white") +
  labs(title="Histogram of Random Walk Chain with b=0.01", x=expression(delta^(t)))





#####---------- Gibbs Sampling ----------#####

x.H <- c(2, 4, 6, 9, 9, 9, 13, 14, 18, 23, 31, 32, 33, 34, 43, 10, 14, 14, 16, 17, 18, 18, 19, 20, 20, 21, 21, 23, 24, 29, 29, 30, 30, 31, 31, 31, 33, 35, 37, 40, 41, 42, 42, 44, 46, 48, 49, 51, 53, 54, 54, 55, 56)
x.C <- c(1, 4, 6, 7, 13, 24, 25, 35, 35, 39, 1, 1, 3, 4, 5, 8, 10, 11, 13, 14, 14, 15, 17, 19, 20, 22, 24, 24, 24, 25, 26, 26, 26, 28, 29, 29, 32, 35, 38, 39, 40, 41, 44, 45, 47, 47, 47, 50, 50, 51)
delta.H <- rep(c(1, 0), c(15, 38)); delta.C <- rep(c(1, 0), c(10, 40))

a=3; b=1; c=60; d=120
a=3/2; b=1/2; c=60/2; d=120/2
a=3*2; b=1*2; c=60*2; d=120*2


length(x.H); length(x.C); mean(delta.H); mean(delta.C)
tapply(x.H, delta.H, mean); tapply(x.C, delta.C, mean)
tapply(x.H, delta.H, sd); tapply(x.C, delta.C, sd)

ggplot() + geom_boxplot(aes(factor(delta.H), x.H)) + ylim(0, 60) + labs(x=expression(delta[i]^H), y=expression(x[i]^H))
ggplot() + geom_boxplot(aes(factor(delta.C), x.C)) + ylim(0, 60) + labs(x=expression(delta[i]^C), y=expression(x[i]^C))

theta <- (a-b+1)/c; tau <- ((b+1)/d)/theta
for(t in 2:n.iter){
  theta[t] <- rgamma(1, a+sum(delta.C)+sum(delta.H)+1, c+sum(x.C)+tau[t-1]*(d+sum(x.H)))
  tau[t] <- rgamma(1, b+sum(delta.H)+1, theta[t]*(d+sum(x.H)))
}

ggplot() + geom_line(aes(1:n.iter, theta)) + labs(x="t", y=expression(theta^(t)))
ggplot() + geom_line(aes(1:n.iter, tau)) + labs(x="t", y=expression(tau^(t)))

acf(theta, lag.max=20, main=expression(theta^(t)))
acf(tau, lag.max=20, main=expression(tau^(t)))

theta <- theta[-(1:500)]; tau <- tau[-(1:500)]
mean(theta); mean(tau); mean(tau*theta); sd(theta); sd(tau); sd(tau*theta)
c(quantile(theta, c(0.025, 0.975))); c(quantile(tau, c(0.025, 0.975)))
c(quantile(tau*theta, c(0.025, 0.975)))

prior <- function(x) x^b*gamma(a+1)/(c+d*x)^(a+1)
c <- 1/integrate(prior, 0, 5)$value; x <- seq(0, 5, length.out=100)
#distribution <- factor(rep(c("prior", "posterior"), c(99, 1)))
ggplot() + geom_density(aes(tau), bw="nrd") + geom_line(aes(x, c*prior(x), lty=distribution)) + labs(x=expression(tau))


c(1/mean(tau*theta), 1/mean(theta), mean(tau))
