###1. Slash Distribution

dslash <- function(y){
  ifelse( y!=0, (1-exp(-y^2/2))/(y^2*sqrt(2*pi)), 1/(2*sqrt(2*pi)) )
}

W <- function(f, g, x){ (f(x)/g(x)) / sum(f(x)/g(x)) }
m <- 100000; n <- 5000

Y1 <- rnorm(m)/runif(m)
X1 <- sample( Y1, n, replace = T, prob = W(dnorm, dslash, Y1) )
X2 <- rnorm(m)
Y2 <- sample( X2, n, replace = T, prob = W(dslash, dnorm, X2) )


library(ggplot2)

x1 <- seq(-3, 3, length.out=n); x2 <- seq(-8, 8, length.out=n)
y1 <- seq(0, 500, 100); y2 <- seq(40, 640, 200)

ggplot() + geom_histogram(aes(X1), color=1, fill="white") +
  geom_line(aes(x1, dnorm(x1)*1250)) +
  labs(title="Normal Samples (Slash Envelope)", x="x", y="Normal Density") +
  scale_y_continuous(breaks=y1, labels=y1/1250)

ggplot() + geom_histogram(aes(Y2), color=1, fill="white") +
  geom_line(aes(x2, dslash(x2)*3200)) +
  labs(title="Slash Samples (Normal Envelope)", x="y", y="Slash Density") +
  scale_y_continuous(breaks=y2, labels=round(y2/3200, 2))



###2. Bayesian Inference

x <- c(8, 3, 4, 3, 1, 7, 2, 6, 2, 7)
lik <- function(lambda){ prod(dpois(x, lambda)) }
n <- 10000; i <- 0; lambda <- list()

for(m in c(10, 20, 50, 100)*10000){
  init <- rlnorm(m, log(4), 0.5); i <- i + 1
  lambda[[i]] <- sample( init, n, replace = T, prob = sapply(init, lik) )
}

post <- function(lambda){ dlnorm(lambda, log(4), 0.5)*prod(dpois(x, lambda)) }

x1 <- seq(0, 10, length.out=100)

ggplot() + geom_histogram(aes(lambda[[1]]), color=1, fill="white") +
  geom_line(aes(x1, sapply(x1, post)*0.81e14)) +
  labs(title="m = 100,000", x=expression(lambda))

ggplot() + geom_histogram(aes(lambda[[2]]), color=1, fill="white") +
  geom_line(aes(x1, sapply(x1, post)*0.81e14)) +
  labs(title="m = 200,000", x=expression(lambda))

ggplot() + geom_histogram(aes(lambda[[3]]), color=1, fill="white") +
  geom_line(aes(x1, sapply(x1, post)*0.81e14)) +
  labs(title="m = 500,000", x=expression(lambda))

ggplot() + geom_histogram(aes(lambda[[4]]), color=1, fill="white") +
  geom_line(aes(x1, sapply(x1, post)*0.81e14)) +
  labs(title="m = 1,000,000", x=expression(lambda))


# Rejection Sampling과 비교

m <- n/0.29
lambda <- rlnorm(m, log(4), 0.5); U <- runif(m)
X <- lambda[ U < sapply(lambda, function(lambda){ lik(lambda)/lik(4.3) }) ]
length(X)

ggplot() + geom_histogram(aes(X), color=1, fill="white") +
  geom_line(aes(x1, sapply(x1, post)*0.8e14)) +
  labs(title="Rejection Sampling", x=expression(lambda))