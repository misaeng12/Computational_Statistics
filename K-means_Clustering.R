MyKmeans <- function(X, K, iter.max=10, nstart) {
  
  WSS <- matrix(nrow=nstart, ncol=K)
  mean.list <- list(); membership.list <- list(); iter <- c()
  
  for(i in 1:nstart){
    
    n <- nrow(X)
    mean <- X[sample(n, K),]
    iter[i] <- 0; error <- 1
    
    while(iter[i] < iter.max & sum(error) > 0){
      which.closest <- function(x){ which.min(dist(rbind(x, mean))[1:K]) }
      membership <- apply(X, 1, which.closest)
      old.mean <- mean
      mean <- aggregate(X, list(membership), mean)[,-1]
      error <- abs(mean - old.mean)
      iter[i] <- iter[i] + 1
    }
    
    for(k in 1:K){
      Ck <- matrix(as.matrix(X)[membership==k,], ncol=ncol(X))
      WSS[i, k] <- sum(apply(Ck, 1, function(x){sum((x-mean[k,])^2)}))
    }
    
    mean.list[[i]] <- mean; membership.list[[i]] <- membership
    
  }
  
  WSS.min <- which.min(rowSums(WSS))
  cluster <- membership.list[[WSS.min]]
  centers <- mean.list[[WSS.min]]
  totss <- sum(apply(X, 1, function(x){sum((x-colMeans(X))^2)}))
  withinss <- WSS[WSS.min,]; tot.withinss <- sum(withinss)
  betweenss <- totss - tot.withinss
  size <- as.vector(table(cluster)); iter <- iter[WSS.min]
  
  print(list(size = size, "Cluster means" = centers,"Clustering vector"=
               cluster, "Within cluster sum of squares by cluster" =
               withinss, "between_SS / total_SS" = betweenss/totss))
  
  return(invisible(list(cluster=cluster, centers=centers, totss=totss,
                        withinss=withinss, tot.withinss=tot.withinss,
                        betweenss=betweenss, size=size, iter=iter)))
  
}


## Example (baseball data)

data <- read.table("baseball.dat.txt", header=T)

kmeans1 <- kmeans(data[,-1], 5, iter.max=10, nstart=10)
kmeans1$tot.withinss
kmeans1$iter

MyKmeans1 <- MyKmeans(data[,-1], 5, iter.max=10, nstart=10)
MyKmeans1$tot.withinss
MyKmeans1$iter