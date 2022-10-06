library(BRISC)
library(matrixStats)

m <- 100
n <- m * m
geom <- cbind(
  kronecker(seq(0, 1, length.out = m), rep(1, m)),
  kronecker(rep(1, m), seq(0, 1, length.out = m))
)

## generate grid

set.seed(123)
nUnknown <- 100
geomTmp <- geom[(geom[, 1] < 1 - 0.9 / (m - 1)) &
                  (geom[, 2] < 1 - 0.9 / (m - 1)), ]

# random locations
geomUnknownRnd <- geomTmp[sample(1:(n - 2 * m + 1), nUnknown, F), ] +
  matrix(runif(nUnknown, 0.2 / m, 0.8 / m), nUnknown, 2)

# grid locations
geomUnknownGrid <- cbind(
  kronecker(seq(0.21, 0.39, length.out = 10), rep(1, 10)),
  kronecker(rep(1, 10), seq(0.41, 0.59, length.out = 10))
)

##

rmvn <- function(n, mu = 0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension not right!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}



grf_gen <- function(geom, alpha) {
  if (ncol(geom) != 2) {
    stop("grf_gen: only works under 2D geometry")
  }
  if (min(geom) < 0 || max(geom) > 1) {
    warning("grf_gen: geom is not in the unit square")
  }
  n <- nrow(geom)
  D <- as.matrix(dist(geom))
  R <- exp(-alpha[2]*D)
  return(rmvn(1, rep(0,n), alpha[1]*R))
}


set.seed(123)
alpha1 <- 1
alpha2 <- sqrt(30)
alpha <- c(alpha1, alpha2)
prTtl <- pnorm(grf_gen(rbind(geom, geomUnknownRnd, geomUnknownGrid), alpha))#this takes time
prUnknownRnd <- prTtl[(n + 1):(n + nUnknown)]
prUnknownGrid <- prTtl[(n + nUnknown + 1):length(prTtl)]
#z <- matrix(prTtl, m, m)



set.seed(123)
yTtl <- rbinom(n = n + 200, size = 1, prob = prTtl)
y <- yTtl[1:n]



mle_BRISC_func_PN <- function(alpha, geom, y, neigh) {
  n <- nrow(geom)
  xi <- rep(0, n)
  ret_total <- Binary_estimation(geom, y, beta = 0, n.neighbors = n.neighbors, sigma.sq = alpha[1], phi = alpha[2], tau.sq = 0, 
                                 neighbor = neigh, verbose = FALSE, mc_iter = 100)
  return(ret_total)
}


# alpha grid
alphaVec1 <- sqrt(seq(15, 45, length.out = 10))/sqrt(30)
alphaVec2 <- sqrt(seq(15, 45, length.out = 10))
alphaPool <- expand.grid(alphaVec1, alphaVec2)


for(mSub in c(15, 25, 50, 100)){
  t1 <- proc.time()
  nSub <- mSub^2
  idx1D <- round(seq(1, m, length.out = mSub))
  idx2D <- c(kronecker(idx1D - 1, rep(m, mSub)) + idx1D)
  
  n.neighbors <- 15
  neigh <- BRISC_neighbor(geom[idx2D, ], n.neighbors = n.neighbors)
  
  set.seed(123)
  lkVecPN_list <- apply(alphaPool, 1, mle_BRISC_func_PN, geom = geom[idx2D, ], y = yTtl[idx2D], neigh = neigh)
  lkVecPN <- sapply(1:length(lkVecPN_list), function(i) mean((colSums(log(lkVecPN_list[[i]]$result$e)))))
  alphaPN <- alphaPool[which.max(lkVecPN), ]
  
  t2 <- proc.time()
  
  t3 <- proc.time()
  spRnd_Grid <- Binary_prediction(lkVecPN_list[[which.max(lkVecPN)]], rbind(geomUnknownRnd, geomUnknownGrid), beta = 0, n.neighbors = n.neighbors, sigma.sq = alphaTN[1], phi = alphaTN[2], tau.sq = 0,
                                       verbose = FALSE, mc_iter = 100)
  
  predicted_out <- exp(colMeans(log(spRnd_Grid$result$predicted_e)))
  t4 <- proc.time()
  
  
  MSERnd <- sum((prTtl[(n + 1):(n + 100)] - predicted_out[1:100])^2) / 100
  MSEGrid <- sum((prTtl[(n + 101):(n + 200)] - predicted_out[101:200])^2) / 100
  
}



