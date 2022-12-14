library(tlrmvnmvt)
library(tlrmvnratio)

for(iter in 1:100){

m <- 100
n <- m * m
geom <- cbind(
  kronecker(seq(0, 1, length.out = m), rep(1, m)),
  kronecker(rep(1, m), seq(0, 1, length.out = m))
)

## generate grid

set.seed(iter+123)
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


set.seed(iter+123)
alpha1 <- 1
alpha2 <- sqrt(30)
alpha <- c(alpha1, alpha2)
prTtl <- pnorm(grf_gen(rbind(geom, geomUnknownRnd, geomUnknownGrid), alpha))#this takes time
prUnknownRnd <- prTtl[(n + 1):(n + nUnknown)]
prUnknownGrid <- prTtl[(n + nUnknown + 1):length(prTtl)]
#z <- matrix(prTtl, m, m)



set.seed(iter+123)
yTtl <- rbinom(n = n + 200, size = 1, prob = prTtl)
y <- yTtl[1:n]



mle_func_TLR <- function(alpha, geom, y) {
  n <- nrow(geom)
  xi <- rep(0, n)
  D <- as.matrix(dist(geom))
  covM <- alpha[1]* exp(-alpha[2]*D)
  xi <- (2 * y - 1) * xi
  covM <- outer(2 * y - 1, 2 * y - 1) * covM
  diag(covM) <- diag(covM) + 1
  ret <- tlrmvnmvt::pmvn(
    lower = rep(-Inf, n), upper = xi, sigma = covM, uselog2 = T, N = 998,
    algorithm = TLRQMC(m = round(sqrt(nrow(covM))), epsl = 1e-4))
  return(ret)
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
  
  set.seed(iter+123)
  lkVecTLR_TLR <- apply(alphaPool, 1, mle_func_TLR, geom = geom[idx2D, ], y = yTtl[idx2D])
  alphaTLR <- alphaPool[which.max(lkVecTLR_TLR), ]
  t2 <- proc.time()
  set.seed(iter+123)
  ySub <- yTtl[idx2D]
  geomSub <- geom[idx2D, , drop = F]
  geomSub[, 1] <- geomSub[, 1] 
  geomSub[, 2] <- geomSub[, 2]
  covM <- matrix(0, nSub + 1, nSub + 1)
  covM[1:nSub, 1:nSub] <- alphaTLR[[1]]*exp(-alphaTLR[[2]]*(as.matrix(dist(geomSub))))
  covM[1:nSub, 1:nSub] <- outer(2 * ySub - 1, 2 * ySub - 1) * covM[1:nSub, 1:nSub]
  diag(covM[1:nSub, 1:nSub]) <- diag(covM[1:nSub, 1:nSub]) + 1
  covM[nSub + 1, nSub + 1] <- 1 + alphaTLR[[1]]
  
  # create the empty vectors of predictive probabilities at random and grid test locations  
  predRnd <- rep(NA, nrow(geomUnknownRnd))
  predGrid <- rep(NA, nrow(geomUnknownGrid))
  
  
  # compute predictive probabilities at the 100 random locations via TLR and save runtime for first prediction
  startTime <- Sys.time()
  for (i in 1:nrow(geomUnknownRnd))
  {
    geomTmp <- geomSub
    geomTmp[, 1] <- geomTmp[, 1] - geomUnknownRnd[i, 1] 
    geomTmp[, 2] <- geomTmp[, 2] - geomUnknownRnd[i, 2] 
    covM[nSub + 1, 1:nSub] <- (2 * ySub - 1) * alphaTLR[[1]]*exp(-alphaTLR[[2]]*sqrt(rowSums(geomTmp^2)))
    covM[1:nSub, nSub + 1] <- covM[nSub + 1, 1:nSub]
    predRnd[i] <- tlrmvnratio::mvn.ratio.tlr(
      upper = 0, sigma = covM, m = mSub,
      N = 998, epsl = 1e-4
    )[[1]]
    if(i == 1)
      endTime <- Sys.time()
  }
  
  # compute predictive probabilities at the 100 grid locations via TLR
  for (i in 1:nrow(geomUnknownGrid))
  {
    geomTmp <- geomSub
    geomTmp[, 1] <- geomTmp[, 1] - geomUnknownGrid[i, 1] 
    geomTmp[, 2] <- geomTmp[, 2] - geomUnknownGrid[i, 2] 
    covM[nSub + 1, 1:nSub] <- (2 * ySub - 1) * alphaTLR[[1]]*exp(-alphaTLR[[2]]*sqrt(rowSums(geomTmp^2)))
    covM[1:nSub, nSub + 1] <- covM[nSub + 1, 1:nSub]
    predGrid[i] <- tlrmvnratio::mvn.ratio.tlr(
      upper = 0, sigma = covM, m = mSub,
      N = 998, epsl = 1e-4
    )[[1]]
  }
  
  # compute and display the runtime and MSEs shown in Table 1 for TLR
  timeCost <- as.numeric(difftime(endTime, startTime, units = "secs"))
  MSERnd_TLR <- sum((prTtl[(n + 1):(n + 100)] - predRnd)^2) / 100
  MSEGrid_TLR <- sum((prTtl[(n + 101):(n + 200)] - predGrid)^2) / 100
  save(MSERnd_TLR, MSEGrid_TLR, file = paste0("../probit_nngp_TLR/TLR_n_",n,"_m_",mSub,"_sigma.sq_",alpha[1],"_phi_",alpha[2],"_iter_",iter,".RData"))
}
}

#collecting the result. 
zip <- list.files(path = "../probit_nngp_TLR", pattern = "TLR_n_10000_m_15", full.names = T)#change m value as required
tot_mat <- matrix(0, 100, 2)
for(j in 1:length(zip)){
    load(zip[j])
    tot_mat[j,] <- c(MSEGrid_TLR, MSERnd_TLR)
}

colMedians(tot_mat)
