
library(data.table)
library(ggplot2)
library(ggpubr)
library(BRISC)
set.seed(1)
idx2D <- sample(1:603, 500, replace = F)
drip <- fread("BerberisMeshoDataBt1.txt")
library(scales)
coord2 <- rescale(drip[[1]])
coord1 <- rescale(drip[[2]])
yTtl <- drip[[3]]
geom <- cbind(coord1, coord2)
save(yTtl, geom, idx2D, file = "Real_data.RData")

###fit truncated normal####
load("Real_data.RData")


mle_func_TN <- function(alpha, geom, y) {
  print(alpha)
  n <- nrow(geom)
  xi <- rep(0, n)
  D <- as.matrix(dist(geom))
  covM <- alpha[1]* exp(-alpha[2]*D)
  xi <- (2 * y - 1) * xi
  covM <- outer(2 * y - 1, 2 * y - 1) * covM
  diag(covM) <- diag(covM) + 1
  ret <- TruncatedNormal::pmvnorm(
    mu = rep(0, length(xi)), sigma = covM,
    ub = xi
  )[[1]]
  return(log(ret))
}
alphaVec1 <- sqrt(seq(0, 200, length.out = 25))
alphaVec2 <- sqrt(seq(0.0, 100, length.out = 10))
alphaPool <- expand.grid(alphaVec1, alphaVec2)

t1 <- proc.time()
set.seed(123)
lkVecTLR_TN <- apply(alphaPool, 1, mle_func_TN, geom = geom[idx2D, ], y = yTtl[idx2D])
alphaTN <- alphaPool[which.max(lkVecTLR_TN), ]
t2 <- proc.time()

set.seed(123)
nSub <- length(idx2D)
ySub <- yTtl[idx2D]
geomSub <- geom[idx2D, , drop = F]
covM <- matrix(0, nSub + 1, nSub + 1)
covM[1:nSub, 1:nSub] <- alphaTN[[1]]*exp(-alphaTN[[2]]*(as.matrix(dist(geomSub))))
covM[1:nSub, 1:nSub] <- outer(2 * ySub - 1, 2 * ySub - 1) * covM[1:nSub, 1:nSub]
diag(covM[1:nSub, 1:nSub]) <- diag(covM[1:nSub, 1:nSub]) + 1
covM[nSub + 1, nSub + 1] <- 1 + alphaTN[[1]]
t1 <- proc.time()
startTime <- Sys.time()
denormTN <- TruncatedNormal::pmvnorm(
  mu = rep(0, nSub),
  sigma = covM[1:nSub, 1:nSub],
  lb = rep(-Inf, nSub),
  ub = rep(0, nSub)
)[[1]]
endTime <- Sys.time()
t2 <- proc.time()
t1 <- proc.time()
geomUnknownRnd <- geom[-idx2D, ]
predRnd <- rep(NA, nrow(geomUnknownRnd))
for (i in 1:nrow(geomUnknownRnd))
{
  geomTmp <- geomSub
  geomTmp[, 1] <- geomTmp[, 1] - geomUnknownRnd[i, 1] 
  geomTmp[, 2] <- geomTmp[, 2] - geomUnknownRnd[i, 2] 
  covM[nSub + 1, 1:nSub] <- (2 * ySub - 1) * alphaTN[[1]]*exp(-alphaTN[[2]]*sqrt(rowSums(geomTmp^2)))
  covM[1:nSub, nSub + 1] <- covM[nSub + 1, 1:nSub]
  predRnd[i] <- TruncatedNormal::pmvnorm(
    mu = rep(0, nSub + 1),
    sigma = covM,
    lb = rep(-Inf, nSub + 1),
    ub = rep(0, nSub + 1)
  )[[1]] / denormTN
}
t2 <- proc.time()

save(alphaTN, predRnd, file = "Real_TN.RData")



#######probit_nngp#####
load(file = "Real_data.RData")
load(file = "Real_TN.RData")

alphaVec1 <- sqrt(seq(0, 200, length.out = 25))
alphaVec2 <- sqrt(seq(0.0, 100, length.out = 10))
alphaPool <- expand.grid(alphaVec1, alphaVec2)



mle_BRISC_func_PN <- function(alpha, geom, y, neigh) {
  n <- nrow(geom)
  xi <- rep(0, n)
  ret_total <- Binary_estimation(geom, y, beta = 0, n.neighbors = n.neighbors, sigma.sq = alpha[1], phi = alpha[2], tau.sq = 0, 
                                 neighbor = neigh, verbose = FALSE, mc_iter = 100)
  return(ret_total)
}

t1 <- proc.time()
n.neighbors <- 15
neigh <- BRISC_neighbor(geom[idx2D, ], n.neighbors = n.neighbors)


set.seed(2)
lkVecPN_list <- apply(alphaPool, 1, mle_BRISC_func_PN, geom = geom[idx2D, ], y = yTtl[idx2D], neigh = neigh)
lkVecPN <- sapply(1:length(lkVecPN_list), function(i) mean((colSums(log(lkVecPN_list[[i]]$result$e)))))
alphaPN <- alphaPool[which.max(lkVecPN), ]
print(alphaPN)



spRnd_Grid <- Binary_prediction(lkVecPN_list[[which.max(lkVecPN)]], geom[-idx2D, ], beta = 0, n.neighbors = n.neighbors, sigma.sq = alphaPN[1], phi = alphaPN[2], tau.sq = 0,
                                verbose = FALSE, mc_iter = 100)

predicted_out <- exp(colMeans(log(spRnd_Grid$result$predicted_e)))


t2 <- proc.time()


#####
mean((predicted_out - predRnd)^2)
library(pROC)
auc(yTtl[-idx2D], as.numeric(predicted_out>0.5))
auc(yTtl[-idx2D], as.numeric(predRnd>0.5))
pred_PN <- yTtl
pred_TN <- yTtl
pred_PN[-idx2D] <- predicted_out
pred_TN[-idx2D] <- predRnd
tot_data <- as.data.frame(cbind(rescale(drip[[2]]), rescale(drip[[1]]), pred_PN, pred_TN))
tot_data$Truth <- as.character(yTtl)
tot_data$Probability <- tot_data$pred_PN
col.br <- colorRampPalette(c("midnightblue", "cyan", "yellow", "red"))
col.pal <- col.br(5)

data1 <- tot_data[idx2D,]
p1 <- ggplot(data1, aes(x=V1, y=V2, color = Probability, shape = Truth)) + geom_point() +    scale_colour_gradientn(colors = col.pal) +
  xlab("X-coord") + ylab("Y-coord") + ggtitle(paste0("(a) Training data")) + theme(plot.title = element_text(hjust = 0.5))+ 
  theme(text = element_text(size = 15))
data1 <- tot_data[-idx2D,]
#data1$Truth <- as.character(yTtl[-idx2D])
p2 <- ggplot(data1, aes(x=V1, y=V2, color = pred_PN, shape = Truth)) + geom_point() +    scale_colour_gradientn(colors = col.pal) +
  xlab("X-coord") + ylab("Y-coord") + ggtitle(paste0("(b) Probit-NNGP prediction")) + theme(plot.title = element_text(hjust = 0.5))+ 
  theme(text = element_text(size = 15))
p3 <- ggplot(data1, aes(x=pred_TN, y=pred_PN)) + geom_point() +
  xlab("TN prediction") + ylab("Probit-NNGP prediction") + ggtitle(paste0("(c) TN vs. probit-NNGP prediction")) + theme(legend.title = element_blank()) + theme(plot.title = element_text(hjust = 0.5))+ 
  theme(text = element_text(size = 15)) + geom_abline(slope = 1, intercept = 0, color = "magenta") + theme(text = element_text(size = 15))



ggsave(
  "complete_data_together.png",
  ggarrange(p1, p2, p3, nrow = 3, common.legend = TRUE, legend="right"),
  width = 6,
  height = 12,
  dpi = 300
)
