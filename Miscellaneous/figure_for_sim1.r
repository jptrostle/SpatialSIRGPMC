# The figure in simulation study #1

rm(list = ls())
library(rTensor)
library(FNN)
library(geoR)
library(truncnorm)
library(stats)
library(Matrix)
library(sp)
library(mvtnorm)
library(tmvtnorm)
set.seed(560788664)
setwd("..")

this_run_code <- "sim1"
this_em_code <- "sim1"

ur_p <- 1
ns <- 25
start_day <- 61
nt <- 80
K <- 10000
lnt <- nt + 1

# SVD kriging nearest-neighbor count
muxnnc <- 10
Sxxnnc <- 10

# SVD kriging for means: how many basis functions to use?
muxsbfc <- 20
muxtbfc <- 10

# SVD kriging for covariances: how many basis functions to use?
Sxxsbfc <- 10
Sxxtbfc <- 10

alphabf <- 5

run_mux_em <- function(theta0, Theta, nnc, nsbf, ntbf, M, U1, U2, tm, ts) {
    
    theta0[1:2] <- (theta0[1:2] - tm) / ts
    theta0 <- matrix(theta0, nrow = 1)
    nn <- c(get.knnx(data = Theta, query = theta0, k = nnc)$nn.index)
    ThetaR <- rbind(Theta[nn, ], theta0)
    distances <- spDists(ThetaR)
    Corr <- matern(distances, 0.002, 2.5)
    
    # setting up kriging details
    Sigma11Inv <- solve(Corr[1:nnc, 1:nnc])
    Sigma12 <- Corr[1:nnc, (nnc + 1)]
    j <- rep(1, nnc)
    tjS11Inv <- t(j) %*% Sigma11Inv
    S12TS11Inv <- t(Sigma12) %*% Sigma11Inv
    this_denom <- c(tjS11Inv %*% j)
    weights <- tjS11Inv / this_denom + S12TS11Inv - c(S12TS11Inv %*% j) * tjS11Inv / this_denom
    
    Mpred <- matrix(0, nrow = nsbf, ncol = ntbf)
    for (si in 1:nsbf) {
        for (ti in 1:ntbf) {
            Mpred[si, ti] <- weights %*% M[si, ti, nn]
        }
    }
    
    mux0 <- U1 %*% Mpred %*% t(U2)
    mux0[mux0 < 0.00001] <- 0.00001
    return(list(mux0 = mux0, nn = nn))
    
}

em_data <- paste("Emulator data/", this_em_code, sep = "")
Theta <- readRDS(paste(em_data, "/ClTheta.rds", sep = ""))
muxM <- readRDS(paste(em_data, "/muxM.rds", sep = ""))
muxSBF <- readRDS(paste(em_data, "/muxSBF.rds", sep = ""))
muxTBF <- readRDS(paste(em_data, "/muxTBF.rds", sep = ""))
SxxM <- readRDS(paste(em_data, "/SxxM.rds", sep = ""))
SxxSBF <- readRDS(paste(em_data, "/SxxSBF.rds", sep = ""))
SxxTBF <- readRDS(paste(em_data, "/SxxTBF.rds", sep = ""))

########################################################

source("Functions/full forward equations v3.r")

Adj <- matrix(0, nrow = ns, ncol = ns)
nsd <- sqrt(ns)
for (i in 1:ns) {
    left <- i - 1
    right <- i + 1
    up <- i - nsd
    down <- i + nsd
    if (i %% nsd != 1) Adj[i, left] <- 1
    if (i %% nsd != 0) Adj[i, right] <- 1
    if (i > nsd) Adj[i, up] <- 1
    if (i <= (ns - nsd)) Adj[i, down] <- 1
}
rm(nsd, i, left, right, up, down)

coords <- make_forward_calc_sequences(ns, Adj)
matrix_of_zeroes <- matrix(0, nrow = ns, ncol = ns)

grid_pop <- rep(100000, ns)
muy0 <- rep(0, ns)
muy0[13] <- 100
mux0 <- grid_pop - muy0

Xparam <- update_XYparam_VECT_full(coords, ns, start_day + nt - 1, mux0, muy0, matrix_of_zeroes, matrix_of_zeroes, matrix_of_zeroes, 0.043, 0.025, 0.019, Adj, grid_pop, num_breaks = 10)
dmux_fe <- Xparam$mux[, 1:(start_day + nt - 1)] - Xparam$mux[, 2:(start_day + nt)]

########################################################

Theta_m <- colMeans(Theta[, 1:2])
Theta_s <- apply(Theta[, 1:2], 2, sd)
for (i in 1:2) Theta[, i] <- (Theta[, i] - Theta_m[i]) / Theta_s[i]

res <- run_mux_em(c(0.043, 0.025, 3, 3), Theta, muxnnc, muxsbfc, muxtbfc, muxM, muxSBF, muxTBF, Theta_m, Theta_s)
mux <- res$mux0
dmux <- ur_p * (mux[, 1:nt] - mux[, 2:(nt + 1)])

s <- 1

date_range <- start_day:(start_day + nt - 1)
y <- readRDS(paste("Simulated data/", this_run_code, "/1.rds", sep = ""))
y[is.nan(y)] <- 0
plot(y[s, date_range], type = "l", xlab = "time", ylab = "New infections")
for (i in 2:200) {
    y <- readRDS(paste("Simulated data/", this_run_code, "/", i, ".rds", sep = ""))
    y[is.nan(y)] <- 0
    lines(y[s, date_range])
}

hmm <- res$nn
CSMat <- readRDS(paste(em_data, "/CSMat.rds", sep = ""))
for (x in hmm) {
    
    tmux <- CSMat[, x]
    tmux <- matrix(tmux, nrow = ns, byrow = FALSE)
    tdmux <- tmux[, 1:nt] - tmux[, 2:(nt + 1)]
    lines(tdmux[s, ], col = "orange")
    
}

lines(dmux[s, ], col = "red", lw = 3)
lines(dmux_fe[s, date_range], col = "blue", lw = 3, lty = 2)
legend("topright", c("Simulated data", "Nearest neighbor", "Forward equation mean", "Emulated mean"), lty = c(1, 1, 2, 1), col = c("black", "orange", "blue", "red"), cex = 0.9)


# eof