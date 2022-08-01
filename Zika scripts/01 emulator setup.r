# Generate emulator bf and weights for Zika analysis
# For the sensitivites in the web appendix, change muxsbfc and/or Sxxsbfc as necessary

rm(list = ls())
library(Matrix)
library(parallel)
library(rTensor)
options(warn = 2)
set.seed(34256679)
setwd("..")

source("Functions/full forward equations v3.r")

#############
# (1) Setup #
#############

load("Data/zika Brazil.rdata")
rm(zika_data, states_names, area_km)

sim_em_data_path <- "Emulator files"

true_gamma <- 1.2
ns <- 27
start_day <- start_t <- 41
nt <- 40
start_inf <- 100
S0 <- 10
K <- 100000

muxsbfc <- 20
muxtbfc <- 10
Sxxsbfc <- 10
Sxxtbfc <- 10


# ************************************************ #

#######################
# (2) Run experiments #
#######################

# setup for my forward-equations code
coords <- make_forward_calc_sequences(ns, Adj)
matrix_of_zeroes <- matrix(0, nrow = ns, ncol = ns)

# make Theta
beta1_vector <- seq(-0.15, 0.00, length.out = K)
beta2_vector <- sample(seq(0.14, 0.24, length.out = K), K, replace = FALSE)
phi_vector <- sample(seq(0.06, 0.125, length.out = K), K, replace = FALSE)
Theta <- cbind(beta1_vector, beta2_vector, phi_vector, S0)

X <- cbind(1, log_density)

# inline function
gen_fe <- function(i, Theta, ns, grid_pop, matrix_of_zeroes, coords, start_t, nt, gamma, Adj, sim_em_data_path, start_inf, X, log_density) {
    
    this_beta <- exp(c(X %*% c(Theta[i, 1:2])))
    this_phi <- Theta[i, 3]
    muy0 <- rep(0, ns)
    muy0[Theta[i, 4]] <- start_inf
    mux0 <- grid_pop - muy0
    
    this_sim <- update_XYparam_VECT_full(coords, ns, start_t + nt - 1, mux0, muy0, matrix_of_zeroes, matrix_of_zeroes, matrix_of_zeroes, this_beta, this_phi, gamma, Adj, grid_pop, num_breaks = 10)
    
    saveRDS(as.vector(this_sim$mux[, start_t:(start_t + nt)]), paste(sim_em_data_path, "/mux", i, ".rds", sep = ""))
    saveRDS(this_sim$sigmaxx[, start_t:(start_t + nt)], paste(sim_em_data_path, "/Sxx", i, ".rds", sep = ""))
    
}

emex.lines <- mclapply(1:K, function(x) {gen_fe(x, Theta, ns, gridpop, matrix_of_zeroes, coords, start_day, nt, true_gamma, Adj, sim_em_data_path, start_inf, X, log_density)}, mc.cores = 10)
saveRDS(Theta, paste(sim_em_data_path, "/Theta_setup.rds", sep = ""))

# numeric checks
reject_list <- c()
CSMatrix_mux <- matrix(0, nrow = ns * (nt + 1), ncol = K)
CSArray_sxx <- array(0, dim = c(ns * (ns + 1) / 2, nt + 1, K))
for (k in 1:K) {
    CSMatrix_mux[, k] <- readRDS(paste(sim_em_data_path, "/mux", k, ".rds", sep = ""))

    if (any(is.na(CSMatrix_mux[, k])) || any(CSMatrix_mux[, k] > max(gridpop))) {
        reject_list <- c(reject_list, k)
    } else {
        this_mux <- matrix(CSMatrix_mux[, k], nrow = ns, byrow = FALSE)
        this_dmux <- this_mux[, 1:40] - this_mux[, 2:41]
        if (any(this_dmux < -1)) reject_list <- c(reject_list, k)
    }

    CSArray_sxx[,, k] <- readRDS(paste(sim_em_data_path, "/Sxx", k, ".rds", sep = ""))
}
rm(this_mux, this_dmux)

# change here to throw out invalid results
keep_list <- setdiff(1:K, reject_list)
orig_K <- K
K <- length(keep_list)
CSMatrix_mux <- CSMatrix_mux[, keep_list]
CSArray_sxx <- CSArray_sxx[,, keep_list]
Theta <- Theta[keep_list, ]

# convert CSArray to a tensor right here
SxxTensor <- as.tensor(array(0, dim = c(ns, ns, nt + 1, K)))
s1 <- s2 <- 1
for (i in 1:(ns * (ns + 1) / 2)) {
    
    if (i > 1) {
        s2 <- s2 + 1
        if (s2 >= (ns + 1)) {
            s1 <- s1 + 1
            s2 <- s1
        }
    }
    SxxTensor[s1, s2,, ] <- SxxTensor[s2, s1,, ] <- CSArray_sxx[i,, ]
}

# save everything out
saveRDS(CSMatrix_mux, paste(sim_em_data_path, "/CSMat.rds", sep = ""))
saveRDS(Theta, paste(sim_em_data_path, "/Theta.rds", sep = ""))
saveRDS(SxxTensor, paste(sim_em_data_path, "/SxxTensor.rds", sep = ""))

for (k in 1:orig_K) {
    file.remove(paste(sim_em_data_path, "/mux", k, ".rds", sep = ""))
    file.remove(paste(sim_em_data_path, "/Sxx", k, ".rds", sep = ""))
}


#############################################
# (3) Calculate weights and basis functions #
#############################################
    
lnt <- nt + 1
Theta <- readRDS(paste(sim_em_data_path, "/Theta.rds", sep = ""))
Orig_Theta <- Theta
np <- dim(Orig_Theta)[2]
Theta <- cbind(Orig_Theta[, 1:(np - 1)], Orig_Theta[, np] %% sqrt(ns), ceiling(Orig_Theta[, np] / sqrt(ns)))
Theta[Theta[, np] == 0, np] <- sqrt(ns)
rm(Orig_Theta)

K <- dim(Theta)[1]

# ----------- #
# Means setup #
# ----------- #

CSMat <- readRDS(paste(sim_em_data_path, "/CSMat.rds", sep = ""))
MuxTensor <- array(0, dim = c(ns, lnt, K))
for (k in 1:K) MuxTensor[,, k] <- matrix(CSMat[, k], nrow = ns, byrow = FALSE)
MuxTensor <- as.tensor(MuxTensor)
rm(CSMat)

Spatial <- matrix(0, nrow = ns, ncol = ns)
Temporal <- matrix(0, nrow = lnt, ncol = lnt)
uhohlist <- c()
for (k in 1:K) {
    thisslice <- MuxTensor[,, k]@data
    Spatial <- Spatial + thisslice %*% t(thisslice)
    Temporal <- Temporal + t(thisslice) %*% thisslice
}		
muxSBF <- eigen(Spatial)$vectors[, 1:muxsbfc]
muxTBF <- eigen(Temporal)$vectors[, 1:muxtbfc]

# calculate mean M matrix
muxM <- array(0, dim = c(muxsbfc, muxtbfc, K))
for (k in 1:K) muxM[,, k] <- t(muxSBF) %*% MuxTensor[,, k]@data %*% muxTBF
rm(k, Spatial, Temporal, MuxTensor)


# ----------------- #
# Covariances setup #
# ----------------- #

    SxxTensor <- readRDS(paste(sim_em_data_path, "/SxxTensor.rds", sep = ""))
    Spatial <- matrix(0, nrow = ns, ncol = ns)
    for (k in 1:K) {
        thisslice <- k_unfold(SxxTensor[,,, k], m = 1)@data
        Spatial <- Spatial + thisslice %*% t(thisslice)
    }
    SxxSBF <- eigen(Spatial)$vectors[, 1:Sxxsbfc]
    
    SxxM3M4 <- ttm(ttm(SxxTensor, t(SxxSBF), m = 1), t(SxxSBF), m = 2)
    SxxM3M4 <- SxxM3M4@data
    rm(SxxTensor)
    gc()
    
    Ch <- array(0, dim = c(Sxxsbfc, Sxxsbfc, lnt, K))
    for (t in 1:lnt) {
        for (k in 1:K) {
            if (min(eigen(SxxM3M4[,, t, k])$values) < 0 ) {
                thiseig <- eigen(SxxM3M4[,, t, k])
                thiseig$values[thiseig$values < 0] <- min(thiseig$values[thiseig$values > 0]) / 10000
                recons_matrix <- thiseig$vectors %*% diag(thiseig$values) %*% t(thiseig$vectors)
                Ch[,, t, k] <- t(chol(recons_matrix))
            } else {
                Ch[,, t, k] <- t(chol(SxxM3M4[,, t, k]))    
            }
        }
    }
    rm(SxxM3M4)
    gc()
    
    Ch <- as.tensor(Ch)
    ChU3 <- svd(k_unfold(Ch, 3)@data)$u[, 1:Sxxtbfc]
    ChM <- ttm(Ch, t(ChU3), m = 3)
    ChM <- ChM@data
    
    #rm(t, k, Spatial, thisslice, Ch, thiseig, recons_matrix, Ch_decomp)
    rm(t, k, Spatial, thisslice, Ch, thiseig, recons_matrix)

      
    SxxM <- ChM
    SxxTBF <- ChU3
    rm(ChM, ChU3)


# -------- #
# Save out #
# -------- #

# Cleaned Theta, muxM, muxSBF, muxTBF, SxxM, SxxSBF, SxxTBF
saveRDS(Theta, paste(sim_em_data_path, "/ClTheta.rds", sep = ""))
saveRDS(muxM, paste(sim_em_data_path, "/muxM.rds", sep = ""))
saveRDS(muxSBF, paste(sim_em_data_path, "/muxSBF.rds", sep = ""))
saveRDS(muxTBF, paste(sim_em_data_path, "/muxTBF.rds", sep = ""))
saveRDS(SxxM, paste(sim_em_data_path, "/SxxM.rds", sep = ""))
saveRDS(SxxSBF, paste(sim_em_data_path, "/SxxSBF.rds", sep = ""))
saveRDS(SxxTBF, paste(sim_em_data_path, "/SxxTBF.rds", sep = ""))



# eof