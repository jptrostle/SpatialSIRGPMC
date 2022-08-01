# Streaming-enabled function for calculating factor matrices and weights for moment-closure emulators

generate_em_bf_weights <- function(sim_em_data_path, Theta, ns, nt, K) {

    lnt <- nt + 1
    Theta <- readRDS(paste(sim_em_data_path, "/Theta.rds", sep = ""))
    Orig_Theta <- Theta
    np <- dim(Orig_Theta)[2]
    Theta <- cbind(Orig_Theta[, 1:(np - 1)], Orig_Theta[, np] %% sqrt(ns), ceiling(Orig_Theta[, np] / sqrt(ns)))
    Theta[Theta[, np] == 0, np] <- sqrt(ns)
    rm(Orig_Theta)
            
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
    Ch_decomp <- hosvd(Ch, ranks = c(Sxxsbfc, Sxxsbfc, Sxxtbfc, K))
    ChU3 <- Ch_decomp$U[[3]]
    ChM <- ttm(Ch, t(ChU3), m = 3)
    ChM <- ChM@data
    
    rm(t, k, Spatial, thisslice, Ch, thiseig, recons_matrix, Ch_decomp)
    
    SxxM <- ChM
    SxxTBF <- ChU3
    rm(ChM, ChU3)
    
    
    # ----------------------- #
    # Now save everything out #
    # ----------------------- #
    
    # Cleaned Theta, muxM, muxSBF, muxTBF, SxxM, SxxSBF, SxxTBF
    saveRDS(Theta, paste(sim_em_data_path, "/ClTheta.rds", sep = ""))
    saveRDS(muxM, paste(sim_em_data_path, "/muxM.rds", sep = ""))
    saveRDS(muxSBF, paste(sim_em_data_path, "/muxSBF.rds", sep = ""))
    saveRDS(muxTBF, paste(sim_em_data_path, "/muxTBF.rds", sep = ""))
    saveRDS(SxxM, paste(sim_em_data_path, "/SxxM.rds", sep = ""))
    saveRDS(SxxSBF, paste(sim_em_data_path, "/SxxSBF.rds", sep = ""))
    saveRDS(SxxTBF, paste(sim_em_data_path, "/SxxTBF.rds", sep = ""))

}

# eof