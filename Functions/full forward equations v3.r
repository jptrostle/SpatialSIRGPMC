# All the code for the moment-closure forward equations
# The first function is to account for memory-efficient storage of the forward equations in the main function
# It's complicated! A simple brute-force series of nested loops will work for the forward equations, but it's slower

make_forward_calc_sequences <- function(ns, Adj) {
    
    library(Matrix)

    # some coordinates to set up for sxx and syy
    coords_same <- c()
    index <- 1
    for (i in 1:ns) {
        for (j in i:ns) {
            coords_same[index] <- paste(i, "-", j, sep = "")
            index <- index + 1
        }
    }
    
    # coordinates to set up for sxy
    coords_xy <- c()
    index <- 1
    for (i in 1:ns) {
        for (j in 1:ns) {
            coords_xy[index] <- paste(i, "-", j, sep = "")
            index <- index + 1
        }
    }
    
    # make hideous indices for storing the off-diagonal elements
    # examples using ns = 4
    # SEQ1 goes c(1, 1, 1, 2, 2, 3)
    # SEQ2 goes c(2, 3, 4, 3, 4, 4)
    # notice the above two line up to go (1,2), (1,3), (1,4), (2,3), and so on. It matches the order of the offdiag sxx storage
    # SEQ3 goes c(5, 9, 13, 10, 14, 15)
    # SEQ4 goes c(2, 3, 4, 7, 8, 12)
    # those are confusing but they are the indices for the sxy interactions that are needed for the offdiag sxx updates
    # for example SEQ3 is for the (j, s) coordinates (2,1), (3,1), and so on
    # then SEQ4 makes more sense because it's (1,2), (1,3), and so on. It's the (s, j) coordinates
    # SEQ5 is for the off-diag coordinates in the whole sxx vector
    SEQ1 <- rep(1:ns, times = seq(ns - 1, 0))
    SEQ2 <- c()
    for (i in 2:ns) SEQ2 <- c(SEQ2, i:ns)
    SEQ3 <- SEQ4 <- c()
    for (i in 1:(ns - 1)) {
        for (j in (i + 1):ns) {
            s4str <- paste(i, "-", j, sep = "")
            s3str <- paste(j, "-", i, sep = "")
            SEQ3 <- c(SEQ3, which(coords_xy == s3str))
            SEQ4 <- c(SEQ4, which(coords_xy == s4str))
        }
    }
    SEQ5 <- c()
    for (i in 1:(ns - 1)) {
        for (j in (i + 1):ns) {
            SEQ5 <- c(SEQ5, which(coords_same == paste(i, "-", j, sep = "")))
        }
    }
    
    # these are the massive adjacency matrices needed for the things like sxy(j, r) where r is adjacent to s (that's AdjMat1)
    # then AdjMat2 is only going to have ns columns because it doesn't need to account for the shifts in sxy
    AdjMat1 <- AdjMat3 <- Matrix(0, nrow = length(SEQ5), ncol = ns^2, sparse = TRUE)
    AdjMat2 <- AdjMat4 <- Matrix(0, nrow = length(SEQ5), ncol = ns, sparse = TRUE)
    for (i in 1:length(SEQ5)) {
        dash_loc <- as.numeric(regexpr("-", coords_same[SEQ5[i]]))
        this_adj_row <- as.numeric(substr(coords_same[SEQ5[i]], 1, dash_loc - 1))
        this_partition_column <- as.numeric(substr(coords_same[SEQ5[i]], dash_loc + 1, nchar(coords_same[SEQ5[i]])))
        column_index <- seq((this_partition_column - 1) * ns + 1, this_partition_column * ns)
        AdjMat1[i, column_index] <- Adj[this_adj_row, ]
        AdjMat2[i, ] <- Adj[this_adj_row, ]
    }
    # yes this could be embedded above, but it's clearer to make the variable names 'right'
    for (i in 1:length(SEQ5)) {
        dash_loc <- as.numeric(regexpr("-", coords_same[SEQ5[i]]))
        this_partition_column <- as.numeric(substr(coords_same[SEQ5[i]], 1, dash_loc - 1))
        this_adj_row <- as.numeric(substr(coords_same[SEQ5[i]], dash_loc + 1, nchar(coords_same[SEQ5[i]])))
        column_index <- seq((this_partition_column - 1) * ns + 1, this_partition_column * ns)
        AdjMat3[i, column_index] <- Adj[this_adj_row, ]
        AdjMat4[i, ] <- Adj[this_adj_row, ]
    }
    rm(i, this_adj_row, this_partition_column)
    
    # AdjMat5 and AdjMat6 are more complicated... they are for (r, j) and (k, s) coordinates but postmultiplied by syy
    AdjMat5 <- AdjMat6 <- Matrix(0, nrow = length(SEQ5), ncol = (ns + 1) * ns / 2, sparse = TRUE)
    for (i in 1:length(SEQ5)) {
        
        this_syy_code <- coords_same[SEQ5[i]]
        dash_loc <- as.numeric(regexpr("-", this_syy_code))
        s_val <- as.numeric(substr(this_syy_code, 1, dash_loc - 1))
        j_val <- as.numeric(substr(this_syy_code, dash_loc + 1, nchar(this_syy_code)))
        
        r_vals <- which(Adj[s_val, ] == 1)
        for (this_r in r_vals) {
            if (this_r <= j_val) {
                this_code <- paste(this_r, "-", j_val, sep = "")
            } else {
                this_code <- paste(j_val, "-", this_r, sep = "")
            }
            
            col_index <- which(coords_same == this_code)
            AdjMat5[i, col_index] <- 1
            
        }
        
        k_vals <- which(Adj[j_val, ] == 1)
        for (this_k in k_vals) {
            if (this_k <= s_val) {
                this_code <- paste(this_k, "-", s_val, sep = "")
            } else {
                this_code <- paste(s_val, "-", this_k, sep = "")
            }
            
            col_index <- which(coords_same == this_code)
            AdjMat6[i, col_index] <- 1
        }
        
    }
    rm(col_index, this_code, this_k, this_r, k_vals, r_vals, s_val, j_val, this_syy_code, dash_loc, i)
    
    # now for more SEQ and AdjMat stuff for the off-diag sxy updates...
    SEQ6 <- SEQ7 <- SEQ8 <- SEQ9 <- c()
    AdjMat7 <- Matrix(0, nrow = (ns - 1) * ns, ncol = ns * (ns + 1) / 2, sparse = TRUE)
    AdjMat8 <- AdjMat10 <- Matrix(0, nrow = (ns - 1) * ns, ncol = ns, sparse = TRUE)
    AdjMat9 <- Matrix(0, nrow = (ns - 1) * ns, ncol = ns^2, sparse = TRUE)
    mat_row_index <- 1
    for (i in 1:(ns^2)) {
        
        this_sxy_code <- coords_xy[i]
        dash_loc <- as.numeric(regexpr("-", this_sxy_code))
        s_val <- as.numeric(substr(this_sxy_code, 1, dash_loc - 1))
        j_val <- as.numeric(substr(this_sxy_code, dash_loc + 1, nchar(this_sxy_code)))
        
        if (s_val == j_val) next
        
        SEQ6 <- c(SEQ6, i)
        SEQ7 <- c(SEQ7, s_val)
        SEQ8 <- c(SEQ8, j_val)
        
        if (s_val < j_val) {
            lookup_code <- paste(s_val, "-", j_val, sep = "")
        } else {
            lookup_code <- paste(j_val, "-", s_val, sep = "")
        }
        
        xx_yy_index <- which(coords_same == lookup_code)
        SEQ9 <- c(SEQ9, xx_yy_index)
        
        # now handle the 'r' neighbors, which interacts with sigmayy
        r_vals <- which(Adj[s_val, ] == 1)
        for (this_r in r_vals) {
            if (this_r <= j_val) {
                syy_code <- paste(this_r, "-", j_val, sep = "")
            } else {
                syy_code <- paste(j_val, "-", this_r, sep = "")
            }
            syy_match <- which(coords_same == syy_code)
            AdjMat7[mat_row_index, syy_match] <- 1
            
            AdjMat8[mat_row_index, this_r] <- 1
        }
        
        # now handle the 'k' neighbors, which confusingly interacts with sigmaxy
        k_vals <- which(Adj[j_val, ] == 1)
        for (this_k in k_vals) {
            
            # somewhat easier here: it goes s-k
            sxy_lookup_code <- paste(s_val, "-", this_k, sep = "")
            sxy_coord_loc <- which(coords_xy == sxy_lookup_code)
            AdjMat9[mat_row_index, sxy_coord_loc] <- 1
            
            AdjMat10[mat_row_index, this_k] <- 1
            
        }
        
        mat_row_index <- mat_row_index + 1
        
    }
    stopifnot(length(SEQ6) == ns * (ns - 1))
    
    # now for stuff for mux / muy / diag covariance updates
    SEQ_DIAG_XY <- c()
    AdjMat_DiagXY <- Matrix(0, nrow = ns, ncol = ns^2, sparse = TRUE)
    for (i in 1:(ns^2)) {
        
        this_sxy_code <- coords_xy[i]
        dash_loc <- as.numeric(regexpr("-", this_sxy_code))
        s_val <- as.numeric(substr(this_sxy_code, 1, dash_loc - 1))
        j_val <- as.numeric(substr(this_sxy_code, dash_loc + 1, nchar(this_sxy_code)))
        
        if (s_val == j_val) SEQ_DIAG_XY <- c(SEQ_DIAG_XY, i)
        
        r_vals <- which(Adj[s_val, ] == 1)
        for (this_r in r_vals) {
            
            sxy_lookup_code <- paste(s_val, "-", this_r, sep = "")
            sxy_coord_loc <- which(coords_xy == sxy_lookup_code)
            AdjMat_DiagXY[s_val, sxy_coord_loc] <- 1
            
        }
    }
    
    AdjMat_DiagSame <- Matrix(0, nrow = ns, ncol = ns * (ns + 1) / 2, sparse = TRUE)
    SEQ_DIAG_XX_XY <- c()
    for (i in 1:(ns * (ns + 1) / 2)) {
        
        this_sxx_code <- coords_same[i]
        dash_loc <- as.numeric(regexpr("-", this_sxx_code))
        s_val <- as.numeric(substr(this_sxx_code, 1, dash_loc - 1))
        j_val <- as.numeric(substr(this_sxx_code, dash_loc + 1, nchar(this_sxx_code)))
        
        if (s_val == j_val) SEQ_DIAG_XX_XY <- c(SEQ_DIAG_XX_XY, i)
        
        r_vals <- which(Adj[s_val, ] == 1)
        for (this_r in r_vals) {
            
            if (s_val <= this_r) {
                sxx_syy_code <- paste(s_val, "-", this_r, sep = "")
            } else {
                sxx_syy_code <- paste(this_r, "-", s_val, sep = "")
            }
            
            sxx_syy_coord_loc <- which(coords_same == sxx_syy_code)
            AdjMat_DiagSame[s_val, sxx_syy_coord_loc] <- 1
        }
        
    }
    
    m_vec <- as.vector(Adj %*% rep(1, ns))
    return(list(SEQ1 = SEQ1, SEQ2 = SEQ2, SEQ3 = SEQ3, SEQ4 = SEQ4, SEQ5 = SEQ5, SEQ6 = SEQ6, SEQ7 = SEQ7, SEQ8 = SEQ8, SEQ9 = SEQ9,
                AdjMat1 = AdjMat1, AdjMat2 = AdjMat2, AdjMat3 = AdjMat3, AdjMat4 = AdjMat4, AdjMat5 = AdjMat5,
                AdjMat6 = AdjMat6, AdjMat7 = AdjMat7, AdjMat8 = AdjMat8, AdjMat9 = AdjMat9, AdjMat10 = AdjMat10,
                AdjMat_DiagXY = AdjMat_DiagXY, AdjMat_DiagSame = AdjMat_DiagSame,
                SEQ_DIAG_XY = SEQ_DIAG_XY, SEQ_DIAG_XX_XY = SEQ_DIAG_XX_XY,
                m_vec = m_vec))
    
    
}




# ************************************************** #
# -------------------------------------------------- #
# ************************************************** #



update_XYparam_VECT_full <- function(coordinates, ns, nt, start_mux, start_muy, start_sigmaxx, start_sigmaxy, start_sigmayy, beta, phi, gamma, adj, grid_pop, num_breaks = 10, return_all = FALSE) {
    
    library(Matrix)

    nt <- nt + 1
    nt_orig <- nt
    nt <- (nt - 1) * num_breaks + 1
    beta <- beta / num_breaks
    phi <- phi / num_breaks
    gamma <- gamma / num_breaks
    
    mux <- muy <- matrix(0, nrow = ns, ncol = nt)
    sigmaxx <- sigmayy <- matrix(0, nrow = (ns * (ns + 1) / 2), ncol = nt)
    sigmaxy <- matrix(0, nrow = ns^2, ncol = nt)
    mux[, 1] <- as.vector(start_mux)
    muy[, 1] <- as.vector(start_muy)
    
    beta <- beta / grid_pop
    phi <- phi / grid_pop
    
    # translate matrix sxx0, sxy0, and syy0 to the vectorized versions
    same_index <- 1
    xy_index <- 1
    for (i in 1:ns) {
        for (j in 1:ns) {
            if (i <= j) {
                sigmaxx[same_index, 1] <- start_sigmaxx[i, j]
                sigmayy[same_index, 1] <- start_sigmayy[i, j]
                same_index <- same_index + 1
            }
            
            sigmaxy[xy_index, 1] <- start_sigmaxy[i, j]
            xy_index <- xy_index + 1
            
        }
    }
    rm(i, j, same_index, xy_index)
    
    # unpack the premade adjacency matrices and sequences
    SEQ1 <- coordinates$SEQ1
    SEQ2 <- coordinates$SEQ2
    SEQ3 <- coordinates$SEQ3
    SEQ4 <- coordinates$SEQ4
    SEQ5 <- coordinates$SEQ5
    SEQ6 <- coordinates$SEQ6
    SEQ7 <- coordinates$SEQ7
    SEQ8 <- coordinates$SEQ8
    SEQ9 <- coordinates$SEQ9
    AdjMat1 <- coordinates$AdjMat1
    AdjMat2 <- coordinates$AdjMat2
    AdjMat3 <- coordinates$AdjMat3
    AdjMat4 <- coordinates$AdjMat4
    AdjMat5 <- coordinates$AdjMat5
    AdjMat6 <- coordinates$AdjMat6
    AdjMat7 <- coordinates$AdjMat7
    AdjMat8 <- coordinates$AdjMat8
    AdjMat9 <- coordinates$AdjMat9
    AdjMat10 <- coordinates$AdjMat10
    SEQ_DIAG_XY <- coordinates$SEQ_DIAG_XY
    SEQ_DIAG_XX_XY <- coordinates$SEQ_DIAG_XX_XY
    AdjMat_DiagXY <- coordinates$AdjMat_DiagXY
    AdjMat_DiagSame <- coordinates$AdjMat_DiagSame
    
    betaSEQ1 <- beta[SEQ1]
    betaSEQ2 <- beta[SEQ2]
    phiSEQ1 <- phi[SEQ1]
    phiSEQ2 <- phi[SEQ2]
    betaSEQ7 <- beta[SEQ7]
    phiSEQ7 <- phi[SEQ7]
    betaSEQ8 <- beta[SEQ8]
    phiSEQ8 <- phi[SEQ8]
    
    for (t in 2:nt) {
        
        prev <- t - 1
        
        # precalculations
        muxprev <- mux[, prev]
        muyprev <- muy[, prev]
        sigmayyprev <- sigmayy[, prev]
        sigmaxyprev <- sigmaxy[, prev]
        adj_mux <- as.vector(adj %*% muxprev)
        adj_muy <- as.vector(adj %*% muyprev)
        muxtimesmuy <- mux[, prev] * muyprev
        diag_sxy <- sigmaxy[SEQ_DIAG_XY, prev]
        muxtimesadjmuy <- muxprev * adj_muy
        AMDXY <- as.vector(AdjMat_DiagXY %*% sigmaxyprev)
        diag_sxx <- sigmaxx[SEQ_DIAG_XX_XY, prev]
        muytimesdiagsxx <- muyprev * diag_sxx
        muxtimesdiagsxy <- muxprev * diag_sxy
        adjmuytimesdiagsxx <- adj_muy * diag_sxx
        diag_syy <- sigmayy[SEQ_DIAG_XX_XY, prev]
        muxAMDSsyy <- muxprev * as.vector(AdjMat_DiagSame %*% sigmayy[, prev])
        muxchange <- beta * (muxtimesmuy + diag_sxy) + phi * (muxtimesadjmuy + AMDXY)
        
        # ----- #
        
        # update mux
        duxdt <- -muxchange
        mux[, t] <- muxprev + as.vector(duxdt)
        mux[mux[, t] < 0, t] <- 0
        
        # update muy
        duydt <- muxchange - gamma * muyprev
        muy[, t] <- muyprev + as.vector(duydt)
        muy[muy[, t] < 0, t] <- 0
        
        # update diagonal sigmaxx
        dsxxdt <- muxchange + beta * (- 2 * muytimesdiagsxx - 2 * muxtimesdiagsxy) +
            phi * (- 2 * adjmuytimesdiagsxx - 2 * muxprev * AMDXY)
        sigmaxx[SEQ_DIAG_XX_XY, t] <- diag_sxx + as.vector(dsxxdt)
        sigmaxx[SEQ_DIAG_XX_XY, t] <- ifelse(sigmaxx[SEQ_DIAG_XX_XY, t] < 0, 0, sigmaxx[SEQ_DIAG_XX_XY, t])
        
        # update diagonal sigmayy
        dsyydt <- muxchange + beta * (2 * (muxprev * diag_syy + muyprev * diag_sxy)) +
            phi * (2 * adj_muy * diag_sxy + 2 * muxAMDSsyy) +
            gamma * (muyprev - 2 * diag_syy)
        sigmayy[SEQ_DIAG_XX_XY, t] <- diag_syy + as.vector(dsyydt)
        sigmaxx[SEQ_DIAG_XX_XY, t] <- ifelse(sigmaxx[SEQ_DIAG_XX_XY, t] < 0, 0, sigmaxx[SEQ_DIAG_XX_XY, t])
        
        # update diagonal sigmaxy
        dsxydt <- -muxchange + beta * (muytimesdiagsxx + muxtimesdiagsxy - muxprev * diag_syy - muyprev * diag_sxy) +
            phi * (adjmuytimesdiagsxx + muxprev * AMDXY - muxAMDSsyy - adj_muy * diag_sxy) -
            gamma * diag_sxy
        sigmaxy[SEQ_DIAG_XY, t] <- diag_sxy + as.vector(dsxydt)
        
        # ------- #
        
        sxyS4 <- sigmaxy[SEQ4, prev]
        muxS1 <- mux[SEQ1, prev]
        muxS2 <- mux[SEQ2, prev]
        sxyS3 <- sigmaxy[SEQ3, prev]
        muyS1 <- muy[SEQ1, prev]
        sxxS5 <- sigmaxx[SEQ5, prev]
        syyS5 <- sigmayy[SEQ5, prev]
        sxyS6 <- sigmaxy[SEQ6, prev]
        muySEQ2 <- muy[SEQ2, prev]
        muxSEQ7 <- mux[SEQ7, prev]
        muxSEQ8 <- mux[SEQ8, prev]
        sigmaxxSEQ9 <- sigmaxx[SEQ9, prev]
        AM2muy <- as.vector(AdjMat2 %*% muyprev)
        AM4muy <- as.vector(AdjMat4 %*% muyprev)
        
        AM1xy <- as.vector(AdjMat1 %*% sigmaxyprev)
        AM3xy <- as.vector(AdjMat3 %*% sigmaxyprev)
        AM5yy <- as.vector(AdjMat5 %*% sigmayyprev)
        AM6yy <- as.vector(AdjMat6 %*% sigmayyprev)
        AM7yy <- as.vector(AdjMat7 %*% sigmayyprev)
        AM8muy <- as.vector(AdjMat8 %*% muyprev)
        AM9xy <- as.vector(AdjMat9 %*% sigmaxyprev)
        AM10muy <- as.vector(AdjMat10 %*% muyprev)
        
        sigmayyS9 <- sigmayy[SEQ9, prev]
        muyS7 <- muy[SEQ7, prev]
        muyS8 <- muy[SEQ8, prev]
        
        # update off-diagonal sigmaxx
        dsigmaxxod_dt <- betaSEQ1 * (-muxS1 * sxyS3 - muyS1 * sxxS5) +
            betaSEQ2 * (-muxS2 * sxyS4 - muySEQ2 * sxxS5) +
            phiSEQ1 * (-muxS1 * AM1xy - AM2muy * sxxS5) +
            phiSEQ2 * (-muxS2 * AM3xy - AM4muy * sxxS5)
        sigmaxx[SEQ5, t] <- as.vector(sxxS5 + dsigmaxxod_dt)
        
        # update off-diagonal sigmayy
        dsigmayyod_dt <- betaSEQ1 * (muxS1 * syyS5 + muyS1 * sxyS4) +
            betaSEQ2 * (muxS2 * syyS5 + muySEQ2 * sxyS3) +
            phiSEQ1 * (muxS1 * AM5yy + AM2muy * sxyS4) +
            phiSEQ2 * (muxS2 * AM6yy + AM4muy * sxyS3) +
            -2 * gamma * syyS5
        sigmayy[SEQ5, t] <- as.vector(syyS5 + dsigmayyod_dt)
        
        # update off-diagonal sigmaxy
        dsigmaxyod_dt <- betaSEQ7 * (-muxSEQ7 * sigmayyS9 - muyS7 * sxyS6) +
            betaSEQ8 * (muxSEQ8 * sxyS6 + muyS8 * sigmaxxSEQ9) +
            phiSEQ7 * (-muxSEQ7 * AM7yy - AM8muy * sxyS6) +
            phiSEQ8 * (muxSEQ8 * AM9xy + AM10muy * sigmaxxSEQ9) +
            -1 * gamma * sxyS6
        sigmaxy[SEQ6, t] <- as.vector(sxyS6 + dsigmaxyod_dt)
        
    } # end of main t loop
    
    # managing the time indices
    this_index <- seq(1, nt, by = num_breaks)
    stopifnot(length(this_index) == nt_orig)
    mux <- mux[, this_index]
    muy <- muy[, this_index]
    sigmaxx <- sigmaxx[, this_index]
    sigmaxy <- sigmaxy[, this_index]
    sigmayy <- sigmayy[, this_index]
    
    if (return_all) {
        return(list(mux = mux, muy = muy, sigmaxx = sigmaxx, sigmaxy = sigmaxy, sigmayy = sigmayy))
    } else {
        return(list(mux = mux, sigmaxx = sigmaxx))
    }
    
}



# eof