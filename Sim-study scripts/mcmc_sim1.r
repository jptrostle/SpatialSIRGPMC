# Same script with p = 0.75: uncomment line in update_lambda function
# Same exact code with nu = 10, no changes needed

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
library(parallel)
RNGkind("L'Ecuyer-CMRG")
set.seed(560788664)
# set as necessary
setwd("..")


# --------- #
# (1) Setup #
# --------- #

this_run_code <- "sim1"
this_em_code <- "sim1"
this_save_code <- "sim1"

ns <- 25
start_day <- 61
nt <- 80
K <- 10000
lnt <- nt + 1

burnin <- 5000
nsamp <- 5000
nruns <- burnin + nsamp

# SVD kriging nearest-neighbor count
muxnnc <- 10
Sxxnnc <- 10

# SVD kriging for means: how many basis functions to use?
muxsbfc <- 20
muxtbfc <- 10

# SVD kriging for covariances: how many basis functions to use?
Sxxsbfc <- 10 # spatial basis-function count
Sxxtbfc <- 10 # temporal basis-function count

alphabf <- 5


# ------------- #
# (2) Functions #
# ------------- #

# emulator function for the means
run_mux_em <- function(theta0, Theta, nnc, nsbf, ntbf, M, U1, U2, tm, ts) {
    
    theta0[1:2] <- (theta0[1:2] - tm) / ts
    theta0 <- matrix(theta0, nrow = 1)
    nn <- c(get.knnx(data = Theta, query = theta0, k = nnc)$nn.index)
    ThetaR <- rbind(Theta[nn, ], theta0)
    distances <- spDists(ThetaR)
    Corr <- matern(distances, 1, 2.5)
    
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
    return(mux0)
    
}

# emulator function for the covariances
run_Sxx_em <- function(theta0, Theta, nnc, nsbf, ntbf, M, U1, U3, lnt, tm, ts) {
    
    theta0[1:2] <- (theta0[1:2] - tm) / ts
    theta0 <- matrix(theta0, nrow = 1)
    nn <- c(get.knnx(data = Theta, query = theta0, k = nnc)$nn.index)
    ThetaR <- rbind(Theta[nn, ], theta0)
    distances <- spDists(ThetaR)
    Corr <- matern(distances, 1, 2.5)
    
    # setting up kriging details
    Sigma11Inv <- solve(Corr[1:nnc, 1:nnc])
    Sigma12 <- Corr[1:nnc, (nnc + 1)]
    j <- rep(1, nnc)
    tjS11Inv <- t(j) %*% Sigma11Inv
    S12TS11Inv <- t(Sigma12) %*% Sigma11Inv
    this_denom <- c(tjS11Inv %*% j)
    weights <- tjS11Inv / this_denom + S12TS11Inv - c(S12TS11Inv %*% j) * tjS11Inv / this_denom
    
    ChMpred <- array(0, dim = c(nsbf, nsbf, ntbf))
    for (si1 in 1:nsbf) {
        for (si2 in 1:nsbf) {
            for (ti in 1:ntbf) {
                ChMpred[si1, si2, ti] <- weights %*% M[si1, si2, ti, nn]
            }
        }
    }
    
    Chpred <- ttm(as.tensor(ChMpred), U3, m = 3)
    Chpred <- Chpred@data
    
    Sbf <- list()
    for (t in 1:lnt) Sbf[[t]] <- U1 %*% Chpred[,, t]
    
    return(Sbf)
    
}

update_lambda <- function(mux, Sbf, alpha, ns, lnt) {
    
    X <- matrix(0, nrow = ns, ncol = lnt)
    mux <- matrix(mux, nrow = ns, ncol = lnt)
    
    for (t in 1:lnt) X[, t] <- mux[, t] + Sbf[[t]] %*% alpha[, t]
    
    lambda <- X[, 1:(lnt - 1)] - X[, 2:lnt]
    #lambda <- 0.75 * lambda
    lambda[lambda < 0.00001] <- 0.00001
    
    return(lambda)
    
}

# main MCMC code
run_mcmc <- function(this_run_id, sim_data_path, sim_results_path, alphabf, burnin, nsamp, ns, nt, start_day, K, Theta, muxnnc, Sxxnnc, muxsbfc, muxtbfc, Sxxsbfc, Sxxtbfc, muxSBF, muxTBF, SxxSBF, SxxTBF, muxM, SxxM) {
    
    Theta_m <- colMeans(Theta[, 1:2])
    Theta_s <- apply(Theta[, 1:2], 2, sd)
    Theta[, 1] <- (Theta[, 1] - Theta_m[1]) / Theta_s[1]
    Theta[, 2] <- (Theta[, 2] - Theta_m[2]) / Theta_s[2]
    
    possible_s <- unique(Theta[, 3:4])
    
    lnt <- nt + 1
    library(splines)
    library(MASS)
    B <- bs(1:lnt, df = alphabf, intercept = TRUE)
    num_std <- diag(diag(B %*% t(B))^-0.5)
    B <- num_std %*% B
    Ba_to_alpha <- function(B, alpha) return(t(B %*% t(alpha)))
    
    y <- readRDS(paste(sim_data_path, "/", this_run_id, ".rds", sep = ""))
    y <- y[, start_day:dim(y)[2]]
    nruns <- burnin + nsamp    
    n <- ns * nt
    
    # buckets for results
    theta <- matrix(0.03, nrow = 2, ncol = nruns)
    nu <- rep(0, nruns)
    s <- matrix(0, nrow = 2, ncol = nruns)
    alpha <- array(0, dim = c(Sxxsbfc, alphabf, nruns)) # bsplines change
    
    # starting values (I change theta later)
    nu[1] <- 2
    theta[1, 1] <- 0.03
    theta[2, 1] <- 0.03
    s[1, 1] <- 2
    s[2, 1] <- 3
    
    # proposals and acceptance rates
    theta_accepts1 <- theta_accepts2 <- s_accepts <- nu_accepts <- 0
    alpha_accepts <- matrix(0, nrow = Sxxsbfc, ncol = alphabf)
    
    # lots of tuning to get starting BP_V for DRAM updates
    BP_V <- matrix(c(1.239760e-05, -3.604449e-06, 8.419056e-05, 0.0001131931, 0.0002180501, -1.781229e-04, -1.961370e-04,
    -3.604449e-06, 1.093984e-06, -5.670947e-05, -0.0001020023, -0.0001489973, -5.781553e-05, 3.323367e-05,
    8.419056e-05, -5.670947e-05, 4.976149e-02, 0.0439270072, 0.0724490666, 7.008686e-02, 3.111200e-02,
    1.131931e-04, -1.020023e-04, 4.392701e-02, 0.1095616917, 0.1243737513, 1.669032e-01, 2.656932e-02,
    2.180501e-04, -1.489973e-04, 7.244907e-02, 0.1243737513, 0.1777368481, 1.913039e-01, 5.862192e-02,
    -1.781229e-04, -5.781553e-05, 7.008686e-02, 0.1669031908, 0.1913038659, 2.836215e-01, 6.412404e-02,
    -1.961370e-04, 3.323367e-05, 3.111200e-02, 0.0265693160, 0.0586219221, 6.412404e-02, 2.244095e-01), nrow = 7, byrow = TRUE)
    BP_V <- (BP_V + t(BP_V)) / 2
    BP_V <- 3.5 * BP_V
    tC_BP_V <- t(chol(BP_V))
    
    # DRAM stuff for alpha
    alpha_sds <- matrix(1, nrow = Sxxsbfc, ncol = alphabf)
    alpha_sds[1, ] <- -3 # this is a placeholder to make sure no problems; it shouldn't ever be used
    
    nu_sd <- 1
    
    # the rest of the precalcs that need to be done now
    mux <- run_mux_em(c(theta[1:2, 1], s[, 1]), Theta, muxnnc, muxsbfc, muxtbfc, muxM, muxSBF, muxTBF, Theta_m, Theta_s)
    Sbf <- run_Sxx_em(c(theta[1:2, 1], s[, 1]), Theta, Sxxnnc, Sxxsbfc, Sxxtbfc, SxxM, SxxSBF, SxxTBF, lnt, Theta_m, Theta_s)
    lambda <- update_lambda(mux, Sbf, Ba_to_alpha(B, alpha[,, 1]), ns, lnt)
    
    # ******************************************************* #
    
    # ------------------------------- #
    # (5) Pick better starting values #
    # ------------------------------- #
    
    best_mse <- sum(y^2)
    pass_accepts <- rep(0, 3)
    for (j in 1:1000) {
        
        this_beta <- runif(1, 0.037, 0.045)
        this_phi <- runif(1, 0.02, 0.03)
        this_s <- sample(1:5, 1)
        this_s1 <- possible_s[this_s, 1]
        this_s2 <- possible_s[this_s, 2]
        
        muxmaybe <- run_mux_em(c(this_beta, this_phi, this_s1, this_s2), Theta, muxnnc, muxsbfc, muxtbfc, muxM, muxSBF, muxTBF, Theta_m, Theta_s)
        Sbfmaybe <- run_Sxx_em(c(this_beta, this_phi, this_s1, this_s2), Theta, Sxxnnc, Sxxsbfc, Sxxtbfc, SxxM, SxxSBF, SxxTBF, lnt, Theta_m, Theta_s)
        lmaybe <- update_lambda(muxmaybe, Sbfmaybe, Ba_to_alpha(B, alpha[,, 1]), ns, lnt)
        
        this_mse <- sum((y - lmaybe)^2)
        
        if (this_mse < best_mse) {
            
            best_mse <- this_mse
            theta[, 1] <- c(this_beta, this_phi)
            s[, 1] <- c(this_s1, this_s2)
            mux <- muxmaybe
            Sbf <- Sbfmaybe
            lambda <- lmaybe
            
            pass_accepts[1] <- pass_accepts[1] + 1
            
        }
    }
    
    # second pass with S fixed
    for (j in 1:1000) {
        
        this_beta <- runif(1, 0.037, 0.045)
        this_phi <- runif(1, 0.02, 0.03)
        
        muxmaybe <- run_mux_em(c(this_beta, this_phi, s[1, 1], s[2, 1]), Theta, muxnnc, muxsbfc, muxtbfc, muxM, muxSBF, muxTBF, Theta_m, Theta_s)
        Sbfmaybe <- run_Sxx_em(c(this_beta, this_phi, s[1, 1], s[2, 1]), Theta, Sxxnnc, Sxxsbfc, Sxxtbfc, SxxM, SxxSBF, SxxTBF, lnt, Theta_m, Theta_s)
        lmaybe <- update_lambda(muxmaybe, Sbfmaybe, Ba_to_alpha(B, alpha[,, 1]), ns, lnt)
        
        this_mse <- sum((y - lmaybe)^2)
        
        if (this_mse < best_mse) {
            
            best_mse <- this_mse
            theta[, 1] <- c(this_beta, this_phi)
            mux <- muxmaybe
            Sbf <- Sbfmaybe
            lambda <- lmaybe
            
            pass_accepts[2] <- pass_accepts[2] + 1
            
        }
        
    }
    
    # now check alphas
    for (pass in 1:10) {
        for (a in 1:Sxxsbfc) {
            for (b in 1:alphabf) {
                for (hmm in seq(-2, 2, 0.05)) {
                    
                    alphamaybe <- alpha[,, 1]
                    alphamaybe[a, b] <- hmm
                    lmaybe <- update_lambda(mux, Sbf, Ba_to_alpha(B, alphamaybe), ns, lnt)
                    this_mse <- sum((y - lmaybe)^2)
                    
                    if (this_mse < best_mse) {
                        
                        best_mse <- this_mse
                        alpha[a, b, 1] <- hmm
                        lambda <- lmaybe
                        
                    }
                }   
            }
        }
    }
    
    # third pass by searching near the last accepted spot
    for (j in 1:1000) {
        
        this_beta <- rnorm(1, theta[1, 1], 0.0001)
        this_phi <- rnorm(1, theta[2, 1], 0.0001)
        
        muxmaybe <- run_mux_em(c(this_beta, this_phi, s[1, 1], s[2, 1]), Theta, muxnnc, muxsbfc, muxtbfc, muxM, muxSBF, muxTBF, Theta_m, Theta_s)
        Sbfmaybe <- run_Sxx_em(c(this_beta, this_phi, s[1, 1], s[2, 1]), Theta, Sxxnnc, Sxxsbfc, Sxxtbfc, SxxM, SxxSBF, SxxTBF, lnt, Theta_m, Theta_s)
        lmaybe <- update_lambda(muxmaybe, Sbfmaybe, Ba_to_alpha(B, alpha[,, 1]), ns, lnt)
        
        this_mse <- sum((y - lmaybe)^2)
        
        if (this_mse < best_mse) {
            
            best_mse <- this_mse
            theta[, 1] <- c(this_beta, this_phi)
            mux <- muxmaybe
            Sbf <- Sbfmaybe
            lambda <- lmaybe
            
            pass_accepts[3] <- pass_accepts[3] + 1
            
        }
        
    }
    
    # now check alphas again
    for (pass in 1:10) {
        for (a in 1:Sxxsbfc) {
            for (b in 1:alphabf) {
                for (hmm in seq(-2, 2, 0.05)) {
                    
                    alphamaybe <- alpha[,, 1]
                    alphamaybe[a, b] <- hmm
                    lmaybe <- update_lambda(mux, Sbf, Ba_to_alpha(B, alphamaybe), ns, lnt)
                    this_mse <- sum((y - lmaybe)^2)
                    
                    if (this_mse < best_mse) {
                        
                        best_mse <- this_mse
                        alpha[a, b, 1] <- hmm
                        lambda <- lmaybe
                        
                    }
                }   
            }
        }
    }
    
    
    # ------------- #
    # (6) Main MCMC #
    # ------------- #
    
    t1 <- Sys.time()
    for (i in 2:nruns) {
        
        if (i %% 100 == 0) print(paste("On iter", i, "at", Sys.time()))
        
        prev <- i - 1
        
        # ----------------- #
        # beta, phi, alpha1 #
        # ----------------- #
        
        first_stage_accepted <- FALSE
        
        cand_vec <- as.vector(tC_BP_V %*% rnorm(7) + c(theta[1:2, prev], alpha[1,, prev]))
        cand_theta <- c(cand_vec[1:2], s[, prev])
        alphamaybe <- alpha[,, prev]
        alphamaybe[1, ] <- cand_vec[3:7]
        muxmaybe <- run_mux_em(cand_theta, Theta, muxnnc, muxsbfc, muxtbfc, muxM, muxSBF, muxTBF, Theta_m, Theta_s)
        Sbfmaybe <- run_Sxx_em(cand_theta, Theta, Sxxnnc, Sxxsbfc, Sxxtbfc, SxxM, SxxSBF, SxxTBF, lnt, Theta_m, Theta_s)
        lmaybe <- update_lambda(muxmaybe, Sbfmaybe,  Ba_to_alpha(B, alphamaybe), ns, lnt)
        
        log_orig <- sum(dnbinom(y, mu = lambda, size = lambda / (nu[prev] - 1), log = TRUE)) + sum(dlnorm(theta[1:2, prev], -2, 1, log = TRUE)) + sum(dnorm(alpha[1, , prev], mean = 0, sd = 1, log = TRUE))
        log_cand <- sum(dnbinom(y, mu = lmaybe, size = lmaybe / (nu[prev] - 1), log = TRUE)) + sum(dlnorm(cand_theta[1:2], -2, 1, log = TRUE)) + sum(dnorm(alphamaybe[1, ], mean = 0, sd = 1, log = TRUE))
        
        log_accep_prob <- log_cand - log_orig
        if (log(runif(1)) < log_accep_prob) {
            theta_accepts1 <- theta_accepts1 + 1
            theta[1:2, i] <- cand_theta[1:2]
            lambda <- lmaybe
            mux <- muxmaybe
            Sbf <- Sbfmaybe
            alpha[1,, i] <- alphamaybe[1, ]
            alpha_accepts[1, ] <- alpha_accepts[1, ] + 1
            first_stage_accepted <- TRUE
        } 
        
        if (!first_stage_accepted) {
            
            cand_vec2 <- as.vector(0.2 * tC_BP_V %*% rnorm(7) + c(theta[1:2, prev], alpha[1, , prev]))
            cand_theta2 <- c(cand_vec2[1:2], s[, prev])
            alphamaybe2 <- alpha[,, prev]
            alphamaybe2[1, ] <- cand_vec2[3:7]
            muxmaybe2 <- run_mux_em(cand_theta2, Theta, muxnnc, muxsbfc, muxtbfc, muxM, muxSBF, muxTBF, Theta_m, Theta_s)
            Sbfmaybe2 <- run_Sxx_em(cand_theta2, Theta, Sxxnnc, Sxxsbfc, Sxxtbfc, SxxM, SxxSBF, SxxTBF, lnt, Theta_m, Theta_s)
            lmaybe2 <- update_lambda(muxmaybe2, Sbfmaybe2, Ba_to_alpha(B, alphamaybe2), ns, lnt)
            
            log_cand2 <- sum(dnbinom(y, mu = lmaybe2, size = lmaybe2 / (nu[prev] - 1), log = TRUE)) + sum(dlnorm(cand_theta2[1:2], -2, 1, log = TRUE)) + sum(dnorm(alphamaybe2[1, ], mean = 0, sd = 1, log = TRUE))
            
            top_of_accep2 <- log_cand2 + dmvnorm(cand_vec, mean = cand_vec2, sigma = BP_V, log = TRUE) + (1 - exp(log_cand - log_cand2))
            bottom_of_accep2 <- log_orig + dmvnorm(cand_vec, mean = c(theta[1:2, prev], alpha[1, , prev]), sigma = BP_V, log = TRUE) + (1 - exp(log_accep_prob))
            log_accep_prob2 <- top_of_accep2 - bottom_of_accep2
            if (log(runif(1)) < log_accep_prob2) {
                
                theta_accepts2 <- theta_accepts2 + 1
                theta[1:2, i] <- cand_theta2[1:2]
                lambda <- lmaybe2
                mux <- muxmaybe2
                Sbf <- Sbfmaybe2
                alpha[1,, i] <- alphamaybe2[1, ]
                alpha_accepts[1, ] <- alpha_accepts[1, ] + 1
                
            } else {
                
                theta[1:2, i] <- theta[1:2, prev]
                alpha[1,, i] <- alpha[1,, prev]
                
            }
            
        }
        
        if (i %% 200 == 0 && i <= burnin && alpha_accepts[1] > 10) {
            adj_index <- 1:i
            BP_V <- 3.5 * cov(t(rbind(theta[, adj_index], alpha[1,, adj_index])))
            tC_BP_V <- t(chol(BP_V))
            
        }        
        
        # ---------- #
        # S1, S2, nu #
        # ---------- #
        
        cand_s_ind <- sample(1:5, 1)
        cand_s1 <- possible_s[cand_s_ind, 1]
        cand_s2 <- possible_s[cand_s_ind, 2]
        cand_theta <- c(theta[1:2, i], cand_s1, cand_s2)
        muxmaybe <- run_mux_em(cand_theta, Theta, muxnnc, muxsbfc, muxtbfc, muxM, muxSBF, muxTBF, Theta_m, Theta_s)
        Sbfmaybe <- run_Sxx_em(cand_theta, Theta, Sxxnnc, Sxxsbfc, Sxxtbfc, SxxM, SxxSBF, SxxTBF, lnt, Theta_m, Theta_s)
        this_alpha <- alpha[,, prev]
        this_alpha[1, ] <- alpha[1,, i]
        lmaybe <- update_lambda(muxmaybe, Sbfmaybe, Ba_to_alpha(B, this_alpha), ns, lnt)
        log_orig <- sum(dnbinom(y, mu = lambda, size = lambda / (nu[prev] - 1), log = TRUE))
        log_cand <- sum(dnbinom(y, mu = lmaybe, size = lmaybe / (nu[prev] - 1), log = TRUE))
        log_accep_prob <- log_cand - log_orig
        if (log(runif(1)) < log_accep_prob) {
            s[, i] <- c(cand_s1, cand_s2)
            lambda <- lmaybe
            mux <- muxmaybe
            Sbf <- Sbfmaybe
            s_accepts <- s_accepts + 1
        } else {
            s[, i] <- s[, prev]
        }
        
        # update nu
        cand_nu <- rtruncnorm(1, mean = nu[prev], sd = nu_sd, a = 1)
        log_orig <- sum(dnbinom(y, mu = lambda, size = lambda / (nu[prev] - 1), log = TRUE)) + log(dtruncnorm(nu[prev], mean = 3, sd = 5, a = 1))
        log_cand <- sum(dnbinom(y, mu = lambda, size = lambda / (cand_nu - 1), log = TRUE)) + log(dtruncnorm(cand_nu, mean = 3, sd = 5, a = 1))
        log_trans <- log(dtruncnorm(cand_nu, mean = nu[prev], sd = nu_sd, a = 1))
        log_back <- log(dtruncnorm(nu[prev], mean = cand_nu, sd = nu_sd, a = 1))
        log_accep_prob <- log_cand - log_orig + log_back - log_trans
        if (log(runif(1)) < log_accep_prob) {
            nu[i] <- cand_nu
            nu_accepts <- nu_accepts + 1
        } else {
            nu[i] <- nu[prev]
        }
        
        # alpha
        for (a in 2:Sxxsbfc) {
            for (b in 1:alphabf) {
                if (a == 1 && b == 1) next
                
                alphamaybe <- alpha[,, prev]
                alphamaybe[a, b] <- rnorm(1, mean = alpha[a, b, prev], sd = alpha_sds[a, b])
                if (a > 1) alphamaybe[1:(a - 1), ] <- alpha[1:(a - 1),, i]
                if (b > 1) alphamaybe[a, 1:(b - 1)] <- alpha[a, 1:(b - 1), i]
                lmaybe <- update_lambda(mux, Sbf, Ba_to_alpha(B, alphamaybe), ns, lnt)
                log_orig <- sum(dnbinom(y, mu = lambda, size = lambda / (nu[i] - 1), log = TRUE)) + dnorm(alpha[a, b, prev], mean = 0, sd = 1, log = TRUE)
                log_cand <- sum(dnbinom(y, mu = lmaybe, size = lmaybe / (nu[i] - 1), log = TRUE)) + dnorm(alphamaybe[a, b], mean = 0, sd = 1, log = TRUE)
                log_accep_prob <- log_cand - log_orig
                if (log(runif(1)) < log_accep_prob) {
                    first_stage_accepted <- 1
                    alpha_accepts[a, b] <- alpha_accepts[a, b] + 1
                    alpha[a, b, i] <- alphamaybe[a, b]
                    lambda <- lmaybe
                } else {
                    alpha[a, b, i] <- alpha[a, b, prev]
                }
            }
        }
        
        
        # ------------------------------- #
        # bookkeeping, tkd sd adjustments #
        # ------------------------------- #
        
        every_num <- 400
        if (i <= burnin && i %% every_num == 0 && i >= every_num) {
            
            # update nu sd
            if (nu_accepts / every_num < 0.2) nu_sd <- nu_sd * 0.9
            if (nu_accepts / every_num > 0.25) nu_sd <- nu_sd * 1.1
            
            for (a in 1:Sxxsbfc) {
                for (b in 1:alphabf) {
                    if (a == 1 && b == 1) next
                    if (alpha_accepts[a, b] / every_num < 0.2) alpha_sds[a, b] <- alpha_sds[a, b] * 0.9
                    if (alpha_accepts[a, b] / every_num > 0.25) alpha_sds[a, b] <- alpha_sds[a, b] * 1.1   
                }
            }
            
            nu_accepts <- 0
            alpha_accepts <- matrix(0, nrow = Sxxsbfc, ncol = alphabf)
            
        }
        
        if (i == burnin) {
            theta_accepts1 <- theta_accepts2 <- s_accepts <- nu_accepts <- 0
            alpha_accepts <- matrix(0, nrow = Sxxsbfc, ncol = alphabf)
            eta_accepts <- 0
        }
        
    }
    t2 <- Sys.time()
    
    results_list <- list(theta = theta, alpha = alpha, nu = nu, s = s,
                         aar = alpha_accepts / nsamp, tar1 = theta_accepts1 / nsamp, tar2 = theta_accepts2 / nsamp,
                         nar = nu_accepts / nsamp, sar = s_accepts / nsamp,
                         rt = difftime(t2, t1, units = "secs"))
    saveRDS(results_list, paste(sim_results_path, "/Results ", this_run_id, ".rds", sep = ""))
    return(1)
    
}


# -------- #
# (3) Prep #
# -------- #

em_data <- paste("Emulator data/", this_em_code, sep = "")
Theta <- readRDS(paste(em_data, "/ClTheta.rds", sep = ""))
muxM <- readRDS(paste(em_data, "/muxM.rds", sep = ""))
muxSBF <- readRDS(paste(em_data, "/muxSBF.rds", sep = ""))
muxTBF <- readRDS(paste(em_data, "/muxTBF.rds", sep = ""))
SxxM <- readRDS(paste(em_data, "/SxxM.rds", sep = ""))
SxxSBF <- readRDS(paste(em_data, "/SxxSBF.rds", sep = ""))
SxxTBF <- readRDS(paste(em_data, "/SxxTBF.rds", sep = ""))


# ----------- # 
# (4) On y va #
# ----------- # 

nsims <- 200
sim_data_path <- paste("Simulated data/", this_run_code, sep = "")
sim_results_path <- paste("Results/", this_save_code, sep = "")
results_codes <- mclapply(1:nsims, function(x) {run_mcmc(x, sim_data_path, sim_results_path, alphabf, burnin, nsamp, ns, nt, start_day, K, Theta, muxnnc, Sxxnnc, muxsbfc, muxtbfc, Sxxsbfc, Sxxtbfc, muxSBF, muxTBF, SxxSBF, SxxTBF, muxM, SxxM)}, mc.cores = 10)
saveRDS(results_codes, paste("Results/", this_save_code, "/results_codes.rds", sep = ""))


# ----------- #
# (5) Results #
# ----------- #

beta_95_cov <- beta_99_cov <- phi_95_cov <- phi_99_cov <- beta_mse <- phi_mse <- 0
not_in_b95 <- not_in_b99 <- not_in_p95 <- not_in_p99 <- c()

for (i in 1:nsims) {
    
    mcmc <- readRDS(paste("Results/", this_save_code, "/Results ", i, ".rds", sep = ""))
    bq <- quantile(mcmc$theta[1, (burnin + 1):nruns], c(0.005, 0.025, 0.975, 0.995))
    pq <- quantile(mcmc$theta[2, (burnin + 1):nruns], c(0.005, 0.025, 0.975, 0.995))
    
    beta_mse <- beta_mse + (mean(mcmc$theta[1, (burnin + 1):nruns]) - 0.043)^2 / nsims
    phi_mse <- phi_mse + (mean(mcmc$theta[2, (burnin + 1):nruns]) - 0.025)^2 / nsims

    if (0.043 > bq[1] && 0.043 < bq[4]) {
        beta_99_cov <- beta_99_cov + 1
    } else {
        not_in_b99 <- c(not_in_b99, i)
    }
    
    if (0.043 > bq[2] && 0.043 < bq[3]) {
        beta_95_cov <- beta_95_cov + 1
    } else {
        not_in_b95 <- c(not_in_b95, i)
    }
    
    if (0.025 > pq[1] && 0.025 < pq[4]) {
        phi_99_cov <- phi_99_cov + 1
    } else {
        not_in_p99 <- c(not_in_p99, i)
    }
    
    if (0.025 > pq[2] && 0.025 < pq[3]) {
        phi_95_cov <- phi_95_cov + 1
    } else {
        not_in_p95 <- c(not_in_p95, i)
    }
    
}

conc_list <- list(b95 = beta_95_cov / nsims, p95 = phi_95_cov / nsims,
                  b99 = beta_99_cov / nsims, p99 = phi_99_cov / nsims,
                  not_in_b95 = not_in_b95, not_in_p95 = not_in_p95,
                  not_in_b99 = not_in_b99, not_in_p99 = not_in_p99,
                  beta_mse = beta_mse, phi_mse = phi_mse)
saveRDS(conc_list, paste("Results/", this_save_code, "/Final results.rds", sep = ""))

# eof