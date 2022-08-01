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

this_run_code <- "sim2"
this_save_code <- "sim2_meanonly"

ur_p <- 1.00
ns <- 25
start_day <- 21
nt <- 80
K <- 60000
lnt <- nt + 1

burnin <- 50000
nsamp <- 10000
nruns <- burnin + nsamp

# SVD kriging nearest-neighbor count
muxnnc <- 10
Sxxnnc <- 10

# SVD kriging for means: how many basis functions to use?
muxsbfc <- 20
muxtbfc <- 10

# SVD kriging for covariances: how many basis functions to use?
Sxxsbfc <- 10
Sxxtbfc <- 10


# ------------- #
# (2) Functions #
# ------------- #

run_mux_em <- function(theta0, Theta, nnc, nsbf, ntbf, M, U1, U2, tm, ts) {
    
    theta0[1:3] <- (theta0[1:3] - tm) / ts
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

update_lambda <- function(mux, ns, lnt) {
    
    mux <- matrix(mux, nrow = ns, ncol = lnt)
    lambda <- mux[, 1:(lnt - 1)] - mux[, 2:lnt]
    lambda[lambda < 0.00001] <- 0.00001
    
    return(lambda)
    
}

# main MCMC code
run_mcmc <- function(this_run_id, sim_data_path, sim_results_path, alphabf, burnin, nsamp, ns, nt, start_day, K, Theta, muxnnc, Sxxnnc, muxsbfc, muxtbfc, Sxxsbfc, Sxxtbfc, muxSBF, muxTBF, SxxSBF, SxxTBF, muxM, SxxM, ur_p) {
    
    Theta_m <- colMeans(Theta[, 1:3])
    Theta_s <- apply(Theta[, 1:3], 2, sd)
    Theta[, 1] <- (Theta[, 1] - Theta_m[1]) / Theta_s[1]
    Theta[, 2] <- (Theta[, 2] - Theta_m[2]) / Theta_s[2]
    Theta[, 3] <- (Theta[, 3] - Theta_m[3]) / Theta_s[3]
    
    possible_s <- unique(Theta[, 4:5])
    
    lnt <- nt + 1
    
    y <- readRDS(paste(sim_data_path, "/", this_run_id, ".rds", sep = ""))
    y <- y[, start_day:dim(y)[2]]
    y[is.nan(y)] <- 0
    nruns <- burnin + nsamp    
    n <- ns * nt
    
    # buckets for results
    theta <- matrix(0.05, nrow = 3, ncol = nruns)
    nu <- rep(0, nruns)
    s <- matrix(0, nrow = 2, ncol = nruns)
    
    # starting values (I change theta later)
    nu[1] <- 2
    theta[1, 1] <- -2.5
    theta[2, 1] <- 0.007
    theta[3, 1] <- 0.04
    
    s[1, ] <- 2
    s[2, ] <- 3
    
    # proposals and acceptance rates
    theta_accepts1 <- theta_accepts2 <- s_accepts <- nu_accepts <- 0
    
    BP_V <- diag(3)
    BP_V[1, 1] <- 1e-3
    BP_V[2, 2] <- 1e-5
    BP_V[3, 3] <- 1e-6
    BP_V <- (BP_V + t(BP_V)) / 2
    tC_BP_V <- t(chol(BP_V))
    
    # remaining TKD stuff
    nu_sd <- 1
    
    # the rest of the precalcs that need to be done now
    mux <- run_mux_em(c(theta[, 1], s[, 1]), Theta, muxnnc, muxsbfc, muxtbfc, muxM, muxSBF, muxTBF, Theta_m, Theta_s)
    lambda <- update_lambda(mux, ns, lnt)
    
    # ******************************************************* #
    
    # ------------------------------- #
    # (5) Pick better starting values #
    # ------------------------------- #
    
    best_mse <- sum(y^2)
    pass_accepts <- rep(0, 3)
    for (j in 1:5000) {
        
        this_beta0 <- runif(1, -3.5, -2.2)
        this_beta1 <- runif(1, 0.08, 0.12)
        this_phi <- runif(1, 0.04, 0.05)
        this_s <- sample(1:5, 1)
        this_s1 <- possible_s[this_s, 1]
        this_s2 <- possible_s[this_s, 2]
        
        muxmaybe <- run_mux_em(c(this_beta0, this_beta1, this_phi, this_s1, this_s2), Theta, muxnnc, muxsbfc, muxtbfc, muxM, muxSBF, muxTBF, Theta_m, Theta_s)
        lmaybe <- update_lambda(muxmaybe, ns, lnt)
        
        this_mse <- sum((y - lmaybe)^2)
        
        if (this_mse < best_mse) {
            
            best_mse <- this_mse
            theta[, 1] <- c(this_beta0, this_beta1, this_phi)
            s[, 1] <- c(this_s1, this_s2)
            mux <- muxmaybe
            lambda <- lmaybe
            
            pass_accepts[1] <- pass_accepts[1] + 1
            
        }
    }
    
    # second pass with S fixed
    for (j in 1:1000) {
        
        this_beta0 <- runif(1, -3.5, -2.2)
        this_beta1 <- runif(1, 0.08, 0.12)
        this_phi <- runif(1, 0.04, 0.05)
        
        muxmaybe <- run_mux_em(c(this_beta0, this_beta1, this_phi, s[1, 1], s[2, 1]), Theta, muxnnc, muxsbfc, muxtbfc, muxM, muxSBF, muxTBF, Theta_m, Theta_s)
        lmaybe <- update_lambda(muxmaybe, ns, lnt)
        
        this_mse <- sum((y - lmaybe)^2)
        
        if (this_mse < best_mse) {
            
            best_mse <- this_mse
            theta[, 1] <- c(this_beta0, this_beta1, this_phi)
            mux <- muxmaybe
            lambda <- lmaybe
            
            pass_accepts[2] <- pass_accepts[2] + 1
            
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
        
        cand_vec <- as.vector(tC_BP_V %*% rnorm(3) + theta[, prev])
        cand_theta <- c(cand_vec[1:3], s[, prev])
        
        muxmaybe <- run_mux_em(cand_theta, Theta, muxnnc, muxsbfc, muxtbfc, muxM, muxSBF, muxTBF, Theta_m, Theta_s)
        lmaybe <- update_lambda(muxmaybe, ns, lnt)
        
        log_orig <- sum(dnbinom(y, mu = lambda, size = lambda / (nu[prev] - 1), log = TRUE)) + dmvnorm(theta[, prev], rep(0, 3), diag(3) * 100, log = TRUE)
        log_cand <- sum(dnbinom(y, mu = lmaybe, size = lmaybe / (nu[prev] - 1), log = TRUE)) + dmvnorm(cand_theta[1:3], rep(0, 3), diag(3) * 100, log = TRUE)
        
        log_accep_prob <- log_cand - log_orig
        if (log(runif(1)) < log_accep_prob) {
            theta_accepts1 <- theta_accepts1 + 1
            theta[, i] <- cand_theta[1:3]
            lambda <- lmaybe
            mux <- muxmaybe
            first_stage_accepted <- TRUE
        } 
        
        if (!first_stage_accepted) {
            
            cand_vec2 <- as.vector(0.2 * tC_BP_V %*% rnorm(3) + theta[, prev])
            cand_theta2 <- c(cand_vec2[1:3], s[, prev])

            muxmaybe2 <- run_mux_em(cand_theta2, Theta, muxnnc, muxsbfc, muxtbfc, muxM, muxSBF, muxTBF, Theta_m, Theta_s)
            lmaybe2 <- update_lambda(muxmaybe2, ns, lnt)
            
            log_cand2 <- sum(dnbinom(y, mu = lmaybe2, size = lmaybe2 / (nu[prev] - 1), log = TRUE)) + dmvnorm(cand_theta2[1:3], rep(0, 3), diag(3) * 100, log = TRUE)
            
            top_of_accep2 <- log_cand2 + dmvnorm(cand_vec, mean = cand_vec2, sigma = BP_V, log = TRUE) + (1 - exp(log_cand - log_cand2))
            bottom_of_accep2 <- log_orig + dmvnorm(cand_vec, mean = theta[, prev], sigma = BP_V, log = TRUE) + (1 - exp(log_accep_prob))
            log_accep_prob2 <- top_of_accep2 - bottom_of_accep2
            if (log(runif(1)) < log_accep_prob2) {
                
                theta_accepts2 <- theta_accepts2 + 1
                theta[, i] <- cand_theta2[1:3]
                lambda <- lmaybe2
                mux <- muxmaybe2
                
            } else {
                
                theta[, i] <- theta[, prev]
                
            }
            
        }
        
        if (i %% 200 == 0 && i <= burnin) {
        
            if (i <= 2000) {
                adj_index <- 1:i
            } else {
                adj_index <- (i - 1999):i
            }
            
            BP_V <- 1.5 * cov(t(theta[, adj_index])) + diag(3) * 1e-8
            tC_BP_V <- t(chol(BP_V))
            
        }        
        
        # ----- #
        # s, nu #
        # ----- #
        
        cand_s_ind <- sample(1:5, 1)
        cand_s1 <- possible_s[cand_s_ind, 1]
        cand_s2 <- possible_s[cand_s_ind, 2]
        cand_theta <- c(theta[, i], cand_s1, cand_s2)
        muxmaybe <- run_mux_em(cand_theta, Theta, muxnnc, muxsbfc, muxtbfc, muxM, muxSBF, muxTBF, Theta_m, Theta_s)
        lmaybe <- update_lambda(muxmaybe, ns, lnt)
        log_orig <- sum(dnbinom(y, mu = lambda, size = lambda / (nu[prev] - 1), log = TRUE))
        log_cand <- sum(dnbinom(y, mu = lmaybe, size = lmaybe / (nu[prev] - 1), log = TRUE))
        log_accep_prob <- log_cand - log_orig
        if (log(runif(1)) < log_accep_prob) {
            s[, i] <- c(cand_s1, cand_s2)
            lambda <- lmaybe
            mux <- muxmaybe
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
        
        # ------------------------------- #
        # bookkeeping, tkd sd adjustments #
        # ------------------------------- #
        
        every_num <- 400
        if (i <= burnin && i %% every_num == 0 && i >= every_num) {
            
            # update nu sds
            if (nu_accepts / every_num < 0.2) nu_sd <- nu_sd * 0.9
            if (nu_accepts / every_num > 0.25) nu_sd <- nu_sd * 1.1
            
            nu_accepts <- 0
            
        }
        
        if (i == burnin) theta_accepts1 <- theta_accepts2 <- s_accepts <- nu_accepts <- 0
              
    }
    t2 <- Sys.time()
    
    
    results_list <- list(theta = theta, nu = nu, s = s,
                         tar1 = theta_accepts1 / nsamp, tar2 = theta_accepts2 / nsamp,
                         nar = nu_accepts / nsamp, sar = s_accepts / nsamp,
                         rt = difftime(t2, t1, units = "secs"))
    saveRDS(results_list, paste(sim_results_path, "/Results ", this_run_id, ".rds", sep = ""))
    return(1)
    
}


# -------- #
# (3) Prep #
# -------- #

em_data <- paste("Emulator data/", this_run_code, sep = "")
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
results_codes <- mclapply(1:nsims, function(x) {run_mcmc(x, sim_data_path, sim_results_path, alphabf, burnin, nsamp, ns, nt, start_day, K, Theta, muxnnc, Sxxnnc, muxsbfc, muxtbfc, Sxxsbfc, Sxxtbfc, muxSBF, muxTBF, SxxSBF, SxxTBF, muxM, SxxM, ur_p)}, mc.cores = 10)
saveRDS(results_codes, paste("Results/", this_save_code, "/results_codes.rds", sep = ""))


# ----------- #
# (5) Results #
# ----------- #

beta0_95_cov <- beta0_99_cov <- beta1_95_cov <- beta1_99_cov <- phi_95_cov <- 
    phi_99_cov <- beta0_mse <- beta1_mse <- phi_mse <- 0
not_in_b095 <- not_in_b195 <- not_in_b099 <- not_in_b199 <- not_in_p95 <- not_in_p99 <- c()
true_beta0 <- -2.83
true_beta1 <- 0.1
true_phi <- 0.045

for (i in 1:nsims) {
    
    mcmc <- readRDS(paste("Results/", this_save_code, "/Results ", i, ".rds", sep = ""))
    b0q <- quantile(mcmc$theta[1, (burnin + 1):nruns], c(0.005, 0.025, 0.975, 0.995))
    b1q <- quantile(mcmc$theta[2, (burnin + 1):nruns], c(0.005, 0.025, 0.975, 0.995))
    pq <- quantile(mcmc$theta[3, (burnin + 1):nruns], c(0.005, 0.025, 0.975, 0.995))
    
    beta0_mse <- beta0_mse + (mean(mcmc$theta[1, (burnin + 1):nruns]) - true_beta0)^2 / nsims
    beta1_mse <- beta1_mse + (mean(mcmc$theta[2, (burnin + 1):nruns]) - true_beta1)^2 / nsims
    phi_mse <- phi_mse + (mean(mcmc$theta[3, (burnin + 1):nruns]) - true_phi)^2 / nsims

    if (true_beta0 > b0q[1] && true_beta0 < b0q[4]) {
        beta0_99_cov <- beta0_99_cov + 1
    } else {
        not_in_b099 <- c(not_in_b099, i)
    }
    
    if (true_beta0 > b0q[2] && true_beta0 < b0q[3]) {
        beta0_95_cov <- beta0_95_cov + 1
    } else {
        not_in_b095 <- c(not_in_b095, i)
    }
    
    if (true_beta1 > b1q[1] && true_beta1 < b1q[4]) {
        beta1_99_cov <- beta1_99_cov + 1
    } else {
        not_in_b199 <- c(not_in_b199, i)
    }
    
    if (true_beta1 > b1q[2] && true_beta1 < b1q[3]) {
        beta1_95_cov <- beta1_95_cov + 1
    } else {
        not_in_b195 <- c(not_in_b195, i)
    }
    
    if (true_phi > pq[1] && true_phi < pq[4]) {
        phi_99_cov <- phi_99_cov + 1
    } else {
        not_in_p99 <- c(not_in_p99, i)
    }
    
    if (true_phi > pq[2] && true_phi < pq[3]) {
        phi_95_cov <- phi_95_cov + 1
    } else {
        not_in_p95 <- c(not_in_p95, i)
    }
    
}

conc_list <- list(b095 = beta0_95_cov / nsims, b195 = beta1_95_cov / nsims, p95 = phi_95_cov / nsims,
                  b099 = beta0_99_cov / nsims, b199 = beta1_99_cov / nsims, p99 = phi_99_cov / nsims,
                  not_in_b095 = not_in_b095, not_in_b195 = not_in_b195, not_in_p95 = not_in_p95,
                  not_in_b099 = not_in_b099, not_in_b199 = not_in_b199, not_in_p99 = not_in_p99,
                  beta0_mse = beta0_mse, beta1_mse = beta1_mse, phi_mse = phi_mse)
saveRDS(conc_list, paste("Results/", this_save_code, "/Final results.rds", sep = ""))

# eof