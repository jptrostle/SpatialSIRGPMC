# Zika model
# For the sensitivites in the web appendix, change muxsbfc and/or Sxxsbfc as necessary

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
library(splines)
library(MASS)
set.seed(19190979)
setwd("..")

# --------- #
# (1) Setup #
# --------- #

load("Data/zika Brazil.rdata")
y <- cbind(zika_data[,, 1], zika_data[,, 2])[, 41:80]
rm(Adj, area_km, states_names, zika_data)

ns <- 27
start_day <- 41
nt <- 40
n <- ns * nt
lnt <- nt + 1
nx <- 2

burnin <- 100000
nsamp <- 250000
nruns <- burnin + nsamp

muxnnc <- 20
Sxxnnc <- 20
muxsbfc <- 20
muxtbfc <- 10
Sxxsbfc <- 10
Sxxtbfc <- 10
alphabf <- 5

Theta <- readRDS("Emulator files/ClTheta.rds")
muxM <- readRDS("Emulator files/muxM.rds")
muxSBF <- readRDS("Emulator files/muxSBF.rds")
muxTBF <- readRDS("Emulator files/muxTBF.rds")
SxxM <- readRDS("Emulator files/SxxM.rds")
SxxSBF <- readRDS("Emulator files/SxxSBF.rds")
SxxTBF <- readRDS("Emulator files/SxxTBF.rds")

# a fix for Theta
Theta <- cbind(Theta[, 1:3], 10)


####################################################################

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

# emulator function for the covariances
run_Sxx_em <- function(theta0, Theta, nnc, nsbf, ntbf, M, U1, U3, lnt, tm, ts) {
    
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

update_lambda <- function(mux, Sbf, alpha, ns, lnt, ur_p) {
    
    X <- matrix(0, nrow = ns, ncol = lnt)
    mux <- matrix(mux, nrow = ns, ncol = lnt)
    
    for (t in 1:lnt) X[, t] <- mux[, t] + Sbf[[t]] %*% alpha[, t]
    
    lambda <- ur_p * (X[, 1:(lnt - 1)] - X[, 2:lnt])
    lambda[lambda < 0.00001] <- 0.00001
    
    return(lambda)
    
}


####################################################################

# -------------- #
# (2) MCMC setup #
# -------------- #

Theta_m <- colMeans(Theta[, 1:(nx + 1)])
Theta_s <- apply(Theta[, 1:(nx + 1)], 2, sd)
for (i in 1:(nx + 1)) Theta[, i] <- (Theta[, i] - Theta_m[i]) / Theta_s[i]

# splines setup
lnt <- nt + 1
B <- bs(1:lnt, df = alphabf, intercept = TRUE)
num_std <- diag(diag(B %*% t(B))^-0.5)
B <- num_std %*% B
Ba_to_alpha <- function(B, alpha) return(t(B %*% t(alpha)))

# buckets for results
theta <- matrix(0.05, nrow = nx + 1, ncol = nruns)
nu <- rep(0, nruns)
s <- rep(10, nruns)
alpha <- array(0, dim = c(Sxxsbfc, alphabf, nruns))

nu[1] <- 2
theta[, 1] <- c(-0.06, 0.18, 0.085)

ur_p <- matrix(0, nrow = ns, ncol = nruns)
# from preliminary analysis
ur_p[, 1] <- c(0.01968, 0.00560, 0.00721, 0.02370, 0.01917, 0.00480, 0.00055, 0.00180, 0.00548, 0.00848, 0.04459, 0.00260, 0.00394, 0.01293, 0.00361, 0.00132, 0.00062, 0.00044, 0.01386, 0.00548, 0.00201, 0.00886, 0.00219, 0.00028, 0.00108, 0.00100, 0.00801)
ur_p_accepts <- rep(0, ns)
ur_p_sd <- rep(0.001, ns)

# proposals and acceptance rates
theta_accepts1 <- theta_accepts2 <- nu_accepts <- 0
alpha_accepts <- matrix(0, nrow = Sxxsbfc, ncol = alphabf)

BP_V <- matrix(c(3.066944e-06, -7.681685e-07, -5.944987e-07, 1.287156e-05, 3.104339e-06, -4.158287e-05, 2.429473e-05, 1.026682e-04,
-7.681685e-07, 5.123048e-07, 2.989016e-08, -1.563372e-06, 3.020410e-06, 1.924299e-05, -4.054162e-05, -6.858674e-05,
-5.944987e-07, 2.989016e-08, 2.492362e-07, -1.818127e-05, 7.971039e-06, -3.197737e-05, 5.186556e-05, 3.040627e-05,
1.287156e-05, -1.563372e-06, -1.818127e-05, 9.657536e-03, -7.717143e-03, 1.917206e-02, -3.206305e-02, -3.594041e-03,
3.104339e-06, 3.020410e-06, 7.971039e-06, -7.717143e-03, 6.884789e-03, -1.659099e-02, 2.742003e-02, 3.973860e-03,
-4.158287e-05, 1.924299e-05, -3.197737e-05, 1.917206e-02, -1.659099e-02, 5.182120e-02, -8.779930e-02, -1.451517e-02,
2.429473e-05, -4.054162e-05, 5.186556e-05, -3.206305e-02, 2.742003e-02, -8.779930e-02, 1.894963e-01, 5.679677e-03,
1.026682e-04, -6.858674e-05, 3.040627e-05, -3.594041e-03, 3.973860e-03, -1.451517e-02, 5.679677e-03, 4.015668e-01), nrow = 8, byrow = TRUE)
BP_V <- (BP_V + t(BP_V)) / 2
tC_BP_V <- t(chol(BP_V))

nu_sd <- 1
alpha_sds <- matrix(1, nrow = Sxxsbfc, ncol = alphabf)
alpha_sds[1, ] <- -3 # this is a placeholder to make sure no problems; it shouldn't ever be used


# ------------- #
# (6) Main MCMC #
# ------------- #

# from preliminary work
alpha[,, 1] <- matrix(c(-0.58868523, -0.6755528, -1.25215420, 0.57448128, -0.50665968,
0.21612563, 1.4247358, -2.83314426, -0.93406372, -0.46252321,
0.30350830, 1.2330043, -0.63719421, -1.54401076, 0.08035355,
0.53038338, 0.7887758, 1.01426561, 0.44656530, 2.33431874,
0.22699163, -0.4147780, -1.02176256, -0.16114805, -1.96875059,
0.09320285, -0.5046840, -0.92270969, -0.24968942, 3.09561268,
-0.88590299, 0.7317008, 1.15946645, 0.29265601, -0.68884741,
0.28287296, 1.3320015, 1.36774416, 0.62915089, -0.48508681,
0.65713037, -0.1392721, -0.18573628, -0.04004736, 0.18996825,
0.57529212, 0.6454127, 0.06808514, -0.12520894, 0.06082921), nrow = Sxxsbfc, byrow = TRUE)

mux <- run_mux_em(c(theta[, 1], s[1]), Theta, muxnnc, muxsbfc, muxtbfc, muxM, muxSBF, muxTBF, Theta_m, Theta_s)
Sbf <- run_Sxx_em(c(theta[, 1], s[1]), Theta, Sxxnnc, Sxxsbfc, Sxxtbfc, SxxM, SxxSBF, SxxTBF, lnt, Theta_m, Theta_s)
lambda <- update_lambda(mux, Sbf, Ba_to_alpha(B, alpha[,, 1]), ns, lnt, ur_p[, 1])

t1 <- Sys.time()
for (i in 2:nruns) {
    
    if (i %% 1000 == 0) print(paste("On iter", i, "at", Sys.time()))
    prev <- i - 1
    
    # ----------------- #
    # beta, phi, alpha1 #
    # ----------------- #
    
    first_stage_accepted <- FALSE
    
    cand_vec <- as.vector(tC_BP_V %*% rnorm(alphabf + nx + 1) + c(theta[, prev], alpha[1,, prev]))
    cand_theta <- c(cand_vec[1:(nx + 1)], s[prev])
    alphamaybe <- alpha[,, prev]
    alphamaybe[1, ] <- cand_vec[(nx + 2):(nx + alphabf + 1)]
    muxmaybe <- run_mux_em(cand_theta, Theta, muxnnc, muxsbfc, muxtbfc, muxM, muxSBF, muxTBF, Theta_m, Theta_s)
    Sbfmaybe <- run_Sxx_em(cand_theta, Theta, Sxxnnc, Sxxsbfc, Sxxtbfc, SxxM, SxxSBF, SxxTBF, lnt, Theta_m, Theta_s)
    lmaybe <- update_lambda(muxmaybe, Sbfmaybe,  Ba_to_alpha(B, alphamaybe), ns, lnt, ur_p[, prev])
    
    log_orig <- sum(dnbinom(y, mu = lambda, size = lambda / (nu[prev] - 1), log = TRUE)) + dmvnorm(theta[, prev], rep(0, nx + 1), diag(nx + 1) * 100, log = TRUE) + sum(dnorm(alpha[1, , prev], mean = 0, sd = 1, log = TRUE))
    log_cand <- sum(dnbinom(y, mu = lmaybe, size = lmaybe / (nu[prev] - 1), log = TRUE)) + dmvnorm(cand_theta[1:(nx + 1)], rep(0, nx + 1), diag(nx + 1) * 100, log = TRUE) + sum(dnorm(alphamaybe[1, ], mean = 0, sd = 1, log = TRUE))
    
    log_accep_prob <- log_cand - log_orig
    if (log(runif(1)) < log_accep_prob) {
        theta_accepts1 <- theta_accepts1 + 1
        theta[, i] <- cand_theta[1:(nx + 1)]
        lambda <- lmaybe
        mux <- muxmaybe
        Sbf <- Sbfmaybe
        alpha[1,, i] <- alphamaybe[1, ]
        alpha_accepts[1, ] <- alpha_accepts[1, ] + 1
        first_stage_accepted <- TRUE
    } 
    
    if (!first_stage_accepted) {
        
        cand_vec2 <- as.vector(0.2 * tC_BP_V %*% rnorm(nx + alphabf + 1) + c(theta[, prev], alpha[1, , prev]))
        cand_theta2 <- c(cand_vec2[1:(nx + 1)], s[prev])
        alphamaybe2 <- alpha[,, prev]
        alphamaybe2[1, ] <- cand_vec2[(nx + 2):(nx + alphabf + 1)]
        muxmaybe2 <- run_mux_em(cand_theta2, Theta, muxnnc, muxsbfc, muxtbfc, muxM, muxSBF, muxTBF, Theta_m, Theta_s)
        Sbfmaybe2 <- run_Sxx_em(cand_theta2, Theta, Sxxnnc, Sxxsbfc, Sxxtbfc, SxxM, SxxSBF, SxxTBF, lnt, Theta_m, Theta_s)
        lmaybe2 <- update_lambda(muxmaybe2, Sbfmaybe2, Ba_to_alpha(B, alphamaybe2), ns, lnt, ur_p[, prev])
        
        log_cand2 <- sum(dnbinom(y, mu = lmaybe2, size = lmaybe2 / (nu[prev] - 1), log = TRUE)) + dmvnorm(cand_theta2[1:(nx + 1)], rep(0, nx + 1), diag(nx + 1) * 100, log = TRUE) + sum(dnorm(alphamaybe2[1, ], mean = 0, sd = 1, log = TRUE))
        
        top_of_accep2 <- log_cand2 + dmvnorm(cand_vec, mean = cand_vec2, sigma = BP_V, log = TRUE) + (1 - exp(log_cand - log_cand2))
        bottom_of_accep2 <- log_orig + dmvnorm(cand_vec, mean = c(theta[, prev], alpha[1, , prev]), sigma = BP_V, log = TRUE) + (1 - exp(log_accep_prob))
        log_accep_prob2 <- top_of_accep2 - bottom_of_accep2
        if (log(runif(1)) < log_accep_prob2) {
            
            theta_accepts2 <- theta_accepts2 + 1
            theta[, i] <- cand_theta2[1:(nx + 1)]
            lambda <- lmaybe2
            mux <- muxmaybe2
            Sbf <- Sbfmaybe2
            alpha[1,, i] <- alphamaybe2[1, ]
            alpha_accepts[1, ] <- alpha_accepts[1, ] + 1
            
        } else {
            
            theta[, i] <- theta[, prev]
            alpha[1,, i] <- alpha[1,, prev]
            
        }
        
    }
    
    # AM step
    if (i %% 200 == 0 && i <= burnin && alpha_accepts[1] > 10) {
        
        adj_index <- 1:i
        
        BP_V <- 1.5 * cov(t(rbind(theta[, adj_index], alpha[1,, adj_index]))) + diag(nx + alphabf + 1) * 1e-8
        tC_BP_V <- t(chol(BP_V))
        
    }        
    
    # ----- #
    # s, nu #
    # ----- #
    
    # assuming S0 corresponds to Maranhao, leave it as S0 = 10
    s[i] <- 10
    
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
            alphamaybe <- alpha[,, prev]
            alphamaybe[a, b] <- rnorm(1, mean = alpha[a, b, prev], sd = alpha_sds[a, b])
            if (a > 1) alphamaybe[1:(a - 1), ] <- alpha[1:(a - 1),, i]
            if (b > 1) alphamaybe[a, 1:(b - 1)] <- alpha[a, 1:(b - 1), i]
            lmaybe <- update_lambda(mux, Sbf, Ba_to_alpha(B, alphamaybe), ns, lnt, ur_p[, prev])
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
    
    ur_p[, i] <- ur_p[, prev]
    
    # ------------------------------- #
    # bookkeeping, tkd sd adjustments #
    # ------------------------------- #
    
    every_num <- 400
    if (i <= burnin && i %% every_num == 0 && i >= every_num) {
        
        # update nu sds
        if (nu_accepts / every_num < 0.2) nu_sd <- nu_sd * 0.9
        if (nu_accepts / every_num > 0.25) nu_sd <- nu_sd * 1.1
        
        for (a in 1:Sxxsbfc) {
            for (b in 1:alphabf) {
                if (a == 1 && b == 1) next
                if (alpha_accepts[a, b] / every_num < 0.2) alpha_sds[a, b] <- alpha_sds[a, b] * 0.9
                if (alpha_accepts[a, b] / every_num > 0.25) alpha_sds[a, b] <- alpha_sds[a, b] * 1.1   
            }
        }
        
        for (x in 1:ns) {
            if (ur_p_accepts[x] / 400 < 0.2) ur_p_sd[x] <- ur_p_sd[x] * 0.9
            if (ur_p_accepts[x] / 400 > 0.3) ur_p_sd[x] <- ur_p_sd[x] * 1.1
        }
        
        nu_accepts <- 0
        alpha_accepts[2:Sxxsbfc, ] <- matrix(0, nrow = Sxxsbfc - 1, ncol = alphabf)
        ur_p_accepts <- rep(0, ns)
    }
    
    if (i == burnin) {
        theta_accepts1 <- theta_accepts2 <- nu_accepts <- 0
        alpha_accepts <- matrix(0, nrow = Sxxsbfc, ncol = alphabf)
        ur_p_accepts <- rep(0, ns)
    }
    
    if (i %% 1000 == 0) {
        par(mfrow = c(2, 2))
        plot(theta[1, 1:i], xlab = "", ylab = "b0 chain", main = "Beta0")
        abline(h = -0.2, col = "red")
        abline(h = 0.05, col = "red")
        plot(theta[2, 1:i], xlab = "", ylab = "b1 chain", main = "Beta1")
        abline(h = 0.13, col = "red")
        abline(h = 0.27, col = "red")
        plot(theta[3, 1:i], xlab = "", ylab = "phi chain", main = "Phi")
        abline(h = 0.055, col = "red")
        abline(h = 0.135, col = "red")
        plot(alpha[1,1,1:i], xlab = "", ylab = "a11 chain", main = "Alpha(1, 1)")
    }
    
}
t2 <- Sys.time()

# eof