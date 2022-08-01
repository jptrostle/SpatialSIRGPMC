# Run on a university cluster using Rmpi

# Change setwd as necessary
setwd("..")
library(Rmpi)
library(snow)
library(rTensor)
options(warn = 2)
set.seed(58207423)
NCORES <- 50

########################################

# ------------------ #
# Main code settings #
# ------------------ #

path <- "sim2"
nsims <- 200

# DGP settings
generate_data <- TRUE
ns <- 25
start_day <- 21
nt <- 80
logdens <- log(runif(ns, 0, 3000))
true_beta0 <- -2.83
true_beta1 <- 0.1
true_beta <- exp(true_beta0 + true_beta1 * logdens)
true_phi <- 0.045
true_gamma <- 0.040
dgp_od <- 2.5
start_inf <- 100
S0 <- 13
underreport <- 1.00

# Emulators: experimental settings
generate_emulator_data <- TRUE
generate_em_bf <- TRUE
K <- 60000
muxnnc <- 10 # nearest neighbors
Sxxnnc <- 10
muxsbfc <- 20 # spatial / temporal basis-function counts
muxtbfc <- 10
Sxxsbfc <- 10
Sxxtbfc <- 10
alphabf <- 5 # including constant

# --------- #
# (0) Setup #
# --------- #

cl <- makeCluster((mpi.universe.size() - 1), type = 'MPI')

# Paths: automatically generate here
sim_data_path <- paste("Data/", path, sep = "")
sim_em_data_path <- paste("Emulator data/", path, sep = "")
sim_results_path <- paste("Results/", path, sep = "")
suppressWarnings(dir.create(sim_data_path))
suppressWarnings(dir.create(sim_em_data_path))
suppressWarnings(dir.create(sim_results_path))

# sanity checks
if (K %% ns != 0) stop("K and ns are not easily compatible")
if (sqrt(ns) != as.integer(sqrt(ns))) stop("Code only works with square spatial domains")

saveRDS(logdens, paste(sim_data_path, "/logdens.rds", sep = ""))

##############################################

# ------------------------------------------ #
# (1) Simulate spatial-SIR jump-process data #
# ------------------------------------------ #

MUX <- rep(100000, ns)
GRID_POP <- MUX
MUY <- rep(0, ns)
MUY[S0] <- start_inf
MUX <- MUX - MUY

if (generate_data) {    

    source("Functions/custom gillespie v4 - logdens.r", echo = FALSE)    
    
    clusterExport(cl, "sim_data_path")
    clusterExport(cl, "ns")
    clusterExport(cl, "start_day")
    clusterExport(cl, "nt")
    clusterExport(cl, "true_beta")
    clusterExport(cl, "true_phi")
    clusterExport(cl, "true_gamma")
    clusterExport(cl, "underreport")
    clusterExport(cl, "dgp_od")
    clusterExport(cl, "my_gillespie")
    clusterExport(cl, "MUX")
    clusterExport(cl, "MUY")
    clusterExport(cl, "GRID_POP")
    
    gill_run <- function(x) my_gillespie(x, sim_data_path, ns, start_day + nt - 1, true_beta, true_phi, true_gamma, under_p = underreport, od = dgp_od, GRID_POP, MUX, MUY, max_obs = 30000000)
    dgp.lines <- clusterApply(cl = cl, 1:nsims, fun = gill_run)
    
}

##############################################

# -------------------------------------- #
# (2) Run experimental forward equations #
# -------------------------------------- #

if (generate_emulator_data) {
    
    source("Functions/full forward equations v3.r")
    
    # adjacency matrix
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
    
    # setup for my forward-equations code
    coords <- make_forward_calc_sequences(ns, Adj)
    matrix_of_zeroes <- matrix(0, nrow = ns, ncol = ns)
    
    # make Theta
    # sim3: a few changes here to account for b0 and b1 + logdens
    beta0_vector <- seq(-3.5375, -2.1225, length.out = K)
    beta1_vector <- sample(seq(-0.05, 0.25, length.out = K), K, replace = FALSE)
    phi_vector <- sample(seq(0.03375, 0.05625, length.out = K), K, replace = FALSE)
    start_locs <- sample(rep(c(8, 12, 13, 14, 18), each = K / 5), K, replace = FALSE)
    Theta <- cbind(beta0_vector, beta1_vector, phi_vector, start_locs)
    
    # inline function
    gen_fe <- function(i, Theta, ns, grid_pop, matrix_of_zeroes, coords, start_t, nt, gamma, Adj, sim_em_data_path, start_inf, logdens) {
        
        this_beta <- exp(Theta[i, 1] + Theta[i, 2] * logdens)
        this_phi <- Theta[i, 3]
        muy0 <- rep(0, ns)
        muy0[Theta[i, 4]] <- start_inf
        mux0 <- grid_pop - muy0
        
        this_sim <- update_XYparam_VECT_full(coords, ns, start_t + nt - 1, mux0, muy0, matrix_of_zeroes, matrix_of_zeroes, matrix_of_zeroes, this_beta, this_phi, gamma, Adj, grid_pop, num_breaks = 10)
        
        saveRDS(as.vector(this_sim$mux[, start_t:(start_t + nt)]), paste(sim_em_data_path, "/mux", i, ".rds", sep = ""))
        saveRDS(this_sim$sigmaxx[, start_t:(start_t + nt)], paste(sim_em_data_path, "/Sxx", i, ".rds", sep = ""))
        
    }
    
	clusterExport(cl, "Theta")
	clusterExport(cl, "ns")
	clusterExport(cl, "GRID_POP")
	clusterExport(cl, "matrix_of_zeroes")
	clusterExport(cl, "coords")
	clusterExport(cl, "start_day")
	clusterExport(cl, "nt")
	clusterExport(cl, "true_gamma")
	clusterExport(cl, "Adj")
	clusterExport(cl, "sim_em_data_path")
	clusterExport(cl, "update_XYparam_VECT_full")
    clusterExport(cl, "gen_fe")
    clusterExport(cl, "start_inf")
    clusterExport(cl, "logdens")
	
    emex_run <- function(x) gen_fe(x, Theta, ns, GRID_POP, matrix_of_zeroes, coords, start_day, nt, true_gamma, Adj, sim_em_data_path, start_inf, logdens)
    emex.lines <- clusterApply(cl = cl, 1:K, fun = emex_run)
    saveRDS(Theta, paste(sim_em_data_path, "Theta_setup.rds"))
    
    CSMatrix_mux <- matrix(0, nrow = ns * (nt + 1), ncol = K)
    CSArray_sxx <- array(0, dim = c(ns * (ns + 1) / 2, nt + 1, K))
    for (k in 1:K) {
        CSMatrix_mux[, k] <- readRDS(paste(sim_em_data_path, "/mux", k, ".rds", sep = ""))
        CSArray_sxx[,, k] <- readRDS(paste(sim_em_data_path, "/Sxx", k, ".rds", sep = ""))
    }
    
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

    for (k in 1:K) {
        file.remove(paste(sim_em_data_path, "/mux", k, ".rds", sep = ""))
        file.remove(paste(sim_em_data_path, "/Sxx", k, ".rds", sep = ""))
    }
    
}

##############################################

# ---------------------- #
# (3) Last emulator step #
# ---------------------- #

source("Functions/mcmc functions -- emulator.r")
generate_em_bf_weights(sim_em_data_path, Theta, ns, nt, K)

stopCluster(cl)
mpi.exit()

# eof