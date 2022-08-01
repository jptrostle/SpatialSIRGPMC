
my_gillespie <- function(this_iter, path, ns, nt_override_value, beta, phi, gamma, under_p = 1, od = 2, ARG_GP, ARG_MUX, ARG_MUY, max_obs = 5000000) {

    RESET_MAX_OBS <- max_obs
    VALID <- FALSE
    while (!VALID) {

        # legacy code
        max_obs <- RESET_MAX_OBS 

        # read starting values defined in the master script
        N <- ARG_GP
        X <- ARG_MUX
        Y <- ARG_MUY

        # save the simulated X and Y as well at key points
        SIM_MUX <- matrix(0, nrow=1, ncol=(ns+1))
        SIM_MUX[, 2:(ns+1)] <- X
        SIM_MUY <- matrix(0, nrow=1, ncol=(ns+1))
        SIM_MUY[, 2:(ns+1)] <- Y

        # expand matrices and make them huge...
        Y <- cbind(Y, matrix(0, nrow=ns, ncol=1))
        X <- cbind(X, matrix(0, nrow=ns, ncol=1))

        # list the adjacencies
        td <- sqrt(ns)
        neighbors <- list("1"=c(2, 1+td))
        for (i in 2:ns) {
            these_neigh <- c()
            if (i > td) these_neigh <- c(these_neigh, i - td)
            if (i <= (ns-td)) these_neigh <- c(these_neigh, i+td)
            if (i %% td != 1) these_neigh <- c(these_neigh, i-1)
            if (i %% td != 0) these_neigh <- c(these_neigh, i+1)
            neighbors[[as.character(i)]] <- these_neigh
        }

        # simulate dataset
        time <- 1
        i <- 2
        all_new_inf <- matrix(0, nrow=max_obs, ncol=3) # the first row is a dummy row to be deleted later
        num_new_infs <- 1 # tracker for which row
        t1 <- Sys.time()
        time_tracker <- rep(0, max_obs)
        last_changed_site <- 0 # to avoid recalculating all the rates over and over and over again
        local_inf_count <- 0
        local_rec_count <- 0
        spatial_count <- 0

        KEEP_GOING <- TRUE
        while (KEEP_GOING) {
            
            prev <- i - 1
            
            # calculate all local infection rates
            if (i == 2) local_inf_events <- rep(0, ns)
            for (j in 1:ns) {
                if (i == 2 || j == last_changed_site) {
                    local_inf_events[j] <- beta / N[j] * X[j,1] * Y[j,1]
                }
            }
            
            # calculate all local recovery rates
            if (i == 2) local_rec_events <- rep(0, ns)
            for (j in 1:ns) {
                if (i == 2 || j == last_changed_site) {
                    local_rec_events[j] <- gamma * Y[j,1]
                }
            }
            
            # calculate all spatial infection rates
            # how to calculate number of spatial events: the four corners have two neighbors
            # the remaining perimeter sites (equals (sqrt(ns)-2)*4) have three neighbors
            # the remaining interior sites have four neighbors
            num_spatial_events <- 4*2 + (sqrt(ns)-2)*4*3 + ((sqrt(ns)-2)^2)*4
            if (i == 2) spatial_events <- rep(0, num_spatial_events)
            if (i == 2) event_list <- list()
            event_tracker <- 1
            for (j in 1:ns) {
                for (k in neighbors[[j]]) {
                    if (i == 2 || (k == last_changed_site || j == last_changed_site)) {
                        # this is the event of site j infecting site k
                        spatial_events[event_tracker] <- phi / N[k] * X[k,1] * Y[j,1]
                        if (i == 2) event_list[[event_tracker]] <- c(j,k)
                    }
                    event_tracker <- event_tracker + 1
                }
            }
            
            # define all_combined_rates
            all_combined_rates <- sum(local_inf_events) + sum(local_rec_events) + sum(spatial_events)
            
            # do the random sampling now
            u1 <- runif(1)
            u2 <- runif(1)
            time_till_event <- -log(u1) / all_combined_rates
            time <- time + time_till_event
            time_tracker[i] <- time
            
            # decide what event happens
            all_events <- c(local_inf_events, local_rec_events, spatial_events)
            if (sum(all_events) == 0) {
                max_obs <- i
                break
            }
            num_events <- length(local_inf_events) + length(local_rec_events) + length(spatial_events)
            which_event <- sample(1:num_events, 1, prob=all_events)
            X[,2] <- X[,1]
            Y[,2] <- Y[,1]
            if (which_event <= length(local_inf_events)) {
                X[which_event,2] <- X[which_event,2] - 1
                Y[which_event,2] <- Y[which_event,2] + 1
                num_new_infs <- num_new_infs + 1
                all_new_inf[num_new_infs,] <- c(time, 1, which_event)
                last_changed_site <- which_event
                local_inf_count <- local_inf_count + 1
            } else if (which_event > length(local_inf_events) && which_event <= (length(local_inf_events) + length(local_rec_events))) {
                j <- which_event - length(local_inf_events)
                Y[j,2] <- Y[j,2] - 1
                last_changed_site <- j
                local_rec_count <- local_rec_count + 1
            } else {
                this_one <- which_event - length(local_inf_events) - length(local_rec_events)
                k <- event_list[[this_one]][2]
                X[k,2] <- X[k,2] - 1
                Y[k,2] <- Y[k,2] + 1
                num_new_infs <- num_new_infs + 1
                all_new_inf[num_new_infs,] <- c(time, 1, k)
                last_changed_site <- k
                spatial_count <- spatial_count + 1
            }

            # if it's a new day, save out the susceptibles and infectious from the end of the previous day
            if (i > 2 && floor(time_tracker[prev]) != floor(time_tracker[i]) && floor(time_tracker[prev]) != 0) {
                save_day <- floor(time_tracker[prev])
                SIM_MUX <- rbind(SIM_MUX, c(save_day, as.vector(X[,1])))
                SIM_MUY <- rbind(SIM_MUY, c(save_day, as.vector(Y[,1])))
            }

            # don't forget to increment!
            i <- i + 1

            # memory-efficiency solution
            X[,1] <- X[,2]
            Y[,1] <- Y[,2]

            # check if we're done running
            if (time > (nt_override_value + 1)) KEEP_GOING <- FALSE

        }
        t2 <- Sys.time()

        # delete the dummy row in the new-infections matrix
        stopifnot(sum(all_new_inf[(num_new_infs+1):(max_obs-1),]) == 0)
        all_new_inf <- all_new_inf[2:num_new_infs,]

        # V2 to V3 change: binom used to be here, but now it's negative binom

        # convert rather tediously the new infections into a matrix of daily new counts
        nt <- ceiling(max(all_new_inf[,1]))
        simdata <- matrix(0, nrow=ns, ncol=nt)
        for (i in 1:ns) {
            these_inf <- all_new_inf[all_new_inf[,3] == i,]
            asdf <- as.data.frame(these_inf)
            asdf$V1 <- floor(asdf$V1)
            aggdata <- by(asdf$V2, asdf$V1, sum)
            for (j in 1:nt) {
                if (j %in% names(aggdata)) {
                    k <- which(names(aggdata) == j)
                    simdata[i,j] <- aggdata[k]
                }
            }
        }

        # the last day is an incomplete day and can mess up parameter inference
        nt <- nt_override_value
        if (dim(simdata)[2] < nt) {
            print("Failed because of not enough time -- probably all recoveries")
            next
        }
        simdata <- simdata[, 1:nt]

        # nbinom step
        # warnings are suppressed so that simdata = 0 doesn't cause an error; I have options(warn = 2) set in the master script
        suppressWarnings(simdata <- matrix(rnbinom(ns * nt, mu = under_p * simdata, size = under_p * simdata / (od - 1)), nrow = ns, ncol = nt, byrow = FALSE))

        rm(event_list, these_inf, all_new_inf, asdf, aggdata, all_combined_rates, all_events, event_tracker, i, j, k, last_changed_site, local_inf_count, local_inf_events, local_rec_count, local_rec_events, num_events, num_new_infs, num_spatial_events, t1, t2, these_neigh, this_one, spatial_count, spatial_events, time_till_event, time, u1, u2, which_event, td, prev)
        rm(X, Y, time_tracker, neighbors)

        # save out just simdata
        rownames(simdata) <- NULL        
        saveRDS(simdata, paste(path, "/", this_iter, ".rds", sep=""))
        VALID <- TRUE

        # also save out simulated MUX and MUY at different time points
        saveRDS(SIM_MUX, paste(path, "/sim_mux_", this_iter, ".rds", sep=""))
        saveRDS(SIM_MUY, paste(path, "/sim_muy_", this_iter, ".rds", sep=""))

    }

    return(1)

}

# eof