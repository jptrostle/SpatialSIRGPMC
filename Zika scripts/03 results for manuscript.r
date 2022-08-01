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
setwd("..")
load("zika results.rdata")

sseq <- 100001:350000

# point estimates
thetahat <- rowMeans(theta[, sseq])
nuhat <- mean(nu[sseq])
alphahat <- apply(alpha[,, sseq], c(1, 2), mean)

# standard errors
thetase <- apply(theta[, sseq], 1, sd)
nuse <- sd(nu[sseq])
alphase <- apply(alpha[,, sseq], c(1, 2), sd)

# 95% credible intervals
theta025 <- apply(theta[, sseq], 1, function(x) quantile(x, 0.025))
theta975 <- apply(theta[, sseq], 1, function(x) quantile(x, 0.975))

# get lambdahat
muxhat <- run_mux_em(c(thetahat, 10), Theta, muxnnc, muxsbfc, muxtbfc, muxM, muxSBF, muxTBF, Theta_m, Theta_s)
Sbfhat <- run_Sxx_em(c(thetahat, 10), Theta, Sxxnnc, Sxxsbfc, Sxxtbfc, SxxM, SxxSBF, SxxTBF, lnt, Theta_m, Theta_s)
lambdahat <- update_lambda(muxhat, Sbfhat, Ba_to_alpha(B, alphahat), ns, lnt, ur_p[, 1])

# plot Rio de Janeiro
par(cex = 1.25)
plot(y[19, ], ylim = c(0, 6000), xlab = "Weeks since 40th week of 2015", ylab = "Rio de Janeiro new infections")
lines(lambdahat[19, ], lw = 2)

# plot Bahia
plot(y[5, ], xlab = "Weeks since 40th week of 2015", ylab = "Bahia new infections")
lines(lambdahat[5, ], lw = 2)

# make easy results table
results <- data.frame(param = c("beta0", "beta1", "phi"), point = thetahat, se = thetase, lb95 = theta025, ub95 = theta975)

# trace plots
par(cex = 1.25)
plot(theta[1, sseq], type = "l", xlab = "", ylab = "beta0 trace plot")
plot(theta[2, sseq], type = "l", xlab = "", ylab = "beta1 trace plot")
plot(theta[3, sseq], type = "l", xlab = "", ylab = "phi trace plot")


##########################

# measuring model discrepancy
md <- rep(0, 27)
for (x in 1:27) md[x] <- norm(y[x, ] - lambdahat[x, ], "2") / norm(y[x, ], "2")
md <- 1 - md

library(units)
library(geobr)
library(sf)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

states <- read_country(year = 2019)
states$name_state <- tolower(states$name_state)
dados <- structure(list(X = 1:27, uf = states$name_state[order(states$name_state)], md = md), class = "data.frame", row.names = c(NA, -27L))
states <- dplyr::left_join(states, dados, by = c("name_state" = "uf"))

p <- states %>% ggplot() + 
    geom_sf(aes(fill = md), size = .15) + 
    scale_fill_gradient(low = "red", high = "blue", limits = c(0.17, 0.84), name = "1 - MD") + 
    xlab("") + 
    ylab("") +
    theme_bw()


p <- p + labs(title = "") +
    theme(plot.caption = element_text(hjust = 0, face = "italic"), 
          plot.title.position = "plot", 
          plot.caption.position =  "plot")

plot(p)


# plot Piaui
plot(y[18, ], xlab = "Weeks since 40th week of 2015", ylab = "Piaui new infections")
lines(lambdahat[18, ], lw = 2)





# eof