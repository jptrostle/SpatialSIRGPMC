# Zika plot, scaled by population size

rm(list = ls())
load("../zika Brazil.rdata")

zika <- cbind(zika_data[,, 1], zika_data[,, 2])
zika <- zika[, 41:80]

zika_norm <- zika
for (s in 1:27) zika_norm[s, ] <- zika_norm[s, ] / gridpop[s]
zika_norm <- zika_norm * 100

par(cex = 1.25)
setEPS()
postscript("brazil.eps")
plot(zika_norm[1,], type = "l", ylim = c(0, 0.085), xlab = "Weeks since 40th week of 2015", ylab = "Weekly new infections (% of population)")
for (s in 2:27) lines(zika_norm[s, ])
dev.off()

