require(deSolve) #loading the ODE solver
require(PBSddesolve)

# Bring in calculated optimal results #
# For switches for both distributions #
switches <- read.csv(file = "~/results/switch_results.csv")

# For polynomials
# "2" represents medium variability, "4" represents very high variability
four.2 <- read.csv(file= "~/results/medium_variability_polynomial_results.csv")
four.4 <- read.csv(file = "~/results/very_high_variability_polynomial_results.csv")

# Separate out by their four parameters of their strategy
four.kappa1 <- c(four.2[,3],four.4[,3])
four.kappa2 <- c(four.2[,4],four.4[,4])
four.kappa3 <- c(four.2[,5],four.4[,5])
four.kappa4 <- c(four.2[,6],four.4[,6])

four.fitness <- c(-four.2[,7],-four.4[,7])

### Make vectors for the strategies so that they can be plotted ###
time.med <- seq(1,628,1)
time.very.high <- seq(1,1152,1)

time.med.10 <- seq(0.1,628,0.1)
time.very.high.10 <- seq(0.1,1152,0.1)


# Just looking at hazards 26 and 29
dist2.four <- c()
dist4.four <- c()


for(i in 1:length(time.med)){
  dist2.four[i] <- ((exp(-exp((four.kappa1[1]*i^3)+(four.kappa2[1]*i^2)+(four.kappa3[1]*i)+four.kappa4[1])))+0.000001)
}

for(i in 1:length(time.very.high)){
  dist4.four[i] <- ((exp(-exp((four.kappa1[2]*i^3)+(four.kappa2[2]*i^2)+(four.kappa3[2]*i)+four.kappa4[2])))+0.000001)
}


dist2.switch <- c()
for(i in 1:length(time.med.10)){
  current.t <- i
  ifelse(current.t>switches$Switch[1],dist2.switch[i]<-switches$Kappa1[1],dist2.switch[i]<-switches$Kappa2[1])
}


dist4.switch <- c()
for(i in 1:length(time.very.high.10)){
  current.t <- i
  ifelse(current.t>switches$Switch[2],dist4.switch[i]<-switches$Kappa1[2],dist4.switch[i]<-switches$Kappa2[2])
}


### Need to compare to fixed results#
results <- read.csv(file="~/results/fixed.results.virulence.00001.csv")

dist2.one <- rep(results[26,3],length(time.med.10))
dist4.one <- rep(results[29,3],length(time.very.high.10))


## Figure 5 in manuscript #
require(RColorBrewer)
col.vec=brewer.pal(n=3, name="Dark2")

jpeg("ecoimm_fig5.tiff", width=2000, height=1600, res=250, family="Arial")
par(mfrow=c(2,1))
par(oma=c(3,0,0,0))
par(mar=c(1,7,1,1))

plot(dist2.four ~ time.med, type="l", col=col.vec[1], lwd=3,xlab=NA,ylab=NA,
     ylim=c(0,1))
lines(dist2.switch ~ time.med.10, col=col.vec[2], lwd=3)
lines(dist2.one ~ time.med.10, lwd=3, col=col.vec[3])
legend(legend=c("Four-parameter polynomial", "Bang-Bang with one switch","Fixed"), col=col.vec,
       lty=c(1,1,1), lwd=c(3,3,3),x=300, y=1, bty="n")
mtext(paste("Optimal", "\n", "immune", "\n", "investment"), las=1, side=2, line=4.5, adj=0.5)
text(labels="a)",x=0,y=0.97)

plot(dist4.four ~ time.very.high, type="l", col=col.vec[1], lwd=3,xlab=NA,ylab=NA,
     ylim=c(0,1))
lines(dist4.switch ~ time.very.high.10, col=col.vec[2], lwd=3)
lines(dist4.one ~ time.very.high.10, lwd=3, col=col.vec[3])
mtext(paste("Optimal", "\n", "immune", "\n", "investment"), las=1, side=2, line=4.5, adj=0.5)
mtext("Time (d)", line=2.5, side=1, cex=1.25)
text(labels="b)",x=0,y=0.97)
dev.off()



