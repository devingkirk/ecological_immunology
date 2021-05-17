require(deSolve) 
require(DEoptim)
require(PBSddesolve)
require(RColorBrewer)


### Figure 1 ###
x <- seq(1,5000,1) # Time vector

# Bring in scale and shape values to make different distributions
values.table <- read.csv("~/results/scale_and_shape_values.csv")

# The hazard function is defined as the Probability Density Function / Survival Function

# Make PDF Matrix #
PDF.matrix <- matrix(nrow=5000,ncol=29)

for(i in 1:29){
  PDF.matrix[,i] <- dgamma(x, shape = values.table[i,3], scale = values.table[i,2])
}

# Make Survival Matrix #
survival.matrix <- matrix(nrow=5000,ncol=29)

for(i in 1:29){
  survival.matrix[,i] <- 1-pgamma(x, shape = values.table[i,3], scale = values.table[i,2])
}


# Make Hazard Matrix (pre - fixed for Inf and NaN) #
hazard.matrix.pre <- matrix(nrow=5000,ncol=29)
hazard.matrix.pre <- PDF.matrix/survival.matrix

# Set NaN to 0 #
hazard.matrix.pre2 <- matrix(nrow=5000,ncol=29)
hazard.matrix.pre2 <- ifelse(is.nan(hazard.matrix.pre),0,hazard.matrix.pre)

# Set Inf to 0 # 
hazard.matrix <- matrix(nrow=5000,ncol=29)
hazard.matrix <- ifelse(is.infinite(hazard.matrix.pre2),0,hazard.matrix.pre2)


##### Create Figure 1 #####
library(wesanderson)
col.vec=wes_palette(name="GrandBudapest2",n=29,type=c("continuous"))

jpeg("ecoimm_fig1.tiff", width=2000, height=1600, res=250, family="Arial")
par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
par(oma=c(3,7,0,0))
plot(survival.matrix[,1] ~ x, xlim=c(0,1000), type="n")
for(i in 1:29){
  lines(survival.matrix[,i] ~ x, lwd=2, col=col.vec[i])
}
mtext(text="Time (d)", side=1, line=2.5, cex=1.25)
mtext(text=paste("Survival", "\n", "probability"),side=2,line=4.5, las=1, cex=1.25, adj=0.5)
text(labels=paste("Most", "\n", "variability"),x=115,y=0.45)
text(labels=paste("Least", "\n", "variability"),x=340,y=0.9)
dev.off()





# For figure 2, bring in the calculated optimal results for each distribution (stored in CSV files)
#  Note that the model used within the optimizer is shown below (labeled: weighted.total.reproduction)

# Here, 'results' is for medium virulence, which we will look at more closely
# 'low' is low virulence, and 'high' is high virulence
results <- read.csv(file="~/results/fixed.results.virulence.00001.csv")
low <- read.csv(file="~/results/fixed.results.virulence.000001.csv")
high <- read.csv(file="~/results/fixed.results.virulence.0001.csv")


# Plot results from intermediate virulence, then all three levels on log scale
## Figure 2 in manuscript #
jpeg("ecoimm_fig2.tiff", width=2000, height=1600, res=250, family="Arial")
par(mfrow=c(2,1))
par(oma=c(3,0,0,0))
par(mar=c(1,7,1,1))

plot(results[,3] ~ results[,2], ylab=NA,
     type="p", ylim=c(0.06,0.15), xlab=NA, main=NA,
     pch=21,col="black",lwd=2, bg="grey",cex=1.5)
mtext(paste("Optimal", "\n", "immune", "\n", "investment"), las=1, side=2, line=4.5, adj=0.5)
text(labels="a)",x=0,y=0.145)

plot(log(results[,3]) ~ results[,2], ylab=NA, type="p", ylim=c(-5,3),
     xlab=NA,main=NA, pch=21,col="black", bg="grey",lwd=2)
points(log(low[,3]) ~ low[,2], col="black", bg="blue",pch=23, lwd=2)
points(log(high[,3]) ~ high[,2], col="black", bg="red", pch=24, lwd=2)
mtext(paste("ln(Optimal", "\n", "immune", "\n", "investment)"), las=1, side=2, line=4.5, adj=0.5)
mtext("Coefficient of variation in background survival", line=2.5, side=1, cex=1.25)
legend(x=.75,y=3.5,legend=c(expression(paste(alpha," = 0.000001")),expression(paste(alpha," = 0.00001")),
                           expression(paste(alpha," = 0.0001"))),pch=c(23,21,24),
       col=c("black","black","black"), pt.bg=c("blue","grey","red"),bty="n", pt.lwd=c(2,2,2))
text(labels="b)",x=0,y=2.5)
dev.off()




# For Figure 3, we need to define the model #
# function to weigh reproduction by probability of being alive ###
weighted.total.reproduction <- function(kappa,nat.hazard){
  TT <- seq(0,length(nat.hazard)-1,1) # time
  kappa<-kappa[1]; # immune investment
  r<-0.4; psi<-0.001; delta<-0.5; gamma<-0.1;  #parameters
  theta <- 0.1; lambda<- 50; #parameters
  P0 <- 1 # initial conditions
  I0 <- 1
  R0 <- 1
  E0 <- 50
  
  nhfunction = splinefun(cbind(TT,nat.hazard))
  
  induced.model<-function(t,x,params){
    r<-params[1]; psi<-params[2]; 
    gamma<-params[3]; delta<-params[4]; kappa<-params[5];theta<-params[6];
    lambda<-params[7]; 
    
    P<-x[1];
    I<-x[2]; 
    R<-x[3];
    E<-x[4];
    Cumu.hazard<-x[5];
    
    dP <- r*P - psi*I*P
    dI <- kappa*P*E-gamma*I
    dR <- delta*E-theta*R
    dE <- lambda - kappa*P*E - delta*E
    
    nh.value = nhfunction(t); # calculate hazard from natural mortality
    for(i in 1:length(nh.value)){
      if(nh.value[i]<0){ nh.value[i]=0 }}
    
    dCumu.hazard <- 0.00001*P + nh.value  
    
    list(c(dP,dI,dR,dE,dCumu.hazard))
  }
  
  out <- as.data.frame(lsoda(y=c(P=P0,I=I0,R=R0,E=E0,Cumu.hazard=0), times=TT, induced.model, parms=c(r,psi,gamma,delta,kappa,theta,lambda)))
  
  out[,7] <- NA
  out[,8] <- NA
  
  
  # Columns: 1 = Time; 2 = Parasites, 3 = Immune System; 4 = Reproduction; 5 = Energy; 6 = Cumulative Hazard; 7 = Survival Function
  
  colnames(out) <- c("time","P","I","R","E","Cumulative Hazard","Survival Function","Fitness")
  
  # Create survival distribution
  out[,7] <- exp(-out[,6])
  
  #weight reproduction by survival
  out[,8]<- out[1:5000,4]*out[1:5000,7]  
   
   return(out)
}




###### For three levels of variability (low = #1, int = #27, high = #29), we want to simulate dynamics for their optimal strategy, a sub-optimal strategy, and a supra-optimal strategy.
# Run simulations #
# Optimal strategy is from the optimizations.
# Setting sub-optimal strategy as Kappa = 0.01 for all.
# Setting supra-optimal strategy as Kappa = 1 for all.

# Low variability, i.e., hazard matrix # 1 
low.optimal <- weighted.total.reproduction(results[1,3],hazard.matrix[,1])
low.sub <- weighted.total.reproduction(0.01,hazard.matrix[,1])
low.supra <- weighted.total.reproduction(1,hazard.matrix[,1])

# Intermediate variability, i.e., hazard matrix # 27
int.optimal <- weighted.total.reproduction(results[27,3],hazard.matrix[,27])
int.sub <- weighted.total.reproduction(0.01,hazard.matrix[,27])
int.supra <- weighted.total.reproduction(1,hazard.matrix[,27])

# High variability, i.e., hazard matrix # 29
high.optimal <- weighted.total.reproduction(results[29,3],hazard.matrix[,29])
high.sub <- weighted.total.reproduction(0.01,hazard.matrix[,29])
high.supra <- weighted.total.reproduction(1,hazard.matrix[,29])



### Also want to calculate relative hazard from natural hazard and parasite-induced hazard ###
## Hazard from parasite at time = t is 0.00001*P(t)
## Hazard from natural mortality at time = t is from the hazard function.
low.optimal.relative.hazards <- hazard.matrix[,1]/(low.optimal[,2]*0.00001)
low.sub.relative.hazards <- hazard.matrix[,1]/(low.sub[,2]*0.00001)
low.supra.relative.hazards <- hazard.matrix[,1]/(low.supra[,2]*0.00001)

int.optimal.relative.hazards <- hazard.matrix[,27]/(int.optimal[,2]*0.00001)
int.sub.relative.hazards <- hazard.matrix[,27]/(int.sub[,2]*0.00001)
int.supra.relative.hazards <- hazard.matrix[,27]/(int.supra[,2]*0.00001)

high.optimal.relative.hazards <- hazard.matrix[,29]/(high.optimal[,2]*0.00001)
high.sub.relative.hazards <- hazard.matrix[,29]/(high.sub[,2]*0.00001)
high.supra.relative.hazards <- hazard.matrix[,29]/(high.supra[,2]*0.00001)




### Now will make 5 x 3 figure.
# Columns 1, 2, and 3 are for low, intermediate and high variability respectively
# For each, plot dynamics for optimal, sub-optimal, and supra-optimal strategies
# Although all run for 5000 time steps, just plot first 1000

jpeg("ecoimm_fig3.tiff", width=3000, height=3000, res=350, family="Arial")
par(mfrow=c(5,3))
par(oma=c(4,10,6,1))
par(mar=c(2,2,1,1))

######### Background hazard rates #
plot(hazard.matrix[,1] ~ low.optimal$time, type="l", col="black", xlim=c(0,1000), ylim=c(0,.01),lwd=2)
arrows(x0=210, x1=210,y0=0.009,y1=0.0103, lwd=1.5, length=0.1) # add arrow to show it goes way beyond y-axis scale
mtext(text="Low", side=3, line=1)
mtext(text=paste("Background", "\n", "hazard rate  "),side=2,line=6.5, las=1, cex=1, adj=0.5)
text(label="a)",x=950, y=0.009, cex=1.5)


plot(hazard.matrix[,27] ~ int.optimal$time, type="l", col="black", xlim=c(0,1000), ylim=c(0,.01),lwd=2)
mtext(text="Intermediate", side=3, line=1)
mtext(text="Variability in host lifespan", side=3, line=4,cex=1.5)
text(label="b)",x=50, y=0.0092, cex=1.5)


plot(hazard.matrix[,29] ~ high.optimal$time, type="l", col="black", xlim=c(0,1000), ylim=c(0,.01),lwd=2)
mtext(text="High", side=3, line=1)
text(label="c)",x=950, y=0.009, cex=1.5)


######### PARASITE DYNAMICS #######3
plot(low.optimal$P ~ low.optimal$time, type="l", col="#468B00", xlim=c(0,1000), ylim=c(0,650),lwd=2)
lines(low.sub$P ~ low.sub$time, col="#00468B", lwd=2)
lines(low.supra$P ~ low.sub$time, col="#8B0000",  lwd=2)
mtext(text=paste("Parasite", "\n", "abundance "),side=2,line=6.5, las=1, cex=1, adj=0.5)
legend(legend=c("Optimal","Sub-optimal","Supra-optimal"),lty=c(1,1,1),lwd=c(2,2,2),
       col=c("#468B00","#00468B","#8B0000"), x=200, y=700, bty="n",cex=1.3)
text(label="d)",x=950, y=585, cex=1.5)


plot(int.optimal$P ~ int.optimal$time, type="l", col="#468B00", xlim=c(0,1000), ylim=c(0,650),lwd=2)
lines(int.sub$P ~ int.sub$time, col="#00468B", lwd=2)
lines(int.supra$P ~ int.sub$time, col="#8B0000",  lwd=2)
text(label="e)",x=950, y=585, cex=1.5)


plot(high.optimal$P ~ high.optimal$time, type="l", col="#468B00", xlim=c(0,1000), ylim=c(0,650),lwd=2)
lines(high.sub$P ~ high.sub$time, col="#00468B", lwd=2)
lines(high.supra$P ~ high.sub$time, col="#8B0000",  lwd=2)
text(label="f)",x=950, y=585, cex=1.5)




#### Survival curves ##
plot(low.optimal$`Survival Function` ~ low.optimal$time, type="l", col="#468B00", xlim=c(0,1000), ylim=c(0,1),lwd=2)
lines(low.sub$`Survival Function` ~ low.sub$time, col="#00468B", lwd=2)
lines(low.supra$`Survival Function` ~ low.sub$time, col="#8B0000",  lwd=2)
mtext(text=paste("Survival"),side=2,line=4.5, las=1, cex=1)
text(label="g)",x=950, y=0.9, cex=1.5)



plot(int.optimal$`Survival Function` ~ int.optimal$time, type="l", col="#468B00", xlim=c(0,1000), ylim=c(0,1),lwd=2)
lines(int.sub$`Survival Function` ~ int.sub$time, col="#00468B", lwd=2)
lines(int.supra$`Survival Function` ~ int.sub$time, col="#8B0000",  lwd=2)
text(label="h)",x=950, y=0.9, cex=1.5)

plot(high.optimal$`Survival Function` ~ high.optimal$time, type="l", col="#468B00", xlim=c(0,1000), ylim=c(0,1),lwd=2)
lines(high.sub$`Survival Function` ~ high.sub$time, col="#00468B", lwd=2)
lines(high.supra$`Survival Function` ~ high.sub$time, col="#8B0000",  lwd=2)
text(label="i)",x=950, y=0.9, cex=1.5)


# Reproduction dynamics #
plot(low.optimal$R ~ low.optimal$time, type="l", col="#468B00", xlim=c(0,1000), ylim=c(0,300),lwd=2)
lines(low.sub$R ~ low.sub$time, col="#00468B", lwd=2)
lines(low.supra$R ~ low.sub$time, col="#8B0000",  lwd=2)
mtext(text=paste("Reproductive", "\n", "cell abundance "),side=2,line=6.5, las=1, cex=1, adj=0.5)
text(label="j)",x=950, y=270, cex=1.5)


plot(int.optimal$R ~ int.optimal$time, type="l", col="#468B00", xlim=c(0,1000), ylim=c(0,300),lwd=2)
lines(int.sub$R ~ int.sub$time, col="#00468B", lwd=2)
lines(int.supra$R ~ int.sub$time, col="#8B0000",  lwd=2)
text(label="k)",x=950, y=270, cex=1.5)

plot(high.optimal$R ~ high.optimal$time, type="l", col="#468B00", xlim=c(0,1000), ylim=c(0,300),lwd=2)
lines(high.sub$R ~ high.sub$time, col="#00468B", lwd=2)
lines(high.supra$R ~ high.sub$time, col="#8B0000",  lwd=2)
text(label="l)",x=950, y=270, cex=1.5)


# Cumulative fitness over time #
plot(cumsum(low.optimal$Fitness) ~ low.optimal$time, type="l", col="#468B00",xlim=c(0,1000),ylim=c(0,23000),lwd=2)
lines(cumsum(low.sub$Fitness) ~ low.sub$time, col="#00468B", lwd=2)
lines(cumsum(low.supra$Fitness) ~ low.sub$time, col="#8B0000",  lwd=2)
mtext(text=paste("Cumulative", "\n", "fitness"),side=2,line=6.5, las=1, cex=1, adj=0.5)
mtext(text="Time (d)", side=1, line=3)
text(label="m)",x=50, y=20700, cex=1.5)


plot(cumsum(int.optimal$Fitness) ~ int.optimal$time, type="l", col="#468B00",xlim=c(0,1000), ylim=c(0,23000),lwd=2)
lines(cumsum(int.sub$Fitness) ~ int.sub$time, col="#00468B", lwd=2)
lines(cumsum(int.supra$Fitness) ~ int.sub$time, col="#8B0000",  lwd=2)
mtext(text="Time (d)", side=1, line=3)
text(label="n)",x=50, y=20700, cex=1.5)


plot(cumsum(high.optimal$Fitness) ~ high.optimal$time, type="l", col="#468B00",xlim=c(0,1000),ylim=c(0,23000), lwd=2)
lines(cumsum(high.sub$Fitness) ~ high.sub$time, col="#00468B", lwd=2)
lines(cumsum(high.supra$Fitness) ~ high.sub$time, col="#8B0000",  lwd=2)
mtext(text="Time (d)", side=1, line=3)
text(label="o)",x=50, y=20700, cex=1.5)
dev.off()