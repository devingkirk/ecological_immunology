require(deSolve) 
require(DEoptim)
require(PBSddesolve)
require(RColorBrewer)


# For figure 4, bring in the calculated optimal results for each distribution (stored in CSV files)

# Here, 'results' is for medium virulence, which we will look at more closely
# 'low' is low virulence, and 'high' is high virulence
results <- read.csv(file="~/results/fixed.results.virulence.00001.csv")
low <- read.csv(file="~/results/fixed.results.virulence.000001.csv")
high <- read.csv(file="~/results/fixed.results.virulence.0001.csv")

# Bring in scale and shape values to make different distributions
values.table <- read.csv("~/results/scale_and_shape_values.csv")

x <- seq(1,5000,1) # Time vector

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


# function to weigh reproduction by probability of being alive ###
weighted.total.reproduction <- function(kappa,nat.hazard){
  TT <- seq(0,length(nat.hazard)-1,1)
  kappa<-kappa[1];
  r<-0.4; psi<-0.001; delta<-0.5; gamma<-0.1;  
  theta <- 0.1; lambda<- 50; 
  P0 <- 1
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
    
    nh.value = nhfunction(t);
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
  
  repro <- out[1:5000,4]*out[1:5000,7]   
  return(-(sum(repro)))
  
}


###### For three levels of variability (low = #1, int = #27, high = #29), we want calculate fitness for a range of kappa values #

kappa.range <- seq(0.005,1,0.005)

low.fitness.surface <- c()
int.fitness.surface <- c()
high.fitness.surface <- c()

for(i in 1:length(kappa.range)){
  low.fitness.surface[i] <- weighted.total.reproduction(kappa.range[i],hazard.matrix[,1])
  int.fitness.surface[i] <- weighted.total.reproduction(kappa.range[i],hazard.matrix[,27])
  high.fitness.surface[i] <- weighted.total.reproduction(kappa.range[i],hazard.matrix[,29])
}



require(wesanderson)
col.vec=wes_palette(name="BottleRocket2",n=3,type=c("discrete"))

jpeg("ecoimm_fig4.tiff", width=2000, height=1400, res=250, family="Arial")
par(mfrow=c(1,1))
par(oma=c(4,8,1,1))
par(mar=c(1,1,1,1))

######### PARASITE DYNAMICS #######3
plot(-low.fitness.surface ~ kappa.range, type="l", col=col.vec[1], xlim=c(0,1), ylim=c(20500,22500),lwd=4)
lines(-int.fitness.surface ~ kappa.range, lwd=4, col=col.vec[2])
lines(-high.fitness.surface ~ kappa.range, lwd=4, col=col.vec[3])
abline(v=results[1,3], col=col.vec[1], lwd=1, lty=2)
abline(v=results[27,3], col=col.vec[2], lwd=1, lty=2)
abline(v=results[29,3], col=col.vec[3], lwd=1, lty=2)
mtext(text=paste("Cumulative", "\n", "fitness"),side=2,line=5.5, las=1, cex=1.25, adj=0.5)
mtext(text=expression(paste("Immune investment (",kappa,")")), side=1, line=3, cex=1.25)
legend(legend=c("Low variability","Intermediate variability","High variability"),
       x= 0.5,y=22470, lwd=c(3,3,3), lty=c(1,1,1), col=col.vec, bty="n")
dev.off()

