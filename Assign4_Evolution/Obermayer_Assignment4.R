drift_ns <- function(w11,w12,w22,p0,n,time){
  p <- p0
  for (t in 1:time){
    a1 <- rbinom(2*n,size=1,prob=p[t])
    pd <- sum(a1)/(2*n) #p after drift
    wbar = pd^2*w11+2*pd*(1-pd)*w12+(1-pd)^2*w22
    pprime = (pd*w11 + (1-pd)*w12)*pd/wbar
    p <- c(p, pprime)
  }
  return(p)
}
#set.seed(1)
p0 <- 0.5                 #initial frequency of A1
time <- 100               #Number of generations
popsize <- c(8,16,32,64)        #Population size, assumed to be constant here
e8 <- NULL
e16 <- NULL
e32 <- NULL
e64 <- NULL
colors<-rainbow(18)


for (i in 1:18){
  e8 <- cbind(e8,drift_ns(w11=0.5,w12=1,w22=0.5,p0=p0,n=8,time=time))
}
e8 <- plot(e8[,1],col=colors[1],lwd=1.5,type="l",ylim=c(0,1),xlab="generation",ylab="proportion A1")
        for (i in 2:18){
          lines(e8[,i],col=colors[i],lwd=1.5)
        }
        
for (i in 1:18){
  e16 <- cbind(e16,drift_ns(w11=0.5,w12=1,w22=0.5,p0=p0,n=16,time=time))
}
e16 <- plot(e16[,1],col=colors[1],lwd=1.5,type="l",ylim=c(0,1),xlab="generation",ylab="proportion A1")
        for (i in 2:18){
          lines(e16[,i],col=colors[i],lwd=1.5)
        }
        
for (i in 1:18){
  e32 <- cbind(e32,drift_ns(w11=0.5,w12=1,w22=0.5,p0=p0,n=32,time=time))
}
e32 <- plot(e32[,1],col=colors[1],lwd=1.5,type="l",ylim=c(0,1),xlab="generation",ylab="proportion A1")
        for (i in 2:18){
          lines(e32[,i],col=colors[i],lwd=1.5)
        }
        
for (i in 1:18){
  e64 <- cbind(e64,drift_ns(w11=0.5,w12=1,w22=0.5,p0=p0,n=64,time=time))
}
e64 <- plot(e64[,1],col=colors[1],lwd=1.5,type="l",ylim=c(0,1),xlab="generation",ylab="proportion A1")
        for (i in 2:18){
          lines(e64[,i],col=colors[i],lwd=1.5)
        }
        
#figure
colors<-rainbow(18)
plot(e[,1],col=colors[1],lwd=1.5,type="l",ylim=c(0,1),xlab="generation",ylab="proportion A1")
for (i in 2:18){
  lines(e[,i],col=colors[i],lwd=1.5)
}

for (n in popsize){
  for (i in 1:18){
    e <- cbind(e,drift_ns(w11=0.5,w12=1,w22=0.5,p0=p0,n=n,time=time))
    colors<-rainbow(18)
    plot(e[,1],col=colors[1],lwd=1.5,type="l",ylim=c(0,1),xlab="generation",ylab="proportion A1")
    for (i in 2:18){
      lines(e[,i],col=colors[i],lwd=1.5)
    }
  }
}

e4
e8
e16
e32
e64
e128