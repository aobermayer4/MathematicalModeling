library(Ryacas)
library(deSolve)

#-----Base SIR-----#

#SIR compartmental model function
SIR <- function (t, n, parms){
  with(as.list(parms), {
    S <- -b*n[1]*n[2]
    I <- b*n[1]*n[2] - g*n[2]
    R <- g*n[2]
    list(c(S, I, R))
  })
}

N <- 472      #Initial Population
b <- 0.000318 #Rate of infection
g <- 0.0175   #Recovery Rate
s0 <- 460     #Initial susceptible
i0 <- 12      #Initial Infected
r0 <- 0       #Initial recovered
initialN <- c(s0,i0,r0) #Initial population vector

#SIR ode function showing number of cases over 90 days
SIR.out <- ode(y = initialN, times = seq(0,90,by=1), func = SIR,
               parms=c(N=N,b=b,g=g),method="ode45")
#Plotting SIR Model
SIRplot <- plot(SIR.out[,2]~SIR.out[,1], #Plotting susceptible
                type="l", lwd=2, col="tomato3",
                xlab="t (days)",ylab="Number of Cases", 
                xaxs="i", yaxs="i",
                xlim=c(0,90), ylim=c(0,500),
                main="Base Numerical SIR Model for Ebola")
lines(SIR.out[,3]~SIR.out[,1], lwd=2, col="turquoise4") #Plotting infected
lines(SIR.out[,4]~SIR.out[,1], lwd=2, col="violetred4") #Plotting recovered


#-----Dimensionless SIR-----#

#SIR dimensionless compartmental model function
SIRd <- function (t, n, parms){
  with(as.list(parms), {
    S <- -b*n[1]*n[2]
    I <- b*n[1]*n[2] - g*n[2]
    R <- g*n[2]
    list(c(S, I, R))
  })
}

N <- 1
b <- 0.2
g <- 0.1
s0 <- 0.95
i0 <- 0.05
r0 <- 0
#Initial parameters changed because model is now based off of percentages of the population
initialN <- c(s0, i0, r0)

SIR.outd <- ode(y = initialN, times = seq(0,90,1), func = SIRd, 
                parms = c(N=N, b=b, g=g), method='ode45')
SIRdplot <- plot(SIR.outd[,2]~SIR.outd[,1], 
                 type='l', lwd=2, col='tomato3', 
                 xlab='t (days)', ylab='Number of Cases', 
                 xaxs="i", yaxs="i",
                 xlim=c(0,90), ylim=c(0,1),
                 main="Dimensionless SIR Model for Ebola")
lines(SIR.outd[,3]~SIR.outd[,1], lwd=2, col='turquoise4')
lines(SIR.outd[,4]~SIR.outd[,1], lwd=2, col='violetred4')


#-----Vaccine SIR-----#

#SIR dimensionless compartmental model function with vaccine introduced
SIRv <- function (t, n, parms){
  with(as.list(parms), {
    S <- -b*n[1]*n[2] - v*n[2] #subtracting vaccinated individuals
    I <- b*n[1]*n[2] - g*n[2]
    R <- g*n[2] + v*n[2] #adding vaccinated individuals
    list(c(S, I, R))
  })
}

N <- 1
b <- 0.2
g <- 0.1
v <- 0 #rate of vaccination, beginning at 0
s0 <- 0.95
i0 <- 0.05
r0 <- 0
initialN <- c(s0, i0, r0)
SIRvaccs <- par(mfrow=c(1,2))

#rate of vaccination = 0
SIR.outv0 <- ode(y = initialN, times = seq(0,90,1), func = SIRv, 
                 parms = c(N=N, b=b, g=g, v=v), method='ode45')
SIRv0 <- plot(SIR.outv0[,2]~SIR.outv0[,1], 
              type='l', lwd=2, col='tomato3', 
              xlab='t (days)', ylab='Number of Cases', 
              xaxs = "i", yaxs = "i",
              xlim=c(0,90), ylim=c(0,1),
              main="Dimensionless SIR Model with Vaccine")
lines(SIR.outv0[,3]~SIR.outv0[,1], lwd=2, col='turquoise4')
lines(SIR.outv0[,4]~SIR.outv0[,1], lwd=2, col='violetred4')

#rate of vaccination = 0.005
v <- 0.005
SIR.outv005 <- ode(y = initialN, times = seq(0,90,1), func = SIRv, 
                   parms = c(N=N, b=b, g=g, v=v), method='ode45')
SIRv005 <- plot(SIR.outv005[,2]~SIR.outv005[,1], 
                type='l', lwd=2, col='tomato3', 
                xlab='t (days)', ylab='Number of Cases', 
                xaxs = "i", yaxs = "i",
                xlim=c(0,90), ylim=c(0,1))
                #main="Dimensionless SIR Model with Vaccine for Ebola")
lines(SIR.outv005[,3]~SIR.outv005[,1], lwd=2, col='turquoise4')
lines(SIR.outv005[,4]~SIR.outv005[,1], lwd=2, col='violetred4')

SIRvaccs2 <- par(mfrow=c(1,2))
#rate of vaccination = 0.01
v <- 0.01
SIR.outv01 <- ode(y = initialN, times = seq(0,90,1), func = SIRv, 
                  parms = c(N=N, b=b, g=g, v=v), method='ode45')
SIRv01 <- plot(SIR.outv01[,2]~SIR.outv01[,1], 
               type='l', lwd=2, col='tomato3', 
               xlab='t (days)', ylab='Number of Cases', 
               xaxs = "i", yaxs = "i",
               xlim=c(0,90), ylim=c(0,1),
               main="Dimensionless SIR Model with Vaccine")
lines(SIR.outv01[,3]~SIR.outv01[,1], lwd=2, col='turquoise4')
lines(SIR.outv01[,4]~SIR.outv01[,1], lwd=2, col='violetred4')

#rate of vaccination = 0.02
v <- 0.02
SIR.outv02 <- ode(y = initialN, times = seq(0,90,1), func = SIRv, 
                  parms = c(N=N, b=b, g=g, v=v), method='ode45')
SIRv02 <- plot(SIR.outv02[,2]~SIR.outv02[,1], 
               type='l', lwd=2, col='tomato3', 
               xlab='t (days)', ylab='Number of Cases', 
               xaxs = "i", yaxs = "i",
               xlim=c(0,90), ylim=c(0,1))
lines(SIR.outv02[,3]~SIR.outv02[,1], lwd=2, col='turquoise4')
lines(SIR.outv02[,4]~SIR.outv02[,1], lwd=2, col='violetred4')

SIRvaccs3 <- par(mfrow=c(1,2))
#rate of vaccination = 0.03
v <- 0.03
SIR.outv03 <- ode(y = initialN, times = seq(0,90,1), func = SIRv, 
                  parms = c(N=N, b=b, g=g, v=v), method='ode45')
SIRv03 <- plot(SIR.outv03[,2]~SIR.outv03[,1], 
               type='l', lwd=2, col='tomato3', 
               xlab='t (days)', ylab='Number of Cases', 
               xaxs = "i", yaxs = "i",
               xlim=c(0,90), ylim=c(0,1),
               main="Dimensionless SIR Model with Vaccine")
lines(SIR.outv03[,3]~SIR.outv03[,1], lwd=2, col='turquoise4')
lines(SIR.outv03[,4]~SIR.outv03[,1], lwd=2, col='violetred4')

#rate of vaccination = 0.06
v <- 0.06
SIR.outv06 <- ode(y = initialN, times = seq(0,90,1), func = SIRv, 
                  parms = c(N=N, b=b, g=g, v=v), method='ode45')
SIRv06 <- plot(SIR.outv06[,2]~SIR.outv06[,1], 
               type='l', lwd=2, col='tomato3', 
               xlab='t (days)', ylab='Number of Cases', 
               xaxs = "i", yaxs = "i",
               xlim=c(0,90), ylim=c(0,1))
lines(SIR.outv06[,3]~SIR.outv06[,1], lwd=2, col='turquoise4')
lines(SIR.outv06[,4]~SIR.outv06[,1], lwd=2, col='violetred4')