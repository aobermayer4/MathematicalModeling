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
v <- c(0,0.005,0.01,0.02,0.03,0.04,0.06) #rate of vaccination, beginning at 0
s0 <- 0.95
i0 <- 0.05
r0 <- 0
initialN <- c(s0, i0, r0)

SIR.outv0 <- ode(y = initialN, times = seq(0,90,1), func = SIRv,parms = c(N=N, b=b, g=g, v=0), method='ode45')
SIR.outv0.005 <- ode(y = initialN, times = seq(0,90,1), func = SIRv,parms = c(N=N, b=b, g=g, v=0.005), method='ode45')
SIR.outv0.01 <- ode(y = initialN, times = seq(0,90,1), func = SIRv,parms = c(N=N, b=b, g=g, v=0.01), method='ode45')
SIR.outv0.02 <- ode(y = initialN, times = seq(0,90,1), func = SIRv,parms = c(N=N, b=b, g=g, v=0.02), method='ode45')
SIR.outv0.03 <- ode(y = initialN, times = seq(0,90,1), func = SIRv,parms = c(N=N, b=b, g=g, v=0.03), method='ode45')
SIR.outv0.04 <- ode(y = initialN, times = seq(0,90,1), func = SIRv,parms = c(N=N, b=b, g=g, v=0.04), method='ode45')
SIR.outv0.06 <- ode(y = initialN, times = seq(0,90,1), func = SIRv,parms = c(N=N, b=b, g=g, v=0.06), method='ode45')


plot(NA,type='l', lwd=2, xlab='t (days)', ylab='Number of Cases',
     xaxs = "i", yaxs = "i",
     xlim=c(0,90), ylim=c(0,1),
     main="Dimensionless SIR Model with Vaccine")
lines(SIR.outv0[,2]~SIR.outv0[,1], lwd=2, col='lightpink3')
lines(SIR.outv0[,3]~SIR.outv0[,1], lwd=2, col='skyblue2')
lines(SIR.outv0[,4]~SIR.outv0[,1], lwd=2, col='lightgoldenrod')
lines(SIR.outv0.005[,2]~SIR.outv0.005[,1], lwd=2, col='hotpink')
lines(SIR.outv0.005[,3]~SIR.outv0.005[,1], lwd=2, col='steelblue1')
lines(SIR.outv0.005[,4]~SIR.outv0.005[,1], lwd=2, col='gold')
lines(SIR.outv0.01[,2]~SIR.outv0.01[,1], lwd=2, col='hotpink4')
lines(SIR.outv0.01[,3]~SIR.outv0.01[,1], lwd=2, col='turquoise2')
lines(SIR.outv0.01[,4]~SIR.outv0.01[,1], lwd=2, col='goldenrod3')
lines(SIR.outv0.02[,2]~SIR.outv0.02[,1], lwd=2, col='firebrick')
lines(SIR.outv0.02[,3]~SIR.outv0.02[,1], lwd=2, col='turquoise4')
lines(SIR.outv0.02[,4]~SIR.outv0.02[,1], lwd=2, col='darkorange1')
lines(SIR.outv0.03[,2]~SIR.outv0.02[,1], lwd=2, col='firebrick4')
lines(SIR.outv0.03[,3]~SIR.outv0.03[,1], lwd=2, col='royalblue3')
lines(SIR.outv0.03[,4]~SIR.outv0.03[,1], lwd=2, col='chocolate')
lines(SIR.outv0.04[,2]~SIR.outv0.04[,1], lwd=2, col='maroon')
lines(SIR.outv0.04[,3]~SIR.outv0.04[,1], lwd=2, col='mediumblue')
lines(SIR.outv0.04[,4]~SIR.outv0.04[,1], lwd=2, col='orange4')
lines(SIR.outv0.06[,2]~SIR.outv0.02[,1], lwd=2, col='darkorchid4')
lines(SIR.outv0.06[,3]~SIR.outv0.06[,1], lwd=2, col='navyblue')
lines(SIR.outv0.06[,4]~SIR.outv0.06[,1], lwd=2, col='saddlebrown')

##In this case, vaccination was incorporated into SIR model to understand the how the S , I, R were affected with introduction of vaccine.
##The individuals that were susceptible to Ebola reduce gradually and plateaus after 45 days, while the cases for infected individuals is highest around 30 days and then decreases.
##The recovery rate increases also increase and about 80% of individuals are recovered when there is no vaccination. But with increase in rate of vaccination the susceptibility of individuals contracting Ebola decreases and recovery rate increases (with vaccination 6% rate the rate of recovery reaches about 90%) 



