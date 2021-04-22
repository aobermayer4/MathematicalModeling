#Question 1
#Part A : What is the stable growth rate of the population as represented in matrix A?

#Generating the Leslie Matrix
A <- matrix(c(0,0.675,0,0,0,0,0.703,0.047,0,0,0,0,0.657,0.019,0,4.665,0,0,0.682,0.061,61.896,0,0,0,0.8091), nrow = 5, ncol = 5)

#Find the eigenvalue and vector of the matrix (A)
EVs <- eigen(A)

#Assign largest eigenvalue to big R (Geometric Growth Factor>Stable population growth rate)
R <- abs(eigen(A)$values[1])
R

#Part B: Which changes the growth rate by a larger proportion, 
  #increasing the survival to juvenile stage to 100% or 
  #decreasing clutch size for the large mature size class by half?

#Juvenile Change
#Generate Leslie Matrix with increased survival of juvenile stage to 100%
AJ <- matrix(c(0, 1.000,0,0,0,0,0.703,0.047,0,0,0,0,0.657,0.019,0,4.665,0,0,0.682,0.061,61.896,0,0,0,0.8091), nrow = 5, ncol = 5)

#Find the eigenvalue and vector of the matrix (AJ)
EVsJ <- eigen(AJ)

#Assign largest eigenvalue to big R (Geometric Growth Factor>Stable population growth rate)
RJ <- abs(eigen(AJ)$values[1])
RJ

#Clutch size change
#Decrease clutch size for large mature size class by half
MCS <- (A[1,5]/2)

#Generate Leslie Matrix with large mature size class by half
AC <- matrix(c(0, 0.675,0,0,0,0,0.703,0.047,0,0,0,0,0.657,0.019,0,4.665,0,0,0.682,0.061,30.948,0,0,0,0.8091), nrow = 5, ncol = 5)

#Find the eigenvalue and vector of the matrix (AC)
EVsC <- eigen(AC)

#Assign largest eigenvalue to big R (Geometric Growth Factor>Stable population growth rate)
RC <- abs(eigen(AC)$values[1])
RC

#Find change in new R values from original
ChangeWithRJ <- eval(abs(R-RJ))
ChangeWithRJ
ChangeWithRC <- eval(abs(R-RC))
ChangeWithRC

#Part C: How much does mortality in stages 4 and 5 need to be reduced to achieve positive growth? 
  #(Assume the change in mortality is equal for the two stages). 

PC <- matrix(c(0,0.675,0,0,0,0,0.703,0.047,0,0,0,0,0.657,0.019,0,4.665,0,0,0.782,0.061,61.896,0,0,0,0.9091), nrow = 5, ncol = 5)
EVsP <- eigen(PC)
PCE <- abs(eigen(PC)$values[1])
PCE


i <- c(0.682, 0.8091)

x <- seq(0,5,0.1)
p <- 0.682
q <- 0.8091

for (i in x) {
  pv <- p + x
  qv <- q + x {
    for (j in pv) {
      A[4,4] <- j
      for (k in qv) {
        A[5,5] <- k {
          for (l in A) {
            
          }
        }
      }
    }
  }
}

for (i in pv) {
  for (j in qv) {
    
  }
}


for (i in x) {
  pv <- p + x
  qv <- q + x {
    for (pv )
    if abs(eigen(A)$values[1]) > 1 {
      break
    }
  }
  print(c(A[4,4], A[5,5]))
}

#Part D: Is there a reduction in mortality in stages 4 and 5 that results in the population doubling in 50 years? If so, what is it?


#Question 2
dem.stoch <- function(p0,p1,n0,times){
  n=n0
  for (i in 1:times){
    new_n = 0 #dependent on how many offspring each current individual has
    if (n[i]==0){
      n <- c(n, new_n)
    }else{
      for (j in 1:n[i]){
        luck = runif(1)
        if (luck < p0){
          new_n = new_n
        }else{
          if (luck < p0 + p1){
            new_n = new_n + 1
          }else{
            new_n = new_n + 2}
        }
      }      
      n <- c(n, new_n)
    }
  }
  return(n)
}

#Makes sure population cannot grow from zero
n0.4<-dem.stoch(0.25,0.5,4,200)
plot(n0.4,type="l",ylab="population",xlab="time",lwd=2)



n0.128<-dem.stoch(0.25,0.5,128,200)
plot(n0.128,ylim=c(0,max(n0.128)),type="l",ylab="population",xlab="time",lwd=2,col="darkorchid2")
lines(dem.stoch(0.25,0.5,4,200),ylab="population",xlab="time",lwd=2,col="black")
lines(dem.stoch(0.25,0.5,8,200),ylab="population",xlab="time",lwd=2,col="firebrick")
lines(dem.stoch(0.25,0.5,16,200),ylab="population",xlab="time",lwd=2,col="orange")
lines(dem.stoch(0.25,0.5,32,200),ylab="population",xlab="time",lwd=2,col="goldenrod3")
lines(dem.stoch(0.25,0.5,64,200),ylab="population",xlab="time",lwd=2,col="forestgreen")
legend(legend=c("4","8","16","32","64","128"),
       text.col=c("black","firebrick","orange","goldenrod3","forestgreen","darkorchid2"),
       lty=1,col=c("black","firebrick","orange","goldenrod3","forestgreen","darkorchid2"),
       lwd=2,x=0,y=180)




#Function for population growth
envr <- function(rbad, rgood, rneut, n0) {
  n <- n0 #Set initial population size outside the loop
  for (i in 1:200) {
    if (runif(1) < 0.5) { #Random switch between good and bad conditions
      R = rbad
    } else {
      if (runif(1) == 0.5) {
        R = rneut
      } else {
        if (runif(1) > 0.5) {
          R = rgood
        }
      }
    }
        n <- c(n, R*n[i]) #Calculate the next population size
      }
      return(n)
    }
    envrpop <- envr(rbad=(4/5),rgood=(5/4),rneut=1,n0=1250)
    
    plot(envrpop,type="l",ylab="Population Size",xlab="Time (years)", main="Logistic Growth with Environmental Variation")
    
#Logistic Model with Environmental Variation    
logistic <- function(rbad, rgood, rneut, k, n0, time){
  n <- n0
    for (t in 1:time){
      if (runif(1) < 0.33) { #Random switch between bad, neutral, and good conditions
        r = rbad
      } else {
        if (0.33 > runif(1) && runif(1) < 0.66) {
          r = rneut
        } else {
          if (runif(1) > 0.66) {
            r = rgood
          }
        }
      }
      nprime <- n[t] + r*n[t]*(k-n[t])/k
      if (nprime < 0) nprime = 0
      n <- c(n, nprime)
    }
  return(n)
}
    
envrlog <- logistic((4/5),(5/4),1,1000,2,20)
    
plot(envrlog,type="l",lwd=2,ylab="Population",xlab="Time (years)")
    
#Basic Logistic Model    
logistic <- function(r, k, n0, time){
  n <- n0
    for (t in 1:time){
      nprime <- n[t] + r*n[t]*(k-n[t])/k
      if (nprime < 0) nprime = 0
      n <- c(n, nprime)
    }
  return(n)
}
    
logR <- logistic(R,1000,2,20)
    
plot(logR,type="l",lwd=2,ylab="Population",xlab="Time (years)")