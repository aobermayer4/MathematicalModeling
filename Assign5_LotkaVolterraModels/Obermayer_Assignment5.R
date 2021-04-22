library(Ryacas)
library(deSolve)

#------------Question 1---------------#

#Species non-interacting, alone, with disturbance
#d-overall disturbance; h-exposure to hazard from disturbance
#The n1 and n2 alone expressions represent the affect of the hazard the disturbance 
#causes in relation to the overall disturbance occurring on the respective population 
#sizes while not interacting with each other.

n1_alone <- expression(k1*(r1-h1*d)/r1)
n1_alone <- ysym(n1_alone)
n2_alone <- expression(k2*(r2-h2*d)/r2)
n2_alone <- ysym(n2_alone)

#Species coexisting with disturbance
#The n1 and n2 coexisting expressions represent the affect of the disturbance, 
#as such described above, though showing the impact on each population while both
#are coexisting in the same space.

n1_coexist <- expression((k1*(r1-h1*d)/r1-a12*k2*(r2-h2*d)/r2)/(1-a12*a21))
n1_coexist <- ysym(n1_coexist)
n2_coexist <- expression((k2*(r2-h2*d)/r2-a21*k1*(r1-h1*d)/r1)/(1-a12*a21))
n2_coexist <- ysym(n2_coexist)

#Species invading when rare with disturbance
#The dn1 and dn2 rare expressions represent the change in population when either 
#species is invading when rare, while disturbance is playing an affect on this invasion.

dn1_rare <- expression(k1*(r1-h1*d)/r1-a12*k2*(r2-h2*d)/r2)
dn1_rare <- ysym(dn1_rare)
dn2_rare <- expression(k2*(r2-h2*d)/r2-a21*k1*(r1-h1*d)/r1)
dn2_rare <- ysym(dn2_rare)

#Solving for the minimum amount of disturbance that allows species 2 to increase when rare
#First, we are solving dn2_rare for the disturbance coefficient to find the minimum
#amount of disturbance that allows species 2 to increase when rare. Prior to this 
#minimum, species 1 is able to dominate and will crowd out species 2 in the absence
#of disturbance or with relatively little disturbance.

d_min <- yac_expr(y_rmvars(solve(dn2_rare, 'd')))
d_min


#Solving for the maximum amount of disturbance that species 1 can withstand and survive
#Next, we are solving our dn1_rare expression for the disturbance coefficient to
#find the maximum amount of disturbance that species 1 can withstand and survive.
#Leading up to this maximum, post-minimum species-1 disturbance, both of the species
#are coexisting with disturbance and species 2 is able to grow in abundance while
#species 1 is on its demise.


d_max <- yac_expr(y_rmvars(solve(dn1_rare, 'd')))
d_max

#Solving for the upper boundary of allowable disturbance
#Lastly, we are solving n2_alone for the disturbance coefficient to find the upper
#bound of disturbance allowable for species 2 to survive while not interacting with 
#species 1. While species 1 was dying off from too much disturbance, species 2 
#was able to grow in abundance until it has also reached its limit of too much disturbance.

d_ub <- yac_expr(y_rmvars(solve(n2_alone, 'd')))
d_ub


#Parameters

r1 <- 0.75
k1 <- 750
a12 <- 0
h1 <- 1
r2 <- 1
k2 <- 1000
a21 <- 2
h2 <- 0.75

#Solving for the d_min, d_max, and d_ub expressions defined above

eval(d_min)
eval(d_max)
eval(d_ub)

#Plot values

plot(NA, type = "l", xlim = c(0,1.4), ylim = c(0,800), xaxs = "i", yaxs = "i", xlab = "Disturbance", ylab = "Abundance", main = "Intermediate Disturbance Principle")
  n1_alone <- expression(k1*(r1-h1*d)/r1)
  d <- seq(0,eval(d_min),0.04)
  lines(eval(n1_alone)~d)
  n1_coexist <- expression((k1*(r1-h1*d)/r1-a12*k2*(r2-h2*d)/r2)/(1-a12*a21))
  n2_coexist <- expression((k2*(r2-h2*d)/r2-a21*k1*(r1-h1*d)/r1)/(1-a12*a21))
  d <- seq(eval(d_min),eval(d_max),0.35)
  lines(eval(n1_coexist)~d)
  lines(eval(n2_coexist)~d)
  n2_alone <- expression(k2*(r2-h2*d)/r2)
  d <- seq(eval(d_max),eval(d_ub),0.55)
  lines(eval(n2_alone)~d)
  
#The essential idea of the intermediate disturbance hypothesis tells us that an 
#intermediate level of disturbance is a good thing because it promotes species 
#through ensuring that one species cannot crowd out another. The plot generated
#is able to support this idea by showing that as the amount of disturbance in the
#environment increases, it was able to generate coexistence for a given period of time
#while the disturbance levels were manageable for each species, between 0.4 and 0.75.
#Although the plot does show the moment in time where both species seize to exist,
#it is able to give a range of disturbance levels that promote species diversity.
  
#When thinking about models of disturbance and competition between more than two species
#there are additional factors that may play a role. 
  
  
#In order for the intermediate
#disturbance hypothesis to contribute, at least one species has to have the ability
#to crowd other species out when not in the presence of disturbance. For instance, 
#that top species must be able to reproduce at a high enough rate to overwhelm its 
#surrounding and use up all the resources so no other competing species can survive
#in its presence. 
  
  
#------------Question 2---------------#
  
r <- 0.7   #Intrinsic rate of increase of prey
a <- 0.01  #Coefficient relating prey-capture rate to pred-prey-collision rate
b <- 0.05  #Coefficient describing the number of prey captures needed to produce a predator
d <- 0.12  #Death rate of predator
c <- 10    #Coefficient for leveling off

#Predator prey model
dn1 <- expression(r*n1*(k-n1)/k-c*(1-exp(-a*n1/c))*n2)
dn2 <- expression(b*c*(1-exp(-a*n1/c))*n2-d*n2)

#Coexistence equilibria
n1hat3 <- expression(-log((b*c-d)/b/c)/a*c)
n2hat3 <- expression(-r*log((b*c-d)/b/c)*c*b*(k*a+log((b*c-d)/b/c)*c)/k/a^2/d)

#k_min, found by solving n2hat3 for k
#Solving for prey because we want to know the min amount of prey to support the pred
#When k is just a little be higher that k_min it is a stable node
kmin <- -log((b*c-d)/b/c)/a*c

#Numerical solution
LVpred_satiate <- function (t, n, parms){
  with(as.list(parms), {
    dn1 <- r * n[1]*(k-n[1])/k - c*(1-exp(-a*n[1]/c))*n[2]
    dn2 <- b*c*(1-exp(-a*n[1]/c))*n[2] - d*n[2]
    list(c(dn1, dn2))
  })
}


j11 <- expression(r*(-log((b*c-d)/b/c)/a*c)*(k-(-log((b*c-d)/b/c)/a*c))/k-c*(1-exp(-a*(-log((b*c-d)/b/c)/a*c)/c))*(-r*log((b*c-d)/b/c)*c*b*(k*a+log((b*c-d)/b/c)*c)/k/a^2/d))
#j11 <- ysym(j11)
j12 <- expression(-log((b*c-d)/b/c)/a*c)
#j12 <- ysym(j12)
j21 <- expression(-r*log((b*c-d)/b/c)*c*b*(k*a+log((b*c-d)/b/c)*c)/k/a^2/d)
#j21 <- ysym(j21)
j22 <- expression(b*c*(1-exp(-a*(-log((b*c-d)/b/c)/a*c)/c))*(-r*log((b*c-d)/b/c)*c*b*(k*a+log((b*c-d)/b/c)*c)/k/a^2/d)-d*(-r*log((b*c-d)/b/c)*c*b*(k*a+log((b*c-d)/b/c)*c)/k/a^2/d))
#j22 <- ysym(j22)


jacob <- matrix(c(j11,j12,j21,j22),nrow = 2, byrow = T)
jacob <- ysym(jacob)

trace <- jacob[1,1] + jacob[2,2]

