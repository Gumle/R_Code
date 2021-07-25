library("profvis")
setwd("~/Documents/Statistics - KU/Computational Statistics/Assignments/Assignment 1 - Density smoothing")
source("Assignment-1-functions.R")
F12 <- log(read.table("infrared.txt", header = T)$F12)
h <- 1
x <- rnorm(10000)
EKernel <- function(x){
  ifelse(-1 <= x & x<= 1,(3/4)*(1-x^2), 0) 
}
f <- cv2(x)

#Density estimations implementation testing:

#Loop
profvis(kdens_loop(x,h))
#Sapply
profvis(kdens_sapply(x,h))
#Bin 
profvis(kdens_bin(x,h))

source("profilingFunctions.R")


profvis(kernDens_bin(x, h = 1))


profvis(cv(x = x,h = 1, kernel = EKernel))
profvis(f2norm(x,h=1))
profvis(fSum(x,h=1))
