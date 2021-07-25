library(splines)
library(ggplot2)

n <- 5000

# creating knots for B-spline basis
knots <- function(inner_knots) {
  sort(c(rep(range(inner_knots), 3), inner_knots))
}


# f(x|beta) = (phi1(x), phi2(x),...,phim(x)) * beta
f <- function(par, x ,inner_knots) {
  
  if(length(par) != length(inner_knots) + 2) {
    stop("none compaterble dimensions")
  }
  phi <- splineDesign(knots(inner_knots), x) # designmatrix
  phi %*% par
}


xx <- seq(0, 1000, length.out = n)
#xx <- rnorm(n)
#xx <- rnorm(10)
#xx <- runif(n, min = 0, max = 500)
inner_knotsxx <- seq(range(xx)[1], range(xx)[2], length.out = 3)
par0 <- rnorm(5)


pvaerdier <- function(x){
  f <- f(par0, x, inner_knotsxx)  #0.1 + 1.2*(x-0.5)^2 + 0.9*x^3  + 0.3*x^4 + 0.2*x^5
  
  exp(f)/(1 + exp(f))
}

yy <- rbinom(n, 1, pvaerdier(xx))

df <- data.frame(x = xx, y = yy, p  = pvaerdier(xx) )
df <- data.frame(x = xx, y = yy)
df$interval.x <- seq(range(df$x)[1],range(df$x)[2], length.out = 20)

library(dplyr)

group_df <- df %>% group_by(y, interval.x)
df2 <- group_df %>% summarise(x = mean(x), count = n())
df2 <- as.data.frame(df2)
str(df2)
str(df)
#My try
data1 <- data.frame(x = df$x[df$y == 1])
data0 <- data.frame(x = df$x[df$y == 0])
xx1 <- hist(data1$x,10,plot = F)$counts
xx0<- hist(data0$x,10,plot = F)$counts
data1 <- data.frame(x = xx1, y = 1)
data0 <- data.frame(x = xx0, y = 0)
df<-merge(data1,data0,all=TRUE)
df <- data.frame(x = xx, y = yy)
#plot(yy ~ xx)

ggplot(df , aes(x = xx, y = yy)) + geom_point() + geom_smooth(method = "gam", formula = y ~ pvaerdier(x), se = F)
ggplot(df2 , aes(x = x, y = y)) + geom_point(aes(size = count)) + geom_smooth(method = "gam", formula = y ~ pvaerdier(x), se = F)

xx.interval <- seq(range(xx)[1], range(xx)[2],length.out = 20)
round(xx.interval, digits = -1)


