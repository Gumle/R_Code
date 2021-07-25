library(splines) # splineDesign (B-splines)
library(LaplacesDemon) # logit and invlogit
library(Rcpp) # c++ 
library(profvis) # profiling
library(numDeriv) # numerical methods for grad and hess
library(ggplot2) 
library(ggpubr)
library(microbenchmark) # benchmarking
library(bench) 
library(gridExtra)
library(dplyr) # data manipulation
library(Matrix) # sparse matrix
setwd("~/Documents/Statistics - KU/Computational Statistics/Assignments/Assignment 3 - Optimization")
source("Debugging_and_tracing.R")
source("Assignment-3-functions.R")
library(simts)

# test of gradiant and hessain:

# Test of grad_H with build in grad as reference:
max(abs(grad_H(par0, x, y, lambda, inner_knots) - grad(function(par) H(par, x, y, lambda, inner_knots), par0 )))
# Test of hessian_H with buildin hessian as reference
max(abs(hessian_H(par0, x, y, lambda, inner_knots) - hessian(function(par) H(par, x, y,  lambda, inner_knots), par0 )))

#inner_knots <- seq(range(x)[1], range(x)[2], length.out = 98)
#par0 <- rnorm(100)

#comparison <- microbenchmark(times = 200,f,f.sparse,grad_H,grad_H.sparse,hessian_H,hessian_H.sparse )
#summary(comparison)


p <- microbenchmark(times = 200, H(par0, x, y, lambda, inner_knots), H.sparse(par0, x, y, lambda, inner_knots))
summary(p)
autoplot(p) + geom_jitter(aes(color = expr), alpha = 0.4) + aes(fill = I("gray"))+
  scale_color_manual(values = c("#354e75", "#884091")) + theme_bw() +
  theme(legend.position = "none")

p <- microbenchmark(times = 200, grad_H(par0, x, y, lambda, inner_knots), grad_H.sparse(par0, x, y, lambda, inner_knots))
summary(p)
autoplot(p) + geom_jitter(aes(color = expr), alpha = 0.4) + aes(fill = I("gray"))+
  scale_color_manual(values = c("#354e75", "#884091")) + theme_bw() +
  theme(legend.position = "none")

p <- microbenchmark(times = 200, hessian_H(par0, x, y, lambda, inner_knots), hessian_H.sparse(par0, x, y, lambda, inner_knots))
summary(p)
autoplot(p) + geom_jitter(aes(color = expr), alpha = 0.4) + aes(fill = I("gray"))+
  scale_color_manual(values = c("#354e75", "#884091")) + theme_bw() +
  theme(legend.position = "none")


Hfun <- function(par0) H(par0, x, y, lambda, inner_knots)
sparseHfun <- function(par0) H.sparse(par0, x, y, lambda, inner_knots)
GradH <- function(par0) grad_H(par0, x, y, lambda, inner_knots)
sparseGradH <- function(par0) grad_H.sparse(par0, x, y, lambda, inner_knots)
HessianH<- function(par0) hessian_H(par0, x, y, lambda, inner_knots)
SparseHessianH <- function(par0) hessian_H.sparse(par0, x, y, lambda, inner_knots)

p <- microbenchmark(times = 200, Hfun(par0),sparseHfun(par0),GradH(par0),sparseGradH(par0),HessianH(par0),SparseHessianH(par0))
library(knitr)
kable(summary(p))
autoplot(p) + geom_jitter(aes(color = expr), alpha = 0.4) + aes(fill = I("gray"))+
  scale_color_manual(values = c("#354e75", "#884091")) + theme_bw() +
  theme(legend.position = "none")
# Compare non-sparse and sparse

length(x)
length(y)
length(inner_knots)
system.time(GD(par0, H, grad_H, hessian_H, maxiter = 5000, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'grad', cb = NULL, x, y, lambda, inner_knots)$par)


system.time(CG(par0, H, grad_H, hessian_H, maxiter = 5000, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'grad', cb = NULL, x, y, lambda, inner_knots)$par)


system.time(Newton(par0, H, grad_H, hessian_H, maxiter = 500, d = 0.1, c = 0.1,gamma0 = 1, epsilon = 1e-5, stop = 'grad', cb = NULL, x, y, lambda, inner_knots)$par)

H_optim(par0)$par

########################################

######################################

# Fix par0, x, y, lamda and inner_knots:

##minimizing with newton 
par_New <- Newton(par0, H, grad_H, hessian_H, maxiter = 500, d = 0.1, c = 0.1,
                 gamma0 = 1, epsilon = 1e-5, stop = 'grad', cb = NULL, x, y, lambda, inner_knots)$par



##minimizing with gradient decent
par_GD <- GD(par0, H, grad_H, hessian_H, maxiter = 5000, d = 0.1, c = 0.1,
            gamma0 = 1, epsilon = 1e-5, stop = 'grad', cb = NULL, x, y, lambda, inner_knots)$par

GD
##minimizing with conjugent gradient decent
par_CG <- CG(par0, H, grad_H, hessian_H, maxiter = 5000, d = 0.1, c = 0.1,
            gamma0 = 1, epsilon = 1e-5, stop = 'grad', cb = NULL, x, y, lambda, inner_knots)$par


##minimizing with optim
par_optim <- H_optim(par0)$par

df_fit <- data.frame(x = x, y = y) 
size <- seq(range(df_fit$x)[1],range(df_fit$x)[2],1)

ggplot(df_fit , aes(x = x, y = y)) +  geom_point() +
  geom_smooth(method = "gam", formula = y ~ pvaerdier_Newton(x), size = 1.5,se = F, col = "dodgerblue4")+
  geom_smooth(method = 'gam', formula = y ~ pvaerdier_GD(x), size = 1.5, se = F, col = "firebrick4") +
  geom_smooth(method = 'gam', formula = y ~ pvaerdier_optim(x), se = F, col = "chocolate1") +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() + 
  ggtitle("Blue: Newton,  Red: GD, Orange: optim. ") +
  xlab("Temperature") + 
  ylab("Probability of death")


###################Convergence plots (tracing)


GD_tracer <- tracer(c("H_val", "grad_norm"), N = 50)
system.time(logit_GD <- GD(par0, H, grad_H, hessian_H, maxiter = 1400, d = 0.1, c = 0.1,
                           gamma0 = 1, epsilon = 1e-5, stop = 'grad', cb = GD_tracer$trace,
                           x, y, lambda, inner_knots))




CG_tracer <- tracer(c("H_val", "grad_norm"), N = 5)
system.time(logit_GD <- CG(par0, H, grad_H, hessian_H, maxiter = 1400, d = 0.1, c = 0.1,
                           gamma0 = 1, epsilon = 1e-5, stop = 'grad', cb = CG_tracer$trace,
                           x, y, lambda, inner_knots))


df_timeGD <- tail(summary(GD_tracer), -100)
df_timeCG <- tail(summary(CG_tracer), -20)


df_timeGD$.time <- df_timeGD$.time*1000 #ms
df_timeCG$.time <- df_timeCG$.time*1000 #ms




p1 <- ggplot(df_timeGD, aes(.time, grad_norm), add = "reg.line") +
  geom_point() +
  stat_smooth(method = "lm", se = F) +
  stat_regline_equation(formula = y ~ x, label.x = 400) +
  scale_y_log10() +
  xlab("Time (ms)") +
  ylab("norm(grad)") +
  ggtitle("GD") +
  theme_bw() 




p2 <- ggplot(df_timeCG, aes(.time, grad_norm)) +
  geom_point() +
  stat_smooth(method = "lm", se = F) +
  stat_regline_equation(formula = y ~ x, label.x = 100) +
  scale_y_log10() +
  xlab("Time (ms)") +
  ylab("norm(grad)") +
  ggtitle("CG") +
  theme_bw() 

grid.arrange(p1, p2, ncol = 2)

###### SIMULATED DATA ######

#####################################################################################################################################
#                                                                                                                                   #
#--------------------------------------------------------- CONVERGE SPEEED ---------------------------------------------------------#
#                                                                                                                                   #
#####################################################################################################################################
n = 10000

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
inner_knotsxx <- seq(range(xx)[1], range(xx)[2], length.out = 3)
par0 <- rnorm(5)
pvaerdier <- function(x)
{
  f <- f(par0, x, inner_knotsxx)  #0.1 + 1.2*(x-0.5)^2 + 0.9*x^3  + 0.3*x^4 + 0.2*x^5
  
  exp(f)/(1 + exp(f))
}
yy <- rbinom(n, 1, pvaerdier(xx))
xx <- sample(xx)
length(xx)
length(yy)
length(inner_knotsxx)
#
yy[1:500]
inner_knotsxx
Newton(par0, H, grad_H, hessian_H, maxiter = 500, d = 0.1, c = 0.1,  gamma0 = 1, epsilon = 1e-5, stop = 'func', cb = NULL, xx, yy, lambda, inner_knotsxx)
GD(par0, H, grad_H, hessian_H, maxiter = 10000, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'func', cb = NULL, xx, yy, lambda, inner_knotsxx)
CG(par0, H, grad_H, hessian_H, maxiter = 10000, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'func', cb = NULL, xx, yy, lambda, inner_knotsxx)

New.Sim$par
Newton.bench <- microbenchmark(Newton(par0, H, grad_H, hessian_H, maxiter = 500, d = 0.1, c = 0.1,  gamma0 = 1, epsilon = 1e-5, stop = 'grad', cb = NULL,
                                      xx[1:500], yy[1:500], lambda, inner_knotsxx),
                               Newton(par0, H, grad_H, hessian_H, maxiter = 500, d = 0.1, c = 0.1,  gamma0 = 1, epsilon = 1e-5, stop = 'grad', cb = NULL,
                                      xx[1:2500], yy[1:2500], lambda, inner_knotsxx),
                               Newton(par0, H, grad_H, hessian_H, maxiter = 500, d = 0.1, c = 0.1,  gamma0 = 1, epsilon = 1e-5, stop = 'grad', cb = NULL,
                                      xx[1:5000], yy[1:5000], lambda, inner_knotsxx),
                               Newton(par0, H, grad_H, hessian_H, maxiter = 500, d = 0.1, c = 0.1,  gamma0 = 1, epsilon = 1e-5, stop = 'grad', cb = NULL,
                                      xx, yy, lambda, inner_knotsxx), times = 2)
summary(Newton.bench)

GD.bench <- microbenchmark(GD(par0, H, grad_H, hessian_H, maxiter = 10000, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'grad',
                              cb = NULL, xx[1:500], yy[1:500], lambda, inner_knotsxx),
                           GD(par0, H, grad_H, hessian_H, maxiter = 10000, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'grad',
                              cb = NULL, xx[1:2500], yy[1:2500], lambda, inner_knotsxx),
                           GD(par0, H, grad_H, hessian_H, maxiter = 10000, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'grad',
                              cb = NULL, xx[1:5000], yy[1:5000], lambda, inner_knotsxx),
                           GD(par0, H, grad_H, hessian_H, maxiter = 10000, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'grad',
                              cb = NULL, xx, yy, lambda, inner_knotsxx), times = 3)

summary(GD.bench)

CG.bench <- microbenchmark(CG(par0, H, grad_H, hessian_H, maxiter = 10000, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'grad',
                              cb = NULL, xx[1:500], yy[1:500], lambda, inner_knotsxx),
                           CG(par0, H, grad_H, hessian_H, maxiter = 10000, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'grad',
                              cb = NULL, xx[1:2500], yy[1:2500], lambda, inner_knotsxx),
                           CG(par0, H, grad_H, hessian_H, maxiter = 10000, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'grad',
                              cb = NULL, xx[1:5000], yy[1:5000], lambda, inner_knotsxx),
                           CG(par0, H, grad_H, hessian_H, maxiter = 10000, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'grad',
                              cb = NULL, xx, yy, lambda, inner_knotsxx), times = 3)

summary(CG.bench)

Newton.time <-summary(Newton.bench)$mean
GD.time <-summary(GD.bench)$mean
CG.time <-summary(CG.bench)$mean



linesNewton<- data.frame("Samplesize" = c(500, 2500, 5000, 10000), 
                           "Time" = Newton.time)

linesGD <- data.frame("Samplesize" = c(500, 2500, 5000, 10000), 
                          "Time" = GD.time)
linesCG <- data.frame("Samplesize" = c(500, 2500, 5000, 10000), 
                          "Time" = CG.time)

ggplot(data = linesNewton, aes(Samplesize, Time)) + geom_line(col = "dodgerblue4") + geom_line(data = linesGD, aes(Samplesize, Time), col = "firebrick4") + geom_line(data = linesCG, aes(Samplesize, Time), col = "darkgreen")+ ggtitle("Blue: SGD,  Red: GD, Green: CG ")
ggplot(data = linesNewton, aes(Samplesize, Time)) + geom_line(col = "dodgerblue4") + geom_line(data = linesGD, aes(Samplesize, Time), col = "firebrick4") + geom_line(data = linesCG, aes(Samplesize, Time), col = "darkgreen")+ ggtitle("Blue: SGD,  Red: GD, Green: CG ")+ ggtitle("Blue: SGD,  Red: GD, Green: EM ")+ ylab("Time [log10]")+ scale_y_log10()

#####################################################################################################################################
#                                                                                                                                   #
#-------------------------------------------------------------- TRACER -------------------------------------------------------------#
#                                                                                                                                   #
#####################################################################################################################################

Newton.tracer <-  tracer(c("H_val", "grad_norm"), N = 5)
GD.tracer <-  tracer(c("H_val", "grad_norm"), N = 5)
CG.tracer <-  tracer(c("H_val", "grad_norm"), N = 5)
Newton(par0, H, grad_H, hessian_H, maxiter = 500, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'grad', xx, yy, lambda, inner_knotsxx, cb = Newton.tracer$trace)
GD(par0, H, grad_H, hessian_H, maxiter = 5000, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'grad', xx, yy, lambda, inner_knotsxx, cb = GD.tracer$trace)
CG(par0, H, grad_H, hessian_H, maxiter = 5000, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'grad', xx, yy, lambda, inner_knotsxx, cb = CG.tracer$trace)
newton_tracer_data <- summary(Newton.tracer)
gd_tracer_data <- summary(GD.tracer)
cg_tracer_data <- summary(CG.tracer)


ggplot(data = newton_tracer_data, aes(.time*1000, H_val)) + geom_point(col = "dodgerblue4") + 
  geom_point(data = gd_tracer_data, aes(.time*1000, H_val), col ="firebrick4") + 
  geom_point(data = cg_tracer_data, aes(.time*1000, H_val), col = "darkgreen") + ggtitle("Blue: Newton,  Red: GD, Green: CG. ")+ labs(x= "Time [milliseconds]", y = "H value") + scale_y_log10()

ggplot(data = newton_tracer_data, aes(.time*1000, grad_norm)) + geom_point(col = "dodgerblue4")  + geom_smooth(data = gd_tracer_data, aes(.time*1000, grad_norm), col = "firebrick4") + 
  geom_smooth(data = cg_tracer_data, aes(.time*1000, grad_norm), col = "darkgreen") + ggtitle("Blue: Newton,  Red: GD, Green: CG.") + labs(x= "Time [milliseconds]") + scale_y_log10()


p1 <- ggplot(data = newton_tracer_data, aes(.time*1000, H_val)) + geom_point(col = "dodgerblue4") +  geom_smooth(data = gd_tracer_data, aes(.time*1000, H_val), col ="firebrick4") + 
  geom_smooth(data = cg_tracer_data, aes(.time*1000, H_val), col = "darkgreen") + ggtitle("Blue: Newton,  Red: GD, Green: CG. ")+ labs(x= "Time [milliseconds]", y = "H value") + scale_y_log10()

p2 <- ggplot(data = newton_tracer_data, aes(.time*1000, grad_norm)) + geom_point(col = "dodgerblue4")  + geom_smooth(data = gd_tracer_data, aes(.time*1000, grad_norm), col = "firebrick4") + 
  geom_smooth(data = cg_tracer_data, aes(.time*1000, grad_norm), col = "darkgreen") + ggtitle("Blue: Newton,  Red: GD, Green: CG.") + labs(x= "Time [milliseconds]", y = "Gradient of H") + scale_y_log10()

grid.arrange(p1, p2, ncol = 2)
# samplesize and runtime (Newton) fixing par, lambda, c, d, and 
#(inner_knots will depend on specific sample)
# but designmatrix has equal size.

#### MORE SIMULATED DATA ####
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


xx <- seq(0, 10, length.out = n)
xx <- sample(1:1000, n, replace=T)
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


pvaerdier_Newtonxx <- function(x) {
  drop(invlogit(f(par_New, x, inner_knotsxx)))
}

#Gradient probabilites
pvaerdier_GDxx <- function(x) {
  drop(invlogit(f(par_GD, x, inner_knotsxx)))
}

#Conjugent gradient probabilites
pvaerdier_CGxx <- function(x) {
  drop(invlogit(f(par_CG, x, inner_knotsxx)))
}

#optim probabilities
pvaerdier_optimxx <- function(x) {
  drop(invlogit(f(par_optim, x, inner_knotsxx)))
}
source("Assignment-3-functions.R")
df <- data.frame(x = xx, y = yy, p  = pvaerdier(xx) )
ggplot(df , aes(x = xx, y = yy)) + geom_point() + geom_smooth(method = "gam", formula = y ~ pvaerdier(x), se = F)
par_optim = H_optim(par0)$par
ggplot(df , aes(x = xx, y = yy)) +  geom_point() +
  geom_smooth(method = "gam", formula = y ~ pvaerdier_Newtonxx(x), size = 1.5,se = F, col = "dodgerblue4")+
  geom_smooth(method = 'gam', formula = y ~ pvaerdier_GDxx(x), size = 1.5, se = F, col = "firebrick4")  +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() + 
  ggtitle("Blue: Newton,  Red: GD ") +
  xlab("Integer") + 
  ylab("Probability")





##minimizing with newton 
par_New <- Newton(par0, H, grad_H, hessian_H, maxiter = 500, d = 0.1, c = 0.1,
                  gamma0 = 1, epsilon = 1e-5, stop = 'grad', cb = NULL, x, y, lambda, inner_knots)$par


##minimizing with gradient decent
par_GD <- GD(par0, H, grad_H, hessian_H, maxiter = 5000, d = 0.1, c = 0.1,
             gamma0 = 1, epsilon = 1e-5, stop = 'grad', cb = NULL, x, y, lambda, inner_knots)$par


##minimizing with conjugent gradient decent
par_CG <- CG(par0, H, grad_H, hessian_H, maxiter = 5000, d = 0.1, c = 0.1,
             gamma0 = 1, epsilon = 1e-5, stop = 'grad', cb = NULL, x, y, lambda, inner_knots)$par

# My tracer function

Newton.tracer <-  tracer(c("H_val", "grad_norm"), N = 1)
GD.tracer <-  tracer(c("H_val", "grad_norm"), N = 10)
CG.tracer <-  tracer(c("H_val", "grad_norm"), N = 10)

Newton(par0, H, grad_H, hessian_H, maxiter = 500, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'grad', x, y, lambda, inner_knots, cb = Newton.tracer$trace)
GD(par0, H, grad_H, hessian_H, maxiter = 5000, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'grad', x, y, lambda, inner_knots, cb = GD.tracer$trace)
CG(par0, H, grad_H, hessian_H, maxiter = 5000, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'grad', x, y, lambda, inner_knots, cb = CG.tracer$trace)


newton_tracer_data <- summary(Newton.tracer)
gd_tracer_data <- summary(GD.tracer)
cg_tracer_data <- summary(CG.tracer)


ggplot(data = newton_tracer_data, aes(.time*1000, H_val)) + geom_point(col = "dodgerblue4") + 
  geom_point(data = gd_tracer_data, aes(.time*1000, H_val), col ="firebrick4") + 
  geom_point(data = cg_tracer_data, aes(.time*1000, H_val), col = "darkgreen") + ggtitle("Blue: Newton,  Red: GD, Green: CG. ")+ labs(x= "Time [milliseconds]", y = "H value") + scale_y_log10()

ggplot(data = newton_tracer_data, aes(.time*1000, grad_norm)) + geom_point(col = "dodgerblue4")  + geom_smooth(data = gd_tracer_data, aes(.time*1000, grad_norm), col = "firebrick4") + 
  geom_smooth(data = cg_tracer_data, aes(.time*1000, grad_norm), col = "darkgreen") + ggtitle("Blue: Newton,  Red: GD, Green: CG.") + labs(x= "Time [milliseconds]") + scale_y_log10()

############ ONE PLOT #################

p1 <- ggplot(data = newton_tracer_data, aes(.time*1000, H_val)) + geom_point(col = "dodgerblue4") +  geom_smooth(data = gd_tracer_data, aes(.time*1000, H_val), col ="firebrick4") + 
  geom_smooth(data = cg_tracer_data, aes(.time*1000, H_val), col = "darkgreen") + ggtitle("Blue: Newton,  Red: GD, Green: CG. ")+ labs(x= "Time [milliseconds]", y = "H value") + scale_y_log10()

p2 <- ggplot(data = newton_tracer_data, aes(.time*1000, grad_norm)) + geom_point(col = "dodgerblue4")  + geom_smooth(data = gd_tracer_data, aes(.time*1000, grad_norm), col = "firebrick4") + 
  geom_smooth(data = cg_tracer_data, aes(.time*1000, grad_norm), col = "darkgreen") + ggtitle("Blue: Newton,  Red: GD, Green: CG.") + labs(x= "Time [milliseconds]", y = "Gradient of H") + scale_y_log10()

grid.arrange(p1, p2, ncol = 2)
