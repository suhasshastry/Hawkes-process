# #Poisson process
# par(mfrow=c(1,1))
# lambda <- 1
# time <- 5
# N <- qpois(1-1e-1, lambda = lambda * time)
# X <- rexp(n = N, rate = lambda)
# S <- c(0, cumsum(X))
# plot(x = S, y = 0:N, type = "s", xlim = c(0, time),xlab = 't',ylab = 'N(t)',
#      main="Homogeneous Poisson process with rate=1")
# grid()
# abline(v = S,col="red",lty=2)
# 
# 
# #Homogeneous Poisson Process paths
# library(tidyverse)
# library(ggplot2)
# lambda <- 1
# time <- 20
# N <- qpois(1-1e-1, lambda = lambda * time)
# X <- rerun(100,rexp(n = N, rate = lambda))
# S <- X %>% map(~cumsum(.))
# df <- matrix(flatten_dbl(S),byrow = F,ncol = 100)
# matplot(x = df, y = 1:N, type = "s",col='black',
#         xlim = c(0, time),
#         main = "Homogeneous Poisson Process paths", xlab = "t", ylab = "N(t)")
# 
# #Homgenoous Poisson Process - Method 2
# simulate_homogenous_poisson <- function(lambda,Time){
#   n <-1
#   t <- array()
#   t[n] <- 0
#   while(TRUE){
#     u <- runif(1)
#     w <- -1/lambda*log(u)
#     t[n+1] <- t[n] + w
#     if(t[n+1] > Time){
#       return(t) #return t 
#     }else{
#       n <- n + 1
#     }
#   }
# }
# 
# plot(simulate_homogenous_poisson(1,10),type="s")
# 
# #Inhomgoneous Poisson Process
# simulate_inhomogenous_poisson <- function(lambda,Time=10){
#   n <- 1
#   m <- 1
#   t <- array()
#   t[n] <- 0
#   s <- array()
#   s[m] <- 0
#   t_seq <- seq(from = 0.01,to = Time,by = 0.01)
#   lambda_bar <- max(lambda(t_seq))
#   while(s[m] < Time){
#     u <- runif(1)
#     w <- -1/lambda_bar*log(u)
#     s[m+1] <- s[m] + w
#     d <- runif(1)
#     if(d<= lambda(round(s[m+1],2))/lambda_bar){
#       t[n+1] <- s[m+1]
#       n <- n+1
#     }
#     m <- m + 1
#   }
#   t
# }
# 
# Time <- 10
# 
# lambda <- function(t){
#   0.5*t
# }
# t <- simulate_inhomogenous_poisson(lambda,Time)
# N <- length(t)
# plot(t,0:(N-1),type = 's',xlab='t',ylab='N(t)',
#      main="Inhomogoneous Poisson Process with rate=t/2")
# abline(v = t,col="red",lty=2)
# 
# #Hawkes Process -thinning procedure
# hawkes <- function(mu,alpha,beta,Time){
#   s <- 0
#   n <- 1
#   t <- array()
#   t[n] <- 0
#   while(s<Time){
#     lambda_bar <- mu + sum(alpha*exp(-beta*(s-t)))
#     u <- runif(1)
#     w <- -1/lambda_bar*log(u)
#     s <- s + w
#     d <- runif(1)
#     if(d*lambda_bar < mu + sum(alpha*exp(-beta*(s-t)))){
#       t[n] <- s
#       n <- n + 1
#     }
#   }
#   t
# }
# hawkes_data <- hawkes(1,2,3,10)
# remove_cum_sum <- function(x){
#   temp <- array()
#   temp[1] <- x[1]
#   for(i in 2:length(x)){
#     temp[i] = x[i]-x[i-1]
#   }
#   temp
# }
# library(PtProcess)
# #mpp(remove_cum_sum(hawkes(1,160,201,20)),simple_gif,marks=list(NULL, NULL),params = p <- c(0.02, 70.77, 0.47, 0.002, 1.25))
# 
# 
# # nro<-10
# # S<-vector(mode="integer",length = nro)
# # S[1]=0
# 
# ##Generation of Exponential random variables with parameter lambda
# S <- round(hawkes(1,2,3,10),2)
# 
# #Plot of the trajectory and add lines in the arrival times
# n_func <- function(t, S) sapply(t, function(t) sum(S <= t))
# t_series <- seq(0, max(S), by = max(S)/100)
# S <- round(hawkes(1,2,3,10),2)
# plot(t_series, n_func(t_series, S),type = "s",
# ylab='N(t)',xlab="t",main="Hawkes Process, lambda = 1, alpha = 2, beta = 3")
# grid()
# abline(v = S,col="red",lty=2)
# 
# 
# loglikelihood_1 <- function(t,init_params){
#   lambda <- init_params[1]
#   alpha <- init_params[2]
#   beta <- init_params[3]
#   t_i <- tail(t,1)
#   #print(t)
#   #print(init_params)
#   temp<- log(lambda + alpha*(sum(exp(-beta*(t_i-t)))-1))
#   #print(paste(temp,'\n\n'))
#   temp
# }
# 
# loglikelihood_2 <- function(t,init_params){
#   lambda <- init_params[1]
#   alpha <- init_params[2]
#   beta <- init_params[3]
#   t_k <- tail(t,1)
#   -lambda*t_k + alpha/beta*sum(exp(-beta*(t_k-t))-1)
# }
# 
# loglikelihood_11 <- function(t,init_params,A){
#   lambda <- init_params[1]
#   alpha <- init_params[2]
#   beta <- init_params[3]
#   t_i <- tail(t,1)
#   A <- exp(-beta*(t_i-t_i-1))
# }
# 
# hawkes_ll <- function(init_params,t){
#   n <- length(t)
#   likelihood1 <- sum(map(1:n,~loglikelihood_1(t[1:.x],init_params))%>% flatten_dbl())
#   likelihood2 <- sum(map(1:n,~loglikelihood_2(t[1:.x],init_params)) %>% flatten_dbl())
#   -1*(likelihood1 + likelihood2)
# }
# t <- hawkes(10,20,30,10)
# (mle1 <- optim(par = c(0.29168892, 0.04131991, 0.79237880), 
#                fn = hawkes_ll, t = data, method = "L-BFGS-B",lower=c(0,0,0)))
# 
# 
# 
# 
# neg.loglik <- function(params, data, opt=TRUE) {
#   mu <- params[1]
#   alpha <- params[2]
#   beta <- params[3]
#   t <- sort(data)
#   r <- rep(0,length(t))
#   for(i in 2:length(t)) {
#     r[i] <- exp(-beta*(t[i]-t[i-1]))*(1+r[i-1])
#   }
#   loglik <- -tail(t,1)*mu
#   loglik <- loglik+alpha/beta*sum(exp(-beta*(tail(t,1)-t))-1)
#   loglik <- loglik+sum(log(mu+alpha*r))
#   if(!opt) {
#     return(list(negloglik=-loglik, mu=mu, alpha=alpha, beta=beta, t=t,
#                 r=r))
#   }
#   else {
#     return(-loglik)
#   }
# }
# 
# chic <- read.csv("/cloud/project/data/chic.csv", sep="")
# data <- chic$x
# S <- data
# t_series <- seq(0, max(S), by = max(S)/100)
# plot(t_series, n_func(t_series, S),type = "s",
#      ylab='N(t)',xlab="t",main="Chicago Burglary data in Beat 423 from 2017 to present")
# grid()
# abline(v = S,col="red",lty=2)
# 
# 
# S <- hawkes(0.29168892, 0.04131991, 0.79237880,700)
# par(mfrow=c(2,1))
# t_series <- seq(0, max(S), by = max(S)/100)
# plot(t_series, n_func(t_series, S),type = "s",
#      ylab='N(t)',xlab="t",main="Simulated Hawkes data with lambda = 0.29, alpha = 0.04,
#      beta = 0.79")
# grid()
# abline(v = S,col="red",lty=2)
# hist(S,100,xlab='t',ylab='N(t)')
# 
# # P values histogram
# S <- rerun(1000,hawkes(0.29168892, 0.04131991, 0.79237880,700))
# p_value <- map(S,~ks.test(data,.x)$p.value) %>% flatten_dbl()
# hist(p_value,20,main="Histogram of p value",xlab="p value")
# abline(v=0.05,col="red",lty=2)
# 
# 
# optim(par=c(2,1,2), fn=neg.loglik, data=data,method = "L-BFGS-B",lower = 0)


library(tidyverse)
library(ggplot2)

#1. Homogenous Poisson Process
lambda <- 1
time <- 5
N <- qpois(1-1e-1, lambda = lambda * time)
X <- rexp(n = N, rate = lambda)
S <- c(0, cumsum(X))
# plot(x = S, y = 0:N, type = "s", xlim = c(0, time),xlab = 't',ylab = 'N(t)',
#      main="Homogeneous Poisson process with rate=1")
# grid()
# abline(v = S,col="red",lty=2)
df <- data.frame(s=S)
ggplot(data=df) + 
  geom_step(mapping = aes(x=s,y=0:N)) +
  coord_cartesian(xlim = c(0, time)) +
  theme_bw() +
  geom_vline(xintercept = df$s,linetype="dashed", color = "red") +
  labs(x='t',y=expression(N[t]))


#2. Inhomgoneous Poisson Process
simulate_inhomogenous_poisson <- function(lambda,Time=10){
  n <- 1
  m <- 1
  t <- array()
  t[n] <- 0
  s <- array()
  s[m] <- 0
  t_seq <- seq(from = 0.01,to = Time,by = 0.01)
  lambda_bar <- max(lambda(t_seq))
  while(s[m] < Time){
    u <- runif(1)
    w <- -1/lambda_bar*log(u)
    s[m+1] <- s[m] + w
    d <- runif(1)
    if(d<= lambda(round(s[m+1],2))/lambda_bar){
      t[n+1] <- s[m+1]
      n <- n+1
    }
    m <- m + 1
  }
  t
}

Time <- 10

lambda <- function(t){
  0.5*t
}
t <- simulate_inhomogenous_poisson(lambda,Time)
N <- length(t)
# plot(t,0:(N-1),type = 's',xlab='t',ylab='N(t)',
#      main="Inhomogoneous Poisson Process with rate=t/2")
# abline(v = t,col="red",lty=2)
df <- data.frame(s=t)
ggplot(data=df) + 
  geom_step(mapping = aes(x=s,y=0:(N-1))) +
  theme_bw() +
  geom_vline(xintercept = df$s,linetype="dashed", color = "red") +
  labs(x='t',y=expression(N[t]))

#Hawkes Process
hawkes <- function(mu,alpha,beta,Time){
  s <- 0
  n <- 1
  t <- array()
  t[n] <- 0
  while(s<Time){
    lambda_bar <- mu + sum(alpha*exp(-beta*(s-t)))
    u <- runif(1)
    w <- -1/lambda_bar*log(u)
    s <- s + w
    d <- runif(1)
    if(d*lambda_bar < mu + sum(alpha*exp(-beta*(s-t)))){
      t[n] <- s
      n <- n + 1
    }
  }
  t
}
S <- round(hawkes(1,2,3,10),2)

n_func <- function(t, S) sapply(t, function(t) sum(S <= t))
t_series <- seq(0, max(S), by = max(S)/100)
S <- round(hawkes(1,2,3,10),2)
# plot(t_series, n_func(t_series, S),type = "s",
#      ylab='N(t)',xlab="t",main="Hawkes Process, lambda = 1, alpha = 2, beta = 3")
# grid()
# abline(v = S,col="red",lty=2)
df <- data.frame(x=t_series,y=n_func(t_series, S))
ggplot(data=df) + 
  geom_step(mapping = aes(x=x,y=y)) +
  theme_bw() +
  geom_vline(xintercept = S,linetype="dashed", color = "red") +
  labs(x='t',y=expression(N[t]))

# MLE optimization

loglikelihood_1 <- function(t,init_params){
  lambda <- init_params[1]
  alpha <- init_params[2]
  beta <- init_params[3]
  t_i <- tail(t,1)
  #print(t)
  #print(init_params)
  temp<- log(lambda + alpha*(sum(exp(-beta*(t_i-t)))-1))
  #print(paste(temp,'\n\n'))
  temp
}

loglikelihood_2 <- function(t,init_params){
  lambda <- init_params[1]
  alpha <- init_params[2]
  beta <- init_params[3]
  t_k <- tail(t,1)
  -lambda*t_k + alpha/beta*sum(exp(-beta*(t_k-t))-1)
}

loglikelihood_11 <- function(t,init_params){
  lambda <- init_params[1]
  alpha <- init_params[2]
  beta <- init_params[3]
  A <- rep(0,length(t))
  for(i in 2:length(t)) {
    A[i] <- exp(-beta*(t[i]-t[i-1]))*(1+A[i-1])
  }
  sum(log(lambda+alpha*A))
}

hawkes_ll <- function(init_params,t){
  n <- length(t)
  #likelihood1 <- sum(map(1:n,~loglikelihood_11(t[1:.x],init_params))%>% flatten_dbl())
  likelihood1 <- loglikelihood_11(t,init_params)
  #likelihood2 <- sum(map(1:n,~loglikelihood_2(t[1:.x],init_params)) %>% flatten_dbl())
  likelihood2 <- loglikelihood_2(t,init_params)
  -1*(likelihood1 + likelihood2)
}
t <- hawkes(10,20,30,10)

(mle1 <- optim(par = c(9, 19, 29), 
               fn = hawkes_ll, t=t, method = "L-BFGS-B",lower=c(0,0,0)))


#Chicago crime data
chic <- read.csv("/cloud/project/data/chic.csv", sep="")
data <- chic$x
S <- data
t_series <- seq(0, max(S), by = max(S)/100)
# plot(t_series, n_func(t_series, S),type = "s",
#      ylab='N(t)',xlab="t",main="Chicago Burglary data in Beat 423 from 2017 to present")
# grid()
# abline(v = S,col="red",lty=2)
# hist(S,100)

df <- data.frame(x=t_series,y=n_func(t_series, S))
plot1 <- ggplot(data=df) +
  geom_step(mapping = aes(x=x,y=y)) +
  theme_bw() +
  coord_cartesian(xlim = c(0, 700)) +
  geom_vline(xintercept = S,linetype="dashed", color = "red") +
  labs(x='t',y=expression(N[t])) +
  ggtitle('Chicago Burglary data in Beat 423 from 2017 to present')

plot2 <- ggplot() + 
  geom_histogram(mapping=aes(x=S),bins = 100) +
  theme_bw() +
  coord_cartesian(xlim = c(0, 700)) +
  labs(x='t',y=expression(N[t])) +
  ggtitle('Chicago Burglary data in Beat 423 from 2017 to present')
grid.arrange(plot1, plot2, nrow=2)

#Optimization of MLE
(mle1 <- optim(par = c(.1,2,3), fn = hawkes_ll, t = data, method = "L-BFGS-B",lower=c(0,0,0))$par)

#Simulated hawkes process with params lm = 0.29,alpha = 0.04, beta=0.79
S <- hawkes(0.29168892, 0.04131991, 0.79237880,700)
df <- data.frame(x=t_series,y=n_func(t_series, S))
plot1 <- ggplot(data=df) +
  geom_step(mapping = aes(x=x,y=y)) +
  theme_bw() +
  coord_cartesian(xlim = c(0, 700)) +
  geom_vline(xintercept = S,linetype="dashed", color = "red") +
  labs(x='t',y=expression(N[t])) +
  ggtitle('Simulated Hawkes data with lambda = 0.29, alpha = 0.04,
     beta = 0.79')

plot2 <- ggplot() + 
  geom_histogram(mapping=aes(x=S),bins = 100) +
  theme_bw() +
  coord_cartesian(xlim = c(0, 700)) +
  labs(x='t',y=expression(N[t])) +
  ggtitle('Simulated Hawkes data with lambda = 0.29, alpha = 0.04,
     beta = 0.79')
library(gridExtra)
grid.arrange(plot1, plot2, nrow=2)


# P values histogram
set.seed(1)
S <- rerun(1000,hawkes(0.29168892, 0.04131991, 0.79237880,700))
p_value <- map(S,~ks.test(data,.x)$p.value) %>% flatten_dbl()
# hist(p_value,20,main="Histogram of p value",xlab="p value")
# abline(v=0.05,col="red",lty=2)
ggplot() +
  geom_histogram(mapping = aes(x=p_value),bins = 21) +
  theme_bw() +
  geom_vline(xintercept = 0.05,linetype="dashed", color = "red")
mean(p_value<0.05)
