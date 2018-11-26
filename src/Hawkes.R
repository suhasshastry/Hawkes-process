#Poisson process
lambda <- 1
time <- 5
N <- qpois(1-1e-1, lambda = lambda * time)
X <- rexp(n = N, rate = lambda)
S <- c(0, cumsum(X))
plot(x = S, y = 0:N, type = "s", xlim = c(0, time),xlab = 't',ylab = 'N(t)',
     main="Homogeneous Poisson process with rate=1")
grid()
abline(v = S,col="red",lty=2)


#Homogeneous Poisson Process paths
library(tidyverse)
library(ggplot2)
lambda <- 1
time <- 20
N <- qpois(1-1e-1, lambda = lambda * time)
X <- rerun(100,rexp(n = N, rate = lambda))
S <- X %>% map(~cumsum(.))
df <- matrix(flatten_dbl(S),byrow = F,ncol = 100)
matplot(x = df, y = 1:N, type = "s",col='black',
        xlim = c(0, time),
        main = "Homogeneous Poisson Process paths", xlab = "t", ylab = "N(t)")

#Homgenoous Poisson Process - Method 2
simulate_homogenous_poisson <- function(lambda,Time){
  n <-1
  t <- array()
  t[n] <- 0
  while(TRUE){
    u <- runif(1)
    w <- -1/lambda*log(u)
    t[n+1] <- t[n] + w
    if(t[n+1] > Time){
      return(t) #return t 
    }else{
      n <- n + 1
    }
  }
}

plot(simulate_homogenous_poisson(1,10),type="s")

#Inhomgoneous Poisson Process
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
plot(t,0:(N-1),type = 's',xlab='t',ylab='N(t)',
     main="Inhomogoneous Poisson Process with rate=t/2")
abline(v = t,col="red",lty=2)

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
hawkes_data <- hawkes(1,2,3,10)
remove_cum_sum <- function(x){
  temp <- array()
  temp[1] <- x[1]
  for(i in 2:length(x)){
    temp[i] = x[i]-x[i-1]
  }
  temp
}
library(PtProcess)
mpp(remove_cum_sum(hawkes(1,160,201,20)),simple_gif,marks=list(NULL, NULL),params = p <- c(0.02, 70.77, 0.47, 0.002, 1.25))


# nro<-10
# S<-vector(mode="integer",length = nro)
# S[1]=0

##Generation of Exponential random variables with parameter lambda
S <- round(hawkes(1,2,3,10),2)

#Plot of the trajectory and add lines in the arrival times
n_func <- function(t, S) sapply(t, function(t) sum(S <= t))
t_series <- seq(0, max(S), by = max(S)/100)
S <- round(hawkes(1,2,3,10),2)
plot(t_series, n_func(t_series, S),type = "s",
ylab='N(t)',xlab="t",main="Hawkes Process, lambda = 1, alpha = 2, beta = 3")
grid()
abline(v = S,col="red",lty=2)


