#Poisson process
lambda <- 1
time <- 5
N <- qpois(1-1e-1, lambda = lambda * time)
X <- rexp(n = N, rate = lambda)
S <- c(0, cumsum(X))
plot(x = S, y = 0:N, type = "s", xlim = c(0, time),xlab = 't',ylab = 'N(t)',main="Poisson process with rate=1") 

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