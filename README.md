
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Hawkes Process

The goal of project is to simulate Hawkes Process. I am considering Hawkes process with exponentially decaying intensity. Mathematically it can be written as 
$$\lambda_t = a + (\lambda_0 -a )e^{-\delta t} + \Sigma_{0 \le T_k <  t}Y_ke^{-\delta (t-T_k)} $$
A function $rHawkes(.)$ takes rate $\lambda$, and decay rate $\delta$ and jump distribution $g(.)$ as the parameters. This function returns random samples from Hawkes process.
