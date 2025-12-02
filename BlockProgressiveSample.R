rm(list=ls(all=TRUE))

n <- 100 #sample slize
N <- c(40, 33, 27) # n1, n2, n3 ...each block sample size
M <- c(35, 30, 25) # s1, s2, s3 ...each block effective sample size
k <- length(M) #Block size

#Initial Value of Parameters
alp <- 1.7; bet <- 1.5;
 

#Defined Censoring Schemes
R <- list()
for(r in 1:k){
	R[[r]] <- c(rep(0,(M[r]-1)),(N[r]-M[r])) # Censoring Schemes
}

# -------------------------
# Weibull quantile function
# F^{-1}(u) = [ -log(1-u) / beta ]^{1/alpha}
# -------------------------
q_weibull <- function(u, alpha, beta) {
  ((-log(1 - u)) / beta)^(1 / alpha)
}

# -------------------------
# Progressive censored sample generator
# qfunc : quantile function F^{-1}(u, alpha, beta)
# alpha, beta : parameters
# N : total sample size (not directly used but included for clarity)
# M : number of observed failures
# R : censoring scheme vector of length M
# -------------------------
sam_generation_progressive <- function(qfunc, alpha, beta, N, M, R) {
  
  if (length(R) != M) {
    stop("Length of R must be equal to M.")
  }
  
  # Step 1: generate uniform variables
  w  <- runif(M)
  v  <- numeric(M)
  u  <- numeric(M)
  
  # Step 2: construct v_i as in the progressive censoring algorithm
  for (i in 1:M) {
    denom   <- i + sum(R[(M - i + 1):M])
    v[i]    <- w[i]^(1 / denom)
  }
  
  # Step 3: construct u_i from v_i
  for (i in 1:M) {
    u[i] <- 1 - prod(v[(M - i + 1):M])
  }
  
  # Step 4: transform uniforms to Weibull via quantile function
  x <- qfunc(u, alpha, beta)
  return(x)
}

# -------------------------
# Generate block progressively censored sample
# -------------------------
set.seed(123)  # for reproducibility

smpl_list <- list()
for(sm in 1:k){
	smpl_list[[sm]] <- sam_generation_progressive(qfunc = q_weibull,
		alpha = alp,
		beta  = bet,
		N     = N[sm],
		M     = M[sm],
		R     = R[[sm]]
	) + (rnorm(M[sm], 0, 0.001))
}
x <- smpl_list; x;
