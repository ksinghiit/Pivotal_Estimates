rm(list=ls(all=T))
rm(list=ls(all=TRUE))

n <- 100 #sample slize
N <- c(40, 33, 27) # n1, n2, n3 ...each block sample size
M <- c(35, 30, 25) # s1, s2, s3 ...each block effective sample size
k <- length(M) #Block size

#Initial Value of Parameters
alp <- 1.7; bet <- 1.5;


NN=2000; #Total generate the pivotal estimates 
N0=1000; #Consider the burn in obtained initial pivotal estimates 

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
ProgressiveSample <- function(qfunc, alpha, beta, N, M, R) {
  
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

smpl_list <- list()
for(sm in 1:k){
  smpl_list[[sm]] <- ProgressiveSample(
    qfunc = q_weibull,
    alpha = alp,
    beta  = bet,
    N     = N[sm],
    M     = M[sm],
    R     = R[[sm]]
  ) + (rnorm(M[sm], 0, 0.001))
}
x <- smpl_list; x;


gtt <- list()
for(i in 1:k){
  gtt[[i]] <- x[[i]];#Weibull distribution
}
gtt;

gx <- function(xt){
  return(xt)#weibull
}

##............. Bisection Method
bisection <- function(f, a, b, n = 1000, tol = 1e-7) {
  # If the signs of the function at the evaluated points, a and b, stop the function and return message.
  #print(c(f(a),f(b)))
  if (!(f(a) < 0) && (f(b) > 0)) {
    stop('signs of f(a) and f(b) differ')
  } else if ((f(a) > 0) && (f(b) < 0)) {
    stop('signs of f(a) and f(b) differ')
  }
  
  for (i in 1:n) {
    c <- (a + b) / 2 # Calculate midpoint
    
    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the 
    # function and return the root.
    if ((f(c) == 0) || ((b - a) / 2) < tol) {
      return(c)
    }
    
    # If another iteration is required, 
    # check the signs of the function at the points c and a and reassign
    # a or b accordingly as the midpoint to be used in the next iteration.
    ifelse(sign(f(c)) == sign(f(a)), 
           a <- c,
           b <- c)
  }
  # If the max number of iterations is reached and no root has been found, 
  # return message and end function.
  print('Too many iterations')
}


##...pivotal based inference

PivotalEstimates <- function(aa, bb, rr, gt, mm, nn){
  
  pivalp <- function(ap){ # pivotal function for alpha
    apstr <- mean(rchisq(5,(2*sum(mm+1))));
    smi <- 0;
    for(i in 1:length(mm)){
      sm21 <- sm22 <- 0;
      for(l in 1:(mm[i]-1)){
        sm21 <- sm21 + ((rr[[i]][l])+1)*((gt[[i]][l])^ap);
        sm22 <- sm22 + ((rr[[i]][l])+1)
      }
      demit <- sm21 + ((nn[i]-sm22)*((gt[[i]][(mm[i])])^ap));#W_imi
      smj <- 0; 
      for(j in 1:(mm[i]-1)){
        sm11 <- sm12 <- 0;
        for(l in 1:(j)){
          sm11 <- sm11 + ((rr[[i]][l])+1)*((gt[[i]][l])^ap);
          sm12 <- sm12 +	((rr[[i]][l])+1)
        }
        numrt <- sm11 + ((nn[i]-sm12)*((gt[[i]][j])^ap)); #W_ij, j=2,3,...,(mm[i]-1)
        smj <- smj + log(numrt/demit);
      }
      smi <- smi + smj;	
    }
    piv_alp <- -2*(smi);
    pp <- piv_alp-apstr;
    return(pp);
  }
  
  #..pivotal_bet_i
  piv_bti <- function(ap, rr, gt, mm, nn){
    pivotbti <-  numeric();
    for(i in 1:length(mm)){
      btstr <- rchisq(1,(2*mm[i]));
      sm21 <- sm22 <- 0;
      for(l in 1:(mm[i]-1)){
        sm21 <- sm21 + ((rr[[i]][l])+1)*((gt[[i]][l])^ap);
        sm22 <- sm22 + ((rr[[i]][l])+1)
      }
      GTi <- 2*(sm21 + ((nn[i]-sm22)*((gt[[i]][mm[i]])^ap)));
      pivotbti[i] <- btstr/GTi;
    }
    return(pivotbti);
  }
  #iteration
  pvt_ap <- pvt_bt <- pvt_sf <- pvt_hf <- pvt_mdtf <- numeric();
  pvt_bti <- matrix(0, nrow=NN, ncol=length(M));
  for(ii in 1:NN){
    pvt_ap_val <- bisection(pivalp, a=aa, b=bb); 
    pvt_ap[ii] <- pvt_ap_val; #pivotal alpha
    pvt_bti_val <- piv_bti(pvt_ap[ii], rr, gt, mm, nn);
    pvt_bti[ii,] <- pvt_bti_val;#pivotal beta_i
    varbti <- etinv <- numeric(); etbti <- 0;
    for(i in 1:length(M)){
      varbti[i] <- var(pvt_bti[,i])
      etinv[i] <- 1/varbti[i] ;
      etbti <- etbti + (etinv[i]*(pvt_bti[ii,i]));
    }
    pvt_bt[ii] <- sum(etbti)/sum(etinv); #pivotal beta
  }
  #pivotal estimates
  pvt_ap_est <- mean(pvt_ap[(N0+1):NN]); #pivt alp
  pvt_bti_est <- numeric()
  for(i in 1:k){
    pvt_bti_est[i] <- mean(pvt_bti[,i][(N0+1):NN]);#pivt beti
  }
  pvt_bt_est <- mean(pvt_bt[(N0+1):NN]); #pivt bet
  
  pivotal_est <- list()
  pivotal_est[[1]] <- pvt_ap_est;
  pivotal_est[[2]] <- pvt_bti_est;
  pivotal_est[[3]] <- pvt_bt_est;
  
  return(pivotal_est)
}
pvt_val <- PivotalEstimates(aa=1, bb=2, rr=R, gt=x, mm=M, nn=N);
pvt_val

