rm(list=ls(all=T))
rm(list=ls(all=TRUE))
start.time <- Sys.time()
library(nleqslv);
library(NLRoot);

## Sample generation for block progressive censored data from UIW Distribution

n <- 100
N <- c(40, 33, 27) # n1, n2, n3
M <- c(35, 30, 25) # s1, s2, s3
  # M <- c(30, 28, 20) # s1, s2, s3

 
x0=2;
N0=1000; NN=2000;los=0.05
times=10
m <- 3
k <- length(M)
R <- list()
for(r in 1:m){
	R[[r]] <- c(rep(0,(M[r]-1)),(N[r]-M[r])) # I
	#R[[r]] <- c((N[r]-M[r]),rep(0,(M[r]-1))) # II
}
alp <- 1.7; bet <- 1.5;
sam_generation_progressive <- function(nn, mm,al, bt,RR){
	ww <- runif(mm);
	vv <- numeric(0);
	ui <- numeric(0);
	for(i1 in 1:mm){
		vv[i1] <- (ww[i1])^(1/(i1 + sum(RR[(mm-i1+1):mm])));
	}
	for(i2 in 1:mm){
		ui[i2] <- 1-prod(vv[(mm-i2+1):mm]);
	}
	#xs <- al*(-1+((1-ui)^(-1/bt))); #quantile of Lomax distribution
	#xs <- al*(1-((1-ui)^(1/bt))); #quantile of Generalized Pareto distribution
	xs <- ((-1/bt)*log(1-ui))^(1/al); #quantile of Weibull distribution
	#xs <- log(1+(((-1/bt)*log(1-ui))^(1/al))); #quantile of Modified Weibull distribution
	return(xs);
}
xsamp <- sam_generation_progressive(N[1],M[1], alp, bet, R[[1]]);xsamp; # sample for progressive censoring

smpl_list <- list()
for(sm in 1:k){
	smpl_list[[sm]] <- sam_generation_progressive(N[sm],M[sm],alp,bet,R[[sm]]) + (rnorm(M[sm], 0, 0.001))
}
x <- smpl_list; x;
gtt <- list()
for(i in 1:k){
	gtt[[i]] <- x[[i]];#Weibull distribution
	#gtt[[i]] <- exp(x[[i]])-1;#Modified Weibull distribution
}
gtt;

gx <- function(xt){
	#return(exp(xt)-1)#modified weibull
	return(xt)#weibull
}

#....... survival function
sf_func <- function(ap,bt,xt){
	return(exp(-bt*((gx(xt))^ap)));
 }
 
  
 
 #..... Hazard rate function
hf_func <- function(ap,bt,xt){
	#return(ap*bt*((gx(xt))^(ap-1))*exp(xt)); #modified weibull
	return(ap*bt*((gx(xt))^(ap-1))); #weibull
 }
 
 
#...Median time to failure
	med_time <- function(ap, bt){
		return((1/(bt^(1/ap)))*gamma((1/ap)+1)); #Weibull distribution
	}



##............. Bisection Method
		bisection <- function(f, a, b, n = 1000, tol = 1e-7) {
			# If the signs of the function at the evaluated points, a and b, stop the function and return message.
			print(c(f(a),f(b)))
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
		


#...HPD interval using M-H algorithm
	
HPDIntv <- function(btb){
	btt <- sort(btb[N0+1:NN]);
	sim <- los*(NN-N0);
	simm <- (1-los)*(NN-N0);
	CI_btt_lwr=CI_btt_upr=CI_btt_lnth=numeric();
	for(s in 1:sim){
		CI_btt_lwr[s] <- btt[s];
		CI_btt_upr[s] <- btt[floor(s+simm)];
		CI_btt_lnth[s] <- CI_btt_upr[s]-CI_btt_lwr[s];
	}
	hpd_intv <- function(lnt_cc, cc_lwr, cc_upr){ ### find the lower and upper bound of the hpd interval
		BB_index <- min(lnt_cc); 
		BB_position <- which(c(lnt_cc)== BB_index);
		BB_int_lwr <- cc_lwr[BB_position[1]];
		BB_int_upr <- cc_upr[BB_position[1]];
		return(c(BB_int_lwr, BB_int_upr));	
	}		
	hpd_CI_btt <- hpd_intv(CI_btt_lnth,CI_btt_lwr,CI_btt_upr);  # hpd credible interval of parameter 
	lnt_hpd_lmd <- hpd_CI_btt[2]-hpd_CI_btt[1]; # hpd interval length of parameter
	return(c(hpd_CI_btt[1],hpd_CI_btt[2],lnt_hpd_lmd));
}
	
##...pivotal based inference

pivotap <- function(aa, bb, rr, gt, mm, nn, xx0){
	
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
					#print(smj)
			}
			# smjj[i] <- log((nn[i]*((gt[[i]][1])^ap))/demit); # first term W_i1
			# # print(smj)
			smi <- smi + smj;	
		}
		#print(smjj)
		# piv_alp <- -2*(sum(smjj)+smi);
		piv_alp <- -2*(smi);
		pp <- piv_alp-apstr;
		
		print(c(piv_alp,apstr))
		return(pp);
	}
	#pppval <- bisection(pivalp, a=aa, b=bb);
	
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
		pvt_ap_val <- bisection(pivalp, a=aa, b=bb); #pivotap(aa=0.2, bb=2, rr=R, gt=x, mm=M, nn=N );
		#pvt_ap_val <- nleqslv(1.5,pivalp)
		#print(pvt_ap_val$message);
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
	
		pvt_sf[ii] <- sf_func(pvt_ap[ii],pvt_bt[ii],xx0); # pivotal SF
		pvt_hf[ii] <- hf_func(pvt_ap[ii],pvt_bt[ii],xx0); # pivotal HF
		pvt_mdtf[ii] <- med_time(pvt_ap[ii],pvt_bt[ii]);  #pivotal MDTF
	}
	# print(pvt_ap[(N0+1):NN])
	# print(sort(pvt_ap[(N0+1):NN]))
	#pivotal estimates
	pvt_ap_est <- mean(pvt_ap[(N0+1):NN]); #pivt alp
	pvt_bti_est1 <- mean(pvt_bti[,1][(N0+1):NN]);#pivt bet1
	pvt_bti_est2 <- mean(pvt_bti[,2][(N0+1):NN]);#pivt bet2
	pvt_bti_est3 <- mean(pvt_bti[,3][(N0+1):NN]);#pivt bet3
	pvt_bt_est <- mean(pvt_bt[(N0+1):NN]); #pivt bet
	pvt_sf_est <- mean(pvt_sf[(N0+1):NN]); #pivt sf
	pvt_hf_est <- mean(pvt_hf[(N0+1):NN]); #pivt hf
	pvt_mdtf_est <- mean(pvt_mdtf[(N0+1):NN]); #pivt mdtf

	#pivotal variance
	pvt_ap_vr <- var(pvt_ap[(N0+1):NN]);#variance alp
	pvt_bti_vr1 <- var(pvt_bti[,1][(N0+1):NN]);#variance bet1
	pvt_bti_vr2 <- var(pvt_bti[,2][(N0+1):NN]); #var beta2
	pvt_bti_vr3 <- var(pvt_bti[,3][(N0+1):NN]);#varbet3
	pvt_bt_vr <- var(pvt_bt[(N0+1):NN]);#var bet
	pvt_sf_vr <- var(pvt_sf[(N0+1):NN]);#var sf
	pvt_hf_vr <- var(pvt_hf[(N0+1):NN]);#var hf
	pvt_mdtf_vr <- var(pvt_mdtf[(N0+1):NN]);#var mdtf

	#Generalized confidence interval;
	hpd_ap <- HPDIntv(pvt_ap);#GIC alp
	hpd_bti1 <- HPDIntv(pvt_bti[,1]);#GIC bet1
	hpd_bti2 <- HPDIntv(pvt_bti[,2]);#GIC bet2
	hpd_bti3 <- HPDIntv(pvt_bti[,3]);#GIC bet3
	hpd_bt <- HPDIntv(pvt_bt);#GIC bet
	hpd_sf <- HPDIntv(pvt_sf);#GIC sf
	hpd_hf <- HPDIntv(pvt_hf);#GIC hf
	hpd_mdtf <- HPDIntv(pvt_mdtf);#GIC mdtf
	
	# estimates
	piv_est <- matrix(0,nrow=k+5, ncol=5);
	piv_est[1,1] <- pvt_ap_est; piv_est[1,2] <- pvt_ap_vr; piv_est[1,3:5] <- hpd_ap; #alpha
	piv_est[2,1] <- pvt_bti_est1; piv_est[2,2] <- pvt_bti_vr1; piv_est[2,3:5] <- hpd_bti1;#beta_1
	piv_est[3,1] <- pvt_bti_est2; piv_est[3,2] <- pvt_bti_vr2; piv_est[3,3:5] <- hpd_bti2;#beta_2
	piv_est[4,1] <- pvt_bti_est3; piv_est[4,2] <- pvt_bti_vr3; piv_est[4,3:5] <- hpd_bti3;#beta_3
	piv_est[5,1] <- pvt_bt_est; piv_est[5,2] <- pvt_bt_vr; piv_est[5,3:5] <- hpd_bt;#beta
	piv_est[6,1] <- pvt_sf_est; piv_est[6,2] <- pvt_sf_vr; piv_est[6,3:5] <- hpd_sf;#SF
	piv_est[7,1] <- pvt_hf_est; piv_est[7,2] <- pvt_hf_vr; piv_est[7,3:5] <- hpd_hf;#HF
	piv_est[8,1] <- pvt_mdtf_est; piv_est[8,2] <- pvt_mdtf_vr; piv_est[8,3:5] <- hpd_mdtf;#MDTF
	return(piv_est);
}
pvt_val <- pivotap(aa=1, bb=2, rr=R, gt=x, mm=M, nn=N, xx0=x0);
pvt_val;

