library(data.table)
library(pracma)

library(matrixStats) #logSumExp
library(genieclust) #adjusted_rand_score

library(ggplot2)
library(ggpubr)


X_df = as.matrix(fread("mixture1.geno"))
Xnew_df = as.matrix(fread("mixture2.geno"))
F_df = as.matrix(fread("mixture1.freq"))
Z_df = as.matrix(fread("mixture1.ganc"))

#X: N by M data matrix
#gamma: P(Z_i=k | X_i=x_i, theta) N by K
M_step = function(X, gamma){
	N = dim(X)[1]; M = dim(X)[2]
	K = dim(gamma)[2]
	
	######### TODO 3a: modify the following to have meaning updates #########
	pis <- colMeans(gamma)
	
	# Frequency matrix F: M by K
	F <- t(X) %*% gamma
	
	# Avoid divide-by-zero error by setting frequency terms associated with zero assignments to 0
	F_colsums <- colSums(gamma)
	F_colsums[F_colsums == 0] <- 1  # Avoid divide-by-zero error
	F <- F / F_colsums
	
	######### end of modification #########
	return(list("pis" = pis, "F" = F))
}

#X: N by M data matrix
#params: a list with two parameters returned from M_step
E_step = function (X, params, thr = 10**(-8)){
	F = params$"F" #Frequency matrix F: M by K
	pis = params$"pis" #Proportion vector pi: length K
	N = dim(X)[1]; M = dim(X)[2]
	K = dim(F)[2] 
	######### TODO 3b: modify the following to have meaning updates #########
	#calculate weighted_log_prob: log(P(X_i=x_i | Z_i=k, theta) * P(Z_i=k | theta))
	#calculate log_prob_sample: log P(Xi=x_i | theta) N by 1. Hint: use logSumExp function
	#calculate log_prob_data: logP(X_1:n=x_1:n | theta) scalar
	#calculate log_gammas: log P(Z_i=k | X_i=x_i, theta) N by K
	log_likelihoods = matrix(0, nrow = N, ncol = K)
	for (k in 1:K) {
	  log_likelihoods[, k] = (X * log(F[, k]) + (1 - X) * log(pmax(1 - F[, k], thr)))
	}
	
	weighted_log_prob = matrix(0, nrow = N, ncol = K)
	for (i in 1:N) {
	  for (k in 1:K) {
	    weighted_log_prob[i, k] = log_likelihoods[i, k] + log(pis[k])
	  }
	}
	
	log_prob_sample = logSumExp(weighted_log_prob)
	log_prob_data = sum(logSumExp(log_likelihoods + log(pis)))
	
	log_gammas = log_likelihoods + log(pis) - matrix(logSumExp(log_likelihoods + log(pis), axis = 2), N, K, byrow = TRUE)
	
	######### end of modification #########
	return(list("log_gammas" = log_gammas, "log_prob_data" = log_prob_data))
}

EM = function(X, K = 2, max_iter = 100, tol = 10**(-4), n_init = 3, debug = FALSE){
	N = dim(X)[1];  M = dim(X)[2]
	res = list()
	best_log_prob_data = -Inf
	converged = FALSE
	
	#loop through different random starting points
	for (init in 1:n_init){
		set.seed(init)
		if(debug){print(paste0("starting EM on random initialization: ", init, " out of ", n_init))}
			
		######### TODO 3c: modify the following to have the full EM updates #########

		#initialize soft assignment 
		gammas = NULL
		
		
		for (n_iter in 1:max_iter){
			prev_log_prob_data = log_prob_data
			
			
			
			
			
		
			######### convergence check #########
			change = (log_prob_data - prev_log_prob_data)/N
			if (abs(change) < tol){
				if(debug){
					print(paste0("random initialization ", init ," converged at iteration ", n_iter))
					print("")
				}
				converged = TRUE
				break
				}
			
		
		######### update on the best initialization #########
		best_init = NULL
		
		
		
		}
	}
	
	
	######### end of modification #########
	res[["converged"]] = converged     
	res[["best_init"]] = best_init
	return(res)
}


#Question 1
GMM<-M_step(X_df, Z_df)
#Question 2
E<-E_step(X_df,GMM)
adjusted_rand_score()

