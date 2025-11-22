library(tidyverse)
library(mgcv)
library(splines)
library(coneproj)

# Multinomial TMLE for m strains (J=0 baseline, J=1,S=s for s=1..m)
# Requires mgcv for smoothing alpha and optionally for p_s(t) estimation
if(!requireNamespace("mgcv", quietly=TRUE)) install.packages("mgcv")
library(mgcv)
# tmle_multinomial.R
library(splines)

safe_clamp <- function(x, lower = 1E-10, upper = 1-1E-10) {
  if (any(lower > upper)) stop("lower must be ≤ upper")
  out <- x
  out[!is.na(x) & x < lower] <- lower
  out[!is.na(x) & x > upper] <- upper
  out
}

tmle_multinomial <- function(dat,
                             T.name="T",
                             J.name="J",
                             V.name = "V",
                             V.early.name="V_early",
                             S.name = "S",
                             psi_delta,
                             m = 2,
                             verbose=TRUE,
                             init_beta = NULL,
                             p_s_fun = predict_GISAID,
                             max_iter = 100,
                             tol = 1e-6,
                             kernel_shape = "epanechikov",
                             bandwidth = NULL,
                             automatic_bw = TRUE, ...) {
  
  ##########
  ### helper: compute Q_{0s}(T_i,V_i) given current beta, alpha and mixture p_s(t)
  # Represent parameters: beta (stacked), alpha_link(T) vector length n, p_s_mixture_fun returns m-vector for t.
  compute_probs <- function(Z_stack, beta_vec, alpha_link_vec, p_s_mixture_vecs) {
    # returns n x (m+1) matrix of probabilities: columns 1..m for strains, last column for baseline J=0
    probs <- matrix(0, nrow = n, ncol = m + 1)
    # compute linear predictors for strain s: for each i and s: lp_is = X_s_i^T beta_s + alpha(T_i)
    # Z_stack row has blocks; extract block dot products
    eta_mat <- matrix(NA, nrow = n, ncol = m)
    for (s in 1:m) {
      start <- (s - 1) * d + 1
      bs <- beta_vec[start:(start + d - 1)]
      Xs <- Z_stack[, start:(start + d - 1), drop = FALSE]
      eta_mat[, s] <- as.numeric(as.numeric(Xs %*% bs) + alpha_link_vec)
    }
    # Now compute numerators: p_s(t) * exp{I(V<=T) f_s + alpha0}
    # We assume p_s_mixture_vecs is n x m matrix giving p_s(t_i)
    numer <- matrix(0, nrow = n, ncol = m)
    for (s in 1:m) {
      numer[, s] <- p_s_mixture_vecs[, s] * exp(eta_mat[, s])
    }
    denom <- 1 + rowSums(numer)
    for (s in 1:m){probs[, s] <- numer[, s] / denom}
    probs[, m + 1] <- 1 / denom    # baseline J=0 prob
    probs
  }
  
  # dat=comp_risk %>% mutate(S = case_when(
  #   type=="Omicron" ~ 2,
  #   type == "Delta" ~ 1,
  #   .default = 0
  # )) %>% 
  #   mutate(J=ifelse(J==0, 0, 1))
  # J.name="J"
  # V.name="V"
  # T.name="T"
  # V.early.name="V_early"
  # S.name="S"
  # #df=6
  # m=2
  # d=3
  # p_s_fun = predict_GISAID
  # init_beta=NULL
  
  # dat must contain: J (0/1), T, V, S (strain index 1..m for J=1, NA or 0 for J=0),
  # S only used when estimating p_s from unvaccinated
  n <- nrow(dat) # number of observations
  J <- dat[[J.name]] # Vaccine preventable/irrelevant event labels
  Tvec <- dat[[T.name]] # Event times
  Vvec <- dat[[V.name]] # Vaccination times
  Svec <- dat[[S.name]] # Strain labels (0 if off-target)
  
  m=m
  
  bases<-vector(mode="list", length=m)
  
  for(i in 1:m){
    if(!is.null(V.early.name)){
      Vvec_early <- dat[[V.early.name]]
      #bases[[i]]<-psi_delta(Tvec, Vvec, Vvec_early)
      bases[[i]]<-psi_delta(Tvec, Vvec, Vvec_early, ...)
    }else{
      bases[[i]]<-psi_delta(Tvec, Vvec, ...)
    }
  }
  
  Z_stack <- do.call(cbind, lapply(1:length(bases),
                                   FUN=function(x){
                                     bases[[x]]$mat
                                   }))
  
  p_tot <- ncol(Z_stack)
  d = p_tot/m
  
  # initial alpha (link) estimate: smooth J~s(T) on link scale and center
  fit_alpha0 <- gam(J ~ s(`T`), dat = bind_cols(T=Tvec, J=J, Z = I(Vvec_early <= Tvec)) %>% filter(Z==0), family = binomial)
  alpha_link <- predict(fit_alpha0, dat, type = "link")
  alpha_link <- alpha_link
  
  if(is.null(init_beta)){beta<-rep(0, p_tot)}else{beta<-init_beta}
  
  probs <- compute_probs(Z_stack, beta, alpha_link, p_s_mixture_vecs = predict_GISAID(dat[[T.name]]))
  Q0s_mat <- probs[,1:m] #Q0s_mat
  Q0 <- rowSums(probs[,1:m])
  
  if(!is.null(bandwidth)){
    bw = bandwidth
  }else if(is.null(bandwidth) & automatic_bw==FALSE){
    message("No bandwidth for nuisance estimation provided")
    message("Using Silverman's rule of thumb")
    bw = 1.06 * sd(dat[[T.name]]) * n^(-1/5)
  }else if(is.null(bandwidth) & automatic_bw==TRUE){
    message("No bandwidth for nuisance estimation provided")
    message("Using LOOCV to select bandwidth")
    # --- helper: gaussian kernel
    K <- function(u, kernel){
      if(kernel=="gaussian"){dnorm(u)}
      if(kernel=="epanechnikov"){
        ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
      }
    }
    
    # --- LOO NW estimator for scalar response y at points T (vectorized)
    # returns vector of length n with leave-one-out NW estimates
    loo_nw_estimates <- function(Tvec, yvec, h, shape) {
      n <- length(Tvec)
      # distance matrix (can be large for big n; consider block or FFT acceleration for very large n)
      # u_ij = (T_j - T_i) / h  (we predict at i using all j != i)
      # We'll compute numerator_i = sum_{j!=i} K(u_ij) * y_j ; denom_i = sum_{j!=i} K(u_ij)
      diffs <- outer(Tvec, Tvec, "-") / h
      W <- K(diffs, shape)
      diag(W) <- 0                     # remove j=i
      denom <- rowSums(W)
      # avoid zero denom
      denom[denom == 0] <- 1e-12
      numer <- W %*% yvec
      as.numeric(numer / denom)
    }
    
    # --- CV score for a single candidate h
    cv_score_for_h <- function(h, Tvec, Vvec, Q0s_mat, Q0vec, psi_list, shape) {
      n <- length(Tvec)
      m <- ncol(Q0s_mat)
      # build responses:
      # denominator z_i
      z <- Q0vec * (1 - Q0vec)
      # numerator components: for each strain s and each basis dimension j
      total_sse <- 0
      # CV error for numerator components
      for (s in seq_len(m)) {
        psi_mat <- psi_list[[s]]$mat  # n x p
        if (is.null(dim(psi_mat))) psi_mat <- matrix(psi_mat, ncol = 1)
        p <- ncol(psi_mat)
        # compute y_{i}^{(s,j)} = I(V<=T) * psi_j * Q0s * (1-Q0)
        # Iv <- as.numeric(Vvec <= Tvec)
        Q0s <- Q0s_mat[, s]
        # for each column j do LOO NW predict and accumulate SSE
        for (j in seq_len(p)) {
          yj <- psi_mat[, j] * Q0s * (1 - Q0vec)
          yhat_loo <- loo_nw_estimates(Tvec, yj, h, shape)
          total_sse <- total_sse + sum( (yj - yhat_loo)^2 )
        }
      }
      # CV error for denominator
      zhat_loo <- loo_nw_estimates(Tvec, z, h, shape)
      total_sse <- total_sse + sum( (z - zhat_loo)^2 )
      total_sse
    }
    
    # --- choose bandwidth on a grid around pilot
    select_bandwidth <- function(Tvec, Vvec, Q0s_mat, Q0vec, psi_list,
                                 grid = NULL, ngrid = 12, shape) {
      n <- length(Tvec)
      pilot <- 1.06 * sd(Tvec) * n^(-1/5)
      if (is.null(grid)) {
        grid <- exp(seq(log(pilot*0.25), log(pilot*2.5), length.out = ngrid))
      }
      scores <- sapply(grid, function(h) cv_score_for_h(h, Tvec, Vvec, Q0s_mat, Q0vec, psi_list, shape))
      best_h <- grid[which.min(scores)]
      list(h = best_h, grid = grid, scores = scores)
    }
    
    bw = select_bandwidth(Tvec, Vvec, Q0s_mat = Q0s_mat, Q0vec = Q0, psi_list = bases, ngrid=25, shape="epanechnikov")$h
  }
  
  for(iter in 1:max_iter){
    
    if(iter==1){
    beta_new <- beta
    }
    probs <- compute_probs(Z_stack, beta_new, alpha_link, p_s_mixture_vecs = predict_GISAID(dat[[T.name]]))
    Q0s_mat <- probs[,1:m] #Q0s_mat
    Q0 <- rowSums(probs[,1:m])
  
    a_hat_list <- vector(mode="list", length=m)
    H_list <- vector(mode="list", length(m))
    offset_list <- vector(mode="list", length(m))
    
    for(s in seq_len(m)){
  
      p0s = Q0s_mat[,s]
    
      # 2) compute offset (eta for strain s) from current probs (use formula below for numerical safety)
      offset_vec <- log( safe_clamp(p0s) / safe_clamp(1 - Q0) )  # n-vector
    
      offset_list[[s]]<-offset_vec
      
      a_mat <- matrix(0, nrow = n, ncol = d)
    
      for (i in seq_len(n)) {
        t0 <- Tvec[i]
        u <- (Tvec - t0) / bw
        if (kernel_shape == "gaussian") {
          w <- dnorm(u)
        } else if (kernel_shape == "epanechnikov") {
          w <- ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
        } else w <- dnorm(u)
        if (sum(w) == 0) w <- rep(1/n, n)
        w <- w / sum(w)
      
        denom <- sum(w * Q0 * (1-Q0))
        # numerator is vector: sum_w [ I(V<=t) * psi(t-V) * Q0s * (1-Q0) ]
        psi_weighted <- bases[[s]]$mat * ((Vvec <= t0) * p0s * (1 - Q0) * w )  # broadcasting per row
        numer_vec <- colSums(psi_weighted)   # length p
        if (denom > 0) a_mat[i, ] <- numer_vec / denom else a_mat[i, ] <- 0
      }
    
      a_hat_list[[s]] <- a_mat
    
      H = bases[[s]]$mat - a_mat
    
      H <- data.frame(H)
      colnames(H)<-paste0("H_", s, "_", seq_len(ncol(H)))
      H_list[[s]]<-H
    }
  
    # fit_init <- gam(list(
    # S~0 + H_1_1 + H_1_2 + H_1_3 + offset(log(p_Delta)),
    # ~0 + H_2_1 + H_2_2 + H_2_3 + offset(log(p_Omicron)),
    # 1+2~ s(`T`)), 
    # data=bind_cols(dat, H_list[[1]], H_list[[2]]), 
    # family=multinom(K=2))
    
    formula_list = vector(mode="list", length=m)
    
    for(i in 1:length(formula_list)){
      
      if(i==1){
        
        formula_list[[i]]<-as.formula(paste0("S~0 + ", 
                          paste0(colnames(H_list[[i]]), collapse="+"),
                          paste0("+offset(off", i, ")"), collapse=""))
        
      }else{
        
        formula_list[[i]]<-as.formula(paste0("~0 + ", 
                                             paste0(colnames(H_list[[i]]), collapse="+"),
                                             paste0("+offset(off", i, ")"), collapse=""))
        
      }
      
    }
    
    for(i in 1:m){
      if(i==1){
        dat_tmp = bind_cols(dat, H_list[[i]])
        dat_tmp[[paste0("off", i)]]<-offset_list[[i]]
      }else{
        dat_tmp = bind_cols(dat_tmp, H_list[[i]])
        dat_tmp[[paste0("off", i)]]<-offset_list[[i]]
      }
    }
    
    fit_init <- gam(formula_list, 
      data=dat_tmp, 
      family=multinom(K=2))
  
    eps_hat <- coef(fit_init)[1:p_tot]
    
    beta_new <- beta_new + eps_hat
    
    if (sqrt(sum(eps_hat^2)) < tol) {
      beta <- beta_new
      if (verbose) cat("Converged at iter", iter, "\n")
      converged <- TRUE
      break
    }
    
    if (verbose) cat(sprintf("iter %d: ||eps|| = %.3e\n", iter, sqrt(sum(eps_hat^2))))
  }
  
  list_out <- list(
    beta = beta
  )
  
  for(s in 1:m){
    
    scores = bases[[s]]$mat * (I(dat$S==s) - Q0s_mat[,s]) - a_hat_list[[s]] * (dat$J - Q0)
    
    Ihat = crossprod(scores)/n
    
    invI <- try(solve(Ihat), silent = TRUE)
    if (inherits(invI, "try-error") || rcond(Ihat) < 1e-12) {
      invI <- try(solve(Ihat + lambda * diag(p)), silent = TRUE)
      if (inherits(invI, "try-error")) invI <- MASS::ginv(Ihat)
    }
    
    EIC_mat <- t(invI %*% t(scores))
    
    var_beta_s = cov(EIC_mat)/n
    
    list_out[[paste0("var_beta", s)]]<-var_beta_s
  }
  
  return(
    list_out
  )
  
}
