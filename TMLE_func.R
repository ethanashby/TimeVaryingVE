# Ready-to-run iterative TMLE for the semiparametric logistic model
# Dependencies: mgcv
if(!requireNamespace("mgcv", quietly=TRUE)) install.packages("mgcv")
library(tidyverse)
library(mgcv)
library(splines)
library(coneproj)

psi_d2 <- function(Tvec, Vvec) {
  # Vmat: n x d
  # default: psi(T-V) = V
  return(cbind(as.numeric(Tvec - Vvec > 0), pmax(0, Tvec-Vvec)))
}

psi_bs <- function(Tvec, Vvec, df=3) {
  # Vmat: n x d
  # default: psi(T-V) = V
  return(cbind(as.numeric(Tvec - Vvec > 0), bs(pmax(0, Tvec-Vvec), df=df, degree = 3)))
}

tmle_iterative <- function(dat, J.name = "J", T.name = "T", V.name="V", psi_delta, maxit = 500, tol = 1e-6,
                           smooth_r = TRUE, smooth_alpha = TRUE, verbose = TRUE) {
  
  #dat=dat_wane_run
  
  # dat: data.frame with columns J (0/1), T, V
  # psi_delta: function(delta) -> numeric vector length d
  # d: dimension of psi
  n <- nrow(dat)
  J <- as.numeric(dat[[J.name]])
  Tvec <- dat[[T.name]]
  Vvec <- dat[[V.name]]
  
  # build X matrix (n x d)
  
  Z<-psi_delta(Tvec, Vvec)
  
  # build_Z <- function() {
  #   X <- matrix(0, nrow = n, ncol = d)
  #   for (i in seq_len(n)) {
  #     if (!is.na(Vvec[i]) && Vvec[i] <= Tvec[i]) {
  #       X[i, ] <- psi_delta(Tvec[i], Vvec[i])
  #     } else {
  #       X[i, ] <- rep(0, d)
  #     }
  #   }
  #   X
  # }
  
  #Z <- build_Z()
  d <- ncol(Z)
  
  # initial estimates
  beta <- rep(0, d)                 # initialize beta to 0
  
  # initial alpha: fit GAM of J ~ s(T) (link scale)
  fit_a0 <- gam(J ~ s(Tvec), family = binomial)
  alpha_hat <- predict(fit_a0, type = "link")
  alpha_hat <- alpha_hat - mean(alpha_hat)  # identifiability: center alpha (mean 0)
  eta <- as.numeric(as.vector(Z %*% beta) + alpha_hat)
  p <- plogis(eta)
  
  iter <- 0
  converged <- FALSE
  for (k in seq_len(maxit)) {
    iter <- k
    w <- p * (1 - p)                  # weight vector length n
    
    # estimate r_j(T) for each column j
    r_mat <- matrix(0, nrow = n, ncol = d)
    for (j in seq_len(d)) {
      # If no smoothing desired (e.g. discrete T with many ties), we can do weighted average by T
      if (!smooth_r) {
        # discrete-level weighted means by T values
        tmp <- tapply(Z[, j] * w, Tvec, sum)
        denom <- tapply(w, Tvec, sum)
        # map back to observations
        denom[is.na(denom)] <- 0
        tmp[is.na(tmp)] <- 0
        r_byT <- numeric(length(denom))
        r_byT[] <- tmp / ifelse(denom == 0, 1, denom)
        # mapping: convert names back to T values
        Tlevels <- as.numeric(names(denom))
        # create named vector
        r_map <- setNames(r_byT, Tlevels)
        r_mat[, j] <- r_map[as.character(Tvec)]
        r_mat[is.na(r_mat[, j]), j] <- 0
      } else {
        # smooth X[,j] on T with weights w using gam (identity link)
        df_tmp <- data.frame(Xcol = Z[, j], T = Tvec)
        # Add a tiny ridge penalty if all weights are zero
        w_use <- w
        if (all(w_use == 0)) w_use <- rep(1e-6, n)
        fit_r <- try(gam(Xcol ~ s(`T`, k = min(20, max(5, length(unique(Tvec))/4))), 
                         data = df_tmp, weights = w_use), silent = TRUE)
        if (inherits(fit_r, "try-error")) {
          # fallback to loess
          r_mat[, j] <- predict(loess(X[, j] ~ Tvec, weights = w_use), newdata = Tvec)
        } else {
          r_mat[, j] <- predict(fit_r, newdata = df_tmp)
        }
        r_mat[is.na(r_mat[, j]), j] <- 0
      }
    }
    
    H <- Z - r_mat   # n x d matrix
    
    # Attempt logistic fluctuation with glm: J ~ -1 + H, offset = qlogis(p)
    # Prepare data frame
    df_fluct <- as.data.frame(H)
    names(df_fluct) <- paste0("H", seq_len(d))
    df_fluct$J <- J
    offset_vec <- qlogis(p)
    
    # Build formula programmatically
    formula_str <- paste("J ~ -1 +", paste(names(df_fluct)[1:d], collapse = "+"))
    fit_eps <- try(glm(as.formula(formula_str), family = binomial, data = df_fluct, offset = offset_vec), silent = TRUE)
    
    if (inherits(fit_eps, "try-error") || any(is.na(coef(fit_eps)))) {
      # Fallback: one-step Newton update (closed form)
      # Solve: eps = (Sum_i w_i H_i H_i^T)^{-1} Sum_i H_i (J_i - p_i)
      W_mat <- 0
      for (i in seq_len(n)) {
        Hi <- matrix(H[i, ], ncol = 1)
        W_mat <- W_mat + (w[i]) * (Hi %*% t(Hi))
      }
      # RHS
      rhs <- colSums(H * (J - p))
      # regularize if near-singular
      reg <- 1e-8 * diag(d)
      eps_hat <- try(as.vector(solve(W_mat + reg, rhs)), silent = TRUE)
      if (inherits(eps_hat, "try-error") || any(is.na(eps_hat))) {
        eps_hat <- rep(0, d)
      }
      if (verbose) cat("Used Newton fallback for fluctuation at iter", k, "\n")
    } else {
      eps_hat <- coef(fit_eps)
      # if names mismatch, coerce
      eps_hat <- as.numeric(eps_hat)
    }
    
    # update: Option A (simple): beta <- beta + eps
    beta_new <- beta + eps_hat
    
    # update logits and alpha profile: eta_new = logit(p) + H eps
    eta_new <- qlogis(p) + as.numeric(H %*% eps_hat)
    alpha_hat_new <- eta_new - as.numeric(Z %*% beta_new)
    
    # smooth/center alpha if requested
    if (smooth_alpha) {
      df_alpha <- data.frame(alpha = alpha_hat_new, T = Tvec)
      fit_alpha <- try(gam(alpha ~ s(`T`, k = min(20, max(5, length(unique(Tvec))/4))), data = df_alpha), silent = TRUE)
      if (inherits(fit_alpha, "try-error")) {
        alpha_hat <- alpha_hat_new - mean(alpha_hat_new)  # at least center
      } else {
        alpha_hat <- predict(fit_alpha, type = "response")
        # ensure link-scale (we used identity for alpha), center for identifiability
        alpha_hat <- alpha_hat - mean(alpha_hat)
      }
    } else {
      alpha_hat <- alpha_hat_new - mean(alpha_hat_new)
    }
    
    # update p
    eta <- as.numeric(as.vector(Z %*% beta_new) + alpha_hat)
    p <- plogis(eta)
    
    # convergence check
    if (sqrt(sum(eps_hat^2)) < tol) {
      beta <- beta_new
      if (verbose) cat("Converged at iter", k, "\n")
      converged <- TRUE
      break
    }
    # assign for next iter
    beta <- beta_new
    if (verbose) cat(sprintf("iter %d: ||eps|| = %.3e\n", k, sqrt(sum(eps_hat^2))))
  } # end for
  
  # Final H and p (recompute)
  w_final <- p * (1 - p)
  # recompute final r_mat (same code as above but single pass)
  r_mat <- matrix(0, nrow = n, ncol = d)
  for (j in seq_len(d)) {
    if (!smooth_r) {
      tmp <- tapply(Z[, j] * w_final, Tvec, sum)
      denom <- tapply(w_final, Tvec, sum)
      tmp[is.na(tmp)] <- 0
      denom[is.na(denom)] <- 0
      r_byT <- tmp / ifelse(denom == 0, 1, denom)
      r_map <- setNames(r_byT, as.numeric(names(denom)))
      r_mat[, j] <- r_map[as.character(Tvec)]
      r_mat[is.na(r_mat[, j]), j] <- 0
    } else {
      df_tmp <- data.frame(Xcol = Z[, j], T = Tvec)
      w_use <- w_final
      if (all(w_use == 0)) w_use <- rep(1e-6, n)
      fit_r <- try(gam(Xcol ~ s(T, k = min(20, max(5, length(unique(Tvec))/4))), data = df_tmp, weights = w_use), silent = TRUE)
      if (inherits(fit_r, "try-error")) {
        r_mat[, j] <- predict(loess(X[, j] ~ Tvec, weights = w_use), newdata = Tvec)
      } else {
        r_mat[, j] <- predict(fit_r, newdata = df_tmp)
      }
      r_mat[is.na(r_mat[, j]), j] <- 0
    }
  }
  H_final <- Z - r_mat
  # estimate information matrix: I_eff_hat = mean_i [ w_final[i] * H_i %o% H_i ] where w_final = p*(1-p)
  I_eff_hat <- matrix(0, nrow = d, ncol = d)
  for (i in seq_len(n)) {
    Hi <- matrix(H_final[i, ], ncol = 1)
    I_eff_hat <- I_eff_hat + w_final[i] * (Hi %*% t(Hi))
  }
  I_eff_hat <- I_eff_hat / n
  
  # variance estimate for beta: var(beta_hat) = (I_eff_hat)^{-1} / n
  reg <- 1e-8 * diag(d)
  I_inv <- try(solve(I_eff_hat + reg), silent = TRUE)
  if (inherits(I_inv, "try-error")) {
    I_inv <- MASS::ginv(I_eff_hat + reg)
  }
  var_beta <- I_inv / n
  se_beta <- sqrt(diag(var_beta))
  
  # return
  list(beta = beta,
       se = se_beta,
       cov = var_beta,
       alpha_link = alpha_hat,
       fitted_p = p,
       H = H_final,
       I_eff = I_eff_hat,
       iterations = iter,
       converged = converged)
}


res<-tmle_iterative(dat=dat_wane_run, psi_delta = psi_d2)
res2<-tmle_iterative(dat=dat_wane_run, psi_delta = psi_bs)

try<-data.frame(A= as.numeric(Vvec <= Tvec), tau = pmax(0, Tvec-Vvec)) 
try<-bind_cols(try, psi_bs(Tvec, Vvec))
try$VE <- as.matrix(try) %*% c(0,0,res2$beta)

ggplot(aes(x=tau, y=VE), data=try)+
  geom_line()+
  geom_abline(intercept=-1, slope=1)
