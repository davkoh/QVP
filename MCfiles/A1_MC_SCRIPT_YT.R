
## ---- Data generator: MC_create_BRW -----------------------------------
MC_create_BRW <- function(T, const, tau, quantspec_cutoff = 0.1) {
  if (quantspec_cutoff > 0.5) quantspec_cutoff <- 1 - quantspec_cutoff
  
  T <- T + 100
  ynum <- 5
  max_beta <- 10
  
  beta_mat <- array(NA_real_, dim = c(max_beta + 1, length(tau), ynum))
  Data_y <- matrix(NA_real_, nrow = T, ncol = ynum)
  Data_x <- array(NA_real_, dim = c(T, max_beta, ynum))
  
  for (ydes in seq_len(ynum)) {
    beta_full <- matrix(NA_real_, nrow = 11, ncol = length(tau))
    beta_full[1, ] <- qnorm(tau) + 1
    eps <- rnorm(T)
    
    if (ydes == 1) {
      for (row in 2:(4 + 1)) beta_full[row, ] <- qnorm(tau) * 0.1 + 1
      X <- matrix(runif(T * 4), nrow = T, ncol = 4)
      beta_use <- c(1, 1, 1, 1)
      theta_use <- c(0.1, 0.1, 0.1, 0.1)
      
    } else if (ydes == 2) {
      for (row in 2:nrow(beta_full)) {
        if (row <= 5) beta_full[row, ] <- qnorm(tau) * 0.1 + 1
        else beta_full[row, ] <- qnorm(tau) * 0
      }
      X <- matrix(runif(T * 10), nrow = T, ncol = 10)
      beta_use <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
      theta_use <- c(0.1, 0.1, 0.1, 0.1, 0, 0, 0, 0, 0, 0)
      
    } else if (ydes == 3) {
      for (row in 2:(7 + 1)) {
        if (row <= 4) beta_full[row, ] <- qnorm(tau) * 0.1 + 1
        else beta_full[row, ] <- qnorm(tau) * 0 + 1
      }
      X <- matrix(runif(T * 7), nrow = T, ncol = 7)
      beta_use <- c(1, 1, 1, 1, 1, 1, 1)
      theta_use <- c(0.1, 0.1, 0.1, 0, 0, 0, 0)
      
    } else if (ydes == 4) {
      tau_r <- round(tau, 3)
      qcut_r <- round(quantspec_cutoff, 3)
      quant_spec_mat <- as.numeric(tau_r <= qcut_r) +
        as.numeric(tau_r >= round(1 - quantspec_cutoff, 3))
      
      for (row in 2:nrow(beta_full)) {
        if (row <= 5) beta_full[row, ] <- qnorm(tau) * 0.1 + 1
        else beta_full[row, ] <- qnorm(tau) * quant_spec_mat
      }
      
      beta_full[(nrow(beta_full) - 1):nrow(beta_full), ] <- 0
      
      X <- matrix(runif(T * 10), nrow = T, ncol = 10)
      beta_use <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
      theta_use <- c(0.1, 0.1, 0.1, 0.1, 1, 1, 1, 1, 0, 0)
      
      qcheck <- as.numeric(eps <= qnorm(quantspec_cutoff)) +
        as.numeric(eps >= qnorm(1 - quantspec_cutoff))
      
    } else if (ydes == 5) {
      for (row in 2:(4 + 1)) {
        if (row <= 3) beta_full[row, ] <- qnorm(tau) * 0.1 + 1
        else beta_full[row, ] <- qnorm(tau) * 1 + 1
      }
      X <- matrix(runif(T * 4), nrow = T, ncol = 4)
      beta_use <- c(1, 1, 1, 1)
      theta_use <- c(0.1, 0.1, 1, 1)
    }
    
    beta_mat[, , ydes] <- beta_full
    
    if (ydes == 4) {
      idx01 <- which(theta_use == 0.1)
      idx1 <- which(theta_use == 1)
      
      part01 <- if (length(idx01) > 0) drop(X[, idx01, drop = FALSE] %*% theta_use[idx01]) else rep(0, T)
      part1 <- if (length(idx1) > 0) drop(X[, idx1, drop = FALSE] %*% theta_use[idx1]) else rep(0, T)
      
      y <- const +
        drop(X %*% beta_use) +
        (1 + part01 + qcheck * part1) * eps
    } else {
      y <- const +
        drop(X %*% beta_use) +
        (1 + drop(X %*% theta_use)) * eps
    }
    
    Data_y[, ydes] <- y
    Data_x[, 1:ncol(X), ydes] <- X
  }
  
  list(Y = Data_y, X = Data_x, coeff_full = beta_mat)
}

## ---- wQS (MATLAB -> R) ------------------------------------
wQS <- function(TrueY, ForecastQuants, Quantiles, weighting = 1L) {
  TrueY <- as.numeric(TrueY)
  T <- length(TrueY)
  qs <- length(Quantiles)
  
  if (!is.matrix(ForecastQuants)) ForecastQuants <- as.matrix(ForecastQuants)
  stopifnot(nrow(ForecastQuants) == T, ncol(ForecastQuants) == qs)
  
  resid <- TrueY - ForecastQuants
  ticklossw <- matrix(Quantiles, nrow = T, ncol = qs, byrow = TRUE) -
    (TrueY <= ForecastQuants)
  QS <- ticklossw * resid
  
  if (weighting == 1L) {
    QStemp <- QS
  } else if (weighting == 2L) {
    w <- Quantiles * (1 - Quantiles)
    QStemp <- sweep(QS, 2, w, `*`)
  } else if (weighting == 3L) {
    w <- (1 - Quantiles)^2
    QStemp <- sweep(QS, 2, w, `*`)
  } else if (weighting == 4L) {
    w <- (Quantiles)^2
    QStemp <- sweep(QS, 2, w, `*`)
  } else {
    stop("weighting must be 1,2,3,4")
  }
  
  output <- rowMeans(QStemp)
  list(output = output, output_full = QS)
}

## ---- psrf (MATLAB -> R) -----------------------------------
psrf <- function(X) {
  if (length(dim(X)) != 3) stop("X must be a 3D array: N x D x M")
  N <- dim(X)[1]; D <- dim(X)[2]; M <- dim(X)[3]
  if (N < 1) stop("Too few samples")
  if (M < 2) stop("psrf requires at least 2 chains (M >= 2)")
  
  # Within-chain variances averaged over chains
  W <- rep(0, D)
  for (m in seq_len(M)) {
    Xm <- X[, , m]
    Xm_centered <- sweep(Xm, 2, colMeans(Xm), `-`)
    W <- W + colSums(Xm_centered * Xm_centered)
  }
  W <- W / ((N - 1) * M)
  
  # Between-chain variance of means (B / N in many references)
  chain_means <- vapply(seq_len(M), function(m) colMeans(X[, , m]), numeric(D))
  chain_means <- t(chain_means)           # M x D
  mbar <- colMeans(chain_means)
  
  Bpn <- rep(0, D)
  for (m in seq_len(M)) {
    diff <- chain_means[m, ] - mbar
    Bpn <- Bpn + diff * diff
  }
  Bpn <- Bpn / (M - 1)
  
  S <- (N - 1) / N * W + Bpn
  R2 <- (M + 1) / M * S / W - (N - 1) / (M * N)
  V <- R2 * W
  R <- sqrt(R2)
  
  B <- Bpn * N
  neff <- pmin(M * N * V / B, M * N)
  
  list(R = R, neff = neff, V = V, W = W, B = B)
}

## ---- QVP_wrapper ----------------------------------
QVP_wrapper_yang_todkar <- function(y, X, tau, iter, burn, chain_i,
                        arg7 = 1, arg8 = 1,
                        incr = 0.05, thin = 1, shrink = TRUE) {
  t0 <- proc.time()[["elapsed"]]
  
  y <- as.numeric(y)
  X <- as.matrix(X)
  
  # Data for formula y ~ .
  df <- data.frame(y = y, as.data.frame(X))
  names(df) <- make.names(names(df), unique = TRUE)
  
  Q <- length(tau)
  draws_keep <- iter - burn
  if (draws_keep <= 0) stop("iter must be > burn")
  
  # Choose nsamp so that after burn-in and thinning we keep exactly (iter-burn) draws
  # Keep indices: (burn + thin), (burn + 2*thin), ..., (burn + draws_keep*thin)
  nsamp <- burn + draws_keep * thin
  
  fit <- qrjoint::qrjoint(
    y ~ .,
    data = df,
    nsamp = nsamp,
    thin = thin,        # we control thinning ourselves so the indexing is explicit
    incr = incr,
    shrink = shrink
  )
  
  betas <- coef(fit, plot = FALSE, nmc = nsamp)
  beta_samp <- betas$beta.samp
  if (is.null(beta_samp)) stop("coef(fit, plot=FALSE)$beta.samp was NULL")
  beta_samp <- beta_samp[2:(dim(beta_samp)[1]-1),,1:(iter-burn)] # dim(beta_samp) = Q x K x iter
  
  # Posterior mean (coef x Q)
  bhat <- t(apply(beta_samp, c(1, 2), mean, na.rm = TRUE))
  
  # full_draw is (coef x Q x draws_keep) which matches:
  full_draw <- aperm(beta_samp, c(2,1,3))
  
  # SAVS not used for qrjoint; return same bhat for compatibility
  bhat_SAVS <- bhat
  
  time <- proc.time()[["elapsed"]] - t0
  list(bhat = bhat, full_draw = full_draw, bhat_SAVS = bhat_SAVS, time = time)
}
## ---- Estimation driver (translated) ----------------------------------
  tau <- seq(0.05, 0.95, by = 0.05)
  iter <- 15000
  burn <- 10000
  MC_T <- 100
  MC_num <- 1
  rho <- 0.0
  Chain_num <- 4
  
  Prior_mat <- c("YT")
  
  set.seed(1998)
  
  # Rearrange Prior_mat so SAVS is always saved for LQVP_HS and LQVP_HS_ASIS
  loc <- !grepl("_SAVS", Prior_mat, fixed = TRUE)
  Prior_temp <- Prior_mat[loc]
  ASIS_flg <- sum(grepl("LQVP_HS_ASIS", Prior_temp, fixed = TRUE))
  LQVP_flg <- sum(grepl("LQVP_HS", Prior_temp, fixed = TRUE)) - ASIS_flg
  Prior_temp <- Prior_temp[!grepl("LQVP_HS", Prior_temp, fixed = TRUE)]
  
  if (LQVP_flg == 1 && ASIS_flg == 1) {
    Prior_mat <- c(Prior_temp, "LQVP_HS","LQVP_HS_SAVS","LQVP_HS_ASIS","LQVP_HS_ASIS_SAVS")
  } else if (LQVP_flg == 1 && ASIS_flg == 0) {
    Prior_mat <- c(Prior_temp, "LQVP_HS","LQVP_HS_SAVS")
  } else if (LQVP_flg == 0 && ASIS_flg == 1) {
    Prior_mat <- c(Prior_temp, "LQVP_HS_ASIS","LQVP_HS_ASIS_SAVS")
  } else {
    Prior_mat <- Prior_temp
  }
  
  results <- vector("list", MC_num)
  
  for (MC_it in seq_len(MC_num)) {
    datause <- MC_create_BRW(MC_T, 1, tau)
    
    True_beta <- datause$coeff_full
    case_tot <- dim(datause$X)[3]
    Prio_mat_length <- length(Prior_mat)
    
    bhat_ALL <- array(NA_real_, dim = c(dim(True_beta), Prio_mat_length))
    bhat_deviation <- array(NA_real_, dim = c(dim(True_beta), Prio_mat_length))
    Rhat <- array(NA_real_, dim = c(dim(True_beta), Prio_mat_length))
    Neff <- array(NA_real_, dim = c(dim(True_beta), Prio_mat_length))
    wQS_ALL <- array(NA_real_, dim = c(4, dim(True_beta)[3], Prio_mat_length))
    Cross_ALL <- array(NA_real_, dim = c(dim(True_beta)[2], dim(True_beta)[3], Prio_mat_length))
    Time_mat <- array(NA_real_, dim = c(case_tot, Prio_mat_length, Chain_num))
    
    Prelocated_y <- matrix(NA_real_, nrow = MC_T, ncol = dim(True_beta)[3])
    Prelocated_x <- array(NA_real_, dim = c(MC_T, dim(True_beta)[1] - 1, dim(True_beta)[3]))
    Forecast_y <- matrix(NA_real_, nrow = 100, ncol = dim(True_beta)[3])
    Forecast_x <- array(NA_real_, dim = c(100, dim(True_beta)[1] - 1, dim(True_beta)[3]))
    
    datause <- MC_create_BRW(MC_T, 1, tau)
    
    Prelocated_y[,] <- datause$Y[1:MC_T, , drop = FALSE]
    Prelocated_x[,,] <- datause$X[1:MC_T, , , drop = FALSE]
    
    Forecast_y[,] <- datause$Y[(MC_T + 1):nrow(datause$Y), , drop = FALSE]
    Forecast_x[,,] <- datause$X[(MC_T + 1):dim(datause$X)[1], , , drop = FALSE]
    
    for (Prior_i in seq_len(Prio_mat_length)) {
      Prior_choice <- Prior_mat[Prior_i]
      
      if (!(Prior_choice %in% c("LQVP_HS_SAVS", "LQVP_HS_ASIS_SAVS"))) {
        
        Rhat_temp <- array(NA_real_, dim = dim(True_beta))
        Neff_temp <- array(NA_real_, dim = dim(True_beta))
        
        bhat_ALL_temp_UNSPARS <- array(NA_real_, dim = dim(True_beta))
        bhat_ALL_temp_SAVS <- array(NA_real_, dim = dim(True_beta))
        bhat_deviation_temp_UNSPARS <- array(NA_real_, dim = dim(True_beta))
        bhat_deviation_temp_SAVS <- array(NA_real_, dim = dim(True_beta))
        
        wQS_temp_UNSPARS <- matrix(NA_real_, nrow = 4, ncol = case_tot)
        wQS_temp_SAVS <- matrix(NA_real_, nrow = 4, ncol = case_tot)
        
        Cross_temp_UNSPARS <- matrix(NA_real_, nrow = dim(True_beta)[2], ncol = case_tot)
        Cross_temp_SAVS <- matrix(NA_real_, nrow = dim(True_beta)[2], ncol = case_tot)
        
        True_beta_temp <- True_beta
        
        for (case_num in seq_len(case_tot)) {
          ysmall <- Prelocated_y[, case_num]
          Xsmall <- Prelocated_x[, , case_num, drop = TRUE]
          Xsmall <- as.matrix(Xsmall)
          Xsmall <- Xsmall[, colSums(is.na(Xsmall)) == 0, drop = FALSE]
          
          p <- ncol(Xsmall)
          Q <- length(tau)
          draws <- iter - burn
          
          bhat_chains <- array(NA_real_, dim = c(p + 1, Q, Chain_num))
          bhat_chains_SAVS <- array(NA_real_, dim = c(p + 1, Q, Chain_num))
          posterior_ALL_temp <- array(NA_real_, dim = c(p + 1, Q, draws, Chain_num))
          time_temp <- rep(NA_real_, Chain_num)
          
          for (chain_i in seq_len(Chain_num)) {
            message(sprintf("Draws for model: %s, Sim-Case: %i, Chain: %i",
                            Prior_choice, case_num, chain_i))
            
            res <- tryCatch(
              QVP_wrapper_yang_todkar(ysmall, Xsmall, tau, iter, burn, chain_i, 1, 1,thin =2),
              error = function(e) NULL
            )
            
            
            if (is.null(res)) {
              bhat <- matrix(NA_real_, nrow = p + 1, ncol = Q)
              bhat_SAVS <- matrix(NA_real_, nrow = p + 1, ncol = Q)
              full_draw <- array(NA_real_, dim = c(p + 1, Q, draws))
              time <- NA_real_
            } else {
              bhat <- res$bhat
              bhat_SAVS <- res$bhat_SAVS
              full_draw <- res$full_draw
              time <- res$time
            }
            
            bhat_chains[, , chain_i] <- bhat
            bhat_chains_SAVS[, , chain_i] <- bhat_SAVS
            posterior_ALL_temp[, , , chain_i] <- full_draw
            time_temp[chain_i] <- time
          }
          
          Time_mat[case_num, Prior_i, ] <- time_temp
          if (Prior_choice %in% c("LQVP_HS", "LQVP_HS_ASIS")) {
            if (Prior_i + 1 <= Prio_mat_length) Time_mat[case_num, Prior_i + 1, ] <- time_temp
          }
          
          # Convergence diagnostics per quantile
          Rhat_mat <- matrix(NA_real_, nrow = dim(bhat_chains)[1], ncol = dim(bhat_chains)[2])
          Neff_mat <- matrix(NA_real_, nrow = dim(bhat_chains)[1], ncol = dim(bhat_chains)[2])
          
          for (quant_i in seq_along(tau)) {
            beta_temp <- posterior_ALL_temp[, quant_i, , , drop = FALSE]
            beta_temp <- array(beta_temp, dim = c(p + 1, draws, Chain_num))
            beta_temp <- aperm(beta_temp, c(3, 1, 2))  # chain x param x draw
            
            diag <- psrf(beta_temp)
            Rhat_mat[, quant_i] <- diag$R
            Neff_mat[, quant_i] <- diag$neff
          }
          
          Rhat_temp[1:(p + 1), , case_num] <- Rhat_mat
          Neff_temp[1:(p + 1), , case_num] <- Neff_mat
          
          # UNSPARS and SAVS variants
          for (bhat_it in 1:2) {
            bhat_chain_use <- if (bhat_it == 1) bhat_chains else bhat_chains_SAVS
            bhat <- apply(bhat_chain_use, c(1, 2), mean, na.rm = TRUE)
            
            bhat_save <- matrix(NA_real_, nrow = dim(True_beta_temp)[1], ncol = dim(True_beta_temp)[2])
            bhat_save[1:nrow(bhat), ] <- bhat
            
            True_beta_use <- True_beta_temp[, , case_num]
            True_beta_use <- True_beta_use[!apply(is.na(True_beta_use), 1, all), , drop = FALSE]
            
            bhat_dev_save <- matrix(NA_real_, nrow = dim(True_beta_temp)[1], ncol = dim(True_beta_temp)[2])
            bhat_dev_save[1:nrow(bhat), ] <- (bhat - True_beta_use)^2
            
            Crossfit <- cbind(1, Xsmall) %*% bhat  # n x Q
            Cross_rate <- rowMeans(apply(Crossfit, 1, function(r) sort(r) != r))
            
            y_OOS <- Forecast_y[, case_num]
            xfit <- Forecast_x[, , case_num, drop = TRUE]
            xfit <- xfit[, colSums(is.na(xfit)) == 0, drop=FALSE]
            yfit <- cbind(1, xfit) %*% bhat        # 100 x Q
            
            wQS_temp <- rep(NA_real_, 4)
            for (w in 1:4) {
              wQS_temp[w] <- mean(wQS(y_OOS, yfit, tau, w)$output, na.rm = TRUE)
            }
            
            if (bhat_it == 1) {
              bhat_ALL_temp_UNSPARS[, , case_num] <- bhat_save
              bhat_deviation_temp_UNSPARS[, , case_num] <- bhat_dev_save
              wQS_temp_UNSPARS[, case_num] <- wQS_temp
              Cross_temp_UNSPARS[, case_num] <- Cross_rate
            } else {
              bhat_ALL_temp_SAVS[, , case_num] <- bhat_save
              bhat_deviation_temp_SAVS[, , case_num] <- bhat_dev_save
              wQS_temp_SAVS[, case_num] <- wQS_temp
              Cross_temp_SAVS[, case_num] <- Cross_rate
            }
          }
        }
        
        # Store UNSPARS
        bhat_ALL[, , , Prior_i] <- bhat_ALL_temp_UNSPARS
        bhat_deviation[, , , Prior_i] <- bhat_deviation_temp_UNSPARS
        Rhat[, , , Prior_i] <- Rhat_temp
        Neff[, , , Prior_i] <- Neff_temp
        wQS_ALL[, , Prior_i] <- wQS_temp_UNSPARS
        Cross_ALL[, , Prior_i] <- Cross_temp_UNSPARS
        
        # Store SAVS next slot if applicable
        if (Prior_choice %in% c("LQVP_HS", "LQVP_HS_ASIS")) {
          
            bhat_ALL[, , , Prior_i] <- bhat_ALL_temp_SAVS
            bhat_deviation[, , , Prior_i ] <- bhat_deviation_temp_SAVS
            Rhat[, , , Prior_i ] <- Rhat_temp
            Neff[, , , Prior_i ] <- Neff_temp
            wQS_ALL[, , Prior_i ] <- wQS_temp_SAVS
            Cross_ALL[, , Prior_i ] <- Cross_temp_SAVS
          
        }
      }
    }
    
    results[[MC_it]] <- list(
      wQS_ALL = wQS_ALL,
      Cross_ALL = Cross_ALL,
      bhat_deviation = bhat_deviation,
      bhat_ALL = bhat_ALL,
      Rhat = Rhat,
      Neff = Neff,
      Prior_mat = Prior_mat,
      tau = tau,
      Time_mat = Time_mat,
      MC_T = MC_T,
      rho = rho
    )

    # Save output (not sure whether we want to save by the iteration again or in one result list)
    file <- sprintf("ZZZ_Results_it=%i_T=%i.RData", MC_it, MC_T)
    save(
      wQS_ALL, Cross_ALL, bhat_deviation, bhat_ALL, Rhat, Neff,
      Prior_mat, tau, Time_mat, MC_T, rho,
      file = file
    )
  }
  
 