#' Calculate the quantile loss based on the given samples, coefficient, and the
#' quantile level
#' @export
#' @param x Input matrix. The convariates of samples.
#' @param y Input vector. The responses of samples.
#' @param beta Input vector. The coefficient.
#' @param u Input scalar. The quantile level.
#' @return A scalar. The quantile loss based on the given inputs.
cdloss <- function(x, y, beta, u) {
  return(mean((y - x %*% beta) * (u - ifelse(y - x %*% beta <= 0, 1, 0))))
}


#' Data generation function to implement the simulation settings in our paper
#' @export
#' @param p Input integer. Default is 150. The dimensionality for the
#' covariate.
#' @param n Input integer. Default is 200. The sample size for each datasets,
#' including all the target and sources.
#' @param s Input integer. Default is 5. The sparsity level in the simulation.
#' @param d Input scalar. Default is 2. Control the scale of the Laplace
#' distribution, larger d leads to smaller similarity between the target and
#' informative sets. Details are in the section for simulation.
#' @param An Input integer. Default is 5. The number of informative set.
#' @param M Input integer. Default is 10. The total number of sources.
#' @param eta Input scalar. Default is 20. The signal-to-noise ratio for
#' residual.
#' @param cov_type Input string. Default is 'auto'. Control the covariance
#' matrix of datasets. The value 'auto' leads to a homogeneous design with
#' auto-covariance; the value 'toeplitz' leads to heterogeneous designs with
#' Toeplitz covariance matrix. No other value is permitted.
#' @param res_type Input string. Default is 'normal'. Control the type of
#' distribution for residual. Only one from 'normal' and 'cauchy' is allowed.
#' @param seed Input integer. The default is 111.  The random seed.
#' @importFrom stats toeplitz rnorm rcauchy qcauchy qnorm
#' @importFrom MASS mvrnorm
#' @importFrom ExtDist rLaplace
#' @return NULL. But will generate correspoding pairs of (X,y) in global
#' enviroment.
data_generation <- function(p = 150, n = 200, s = 20, d = 2, An = 5, M = 10,
                            eta = 20, cov_type = 'auto', res_type = 'normal',
                            seed = 111){
  set.seed(seed)
  # Construct beta
  beta_0 <- rep(0, p + 1)
  beta_0[2:(s + 1)] <- 1
  assign('beta_0', beta_0, envir = .GlobalEnv)
  A_index <- c(rep(1, An), rep(0, M - An))
  for (i in 1:length(A_index)) {
    if (A_index[i] == 1) {
      Hk <- sample(seq(from = 1, to = p, by = 1), size = p / 2,
                   replace = FALSE)
      ksi <- rLaplace(n = length(Hk), mu = 0, b = 2 * d / p)
      ksi_add <- rep(0, p)
      ksi_add[Hk] <- ksi
      nam <- paste0("beta_", i)
      assign(nam, beta_0 + c(0, ksi_add), envir = .GlobalEnv)
    } else{
      Hk <- sample(seq(from = 1, to = p, by = 1), size = p / 2,
                   replace = FALSE)
      ksi <- rLaplace(n = length(Hk), mu = 0, b = 2 * 70 / p)
      ksi_add <- rep(0, p)
      ksi_add[Hk] <- ksi
      nam <- paste0("beta_", i)
      assign(nam, beta_0 + c(0, ksi_add), envir = .GlobalEnv)
    }
  }
  # Construct covariates:
  if (cov_type == 'auto') {
    sigma <- matrix(rep(0, p * p), nrow = p, byrow = TRUE)
    for (i in 1:p) {
      for (j in 1:i) {
        sigma[i, j] <- 0.5 ^ abs(i - j)
        sigma[j, i] <- sigma[i, j]
      }
    }
    for (i in 0:M) {
      nam <- paste0("X_", i)
      if (i == 0) {
        assign(nam, cbind(rep(1, floor(n / 0.8) + 1),
                          mvrnorm(n = floor(n / 0.8) + 1, rep(0, p), sigma,
                          tol = 1e-10))
               ,envir = .GlobalEnv)
      } else {
        assign(nam, cbind(rep(1, n),
                          mvrnorm(n = n, rep(0, p), sigma, tol = 1e-10)),
               envir = .GlobalEnv)
      }
      assign(paste0("sigma_", i), sigma)
    }
  } else if (cov_type == 'toeplitz') {
    for (i in 0:M) {
      if (i != 0) {
        row1 <- c(1, rep(1 / (i + 1), 2 * i - 1), rep(0, p - 2 * i))
        sigma <- toeplitz(row1)
      } else{
        sigma <- diag(rep(1, p))
      }
      nam <- paste0("X_", i)
      if (i == 0) {
        assign(nam, cbind(rep(1, floor(n / 0.8) + 1),
                          mvrnorm(n = floor(n / 0.8) + 1, rep(0, p), sigma,
                                  tol = 1e-10))
               ,envir = .GlobalEnv)
      } else{
        assign(nam, cbind(rep(1, n),
                          mvrnorm(n = n, rep(0, p), sigma, tol = 1e-10)),
               envir = .GlobalEnv)
      }
      assign(paste0("sigma_", i), sigma)
    }
  } else{
    print('wrong covariance type!')
  }
  # Construct response:
  if (res_type == 'normal'){
    btild0 <- get(paste0('beta_', 0))
    bt0 <- btild0[2:length(btild0)]
    assign("y_0", get(paste0('X_', 0)) %*% btild0 +
          rnorm(n = floor(n / 0.8) + 1,
                qnorm(0.2, mean = 0, sd = t(bt0) %*%
                      get(paste0('sigma_', 0)) %*% bt0 / eta),
                sd = t(bt0) %*% get(paste0('sigma_', 0)) %*% bt0 / eta),
          envir = .GlobalEnv)
    for (i in 1:M){
      btildi <- get(paste0('beta_', i))
      bti <- btildi[2:length(btildi)]
      nam <- paste0("y_", i)
      assign(nam, get(paste0('X_', i)) %*% btildi +
            rnorm(n = n,
                  qnorm(0.2, mean = 0, sd = t(bti) %*%
                        get(paste0('sigma_', i)) %*% bti / eta),
                  sd = t(bti) %*% get(paste0('sigma_', i)) %*% bti / eta),
            envir = .GlobalEnv)
    }
  } else if (res_type == 'cauchy'){
    btild0 <- get(paste0('beta_', 0))
    bt0 <- btild0[2:length(btild0)]
    assign("y_0", get(paste0('X_', 0)) %*% btild0 +
           rcauchy(n = floor(n / 0.8) + 1,
                  qcauchy(0.2, location = 0,
                          scale =  t(bt0) %*%
                            get(paste0('sigma_', 0)) %*% bt0 / eta),
                  scale = t(bt0) %*%
                    get(paste0('sigma_', 0)) %*% bt0 / eta),
           envir = .GlobalEnv)
    for (i in 1:M) {
      btildi <- get(paste0('beta_', i))
      bti <- btildi[2:length(btildi)]
      nam <- paste0("y_", i)
      assign(nam, get(paste0('X_', i)) %*% btildi +
             rcauchy(n = n,
                     qcauchy(0.2, location = 0,
                     scale =  t(bti) %*%
                       get(paste0('sigma_', i)) %*% bti / eta),
                     scale = t(bti) %*%
                       get(paste0('sigma_', i)) %*% bti / eta),
             envir = .GlobalEnv)
    }
  } else {
    print('Error on residual setting!')
  }
}

#' Build the smoothing kernel for estimation of density
#' @param x Numeric variable
#' @return Numeric scalar
K <- function(x) {
  if (x <= -1) {
    return(0)
  } else if (x >= 1) {
    return(0)
  } else{
    return((-315 / 64) * x ^ 6 + (735 / 64) * x ^ 4 - (525 / 64) * x ^ 2 +
             (105 / 64))
  }
}


#' Select bandwidth h by cross validation
#' @export
#' @param x design matrix which is in size of n*p.
#' @param y response vector.
#' @param u A scalar for the quantile level.
#' @importFrom quantreg rq.fit.lasso
#' @importFrom glmnet cv.glmnet glmnet
hchoose <- function(x, y, u) {
  n <- nrow(x)
  beta_hat <- quantreg::rq.fit.lasso(x, y, tau = u)$coefficients
  shat <- sum(abs(beta_hat) > (max(abs(beta_hat[2:length(beta_hat)])) * 1e-10))
  clist <- exp(seq(-5,3,length.out=25))
  hlist <- sapply(clist, function(x) {
     sqrt(shat * log(n) / n) + x * shat ^ (3 / 2) * log(n) / n
  })
  qlist <- sapply(hlist, function(h) {
    f_hat <- mean(apply(((y - x %*% beta_hat) / h), 1, K)) / h
    if (f_hat != 0) {
      y_tild <- x %*% beta_hat - f_hat ^ -1 * (ifelse(y - x %*% beta_hat <= 0,
                                                      1, 0) - u)
      cvfit <- cv.glmnet(x = x[, 2:ncol(x)], y = y_tild, alpha = 1,
                         standardize = FALSE, maxit = 1e6)
      bfit <- glmnet(x[, 2:ncol(x)], y_tild, lambda = cvfit$lambda.min,
                     standardize = FALSE)
      b_hat <- rbind(bfit$a0, bfit$beta)
      return(cdloss(x = x, y = y, beta = b_hat, u = u))
    } else{
      return(1e8)
    }
  })
  h <- sqrt(shat * log(n) / n) + clist[which.min(qlist)] * shat ^ (3 / 2) *
    log(n) / n
  return(h)
}



#' Oracle Trans-Lasso QR with a given informative sources set
#' @export
#' @param x_target A n*p matrix as the covariate from the target population.
#' @param y_target A vector with length n as the response from the target
#' population.
#' @param x_aux_bd A list object that each element is a n*p matrix as the
#' covariate from informative sources.
#' @param y_aux_bd A list object that each element is a vector with length n
#' as the covariate from informative sources.
#' @param u_target A scalar in (0,1). The quantile level for model of target
#' population.
#' @param u_aux_bd A vector that contains the quantile level for informative
#' populations.
#' @param h A vector with length (K+1), where K is the number of auxiliary
#' datasets. It pre-determine the bandwidth. The default is NULL, if NULL,
#' cross validation will be carried out.
#' @param lambda1 A scalar for the penalty strength in fusion step. Default
#' is null, if null, a cross validation will be carried for searching.
#' @param lambda2 A scalar for the penalty strength in debiasing step. Default
#' is null, if null, lambda2 will be calculated based on lambda1 according to
#' the proportion of sample sizes as described in our paper.
#' @param parallel A logic variable, default is FALSE.
#' @param ncore The integer, the number of cores used for parallel computation.
#' @param seed An integer for random seed, default is 111.
#' @importFrom quantreg rq.fit.lasso
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach `%dopar%`
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @return A vector. The corresponding transfer learning estimator based on
#' the given auxiliary sources.
rq.transfer <- function(x_target, y_target, x_aux_bd, y_aux_bd, u_target,
                        u_aux_bd, h = NULL, lambda1 = NULL, lambda2 = NULL,
                        parallel = FALSE, ncore = 20, seed = 111) {
  for (k in 1:length(x_aux_bd)) {
    if (ncol(x_target) != ncol(x_aux_bd[[k]])) {
      print("Two numbers of features don't match!")
      break
    }
  }
  set.seed(seed)
  nt <- nrow(x_target)
  # Step 1:
  # For target
  beta_target_hat <- rq.fit.lasso(x_target, y_target,
                                  tau = u_target)$coefficients
  # CV for h selection if NULL
  if (is.null(h)) {
    h_target <- hchoose(x_target, y_target, u_target)
  } else{
    h_target <- h[1]
  }
  f_target_hat <- mean(apply(((y_target - x_target %*% beta_target_hat) /
                      h_target), 1, K)) / h_target
  y_target_tild <- x_target %*% beta_target_hat - f_target_hat ^ -1 *
    (ifelse(y_target - x_target %*% beta_target_hat <= 0, 1, 0) - u_target)
  # For source
  if (parallel == TRUE){
    i <- NULL # pass R-CMD check
    cl <- makeCluster(min(detectCores() - 3, ncore))
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = length(x_aux_bd), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    y_aux_tild_bd <- foreach(i = 1:length(x_aux_bd),.options.snow = opts,
                             .combine = 'c', .multicombine=T,
                             .packages = c('quantreg'),
                             .inorder = TRUE) %dopar% {
      set.seed(seed)
      x_source <- x_aux_bd[[i]]
      y_source <- y_aux_bd[[i]]
      u_source <- u_aux_bd[i]
      ns <- nrow(x_source)
      beta_aux_hat <- rq.fit.lasso(x_source, y_source,
                                   tau = u_source)$coefficients
      if (is.null(h)) {
        h_source <- hchoose(x_source, y_source, u_source)
      } else{
        h_source <- h[k + 1]
      }
      f_aux_hat <- mean(apply(((y_source - x_source %*% beta_aux_hat) /
                       h_source), 1, K)) / h_source
      list(x_source %*% beta_aux_hat - f_aux_hat ^ -1 *
          (ifelse(y_source - x_source %*% beta_aux_hat <= 0, 1, 0) -
          u_source))
    }
    close(pb)
    stopCluster(cl)
  } else {
    y_aux_tild_bd <- list()
    for (k in 1:length(x_aux_bd)) {
      set.seed(seed)
      x_source <- x_aux_bd[[k]]
      y_source <- y_aux_bd[[k]]
      u_source <- u_aux_bd[k]
      ns <- nrow(x_source)
      beta_aux_hat <- rq.fit.lasso(x_source, y_source,
                                   tau = u_source)$coefficients
      if (is.null(h)) {
        h_source <- hchoose(x_source, y_source, u_source)
      } else{
        h_source <- h[k + 1]
      }
      f_aux_hat <- mean(apply(((y_source - x_source %*% beta_aux_hat) /
                                 h_source), 1, K)) / h_source
      y_aux_tild_bd[[k]] <- x_source %*% beta_aux_hat - f_aux_hat ^ -1 *
        (ifelse(y_source - x_source %*% beta_aux_hat <= 0, 1, 0) - u_source)
    }
  }
  # Step 2:
  x_comb <- x_target
  for (k in 1:length(x_aux_bd)) {
    x_comb <- rbind(x_comb, x_aux_bd[[k]])
  }
  y_comb <- y_target_tild
  for (k in 1:length(x_aux_bd)) {
    y_comb <- c(y_comb, y_aux_tild_bd[[k]])
  }
  if (is.null(lambda1)) {
    lambda1 <- cv.glmnet(x_comb[, 2:ncol(x_comb)], y_comb, alpha = 1,
                         standardize = FALSE)$lambda.min
  }
  wfit <- glmnet(x_comb[, 2:ncol(x_comb)], y_comb, lambda = lambda1,
                  standardize = FALSE, alpha = 1)
  w_hat <- rbind(wfit$a0, wfit$beta)
  # Step 3:
  if (is.null(lambda2)) {
    lambda2 <- lambda1 * sqrt(nrow(x_comb) / nrow(x_target))
  }
  deltafit <- glmnet(x_target[, 2:ncol(x_target)],
                     y_target_tild - (x_target %*% w_hat),
                     lambda = lambda2, standardize = FALSE, alpha = 1)
  delta_hat <- rbind(deltafit$a0, deltafit$beta)
  return(w_hat + delta_hat)
}


#' Build fusion model between datasets without debias procedure
#' @export
#' @param x_target A n*p matrix as the covariate from the target population.
#' @param y_target A vector with length n as the response from the target
#' population.
#' @param x_aux_bd A list object that each element is a n*p matrix as the
#' covariate from informative sources.
#' @param y_aux_bd A list object that each element is a vector with length n
#' as the covariate from informative sources.
#' @param u_target A scalar in (0,1). The quantile level for model of target
#' population.
#' @param u_aux_bd A vector that contains the quantile level for informative
#' populations.
#' @param h A vector with length (K+1), where K is the number of auxiliary
#' datasets. It pre-determine the bandwidth. The default is NULL, if NULL,
#' cross validation will be carried out.
#' @param lambda1 A scalar for the penalty strength in fusion step. Default
#' is null, if null, a cross validation will be carried for searching.
#' @param parallel A logic variable, default is FALSE.
#' @param ncore The integer, the number of cores used for parallel computation.
#' @param seed An integer for random seed, default is 111.
#' @importFrom quantreg rq.fit.lasso
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach `%dopar%`
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @return A vector. The fusion estimator withou debias procedure based on the
#' given auxiliary sources.
rq.fusion <- function(x_target, y_target, x_aux_bd, y_aux_bd, u_target,
                      u_aux_bd, h = NULL, lambda1 = NULL, parallel = FALSE,
                      ncore = 20, seed = 111) {
  for (k in 1:length(x_aux_bd)) {
    if (ncol(x_target) != ncol(x_aux_bd[[k]])) {
      print("Two numbers of features don't match!")
      break
    }
  }
  set.seed(seed)
  nt <- nrow(x_target)
  # Step 1:
  # For target
  beta_target_hat <- rq.fit.lasso(x_target, y_target,
                                  tau = u_target)$coefficients
  if (is.null(h)) {
    h_target <- hchoose(x_target, y_target, u_target)
  } else{
    h_target <- h[1]
  }
  f_target_hat <- mean(apply(((y_target - x_target %*% beta_target_hat) /
                      h_target), 1, K)) / h_target
  y_target_tild <- x_target %*% beta_target_hat - f_target_hat ^ -1 *
    (ifelse(y_target - x_target %*% beta_target_hat <= 0, 1, 0) - u_target)
  # For source
  if (parallel == TRUE){
    i <- NULL # pass R-CMD check
    cl <- makeCluster(min(detectCores() - 3, ncore))
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = length(x_aux_bd), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    y_aux_tild_bd <- foreach(i = 1:length(x_aux_bd),.options.snow = opts,
                             .combine = 'c', .multicombine=T,
                             .packages = c('quantreg'),
                             .inorder = TRUE) %dopar% {
      set.seed(seed)
      x_source <- x_aux_bd[[i]]
      y_source <- y_aux_bd[[i]]
      u_source <- u_aux_bd[i]
      ns <- nrow(x_source)
      beta_aux_hat <- rq.fit.lasso(x_source, y_source,
                                   tau = u_source)$coefficients
      if (is.null(h)) {
        h_source <- hchoose(x_source, y_source, u_source)
      } else{
        h_source <- h[k + 1]
      }
      f_aux_hat <- mean(apply(((y_source - x_source %*% beta_aux_hat) /
                       h_source), 1, K)) / h_source
      list(x_source %*% beta_aux_hat - f_aux_hat ^ -1 *
          (ifelse(y_source - x_source %*% beta_aux_hat <= 0, 1, 0) - u_source))
    }
    close(pb)
    stopCluster(cl)
  } else {
    y_aux_tild_bd <- list()
    for (k in 1:length(x_aux_bd)) {
      set.seed(seed)
      x_source <- x_aux_bd[[k]]
      y_source <- y_aux_bd[[k]]
      u_source <- u_aux_bd[k]
      ns <- nrow(x_source)
      beta_aux_hat <- rq.fit.lasso(x_source, y_source,
                                   tau = u_source)$coefficients
      if (is.null(h)) {
        h_source <- hchoose(x_source, y_source, u_source)
      } else{
        h_source <- h[k + 1]
      }
      f_aux_hat <- mean(apply(((y_source - x_source %*% beta_aux_hat) /
                       h_source), 1, K)) / h_source
      y_aux_tild_bd[[k]] <- x_source %*% beta_aux_hat - f_aux_hat ^ -1 *
        (ifelse(y_source - x_source %*% beta_aux_hat <= 0, 1, 0) - u_source)
    }
  }
  # Step 2:
  x_comb <- x_target
  for (k in 1:length(x_aux_bd)) {
    x_comb <- rbind(x_comb, x_aux_bd[[k]])
  }
  y_comb <- y_target_tild
  for (k in 1:length(x_aux_bd)) {
    y_comb <- c(y_comb, y_aux_tild_bd[[k]])
  }
  if (is.null(lambda1)) {
    lambda1 <- cv.glmnet(x_comb[, 2:ncol(x_comb)], y_comb, alpha = 1,
                         standardize = FALSE)$lambda.min
  }
  wfit <- glmnet(x_comb[, 2:ncol(x_comb)], y_comb, lambda = lambda1,
                  standardize = FALSE, alpha = 1)
  w_hat <- rbind(wfit$a0, wfit$beta)
  return(w_hat)
}

#' Build the auxiliary function to compute quantile loss for informative
#' sources detection
#' @param betahat Estimator for target only from training part.
#' @param betaic Estimator for target only from testing part.
#' @param X_measure A n*p matrix for covariate of population that is waiting
#' for the test.
#' @param y_measure A vector as the response from the population that is
#' waiting for the test.
#' @param u_target A scalar in (0,1). The quantile level for informative
#' source detection testing.
#' @param hic A scalar >0. The smoothing bandwidth for calculate the loss. If
#' not provided, then the cross validation will be carried out.
#' @return A scalar. The quantile loss for given estimator on the given
#' testing samples.
Q_loss <- function(betahat, betaic, X_measure, y_measure, u_target,
                   hic = NULL) {
  if (is.null(hic)) {
    hic <- hchoose(X_measure, y_measure, u_target)
  }
  f_hat <- mean(apply((y_measure - X_measure %*% betaic / hic), 1, K)) / hic
  y_measure_tilde <- X_measure %*% betaic - f_hat ^ (-1) *
    (ifelse(y_measure - X_measure %*% betaic <= 0, 1, 0) - u_target)
  return(mean((y_measure_tilde - X_measure %*% betahat) ^ 2))
}

#' Main function for informative sources detection among all dataset
#' @export
#' @param x_target A n*p matrix as the covariate from the target population.
#' @param y_target A vector with length n as the response from the target
#' population.
#' @param x_aux_bd A list object that each element is a n*p matrix as the
#' covariate from informative sources.
#' @param y_aux_bd A list object that each element is a vector with length n
#' as the covariate from informative sources.
#' @param u_target A scalar in (0,1). The quantile level for model of target
#' population.
#' @param u_aux_bd A vector that contains the quantile level for informative
#' populations.
#' @param epsilon A scalar. Default is 0.01. The strict level for informative
#' sources detectiond, larger the value is, less strict the procedure is.
#' @param h A vector with length (K+1), where K is the number of auxiliary
#' datasets. It pre-determine the bandwidth. The default is NULL, if NULL,
#' cross validation will be carried out.
#' @param lambda1 A scalar for the penalty strength in fusion step. Default
#' is null, if null, a cross validation will be carried for searching.
#' @param parallel A logic variable, default is TRUE
#' @param ncore The integer, the number of cores used for parallel computation.
#' @param info_num An integar. Default is NULL. The given number of informative
#' sources under pseduo running.
#' @param verbose A logic variable, default is FALSE.
#' @param seed A integer variable for the random seed, default is 111.
#' @importFrom quantreg rq.fit.lasso
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom foreach foreach `%dopar%`
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils head setTxtProgressBar txtProgressBar
#' @return A vector of index. The index of the informative sources.
info_detect <- function(x_target, y_target, x_aux_bd, y_aux_bd, u_target,
                        u_aux_bd, epsilon = 0.01, h = NULL, lambda1 = NULL,
                        parallel = TRUE, ncore = 10, info_num = NULL,
                        seed = 111, verbose = FALSE){
  if (length(x_aux_bd) != length(y_aux_bd)) {
    print('the # of datasets for x and y are not agreed!')
  } else {
    set.seed(seed)
    M <- length(x_aux_bd)
    total_index <- 1:length(x_aux_bd)
    # Step 1.1: cut the target half and half: (remark: 1 for I, 2 for Ic)
    target_index <- 1:nrow(x_target)
    X0_index <- sample(target_index, size = round(nrow(x_target) / 2),
                       replace = FALSE)
    X0_cut1 <- x_target[X0_index,]
    X0_cut2 <- x_target[-X0_index,]
    y0_cut1 <- y_target[X0_index]
    y0_cut2 <- y_target[-X0_index]
    # CV1:
    # Step 1.2: train sparse quantile regression for I for X0, y0:
    beta_0_hat <- rq.fit.lasso(X0_cut1, y0_cut1, tau = u_target)$coefficients
    # Step 1.3: train sparse quantile regression for I for X0, y0 + k-th
    # source:
    if (parallel == FALSE){
      for (k in 1:length(x_aux_bd)) {
        if (is.null(h)) {
          assign(paste0(paste0('beta_', k), '_hat'),
                 rq.fusion(X0_cut1, y0_cut1, list(x_aux_bd[[k]]),
                           list(y_aux_bd[[k]]), u_target, c(u_aux_bd[k]),
                           lambda1 = lambda1))
        } else{
          assign(paste0(paste0('beta_', k), '_hat'),
                 rq.fusion(X0_cut1, y0_cut1, list(x_aux_bd[[k]]),
                           list(y_aux_bd[[k]]), u_target, c(u_aux_bd[k]),
                           h = c(h[1], h[k + 1]), lambda1 = lambda1))
        }
      }
    } else{
      i <- NULL # pass R-CMD check
      cl <- makeCluster(min(detectCores() - 3, ncore))
      registerDoSNOW(cl)
      pb <- txtProgressBar(max = M, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      betasources <- foreach(i = 1:M, .options.snow = opts, .combine = cbind,
                             .packages = c('quantreg', 'glmnet'),
                             .inorder = TRUE) %dopar% {
        if (is.null(h)) {
          rq.fusion(X0_cut1, y0_cut1, list(x_aux_bd[[i]]), list(y_aux_bd[[i]]),
                    u_target, c(u_aux_bd[i]), lambda1 = lambda1)
        } else{
          rq.fusion(X0_cut1, y0_cut1, list(x_aux_bd[[i]]), list(y_aux_bd[[i]]),
                    u_target, c(u_aux_bd[i]), h = c(h[1], h[i + 1]),
                    lambda1 = lambda1)
        }
      }
      close(pb)
      stopCluster(cl)
      for (k in 1:length(x_aux_bd)) {
        assign(paste0(paste0('beta_', k), '_hat'), betasources[, k])
      }
    }
    # Step 2: Calculate the loss based on the loss function:
    beta_0_hatic <- rq.fit.lasso(X0_cut2, y0_cut2, tau = u_target)$coefficients
    #hic <- hchoose(X0_cut2, y0_cut2, u_target)
    hic <- h[1]
    for (k in c(0, (1:length(x_aux_bd)))) {
      assign(paste0('Q', k),
             Q_loss(get(paste0(paste0('beta_', k), '_hat')), beta_0_hatic,
                    X0_cut2, y0_cut2, u_target, hic = hic))
    }
    # Step 3: Get informative sets:
    # Step 3.1 true
    index_true1 <- c()
    if (verbose == TRUE) {
      print(get(paste0('Q', 0)))
    }
    for (k in 1:length(x_aux_bd)) {
      if (verbose == TRUE) {
        print(get(paste0('Q', k)))
      }
      if (get(paste0('Q', k)) <= (1 + epsilon) * get(paste0('Q', 0))) {
        index_true1 <- c(index_true1, k)
      } else{
        next
      }
    }
    # Step 3.2 psd
    index_psd1 <- NULL
    if (!is.null(info_num)) {
      Q_K_vec <- c()
      for (k in 1:length(x_aux_bd)) {
        Q_K_vec <- c(Q_K_vec, get(paste0('Q', k)))
      }
      index_psd1 <- match(head(sort(Q_K_vec, decreasing = FALSE), info_num),
                          Q_K_vec)
    } else{
      index_psd1 <- c()
    }
    # CV2
    # Step 1.2: train sparse quantile regression for I for X0, y0:
    beta_0_hat <- rq.fit.lasso(X0_cut2, y0_cut2, tau = u_target)$coefficients
    # Step 1.3: train sparse quantile regression for I for X0, y0 + k-th
    # source:
    if (parallel == FALSE){
      for (k in 1:length(x_aux_bd)) {
        if (is.null(h)) {
          assign(paste0(paste0('beta_', k), '_hat'),
                 rq.fusion(X0_cut2, y0_cut2, list(x_aux_bd[[k]]),
                           list(y_aux_bd[[k]]), u_target, c(u_aux_bd[k]),
                           lambda1 = lambda1))
        } else{
          assign(paste0(paste0('beta_', k), '_hat'),
                 rq.fusion(X0_cut2, y0_cut2, list(x_aux_bd[[k]]),
                           list(y_aux_bd[[k]]), u_target, c(u_aux_bd[k]),
                           h = c(h[1], h[k + 1]), lambda1 = lambda1))
        }
      }
    } else{
      i <- NULL # pass R-CMD check
      cl <- makeCluster(min(detectCores() - 3, ncore))
      registerDoSNOW(cl)
      pb <- txtProgressBar(max = M, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      betasources <- foreach(i = 1:M, .options.snow = opts, .combine = cbind,
                             .packages = c('quantreg', 'glmnet'),
                             .inorder = TRUE) %dopar% {
        if (is.null(h)) {
          rq.fusion(X0_cut2, y0_cut2, list(x_aux_bd[[i]]), list(y_aux_bd[[i]]),
                    u_target, c(u_aux_bd[i]), lambda1 = lambda1)
        } else {
          rq.fusion(X0_cut2, y0_cut2, list(x_aux_bd[[i]]), list(y_aux_bd[[i]]),
                    u_target, c(u_aux_bd[i]), h = c(h[1], h[i + 1]),
                    lambda1 = lambda1)
          }
      }
      close(pb)
      stopCluster(cl)
      for (k in 1:length(x_aux_bd)) {
        assign(paste0(paste0('beta_', k), '_hat'), betasources[, k])
      }
    }
    # Step 2: Calculate the loss based on the loss function:
    beta_0_hatic <- rq.fit.lasso(X0_cut1, y0_cut1, tau = u_target)$coefficients
    #hic <- hchoose(X0_cut1,y0_cut1,u_target)
    hic <- h[1]
    for (k in c(0, (1:length(x_aux_bd)))) {
      assign(paste0('Q', k),
             Q_loss(get(paste0(paste0('beta_', k), '_hat')),
                    beta_0_hatic, X0_cut1, y0_cut1, u_target, hic = hic))
    }
    # Step 3: Get informative sets:
    # Step 3.1 true
    index_true2 <- c()
    if (verbose == TRUE) {
      print(get(paste0('Q', 0)))
    }
    for (k in 1:length(x_aux_bd)) {
      if (verbose == TRUE) {
        print(get(paste0('Q', k)))
      }
      if (get(paste0('Q', k)) <= (1 + epsilon) * get(paste0('Q', 0))) {
        index_true2 <- c(index_true2, k)
      } else{
        next
      }
    }
    # Step 3.2 psd
    index_psd2 <- NULL
    if (!is.null(info_num)) {
      Q_K_vec <- c()
      for (k in 1:length(x_aux_bd)) {
        Q_K_vec <- c(Q_K_vec, get(paste0('Q', k)))
      }
      index_psd2 <- match(head(sort(Q_K_vec, decreasing = FALSE), info_num),
                          Q_K_vec)
    } else{
      index_psd2 <- c()
    }
    index_true <- intersect(index_true1, index_true2)
    index_psd <- intersect(index_psd1, index_psd2)
    return(list('true' = index_true, 'psd' = index_psd))
  }
}

#' Oracle Trans-pooling QR with a given informative sources set
#' @export
#' @param x_target A n*p matrix as the covariate from the target population.
#' @param y_target A vector with length n as the response from the target
#' population.
#' @param x_aux_bd A list object that each element is a n*p matrix as the
#' covariate from informative sources.
#' @param y_aux_bd A list object that each element is a vector with length n
#' as the covariate from informative sources.
#' @param u A scalar in (0,1). The quantile level for all models.
#' @importFrom quantreg rq.fit.lasso
#' @return A vector. The corresponding transfer learning estimator based on
#' the given auxiliary sources based on Huang 2022.
rq.transfer.pool <- function(x_target, y_target, x_aux_bd, y_aux_bd, u){
  for (k in 1:length(x_aux_bd)) {
    if (ncol(x_target) != ncol(x_aux_bd[[k]])) {
      print("Two numbers of features don't match!")
      break
    }
  }
  # fusion
  x_comb <- x_target
  for (k in 1:length(x_aux_bd)){
    x_comb <- rbind(x_comb, x_aux_bd[[k]])
  }
  y_comb <- y_target
  for (k in 1:length(x_aux_bd)){
    y_comb <- c(y_comb, y_aux_bd[[k]])
  }
  beta_fusion <- rq.fit.lasso(x_comb, y_comb, tau = u)$coefficients
  # debias
  beta_db <- rq.fit.lasso(x_target,
                          y_target - (x_target %*% beta_fusion),
                          tau = u)$coefficients
  return(beta_fusion + beta_db)
}

#' Build fusion model for Trans-pooling QR between datasets without debias
#' procedure
#' @export
#' @param x_target A n*p matrix as the covariate from the target population.
#' @param y_target A vector with length n as the response from the target
#' population.
#' @param x_aux_bd A list object that each element is a n*p matrix as the
#' covariate from informative sources.
#' @param y_aux_bd A list object that each element is a vector with length n
#' as the covariate from informative sources.
#' @param u A scalar in (0,1). The quantile level for all models.
#' @importFrom quantreg rq.fit.lasso
#' @return A vector. The fusion estimator withou debias procedure based on the
#' given auxiliary sources.
rq.fusion.pool <- function(x_target, y_target, x_aux_bd, y_aux_bd, u) {
  for (k in 1:length(x_aux_bd)) {
    if (ncol(x_target) != ncol(x_aux_bd[[k]])) {
      print("Two numbers of features don't match!")
      break
    }
  }
  # fusion
  x_comb <- x_target
  for (k in 1:length(x_aux_bd)){
    x_comb <- rbind(x_comb, x_aux_bd[[k]])
  }
  y_comb <- y_target
  for (k in 1:length(x_aux_bd)){
    y_comb <- c(y_comb, y_aux_bd[[k]])
  }
  beta_fusion <- rq.fit.lasso(x_comb, y_comb, tau = u)$coefficients
  return(beta_fusion)
}


#' Main function for informative sources detection among all dataset for
#' Trans-pooling QR
#' @export
#' @param x_target A n*p matrix as the covariate from the target population.
#' @param y_target A vector with length n as the response from the target
#' population.
#' @param x_aux_bd A list object that each element is a n*p matrix as the
#' covariate from informative sources.
#' @param y_aux_bd A list object that each element is a vector with length n
#' as the covariate from informative sources.
#' @param u A scalar in (0,1). The quantile level for all models.
#' @param epsilon A scalar. Default is 0.01. The strict level for informative
#' sources detectiond, larger the value is, less strict the procedure is.
#' @param parallel A logic variable, default is TRUE
#' @param ncore The integer, the number of cores used for parallel computation.
#' @param info_num An integar. Default is NULL. The given number of informative
#' sources under pseduo running.
#' @param verbose A logic variable, default is FALSE.
#' @param seed A integer variable for the random seed, default is 111.
#' @importFrom quantreg rq.fit.lasso
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom foreach foreach `%dopar%`
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils head setTxtProgressBar txtProgressBar
#' @return A vector of index. The index of the informative sources.
info_detect.pool <- function(x_target, y_target, x_aux_bd, y_aux_bd, u,
                        epsilon = 0.01, parallel = TRUE, ncore = 10,
                        info_num = NULL, seed = 111, verbose = FALSE){
  if (length(x_aux_bd) != length(y_aux_bd)) {
    print('the # of datasets for x and y are not agreed!')
  } else {
    set.seed(seed)
    M <- length(x_aux_bd)
    total_index <- 1:length(x_aux_bd)
    # Step 1.1: cut the target half and half: (remark: 1 for I, 2 for Ic)
    target_index <- 1:nrow(x_target)
    X0_index <- sample(target_index, size = round(nrow(x_target) / 2),
                       replace = FALSE)
    X0_cut1 <- x_target[X0_index,]
    X0_cut2 <- x_target[-X0_index,]
    y0_cut1 <- y_target[X0_index]
    y0_cut2 <- y_target[-X0_index]
    # CV1:
    # Step 1.2: train sparse quantile regression for I for X0, y0:
    beta_0_hat <- rq.fit.lasso(X0_cut1, y0_cut1, tau=u)$coefficients
    # Step 1.3: train sparse quantile regression for I for X0, y0 + k-th
    # source:
    if (parallel == FALSE){
      for (k in 1:length(x_aux_bd)) {
        assign(paste0(paste0('beta_', k), '_hat'),
               rq.fusion.pool(X0_cut1, y0_cut1, list(x_aux_bd[[k]]),
                              list(y_aux_bd[[k]]), u))
      }
    } else{
      i <- NULL # pass R-CMD check
      cl <- makeCluster(min(detectCores() - 3, ncore))
      registerDoSNOW(cl)
      pb <- txtProgressBar(max = M, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      betasources <- foreach(i = 1:M, .options.snow = opts, .combine = cbind,
                             .packages = c('quantreg'),
                             .inorder = TRUE) %dopar% {
        rq.fusion.pool(X0_cut1, y0_cut1, list(x_aux_bd[[i]]),
                       list(y_aux_bd[[i]]),u)
      }
      close(pb)
      stopCluster(cl)
      for (k in 1:length(x_aux_bd)) {
        assign(paste0(paste0('beta_', k), '_hat'), betasources[, k])
      }
    }
    # Step 2: Calculate the loss based on the loss function:
    for (k in c(0, (1:length(x_aux_bd)))) {
      assign(paste0('Q', k),
             cdloss(X0_cut2, y0_cut2, get(paste0(paste0('beta_', k), '_hat')),
                    u))
    }
    # Step 3: Get informative sets:
    # Step 3.1 true
    index_true1 <- c()
    if (verbose == TRUE) {
      print(get(paste0('Q', 0)))
    }
    for (k in 1:length(x_aux_bd)) {
      if (verbose == TRUE) {
        print(get(paste0('Q', k)))
      }
      if (get(paste0('Q', k)) <= (1 + epsilon) * get(paste0('Q', 0))) {
        index_true1 <- c(index_true1, k)
      } else{
        next
      }
    }
    # Step 3.2 psd
    index_psd1 <- NULL
    if (!is.null(info_num)) {
      Q_K_vec <- c()
      for (k in 1:length(x_aux_bd)) {
        Q_K_vec <- c(Q_K_vec, get(paste0('Q', k)))
      }
      index_psd1 <- match(head(sort(Q_K_vec, decreasing = FALSE), info_num),
                          Q_K_vec)
    } else{
      index_psd1 <- c()
    }
    # CV2
    # Step 1.2: train sparse quantile regression for I for X0, y0:
    beta_0_hat <- rq.fit.lasso(X0_cut2, y0_cut2, tau = u)$coefficients
    # Step 1.3: train sparse quantile regression for I for X0, y0 + k-th
    if (parallel == FALSE){
      for (k in 1:length(x_aux_bd)) {
        assign(paste0(paste0('beta_', k), '_hat'),
               rq.fusion.pool(X0_cut2, y0_cut2, list(x_aux_bd[[k]]),
                              list(y_aux_bd[[k]]), u))
      }
    } else{
      i <- NULL # pass R-CMD check
      cl <- makeCluster(min(detectCores() - 3, ncore))
      registerDoSNOW(cl)
      pb <- txtProgressBar(max = M, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      betasources <- foreach(i = 1:M, .options.snow = opts, .combine = cbind,
                             .packages = c('quantreg'),
                             .inorder = TRUE) %dopar% {
        rq.fusion.pool(X0_cut2, y0_cut2, list(x_aux_bd[[i]]),
                       list(y_aux_bd[[i]]),u)
      }
      close(pb)
      stopCluster(cl)
      for (k in 1:length(x_aux_bd)) {
        assign(paste0(paste0('beta_', k), '_hat'), betasources[, k])
      }
    }
    # Step 2: Calculate the loss based on the loss function:
    for (k in c(0, (1:length(x_aux_bd)))) {
      assign(paste0('Q', k),
             cdloss(X0_cut1, y0_cut1, get(paste0(paste0('beta_', k), '_hat')),
                    u))
    }
    # Step 3: Get informative sets:
    # Step 3.1 true
    index_true2 <- c()
    if (verbose == TRUE) {
      print(get(paste0('Q', 0)))
    }
    for (k in 1:length(x_aux_bd)) {
      if (verbose == TRUE) {
        print(get(paste0('Q', k)))
      }
      if (get(paste0('Q', k)) <= (1 + epsilon) * get(paste0('Q', 0))) {
        index_true2 <- c(index_true2, k)
      } else{
        next
      }
    }
    # Step 3.2 psd
    index_psd2 <- NULL
    if (!is.null(info_num)) {
      Q_K_vec <- c()
      for (k in 1:length(x_aux_bd)) {
        Q_K_vec <- c(Q_K_vec, get(paste0('Q', k)))
      }
      index_psd2 <- match(head(sort(Q_K_vec, decreasing = FALSE), info_num),
                          Q_K_vec)
    } else{
      index_psd2 <- c()
    }
    index_true <- intersect(index_true1, index_true2)
    index_psd <- intersect(index_psd1, index_psd2)
    return(list('true' = index_true, 'psd' = index_psd))
  }
}

