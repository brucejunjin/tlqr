#' Calculate the quantile loss based on the given samples, coefficient, and the
#' quantile level
#' @export
#' @param x Input design matrix. The convariates of samples.
#' @param y Input response vector. The responses of samples.
#' @param beta Input vector. The coefficient.
#' @param u Input scalar. The quantile level.
#' @return A scalar. The quantile loss based on the given inputs.
QuantileLoss <- function(x, y, beta, u) {
  xtilde <- cbind(rep(1, nrow(x)), x)
  return(mean((y - xtilde %*% beta) *
              (u - ifelse(y - xtilde %*% beta <= 0, 1, 0))))
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
#' @param seed An integer for random seed, default is 111.
#' @importFrom stats toeplitz rnorm rcauchy qcauchy qnorm
#' @importFrom MASS mvrnorm
#' @importFrom ExtDist rLaplace
#' @return NULL. But will generate correspoding pairs of (X,y) in global
#' enviroment.
DataDemo <- function(p = 150, n = 200, s = 20, d = 2, An = 5, M = 10,
                     eta = 20, cov_type = 'auto', res_type = 'normal',
                     seed = 111)
{
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
        assign(nam, mvrnorm(n = floor(n / 0.8) + 1, rep(0, p), sigma,
                                  tol = 1e-10)
               ,envir = .GlobalEnv)
      } else {
        assign(nam, mvrnorm(n = n, rep(0, p), sigma, tol = 1e-10),
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
        assign(nam, mvrnorm(n = floor(n / 0.8) + 1, rep(0, p), sigma,
                                  tol = 1e-10)
               ,envir = .GlobalEnv)
      } else{
        assign(nam, mvrnorm(n = n, rep(0, p), sigma, tol = 1e-10),
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
    assign("y_0", get(paste0('X_', 0)) %*% bt0 +
             rnorm(n = floor(n / 0.8) + 1,
                   qnorm(0.2, mean = 0, sd = t(bt0) %*%
                           get(paste0('sigma_', 0)) %*% bt0 / eta),
                   sd = t(bt0) %*% get(paste0('sigma_', 0)) %*% bt0 / eta),
           envir = .GlobalEnv)
    for (i in 1:M){
      btildi <- get(paste0('beta_', i))
      bti <- btildi[2:length(btildi)]
      nam <- paste0("y_", i)
      assign(nam, get(paste0('X_', i)) %*% bti +
               rnorm(n = n,
                     qnorm(0.2, mean = 0, sd = t(bti) %*%
                             get(paste0('sigma_', i)) %*% bti / eta),
                     sd = t(bti) %*% get(paste0('sigma_', i)) %*% bti / eta),
             envir = .GlobalEnv)
    }
  } else if (res_type == 'cauchy'){
    btild0 <- get(paste0('beta_', 0))
    bt0 <- btild0[2:length(btild0)]
    assign("y_0", get(paste0('X_', 0)) %*% bt0 +
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
      assign(nam, get(paste0('X_', i)) %*% bti +
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
#' @param control The list object for control, see \code{tlqr.control}.
#' @importFrom quantreg rq.fit.lasso
#' @importFrom glmnet cv.glmnet glmnet
ChooseBandwidth <- function(x, y, u, control = list()) {
  internal.params <- do.call("tlqr.control", control)
  lowc <- internal.params$lowc
  upc <- internal.params$upc
  thresh <- internal.params$thresh
  length.out <- internal.params$length.out
  max.it <- internal.params$max.it
  max.round <- internal.params$max.round
  iter.round <- 1
  n <- nrow(x)
  xtilde <- cbind(rep(1, n), x)
  beta_hat <- quantreg::rq.fit.lasso(xtilde, y, tau = u)$coefficients
  shat <- sum(abs(beta_hat) > (max(abs(beta_hat[2:length(beta_hat)]))
                               * thresh))
  if (shat >= 0.15 * n){
    length.out = 2 * length.out
  }
  lowb <- lowc * (sqrt(shat * log(n) / n) + shat ^ (3 / 2) * log(n) / n)
  upb <- upc * (sqrt(shat * log(n) / n) + shat ^ (3 / 2) * log(n) / n)
  hlist <- seq(from = lowb, to = upb, length.out = length.out)
  qlist <- sapply(hlist, function(h) {
    f_hat <- mean(apply(((y - xtilde %*% beta_hat) / h), 1, K)) / h
    if (f_hat != 0) {
      y_tild <- xtilde %*% beta_hat - f_hat ^ -1 *
        (ifelse(y - xtilde %*% beta_hat <= 0, 1, 0) - u)
      cvfit <- cv.glmnet(x = x, y = y_tild, alpha = 1,
                         standardize = FALSE, maxit = max.it)
      bfit <- glmnet(x = x, y = y_tild, lambda = cvfit$lambda.min,
                     standardize = FALSE)
      b_hat <- rbind(bfit$a0, bfit$beta)
      return(QuantileLoss(x = x, y = y, beta = as.vector(b_hat), u = u))
    } else{
      return(1e8)
    }
  })
  if (which.min(qlist) == 1) {
    if (iter.round <= max.round) {
      upc <- hlist[2] / (sqrt(shat * log(n) / n) + shat ^ (3 / 2) * log(n) / n)
      return(ChooseBandwidth(x, y, u,
                             control = list(lowc = 0.5 * lowc, upc = upc)))
    } else{
      warning('The minimum bandwidth in range is chosen!')
    }
  } else if (which.min(qlist) == length.out) {
    if (iter.round <= max.round) {
      lowc <- hlist[length.out - 1] /
        (sqrt(shat * log(n) / n) + shat ^ (3 / 2) * log(n) / n)
      return(ChooseBandwidth(x, y, u,
                             control = list(lowc = lowc, upc = 1.5 * upc)))
    } else{
      warning('The maximum bandwidth in range is chosen!')
    }
  }
  h <- hlist[which.min(qlist)]
  return(h)
}


#' Internal tlqr parameters
#' @param h A vector with length (K+1), where K is the number of auxiliary
#' datasets. It pre-determine the bandwidth. The default is \code{NULL},
#' if \code{NULL}, cross validation will be carried out.
#' @param hic A scalar >0. The smoothing bandwidth for calculate the loss.
#' Default is \code{NULL}. If \code{NULL}, then the cross validation will be
#' carried out.
#' @param parallel A logic variable to determine whether the Step 1 in our
#' proposed method run with parallel computation, default is \code{FALSE}.
#' @param ncore The integer, the number of cores used for parallel computation.
#' @param seed An integer for random seed for data split in infotmative set
#' detection procedure, default is 111.
#' @param lowc Constant used for calculating lower bound of bandwidth according
#' to our paper, default is 0.05.
#' @param upc Constant used for calculating upper bound of bandwidth according
#' to our paper, default is 1.
#' @param thresh The threshold used to determine the sparsity of initial,
#' default is 1e-10.
#' @param length.out The number of grid search points for searching the best
#' bandwidth, default is 25.
#' @param max.it The maximum iteration number for lasso problem via glmnet,
#' default is 1e6.
#' @param max.round The maximum round of expanding the range of bandwidth
#' selection before warning if provided is not suitable. Default is 5.
#' @param verbose A logic variable to control the print of loss during
#' informative set detection, default is FALSE.
#' @export
tlqr.control <- function(h = NULL, hic = NULL, parallel = FALSE, ncore = 20,
                         seed = 111, lowc = 0.05, upc = 1, thresh = 1e-10,
                         length.out = 25, max.it = 1e6, max.round = 5,
                         verbose = FALSE) {
  value <- list('h' = h, 'hic' = hic, 'parallel' = parallel, 'ncore' = ncore,
                'seed' = seed, 'lowc' = lowc, 'upc' = upc, 'thresh' = thresh,
                'length.out' = length.out, 'max.it' = max.it,
                'max.round' = max.round, 'verbose' = verbose)
  value
}



#' Oracle Trans-Lasso QR with a given informative sources set
#' @export
#' @param x.target A n*p matrix as the covariate from the target population.
#' @param y.target A vector with length n as the response from the target
#' population.
#' @param x.source A list object that each element is a n*p matrix as the
#' covariate from informative sources.
#' @param y.source A list object that each element is a vector with length n
#' as the covariate from informative sources.
#' @param u.target A scalar in (0,1). The quantile level for model of target
#' population.
#' @param u.source A vector that contains the quantile level for informative
#' populations.
#' @param lambda1 A scalar for the penalty strength in fusion step. Default
#' is null, if null, a cross validation will be carried for searching.
#' @param lambda2 A scalar for the penalty strength in debiasing step. Default
#' is null, if null, lambda2 will be calculated based on lambda1 according to
#' the proportion of sample sizes as described in our paper.
#' @param control The list object for control, see \code{tlqr.control}.
#' @importFrom quantreg rq.fit.lasso
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach `%dopar%`
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @return A vector. The corresponding transfer learning estimator based on
#' the given auxiliary sources.
TransferQR <- function(x.target, y.target, x.source, y.source, u.target,
                       u.source, lambda1 = NULL, lambda2 = NULL,
                       control = list()) {
  for (k in 1:length(x.source)) {
    if (ncol(x.target) != ncol(x.source[[k]])) {
      print("Two numbers of features don't match!")
      break
    }
  }
  internal.params <- do.call("tlqr.control", control)
  seed <- internal.params$seed
  h <- internal.params$h
  parallel <- internal.params$parallel
  ncore <- internal.params$ncore
  set.seed(seed)
  nt <- nrow(x.target)
  # Step 1:
  # For target
  x.target.tilde <- cbind(rep(1, nt), x.target)
  beta_target_hat <- rq.fit.lasso(x.target.tilde, y.target,
                                  tau = u.target)$coefficients
  # CV for h selection if NULL
  if (is.null(h)) {
    h_target <- ChooseBandwidth(x.target, y.target, u.target, control)
  } else{
    h_target <- h[1]
  }
  f_target_hat <- mean(apply(((y.target - x.target.tilde %*% beta_target_hat) /
                      h_target), 1, K)) / h_target
  y.target.tilde <- x.target.tilde %*% beta_target_hat - f_target_hat ^ -1 *
    (ifelse(y.target - x.target.tilde %*% beta_target_hat <= 0, 1, 0) -
       u.target)
  # For source
  if (parallel == TRUE){
    i <- NULL # pass R-CMD check
    ncore <- min(ncore, length(x.source))
    cl <- makeCluster(min(detectCores() - 3, ncore))
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = length(x.source), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    y_aux_tild_bd <- foreach(i = 1:length(x.source),.options.snow = opts,
                             .combine = 'c', .multicombine=T,
                             .packages = c('quantreg'),
                             .inorder = TRUE) %dopar% {
      set.seed(seed)
      x_source <- x.source[[i]]
      ns <- nrow(x_source)
      x_source_tilde <- cbind(rep(1, ns), x_source)
      y_source <- y.source[[i]]
      u_source <- u.source[i]
      beta_aux_hat <- rq.fit.lasso(x_source_tilde, y_source,
                                   tau = u_source)$coefficients
      if (is.null(h)) {
        h_source <- ChooseBandwidth(x_source, y_source, u_source, control)
      } else{
        h_source <- h[k + 1]
      }
      f_aux_hat <- mean(apply(((y_source - x_source_tilde %*% beta_aux_hat) /
                       h_source), 1, K)) / h_source
      list(x_source_tilde %*% beta_aux_hat - f_aux_hat ^ -1 *
          (ifelse(y_source - x_source_tilde %*% beta_aux_hat <= 0, 1, 0) -
             u_source))
    }
    close(pb)
    stopCluster(cl)
  } else {
    y_aux_tild_bd <- list()
    for (k in 1:length(x.source)) {
      set.seed(seed)
      x_source <- x.source[[k]]
      ns <- nrow(x_source)
      x_source_tilde <- cbind(rep(1, ns), x_source)
      y_source <- y.source[[k]]
      u_source <- u.source[k]
      beta_aux_hat <- rq.fit.lasso(x_source_tilde, y_source,
                                   tau = u_source)$coefficients
      if (is.null(h)) {
        h_source <- ChooseBandwidth(x_source, y_source, u_source, control)
      } else{
        h_source <- h[k + 1]
      }
      f_aux_hat <- mean(apply(((y_source - x_source_tilde %*% beta_aux_hat) /
                                 h_source), 1, K)) / h_source
      y_aux_tild_bd[[k]] <- x_source_tilde %*% beta_aux_hat - f_aux_hat ^ -1 *
        (ifelse(y_source - x_source_tilde %*% beta_aux_hat <= 0, 1, 0) -
           u_source)
    }
  }
  # Step 2:
  x_comb <- x.target
  for (k in 1:length(x.source)) {
    x_comb <- rbind(x_comb, x.source[[k]])
  }
  y_comb <- y.target.tilde
  for (k in 1:length(x.source)) {
    y_comb <- c(y_comb, y_aux_tild_bd[[k]])
  }
  if (is.null(lambda1)) {
    lambda1 <- cv.glmnet(x = x_comb, y = y_comb, alpha = 1,
                         standardize = FALSE)$lambda.min
  }
  wfit <- glmnet(x = x_comb, y = y_comb, lambda = lambda1,
                  standardize = FALSE, alpha = 1)
  w_hat <- as.vector(rbind(wfit$a0, wfit$beta))
  # Step 3:
  if (is.null(lambda2)) {
    lambda2 <- lambda1 * sqrt(nrow(x_comb) / nrow(x.target))
  }
  deltafit <- glmnet(x = x.target, y = y.target.tilde -
                       (x.target.tilde %*% w_hat), lambda = lambda2,
                     standardize = FALSE, alpha = 1)
  delta_hat <- as.vector(rbind(deltafit$a0, deltafit$beta))
  return(w_hat + delta_hat)
}


#' Build fusion model between datasets without debias procedure
#' @export
#' @param x.target A n*p matrix as the covariate from the target population.
#' @param y.target A vector with length n as the response from the target
#' population.
#' @param x.source A list object that each element is a n*p matrix as the
#' covariate from informative sources.
#' @param y.source A list object that each element is a vector with length n
#' as the covariate from informative sources.
#' @param u.target A scalar in (0,1). The quantile level for model of target
#' population.
#' @param u.source A vector that contains the quantile level for informative
#' populations.
#' @param lambda1 A scalar for the penalty strength in fusion step. Default
#' is null, if null, a cross validation will be carried for searching.
#' @param control The list object for control, see \code{tlqr.control}.
#' @importFrom quantreg rq.fit.lasso
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach `%dopar%`
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @return A vector. The fusion estimator withou debias procedure based on the
#' given auxiliary sources.
FusionQR <- function(x.target, y.target, x.source, y.source, u.target,
                     u.source, lambda1 = NULL, control = list()) {
  for (k in 1:length(x.source)) {
    if (ncol(x.target) != ncol(x.source[[k]])) {
      print("Two numbers of features don't match!")
      break
    }
  }
  internal.params <- do.call("tlqr.control", control)
  seed <- internal.params$seed
  h <- internal.params$h
  parallel <- internal.params$parallel
  ncore <- internal.params$ncore
  set.seed(seed)
  nt <- nrow(x.target)
  # Step 1:
  # For target
  x.target.tilde <- cbind(rep(1, nt), x.target)
  beta_target_hat <- rq.fit.lasso(x.target.tilde, y.target,
                                  tau = u.target)$coefficients
  if (is.null(h)) {
    h_target <- ChooseBandwidth(x.target, y.target, u.target, control)
  } else{
    h_target <- h[1]
  }
  f_target_hat <- mean(apply(((y.target - x.target.tilde %*% beta_target_hat) /
                      h_target), 1, K)) / h_target
  y.target.tilde <- x.target.tilde %*% beta_target_hat - f_target_hat ^ -1 *
    (ifelse(y.target - x.target.tilde %*% beta_target_hat <= 0, 1, 0) -
       u.target)
  # For source
  if (parallel == TRUE){
    i <- NULL # pass R-CMD check
    ncore <- min(ncore, length(x.source))
    cl <- makeCluster(min(detectCores() - 3, ncore))
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = length(x.source), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    y_aux_tild_bd <- foreach(i = 1:length(x.source),.options.snow = opts,
                             .combine = 'c', .multicombine=T,
                             .packages = c('quantreg'),
                             .inorder = TRUE) %dopar% {
      set.seed(seed)
      x_source <- x.source[[i]]
      ns <- nrow(x_source)
      x_source_tilde <- cbind(rep(1, ns), x_source)
      y_source <- y.source[[i]]
      u_source <- u.source[i]
      beta_aux_hat <- rq.fit.lasso(x_source_tilde, y_source,
                                   tau = u_source)$coefficients
      if (is.null(h)) {
        h_source <- ChooseBandwidth(x_source, y_source, u_source, control)
      } else{
        h_source <- h[k + 1]
      }
      f_aux_hat <- mean(apply(((y_source - x_source_tilde %*% beta_aux_hat) /
                       h_source), 1, K)) / h_source
      list(x_source_tilde %*% beta_aux_hat - f_aux_hat ^ -1 *
          (ifelse(y_source - x_source_tilde %*% beta_aux_hat <= 0, 1, 0) -
             u_source))
    }
    close(pb)
    stopCluster(cl)
  } else {
    y_aux_tild_bd <- list()
    for (k in 1:length(x.source)) {
      set.seed(seed)
      x_source <- x.source[[k]]
      ns <- nrow(x_source)
      x_source_tilde <- cbind(rep(1, ns), x_source)
      y_source <- y.source[[k]]
      u_source <- u.source[k]
      beta_aux_hat <- rq.fit.lasso(x_source_tilde, y_source,
                                   tau = u_source)$coefficients
      if (is.null(h)) {
        h_source <- ChooseBandwidth(x_source, y_source, u_source, control)
      } else{
        h_source <- h[k + 1]
      }
      f_aux_hat <- mean(apply(((y_source - x_source_tilde %*% beta_aux_hat) /
                       h_source), 1, K)) / h_source
      y_aux_tild_bd[[k]] <- x_source_tilde %*% beta_aux_hat - f_aux_hat ^ -1 *
        (ifelse(y_source - x_source_tilde %*% beta_aux_hat <= 0, 1, 0) -
           u_source)
    }
  }
  # Step 2:
  x_comb <- x.target
  for (k in 1:length(x.source)) {
    x_comb <- rbind(x_comb, x.source[[k]])
  }
  y_comb <- y.target.tilde
  for (k in 1:length(x.source)) {
    y_comb <- c(y_comb, y_aux_tild_bd[[k]])
  }
  if (is.null(lambda1)) {
    lambda1 <- cv.glmnet(x = x_comb, y = y_comb, alpha = 1,
                         standardize = FALSE)$lambda.min
  }
  wfit <- glmnet(x = x_comb, y = y_comb, lambda = lambda1,
                  standardize = FALSE, alpha = 1)
  w_hat <- as.vector(rbind(wfit$a0, wfit$beta))
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
#' @param u.target A scalar in (0,1). The quantile level for informative
#' source detection testing.
#' @param control The list object for control, see \code{tlqr.control}.
#' @return A scalar. The quantile loss for given estimator on the given
#' testing samples.
Q_loss <- function(betahat, betaic, X_measure, y_measure, u.target,
                   control = list()) {
  internal.params <- do.call("tlqr.control", control)
  hic <- internal.params$hic
  if (is.null(hic)) {
    hic <- ChooseBandwidth(X_measure, y_measure, u.target, control)
  }
  X_measure_tilde <- cbind(rep(1, nrow(X_measure)), X_measure)
  f_hat <- mean(apply((y_measure - X_measure_tilde %*% betaic / hic), 1, K)) /
    hic
  y_measure_tilde <- X_measure_tilde %*% betaic - f_hat ^ (-1) *
    (ifelse(y_measure - X_measure_tilde %*% betaic <= 0, 1, 0) - u.target)
  return(mean((y_measure_tilde - X_measure_tilde %*% betahat) ^ 2))
}


#' Main function for informative sources detection among all dataset
#' @export
#' @param x.target A n*p matrix as the covariate from the target population.
#' @param y.target A vector with length n as the response from the target
#' population.
#' @param x.source A list object that each element is a n*p matrix as the
#' covariate from informative sources.
#' @param y.source A list object that each element is a vector with length n
#' as the covariate from informative sources.
#' @param u.target A scalar in (0,1). The quantile level for model of target
#' population.
#' @param u.source A vector that contains the quantile level for informative
#' populations.
#' @param epsilon A scalar. Default is 0.01. The strict level for informative
#' sources detectiond, larger the value is, less strict the procedure is.
#' @param lambda1 A scalar for the penalty strength in fusion step. Default
#' is null, if null, a cross validation will be carried for searching.
#' @param ninfo An integar. Default is NULL. The given number of informative
#' sources under pseduo running.
#' @param control The list object for control, see \code{tlqr.control}.
#' @importFrom quantreg rq.fit.lasso
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom foreach foreach `%dopar%`
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils head setTxtProgressBar txtProgressBar
#' @return A vector of index. The index of the informative sources.
DetectQR <- function(x.target, y.target, x.source, y.source, u.target,
                     u.source, epsilon = 0.01, lambda1 = NULL, ninfo = NULL,
                     control = list()){
  if (length(x.source) != length(y.source)) {
    print('the # of datasets for x and y are not agreed!')
  } else {
    internal.params <- do.call("tlqr.control", control)
    seed <- internal.params$seed
    h <- internal.params$h
    parallel <- internal.params$parallel
    ncore <- internal.params$ncore
    verbose <- internal.params$verbose
    set.seed(seed)
    M <- length(x.source)
    # Step 1.1: cut the target half and half: (remark: 1 for I, 2 for Ic)
    target_index <- 1:nrow(x.target)
    X0_index <- sample(target_index, size = round(nrow(x.target) / 2),
                       replace = FALSE)
    X0_cut1 <- x.target[X0_index,]
    X0_cut2 <- x.target[-X0_index,]
    y0_cut1 <- y.target[X0_index]
    y0_cut2 <- y.target[-X0_index]
    # CV1:
    X0_cut1_tilde <- cbind(rep(1, nrow(X0_cut1)), X0_cut1)
    # Step 1.2: train sparse quantile regression for I for X0, y0:
    beta_0_hat <- rq.fit.lasso(X0_cut1_tilde, y0_cut1,
                               tau = u.target)$coefficients
    # Step 1.3: train sparse quantile regression for I for X0, y0 + k-th
    # source:
    if (parallel == FALSE){
      for (k in 1:length(x.source)) {
        assign(paste0(paste0('beta_', k), '_hat'),
               FusionQR(X0_cut1, y0_cut1, list(x.source[[k]]),
                        list(y.source[[k]]), u.target, c(u.source[k]),
                        lambda1 = lambda1, control))
      }
    } else{
      i <- NULL # pass R-CMD check
      ncore <- min(ncore, M)
      cl <- makeCluster(min(detectCores() - 3, ncore))
      registerDoSNOW(cl)
      pb <- txtProgressBar(max = M, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      betasources <- foreach(i = 1:M, .options.snow = opts, .combine = cbind,
                             .packages = c('quantreg', 'glmnet'),
                             .inorder = TRUE) %dopar% {
        FusionQR(X0_cut1, y0_cut1, list(x.source[[i]]), list(y.source[[i]]),
                 u.target, c(u.source[i]), lambda1 = lambda1, control)
      }
      close(pb)
      stopCluster(cl)
      for (k in 1:length(x.source)) {
        assign(paste0(paste0('beta_', k), '_hat'), betasources[, k])
      }
    }
    ## Step 2: Calculate the loss based on the loss function:
    X0_cut2_tilde <- cbind(rep(1, nrow(X0_cut2)), X0_cut2)
    beta_0_hatic <- rq.fit.lasso(X0_cut2_tilde, y0_cut2,
                                 tau = u.target)$coefficients
    ## hic <- ChooseBandwidth(X0_cut2, y0_cut2, u.target)
    hic <- h[1]
    for (k in c(0, (1:length(x.source)))) {
      assign(paste0('Q', k),
             Q_loss(get(paste0(paste0('beta_', k), '_hat')), beta_0_hatic,
                    X0_cut2, y0_cut2, u.target, control = list(hic = hic)))
    }
    ## Step 3: Get informative sets:
    ## Step 3.1 true
    index_true1 <- c()
    if (verbose == TRUE) {
      print(get(paste0('Q', 0)))
    }
    for (k in 1:length(x.source)) {
      if (verbose == TRUE) {
        print(get(paste0('Q', k)))
      }
      if (get(paste0('Q', k)) <= (1 + epsilon) * get(paste0('Q', 0))) {
        index_true1 <- c(index_true1, k)
      } else{
        next
      }
    }
    ## Step 3.2 psd
    index_psd1 <- NULL
    if (!is.null(ninfo)) {
      Q_K_vec <- c()
      for (k in 1:length(x.source)) {
        Q_K_vec <- c(Q_K_vec, get(paste0('Q', k)))
      }
      index_psd1 <- match(head(sort(Q_K_vec, decreasing = FALSE), ninfo),
                          Q_K_vec)
    } else{
      index_psd1 <- c()
    }
    ## CV2
    ## Step 1.2: train sparse quantile regression for I for X0, y0:
    beta_0_hat <- rq.fit.lasso(X0_cut2_tilde, y0_cut2,
                               tau = u.target)$coefficients
    # Step 1.3: train sparse quantile regression for I for X0, y0 + k-th
    # source:
    if (parallel == FALSE){
      for (k in 1:length(x.source)) {
        assign(paste0(paste0('beta_', k), '_hat'),
               FusionQR(X0_cut2, y0_cut2, list(x.source[[k]]),
                        list(y.source[[k]]), u.target, c(u.source[k]),
                        lambda1 = lambda1, control))
      }
    } else{
      i <- NULL # pass R-CMD check
      ncore <- min(ncore, M)
      cl <- makeCluster(min(detectCores() - 3, ncore))
      registerDoSNOW(cl)
      pb <- txtProgressBar(max = M, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      betasources <- foreach(i = 1:M, .options.snow = opts, .combine = cbind,
                             .packages = c('quantreg', 'glmnet'),
                             .inorder = TRUE) %dopar% {
        FusionQR(X0_cut2, y0_cut2, list(x.source[[i]]), list(y.source[[i]]),
                 u.target, c(u.source[i]), lambda1 = lambda1, control)
      }
      close(pb)
      stopCluster(cl)
      for (k in 1:length(x.source)) {
        assign(paste0(paste0('beta_', k), '_hat'), betasources[, k])
      }
    }
    ## Step 2: Calculate the loss based on the loss function:
    beta_0_hatic <- rq.fit.lasso(X0_cut1_tilde, y0_cut1,
                                 tau = u.target)$coefficients
    ## hic <- ChooseBandwidth(X0_cut1,y0_cut1,u.target)
    hic <- h[1]
    for (k in c(0, (1:length(x.source)))) {
      assign(paste0('Q', k),
             Q_loss(get(paste0(paste0('beta_', k), '_hat')), beta_0_hatic,
                    X0_cut1, y0_cut1, u.target, control = list(hic = hic)))
    }
    # Step 3: Get informative sets:
    # Step 3.1 true
    index_true2 <- c()
    if (verbose == TRUE) {
      print(get(paste0('Q', 0)))
    }
    for (k in 1:length(x.source)) {
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
    if (!is.null(ninfo)) {
      Q_K_vec <- c()
      for (k in 1:length(x.source)) {
        Q_K_vec <- c(Q_K_vec, get(paste0('Q', k)))
      }
      index_psd2 <- match(head(sort(Q_K_vec, decreasing = FALSE), ninfo),
                          Q_K_vec)
    } else{
      index_psd2 <- c()
    }
    index_true <- intersect(index_true1, index_true2)
    index_psd <- intersect(index_psd1, index_psd2)
    return(list('true' = index_true, 'psd' = index_psd))
  }
}

#' Trans-Lasso QR with informative set detection
#' @export
#' @param x.target A n*p matrix as the covariate from the target population.
#' @param y.target A vector with length n as the response from the target
#' population.
#' @param x.source A list object that each element is a n*p matrix as the
#' covariate from informative sources.
#' @param y.source A list object that each element is a vector with length n
#' as the covariate from informative sources.
#' @param u.target A scalar in (0,1). The quantile level for model of target
#' population.
#' @param u.source A vector that contains the quantile level for informative
#' populations.
#' @param epsilon A scalar. Default is 0.01. The strict level for informative
#' sources detectiond, larger the value is, less strict the procedure is.
#' @param ninfo An integar. Default is NULL. The given number of informative
#' sources under pseduo running.
#' @param lambda.det A scalar for the penalty strength in detection procedure.
#' Default is null, if null, a cross validation will be carried for searching.
#' @param lambda1 A scalar for the penalty strength in fusion step. Default
#' is null, if null, a cross validation will be carried for searching.
#' @param lambda2 A scalar for the penalty strength in debiasing step. Default
#' is null, if null, lambda2 will be calculated based on lambda1 according to
#' the proportion of sample sizes as described in our paper.
#' @param control The list object for control, see \code{tlqr.control}.
#' @importFrom quantreg rq.fit.lasso
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach `%dopar%`
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @return A vector. The corresponding transfer learning estimator based on
#' the given auxiliary sources.
DetectTransferQR <- function(x.target, y.target, x.source, y.source, u.target,
                             u.source, epsilon = 0.01, ninfo = NULL,
                             lambda.det = NULL, lambda1 = NULL, lambda2 = NULL,
                             control = list()) {
  for (k in 1:length(x.source)) {
    if (ncol(x.target) != ncol(x.source[[k]])) {
      print("Two numbers of features don't match!")
      break
    }
  }
  info_det <- DetectQR(x.target, y.target, x.source, y.source, u.target,
                       u.source, epsilon, lambda1 = lambda.det, ninfo, control)
  if (is.null(ninfo)){
    info_index <- info_det$true
    if (length(info_index) != 0) {
      x.source <- lapply(1:length(info_index), function(k) {
        return(x.source[[info_index[k]]])})
      y.source <- lapply(1:length(info_index), function(k) {
        return(y.source[[info_index[k]]])})
      u.source <- sapply(1:length(info_index), function(k) {
        return(u.source[info_index[k]])})
      beta_hat <- TransferQR(x.target, y.target, x.source, y.source, u.target,
                             u.source, lambda1, lambda2, control)
    } else{
      beta_hat <- rq.fit.lasso(cbind(rep(1, nrow(x.target)), x.target),
                               y.target, tau = u.target)$coefficients
    }
  } else {
    info_index <- info_det$psd
    if (length(info_index) != 0) {
      x.source <- lapply(1:length(info_index), function(k) {
        return(x.source[[info_index[k]]])})
      y.source <- lapply(1:length(info_index), function(k) {
        return(y.source[[info_index[k]]])})
      u.source <- sapply(1:length(info_index), function(k) {
        return(u.source[info_index[k]])})
      beta_hat <- TransferQR(x.target, y.target, x.source, y.source, u.target,
                             u.source, lambda1, lambda2, control)
    } else{
      beta_hat <- rq.fit.lasso(cbind(rep(1, nrow(x.target)), x.target),
                               y.target, tau = u.target)$coefficients
    }
  }
  return(beta_hat)
}



#' Oracle Trans-pooling QR with a given informative sources set
#' @export
#' @param x.target A n*p matrix as the covariate from the target population.
#' @param y.target A vector with length n as the response from the target
#' population.
#' @param x.source A list object that each element is a n*p matrix as the
#' covariate from informative sources.
#' @param y.source A list object that each element is a vector with length n
#' as the covariate from informative sources.
#' @param u A scalar in (0,1). The quantile level for all models.
#' @importFrom quantreg rq.fit.lasso
#' @return A vector. The corresponding transfer learning estimator based on
#' the given auxiliary sources based on Huang 2022.
TransferPoolQR <- function(x.target, y.target, x.source, y.source, u){
  for (k in 1:length(x.source)) {
    if (ncol(x.target) != ncol(x.source[[k]])) {
      print("Two numbers of features don't match!")
      break
    }
  }
  # fusion
  x.target.tilde <- cbind(rep(1, nrow(x.target)), x.target)
  x_comb <- x.target
  for (k in 1:length(x.source)){
    x_comb <- rbind(x_comb, x.source[[k]])
  }
  x_comb_tilde <- cbind(rep(1, nrow(x_comb)), x_comb)
  y_comb <- y.target
  for (k in 1:length(x.source)){
    y_comb <- c(y_comb, y.source[[k]])
  }
  beta_fusion <- rq.fit.lasso(x_comb_tilde, y_comb, tau = u)$coefficients
  # debias
  beta_db <- rq.fit.lasso(x.target.tilde,
                          y.target - (x.target.tilde %*% beta_fusion),
                          tau = u)$coefficients
  return(beta_fusion + beta_db)
}

#' Build fusion model for Trans-pooling QR between datasets without debias
#' procedure
#' @export
#' @param x.target A n*p matrix as the covariate from the target population.
#' @param y.target A vector with length n as the response from the target
#' population.
#' @param x.source A list object that each element is a n*p matrix as the
#' covariate from informative sources.
#' @param y.source A list object that each element is a vector with length n
#' as the covariate from informative sources.
#' @param u A scalar in (0,1). The quantile level for all models.
#' @importFrom quantreg rq.fit.lasso
#' @return A vector. The fusion estimator withou debias procedure based on the
#' given auxiliary sources.
FusionPoolQR <- function(x.target, y.target, x.source, y.source, u) {
  for (k in 1:length(x.source)) {
    if (ncol(x.target) != ncol(x.source[[k]])) {
      print("Two numbers of features don't match!")
      break
    }
  }
  # fusion
  x_comb <- x.target
  for (k in 1:length(x.source)){
    x_comb <- rbind(x_comb, x.source[[k]])
  }
  x_comb_tilde <- cbind(rep(1, nrow(x_comb)), x_comb)
  y_comb <- y.target
  for (k in 1:length(x.source)){
    y_comb <- c(y_comb, y.source[[k]])
  }
  beta_fusion <- rq.fit.lasso(x_comb_tilde, y_comb, tau = u)$coefficients
  return(beta_fusion)
}


#' Main function for informative sources detection among all dataset for
#' Trans-pooling QR
#' @export
#' @param x.target A n*p matrix as the covariate from the target population.
#' @param y.target A vector with length n as the response from the target
#' population.
#' @param x.source A list object that each element is a n*p matrix as the
#' covariate from informative sources.
#' @param y.source A list object that each element is a vector with length n
#' as the covariate from informative sources.
#' @param u A scalar in (0,1). The quantile level for all models.
#' @param epsilon A scalar. Default is 0.01. The strict level for informative
#' sources detectiond, larger the value is, less strict the procedure is.
#' @param ninfo An integar. Default is NULL. The given number of informative
#' sources under pseduo running.
#' @param control The list object for control, see \code{tlqr.control}.
#' @importFrom quantreg rq.fit.lasso
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom foreach foreach `%dopar%`
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils head setTxtProgressBar txtProgressBar
#' @return A vector of index. The index of the informative sources.
DetectPoolQR <- function(x.target, y.target, x.source, y.source, u,
                         epsilon = 0.01, ninfo = NULL, control = list()){
  if (length(x.source) != length(y.source)) {
    print('the # of datasets for x and y are not agreed!')
  } else {
    internal.params <- do.call("tlqr.control", control)
    seed <- internal.params$seed
    parallel <- internal.params$parallel
    ncore <- internal.params$ncore
    verbose <- internal.params$verbose
    set.seed(seed)
    M <- length(x.source)
    # Step 1.1: cut the target half and half: (remark: 1 for I, 2 for Ic)
    target_index <- 1:nrow(x.target)
    X0_index <- sample(target_index, size = round(nrow(x.target) / 2),
                       replace = FALSE)
    X0_cut1 <- x.target[X0_index,]
    X0_cut2 <- x.target[-X0_index,]
    y0_cut1 <- y.target[X0_index]
    y0_cut2 <- y.target[-X0_index]
    # CV1:
    # Step 1.2: train sparse quantile regression for I for X0, y0:
    X0_cut1_tilde <- cbind(rep(1, nrow(X0_cut1)), X0_cut1)
    beta_0_hat <- rq.fit.lasso(X0_cut1_tilde, y0_cut1, tau=u)$coefficients
    # Step 1.3: train sparse quantile regression for I for X0, y0 + k-th
    # source:
    if (parallel == FALSE){
      for (k in 1:length(x.source)) {
        assign(paste0(paste0('beta_', k), '_hat'),
               FusionPoolQR(X0_cut1, y0_cut1, list(x.source[[k]]),
                            list(y.source[[k]]), u))
      }
    } else{
      i <- NULL # pass R-CMD check
      ncore <- min(ncore, M)
      cl <- makeCluster(min(detectCores() - 3, ncore))
      registerDoSNOW(cl)
      pb <- txtProgressBar(max = M, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      betasources <- foreach(i = 1:M, .options.snow = opts, .combine = cbind,
                             .packages = c('quantreg'),
                             .inorder = TRUE) %dopar% {
        FusionPoolQR(X0_cut1, y0_cut1, list(x.source[[i]]),
                     list(y.source[[i]]),u)
      }
      close(pb)
      stopCluster(cl)
      for (k in 1:length(x.source)) {
        assign(paste0(paste0('beta_', k), '_hat'), betasources[, k])
      }
    }
    # Step 2: Calculate the loss based on the loss function:
    for (k in c(0, (1:length(x.source)))) {
      assign(paste0('Q', k),
             QuantileLoss(X0_cut2, y0_cut2,
                          get(paste0(paste0('beta_', k), '_hat')), u))
    }
    # Step 3: Get informative sets:
    # Step 3.1 true
    index_true1 <- c()
    if (verbose == TRUE) {
      print(get(paste0('Q', 0)))
    }
    for (k in 1:length(x.source)) {
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
    if (!is.null(ninfo)) {
      Q_K_vec <- c()
      for (k in 1:length(x.source)) {
        Q_K_vec <- c(Q_K_vec, get(paste0('Q', k)))
      }
      index_psd1 <- match(head(sort(Q_K_vec, decreasing = FALSE), ninfo),
                          Q_K_vec)
    } else{
      index_psd1 <- c()
    }
    # CV2
    # Step 1.2: train sparse quantile regression for I for X0, y0:
    X0_cut2_tilde <- cbind(rep(1, nrow(X0_cut2)), X0_cut2)
    beta_0_hat <- rq.fit.lasso(X0_cut2_tilde, y0_cut2, tau = u)$coefficients
    # Step 1.3: train sparse quantile regression for I for X0, y0 + k-th
    if (parallel == FALSE){
      for (k in 1:length(x.source)) {
        assign(paste0(paste0('beta_', k), '_hat'),
               FusionPoolQR(X0_cut2, y0_cut2, list(x.source[[k]]),
                            list(y.source[[k]]), u))
      }
    } else{
      i <- NULL # pass R-CMD check
      ncore <- min(ncore, M)
      cl <- makeCluster(min(detectCores() - 3, ncore))
      registerDoSNOW(cl)
      pb <- txtProgressBar(max = M, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      betasources <- foreach(i = 1:M, .options.snow = opts, .combine = cbind,
                             .packages = c('quantreg'),
                             .inorder = TRUE) %dopar% {
        FusionPoolQR(X0_cut2, y0_cut2, list(x.source[[i]]),
                     list(y.source[[i]]),u)
      }
      close(pb)
      stopCluster(cl)
      for (k in 1:length(x.source)) {
        assign(paste0(paste0('beta_', k), '_hat'), betasources[, k])
      }
    }
    # Step 2: Calculate the loss based on the loss function:
    for (k in c(0, (1:length(x.source)))) {
      assign(paste0('Q', k),
             QuantileLoss(X0_cut1, y0_cut1,
                          get(paste0(paste0('beta_', k), '_hat')), u))
    }
    # Step 3: Get informative sets:
    # Step 3.1 true
    index_true2 <- c()
    if (verbose == TRUE) {
      print(get(paste0('Q', 0)))
    }
    for (k in 1:length(x.source)) {
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
    if (!is.null(ninfo)) {
      Q_K_vec <- c()
      for (k in 1:length(x.source)) {
        Q_K_vec <- c(Q_K_vec, get(paste0('Q', k)))
      }
      index_psd2 <- match(head(sort(Q_K_vec, decreasing = FALSE), ninfo),
                          Q_K_vec)
    } else{
      index_psd2 <- c()
    }
    index_true <- intersect(index_true1, index_true2)
    index_psd <- intersect(index_psd1, index_psd2)
    return(list('true' = index_true, 'psd' = index_psd))
  }
}


#' Trans-Lasso QR with informative set detection
#' @export
#' @param x.target A n*p matrix as the covariate from the target population.
#' @param y.target A vector with length n as the response from the target
#' population.
#' @param x.source A list object that each element is a n*p matrix as the
#' covariate from informative sources.
#' @param y.source A list object that each element is a vector with length n
#' as the covariate from informative sources.
#' @param u A scalar in (0,1). The quantile level for all models.
#' @param epsilon A scalar. Default is 0.01. The strict level for informative
#' sources detectiond, larger the value is, less strict the procedure is.
#' @param ninfo An integar. Default is NULL. The given number of informative
#' sources under pseduo running.
#' @param control The list object for control, see \code{tlqr.control}.
#' @importFrom quantreg rq.fit.lasso
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach `%dopar%`
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @return A vector. The corresponding transfer learning estimator based on
#' the given auxiliary sources.
DetectTransferPoolQR <- function(x.target, y.target, x.source, y.source,
                                 u, epsilon = 0.01, ninfo = NULL,
                                 control = list()) {
  for (k in 1:length(x.source)) {
    if (ncol(x.target) != ncol(x.source[[k]])) {
      print("Two numbers of features don't match!")
      break
    }
  }
  info_det <- DetectPoolQR(x.target, y.target, x.source, y.source, u, epsilon,
                           ninfo, control)
  if (is.null(ninfo)){
    info_index <- info_det$true
    if (length(info_index) != 0) {
      x.source <- lapply(1:length(info_index), function(k) {
        return(x.source[[info_index[k]]])})
      y.source <- lapply(1:length(info_index), function(k) {
        return(y.source[[info_index[k]]])})
      beta_hat <- TransferPoolQR(x.target, y.target, x.source, y.source, u)
    } else{
      beta_hat <- rq.fit.lasso(cbind(rep(1, nrow(x.target)), x.target),
                               y.target, tau = u)$coefficients
    }
  } else {
    info_index <- info_det$psd
    if (length(info_index) != 0) {
      x.source <- lapply(1:length(info_index), function(k) {
        return(x.source[[info_index[k]]])})
      y.source <- lapply(1:length(info_index), function(k) {
        return(y.source[[info_index[k]]])})
      beta_hat <- TransferPoolQR(x.target, y.target, x.source, y.source, u)
    } else{
      beta_hat <- rq.fit.lasso(cbind(rep(1, nrow(x.target)), x.target),
                               y.target, tau = u)$coefficients
    }
  }
  return(beta_hat)
}


