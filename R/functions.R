# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Quantile loss
#' @description Calculate the quantile loss based on the given samples, coefficient, and quantile level.
#' @export
#' @param x Input matrix. The convariates of samples.
#' @param y Input vector. The responses of samples.
#' @param beta Input vector. The coefficient.
#' @param u Input scalar. The quantile level.
#' @return A scalar. The quantile loss based on the given inputs.
cdloss <- function(x, y, beta, u) {
  return(sum((y - x %*% beta) * (u - ifelse(y - x %*% beta <= 0, 1, 0))) / nrow(x))
}

#' Simulated data generation
#' @description Data generation function to implement the simulation settings in our paper.
#' @export
#' @param p Input integer. Default is 150. The dimensionality for the covariate.
#' @param n Input integer. Default is 200. The sample size for each datasets, including all the target and sources.
#' @param s Input integer. Default is 5. The sparsity level in the simulation.
#' @param d Input scalar. Default is 2. Control the scale of the Laplace distribution, larger d leads to smaller similarity between the target and informative sets. Details are in the section for simulation.
#' @param An Input integer. Default is 5. The number of informative set.
#' @param M Input integer. Default is 10. The total number of sources.
#' @param eta Input scalar. Default is 20. The signal-to-noise ratio for residual.
#' @param cov_type Input string. Default is 'auto'. Control the covariance matrix of datasets. The value 'auto' leads to a homogeneous design with auto-covariance; the value 'toeplitz' leads to heterogeneous designs with Toeplitz covariance matrix. No other value is permitted.
#' @param res_type Input string. Default is 'normal'. Control the type of distribution for residual. Only one from 'normal' and 'cauchy' is allowed.
#' @param seed Input integer. The default is 111.  The random seed.
#' @import stats
#' @import MASS
#' @importFrom ExtDist rLaplace
#' @return NULL. But will generate correspoding pairs of (X,y) in global enviroment.
data_generation <- function(p=150,n=200,s=5,d=2,An=5,M=10,eta=20,cov_type='auto',res_type='normal',seed=111){
  set.seed(seed)
  ## Construct beta
  beta_0 <- rep(0, p)
  beta_0[1:s] <- 1
  assign('beta_0', beta_0, envir = .GlobalEnv)
  A_index <- c(rep(1, An), rep(0, M - An))
  for (i in 1:length(A_index)) {
    if (A_index[i] == 1) {
      Hk <- sample(seq(from = 1, to = p, by = 1),
                   size = p / 2,
                   replace = F)
      ksi <- rLaplace(n = length(Hk),
                      mu = 0,
                      b = 2 * d / p)
      ksi_add <- rep(0, p)
      ksi_add[Hk] <- ksi
      nam <- paste0("beta_", i)
      assign(nam, beta_0 + ksi_add, envir = .GlobalEnv)
    } else{
      Hk <- sample(seq(from = 1, to = p, by = 1),
                   size = p / 2,
                   replace = F)
      ksi <- rLaplace(n = length(Hk),
                      mu = 0,
                      b = 140 / p)
      ksi_add <- rep(0, p)
      ksi_add[Hk] <- ksi
      nam <- paste0("beta_", i)
      assign(nam, beta_0 + ksi_add, envir = .GlobalEnv)
    }
  }
  ## Construct covariates:
  if (cov_type == 'auto') {
    sigma <- matrix(rep(0, p * p), nrow = p, byrow = T)
    for (i in 1:p) {
      for (j in 1:i) {
        sigma[i, j] <- 0.5 ^ abs(i - j)
        sigma[j, i] <- sigma[i, j]
      }
    }
    for (i in 0:M) {
      nam <- paste0("X_", i)
      if (i == 0) {
        assign(nam, MASS::mvrnorm(n = floor(n / 0.8) + 1, rep(0, p), sigma, tol = 1e-10),envir = .GlobalEnv)
      } else {
        assign(nam, MASS::mvrnorm(n = n, rep(0, p), sigma, tol = 1e-10),envir = .GlobalEnv)
      }
      assign(paste0("sigma_", i), sigma)
    }
  } else if (cov_type == 'toeplitz') {
    for (i in 0:M) {
      if (i != 0) {
        row1 <- c(1, rep(1 / (i + 1), 2 * i - 1), rep(0, p - 2 * i))
        sigma <- stats::toeplitz(row1)
      } else{
        sigma <- diag(rep(1, p))
      }
      nam <- paste0("X_", i)
      if (i == 0) {
        assign(nam, MASS::mvrnorm(n = floor(n / 0.8) + 1, rep(0, p), sigma, tol = 1e-10),envir = .GlobalEnv)
      } else{
        assign(nam, MASS::mvrnorm(n = n, rep(0, p), sigma, tol = 1e-10),envir = .GlobalEnv)
      }
      assign(paste0("sigma_", i), sigma)
    }
  } else{
    print('wrong covariance type!')
  }
  ## Construct response:
  if (res_type=='normal'){
    assign("y_0", get(paste0('X_',0)) %*% get(paste0('beta_',0)) +
             stats::rnorm(n = floor(n/0.8)+1,qnorm(0.2,mean=0,sd= t(get(paste0('beta_',0)))%*% get(paste0('sigma_',0)) %*% get(paste0('beta_',0))/eta),
                   sd =t(get(paste0('beta_',0)))%*% get(paste0('sigma_',0)) %*% get(paste0('beta_',0))/eta),envir = .GlobalEnv)
    for (i in 1:M){
      nam <- paste0("y_", i)
      assign(nam, get(paste0('X_',i)) %*% get(paste0('beta_',i)) +
               stats::rnorm(n = n,qnorm(0.2,mean=0,sd= t(get(paste0('beta_',i)))%*% get(paste0('sigma_',i)) %*% get(paste0('beta_',i))/eta),
                     sd =t(get(paste0('beta_',i)))%*% get(paste0('sigma_',i)) %*% get(paste0('beta_',i))/eta),envir = .GlobalEnv)
    }
  } else if (res_type=='cauchy'){
    assign("y_0", get(paste0('X_',0)) %*% get(paste0('beta_',0)) +
             stats::rcauchy(n = n,qcauchy(0.2,location=0,scale=  t(get(paste0('beta_',0)))%*% get(paste0('sigma_',0)) %*% get(paste0('beta_',0))/eta),
                     scale= t(get(paste0('beta_',0)))%*% get(paste0('sigma_',0)) %*% get(paste0('beta_',0))/eta),envir = .GlobalEnv)
    for (i in 1:M){
      nam <- paste0("y_", i)
      assign(nam, get(paste0('X_',i)) %*% get(paste0('beta_',i)) +
               stats::rcauchy(n = n,qcauchy(0.2,location=0,scale=  t(get(paste0('beta_',i)))%*% get(paste0('sigma_',i)) %*% get(paste0('beta_',i))/eta),
                       scale= t(get(paste0('beta_',i)))%*% get(paste0('sigma_',i)) %*% get(paste0('beta_',i))/eta),envir = .GlobalEnv)
    }
  } else {
    print('Error on residual setting!')
  }
}

#' Smoothing kernel
#' @description Built-in function for smoothing kernel. Not exported to users.
#' @param x Numeric variable
#' @return Numeric scalar
K <- function(x) {
  if (x <= -1) {
    return(0)
  } else if (x >= 1) {
    return(0)
  } else{
    return((-315 / 64) * x ^ 6 + (735 / 64) * x ^ 4 - (525 / 64) * x ^ 2 + (105 / 64))
  }
}

#' Tunning the hyper-parameter for high-dimensional model
#' @description Build the function for tunning the parameter based on given data. In this function, the given samples will be split into training data and validation data according to the given proportion variable 20 times. Then, for each of the parameters generated by the grid from e^lower to the e^upper with a fixed step length, a high-dimensional model will be trained 20 times, and the function will output the best parameter that has the minimum average loss from 20 times repetition.
#' @export
#' @param x A n*p matrix as the training and validation covariate.
#' @param y A vector with length n as the response from training and validation.
#' @param u A scalar. The quantile level for model.
#' @param size A scalar in [0,1]. default is 0.05, the proportion of validation sample.
#' @param lower A scalar. The default is -2, the lower bound for grid search, established by exp(lower).
#' @param upper A scalar larget than lower. The default is +2, the upper bound for grid search, established by exp(upper).
#' @param step A scalar. The default is 0.05, the length of step for the grid search.
#' @param machines An integer. The default is 10, the number of machines for parellel computation.
#' @param seed An integer. The default is 111. The random seed.
#' @import doParallel
#' @import foreach
#' @import glmnet
#' @import quantreg
#' @import parallel
#' @return A scalar. The best lambda value for penalized model.
tune_param <- function(x,y,u,size=0.05,lower=-2,upper=2,step=0.05,machines=10,seed=111){
  n <- nrow(x)
  lambda_seq = exp(seq(lower,upper,by=step))
  tryCatch({
    cl <- parallel::makeCluster(machines)
    doParallel::registerDoParallel(cl)
    pd_error <- foreach::foreach(i = 1:length(lambda_seq), .combine=rbind,.inorder = FALSE, .packages =c('quantreg')) %dopar% {
      lam = lambda_seq[i]
      iter_error = c()
      set.seed(111)
      for (T in 1:20){
        # train-test-split
        train_index <- sample(1:n,size=(1-size)*n,replace=F)
        X_0_train <- x[train_index,]
        X_0_test <- x[-train_index,]
        y_0_train <- y[train_index]
        y_0_test <- y[-train_index]
        duplicated.columns <- duplicated(t(X_0_train))
        X_0_train <- X_0_train[,!duplicated.columns]
        X_0_test <- X_0_test[,!duplicated.columns]
        beta_hat <-quantreg::rq.fit.lasso(X_0_train,y_0_train,tau=u,lambda = lam)$coefficients
        iter_error <- c(iter_error,sum((y_0_test- X_0_test %*% beta_hat)*(u- ifelse(y_0_test-X_0_test %*% beta_hat<=0,1,0))))
      }
      df_add <- data.frame(lambda=c(lam),performance = c(mean(iter_error)))
      df_add
    }
    parallel::stopCluster(cl)
  }, error = function(cond){
    return(tune_param(x,y,u,size=0.05,lower=-2,upper=2,step=0.05,machines=10,seed=seed+1))
  })
  return(pd_error[which.min(pd_error$performance),]$lambda)
}


#' Get the smoothened response
#' @description Get the smoothened response (from the quantile model to least square model) variable by iterations.
#' @export
#' @param x A n*p matrix as the covariate for response smoothening.
#' @param y A vector with length n as the response that needs to be smoothened.
#' @param u A scalar. The quantile level.
#' @param maxit A integer. The maximum iteration number.
#' @param reltol A scalar. Default is 1e-2. The relative tolerance for the stopping criterion.
#' @import quantreg
#' @import glmnet
#' @return y: A vector as the smoothed response. h: A scalar as the selected bandwidth. converge: The convergence status, 0 is succuessful convergence.
smoothy <- function(x,y,u,maxit=100,reltol=5e-2){
  beta_old <- quantreg::rq.fit.lasso(x,y,u)$coefficients
  converge <- 0
  for (g in 1:maxit){
    s0 <- sum(abs(beta_old)>=(max(beta_old)/5))
    h_target <- sqrt(s0*log(nrow(x))/nrow(x))+(s0^((2*(g-1)+1)/2))*(log(nrow(x))/nrow(x))^(g/2)
    f_target_hat <- 0
    while (f_target_hat <=0){
      f_target_hat <- sum(apply(((y - x %*% beta_old) / h_target),1,K))/(nrow(x)*h_target)
      h_target <- h_target * 1.05
    }
    y_tild <- x %*% beta_old - f_target_hat^-1 *(ifelse(y-x %*% beta_old<=0,1,0)-u)
    fit <- glmnet::cv.glmnet(x=x,y = y_tild,alpha=1)
    beta_new <- as.vector(coef(fit, s = "lambda.min"))[2:length(coef(fit, s = "lambda.min"))]
    if (sum(abs(beta_new-beta_old))>reltol*sum(abs(beta_old))){
      beta_old <- beta_new
    } else {
      if (f_target_hat^-1 >= 1e5){
        converge <- 1
      }
      return(list(y=y_tild,h=h_target,converge=converge))
    }
  }
  if (f_target_hat^-1 >= 1e5){
    converge <- 1
  }
  return(list(y=y_tild,h=h_target,converge=converge))
}

#' Oracle Trans-Lasso QR
#' @description Implement Oracle Trans-Lasso QR with a given informative sources set.
#' @export
#' @param x_target A n*p matrix as the covariate from the target population.
#' @param y_target A vector with length n as the response from the target population.
#' @param x_aux_bd A list object that each element is a n*p matrix as the covariate from informative sources.
#' @param y_aux_bd A list object that each element is a vector with length n as the covariate from informative sources.
#' @param u_target A scalar in (0,1). The quantile level for model of target population.
#' @param u_aux_bd A vector that contains the quantile level for informative populations.
#' @param mode A string. Default is 'real'. Only 'simulation' or 'real' is allowed. The difference of mode leads to different operations for the selection of strength of penalty in single-source modeling. 'real' mode takes longer time.
#' @import quantreg
#' @import glmnet
#' @return A vector. The corresponding transfer learning estimator based on the given auxiliary sources.
rq.transfer <- function(x_target,y_target,x_aux_bd,y_aux_bd,u_target,u_aux_bd,mode='real'){
  for (k in 1:length(x_aux_bd)){
    if (ncol(x_target)!=ncol(x_aux_bd[[k]])){
      print("Two numbers of features don't match!")
      break
    }
  }
  n <- nrow(x_target)
  if (mode=='simulation'){
    lambda_1 <- tune_param(x_target,y_target,u_target)
    # Step 1:
    ## For target
    beta_target_hat <- quantreg::rq.fit.lasso(x_target,y_target,u_target,lambda_1)$coefficients
    h_target <- sqrt(6*log(nrow(x_target))/nrow(x_target))+sqrt(6)^-1*(6^2*log(nrow(x_target))/nrow(x_target))^0.5
    f_target_hat <- sum(apply((y_target - x_target %*% beta_target_hat / h_target),1,K))/(nrow(x_target)*h_target)
    y_target_tild <- x_target %*% beta_target_hat - f_target_hat^-1 *(ifelse(y_target-x_target %*% beta_target_hat<=0,1,0)-u_target)
    ## For source
    y_aux_tild_bd <- list()
    for (k in 1:length(x_aux_bd)){
      beta_aux_hat <-  quantreg::rq.fit.lasso(x_aux_bd[[k]],y_aux_bd[[k]],u_aux_bd[k],lambda_1)$coefficients
      h_source <- sqrt(6*log(nrow(x_aux_bd[[k]]))/nrow(x_aux_bd[[k]]))+sqrt(6)^-1*(6^2*log(nrow(x_aux_bd[[k]]))/nrow(x_aux_bd[[k]]))^0.5
      f_aux_hat <- sum(apply((y_aux_bd[[k]] - x_aux_bd[[k]] %*% beta_aux_hat / h_source),1,K))/(nrow(x_aux_bd[[k]])*h_source)
      y_aux_tild_bd[[k]] <- x_aux_bd[[k]] %*% beta_aux_hat - f_aux_hat^-1 *(ifelse(y_aux_bd[[k]]-x_aux_bd[[k]] %*% beta_aux_hat<=0,1,0)-u_aux_bd[k])
    }
    # Step 2:
    x_comb <- x_target
    for (k in 1:length(x_aux_bd)){
      x_comb <- rbind(x_comb,x_aux_bd[[k]])
    }
    y_comb <- y_target_tild
    for (k in 1:length(x_aux_bd)){
      y_comb <- c(y_comb,y_aux_tild_bd[[k]])
    }
    lambda_2 <- sqrt(2*log(nrow(x_target))/nrow(x_comb))
    w_hat <- glmnet::glmnet(x_comb,y_comb,lambda=lambda_2)$beta
    # Step 2:
    lambda_3 <- sqrt(2*log(nrow(x_target))/nrow(x_target))
    delta_hat <- glmnet::glmnet(x_target,y_target_tild-(x_target %*% w_hat),lambda=lambda_3)$beta
    return(w_hat+delta_hat)
  }
  if (mode=='real'){
    # Step 1:
    ## For target
    smoothfit0 <- smoothy(x_target,y_target,u_target)
    y_target_tild <- smoothfit0$y
    conv <- smoothfit0$converge
    if (conv!=0){
      print('We encourage a larger sample size of the target \nto facilate the sparsity condition required.')
    }
    ## For source
    y_aux_tild_bd = list()
    for (k in 1:length(x_aux_bd)){
      smoothfit1 <- smoothy(x_aux_bd[[k]],y_aux_bd[[k]],u_aux_bd[k])
      y_aux_tild_bd[[k]] <- smoothfit1$y
      conv <- smoothfit1$converge
      if (conv!=0){
        print(paste0(paste0('We encourage a larger sample size of the number',k),' source \nto facilate the sparsity condition required.'))
      }
    }
    # Step 2:
    x_comb <- x_target
    for (k in 1:length(x_aux_bd)){
      x_comb <- rbind(x_comb,x_aux_bd[[k]])
    }
    y_comb <- y_target_tild
    for (k in 1:length(x_aux_bd)){
      y_comb <- c(y_comb,y_aux_tild_bd[[k]])
    }
    fit <- glmnet::cv.glmnet(x=x_comb,y = y_comb,alpha=1)
    w_hat <- as.vector(coef(fit, s = "lambda.min"))[2:length(coef(fit, s = "lambda.min"))]
    # Step 2:
    fit <- glmnet::cv.glmnet(x=x_target,y = y_target_tild-(x_target %*% w_hat),alpha=1)
    delta_hat <- as.vector(coef(fit, s = "lambda.min"))[2:length(coef(fit, s = "lambda.min"))]
    return(w_hat+delta_hat)
  }
}

#' Fusion model
#' @description Function for building fusion model between provided target and sources without the final debias procedure.
#' @export
#' @param x_target A n*p matrix as the covariate from the target population.
#' @param y_target A vector with length n as the response from the target population.
#' @param x_aux_bd A list object that each element is a n*p matrix as the covariate from informative sources.
#' @param y_aux_bd A list object that each element is a vector with length n as the covariate from informative sources.
#' @param u_target A scalar in (0,1). The quantile level for model of target population.
#' @param u_aux_bd A vector that contains the quantile level for informative populations.
#' @param mode A string. Default is 'real'. Only 'simulation' or 'real' is allowed. The difference of mode leads to different operations for the selection of strength of penalty in single-source modeling. 'real' mode takes longer time.
#' @import quantreg
#' @import glmnet
#' @return A vector. The fusion estimator withou debias procedure based on the given auxiliary sources.
rq.fusion <- function(x_target,y_target,x_aux_bd,y_aux_bd,u_target,u_aux_bd,mode='real'){
  for (k in 1:length(x_aux_bd)){
    if (ncol(x_target)!=ncol(x_aux_bd[[k]])){
      print("Two numbers of features don't match!")
      break
    }
  }
  n <- nrow(x_target)
  if (mode=='simulation'){
    lambda_1 <- tune_param(x_target,y_target,u_target)
    # Step 1:
    ## For target
    beta_target_hat <- quantreg::rq.fit.lasso(x_target,y_target,u_target,lambda_1)$coefficients
    h_target <- sqrt(6*log(nrow(x_target))/nrow(x_target))+sqrt(6)^-1*(6^2*log(nrow(x_target))/nrow(x_target))^0.5
    f_target_hat <- sum(apply((y_target - x_target %*% beta_target_hat / h_target),1,K))/(nrow(x_target)*h_target)
    y_target_tild <- x_target %*% beta_target_hat - f_target_hat^-1 *(ifelse(y_target-x_target %*% beta_target_hat<=0,1,0)-u_target)
    ## For source
    y_aux_tild_bd <- list()
    for (k in 1:length(x_aux_bd)){
      beta_aux_hat <-  quantreg::rq.fit.lasso(x_aux_bd[[k]],y_aux_bd[[k]],u_aux_bd[k],lambda_1)$coefficients
      h_source <- sqrt(6*log(nrow(x_aux_bd[[k]]))/nrow(x_aux_bd[[k]]))+sqrt(6)^-1*(6^2*log(nrow(x_aux_bd[[k]]))/nrow(x_aux_bd[[k]]))^0.5
      f_aux_hat <- sum(apply((y_aux_bd[[k]] - x_aux_bd[[k]] %*% beta_aux_hat / h_source),1,K))/(nrow(x_aux_bd[[k]])*h_source)
      y_aux_tild_bd[[k]] <- x_aux_bd[[k]] %*% beta_aux_hat - f_aux_hat^-1 *(ifelse(y_aux_bd[[k]]-x_aux_bd[[k]] %*% beta_aux_hat<=0,1,0)-u_aux_bd[k])
    }
    # Step 2:
    x_comb <- x_target
    for (k in 1:length(x_aux_bd)){
      x_comb <- rbind(x_comb,x_aux_bd[[k]])
    }
    y_comb <- y_target_tild
    for (k in 1:length(x_aux_bd)){
      y_comb <- c(y_comb,y_aux_tild_bd[[k]])
    }
    lambda_2 <- sqrt(2*log(nrow(x_target))/nrow(x_comb))
    w_hat <- glmnet::glmnet(x_comb,y_comb,lambda=lambda_2)$beta
    return(w_hat)
  }
  if (mode=='real'){
    # Step 1:
    ## For target
    smoothfit0 <- smoothy(x_target,y_target,u_target)
    y_target_tild <- smoothfit0$y
    conv <- smoothfit0$converge
    if (conv!=0){
      print('We encourage a larger sample size of the target \nto facilate the sparsity condition required.')
    }
    ## For source
    y_aux_tild_bd = list()
    for (k in 1:length(x_aux_bd)){
      smoothfit1 <- smoothy(x_aux_bd[[k]],y_aux_bd[[k]],u_aux_bd[k])
      y_aux_tild_bd[[k]] <- smoothfit1$y
      conv <- smoothfit1$converge
      if (conv!=0){
        print(paste0(paste0('We encourage a larger sample size of the number',k),' source \nto facilate the sparsity condition required.'))
      }
    }
    # Step 2:
    x_comb <- x_target
    for (k in 1:length(x_aux_bd)){
      x_comb <- rbind(x_comb,x_aux_bd[[k]])
    }
    y_comb <- y_target_tild
    for (k in 1:length(x_aux_bd)){
      y_comb <- c(y_comb,y_aux_tild_bd[[k]])
    }
    fit <- glmnet::cv.glmnet(x=x_comb,y = y_comb,alpha=1)
    w_hat <- as.vector(coef(fit, s = "lambda.min"))[2:length(coef(fit, s = "lambda.min"))]
    return(w_hat)
  }
}

#' Quantile loss for informative sources detection
#' @description Build the auxiliary function to compute quantile loss for informative sources detection. This is different with "cdloss" function, here, we use the MSE loss on the transformed least square problem after smoothening. Not exported to users.
#' @param betahat Estimator for target only from training part.
#' @param betaic Estimator for target only from testing part.
#' @param X_measure A n*p matrix for covariate of population that is waiting for the test.
#' @param y_measure A vector as the response from the population that is waiting for the test.
#' @param u_target A scalar in (0,1). The quantile level for informative source detection testing.
#' @param hic A scalar >0. The smoothing bandwidth for calculate the loss.
#' @return A scalar. The quantile loss for given estimator on the given testing samples.
Q_loss <- function(betahat,betaic,X_measure,y_measure,u_target, hic=NULL){
  if (is.null(hic)){
    hic = smoothy(X_measure,y_measure,u_target)$h
  }
  f_hat <- sum(apply((y_measure - X_measure %*% betaic / hic),1,K))/(round(nrow(X_measure))*hic)
  y_measure_tilde <- X_measure %*% betaic - f_hat^(-1)*(ifelse(y_measure-X_measure %*% betaic<=0,1,0)-u_target)
  return(mean((y_measure_tilde-X_measure %*% betahat)^2))
}

#' Informative sources detection
#' @description The main function for informative sources detection among all datasets. This is a general function that can implement both Pseudo and normal informative sources detection, it depends on the input variables explained following.
#' @export
#' @param x_target A n*p matrix as the covariate from the target population.
#' @param y_target A vector with length n as the response from the target population.
#' @param x_aux_bd A list object that each element is a n*p matrix as the covariate from informative sources.
#' @param y_aux_bd A list object that each element is a vector with length n as the covariate from informative sources.
#' @param u_target A scalar in (0,1). The quantile level for model of target population.
#' @param u_aux_bd A vector that contains the quantile level for informative populations.
#' @param epsilon A scalar. Default is NULL. The strict level for informative sources detectiond, larger the value is, less strict the procedure is.
#' @param mode A string. Default is 'real'. Only 'simulation' or 'real' is allowed. The difference of mode leads to different operations for the selection of strength of penalty in single-source modeling. 'real' mode takes longer time.
#' @param psd A boolean variable. Default is true. Whether the procedure is pseudo.
#' @param info_num A integar. Default is NULL. The given number of informative sources under pseduo running.
#' @import quantreg
#' @import utils
#' @return A vector of index. The index of the informative sources.
info_detect <- function(x_target,y_target,x_aux_bd,y_aux_bd,u_target,u_aux_bd,mode='real',psd=FALSE,info_num=NULL,epsilon=NULL){
  if (psd==TRUE){
    if (is.null(info_num)){
      print('The argument info_num should be given!')
      return()
    }
  } else{
    if (is.null(epsilon)){
      print('The argument epsilon should be given!')
      return()
    }
  }
  if (length(x_aux_bd)!=length(y_aux_bd)){
    print('the # of datasets for x and y are not agreed!')
  } else {
    total_index <- 1:length(x_aux_bd)
    # Step 1.1: cut the target half and half: (remark: 1 for I, 2 for Ic)
    target_index <- 1:nrow(x_target)
    X0_index <- sample(target_index,size = round(nrow(x_target)/2),replace=F)
    X0_cut1 <- x_target[X0_index,]
    X0_cut2 <- x_target[-X0_index,]
    y0_cut1 <- y_target[X0_index]
    y0_cut2 <- y_target[-X0_index]
    if (mode=='simulation'){
      # Step 1.2: train sparse quantile regression for I for X0, y0:
      lambda_0 <- tune_param(X0_cut1,y0_cut1,u=u_target)
      beta_0_hat <- quantreg::rq.fit.lasso(X0_cut1,y0_cut1,tau=u_target,lambda=lambda_0)$coefficients
      # Step 1.3: train sparse quantile regression for I for X0, y0 + k-th source:
      for (k in 1:length(x_aux_bd)){
        assign(paste0(paste0('beta_',k),'_hat'), rq.fusion(X0_cut1,y0_cut1,list(get(paste0('X_',k))),list(get(paste0('y_',k))),u_target,c(u_aux_bd[k]),mode=mode))
      }
      # Step 2: Calculate the loss based on the loss function:
      beta_0_hatic <- quantreg::rq.fit.lasso(X0_cut2,y0_cut2,tau=u_target,lambda=lambda_0)$coefficients
      for (k in c(0,(1:length(x_aux_bd)))){
        assign(paste0('Q',k), Q_loss(get(paste0(paste0('beta_',k),'_hat')),beta_0_hatic,X0_cut2,y0_cut2,u_target,hic=sqrt(6*log(round(nrow(X0_cut2)))/round(nrow(X0_cut2)))+sqrt(6)^-1*(6^2*log(round(nrow(X0_cut2)))/round(nrow(X0_cut2)))^0.5))
      }
    } else {
      # Step 1.2: train sparse quantile regression for I for X0, y0:
      beta_0_hat <- quantreg::rq.fit.lasso(X0_cut1,y0_cut1,tau=u_target)$coefficients
      # Step 1.3: train sparse quantile regression for I for X0, y0 + k-th source:
      for (k in 1:length(x_aux_bd)){
        assign(paste0(paste0('beta_',k),'_hat'), rq.fusion(X0_cut1,y0_cut1,list(get(paste0('X_',k))),list(get(paste0('y_',k))),u_target,c(u_aux_bd[k]),mode=mode))
      }
      # Step 2: Calculate the loss based on the loss function:
      beta_0_hatic <- quantreg::rq.fit.lasso(X0_cut2,y0_cut2,tau=u_target)$coefficients
      for (k in c(0,(1:length(x_aux_bd)))){
        assign(paste0('Q',k), Q_loss(get(paste0(paste0('beta_',k),'_hat')),beta_0_hatic,X0_cut2,y0_cut2,u_target))
      }
    }
    # Step 3: Get informative sets:
    if (psd==FALSE){
      inform_index <- c()
      print(get(paste0('Q',0)))
      for (k in 1:length(x_aux_bd)){
        print(get(paste0('Q',k)))
        if (get(paste0('Q',k))<=(1+epsilon)*get(paste0('Q',0))){
          inform_index <- c(inform_index,k)
        } else{
          next
        }
      }
      return(inform_index)
    }
    if (psd==TRUE){
      print(get(paste0('Q',0)))
      Q_K_vec <- c()
      for (k in 1:length(x_aux_bd)){
        print(get(paste0('Q',k)))
        Q_K_vec <- c(Q_K_vec,get(paste0('Q',k)))
      }
      return(match(utils::head(sort(Q_K_vec,decreasing = F),info_num),Q_K_vec))
    }
  }
}
