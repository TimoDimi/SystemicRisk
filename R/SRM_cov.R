
vcovA <- function(SRM_object, sparsity="kernel", bandwidth="MachadoSilva",...) {

  risk_measure <- SRM_object$risk_measure
  data_trunc <- slice(SRM_object$data,-1) # Cut off first line of data as there are no predictions available for this line/date!

  TT <- dim(data_trunc)[1]
  theta_info <- theta_fun(model=SRM_object$model, theta=SRM_object$theta, df=data_trunc %>% dplyr::select(-c("VaR", "risk_measure")))
  q1 <- theta_info$length_theta1
  q2 <- theta_info$length_theta2
  y <- data_trunc$y
  x <- data_trunc$x
  v <- data_trunc$VaR
  c <- data_trunc$risk_measure
  alpha <- SRM_object$prob_level$alpha
  beta <- SRM_object$prob_level$beta

  # m <- model_fun(theta=SRM_object$theta, df=data_trunc %>% dplyr::select(-c("VaR", "risk_measure")), prob_level=SRM_object$prob_level, model=SRM_object$model, risk_measure=risk_measure)
  # v <- m$m1
  # c <- m$m2

  nabla_m <- nabla_fun(theta=SRM_object$theta, df=data_trunc %>% dplyr::select(-c("VaR", "risk_measure")), prob_level=SRM_object$prob_level, model=SRM_object$model, risk_measure=risk_measure)
  nabla_v <- nabla_m$nabla_m1[,1:q1]
  nabla_c <- nabla_m$nabla_m2[,(q1+1):(q1+q2)]

  # Xi: Outer product of the gradients
  Xi_vv_vec <- array(NA, dim=c(TT, q1, q1))
  Xi_cc_vec <- array(NA, dim=c(TT, q2, q2))
  Xi_vc_vec <- array(NA, dim=c(TT, q1, q2))
  Xi_cv_vec <- array(NA, dim=c(TT, q2, q1))
  for (tt in 1:TT){
    Xi_vv_vec[tt,,] <- outer(nabla_v[tt,], nabla_v[tt,])
    Xi_cc_vec[tt,,] <- outer(nabla_c[tt,], nabla_c[tt,])
    Xi_vc_vec[tt,,] <- outer(nabla_v[tt,], nabla_c[tt,])
    Xi_cv_vec[tt,,] <- t(Xi_vc_vec[tt,,])
  }

  # Estimate the conditional distributions/densities
  if (sparsity=="kernel"){
    cond_dist <- sparsity_kernel(x=x, y=y, v=v, c=c, beta=beta, alpha=alpha, TT=TT, bandwidth=bandwidth)
  } else if (sparsity=="ARB"){
    cond_dist <- sparsity_ARB(x=x, y=y, v=v, c=c, beta=beta, alpha=alpha, TT=TT, nabla_v=nabla_v, Xi_vv_vec=Xi_vv_vec)
  } else {
    error("Enter a correct value for 'sparsity'.")
  }

  fX_est <- cond_dist$fX_est
  F_delta1_est <- cond_dist$F_delta1_est
  fY_est <- cond_dist$fY_est
  F_delta2_est <- cond_dist$F_delta2_est

  if (risk_measure=="CoVaR"){
    # V, C and D Matrix
    V <- beta*(1-beta) * apply(Xi_vv_vec, c(2,3), mean, na.rm=TRUE)
    D <- alpha*(1-alpha)*(1-beta) * apply(Xi_cc_vec, c(2,3), mean, na.rm=TRUE)
    C <- rbind(cbind(V, matrix(0, nrow=q1, ncol=q2)),
               cbind(matrix(0, nrow=q2, ncol=q1), D))

    # Lambda matrices
    Lambda_vec <- Xi_vv_vec * array(fX_est, dim=c(TT,q1,q1))
    Lambda1_vec <- Xi_cc_vec * array(fY_est - F_delta2_est, dim=c(TT,q2,q2))
    Lambda2_vec <- Xi_cv_vec * array(alpha*fX_est - F_delta1_est, dim=c(TT,q2,q1))
  } else if (risk_measure=="MES"){
    # V, C and M Matrix
    V <- beta*(1-beta) * apply(Xi_vv_vec, c(2,3), mean, na.rm=TRUE)
    M_YH <- apply( array(((x>v) * y^2 - (1-beta) * c^2), dim=c(TT,q2,q2)) * Xi_cc_vec, c(2,3), mean, na.rm=TRUE)  # IS THIS CORRECT???
    M <- apply( array(((x>v) * (y-c)^2), dim=c(TT,q2,q2)) * Xi_cc_vec, c(2,3), mean, na.rm=TRUE)

    C <- rbind(cbind(V, matrix(0, nrow=q1, ncol=q2)),
               cbind(matrix(0, nrow=q2, ncol=q1), M))

    # Lambda matrices
    Lambda_vec <- Xi_vv_vec * array(fX_est, dim=c(TT,q1,q1))
    Lambda1_vec <- Xi_cc_vec * array((1-beta), dim=c(TT,q2,q2))
    Lambda2_vec <- Xi_cv_vec * array((y-c)*fX_est, dim=c(TT,q2,q1))
  }

  Lambda <- apply(Lambda_vec, c(2,3), mean, na.rm=TRUE)
  Lambda1 <- apply(Lambda1_vec, c(2,3), mean, na.rm=TRUE)
  Lambda2 <- apply(Lambda2_vec, c(2,3), mean, na.rm=TRUE)
  Lambda_inv <- tryCatch(solve(Lambda), error=function(e){MASS::ginv(Lambda)})
  Lambda1_inv <- tryCatch(solve(Lambda1), error=function(e){MASS::ginv(Lambda1)})


  # Partial covariance matrices:
  Gamma <-  cbind(Lambda1_inv %*% Lambda2 %*% Lambda_inv,  -Lambda1_inv)
  cov_theta_v <- 1/TT * (Lambda_inv %*% V %*% Lambda_inv)
  cov_theta_c <- 1/TT * (Gamma %*% C %*% t(Gamma))

  # This is only the block-diagonal covariance matrix!
  cov_BlockDiag <- rbind(cbind(cov_theta_v, matrix(0, nrow=q1, ncol=q2)),
                         cbind(matrix(0, nrow=q2, ncol=q1), cov_theta_c))


  # Asymptotic covariance from Remark 4
  barGamma <- rbind(cbind(-Lambda_inv,  matrix(0, nrow=q1, ncol=q2)),
                    cbind(Lambda1_inv %*% Lambda2 %*% Lambda_inv,  -Lambda1_inv))
  cov_theta <- 1/TT * (barGamma %*% C %*% t(barGamma))

  # rownames(cov_theta_v) <- colnames(cov_theta_v) <- names(stats::coef(object))

  cov_theta
}



sparsity_kernel <- function(x, y, v, c, beta, alpha, TT, bandwidth){
  if (bandwidth=="MachadoSilva"){
    # Koenker (2005), Machado and Silva (2013) bandwidth selection
    m_T <- function(tau){(qnorm(0.975))^(2/3) * ((1.5*(dnorm(qnorm(tau)))^2)/(2*(qnorm(tau))^2 + 1))^(1/3)}
    eps <- 0.0001
    b_T <- stats::mad(x-v) * TT^(-1/3) * (qnorm( min(1-eps, beta + m_T(beta))) - qnorm( max(eps, beta - m_T(beta))))
    c_T <- stats::mad(y-c) * ((1-beta)*TT)^(-1/3) * (qnorm( min(1-eps, alpha + m_T(alpha))) - qnorm( max(eps, alpha - m_T(alpha))))
  } else if (bandwidth=="own"){
    b_T <- 0.1/sqrt(beta*(1-beta)) * stats::mad(x-v) * TT^(-1/3)
    # Larger bandwidth as c is (at least for beta >= 0.5) further in the tail of y than v of x
    c_T <- 0.1/sqrt((1-beta) * alpha * (1-alpha)) * stats::mad(y-c) * TT^(-1/3)
  } else {
    b_T <- c_T <- TT^{-1/3}
  }

  fX_est <- 1/(2*b_T) * (abs(x-v) <= b_T)
  F_delta1_est <- 1/(2*b_T) * (abs(x-v) <= b_T) * (y<=c)

  fY_est <- 1/(2*c_T) * (abs(y-c) <= c_T)
  F_delta2_est <- 1/(2*c_T) * (abs(y-c) <= c_T) * (x<=v)

  return(list(fX_est=fX_est, F_delta1_est=F_delta1_est, fY_est=fY_est, F_delta2_est=F_delta2_est))
}


sparsity_ARB <- function(x, y, v, c, beta, alpha, TT, nabla_v, Xi_vv_vec){
  # Adaptiv Random Bandwidth (ARB) method
  # So far, this only works for f_X(v), but not for the other partial derivatives
  # Should be generalized for the other derivatives of the CDF...
  cond_dist <- sparsity_kernel(x=x, y=y, v=v, c=c, beta=beta, alpha=alpha, TT=TT, bandwidth="MachadoSilva")
  V <- beta*(1-beta) * apply(Xi_vv_vec, c(2,3), mean, na.rm=TRUE)

  n <- 1000
  m_update <- 3
  q1 <- dim(nabla_v)[2]
  Vd <- 10*diag(q1)
  for (m in 1:m_update){
    fX_est <- sapply(1:n, function(i){
      rand_bw <- as.numeric(nabla_v %*% t(mvtnorm::rmvnorm(n=1, sigma=Vd)/sqrt(TT)))
      -((x <= v) - (x <= v+rand_bw)) / rand_bw
    }) %>% rowMeans()

    Lambda <- apply(Xi_vv_vec * array(fX_est, dim=c(TT,q1,q1)), c(2,3), mean, na.rm=TRUE)
    Lambda_inv <- solve(Lambda)
    Vd <- Lambda_inv %*% V %*% Lambda_inv
  }

  cond_dist$fX_est <- fX_est
  return(cond_dist)
}



# Helper function that samples moving block bootstrap (MBB) indices
MBB_indices <- function(l, bl){
  if (bl > l) bl <- l

  Index_set <- 1:(l-bl+1)
  k <- ceiling(l/bl)
  B_l <- sample(Index_set, k, replace=TRUE) # Sample the block starting points
  B_Index <- sapply(B_l, function(i) i:(i+bl-1)) %>% c() %>% .[1:l] # Generate the entire series of
  B_Index
}



vcovB <- function(SRM_object, block_length='default', B=500) {
  if (B < 500)
    warning("The number of bootstrap iterations is small!")

  data <- dplyr::select(SRM_object$data, -c("VaR", "risk_measure"))

  # Draw the bootstrap indices
  TT <- dim(SRM_object$data)[1]
  if (block_length=='default') block_length <- floor(TT^(1/3))

  # Split the starting value (the original parameter estimate)
  theta0_list <- theta_fun(model=SRM_object$model, theta=SRM_object$theta, df=data)
  theta01 <- theta0_list$theta1
  theta02 <- theta0_list$theta2

  theta_B <- matrix(NA, nrow=B, ncol=theta0_list$length_theta)

  for (b in 1:B) {
    B_idx <- MBB_indices(TT, block_length)

    # Bootstrap the loss functions according to Goncalves et al (2022, JBES)
    # Moving Block Bootstrapped VaR Optimization:
    # Use the bootstrap indices for the result of the "loss_model_VaR" function
    opt_v <- optim(par=theta01,
                   fn=function(theta,...){
                     mean(loss_model_VaR(theta,...) %>% .[B_idx], na.rm=TRUE)}, # Does the "[B_idx]" work?
                   df=data, model=SRM_object$model, prob_level=SRM_object$prob_level)

    theta_est_v <- opt_v$par
    m1_est <- model_fun(theta_est_v, df=data, prob_level=SRM_object$prob_level, model=SRM_object$model, model_type="first")$m1

    # Second bootstrap step:
    # Important: Only select the bootstrap losses, instead of directly resamlping the data!
    opt_m2 <- optim(par=theta02,
                    fn=function(theta,...){
                      switch(SRM_object$risk_measure,
                             MES = {mean(loss_model_MES(theta,...) %>% .[B_idx],na.rm=TRUE)},
                             CoVaR = {mean(loss_model_CoVaR(theta,...) %>% .[B_idx],na.rm=TRUE)})},
                    df=data, m1=m1_est, model=SRM_object$model, prob_level=SRM_object$prob_level)
    theta_est_m2 <- opt_m2$par

    theta_B[b,] <- c(theta_est_v, theta_est_m2)
  }


  cov_B <- cov(theta_B)
  # rownames(cov_B) <- colnames(cov_B) <- names(stats::coef(object))
  cov_B
}


# vcovBfast <- function(SRM_object, block_length='default', B=100) {
#   if (B < 500)
#     warning("The number of bootstrap iterations is small!")
#
#   # Draw the bootstrap indices
#   TT <- dim(SRM_object$data)[1]
#   if (block_length=='default') {block_length <- floor(TT^(1/3))}
#
#   # # Split the starting value (the original parameter estimate)
#   theta0_list <- theta_fun(model=SRM_object$model, theta=SRM_object$theta, df=SRM_object$data)
#   # theta01 <- theta0_list$theta1
#   # theta02 <- theta0_list$theta2
#
#   TT <- dim(SRM_object$data)[1]
#   theta_info <- theta_fun(model=SRM_object$model, theta=SRM_object$theta)
#   q1 <- theta_info$length_theta1
#   q2 <- theta_info$length_theta2
#   y <- SRM_object$data$y
#   x <- SRM_object$data$x
#   alpha = SRM_object$prob_level$alpha
#   beta = SRM_object$prob_level$beta
#
#   m <- model_fun(theta=SRM_object$theta, df=SRM_object$data, prob_level=SRM_object$prob_level, model=SRM_object$model, risk_measure="CoVaR")
#   v <- m$m1
#   c <- m$m2
#
#   nabla_m <- nabla_fun(theta=SRM_object$theta, df=SRM_object$data, prob_level=SRM_object$prob_level, model=SRM_object$model, risk_measure="CoVaR")
#   nabla_v <- nabla_m$nabla_m1[,1:q1]
#   nabla_c <- nabla_m$nabla_m2[,(q1+1):(q1+q2)]
#
#   # Xi: Outer product of the gradients
#   Xi_vv_vec <- array(NA, dim=c(TT, q1, q1))
#   Xi_vc_vec <- array(NA, dim=c(TT, q1, q2))
#   Xi_cc_vec <- array(NA, dim=c(TT, q2, q2))
#   for (tt in 1:TT){
#     Xi_vv_vec[tt,,] <- outer(nabla_v[tt,], nabla_v[tt,]) # t(nabla_v[tt,]) %*% nabla_v[tt,]
#     Xi_vc_vec[tt,,] <- outer(nabla_v[tt,], nabla_c[tt,]) # t(nabla_v[tt,]) %*% nabla_c[tt,]
#     Xi_cc_vec[tt,,] <- outer(nabla_c[tt,], nabla_c[tt,]) # t(nabla_c[tt,]) %*% nabla_c[tt,]
#   }
#
#   ### Estimate the nuisance quantities for the Lambdas through indicator functions
#   b_T <- TT^(-1/3)
#   fX_est <- 1/(2*b_T) * (abs(x-v) < b_T)
#   F_delta1_est <- 1/(2*b_T) * (abs(x-v) < b_T) * (y<=c)
#
#   # Larger bandwidth as c is further in the tail of y than v of x
#   c_T <- 5*TT^(-1/3)
#   fY_est <- 1/(2*c_T) * (abs(y-c) < c_T)
#   F_delta2_est <- 1/(2*c_T) * (abs(y-c) < c_T) * (x<=v)
#
#   Lambda_vec <- Xi_vv_vec * array(fX_est, dim=c(TT,q1,q1))
#   Lambda1_vec <- Xi_cc_vec * array(fY_est - F_delta2_est, dim=c(TT,q2,q2))
#   Lambda2_vec <- Xi_vc_vec * array(alpha*fX_est - F_delta1_est, dim=c(TT,q1,q2))
#
#   Lambda <- apply(Lambda_vec, c(2,3), mean, na.rm=TRUE)
#   Lambda1 <- apply(Lambda1_vec, c(2,3), mean, na.rm=TRUE)
#   Lambda2 <- apply(Lambda2_vec, c(2,3), mean, na.rm=TRUE)
#   Lambda_inv <- tryCatch(solve(Lambda), error=function(e){MASS::ginv(Lambda)})
#   Lambda1_inv <- tryCatch(solve(Lambda1), error=function(e){MASS::ginv(Lambda1)})
#
#
#   # Compute the moment conditions which we resample from below
#   A <- abind::abind(nabla_m$nabla_m1, nabla_m$nabla_m2, along=3)
#   momcond <- momcond_CoVaR(SRM_object$theta, SRM_object$data, A, SRM_object$model, SRM_object$prob_level)
#   momcond_trans <- cbind(momcond[,1:q1],
#                          t(t(momcond[,(q1+1):(q1+q2)]) - Lambda2 %*% Lambda_inv %*% t(momcond[,1:q1])))
#
#
#   theta_B <- matrix(NA, nrow=B, ncol=theta0_list$length_theta)
#   for (b in 1:B) {
#     B_idx <- MBB_indices(TT, block_length)
#     theta_B[b,] <- SRM_object$theta -
#       rbind(cbind(Lambda_inv, matrix(0, nrow=q1, ncol=q2)),
#             cbind(matrix(0, nrow=q2, ncol=q1), Lambda1_inv)) %*%
#       colMeans(momcond_trans[B_idx,], na.rm=TRUE)
#   }
#
#   cov_B <- cov(theta_B)
#   # rownames(cov_B) <- colnames(cov_B) <- names(stats::coef(object))
#   cov_B
# }






