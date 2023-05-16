
################################################################################
###    General model, nabla and split_theta functions
################################################################################

# A model should have three components, i.e., it is a list of length 3 or a (pre-implemented) string
# - m(theta)
# - nabla_m(theta)
# - theta_split which splits the parameter vector into theta1 and theta2 of length q1 and q2
# Simply call theta_split(theta)$theta1 and theta_split(theta)$theta2 to get the subvectors


#' SRM Model Function
#'
#' Computes the joint VaR/CoVaR model depending on either a cross-sectional, CoCAViaR or user-written model type. See also \link{SRM}
#'
#' @details
#' 'model_fun' computes the model, 'nabla_fun' its gradient, and 'theta_fun' some generic information about the respective model parameters.
#' Currently implemented are the following models:
#' \itemize{
#'   \item "joint_linear": Linear models with covariates z for both, the VaR and CoVaR
#'   \item "CoCAViaR_SAV_diag": Symmetric Absolute Value (SAV) CoCAViaR model with diagonal matrices A and B
#'   \item "CoCAViaR_SAV_fullA": Symmetric Absolute Value (SAV) CoCAViaR model with full matrix A and and diagonal B
#'   \item "CoCAViaR_SAV_full": Symmetric Absolute Value (SAV) CoCAViaR model with full matrix A and and upper triangular matrix B
#'   \item "CoCAViaR_AS_pos": Asymmetric Slope (AS) CoCAViaR model where only the positive losses (negative returns) act as forcing variables
#'   \item "CoCAViaR_AS_signs": Asymmetric Slope (AS) CoCAViaR model where only the positive and negative losses act as different forcing variables
#'   \item "CoCAViaR_AS_mixed": Asymmetric Slope (AS) CoCAViaR model with 'mixed' forcing variables
#' }
#'
#' Own model functions can be written as follows:
#' ToDo: Specify details!
#'
#' @param theta Parameter vector theta
#' @param df data.frame/tibble containing the collected data
#' @param prob_level list of probability levels beta and alpha
#' @param model string or list indicating the model, see Details for more
#' @param risk_measure The system risk measure; currently only implemented for the "CoVaR"
#' @param forecast logical variable whether the forecast should be produced in the call or not
#' @param model_type either "joint", "first" or "second" indicating which model (both, VaR, or CoVaR) should be returned
#' @param m1 Fitted values for the first (VaR) model; required for model_type="second"
#'
#' @return a list with entries m1 and m2 containing the VaR and CoVaR values as entries
#' @export
model_fun <- function(theta, df, prob_level, model="joint_linear", risk_measure="CoVaR", forecast=FALSE, model_type="joint", m1=NULL){
  if (!is.list(model)) {
    m <- do.call(paste0("model_", model), args=list(theta=theta, df=df, prob_level=prob_level, risk_measure=risk_measure, forecast=forecast, model_type=model_type, m1=m1))
  } else {
    m <- model$fun(theta=theta, df=df, prob_level=prob_level, risk_measure=risk_measure, forecast=forecast, model_type=model_type, m1=m1)
  }
  return(m)
}


nabla_fun <- function(theta, df, prob_level, model="joint_linear", risk_measure="CoVaR"){
  if (!is.list(model)) {
    nabla_m <- do.call(paste0("nabla_", model), args=list(theta=theta, df=df, prob_level=prob_level, risk_measure=risk_measure))
  } else {
    nabla_m <- model$nabla(theta, df)
  }
  nabla_m
}


theta_fun <- function(model, theta=NULL, df=NULL){
  if (!is.list(model)) {
    theta_list <- do.call(paste0("theta_", model), args=list(theta=theta, df=df))
  } else {
    theta_list <- model$theta_fun(theta, df)
  }
  theta_list
}




################################################################################
###    Joint linear model functions
################################################################################

# Careful: I changed it such that the covariates z[tt-1] "belong to" x[tt] and y[tt],
# i.e., they are in the row above!!!!!
#

model_joint_linear <- function(theta, df, prob_level, risk_measure, forecast=FALSE, model_type="joint", m1=NULL){

  TT <- dim(df)[1]
  if (forecast){TT <- TT + 1} # To produce a model forecast, simply add 1 to TT

  # Covariates
  z <- df %>%
    tibble::as_tibble() %>%
    dplyr::select(-c("Date", "x", "y")) %>%
    as.matrix()

  # If forecasts are produced, we need the last line of covariates, otherwise we don't!
  if (forecast){z <- rbind(lag(z),tail(z,1))}
  else{z <- lag(z) } # to capture the time lag between x[tt], y[tt] and z[tt-1]!

  q <- length(theta)
  if (model_type=="joint"){
    m1 <- z %*% theta[1:(q/2)]
    m2 <- z %*% theta[(q/2+1):q]
    return(list(m1=m1, m2=m2))
  } else if (model_type=="first"){
    m1 <- z %*% theta
    return(list(m1=m1))
  } else if (model_type=="second"){
    m2 <- z %*% theta
    return(list(m1=m1, m2=m2))
  }
}



#' Gradient of the joint linear VaR/CoVaR model
#'
#' @param theta Parameter vector
#' @param df Data
#' @param prob_level List of the probability levels alpha and beta
#' @param risk_measure Systemic risk measure, currently only "CoVaR"
#'
#' @return list of model gradients
nabla_joint_linear <- function(theta, df, prob_level, risk_measure){
  z <- df %>%
    tibble::as_tibble() %>%
    dplyr::select(-c("Date", "x", "y")) %>%
    as.matrix() %>%
    lag() # to capture the time lag between x[tt], y[tt] and z[tt-1]!

  nabla_m1 <- cbind(z, array(0,dim=dim(z)))
  nabla_m2 <- cbind(array(0,dim=dim(z)), z)

  list(nabla_m1=nabla_m1, nabla_m2=nabla_m2)
}



theta_joint_linear <- function(theta, df){
  # We get a result for q when either theta or df are specified
  if (!is.null(df)){
    # Automatically get the dimension
    q <- 2 * df %>%
      tibble::as_tibble() %>%
      dplyr::select(-c("Date", "x", "y")) %>%
      dim() %>%
      .[2]
  } else {
    q <- length(theta)
    # q <- length(theta)/2
  }

  if (!(length(theta)%% 2) == 0){stop("Length of theta is not even!")}
  theta1 <- theta[1:(q/2)]
  theta2 <- theta[-(1:(q/2))]
  theta_start_default <- rep(0,q)
  return(list(theta1=theta1, theta2=theta2, length_theta=q, length_theta1=q/2, length_theta2=q/2, theta_start_default=theta_start_default, theta_names=NULL))
}






################################################################################
###   6 Parameter "CoCAViaR_SAV_diag" model of the symmetric absolute value class with diagonal A and B matrices.
###   This model more or less corresponds to an absolute value CCC-GARCH model with diagonal A and B.
################################################################################

model_CoCAViaR_SAV_diag <- function(theta, df, prob_level, risk_measure, forecast=FALSE, model_type="joint", m1=NULL){
  q <- length(theta)
  alpha <- prob_level$alpha
  beta <- prob_level$beta

  TT <- dim(df)[1]
  if (forecast){TT <- TT + 1} # To produce a model forecast, simply add 1 to TT

  if (model_type=="joint"){
    m1 <- m2 <- rep(NA,TT)
    m1[1] <- quantile(df$x, beta)
    # Initialize m2 depending on the risk_measure we estimate
    if (risk_measure=="CoVaR"){
      m2[1] <- quantile(df$y[df$x >= m1[1]] , alpha)
    } else if (risk_measure=="MES"){
      m2[1] <- mean(df$y[df$x >= m1[1]])
    }
    for (tt in 2:TT){
      m1[tt] <- theta[1] + theta[2]*abs(df$x[tt-1]) + theta[3]*m1[tt-1]
      m2[tt] <- theta[4] + theta[5]*abs(df$y[tt-1]) + theta[6]*m2[tt-1]
    }
    return(list(m1=m1, m2=m2))
  } else if (model_type=="first"){
    m1 <- rep(NA,TT)
    m1[1] <- quantile(df$x, beta)
    for (tt in 2:TT){
      m1[tt] <- theta[1] + theta[2]*abs(df$x[tt-1]) + theta[3]*m1[tt-1]
    }
    return(list(m1=m1))
  } else if (model_type=="second"){
    stopifnot(!is.null(m1)) # A vector of m1 must be given!
    m2 <- rep(NA,TT)
    # Initialize m2 depending on the risk_measure we estimate
    if (risk_measure=="CoVaR"){
      m2[1] <- quantile(df$y[df$x >= m1[1]] , alpha)
    } else if (risk_measure=="MES"){
      m2[1] <- mean(df$y[df$x >= m1[1]])
    }
    for (tt in 2:TT){
      m2[tt] <- theta[1] + theta[2]*abs(df$y[tt-1]) + theta[3]*m2[tt-1]
    }
    return(list(m1=m1, m2=m2))
  }
}


nabla_CoCAViaR_SAV_diag <- function(theta, df, prob_level, risk_measure){
  TT <- dim(df)[1]
  m <- model_CoCAViaR_SAV_diag(theta, df, prob_level, risk_measure)

  # Manually implement the model gradient
  nabla_m1 <- cbind(1, lag(abs(df$x)), lag(m$m1), array(0, dim=c(TT,3)))
  nabla_m2 <- cbind(array(0, dim=c(TT,3)), 1, lag(abs(df$y)), lag(m$m2))

  list(nabla_m1=nabla_m1, nabla_m2=nabla_m2)
}


theta_CoCAViaR_SAV_diag <- function(theta, df=NULL){
  length_theta <- 6
  length_theta1 <- 3
  length_theta2 <- 3
  theta1 <- theta[1:length_theta1]
  theta2 <- theta[(length_theta1+1):(length_theta1+length_theta2)]
  theta_start_default <- c(0.01, 0.1,      0.75,
                           0.01,      0.1,      0.75)
  theta_names <- list(VaR=c("(Intercept)", "lag |X|", "lag VaR"),
                      CoVaR=c("(Intercept)", "lag |Y|", "lag CoVaR"),
                      MES=c("(Intercept)", "lag |Y|", "lag MES"))
  return(list(theta1=theta1, theta2=theta2, length_theta=length_theta, length_theta1=length_theta1, length_theta2=length_theta2, theta_start_default=theta_start_default, theta_names=theta_names))
}



################################################################################
###    8 Parameter "CoCAViaR_SAV_fullA" model with the restriction B_{12} = B_{21} = 0
################################################################################

model_CoCAViaR_SAV_fullA <- function(theta, df, prob_level, risk_measure, forecast=FALSE, model_type="joint", m1=NULL){
  q <- length(theta)
  alpha <- prob_level$alpha
  beta <- prob_level$beta

  TT <- dim(df)[1]
  if (forecast){TT <- TT + 1} # To produce a model forecast, simply add 1 to TT

  if (model_type=="joint"){
    m1 <- m2 <- rep(NA,TT)
    m1[1] <- quantile(df$x, beta)
    # Initialize m2 depending on the risk_measure we estimate
    if (risk_measure=="CoVaR"){
      m2[1] <- quantile(df$y[df$x >= m1[1]] , alpha)
    } else if (risk_measure=="MES"){
      m2[1] <- mean(df$y[df$x >= m1[1]])
    }
    for (tt in 2:TT){
      m1[tt] <- theta[1] + theta[2]*abs(df$x[tt-1]) + theta[3]*abs(df$y[tt-1]) + theta[4]*m1[tt-1]
      m2[tt] <- theta[5] + theta[6]*abs(df$x[tt-1]) + theta[7]*abs(df$y[tt-1]) + theta[8]*m2[tt-1]
    }
    return(list(m1=m1, m2=m2))
  } else if (model_type=="first"){
    m1 <- rep(NA,TT)
    m1[1] <- quantile(df$x, beta)
    for (tt in 2:TT){
      m1[tt] <- theta[1] + theta[2]*abs(df$x[tt-1]) + theta[3]*abs(df$y[tt-1]) + theta[4]*m1[tt-1]
    }
    return(list(m1=m1))
  } else if (model_type=="second"){
    stopifnot(!is.null(m1)) # A vector of m1 must be given!
    m2 <- rep(NA,TT)
    # Initialize m2 depending on the risk_measure we estimate
    if (risk_measure=="CoVaR"){
      m2[1] <- quantile(df$y[df$x >= m1[1]] , alpha)
    } else if (risk_measure=="MES"){
      m2[1] <- mean(df$y[df$x >= m1[1]])
    }
    for (tt in 2:TT){
      m2[tt] <- theta[1] + theta[2]*abs(df$x[tt-1]) + theta[3]*abs(df$y[tt-1]) + theta[4]*m2[tt-1]
    }
    return(list(m1=m1, m2=m2))
  }
}


nabla_CoCAViaR_SAV_fullA <- function(theta, df, prob_level, risk_measure){
  q <- length(theta)
  TT <- dim(df)[1]

  m <- model_CoCAViaR_SAV_fullA(theta, df, prob_level, risk_measure)

  # Manually implement the gradient
  nabla_m1 <- cbind(1, lag(abs(df$x)), lag(abs(df$y)), lag(m$m1), array(0, dim=c(TT,4)))
  nabla_m2 <- cbind(array(0, dim=c(TT,4)), 1, lag(abs(df$x)), lag(abs(df$y)), lag(m$m2))

  list(nabla_m1=nabla_m1, nabla_m2=nabla_m2)
}


theta_CoCAViaR_SAV_fullA <- function(theta, df=NULL){
  length_theta <- 8
  length_theta1 <- 4
  length_theta2 <- 4
  theta1 <- theta[1:length_theta1]
  theta2 <- theta[(length_theta1+1):(length_theta1+length_theta2)]
  theta_start_default <- c(0.01, 0.1, 0.01, 0.8,
                           0.01, 0.01, 0.1,      0.8)
  theta_names <- list(VaR=c("(Intercept)", "lag |X|",  "lag |Y|", "lag VaR"),
                      CoVaR=c("(Intercept)", "lag |X|",  "lag |Y|", "lag CoVaR"),
                      MES=c("(Intercept)", "lag |X|",  "lag |Y|", "lag MES"))
  return(list(theta1=theta1, theta2=theta2, length_theta=length_theta, length_theta1=length_theta1, length_theta2=length_theta2, theta_start_default=theta_start_default, theta_names=theta_names))
}





################################################################################
###    9 Parameter "CoCAViaR_SAV_full" model with the restriction: B_{21}=0
################################################################################

model_CoCAViaR_SAV_full <- function(theta, df, prob_level, risk_measure, forecast=FALSE, model_type="joint", m1=NULL){
  q <- length(theta)
  alpha <- prob_level$alpha
  beta <- prob_level$beta

  TT <- dim(df)[1]
  if (forecast){TT <- TT + 1} # To produce a model forecast, simply add 1 to TT

  if (model_type=="joint"){
    m1 <- m2 <- rep(NA,TT)
    m1[1] <- quantile(df$x, beta)
    # Initialize m2 depending on the risk_measure we estimate
    if (risk_measure=="CoVaR"){
      m2[1] <- quantile(df$y[df$x >= m1[1]] , alpha)
    } else if (risk_measure=="MES"){
      m2[1] <- mean(df$y[df$x >= m1[1]])
    }
    for (tt in 2:TT){
      m1[tt] <- theta[1] + theta[2]*abs(df$x[tt-1]) + theta[3]*abs(df$y[tt-1])  + theta[4]*m1[tt-1]
      m2[tt] <- theta[5] + theta[6]*abs(df$x[tt-1]) + theta[7]*abs(df$y[tt-1]) + theta[8]*m1[tt-1] + theta[9]*m2[tt-1]
    }
    return(list(m1=m1, m2=m2))
  } else if (model_type=="first"){
    m1 <- rep(NA,TT)
    m1[1] <- quantile(df$x, beta)
    for (tt in 2:TT){
      m1[tt] <- theta[1] + theta[2]*abs(df$x[tt-1]) + theta[3]*abs(df$y[tt-1]) + theta[4]*m1[tt-1]
    }
    return(list(m1=m1))
  } else if (model_type=="second"){
    stopifnot(!is.null(m1)) # A vector of m1 must be given!
    m2 <- rep(NA,TT)
    # Initialize m2 depending on the risk_measure we estimate
    if (risk_measure=="CoVaR"){
      m2[1] <- quantile(df$y[df$x >= m1[1]] , alpha)
    } else if (risk_measure=="MES"){
      m2[1] <- mean(df$y[df$x >= m1[1]])
    }
    for (tt in 2:TT){
      m2[tt] <- theta[1] + theta[2]*abs(df$x[tt-1]) + theta[3]*abs(df$y[tt-1]) + theta[4]*m1[tt-1] + theta[5]*m2[tt-1]
    }
    return(list(m1=m1, m2=m2))
  }
}


nabla_CoCAViaR_SAV_full <- function(theta, df, prob_level, risk_measure){
  q <- length(theta)
  TT <- dim(df)[1]

  m <- model_CoCAViaR_SAV_full(theta, df, prob_level, risk_measure)

  # Manually implement the gradient
  nabla_m1 <- cbind(1, lag(abs(df$x)), lag(abs(df$y)), lag(m$m1), array(0, dim=c(TT,5)))
  nabla_m2 <- cbind(array(0, dim=c(TT,4)), 1, lag(abs(df$x)), lag(abs(df$y)), lag(m$m1), lag(m$m2))

  list(nabla_m1=nabla_m1, nabla_m2=nabla_m2)
}


theta_CoCAViaR_SAV_full <- function(theta, df=NULL){
  length_theta <- 9
  length_theta1 <- 4
  length_theta2 <- 5
  theta1 <- theta[1:length_theta1]
  theta2 <- theta[(length_theta1+1):(length_theta1+length_theta2)]
  theta_start_default <- c(0.01, 0.1, 0.01, 0.8,
                           0.01, 0.01, 0.1, 0.05, 0.8)
  theta_names <- list(VaR=c("(Intercept)", "lag |X|",  "lag |Y|", "lag VaR"),
                      CoVaR=c("(Intercept)", "lag |X|",  "lag |Y|", "lag VaR", "lag CoVaR"),
                      CoVaR=c("(Intercept)", "lag |X|",  "lag |Y|", "lag VaR", "lag MES"))
  return(list(theta1=theta1, theta2=theta2, length_theta=length_theta, length_theta1=length_theta1, length_theta2=length_theta2, theta_start_default=theta_start_default, theta_names=theta_names))
}









################################################################################
###    8 Parameter asymmetric slopde "CoCAViaR_AS_pos" model with the restriction B_{12} = B_{21} = 0 and only the positive componentes of X_t and Y_t enter the equations
################################################################################

model_CoCAViaR_AS_pos <- function(theta, df, prob_level, risk_measure, forecast=FALSE, model_type="joint", m1=NULL){
  q <- length(theta)
  alpha <- prob_level$alpha
  beta <- prob_level$beta

  TT <- dim(df)[1]
  if (forecast){TT <- TT + 1} # To produce a model forecast, simply add 1 to TT

  if (model_type=="joint"){
    m1 <- m2 <- rep(NA,TT)
    m1[1] <- quantile(df$x, beta)
    # Initialize m2 depending on the risk_measure we estimate
    if (risk_measure=="CoVaR"){
      m2[1] <- quantile(df$y[df$x >= m1[1]] , alpha)
    } else if (risk_measure=="MES"){
      m2[1] <- mean(df$y[df$x >= m1[1]])
    }
    for (tt in 2:TT){
      m1[tt] <- theta[1] + theta[2]*pmax(df$x[tt-1],0) + theta[3]*pmax(df$y[tt-1],0) + theta[4]*m1[tt-1]
      m2[tt] <- theta[5] + theta[6]*pmax(df$x[tt-1],0) + theta[7]*pmax(df$y[tt-1],0) + theta[8]*m2[tt-1]
    }
    return(list(m1=m1, m2=m2))
  } else if (model_type=="first"){
    m1 <- rep(NA,TT)
    m1[1] <- quantile(df$x, beta)
    for (tt in 2:TT){
      m1[tt] <- theta[1] + theta[2]*pmax(df$x[tt-1],0) + theta[3]*pmax(df$y[tt-1],0) + theta[4]*m1[tt-1]
    }
    return(list(m1=m1))
  } else if (model_type=="second"){
    stopifnot(!is.null(m1)) # A vector of m1 must be given!
    m2 <- rep(NA,TT)
    # Initialize m2 depending on the risk_measure we estimate
    if (risk_measure=="CoVaR"){
      m2[1] <- quantile(df$y[df$x >= m1[1]] , alpha)
    } else if (risk_measure=="MES"){
      m2[1] <- mean(df$y[df$x >= m1[1]])
    }
    for (tt in 2:TT){
      m2[tt] <- theta[1] + theta[2]*pmax(df$x[tt-1],0) + theta[3]*pmax(df$y[tt-1],0) + theta[4]*m2[tt-1]
    }
    return(list(m1=m1, m2=m2))
  }
}


nabla_CoCAViaR_AS_pos <- function(theta, df, prob_level, risk_measure){
  q <- length(theta)
  TT <- dim(df)[1]

  m <- model_CoCAViaR_AS_pos(theta, df, prob_level, risk_measure)

  # Manually implement the gradient
  nabla_m1 <- cbind(1, lag(pmax(df$x,0)), lag(pmax(df$y,0)), lag(m$m1), array(0, dim=c(TT,4)))
  nabla_m2 <- cbind(array(0, dim=c(TT,4)), 1, lag(pmax(df$x,0)), lag(pmax(df$y,0)), lag(m$m2))

  list(nabla_m1=nabla_m1, nabla_m2=nabla_m2)
}


theta_CoCAViaR_AS_pos <- function(theta, df=NULL){
  length_theta <- 8
  length_theta1 <- 4
  length_theta2 <- 4
  theta1 <- theta[1:length_theta1]
  theta2 <- theta[(length_theta1+1):(length_theta1+length_theta2)]
  theta_start_default <- c(0.01, 0.1, 0.05, 0.75,
                           0.01, 0.05, 1,       0.75)
  theta_names <- list(VaR=c("(Intercept)", "lag X+", "lag Y+", "lag VaR"),
                      CoVaR=c("(Intercept)", "lag X+","lag Y+", "lag CoVaR"),
                      MES=c("(Intercept)", "lag X+","lag Y+", "lag MES"),)
  return(list(theta1=theta1, theta2=theta2, length_theta=length_theta, length_theta1=length_theta1, length_theta2=length_theta2, theta_start_default=theta_start_default, theta_names=theta_names))
}



################################################################################
###    10 Parameter "CoCAViaR_AS_signs" model where X_t and Y_t can enter in their signed components
################################################################################

model_CoCAViaR_AS_signs <- function(theta, df, prob_level, risk_measure, forecast=FALSE, model_type="joint", m1=NULL){
  q <- length(theta)
  alpha <- prob_level$alpha
  beta <- prob_level$beta

  TT <- dim(df)[1]
  if (forecast){TT <- TT + 1} # To produce a model forecast, simply add 1 to TT

  if (model_type=="joint"){
    m1 <- m2 <- rep(NA,TT)
    m1[1] <- quantile(df$x, beta)
    # Initialize m2 depending on the risk_measure we estimate
    if (risk_measure=="CoVaR"){
      m2[1] <- quantile(df$y[df$x >= m1[1]] , alpha)
    } else if (risk_measure=="MES"){
      m2[1] <- mean(df$y[df$x >= m1[1]])
    }
    for (tt in 2:TT){
      m1[tt] <- theta[1] + theta[2]*pmax(df$x[tt-1],0) + theta[3]*(-pmin(df$x[tt-1],0)) + theta[4]*m1[tt-1]
      m2[tt] <- theta[5] + theta[6]*pmax(df$x[tt-1],0) + theta[7]*(-pmin(df$x[tt-1],0)) + theta[8]*pmax(df$y[tt-1],0) + theta[9]*(-pmin(df$y[tt-1],0)) + theta[10]*m2[tt-1]
    }
    return(list(m1=m1, m2=m2))
  } else if (model_type=="first"){
    m1 <- rep(NA,TT)
    m1[1] <- quantile(df$x, beta)
    for (tt in 2:TT){
      m1[tt] <- theta[1] + theta[2]*pmax(df$x[tt-1],0) + theta[3]*(-pmin(df$x[tt-1],0)) + theta[4]*m1[tt-1]
    }
    return(list(m1=m1))
  } else if (model_type=="second"){
    stopifnot(!is.null(m1)) # A vector of m1 must be given!
    m2 <- rep(NA,TT)
    # Initialize m2 depending on the risk_measure we estimate
    if (risk_measure=="CoVaR"){
      m2[1] <- quantile(df$y[df$x >= m1[1]] , alpha)
    } else if (risk_measure=="MES"){
      m2[1] <- mean(df$y[df$x >= m1[1]])
    }
    for (tt in 2:TT){
      m2[tt] <- theta[1] + theta[2]*pmax(df$x[tt-1],0) + theta[3]*(-pmin(df$x[tt-1],0)) + theta[4]*pmax(df$y[tt-1],0) + theta[5]*(-pmin(df$y[tt-1],0)) + theta[6]*m2[tt-1]
    }
    return(list(m1=m1, m2=m2))
  }
}


nabla_CoCAViaR_AS_signs <- function(theta, df, prob_level, risk_measure){
  q <- length(theta)
  TT <- dim(df)[1]

  m <- model_CoCAViaR_AS_signs(theta, df, prob_level, risk_measure)

  # Manually implement the gradient
  nabla_m1 <- cbind(1, lag(pmax(df$x,0)), -lag(pmin(df$x,0)), lag(m$m1), array(0, dim=c(TT,6)))
  nabla_m2 <- cbind(array(0, dim=c(TT,4)), 1, lag(pmax(df$x,0)), -lag(pmin(df$x,0)), lag(pmax(df$y,0)), -lag(pmin(df$y,0)), lag(m$m2))

  list(nabla_m1=nabla_m1, nabla_m2=nabla_m2)
}


theta_CoCAViaR_AS_signs <- function(theta, df=NULL){
  length_theta <- 10
  length_theta1 <- 4
  length_theta2 <- 6
  theta1 <- theta[1:length_theta1]
  theta2 <- theta[(length_theta1+1):(length_theta1+length_theta2)]
  theta_start_default <- c(0.01, 0.03, 0.08, 0.75,
                           0.01, 0.01, 0.01, 0.03, 0.08,       0.75)
  theta_names <- list(VaR=c("(Intercept)", "lag X+", "lag X-", "lag VaR"),
                      CoVaR=c("(Intercept)", "lag X+", "lag X-", "lag Y+", "lag Y-", "lag CoVaR"),
                      MES=c("(Intercept)", "lag X+", "lag X-", "lag Y+", "lag Y-", "lag MES"))
  return(list(theta1=theta1, theta2=theta2, length_theta=length_theta, length_theta1=length_theta1, length_theta2=length_theta2, theta_start_default=theta_start_default, theta_names=theta_names))
}




################################################################################
###    10 Parameter "CoCAViaR_AS_mixed" model where X_t and Y_t enter in their signed components, and cross terms in their absolute value
################################################################################

model_CoCAViaR_AS_mixed <- function(theta, df, prob_level, risk_measure, forecast=FALSE, model_type="joint", m1=NULL){
  q <- length(theta)
  alpha <- prob_level$alpha
  beta <- prob_level$beta

  TT <- dim(df)[1]
  if (forecast){TT <- TT + 1} # To produce a model forecast, simply add 1 to TT

  if (model_type=="joint"){
    m1 <- m2 <- rep(NA,TT)
    m1[1] <- quantile(df$x, beta)
    # Initialize m2 depending on the risk_measure we estimate
    if (risk_measure=="CoVaR"){
      m2[1] <- quantile(df$y[df$x >= m1[1]] , alpha)
    } else if (risk_measure=="MES"){
      m2[1] <- mean(df$y[df$x >= m1[1]])
    }
    for (tt in 2:TT){
      m1[tt] <- theta[1] + theta[2]*pmax(df$x[tt-1],0) + theta[3]*(-pmin(df$x[tt-1],0)) + theta[4]*abs(df$y[tt-1]) + theta[5]*m1[tt-1]
      m2[tt] <- theta[6] + theta[7]*abs(df$x[tt-1]) + theta[8]*pmax(df$y[tt-1],0) + theta[9]*(-pmin(df$y[tt-1],0)) + theta[10]*m2[tt-1]
    }
    return(list(m1=m1, m2=m2))
  } else if (model_type=="first"){
    m1 <- rep(NA,TT)
    m1[1] <- quantile(df$x, beta)
    for (tt in 2:TT){
      m1[tt] <- theta[1] + theta[2]*pmax(df$x[tt-1],0) + theta[3]*(-pmin(df$x[tt-1],0)) + theta[4]*abs(df$y[tt-1]) + theta[5]*m1[tt-1]
    }
    return(list(m1=m1))
  } else if (model_type=="second"){
    stopifnot(!is.null(m1)) # A vector of m1 must be given!
    m2 <- rep(NA,TT)
    # Initialize m2 depending on the risk_measure we estimate
    if (risk_measure=="CoVaR"){
      m2[1] <- quantile(df$y[df$x >= m1[1]] , alpha)
    } else if (risk_measure=="MES"){
      m2[1] <- mean(df$y[df$x >= m1[1]])
    }
    for (tt in 2:TT){
      m2[tt] <- theta[1] + theta[2]*abs(df$x[tt-1]) + theta[3]*pmax(df$y[tt-1],0) + theta[4]*(-pmin(df$y[tt-1],0)) + theta[5]*m2[tt-1]
    }
    return(list(m1=m1, m2=m2))
  }
}


nabla_CoCAViaR_AS_mixed <- function(theta, df, prob_level, risk_measure){
  q <- length(theta)
  TT <- dim(df)[1]

  m <- model_CoCAViaR_AS_mixed(theta, df, prob_level, risk_measure)

  # Manually implement the gradient
  nabla_m1 <- cbind(1, lag(pmax(df$x,0)), -lag(pmin(df$x,0)), lag(abs(df$y)), lag(m$m1), array(0, dim=c(TT,5)))
  nabla_m2 <- cbind(array(0, dim=c(TT,5)), 1, lag(abs(df$x)), lag(pmax(df$y,0)), -lag(pmin(df$y,0)), lag(m$m2))

  list(nabla_m1=nabla_m1, nabla_m2=nabla_m2)
}


theta_CoCAViaR_AS_mixed <- function(theta, df=NULL){
  length_theta <- 10
  length_theta1 <- 5
  length_theta2 <- 5
  theta1 <- theta[1:length_theta1]
  theta2 <- theta[(length_theta1+1):(length_theta1+length_theta2)]
  theta_start_default <- c(0.01, 0.03, 0.08, 0.01, 0.75,
                           0.01, 0.01, 0.03, 0.08,       0.75)
  theta_names <- list(VaR=c("(Intercept)", "lag X+", "lag X-", "lag |Y|", "lag VaR"),
                      CoVaR=c("(Intercept)", "lag |X|", "lag Y+", "lag Y-","lag CoVaR"),
                      MES=c("(Intercept)", "lag |X|", "lag Y+", "lag Y-","lag MES"))
  return(list(theta1=theta1, theta2=theta2, length_theta=length_theta, length_theta1=length_theta1, length_theta2=length_theta2, theta_start_default=theta_start_default, theta_names=theta_names))
}







################################################################################
###   "CoCAViaR_Z" model corresponding to the "CoCAViaR_SAV_diag" model plus added linear covariates Z.
################################################################################

model_CoCAViaR_Z <- function(theta, df, prob_level, risk_measure, forecast=FALSE, model_type="joint", m1=NULL){
  q <- length(theta)
  alpha <- prob_level$alpha
  beta <- prob_level$beta

  # Extract covariates Z:
  # No lag here as this is captured in the [tt-1,] below!
  z <- df %>%
    tibble::as_tibble() %>%
    dplyr::select(-c("Date", "x", "y")) %>%
    as.matrix()

  # Forecasts?
  TT <- dim(df)[1]
  if (forecast){TT <- TT + 1} # To produce a model forecast, simply add 1 to TT

  if (model_type=="joint"){
    m1 <- m2 <- rep(NA,TT)
    m1[1] <- quantile(df$x, beta)
    # Initialize m2 depending on the risk_measure we estimate
    if (risk_measure=="CoVaR"){
      m2[1] <- quantile(df$y[df$x >= m1[1]] , alpha)
    } else if (risk_measure=="MES"){
      m2[1] <- mean(df$y[df$x >= m1[1]])
    }
    for (tt in 2:TT){
      m1[tt] <- theta[1] + theta[2]*abs(df$x[tt-1]) + theta[3]*m1[tt-1] + as.numeric(z[tt-1,] %*% theta[4:(q/2)])
      m2[tt] <- theta[(q/2)+1] + theta[(q/2)+2]*abs(df$y[tt-1]) + theta[(q/2)+3]*m2[tt-1] + as.numeric(z[tt-1,] %*% theta[((q/2)+4):q])
    }
    return(list(m1=m1, m2=m2))
  } else if (model_type=="first"){
    m1 <- rep(NA,TT)
    m1[1] <- quantile(df$x, beta)
    for (tt in 2:TT){
      m1[tt] <- theta[1] + theta[2]*abs(df$x[tt-1]) + theta[3]*m1[tt-1] + as.numeric(z[tt-1,] %*% theta[4:q])
    }
    return(list(m1=m1))
  } else if (model_type=="second"){
    stopifnot(!is.null(m1)) # A vector of m1 must be given!
    m2 <- rep(NA,TT)
    # Initialize m2 depending on the risk_measure we estimate
    if (risk_measure=="CoVaR"){
      m2[1] <- quantile(df$y[df$x >= m1[1]] , alpha)
    } else if (risk_measure=="MES"){
      m2[1] <- mean(df$y[df$x >= m1[1]])
    }
    for (tt in 2:TT){
      m2[tt] <- theta[1] + theta[2]*abs(df$y[tt-1]) + theta[3]*m2[tt-1] + as.numeric(z[tt-1,] %*% theta[4:q])
    }
    return(list(m1=m1, m2=m2))
  }
}


nabla_CoCAViaR_Z <- function(theta, df, prob_level, risk_measure){
  TT <- dim(df)[1]
  m <- model_CoCAViaR_SAV_diag(theta, df, prob_level, risk_measure)

  # Matrix of covariates
  z <- df %>%
    tibble::as_tibble() %>%
    dplyr::select(-c("Date", "x", "y")) %>%
    as.matrix() %>%
    lag() # to capture the time lag between x[tt], y[tt] and z[tt-1]!

  # Manually implement the model gradient
  nabla_m1 <- cbind(1, lag(abs(df$x)), lag(m$m1), z, array(0, dim=c(TT,3+dim(z)[2])))
  nabla_m2 <- cbind(array(0, dim=c(TT,3+dim(z)[2])), 1, lag(abs(df$y)), lag(m$m2), z)

  list(nabla_m1=nabla_m1, nabla_m2=nabla_m2)
}


theta_CoCAViaR_Z <- function(theta, df=NULL){
  # We get a result for q when either theta or df are specified
  if (!is.null(df)){

    df_z <- df %>%
      tibble::as_tibble() %>%
      dplyr::select(-c("Date", "x", "y"))

    # Automatically get the dimension
    q <- 6 + 2*df_z %>% dim() %>% .[2]

    # Variable names
    theta_names <- list(VaR=c("(Intercept)", "lag |X|", "lag VaR", colnames(df_z)),
                        CoVaR=c("(Intercept)", "lag |Y|", "lag CoVaR", colnames(df_z)),
                        MES=c("(Intercept)", "lag |Y|", "lag MES", colnames(df_z)))
  } else {
    q <- length(theta)
  }

  if (!(length(theta)%% 2) == 0){stop("Length of theta is not even!")}
  theta1 <- theta[1:(q/2)]
  theta2 <- theta[-(1:(q/2))]
  theta_start_default <- rep(0,q)

  theta_start_default <- c(0.01, 0.1,      0.75,      rep(0,q/2-3),
                           0.01,      0.1,      0.75, rep(0,q/2-3))


  return(list(theta1=theta1, theta2=theta2, length_theta=q, length_theta1=q/2, length_theta2=q/2, theta_start_default=theta_start_default, theta_names=theta_names))
}


