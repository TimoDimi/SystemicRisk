

loss_model_VaR <- function(theta, df, model, prob_level){
  # The parameter theta only contains the VaR parameters here for the M-estimator!
  beta <- prob_level$beta

  m <- model_fun(theta=theta, df=df, prob_level=prob_level, model=model, model_type="first")
  v <- m$m1
  loss <- (v - df$x) * ((df$x <= v) - beta)
  return(loss)
}



loss_model_CoVaR <- function(theta, df, m1, model, prob_level){
  # The parameter theta only contains the CoVaR parameters here for the second step M-estimator!
  # One has to pass predictions m1 here!
  alpha <- prob_level$alpha

  m <- model_fun(theta=theta, df=df, prob_level=prob_level, model=model, risk_measure="CoVaR", model_type="second", m1=m1)
  v <- m1
  c <- m$m2

  loss_c <- (df$x > v) * ((c - df$y) * ((df$y <= c) - alpha))

  return(loss_c)
}



loss_model_MES <- function(theta, df, m1, model, prob_level){
  # The parameter theta only contains the MES parameters here for the second step M-estimator!
  # One has to pass predictions m1 here!

  m <- model_fun(theta=theta, df=df, prob_level=prob_level, model=model, risk_measure="MES", model_type="second", m1=m1)
  v <- m1
  m <- m$m2

  loss_m <- (df$x > v) * (m - df$y)^2

  return(loss_m)
}


collect_data <- function(data, x, y, z=NULL){
  # This function collects the data from either data, or x,y,z and saves it as a data frame

  # Checks if everything is entered correctly!
  if (is.null(data) & (is.null(x) | is.null(y) )){stop("data and (x or y) are not specified")}

  # if (!is.numeric(y) | !is.numeric(x) | !is.numeric(z))
  #   stop("There is a value in x,y, or z that is not numeric!")

  if (any(is.na(y)) | any(is.na(x)) | any(is.na(z)))
    stop("Data contains NAs!")

  if (!(all(is.finite(y)) & all(is.finite(x)) & all(is.finite(z))))
    stop("Not all values are finite!")


  # check dimensions of x,y and z
  # if ( (length(x)!=length(y)) || (length(x)!=dim(z)[1]) ) {stop("Dimensions of x,y, or z do not match!")}

  # Stop if columns of z are multicolinear
  # if (det(t(z) %*% z) < 10^(-8)) {stop("The matrix z is (close to) perfect multicolinearity!")}

  # Collect data
  if (is.null(data)){
    # Case that x,y (and z) are specified
    if (is.null(z)){
      if ( (length(x)!=length(y))) {stop("Dimensions of x and y do not match!")}
      data <- data.frame(x=x, y=y)
    } else {
      # If z is a vector, convert to a nx1 matrix
      if (!is.matrix(z)) { z <- matrix(z, ncol=1, dimnames=list(NULL, c("z1"))) }

      # Delete any constant predictors (i.e., intercepts) of z; add later on
      if (any(apply(z, 2, var)==0)) {
        intercept_index <- which(apply(z, 2, var)==0)
        z <- z[,-intercept_index]
      }

      # Empty colnames of z?
      if(is.null(colnames(z))){
        colnames(z) <- paste0("z", 1:(dim(z)[2]))
      }

      # Add intercept
      z <- cbind("Intercept"=1, z)

      data <- data.frame(x=x, y=y, z) %>% as_tibble()
    }
  } else {
  data <- data %>%
    # dplyr::select(x,y, contains("z"), contains("Date"), contains("Date_index")) %>%
    tibble::as_tibble()
  }

  # Treat empty Date_index and Date columns
  # if(!("Date_index" %in% colnames(data))) {data <- data %>% tibble::add_column(Date_index=1:nrow(data))}
  if(!("Date" %in% colnames(data))) {data <- data %>% dplyr::mutate(Date=lubridate::as_date(1:nrow(data)))}

  data %>% dplyr::select(Date, x, y, everything())
}









#' Co Value at Risk of a bivariate t distribution
#'
#' @param alpha Probability level for the CoVaR
#' @param beta Probability level for the VaR
#' @param nu Degrees of freedom of the bivariate t-distribution
#' @param H Covariance matrix of the bivariate t-distribution (2x2)
#'
#' @return
#' @export
#'
#' @examples
CoVaR_tdist <- function(alpha, beta, nu, H){
  root.CoVaR.t <- function(x){
    sigma <-  (nu - 2) / nu * H
    sd.wt <- sqrt(nu / (nu-2) )
    VaR   <- qt(p = beta, df = nu) / sd.wt * sqrt( H[1,1] )
    prob  <- mvtnorm::pmvt(lower=c(VaR, x), upper=c(Inf, Inf), df=nu, sigma=sigma)

    return( prob - (1-alpha) * (1-beta) )
  }
  uniroot(root.CoVaR.t, interval = c(0, 10), extendInt = "yes")$root
}




#' Marginal Expected Shortfall of a bivariate t distribution
#'
#' @param beta Probability level for the VaR
#' @param nu Degrees of freedom of the bivariate t-distribution
#' @param H Covariance matrix of the bivariate t-distribution (2x2)
#'
#' @return
#' @export
#'
#' @examples
MES_tdist <- function(beta, nu, H){
  Corr <- stats::cov2cor(H)
  rho <- Corr[1,2]

  if(rho==0){
    MES <- 0
  }else{
    # sd.wt <- sqrt(nu / (nu-2) )
    f.MES <- function(x, nu, H){ x[2] * mvtnorm::dmvt(x, sigma = H * (nu - 2) / nu, df = nu, log = FALSE) } # "x" is vector
    VaR   <- qt(p = beta, df = nu) * sqrt((nu-2)/nu * H[1,1])
    MES   <- 1/(1 - beta) * cubature::adaptIntegrate(f.MES, lowerLimit = c(VaR, -Inf), upperLimit = c(Inf, Inf), nu, H)$integral
  }
  return(MES)
}



# ===========================================================
# Function for MES calculation for bivariate t-distribution
# ===========================================================

# Inputs:
# nu   = degrees of freedom
# rho  = correlation
# beta = risk level

MES.true.t.rho <- function(nu, rho, beta){
  if(rho==0){
    MES <- 0
  }
  else{
    sd.wt <- sqrt(nu / (nu-2) )
    f.MES <- function(x, nu, rho) { x[2] * mvtnorm::dmvt(x, sigma = (nu - 2) / nu * matrix(c(1, rho, rho, 1), nrow=2), df = nu, log = FALSE) } # "x" is vector
    VaR   <- qt(p = beta, df = nu) / sd.wt
    MES   <- 1/(1 - beta) * cubature::adaptIntegrate(f.MES, lowerLimit = c(VaR, -Inf), upperLimit = c(Inf, Inf), nu, rho)$integral
  }
  return( MES )
}










