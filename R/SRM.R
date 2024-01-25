

#' Joint Dynamic Models for the (VaR, CoVaR) or (VaR, MES)
#'
#' Estimates a joint dynamic semiparametric 'SRM' model for the pair (VaR, CoVaR):
#' \deqn{VaR_\beta(X_t | Z_t) = v_t(\theta^v)
#' \deqn{CoVaR_\alpha(Y_t | Z_t) = c_t(\theta^c)
#'
#' ToDo: Detailed description
#'
#' @param data A tsibble that contains the variables x, y, possibly covariates and an index columns with the name Date.
#' @param model Specify the model type; see \link{model_fun} for details.
#' @param risk_measure The systemic risk measure under consideration; currently only the "CoVaR" is implemented.
#' @param beta Probability level for the VaR.
#' @param alpha Probability level for the CoVaR.
#' @param theta0 Starting values for the model parameters. If NULL, then standard values are used.
#' @param optim_replications A vector with two integer entries indicating how often the M-estimator of the (VaR, CoVaR) model will be restarted. The default is c(1,3).
#' @param ... Further arguments (does not apply here).
#'
#' @return A 'SRM' object
#'
#' @seealso \code{\link{SRMroll}} for a class of rolling window SRM models,
#' \code{\link{model_functions}} for the specification of model functions
#'
#' @examples
#' ### (1) A predictive SRM
#' library(tsibble)
#'
#' # (1) Simulate bivariate data (x,y) with covariate z
#' eps <- mvtnorm::rmvt(n=1000, sigma=matrix(c(1,0.5,0.5,1),2), df = 8)
#' z <- rnorm(1000)
#' xy <- c(1,1) + cbind(2*z, 2.5*z) + eps
#'
#' # Collect data as tsibble
#' data <- tsibble(Date=1:length(z),
#'                 x=xy[,1],
#'                 y=xy[,2],
#'                 Intercept=1,
#'                 z=z,
#'                 index=Date)
#'
#' # Estimate the 'joint_linear' SRM regression model
#' obj_fit <- SRM(data=data,
#'                 model = "joint_linear",
#'                 beta=0.95,
#'                 alpha=0.95)
#'
#' # Estimate the standard errors of the parameters
#' summary(obj_fit)
#'
#'
#'
#' ### (1) A Dynamic SRM Forecasting Model
#'
#'  # Get financial data and generate the data tsibble
#' library(rmgarch)
#' library(dplyr)
#' library(lubridate)
#' data(dji30retw)
#'
#' data <- dji30retw %>%
#'   dplyr::mutate(DJ30=rowMeans(.)) %>%
#'   tibble::rownames_to_column(var = "Date") %>%
#'   dplyr::mutate(Date=lubridate::as_date((Date))) %>%
#'   tsibble::as_tsibble(index="Date") %>%
#'   dplyr::select(Date, x=JPM, y=DJ30) %>%
#'   dplyr::mutate(x=-100*x, y=-100*y)

#' # Estimate the "CoCAViaR_SAV_diag" model on the negative percentage log-returns
#' obj <- SRM(data=data,
#'            model="CoCAViaR_SAV_diag")
#'
#' # Covariance estimation and display the parameter estimates with standard errors
#' summary(obj)
#'
#' # Plot the estimated time series
#' plot(obj)
#'
#'
#' @references \href{https://arxiv.org/abs/2206.14275}{Dynamic CoVaR Modeling} and \href{https://arxiv.org/abs/2311.13327}{Regressions Under Adverse Conditions}
#' @rdname SRM
#' @export
SRM <- function(...) {
  UseMethod('SRM')
}



#' @rdname SRM
#' @export
SRM.default <- function(data=NULL,
                         model="CoCAViaR_SAV_diag", risk_measure="CoVaR", beta=0.95, alpha=0.95,
                         theta0=NULL, init_method="omega", optim_replications=c(1,3)){

  fit <- SRM.fit(data=data,
                  model=model, risk_measure=risk_measure, beta=beta, alpha=alpha,
                  theta0=theta0, init_method=init_method, optim_replications=optim_replications)

  fit$call <- match.call()
  fit
}




#' @rdname SRM
#' @export
SRM.fit <- function(data,
                    model, risk_measure, beta, alpha,
                    theta0, init_method, optim_replications){

  # Collect input data as a tibble
  # data <- collect_data(data=data, x=x, y=y, z=z)

  # Check inputs:
  if (!is_tsibble(data)) stop("Error: Please enter a 'tsibble' object for the argument 'data'.")


  models_implemented <- c("joint_linear",
                          "CoCAViaR_SAV_diag", "CoCAViaR_SAV_fullA", "CoCAViaR_SAV_full",
                          "CoCAViaR_AS_pos", "CoCAViaR_AS_signs", "CoCAViaR_AS_mixed",
                          "CoCAViaR_8pCrossCoVaR", "CoCAViaR_10p", "CoCAViaR_7pVaRViolation",
                          "CoCAViaR_Z")

  if ( !( is.list(model) | (is.character(model) & (model %in% models_implemented)) ) ){
    stop("Please provide for 'model' either a list of functions or a string matching one of the pre-implemented models!")
  }

  TT <- dim(data)[1]
  risk_measure <- match.arg(risk_measure, c("MES","CoVaR"))
  prob_level <- switch(risk_measure,
                       MES = {list(beta=beta, alpha=NA)},
                       CoVaR = {list(beta=beta, alpha=alpha)})

  # Assign default implemented starting values if theta0==NULL
  if (is.null(theta0)){
    theta0_list <- theta_fun(model=model, theta=theta0, df=data)
    theta0 <- theta0_list$theta_start_default
  }
  # Split the starting value!
  theta0_list <- theta_fun(model=model, theta=theta0, df=data)
  theta01 <- theta0_list$theta1
  theta02 <- theta0_list$theta2


  ### VaR Optimization
  thetav_est_rep <- matrix(NA, nrow=optim_replications[1], ncol=length(theta01))
  Mest_obj1 <- rep(NA, optim_replications[1])

  for (rep in 1:optim_replications[1]){
    opt_v <- optim(par=theta01,
                   fn=function(theta,...){
                     mean(loss_model_VaR(theta,...),na.rm=TRUE)},
                   df=data, model=model, prob_level=prob_level, init_method=init_method)

    thetav_est_rep[rep, ] <- opt_v$par
    Mest_obj1[rep] <- opt_v$value

    theta01 <- thetav_est_rep[rep, ] +
      MASS::mvrnorm(n=1, mu=rep(0,length(theta01)), Sigma=0.1*diag(length(theta01)))
    while(!is.finite(mean(loss_model_VaR(theta01,df=data, model=model, prob_level=prob_level, init_method=init_method), na.rm=TRUE))){
      theta01 <- thetav_est_rep[rep, ] +
        MASS::mvrnorm(n=1, mu=rep(0,length(theta01)), Sigma=0.1*diag(length(theta01)))
    }
  }
  thetav_est <- thetav_est_rep[which.min(Mest_obj1),]
  m1_est <- model_fun(thetav_est, df=data, prob_level=prob_level, model=model, model_type="first", init_method=init_method)$m1


  ### Second step: CoVaR/MES optimization
  theta2_est_rep <- matrix(NA, nrow=optim_replications[2], ncol=length(theta02))
  Mest_obj2 <- rep(NA, optim_replications[2])

  for (rep in 1:optim_replications[2]){
    opt_m2 <- optim(par=theta02,
                    fn=function(theta,...){
                      switch(risk_measure,
                             MES = {mean(loss_model_MES(theta,...),na.rm=TRUE)},
                             CoVaR = {mean(loss_model_CoVaR(theta,...),na.rm=TRUE)})},
                    df=data, m1=m1_est, model=model, prob_level=prob_level, init_method=init_method)

    theta2_est_rep[rep, ] <- opt_m2$par
    Mest_obj2[rep] <- opt_m2$value


    theta02 <- theta2_est_rep[rep, ] +
      MASS::mvrnorm(n=1, mu=rep(0,length(theta02)), Sigma=0.1*diag(length(theta02)))
    while(!is.finite(switch(risk_measure,
                            MES = {mean(loss_model_MES(theta02, df=data, m1=m1_est, model=model, prob_level=prob_level, init_method=init_method), na.rm=TRUE)},
                            CoVaR = {mean(loss_model_CoVaR(theta02, df=data, m1=m1_est, model=model, prob_level=prob_level, init_method=init_method), na.rm=TRUE)}))){
      theta02 <- theta2_est_rep[rep, ] +
        MASS::mvrnorm(n=1, mu=rep(0,length(theta02)), Sigma=0.1*diag(length(theta02)))
    }
  }
  theta2_est <- theta2_est_rep[which.min(Mest_obj2),]

  theta_est <- c(thetav_est, theta2_est)

  # Create a data.frame with in sample predictions
  m_est  <- model_fun(theta_est, data, prob_level, model, risk_measure, init_method=init_method)
  m_est$m1[1] <- m_est$m2[1] <- NA # The first values are NAs, just for the estimation of dynamic models we have to select starting values!
  data <- data %>%
    dplyr::mutate(VaR=as.numeric(m_est$m1), risk_measure:=as.numeric(m_est$m2)) %>%
    dplyr::select(Date, x, y, tidyselect::everything(), VaR, risk_measure)

  # create a variable with colnames
  if (model == "joint_linear"){
    colnames_help <- data %>%
      as_tibble() %>%
      dplyr::select(-c(Date, x, y, VaR, risk_measure)) %>%
      colnames()
    colnames <- list(VaR=colnames_help, CoVaR=colnames_help, MES=colnames_help)
  } else {
    colnames <- theta0_list$theta_names
  }

  # Return an object of class "SRM"
  fit <- list(
    theta=theta_est,
    data=data,
    colnames=colnames,
    risk_measure=risk_measure,
    model=model,
    prob_level=prob_level,
    call=match.call()
  )

  class(fit) <- "SRM"
  fit
}


#' SRM summary method
#'
#' @param SRM_object A SRM object
#' @param method The method to compute the standard errors. Either "asymptotic" or "boot"
#' @param ... Further input parameters
#'
#' @return A summary.SRM object
#' @export
summary.SRM <- function(SRM_object, method='asymptotic',...){
  cov <- cov.SRM(SRM_object, method,...)
  se <- sqrt(diag(cov))
  tval <- SRM_object$theta/se
  coef_mat <- cbind(
    Estimate     = SRM_object$theta,
    `Std. Error` = se,
    `t value`    = tval,
    `Pr(>|t|)`   = 2 * stats::pt(abs(tval), dim(SRM_object$data)[1] - length(SRM_object$theta) , lower.tail = FALSE)
  )

  # VaR and CoVaR coef mats
  coef_mat_VaR <- coef_mat[1:length(SRM_object$colnames$VaR),]
  rownames(coef_mat_VaR) <- SRM_object$colnames$VaR

  coef_mat_risk_measure <- coef_mat[-(1:length(SRM_object$colnames[[SRM_object$risk_measure]])),]
  rownames(coef_mat_risk_measure) <- SRM_object$colnames[[SRM_object$risk_measure]]

  # Return an object of class "summary.SRM"
  object <- list(cov=cov,
                 se=se,
                 coef_mat_VaR=coef_mat_VaR,
                 coef_mat_risk_measure=coef_mat_risk_measure,
                 risk_measure=SRM_object$risk_measure)
  class(object) <- "summary.SRM"
  object
}


#' SRM p arameter covariance estimation
#'
#' @param SRM_object A SRM object
#' @param method The method to compute the standard errors. Either "asymptotic" or "boot"
#' @param ... Further input parameters
#'
#' @return covariance matrix
#' @export
cov.SRM <- function(SRM_object, method, ...) {
  if (method == 'asymptotic') {
    cov <- vcovA(SRM_object, ...)
  } else if(method == 'boot') {
    cov <- vcovB(SRM_object, ...)
  } else if(method == 'boot_fast') {
    cov <- vcovBfast(SRM_object, ...)
  } else {
    stop('method can be asymptotic, boot, or boot_fast')
  }
  cov
}



#' print method for the sum.SRM object
#'
#' @param sum.SRM_object A sum.SRM object
#'
#' @return Nothing
#' @export
print.summary.SRM <- function(sum.SRM_object){
  # Print the VaR and CoVaR coefficients in separate calls:
  cat("\nVaR Coefficients:\n")
  stats::printCoefmat(sum.SRM_object$coef_mat_VaR,
                      signif.legend=FALSE)
  # cat("\nCoVaR Coefficients:\n")
  cat("\n",sum.SRM_object$risk_measure,"Parameter Estimates:\n")
  stats::printCoefmat(sum.SRM_object$coef_mat_risk_measure)
}



#' print method for the SRM class
#'
#' @param obj SRM object
#' @param digits printed digits after the comma
#'
#' @export
print.SRM <- function(obj, digits=4){
  theta_info <- theta_fun(model=obj$model, theta=obj$theta, df=obj$data %>% dplyr::select(-c("VaR", "risk_measure")))
  q1 <- theta_info$length_theta1
  q2 <- theta_info$length_theta2

  cat("Call:\n")
  cat(deparse(obj$call), "\n")
  cat("\nVaR Parameter Estimates:\n")
  print(format(obj$theta[1:q1], digits = digits), quote = FALSE)
  cat("\n",obj$risk_measure,"Parameter Estimates:\n")
  print(format(obj$theta[(q1+1):(q1+q2)], digits = digits), quote = FALSE)
}


#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot



#' Autoplot a SRM object
#'
#' @param obj SRM object
#' @param facet_names Titles of the facets
#'
#' @return A ggplot object
#' @export
#'
#' @import ggplot2
#' @importFrom magrittr `%>%`
autoplot.SRM <- function(obj, facet_names=NULL){

  if (is.null(facet_names)){
    if (obj$risk_measure=="CoVaR"){facet_names <- c("X / VaR","Y / CoVaR")}
    else{facet_names <- c("X / VaR","Y / MES")}
  }

  df_tmp <- obj$data %>%
    dplyr::mutate(VaR_violation=(x > VaR))

  # ToDo: Use dplyr functions instead of reshape2!
  df_long <- dplyr::left_join(
    df_tmp %>%
      dplyr::select(Date, VaR_violation, x, y) %>%
      reshape2::melt(id.vars=c("Date","VaR_violation"), variable.name="Symbol", value.name="NegativeReturns"),
    df_tmp %>%
      dplyr::select(Date, VaR_violation, VaR, risk_measure) %>%
      dplyr::rename(x=VaR, y=risk_measure) %>%
      reshape2::melt(id.vars=c("Date","VaR_violation"), variable.name="Symbol", value.name="risk_measureForecasts"),
    by=c("Date", "VaR_violation", "Symbol")) %>%
    tibble::as_tibble() %>%
    dplyr::rename("VaR Exceedance" = "VaR_violation")

  # Rename for the facet names
  levels(df_long$Symbol) <- facet_names

  # ggplot2
  p <- ggplot(df_long %>% arrange(`VaR Exceedance`)) +
    ggplot2::geom_point(aes(x=Date, y=NegativeReturns, color=`VaR Exceedance`)) +
    ggplot2::scale_colour_manual(values = c("grey", "black")) +
    ggnewscale::new_scale_color() + # For two color scales!!!
    ggplot2::geom_line(aes(x=Date, y=risk_measureForecasts, color=Symbol)) +
    ggplot2::scale_colour_manual(values = c("red", "blue")) +
    ggplot2::facet_wrap(~Symbol, ncol=1) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="bottom") +
    ggplot2::ylab("Negative Returns")

  p
}



#' SRM plot method
#'
#' @param obj
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot.SRM <- function(obj, ...){
  p <- autoplot(obj, ...)
  print(p)
}


#' SRM fitted method
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
fitted.SRM <- function(obj){
  m <- model_fun(obj$theta, obj$data, obj$prob_level, obj$model, obj$risk_measure)
  cbind(m$m1, m$m2)
}



#' SRM residuals
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
residuals.SRM <- function(obj){
  cbind(obj$data$x, obj$data$y) - fitted(obj)
}

# Not sure which one to use. Maybe both and use predict for cross-sectional predictions, and forecast for time-series forecasting!
predict.SRM <- function(obj){
}



#' Generic forecast method
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
forecast <- function(...) {
  UseMethod('forecast')
}


#' SRM forecast method
#'
#' @param obj SRM object
#' @param newdata new data set to use for forecasting (for models with regressors or in dynamic models if the model is no re-estimated in this specific time step)
#'
#' @return Vector of forecasts
#' @export
forecast.SRM <- function(obj, newdata=NULL){
  if (is.null(newdata)) {
    df <- obj$data %>%
      dplyr::select(-c(VaR,risk_measure))
  } else {
    # ToDo: Check if newdata is reasonable!
    df <- newdata
  }

  m <- model_fun(obj$theta, df=df, prob_level=obj$prob_level, model=obj$model, risk_measure=obj$risk_measure, forecast=TRUE)
  forecast <- c(m1_FC=tail(m$m1,1), m2_FC=tail(m$m2,1))
  forecast
}



#' Generic PIforecast method
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
PIforecast <- function(...) {
  UseMethod('PIforecast')
}




#' SRM PIforecast method
#'
#' @param obj
#' @param level
#' @param newdata
#'
#' @return
#' @export
#'
#' @examples
PIforecast.SRM <- function(obj, level=0.9, newdata=NULL){
  se <- summary(obj)$se

  if (is.null(newdata)) {
    df <- obj$data %>%
      dplyr::select(-c(VaR,risk_measure))
  } else {
    # ToDo: Check if newdata is reasonable!
    df <- newdata
  }


  # Very naive implementation for now!
  m <- model_fun(obj$theta, df=df, prob_level=obj$prob_level, model=obj$model, risk_measure=obj$risk_measure, forecast=TRUE)
  # No!!!!! Only apply the recursion in the last iteration!
  theta_low <- obj$theta + qnorm((1-level)/2)*se
  theta_up <- obj$theta - qnorm((1-level)/2)*se
  TT <- dim(df)[1]

  ##################### This is incorrect! We need to take into account the whole past!
  # m1_low <- theta_low[1] + theta_low[2]*abs(df$x[TT]) + theta_low[3]*m$m1[TT]
  # m2_low <- theta_low[4] + theta_low[5]*abs(df$y[TT]) + theta_low[6]*m$m2[TT]
  #
  # m1_up <- theta_up[1] + theta_up[2]*abs(df$x[TT]) + theta_up[3]*m$m1[TT]
  # m2_up <- theta_up[4] + theta_up[5]*abs(df$y[TT]) + theta_up[6]*m$m2[TT]


  # This is correct I think: It takes into account the different theta's on the whole past!
  m_low <- model_fun(theta_low, df=df, prob_level=obj$prob_level, model=obj$model, risk_measure=obj$risk_measure, forecast=TRUE)
  m_up <- model_fun(theta_up, df=df, prob_level=obj$prob_level, model=obj$model, risk_measure=obj$risk_measure, forecast=TRUE)

  PIforecast <- list(m1_FC=c(tail(m_low$m1,1), tail(m_up$m1,1)),
                     m2_FC=c(tail(m_low$m2,1), tail(m_up$m2,1)))
  PIforecast
}
