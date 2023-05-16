

#' Forecast Evaluation for the Systemic Risk Measures CoVaR and MES
#'
#' @param data A tsibble containing columns for Date, x,y, VaR1, VaR2, risk_measure1 and risk_measure2, where the latter two are either CoVaR or MES forecasts, as specified in the argument 'risk_measure'.
#' @param x Observation (corresponding to the VaR series)
#' @param y Observation (corresponding to the CoVaR series)
#' @param VaR1 Baseline VaR Forecasts
#' @param VaR2 Competitor VaR Forecasts
#' @param CoVaR1 Baseline CoVaR Forecasts
#' @param CoVaR2 Competitor CoVaR Forecasts
#' @param MES1 Baseline MES Forecasts
#' @param MES2 Competitor MES Forecasts
#' @param risk_measure The Systemic Risk Measure to be evaluated; either CoVaR or MES
#' @param beta Probability level of the VaR
#' @param alpha Probability level of the CoVaR
#' @param sided "onehalf" or "two"; see Fissler and Hoga (2021, arXiv) for details
#' @param cov_method covariance estimation method; "HAC" or "iid"
#' @param sig_level Significance level for the illustrative plots
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom magrittr `%>%`
SystemicRiskFCeval <- function(data=NULL,
                               x=NULL, y=NULL,
                               VaR1=NULL, VaR2=NULL,
                               CoVaR1=NULL, CoVaR2=NULL,
                               MES1=NULL, MES2=NULL,
                               risk_measure="CoVaR", beta=0.95, alpha=0.95,
                               sided="onehalf", cov_method="HAC", sig_level=0.05){


  # ToDo: Somehow deal elegantly with CoVaR and MES input data!
  if (is.null(data)){
    if (is.null(CoVaR1)){
      risk_measure <- "MES"
      data <- tsibble::tsibble(Date=1:length(x), x=x, y=y, VaR1=VaR1, VaR2=VaR2, risk_measure1=MES1, risk_measure2=MES2, index=Date)
    } else {
      risk_measure <- "CoVaR"
      data <- tsibble::tsibble(Date=1:length(x), x=x, y=y, VaR1=VaR1, VaR2=VaR2, risk_measure1=CoVaR1, risk_measure2=CoVaR2, index=Date)
    }
  } else{
    if (!is_tsibble(data)) stop("Error: Please enter a 'tsibble' object for the argument 'data'.")
  }

  data <- data %>% stats::na.omit()

  LossDiffVaR <- with(data, loss_VaR(x=x, VaR=VaR1, beta=beta)) - with(data, loss_VaR(x=x, VaR=VaR2, beta=beta))
  LossDiffrisk_measure <- switch(risk_measure,
                        MES = {with(data, loss_MES(x=x, y=y, VaR=VaR1, MES=risk_measure1)) - with(data, loss_MES(x=x, y=y, VaR=VaR2, MES=risk_measure2))},
                        CoVaR = {with(data, loss_CoVaR(x=x, y=y, VaR=VaR1, CoVaR=risk_measure1, alpha=alpha)) - with(data, loss_CoVaR(x=x, y=y, VaR=VaR2, CoVaR=risk_measure2, alpha=alpha))}
  )
  LossDiff <- cbind(LossDiffVaR, LossDiffrisk_measure)
  MeanLossDiff <- colMeans(LossDiff, na.rm = T)

  # Wald type tests
  n <- dim(data)[1]

  # HAC or iid covariance estimator
  if (cov_method=="HAC"){
    Omega <- n * sandwich::vcovHAC(lm(LossDiff~1))
  } else {
    Omega <- cov(LossDiff)
  }
  Omega_inv = tryCatch(solve(Omega), error = function(e) {MASS::ginv(Omega)})

  if (sided=="onehalf"){
    b          <- Omega_inv[1, 2] * n
    c          <- Omega_inv[2, 2] * n
    d_1n       <- mean(LossDiff[,1])
    d_2n       <- mean(LossDiff[,2])
    d_2n_tilde <- max(d_2n, -(b/c) * d_1n)

    TestStat <- n * (t(c(d_1n, d_2n_tilde))) %*% Omega_inv %*% t((t(c(d_1n, d_2n_tilde))))     # Wald test statistic
    pval <- 0.5 + 0.5 * pchisq(q=TestStat, df=2, lower.tail = FALSE) - 0.5*pchisq(q=TestStat, df=1)    # p-value of onehalf sided test
  } else {
    TestStat <- n * (t(colMeans(LossDiff))) %*% Omega_inv %*% t((t(colMeans(LossDiff))))  # Wald test statistic
    pval <- pchisq(TestStat, df=2, lower.tail=FALSE)                    # p-value
  }


  # Contour data for the ellipse plot
  # This only applies to the ONEHALF-sided test
  # Calculate and plot something different for the TWO-sided test
  new_sig <- function(x, sig_level){
    ret <- sig_level - 0.5 * (1 + x - pchisq(qchisq(1-x, df=2) , df=1))
  }
  nu_tilde = stats::uniroot(new_sig, interval = c(0,1), sig_level=sig_level)$root   # calculate "new" significance level

  npoints <- 1000
  # Get points for a central ellipse of acceptance region
  contour_data <- mixtools::ellipse(mu=c(0,0),
                                    sigma = Omega/dim(data)[1],
                                    alpha = nu_tilde,
                                    npoints=npoints,
                                    draw=FALSE) %>%
    as.data.frame()
  colnames(contour_data) <- c("x","y")
  ellipse_data <- bind_rows(tail(contour_data, -which.min(contour_data$x)),
                            head(contour_data, which.min(contour_data$x))) %>%
    dplyr::mutate(is_upper = (row_number() > dim(.)[1]/2))


  # Calculate the zone!
  if (pval > sig_level) {zone <- "yellow"}
  else if (MeanLossDiff[1] < min(ellipse_data$x)) {zone <- "red"}
  else if (MeanLossDiff[1] > max(ellipse_data$x)) {zone <- "grey"}
  else if (MeanLossDiff[2] >= 0) {zone <- "green"}
  else {zone <- "orange"}

  # Return the object
  obj <- list(MeanLossDiff=MeanLossDiff,
              LossDiffs=LossDiff,
              TestStat=TestStat,
              pval=pval,
              zone=zone,
              data=data,
              ellipse_data=ellipse_data,
              sig_level=sig_level,
              Omega=Omega,
              sided=sided,
              risk_measure=risk_measure,
              beta=beta,
              alpha=alpha
              )
  class(obj) <- "SystemicRiskFCeval"
  return(obj)
}




#'  Plot method for SystemicRiskFCeval
#'
#' @param SystemicRiskFCeval object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot.SystemicRiskFCeval <- function(obj, ...){
  p <- autoplot(obj, ...)
  print(p)
}




#' Autoplot method for SystemicRiskFCeval
#'
#' @param obj SystemicRiskFCeval object
#'
#' @return
#' @export
#'
#' @examples
autoplot.SystemicRiskFCeval <- function(obj){
  x_lim <- 1.25 * max(abs(obj$ellipse_data$x), abs(obj$MeanLossDiff[1]))
  y_lim <- 1.25 * max(abs(obj$ellipse_data$y), abs(obj$MeanLossDiff[2]))

  p <- ggplot2::ggplot(obj$ellipse_data, aes(x,y)) +
    ggplot2::geom_rect(ggplot2::aes(xmin=-Inf, xmax=min(x), ymin=-Inf, ymax=Inf), fill="red") +
    ggplot2::geom_rect(ggplot2::aes(xmin=max(x), xmax=Inf, ymin=-Inf, ymax=Inf), fill="grey") +
    ggplot2::geom_ribbon(data=obj$ellipse_data,
                         ggplot2::aes(x=x, ymin=min(y), ymax=max(y)), fill="yellow") +  # This can be done better!
    ggplot2::geom_ribbon(data=obj$ellipse_data %>% dplyr::filter(is_upper==FALSE),
                         ggplot2::aes(x=x, ymin=-Inf, ymax=y), fill="orange") +
    ggplot2::geom_ribbon(data=obj$ellipse_data  %>% dplyr::filter(is_upper==TRUE),
                         ggplot2::aes(x=x, ymin=y, ymax=Inf), fill="green") +
    ggplot2::geom_path(color="grey20") +
    ggplot2::geom_vline(xintercept=0, color="black", linetype="dotted") +
    ggplot2::geom_hline(yintercept=0, color="black", linetype="dotted") +
    ggplot2::geom_vline(aes(xintercept=min(x)), color="grey20") +
    ggplot2::geom_vline(aes(xintercept=max(x)), color="grey20") +
    ggplot2::geom_point(aes(x=obj$MeanLossDiff[1], y=obj$MeanLossDiff[2]), shape=4, color="black", size=3) +
    ggplot2::theme_bw() +
    # theme(panel.ontop = TRUE, panel.background = element_rect(color = NA, fill = NA)) +
    ggplot2::scale_x_continuous(limits = c(-x_lim,x_lim)) +
    ggplot2::scale_y_continuous(limits = c(-y_lim,y_lim)) +
    ggplot2::xlab("Average VaR Loss Difference") +
    ggplot2::ylab("Average CoVaR Loss Difference")

  p
}





#' Print method for SystemicRiskFCeval
#'
#' @param obj SystemicRiskFCeval object
#' @param digits Printed digits after the comma
#'
#' @return
#' @export
#'
#' @examples
print.SystemicRiskFCeval <- function(obj, digits=4){
  cat("Mean loss difference:\n")
  print(format(obj$MeanLossDiff, digits = digits), quote = FALSE)
  cat("Loss difference covariance:\n")
  print(format(obj$Omega, digits = digits), quote = FALSE)
}





#' VaR part of the joint (VaR, CoVaR) loss function
#'
#' @param x Observation
#' @param VaR VaR Forecast
#' @param beta Quantile level beta
#' @param b degree of homogeneity (b>=0)
#'
#' @return
#' @export
#'
#' @examples
loss_VaR <- function(x, VaR, beta, b=1){
  if(b==0){   # requires positive VaR forecasts (here h(x)=log(x))
    S <- ( 1*(x <= VaR) - beta ) * log( pmax.int(VaR, 1e-10) ) +  1*(x > VaR) * log( pmax.int(x, 1e-10) )
  }
  else{       # h(x) = sign(x) * x^b
    S <- ( 1*(x <= VaR) - beta ) * ( sign(VaR) * abs(VaR)^b - sign(x) * abs(x)^b )
  }
  return(S)
}




#' CoVaR part of the joint (VaR, CoVaR) loss function
#'
#' @param x Observation (corresponding to the VaR series)
#' @param y Observation (corresponding to the CoVaR series)
#' @param VaR VaR forecasts
#' @param CoVaR CoVaR forecasts
#' @param alpha CoVaR quantile level alpha
#' @param b degree of homogeneity (b>=0)
#'
#' @return
#' @export
#'
#' @examples
loss_CoVaR <- function(x, y, VaR, CoVaR, alpha, b=1){
  if(b==0){  # requires positive CoVaR forecasts (here g(x)=log(x))
    S <- 1*(x>VaR) * ( ( 1*(y <= CoVaR) - alpha ) * log( pmax.int(CoVaR, 1e-10) ) +  1*(y > CoVaR) * log( pmax.int(y, 1e-10) ) )
  }
  else{      # g(x) = sign(x) * x^b
    S <- 1*(x>VaR) * ( ( 1*(y <= CoVaR) - alpha ) * ( sign(CoVaR) * abs(CoVaR)^b - sign(y) * abs(y)^b) )
  }
  return(S)
}



#'  MES part of the joint (VaR, MES) loss function
#'
#' @param x Observation (corresponding to the VaR series)
#' @param y Observation (corresponding to the MES series)
#' @param VaR VaR forecasts
#' @param MES MES forecasts
#'
#' @return
#' @export
#'
#' @examples
loss_MES <- function(x, y, VaR, MES){
  S <- 1*(x>VaR) * (MES - y)^2
  return(S)
}
