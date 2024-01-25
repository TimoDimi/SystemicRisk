
#' Rolling Forecasts for Dynamic (VaR, CoVaR) or (VaR, MES) Models
#'
#' @param data a tsibble object that contains the variables x,y, possibly covariates and the index columns with names Date and Date_index.
#' @param model Specify the model type; see \link{model_fun} for details.
#' @param length_IS Length of the in-sample estimation length.
#' @param refit_freq Frequency of model refitting.
#' @param risk_measure The systemic risk measure under consideration; currently only the "CoVaR" is implemented.
#' @param beta Probability level for the VaR.
#' @param alpha Probability level for the CoVaR.
#' @param theta0 Starting values for the model parameters. If NULL, then standard values are used.
#' @param optim_replications A vector with two integer entries indicating how often the M-estimator of the (VaR, CoVaR) model will be restarted. The default is c(1,3).
#'
#' @return A 'SRMroll' object
#'
#' @export

#'
#' @examples
#' # Get financial data and generate the data tsibble
#' library(rmgarch)
#' data(dji30retw)
#'
#' data <- dji30retw %>%
#'   dplyr::mutate(DJ30=rowMeans(.)) %>%
#'   tibble::rownames_to_column(var = "Date") %>%
#'   dplyr::mutate(Date=lubridate::as_date((Date))) %>%
#'   tsibble::as_tsibble(index="Date") %>%
#'   dplyr::select(Date, x=JPM, y=DJ30) %>%
#'   dplyr::mutate(x=-100*x, y=-100*y)
#'
#' # Estimate the "CoCAViaR_SAV_diag" model on the negative percentage log-returns
#' obj_roll <- SRMroll(data=data,
#'                      model="CoCAViaR_SAV_diag",
#'                      length_IS=700, refit_freq=100)
#'
#' # Plot the forecasts
#' plot(obj_roll)
#'
#' @references \href{https://arxiv.org/abs/2206.14275}{Dynamic Co-Quantile Regression}
#' @importFrom magrittr `%>%`
SRMroll <- function(data=NULL,
                     model="CoCAViaR_SAV_diag",
                     length_IS=1000, refit_freq=100,
                     risk_measure="CoVaR", beta=0.95, alpha=0.95,
                     theta0=NULL, init_method="omega", optim_replications=c(1,3)){

  # data must be a tsibble object with columns Date, x, y and possibly covariates
  TT <- dim(data)[1]
  length_OOS <- TT - length_IS
  refit_points <- seq(length_IS+1,TT,by=refit_freq)

  FC_df <- data %>%
    dplyr::slice((length_IS+1):TT) %>%
    dplyr::mutate(m1_FC=NA, m2_FC=NA)

  SRM_objects <- list()

  # Loop over all OOS days
  for (tt in (length_IS+1):TT){
    data_tt <- data[(tt-length_IS):(tt-1),]

    # Only refit at certain points!
    if (tt %in% refit_points){

      # Iterative starting values from the previous fit
      if (tt==refit_points[1]) theta_start <- theta0 else theta_start <- SRM_obj$theta

      # Possible extension: Include some error handling if this fit fails?
      SRM_obj <- SRM(data=data_tt, model=model,
                       risk_measure=risk_measure, beta=beta, alpha=alpha,
                       theta0=theta_start, init_method=init_method, optim_replications=optim_replications)
      SRM_objects <- append(SRM_objects, list(SRM_obj))
    }
    FCs <- forecast(SRM_obj,
                    newdata=data_tt)
                      # dplyr::select(x,y))

    FC_df[tt-length_IS,] <- FC_df %>%
      dplyr::slice(tt-length_IS) %>%
      dplyr::mutate(m1_FC=FCs[1], m2_FC=FCs[2])
  }

  FC_df <- FC_df %>%
    dplyr::rename(VaR=m1_FC, risk_measure=m2_FC) %>%
    dplyr::select(Date, x, y, VaR, risk_measure)


  # Return an object of class "SRMroll"
  obj <- list(
    FC_df=FC_df,
    data=data,
    risk_measure=risk_measure,
    alpha=alpha,
    beta=beta,
    SRM_objects=SRM_objects)

  class(obj) <- "SRMroll"

  obj
}


#' SRMroll print method
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
print.SRMroll <- function(obj){
  print(obj$FC_df)
}





#' SRMroll plot method
#'
#' @param obj
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot.SRMroll <- function(obj, ...){
  p <- autoplot(obj, ...)
  print(p)
}





#' SRMroll autoplot method
#'
#' @param obj
#' @param facet_names
#'
#' @return
#' @export
#'
#' @examples
#' @import ggplot2
#' @importFrom magrittr `%>%`
autoplot.SRMroll <- function(obj, facet_names=NULL){
  # Add options etc! this is a little lame so far...
  # Shall we include the estimates into the df?
  # ToDo: Improve facet and legend labels

  if (is.null(facet_names)){
    if (obj$risk_measure=="CoVaR"){facet_names <- c("X / VaR","Y / CoVaR")}
    else{facet_names <- c("X / VaR","Y / MES")}
  }

  df_tmp <- obj$FC_df %>%
    dplyr::mutate(VaR_violation=(x > VaR))

  df_long <- dplyr::left_join(
    df_tmp %>%
      dplyr::select(Date, VaR_violation, x, y) %>%
      reshape2::melt(id.vars=c("Date","VaR_violation"), variable.name="Symbol", value.name="NegativeReturns"),
    df_tmp %>%
      dplyr::select(Date, VaR_violation, VaR, risk_measure) %>%
      dplyr::rename(x=VaR, y=risk_measure) %>%
      reshape2::melt(id.vars=c("Date","VaR_violation"), variable.name="Symbol", value.name="risk_measureForecasts"),
    by=c("Date", "VaR_violation", "Symbol")) %>%
    tibble::as_tibble()

  levels(df_long$Symbol) <- facet_names

  p <- ggplot2::ggplot(df_long %>% arrange(VaR_violation)) +
    ggplot2::geom_point(aes(x=Date, y=NegativeReturns, color=VaR_violation)) +
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


