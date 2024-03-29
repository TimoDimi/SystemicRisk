% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SRMroll.R
\name{SRMroll}
\alias{SRMroll}
\title{Rolling Forecasts for Dynamic (VaR, CoVaR) or (VaR, MES) Models}
\usage{
SRMroll(
  data = NULL,
  model = "CoCAViaR_SAV_diag",
  length_IS = 1000,
  refit_freq = 100,
  risk_measure = "CoVaR",
  beta = 0.95,
  alpha = 0.95,
  theta0 = NULL,
  init_method = "omega",
  optim_replications = c(1, 3)
)
}
\arguments{
\item{data}{a tsibble object that contains the variables x,y, possibly covariates and the index columns with names Date and Date_index.}

\item{model}{Specify the model type; see \link{model_fun} for details.}

\item{length_IS}{Length of the in-sample estimation length.}

\item{refit_freq}{Frequency of model refitting.}

\item{risk_measure}{The systemic risk measure under consideration; currently only the "CoVaR" is implemented.}

\item{beta}{Probability level for the VaR.}

\item{alpha}{Probability level for the CoVaR.}

\item{theta0}{Starting values for the model parameters. If NULL, then standard values are used.}

\item{optim_replications}{A vector with two integer entries indicating how often the M-estimator of the (VaR, CoVaR) model will be restarted. The default is c(1,3).}
}
\value{
A 'SRMroll' object
}
\description{
Rolling Forecasts for Dynamic (VaR, CoVaR) or (VaR, MES) Models
}
\examples{
# Get financial data and generate the data tsibble
library(rmgarch)
data(dji30retw)

data <- dji30retw \%>\%
  dplyr::mutate(DJ30=rowMeans(.)) \%>\%
  tibble::rownames_to_column(var = "Date") \%>\%
  dplyr::mutate(Date=lubridate::as_date((Date))) \%>\%
  tsibble::as_tsibble(index="Date") \%>\%
  dplyr::select(Date, x=JPM, y=DJ30) \%>\%
  dplyr::mutate(x=-100*x, y=-100*y)

# Estimate the "CoCAViaR_SAV_diag" model on the negative percentage log-returns
obj_roll <- SRMroll(data=data,
                     model="CoCAViaR_SAV_diag",
                     length_IS=700, refit_freq=100)

# Plot the forecasts
plot(obj_roll)

}
\references{
\href{https://arxiv.org/abs/2206.14275}{Dynamic Co-Quantile Regression}
}
