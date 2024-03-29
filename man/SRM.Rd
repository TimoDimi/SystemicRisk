% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SRM.R
\name{SRM}
\alias{SRM}
\alias{SRM.default}
\alias{SRM.fit}
\title{Joint Dynamic Models for the (VaR, CoVaR) or (VaR, MES)}
\usage{
SRM(...)

\method{SRM}{default}(
  data = NULL,
  model = "CoCAViaR_SAV_diag",
  risk_measure = "CoVaR",
  beta = 0.95,
  alpha = 0.95,
  theta0 = NULL,
  init_method = "omega",
  optim_replications = c(1, 3)
)

\method{SRM}{fit}(
  data,
  model,
  risk_measure,
  beta,
  alpha,
  theta0,
  init_method,
  optim_replications
)
}
\arguments{
\item{...}{Further arguments (does not apply here).}

\item{data}{A tsibble that contains the variables x, y, possibly covariates and an index columns with the name Date.}

\item{model}{Specify the model type; see \link{model_fun} for details.}

\item{risk_measure}{The systemic risk measure under consideration; currently only the "CoVaR" is implemented.}

\item{beta}{Probability level for the VaR.}

\item{alpha}{Probability level for the CoVaR.}

\item{theta0}{Starting values for the model parameters. If NULL, then standard values are used.}

\item{optim_replications}{A vector with two integer entries indicating how often the M-estimator of the (VaR, CoVaR) model will be restarted. The default is c(1,3).}
}
\value{
A 'SRM' object
}
\description{

}
\details{
ToDo: Detailed description
}
\examples{
### (1) A predictive SRM
library(tsibble)

# (1) Simulate bivariate data (x,y) with covariate z
eps <- mvtnorm::rmvt(n=1000, sigma=matrix(c(1,0.5,0.5,1),2), df = 8)
z <- rnorm(1000)
xy <- c(1,1) + cbind(2*z, 2.5*z) + eps

# Collect data as tsibble
data <- tsibble(Date=1:length(z),
                x=xy[,1],
                y=xy[,2],
                Intercept=1,
                z=z,
                index=Date)

# Estimate the 'joint_linear' SRM regression model
obj_fit <- SRM(data=data,
                model = "joint_linear",
                beta=0.95,
                alpha=0.95)

# Estimate the standard errors of the parameters
summary(obj_fit)



### (1) A Dynamic SRM Forecasting Model

 # Get financial data and generate the data tsibble
library(rmgarch)
library(dplyr)
library(lubridate)
data(dji30retw)

data <- dji30retw \%>\%
  dplyr::mutate(DJ30=rowMeans(.)) \%>\%
  tibble::rownames_to_column(var = "Date") \%>\%
  dplyr::mutate(Date=lubridate::as_date((Date))) \%>\%
  tsibble::as_tsibble(index="Date") \%>\%
  dplyr::select(Date, x=JPM, y=DJ30) \%>\%
  dplyr::mutate(x=-100*x, y=-100*y)
# Estimate the "CoCAViaR_SAV_diag" model on the negative percentage log-returns
obj <- SRM(data=data,
           model="CoCAViaR_SAV_diag")

# Covariance estimation and display the parameter estimates with standard errors
summary(obj)

# Plot the estimated time series
plot(obj)


}
\references{
\href{https://arxiv.org/abs/2206.14275}{Dynamic CoVaR Modeling} and \href{https://arxiv.org/abs/2311.13327}{Regressions Under Adverse Conditions}
}
\seealso{
\code{\link{SRMroll}} for a class of rolling window SRM models,
\code{\link{model_functions}} for the specification of model functions
}
