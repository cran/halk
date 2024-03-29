% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assess_age_models.R
\name{integral_quotient}
\alias{integral_quotient}
\title{Compute the quotient of integrals as a measure of percent error between two
curves}
\usage{
integral_quotient(
  ref_curve_params,
  comp_curve_params,
  min_x,
  max_x,
  curve_fun = function(x, linf, k, t0) {
     out <- linf * (1 - exp(-k * (x - t0)))
    
    return(out)
 }
)
}
\arguments{
\item{ref_curve_params}{A list of named parameters for the reference curve
(i.e. the standard that is being compared to)}

\item{comp_curve_params}{A list of named parameters for the curve that is
being compared}

\item{min_x}{The minimum value across which to integrate}

\item{max_x}{The maximum value across which to integrate}

\item{curve_fun}{The function that is being compared. Defaults to an
anonymous function that is the von Bertalanffy growth function.}
}
\value{
A value of the area between curves divided by the area under the
  reference curve
}
\description{
This is a method for comparing how "close" or "accurate" one curve is to
another (reference) curve. The method works by dividing the area between the
curves by the area under the reference curve. See Details for more
information
}
\details{
The integral quotient method provides a basis for comparison between two
curves by dividing the area between the curves by the area under the
reference curve (i.e. the quotient of integrals)
}
\examples{
ref_curve_params <- list(linf = 60, k = 0.25, t0 = -0.5)
comp_curve_params <- list(linf = 62, k = 0.25, t0 = -0.4)
comp_curve2_params <- list(linf = 65, k = 0.25, t0 = -1)
comp_curve_iq <-
 integral_quotient(ref_curve_params, comp_curve_params, 0, 10)
comp_curve2_iq <-
  integral_quotient(ref_curve_params, comp_curve2_params, 0, 10)
vbgf <- function (x, linf, k, t0) {linf * (1 - exp(-k * (x - t0)))}
curve(
  vbgf(x, ref_curve_params$linf, ref_curve_params$k, ref_curve_params$t0),
  from = 0,
  to = 10,
  ylim = c(0, 60),
  xlab = "Age", ylab = "Length"
)
curve(
  vbgf(x, comp_curve_params$linf, comp_curve_params$k, comp_curve_params$t0),
  add = TRUE,
  col = "blue"
)
curve(
  vbgf(x, comp_curve2_params$linf, comp_curve2_params$k, comp_curve2_params$t0),
  add = TRUE,
  col = "red"
)
text(9, 40, labels = paste0(comp_curve_iq, "\%"), col = "blue")
text(9, 43, labels = paste0(comp_curve2_iq, "\%"), col = "red")
}
