% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{assign_alk_attributes}
\alias{assign_alk_attributes}
\title{Assign associated age-length key attributes to a data.frame}
\usage{
assign_alk_attributes(
  data,
  size_col = "length",
  age_col = "age",
  autobin = TRUE,
  size_bin = 1,
  min_age = NULL,
  plus_group = NULL,
  alk_n = NULL,
  classes = "alk",
  dnorm_params = NULL,
  levels = NULL
)
}
\arguments{
\item{data}{A data.frame}

\item{size_col}{Character. Name of the column representing sizes}

\item{age_col}{Character. Name of the column representing ages}

\item{autobin}{Logical to set the attribute of autobin}

\item{size_bin}{Numeric. What is the width of size bins}

\item{min_age}{Numeric. The minimum age that was included in the alk}

\item{plus_group}{Numeric. The age that represents the plus group}

\item{alk_n}{Numeric. The number of samples that went into creating the alk}

\item{classes}{Character. The class that should get prepended to the
data.frame class(es)}

\item{dnorm_params}{The value of parameters that went into creating the
normal distributions on the age groups}

\item{levels}{Character vector of the levels used. This creates the "levels"
attribute if present}
}
\value{
A data.frame with associated attributes assigned
}
\description{
This is just a helper function to assign the needed attributes and classes
to a data.frame that is produced by either \code{\link{make_alk}} or
\code{\link{make_halk}}.
}
