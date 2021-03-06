% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bopr.R
\name{online_leraning.matrix}
\alias{online_leraning.matrix}
\title{Online learning with stream of samples.}
\usage{
online_leraning.matrix(object, x = NULL, y = NULL, allow.new = TRUE)
}
\arguments{
\item{object}{\code{BOPR} object}

\item{x}{a  matrix of predictors.}

\item{y}{a factor vector with 2 level}

\item{allow.new}{allow new levels of current variable}
}
\value{
S3 \code{BOPR} object; a list of consisting of
 \item{beta_matrix}{beta matrix with mean and variance}
 \item{beta}{scaling parameter}
 \item{prior_prob}{prior of initial parametes}
 \item{epsilon}{parameter to apply dynamics}
}
\description{
Online learning with stream of samples.
}

