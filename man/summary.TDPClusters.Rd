% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TDPClusters.R
\name{summary.TDPClusters}
\alias{summary.TDPClusters}
\title{Summarizing cluster information and generating summary table for found clusters.}
\usage{
summary(object, ..., rest = FALSE)
}
\value{
Returns a summary table reporting Size, TDN (lower bound), #TrueNull (upper bound), TDP (lower bound), maximum statistic and the corresponding ID for each found cluster.
}
\description{
\code{summary} method for class \code{\link{TDPClusters}}.
}
\details{
If the output is not assigned, the summary table will be printed on console.
}
\section{Slots}{

\describe{
\item{\code{object}}{Object of class "TDPClusters". Stores the \code{\link{TDPClusters}} object, usually, a result of a call to \code{\link{TDPQuery}}.}

\item{\code{rest}}{Object of class "logical". By default, \code{rest = FALSE} indicates that information on the rest nodes will not be shown in the summary table.}
}}

\examples{

table <- summary(tdpclusters)      # get summary table (without rest)
summary(tdpclusters, rest = TRUE)  # print summary table (with rest)

}
\author{
Xu Chen, Thijmen Krebs, Wouter Weeda.
}
