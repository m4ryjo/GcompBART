\name{GcompBART}
\alias{GcompBART}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
G-computation based on BART
}
\description{
Implements G-computation based on Bayesian Additive Regression Trees (BART). Effect estimates of a (time-varying) treatment/exposure can be obtained using both cross-sectional and longitudinal data with possibly time-varying confounding. BART is  a Bayesian "sum-of-trees" model and provide flexible nonparametric modeling of covariates for continuous and binary outcomes. The package accommodates binary or continuous outcomes possibly at multiple time points. The package also provide mortal-cohort inferences on those who are still alive at a specific time point (partly conditional inference). For details see Josefsson and Daniels, 2021 and Josefsson, Daniels and Pudas, 2021.
}
\usage{
GcompBART(
  data, var.type, exp.reg, cont.reg, J=2000, Ndraws=200, Suppress = TRUE, By = Rep/10, ...
  )
}

\arguments{
  \item{data}{
  A matrix or data frame possibly containing confounder(s), exposure(s), outcome(s) and mortality indicator(s) in temporal order from left to right.
  }

  \item{var.type}{
  Vector of variable specifications for data. C=confounder, Fi=fixed, O=outcome, S=survival
  }

  \item{exp.reg}{
  A vector specifying the regime for exposed/treated.
  }

  \item{cont.reg}{
  A vector specifying the contrasting regime for unexposed/untreated.
  }

  \item{J}{
  Size of pseudo data. Default is set to 2,000.
  }

  \item{Ndraws}{
  Number of posterior draws. Default is set to 200.
  }

  \item{Suppress}{
  Indicate if the output should be suppressed. Default is TRUE
  }

  \item{By}{
  If Suppress is set to FALSE, output is provided for the By:th iteration.
  }

  \item{...}{
  Additional arguments passed to pbart or wbart.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\code{GcompBART} returns a list, with the following components:

  \item{summary}{
    A matrix with summary statistics, i.e. means and 95 \% creadible intervals, where each row corresponds to a specific regime and a time point.
   }
  \item{y_hat }{A matrix with (2 * the number of outcome measures) rows and Ndraw columns. Each column corresponds to a draw from the posterior distribution averged over the pseudo sample and each row corresponds to a specific regime and a time point. Burn-in is dropped.}

}
\references{
Chipman, H., George, E., and McCulloch R. (2010)
   Bayesian Additive Regression Trees.
   \emph{The Annals of Applied Statistics}, \bold{4,1}, 266-298 <doi:10.1214/09-AOAS285>.

Chipman, H., George, E., and McCulloch R. (2006)
   Bayesian Ensemble Learning.
   Advances in Neural Information Processing Systems 19,
   Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.

Friedman, J.H. (1991)
   Multivariate adaptive regression splines.
   \emph{The Annals of Statistics}, \bold{19}, 1--67.

Josefsson, M., Daniels, M.J., and Pudas, S. (2021)
  A Bayesian semi-parametric approach for inference on the population partly conditional mean from longitudinal data with dropout.
  \emph{Biostatistics}, Epub ahead of print.

Josefsson, M., and Daniels, M.J. (2021)
  Bayesian semi-parametric G-computation for causal inference in a cohort study with non-ignorable dropout and death.
  \emph{Journal of the Royal Statistical Society - Series C}, 70(2), 398-414
}
\author{
Maria Josefsson
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##simulate data (modified example from Friedman, 1991)
f <- function(x){
  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma <- 1.0  #y = f(x) + sigma*z , z~N(0,1)
n <- 100      #number of observations
set.seed(99)
x <- matrix(c(runif(n*4), rbinom(n, 1, prob = 0.5), runif(n*5)), n, 10) #10 variables, only first 5 matter
y <- f(x) + rnorm(n)

##test GcompBART
set.seed(99)
gcbfit <- GcompBART(data = cbind(x, y), var.type = c(rep("C", 4), "Fi", rep("C", 5), "O"), exp.reg = 1, cont.reg = 0, J = 100, Ndraws = 10, Suppress = TRUE)
gcbfit$summary

\dontrun{
##run GcompBART
  gcbfit <- GcompBART(data = cbind(x, y), var.type = c(rep("C", 4), "Fi", rep("C", 5), "O"), exp.reg = 1, cont.reg = 0, J = 100, Ndraws = 100, Suppress = TRUE)
  gcbfit$summary
}
}
\keyword{gcomputation}
\keyword{regression standardization}
\keyword{nonparametric}
