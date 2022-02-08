\name{pwbart_it}
\alias{pwbart_it}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predicting new observations from one posterior draw with a previously fitted BART model
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pwbart_it(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
x.test,	#
   treedraws,	# treedraws from bart
   mu=0,	#mean to add on
   it, #it
   transposed = FALSE

  \item{x.test}{
x matrix to predict at
}

\item{treedraws}{
treedraws from bart model
}

\item{mu}{
mean to add on
}

\item{it}{
which posterior draw to predict}

\item{transposed}{
If x.train is transposed, transposed=TRUE is default
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
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
f <- function(x){
  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma <- 1.0  #y = f(x) + sigma*z , z~N(0,1)
n <- 100      #number of observations
set.seed(99)
x <- matrix(c(runif(n*4), rbinom(n, 1, prob = 0.5), runif(n*5)), n, 10) #10 variables, only first 5 matter
y <- f(x) + rnorm(n)

set.seed(123)
x.test = matrix(c(runif(10*10)),10,10)
post <- wbart(x.train = x, y.train = y, ndpost = 5)
pwbart_it(x.test, post$treedraws, mu=mean(y), it = 2)

