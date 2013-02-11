\name{predict.bmrf.glmnet}
\alias{predict.bmrf.glmnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 elastic net prediction
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.bmrf.glmnet(bmrf, go.idx, dfmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bmrf}{
%%     ~~Describe \code{bmrf} here~~
}
  \item{go.idx}{
%%     ~~Describe \code{go.idx} here~~
}
  \item{dfmax}{
%%     ~~Describe \code{dfmax} here~~
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (bmrf, go.idx, dfmax) 
{
    yf <- bmrf@go[-bmrf@unknown.idcs, go.idx]
    xf <- bmrf@fd[-bmrf@unknown.idcs, ]
    f = try(cv.glmnet(y = yf, x = xf, family = "binomial", dfmax = dfmax, 
        alpha = 0.5, maxit = 1000, type.measure = "auc"))
    if (class(f) == "try-error") 
        return(NULL)
    yhat = try(predict(f, bmrf@fd, s = "lambda.min"))
    names(Yhat) = rownames(bmrf@fd)
    if (class(yhat) == "try-error") 
        return(NULL)
    return(yhat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
