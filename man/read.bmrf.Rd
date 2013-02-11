\name{read.bmrf}
\alias{read.bmrf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
read bmrf from files
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.bmrf(net.file, go.file, fd.file, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net.file}{
%%     ~~Describe \code{net.file} here~~
}
  \item{go.file}{
%%     ~~Describe \code{go.file} here~~
}
  \item{fd.file}{
%%     ~~Describe \code{fd.file} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (net.file, go.file, fd.file, verbose = FALSE, ...) 
{
    neti <- read.net(net.file, verbose = verbose)
    go <- read.terms(go.file, neti$protein.idcs, verbose = verbose)
    fd <- read.terms(fd.file, neti$protein.idcs, verbose = verbose)
    return(bmrf(neti$net, go$mat, fd$mat, ...))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
