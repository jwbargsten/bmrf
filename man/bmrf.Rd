\name{bmrf}
\alias{bmrf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bmrf(net, go, fd, burnin = 20, niter = 20, minGOsize = 20, maxGOsize = (0.9 * nrow(net)), minFDsize = 20, maxFDsize = (0.9 * nrow(net)), verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{
%%     ~~Describe \code{net} here~~
}
  \item{go}{
%%     ~~Describe \code{go} here~~
}
  \item{fd}{
%%     ~~Describe \code{fd} here~~
}
  \item{burnin}{
%%     ~~Describe \code{burnin} here~~
}
  \item{niter}{
%%     ~~Describe \code{niter} here~~
}
  \item{minGOsize}{
%%     ~~Describe \code{minGOsize} here~~
}
  \item{maxGOsize}{
%%     ~~Describe \code{maxGOsize} here~~
}
  \item{minFDsize}{
%%     ~~Describe \code{minFDsize} here~~
}
  \item{maxFDsize}{
%%     ~~Describe \code{maxFDsize} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (net, go, fd, burnin = 20, niter = 20, minGOsize = 20, 
    maxGOsize = (0.9 * nrow(net)), minFDsize = 20, maxFDsize = (0.9 * 
        nrow(net)), verbose = FALSE) 
{
    go.cs <- colSums(go)
    go <- go[, go.cs >= minGOsize & go.cs <= maxGOsize]
    fd.cs <- colSums(fd)
    fd <- fd[, fd.cs >= minFDsize & fd.cs <= maxFDsize]
    u = which(rowSums(go) == 0)
    new("BMRF", net = net, go = go, fd = fd, burnin = as.integer(burnin), 
        niter = as.integer(niter), maxGOsize = maxGOsize, minGOsize = minGOsize, 
        maxFDsize = maxFDsize, minFDsize = minFDsize, unknown.idcs = u, 
        verbose = verbose)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
