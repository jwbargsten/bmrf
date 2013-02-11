\name{read.terms}
\alias{read.terms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
read terms from file (ipr or go)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.terms(f, p.idcs, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{
%%     ~~Describe \code{f} here~~
}
  \item{p.idcs}{
%%     ~~Describe \code{p.idcs} here~~
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
function (f, p.idcs, verbose = FALSE) 
{
    data = read.table(f, sep = "\t")
    terms.uniq = sort(levels(data[, 2]))
    t.idcs = seq(along = terms.uniq)
    names(t.idcs) = terms.uniq
    protein.ids = as.character(data[, 1])
    terms = as.character(data[, 2])
    L = sparseMatrix(i = p.idcs[protein.ids], j = t.idcs[terms], 
        x = 1, dims = c(length(p.idcs), length(t.idcs)))
    L@x[] = 1
    rownames(L) = names(p.idcs)
    colnames(L) = names(t.idcs)
    if (verbose) 
        cat("reading terms OK\n")
    return(list(mat = L, term.idcs = t.idcs))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
