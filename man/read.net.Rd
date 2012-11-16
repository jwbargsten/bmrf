\name{read.net}
\alias{read.net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.net(f, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{
%%     ~~Describe \code{f} here~~
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
function (f, verbose = FALSE) 
{
    data = read.table(f, sep = "\t")
    protein.ids.uniq = sort(unique(c(levels(data[, 1]), levels(data[, 
        2]))))
    idcs = seq(along = protein.ids.uniq)
    names(idcs) = protein.ids.uniq
    l = length(protein.ids.uniq)
    protein.ids.a = as.character(data[, 1])
    protein.ids.b = as.character(data[, 2])
    net = sparseMatrix(i = idcs[protein.ids.a], j = idcs[protein.ids.b], 
        x = 1, dims = c(l, l))
    net = net + t(net)
    rownames(net) = protein.ids.uniq
    colnames(net) = protein.ids.uniq
    if (verbose) 
        cat("reading network OK\n")
    return(list(net = net, protein.idcs = idcs))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
