\name{BMRF-class}
\Rdversion{1.1}
\docType{class}
\alias{BMRF-class}
\alias{show,BMRF-method}

\title{Class \code{"BMRF"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("BMRF", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{net}:}{Object of class \code{"Matrix"} ~~ }
    \item{\code{go}:}{Object of class \code{"Matrix"} ~~ }
    \item{\code{fd}:}{Object of class \code{"Matrix"} ~~ }
    \item{\code{burnin}:}{Object of class \code{"integer"} ~~ }
    \item{\code{niter}:}{Object of class \code{"integer"} ~~ }
    \item{\code{minGOsize}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{maxGOsize}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{minFDsize}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{maxFDsize}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{unknown.idcs}:}{Object of class \code{"vector"} ~~ }
    \item{\code{verbose}:}{Object of class \code{"logical"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "BMRF")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("BMRF")
}
\keyword{classes}
