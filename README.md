# BMRF

Implementation of Bayesian Markov Random Field-based protein function prediction algorithm.

The package is pre-alpha-stage, so take care.

You can find more information on <http://www.ab.wur.nl/bmrf>.

Installation goes like this:

in R:
```
> install.packages(c("glmnet", "mvtnorm", "brglm"))

in the shell:

$ g clone https://github.com/jwbargsten/bmrf.git
$ cd bmrf
$ R CMD INSTALL .

```
An example is in:
`tests/test_basic.R`
