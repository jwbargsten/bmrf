setClass("BMRF", representation(
    net = "Matrix",
    go = "Matrix",
    fd = "Matrix",
    burnin = "integer",
    niter = "integer",
    minGOsize = "numeric",
    maxGOsize = "numeric",
    minFDsize = "numeric",
    maxFDsize = "numeric",
    unknown.idcs = "vector",
    verbose = "logical"
  )
)

     

