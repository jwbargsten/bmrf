proteins <- function(bmrf) {
  names = rownames(bmrf@net)
  uidx = bmrf@unknown.idcs

  return(list(names=names, uidx=uidx))
}
