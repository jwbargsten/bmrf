#!/usr/bin/env Rscript

tpl_id <- "AT1G15750"

ai <- read.table("AI_interactions.txt", header=TRUE, sep="\t")

direct_int <- ai[ai[,"ida"] == tpl_id | ai[,"idb"] == tpl_id, ]
direct_ids <- unique(as.character(direct_int[,"ida"]), as.character(direct_int[,"idb"]))

two_lvl <- lapply(direct_ids, function(id, ai) {
  int <- ai[ai[,"ida"] == id | ai[,"idb"] == id, ]
  int
}, ai)

z <- Reduce(function(x, y) merge(x, y, all=T), two_lvl, accumulate=F)
write.table(cbind(z[,c("ida", "idb")],1), "topless_2lvl_int.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
