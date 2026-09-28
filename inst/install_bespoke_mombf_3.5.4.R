try(detach("package:mombf", unload=TRUE))
try(detach("package:cli", unload=TRUE))
try(detach("package:sparseMatrixStats", unload = TRUE))
try(remove.packages("mombf"))
# .rs.restartR() # can help

# dependency "sparseMatrixStats" available from Bioconductor
remotes::install_bioc("sparseMatrixStats")

Rcpp::compileAttributes(pkgdir = "./inst/mombf_3.5.4/mombf")
remotes::install_local("./inst/mombf_3.5.4/mombf", force = TRUE)
library(mombf)
