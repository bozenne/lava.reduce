# setwd(file.path(butils:::path_gitHub(),"lavaReduce","tests"))


library("lavaReduce")
suppressPackageStartupMessages(library("testthat"))
source("FCT.R")
test_check("lavaReduce")


