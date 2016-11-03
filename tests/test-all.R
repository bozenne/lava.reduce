# setwd(file.path(butils:::path_gitHub(),"lava.reduce","tests"))


library("lava.reduce")
suppressPackageStartupMessages(library("testthat"))
source("FCT.R")
test_check("lava.reduce")


