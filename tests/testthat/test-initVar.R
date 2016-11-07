# library(testthat)
# library(lava.reduce)
# library(lava)

context("#### initVar #### \n")

initVar_link <-  lava.reduce:::initVar_link

lava.options(symbol = c("~",","))
initVar_link(var1 = a~b)
initVar_link(var1 = a ~ b)
initVar_link(var1 = "a ~ b", format = "formula")
initVar_link(var1 = a ~ b+c+d*e, format = "list")
initVar_link(var1 = a ~ b+c+d*e, format = "txt.formula")
initVar_link(var1 = a ~ b+c+d*e, format = "formula")

initVar_link(var1 = "a,b")
initVar_link(var1 = "a", var2 = "b")

initVar_link(var1 = Y~X1+X2)
initVar_link(var1 = Y~X1+X2, repVar1 = TRUE)
initVar_link(var1 = Y~X1+X2, format = "formula")
initVar_link(var1 = Y~X1+X2, format = "txt.formula")

lava.options(symbols = c("<-","<->"))
initVar_link(var1 = "Y<-X1+X2", repVar1 = TRUE)
initVar_link(var1 = "Y<-X1+X2", format = "formula")
initVar_link(var1 = "Y<-X1+X2", format = "txt.formula")

initVar_link(var1 = "X1")
initVar_link(var1 = ~X1)
initVar_link(var1 = ~X1+X2)
