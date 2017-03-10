# library(testthat)
# library(lava.reduce)
# library(lava)

context("#### destructors #### \n")

lava.options(symbols = c("~","~~"))

m <- lvm.reduced()
m <- regression(m, x=paste0("x",1:10),y="y1", reduce = TRUE)
m <- regression(m, x=paste0("x",1:10),y="y2", reduce = TRUE)
m <- regression(m, x=c("z1","z2"),y="y1")
m <- regression(m, x=c("z3"),y="y2")

# {{{ cancel
test_that("cancel coefficient not in lp (compatibility with lava)", {
    coef2RM <- c("y1~z1","y1~z2")
    formula2RM <- combine.formula(coef2RM)[[1]]
    
    mc <- cancel(m, formula2RM)    
    expect_false(any(coef2RM %in% coef(mc)))
    expect_true(all(setdiff(coef(m),coef2RM) %in% coef(mc)))

    mc2 <- cancel(m, coef2RM)
    mc3 <- m; cancel(mc3) <- formula2RM
    mc4 <- m; cancel(mc4) <- coef2RM    
    expect_equal(mc,mc2)
    expect_equal(mc,mc3)
    expect_equal(mc,mc4)  
})

test_that("cancel coefficient in lp", {
    coef2RM <- c("y1~x1","y1~x2")
    formula2RM <- combine.formula(coef2RM)[[1]]

    mc <- cancel(m, formula2RM)    
    expect_false(any(coef2RM %in% coef(mc)))
    expect_true(all(setdiff(coef(m),coef2RM) %in% coef(mc)))

    mc2 <- cancel(m, coef2RM)
    mc3 <- m; cancel(mc3) <- formula2RM
    mc4 <- m; cancel(mc4) <- coef2RM    
    expect_equal(mc,mc2)
    expect_equal(mc,mc3)
    expect_equal(mc,mc4)
})

test_that("move coefficient from lp to normal", {
    coef2RM <- c("y1~x1","y1~x2")
    formula2RM <- combine.formula(coef2RM)[[1]]

    mc <- cancel(m, formula2RM, restaure = TRUE)    
    expect_true(all(coef(m) %in% coef(mc)))
    expect_false(any(coef2RM %in% lp(mc, type = "link")))

    mc2 <- cancel(m, coef2RM, restaure = TRUE)
    mc3 <- m; cancel(mc3, restaure = TRUE) <- formula2RM
    mc4 <- m; cancel(mc4, restaure = TRUE) <- coef2RM    
    expect_equal(mc,mc2)
    expect_equal(mc,mc3)
    expect_equal(mc,mc4)
})

test_that("remove a complete lp", {
    lpName <-  lp(m, type = "name")[1]
    coef2RM <- lp(m, lp = lpName, type = "link")
    formula2RM <- combine.formula(coef2RM)[[1]]

    mc <- cancel(m, formula2RM, clean = TRUE)
    expect_false(lpName %in% lp(mc, type = "name"))

    mc2 <- cancel(m, coef2RM, clean = TRUE)
    mc3 <- m; cancel(mc3, clean = TRUE) <- formula2RM
    mc4 <- m; cancel(mc4, clean = TRUE) <- coef2RM    
    expect_equal(mc,mc2)
    expect_equal(mc,mc3)
    expect_equal(mc,mc4)

    lpName <-  lp(m, type = "name")
    coef2RM <- lp(m, lp = lpName, type = "link")

    mc <- cancel(m, coef2RM, clean = TRUE)    
    expect_false("lvm.reduced" %in% class(mc))
})
# }}}

# {{{ kill
test_that("kill variables outside lp (compatibility with lava)",{
    var2RM <- exogenous(m, lp = FALSE)
    formula2RM <- as.formula(paste0("~",paste(var2RM,collapse ="+")))
    
    m1 <- m
    kill(m1) <- formula2RM
    expect_false(any(var2RM %in% vars(m1)))

    m2 <- kill(m, var2RM)
    m3 <- m; kill(m3) <- formula2RM
    m4 <- m; kill(m4) <- var2RM    
    expect_equal(m1,m2)
    expect_equal(m1,m3)
    expect_equal(m1,m4)
})

test_that("kill variables in linear predictors",{
    var2RM <- lp(m, type = "x")[1:3]
    formula2RM <- as.formula(paste0("~",paste(var2RM,collapse ="+")))
    
    m1 <- m
    kill(m1) <- formula2RM
    expect_false(any(var2RM %in% lp(m1, type = "x")))
    expect_false(any(var2RM %in% vars(m1)))

    m2 <- kill(m, var2RM)
    m3 <- m; kill(m3) <- formula2RM
    m4 <- m; kill(m4) <- var2RM    
    expect_equal(m1,m2)
    expect_equal(m1,m3)
    expect_equal(m1,m4)

    m1 <- kill(m, unique(lp(m, type = "x")))
    expect_false("lvm.reduced" %in% class(m1))
})

test_that("kill linear predictors",{
    var2RM <- lp(m, type = "name")[1]
    formula2RM <- as.formula(paste0("~",var2RM))
    
    m1 <- m
    kill(m1) <- formula2RM
    expect_false(var2RM %in% lp(m1, type = "name"))
    expect_false(any(lp(m, lp = var2RM, type = "link") %in% coef(m1)))
    expect_true(all(setdiff(coef(m),lp(m, lp = var2RM, type = "link")) %in% coef(m1)))
    
    m2 <- kill(m, var2RM)
    m3 <- m; kill(m3) <- formula2RM
    m4 <- m; kill(m4) <- var2RM    
    expect_equal(m1,m2)
    expect_equal(m1,m3)
    expect_equal(m1,m4)

    m1 <- kill(m, lp(m, type = "name"))
    expect_false("lvm.reduced" %in% class(m1))
    expect_true(identical(vars(m1),vars(m, lp = FALSE)))    
})

# }}}

# {{{ clean
m1 <- lvm()
m1 <- regression(m1, x=paste0("x",1:5),y="y1")
m1 <- regression(m1, x=paste0("x",1:5),y="y2")
covariance(m1) <- y1~y2

m2 <- lvm(y1 ~ eta + x1, y2 ~ eta, y3 ~ eta + x2)
latent(m2) <- ~eta

test_that("do not kill when it should not", {
    expect_equal(clean(m1),m1)
    expect_equal(clean(m2),m2)
})

test_that("remove exogenous variables", {
 cancel(m1) <- y1 ~ x1
 cancel(m1) <- y2 ~ x1
 m10 <- clean(m1)
 expect_false("x1" %in% vars(m10))
 expect_true(all(setdiff(vars(m1),"x1") %in% vars(m10)))
})

test_that("remove latent and endogenous variables", {
    cancel(m2) <- y1 ~ eta
    cancel(m2) <- y2 ~ eta
    cancel(m2) <- y3 ~ eta
    m20 <- clean(m2)
    expect_false("eta" %in% vars(m20))
    expect_false("y2" %in% vars(m20))
    expect_true(all(setdiff(vars(m2),c("y2","eta")) %in% vars(m20)))
})
# }}}
