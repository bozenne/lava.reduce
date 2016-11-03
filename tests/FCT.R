checkMoment <- function(lvmRed, elvm, param = coef(elvm), data = elvm$data$model.frame){
  
  indexRED <- na.omit(match(coef(lvmRed),names(param)))
  
  ## test logLik
  expect_equal(gaussianLP_logLik.lvm(lvmRed, p = param[indexRED], data = data),
               lava:::gaussian_logLik.lvm(object = elvm, data=data, p = param) 
  )
  
  ## test gradient
  expect_equal(unname(gaussianLP_gradient.lvm(lvmRed, p = param[indexRED], data = data, indiv = FALSE)),
               lava:::gaussian_gradient.lvm(x = elvm$model, data= data, p=param, S = elvm$S, n = elvm$data$n, mu = elvm$mu)[indexRED]
  )
  
  ## test score
  expect_equal(unname(gaussianLP_score.lvm(lvmRed, p = param[indexRED], data = data, indiv = FALSE)),
               unname(lava:::gaussian_score.lvm(x = elvm$model, data=data, p=param, S = elvm$S, n = elvm$data$n, mu = elvm$mu)[,indexRED,drop=FALSE])
  )  
  expect_equal(unname(gaussianLP_score.lvm(lvmRed, p = param[indexRED], data = data, indiv = TRUE)),
               unname(score(elvm$model,data=data,p=param,indiv=TRUE)[,indexRED,drop=FALSE])
  )
  expect_equal(unname(gaussianLP_score.lvm(lvmRed, p = param[indexRED], data = data, indiv = FALSE)),
               unname(score(elvm$model,data=data, p = param,indiv=FALSE)[,indexRED,drop=FALSE])
  )  
  
  ## test hessian
  Hred_E <- gaussian2LP_hessian.lvm(x = lvmRed, data=data, p=param[indexRED])
  Hlava <- lava:::gaussian2_hessian.lvm(x = elvm$model, data=data, p=param, n = elvm$data$n, mu = elvm$mu, S = elvm$S)
  expect_equal(attr(Hred_E,"grad"),
               attr(Hlava,"grad")[indexRED]
  ) 
  attr(Hred_E,"grad") <- NULL
  attr(Hlava,"grad") <- NULL
  expect_equal(unname(Hred_E), unname(Hlava[indexRED,indexRED]))
  
  Hred_num <- gaussian1LP_hessian.lvm(x = lvmRed, data=data, p=param[indexRED])
  Hlava <- lava:::gaussian1_hessian.lvm(x = elvm$model,  p=param, n = elvm$data$n, mu = elvm$mu, S = elvm$S)
  expect_equal(unname(Hred_num),unname(Hlava[indexRED,indexRED]), tolerance = 1e-6)
  
  Hred_I <- gaussianLP_hessian.lvm(x = lvmRed, data=data, p=param[indexRED], type = "information")
  Hlava <- lava:::gaussian_hessian.lvm(x = elvm$model, data=data, p=param, n = elvm$data$n, mu = elvm$mu, S = elvm$S)
  # Hred_I-Hlava[indexRED,indexRED]
  
  return(invisible(list(Hred_E = Hred_E,
                        Hred_num = Hred_num,
                        Hred_I = Hred_I,
                        Hlava = Hlava)))
}