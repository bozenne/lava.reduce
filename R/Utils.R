#' @title Normalize var1 and var2
#' @description Convert var1 and var2 from formula or covariance to character
#' 
#' @param var1 a character indicating the endogeneous variable or a formula
#' @param var2 an optional character indicating the exogeneous variable
#' @param repVar1 should var1 be duplicated to match var2 length. Only active if format = "list".
#' @param format should the name of the variable be return (format = "list"), a vector of character formula ("txt.formula") or a list of formula ("formula")
#' 
#' @details See test file test/testthat/test-Reduce.R for examples
#' @export
initVar_link <- function(var1, var2, repVar1 = FALSE, format = "list"){
  
  Slink <- c(lava.options()$symbol[1],lava.options()$symbols[1])
  pSlink <- paste(c(Slink,"~"), collapse = "|")
  Scov <- c(lava.options()$symbol[2],lava.options()$symbols[2])
  pScov <- paste(Scov, collapse = "|")

  if(missing(var2) && is.character(var1)){ 
    if(grepl(pScov,var1)==TRUE){ ## covariance
      Scov <- names(unlist(sapply(Scov, grep, x = var1)))[1]
      var1 <- gsub(Scov,"~", x = var1)
      sep <- Scov
    }
    if(grepl(pSlink,var1)==TRUE){ # regression
      Slink <- names(unlist(sapply(c(Slink,"~"), grep, x = var1)))[1]
      var1 <- stats::as.formula(gsub(Slink,"~",var1))
      sep <- if(format == "formula"){"~"}else{Slink}
    }
  }else{
    sep <- if(format == "formula"){"~"}else{Slink}
  }
  
  if(class(var1) == "formula"){
    var2 <- select.regressor(var1, type = "vars")
    n.var2 <- length(var2)
    var1 <- select.response(var1, type = "vars")
  }
  
  #### convert to format
  if(format == "formula"){
    n.var2 <- length(var2)
    var1 <- rep(var1, times = n.var2)
    res <- sapply(1:n.var2, function(i){stats::as.formula(paste(var1[i], var2[i], sep = sep))})
    
  }else if(format == "txt.formula"){
    n.var2 <- length(var2)
    var1 <- rep(var1, times = n.var2)
    res <- sapply(1:n.var2, function(i){paste(var1[i], var2[i], sep = sep)})
    
  }else if(format == "list"){
    if(repVar1 && !missing(var1)){var1 <- rep(var1, length(var2))}
    
    res <- list(var1 = var1,
                var2 = if(!missing(var2)){var2}else{NULL} 
    )
  }
 
  ## export 
  return(res)
}



#' @title Response variable of a formula
#' @description Return the reponse variable contained in the formula
#' @name select.response
#' 
#' @param x a formula
#' @param type either return an object of type call (\code{"call"}) or the names of the variables (\code{"vars"})
#' @param ... additional arguments to be passed to the low level functions
#' 
#' @examples 
#' select.response(Y1~X1+X2)
#' select.response(Y1~X1+X2, type = "vars")
#' 
#' select.response(Y1~X1+Y1)
#' select.response(Y1+Y2~X1+Y1, type = "vars")
#' 
#' select.response(~X1+X2)
#' select.response(~X1+X2, type = "vars")

#' @rdname select.response
#' @export
`select.response` <-  function(x,...) UseMethod("select.response")

#' @rdname select.response
#' @export
select.response.formula <- function(x, type = "call", ...){
  
  match.arg(type, c("call","vars"))
  
  if(length(x)==3){
    res <- x[[2]]
    if(type == "vars"){
      res <- all.vars(res)
    }
  }else{
    res <- NULL
  }
  
  return(res)
}

#' @title Regressor of a formula
#' @description Return the regressor variables contained in the formula
#' @name select.regressor
#' 
#' @param x a formula
#' @param type either return an object of type call (\code{"call"}) or the names of the variables (\code{"vars"})
#' @param ... additional arguments to be passed to the low level functions
#'   
#' @examples 
#' select.regressor(Y1~X1+X2)
#' select.regressor(Y1~X1+X2, type = "vars")
#' 
#' select.regressor(Y1~X1+Y1)
#' select.regressor(Y1+Y2~X1+Y1, type = "vars")
#' 
#' select.regressor(~X1+X2)
#' select.regressor(~X1+X2, type = "vars")


#' @rdname select.regressor
#' @export
`select.regressor` <-  function(x,...) UseMethod("select.regressor")

#' @rdname select.regressor
#' @export
select.regressor.formula <- function(x, type = "call", ...){
  
  match.arg(type, c("call","vars"))
  
  if(length(x)==3){
    res <- x[[3]]
    
  }else if(length(x)==2){
    res <- x[[2]]
  }else{
    res <- NULL
  }
  if(type == "vars"){
    res <- all.vars(res)
  }
  
  return(res)
}

#' @title Combine formula
#' @description Combine formula by outcome
#' 
#' @param ls.formula a list of formula
#' @param as.formula should as.formula be applied to each element of the list
#' @param as.unique should regressors appears at most once in the formula
#' 
#' @examples
#' combine.formula(list(Y~X1,Y~X3+X5,Y1~X2))
#' lava.options(symbol = c("~",","))
#' combine.formula(list("Y~X1","Y~X3+X5","Y1~X2"))
#' lava.options(symbol = c("<-","<->"))
#' combine.formula(list("Y<-X1","Y<-X3+X5","Y1<-X2"))
#' 
#' combine.formula(list(Y~X1,Y~X3+X1,Y1~X2))
#' combine.formula(list(Y~X1,Y~X3+X1,Y1~X2), as.unique = TRUE)
#' 
#' @export
combine.formula <- function(ls.formula, as.formula = TRUE, as.unique = FALSE){
  
  if(length(ls.formula)==0){return(NULL)}
  if(class(ls.formula)=="formula"){ls.formula <- list(ls.formula)}
  
  ls.Vars <- lapply(ls.formula, initVar_link)
  
  ls.endogeneous <- unlist(lapply(ls.Vars, "[[", 1))
  ls.X <- lapply(ls.Vars, "[[", 2)
  endogenous <- unique(ls.endogeneous)
  n.endogeneous <- length(endogenous)
  
  ls.formula2 <- vector(n.endogeneous, mode = "list")
  for(iterE in 1:n.endogeneous){
    X <- unlist(ls.X[which(ls.endogeneous==endogenous[iterE])])
    if(as.unique){X <- unique(X)}
    txt <- paste(endogenous[iterE],"~",paste(X, collapse = " + "))
    ls.formula2[[iterE]] <- as.formula(txt)
  }
  
  return(ls.formula2)
}


#### NOT USED ####

#' @title formula character conversion
#' @description Conversion of formula into character string or vice versa
#' @name convFormulaCharacter
#' 
#' @param f a formula.
#' @param txt a character string.
#' @param type should the normal formula operator be used (\code{"formula"}) or the one of lava.option (\code{"symbols"} or \code{"symbol"}).
#' 
#' @examples
#' formula2character(Y1~X1+X2)
#' formula2character(Y1~X1+X2, type = "symbol")
#' formula2character(Y1~X1+X2, type = "symbols")
#' 
#' character2formula("Y1~X1+X2")
#' 
#' m <- lvm(Y~X)
#' character2formula(coef(m)[2])

#' @rdname convFormulaCharacter
#' @export
formula2character <- function(f, type = "formula"){
  
  match.arg(type, choices = c("formula", "symbols", "symbol"))
  
  if(type == "formula"){
    txt <- paste(deparse(f), collapse = "+")
  }else {
    txt <- as.character(f)
    txt[1] <- lava.options()[[type]][1]
    txt <- paste(txt[2],txt[1],txt[3], sep = "")
  }
  
  return(gsub("[[:blank:]]","",txt))
  
}

#' @rdname convFormulaCharacter
#' @export
character2formula <- function(txt){
  
  txt <- gsub(paste(lava.options()$symbols[1],lava.options()$symbol[1],sep="|"),"~",txt)
  return(stats::as.formula(txt))
  
}

#' @title Find methods belonging to other classes
#' @description Find methods that apply to the object for each class of the object
#' 
#' @param object an object
#' @param FUN the name of a generic method
#' @param class the level of class
#' @param export should all methods be output (\code{"all"}) or only the first one (\code{"first"})
#' @param value should the method be returned. Else its name will be returned.
#' 
#' @examples
#' m <- lvm.reduced(Y~X)
#' getS3methodParent(m, "coef")
#' getS3methodParent(m, "coef", value = TRUE)
#' 
#' @export
getS3methodParent <- function(object, FUN, class = 1, export = "all", value = FALSE){
  
  match.arg(export, c("first","all"))
  
   object.class <- class(object)
   n.class <- length(object.class)
   if(n.class == 0){
     stop("\'object\' inherit from no class \n")
   }
   
   
   if(is.character(class)){
     if(class %in% object.class == FALSE){
       stop("\'object\' does not inherit of  class ",class," \n")
     }else{
       class <- object.class[-(1:which(object.class == class))]
     }
   }else{
     if(class %in% 0:(n.class-1) == FALSE){
       stop("\'object\' inherit from ",n.class," class(es)\n",
            "\'class\' cannot be ",class,"\n")
     }else{
       class <- object.class[-(1:class)]
     }
   }
   ls.method <- sapply(class, function(c){getS3method(f = FUN, class = c, optional = TRUE)})
   test <- unlist(lapply(ls.method, length))
   if(all(test == 0)){
     return(NULL)
   }else if(any(test == 0)){
     ls.method[test == 0] <- NULL
   }
   
   if(value == FALSE){
     available.method <- paste(FUN, class[test>0], sep = ".")
     if(export == "first"){
       return(available.method[1])
     }else if(export == "all"){
       return(available.method)
     }
   }else if(export == "first"){
     return(ls.method[[1]])
   }else if(export == "all"){
     return(ls.method)
   }
}

#' @title Call the first method inhereted by the object
#' @description Call the first method inhereted by the object
#' 
#' @param object an object
#' @param FUN the name of a generic method
#' @param class the level of class
#' @param ... additional arguments to be passed to the method that will be called
#' 
#' @examples
#' m <- lvm.reduced(Y~X)
#' callS3methodParent(m, "coef")
#' callS3methodParent(m, FUN = "coef", class = "lvm.reduced")
#' 
#' @export
callS3methodParent <- function(object, FUN, class = 1,...){
  
  return(getS3methodParent(object, FUN = FUN, class = class, export = "first", value = TRUE)(object, ...))
  
}