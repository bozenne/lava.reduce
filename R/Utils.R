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
    var2 <- all.vars(stats::delete.response(stats::terms(var1)))
    n.var2 <- length(var2)
    var1 <- setdiff(all.vars(var1),var2)
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
 
  #### export ####
  return(res)
}



#' @title Response variable of a formula
#' @description Return the reponse variable contained in the formula
#' 
#' @param formula a formula
#' 
#' @examples 
#' select.response(Y1~X1+X2)
#' 
#' @export
select.response <- function(formula){
  return(
    setdiff(all.vars(formula),
            all.vars(delete.response(terms(formula))))
  )
}

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

