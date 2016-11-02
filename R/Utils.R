#' @title Normalize var1 and var2
#' @description Convert var1 and var2 from formula or covariance to character
#' 
#' @param var1 a character indicating the endogeneous variable or a formula
#' @param var2 an optional character indicating the exogeneous variable
#' @param repVar1 should var1 be duplicated to match var2 length. Only active if format = "list".
#' @param format should the name of the variable be return (format = "list"), a vector of character formula ("txt.formula") or a list of formula ("formula")
#' 
#' @details See test file test/testthat/test-Reduce.R for examples
#' 
initVar_link <- function(var1, var2, repVar1 = FALSE, format = "list"){
  
  Slink <- lava.options()$symbol[1]
  Scov <- lava.options()$symbol[2]
  
  if(missing(var2) && is.character(var1)){
    if(grepl(Scov,var1)==TRUE){
      var1 <- gsub(Scov,"~", x = var1)
      sep <- Scov
    }
    if(grepl(Slink,var1)==TRUE){
      var1 <- as.formula(gsub(Slink,"~",var1))
      sep <- if(format == "formula"){"~"}else{Slink}
    }
  }else{
    sep <- if(format == "formula"){"~"}else{Slink}
  }
  
  if(class(var1) == "formula"){
    var2 <- all.vars(delete.response(terms(var1)))
    n.var2 <- length(var2)
    var1 <- setdiff(all.vars(var1),var2)
  }
  
  #### convert to format
  if(format == "formula"){
    n.var2 <- length(var2)
    var1 <- rep(var1, times = n.var2)
    res <- sapply(1:n.var2, function(i){as.formula(paste(var1[i], var2[i], sep = sep))})
    
  }else if(format == "txt.formula"){
    n.var2 <- length(var2)
    var1 <- rep(var1, times = n.var2)
    res <- sapply(1:n.var2, function(i){paste(var1[i], var2[i], sep = sep)})
    
  }else if(format == "list"){
    if(repVar1){var1 <- rep(var1, length(var2))}
    
    res <- list(var1 = var1,
                var2 = var2)
  }
 
  #### export ####
  return(res)
}