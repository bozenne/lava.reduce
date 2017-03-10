#' @title Remove variables from linear predictors
#' @description Remove one or several variable from the linear predictors
#' @name killLP
#' 
#' @param x \code{lvm}-object
#' @param value the name of the links that should be removed
#' @param lp should the variables be removed from the linear predictor
#' @param expar should the external parameters be also removed
#' @param restaure should the link be kept while removing the linear predictor
#' @param clean should the lvm object be simplified using the \code{clean} function
#' @param ... argument passed to \code{clean}
#' 
#' @examples
#'
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y")
#' rm <- reduce(m, rm.exo = FALSE)
#'
#' kill(rm, lp = FALSE, value = "x1")
#' kill(rm, lp = FALSE, value = c("x2","x3"))
#' kill(rm, lp = FALSE, value = ~x4+x5)
#'
#'
#' m <- lvm()
#' regression(m) <- y ~ x1+x2+x3
#' cancel(m) <- ~x1+x2+x3
#' cancel(m, c("x1","x2","x3"))
#' cancel(m) <- y~x1+x2+x3
#' m
#'
#' 
#' # see test/testthat/test-cancel.R

# {{{ kill.lvm.reduce
#' @rdname killLP
#' @export
kill.lvm.reduced  <- function(x, value = NULL, lp = TRUE, expar = TRUE, restaure = FALSE,
                              clean = TRUE, ...){

    ## normalize value
    if("formula" %in% class(value)){
        value <- select.regressor(value, type = "vars")
    }

    ## compatibility with lava
    if(is.null(lp(x))){
        
        return(callS3methodParent(x, FUN = "kill", class = "lvm.reduced", value = value))
        
    }else{  ## lava.reduce
        
        if(any(value %in% lp(x))){ ## remove linear predictors
            lp.rm <- value[value %in% lp(x)]
            x <- cancel(x, value = lp(x, lp = lp.rm, type = "link"), expar = expar, restaure = restaure, clean = FALSE)    
        }

        ## if the object is still lvm.reduce (could not be the case if all lp has been removed and simplify is TRUE)
        ## remove variables in the linear predictor
        if(lp && any(value %in% lp(x, type = "x"))){ 
            var.rm <- value[value %in% lp(x, type = "x")]
            index.x <- lapply(lp(x, type = "x", format = "list2"), function(var){which(var$x %in% var.rm)})
            name.lp <- lp(x, type = "name")
    
            for(iterLP in 1:length(index.x)){
                if(length(index.x)==0){next}
                x <- cancel(x, value = lp(x, lp = name.lp[iterLP], type = "link")[index.x[[iterLP]]], expar = expar,
                            restaure = restaure, clean = FALSE)
            }
    
        }

        ## remove variables outside the linear predictor
        if(any(value %in% vars(x, lp = FALSE, xlp = FALSE))){ 
            var.rm <- value[value %in% vars(x, lp = FALSE, xlp = FALSE)]
    
            x <- callS3methodParent(x, FUN = "kill", class = "lvm.reduced", value = var.rm)
    
        }
    }

    ## clean object
    if(clean){
        x <- clean(x, ...)
    }   
  
  return(x)
}

# }}}

# {{{ cancel.lvm.reduce
#' @name killLP
#' @export
cancel.lvm.reduced  <- function(x, value = NULL, expar = TRUE, restaure = FALSE,
                                clean = FALSE, ...){

    ## normalize value
    allCoef <- initVar_links(value, format = "txt.formula")
    
    ## order by linear predictor  
    if(any(allCoef %in% coef(x) == FALSE)){
        stop("unknown link: \"",paste(allCoef[allCoef %in% coef(x) == FALSE], collapse = "\" \""),"\" \n",
             "possible link: \"",paste(coef(x)[coef(x) %in%allCoef == FALSE], collapse = "\" \""),"\" \n")
    }
  
    allCoef.lp <- intersect(allCoef, lp(x, "link"))
    allCoef.nlp <- setdiff(allCoef, lp(x, "link"))
  
    ## compatibility with lava: links outside the linear predictor
    if(length(allCoef.nlp)>0){
        list.nlp <- initVar_links(allCoef.nlp, format = "list")
        
        ## not working because cancel.lvm calls cancel
        # f.nlp <- combine.formula(allCoef.nlp)
        for(iterF in 1:length(list.nlp$var2)){
            #   x <- callS3methodParent(x, FUN = "cancel", class = "lvm.reduced", 
            #                           value = f.nlp[[iterF]])Q
            x <- callS3methodParent(x, FUN = "cancel", class = "lvm.reduced", 
                                       value = c(list.nlp$var1[iterF],list.nlp$var2[iterF]))
        }
        
    }
  
    ## coefficients in LP
    if(length(allCoef.lp)>0){
        linkLP <- lp(x, type = "link", format = "list")        
        name.lp <- lp(x, type = "name")        

        for(iterLP in 1:length(linkLP)){ # iterLP <- 1
            iCoef.lp <- allCoef.lp[allCoef.lp %in% linkLP[[iterLP]]]
            if(length(iCoef.lp)==0){next} 
            
            ## remove external parameters from the LVM 
            if(expar){
                parameter(x, remove = TRUE) <- unlist(iCoef.lp) #must be a character
                # x$expar <- x$expar[setdiff(names(x$expar),coefLP)]
                # x$exfix <- x$exfix[setdiff(names(x$exfix),coefLP)]
                # if(length(x$exfix)==0){x$exfix <- NULL}
                # x$attributes$parameter <- x$attributes$parameter[setdiff(names(x$attributes$parameter),coefLP)]
                # x$index$npar.ex <- x$index$npar.ex[setdiff(names(x$index$npar.ex),coefLP)]
                # x$index$parval <- x$index$parval[setdiff(names(x$index$parval),coefLP)]
            }
      
            newlp <- lp(x, type = c("x","con","name","link"), lp = name.lp[iterLP], format = "list2")
            newlp.link <- newlp[[1]]$link
            newlp.x <- newlp[[1]]$x
      
            indexRM <- which(newlp.link %in% iCoef.lp)

            ## restaure the links removed from the lp in the lvm
            if(restaure){
                f <- as.formula(paste( lp(x, type = "endo", lp = iterLP), "~",  paste(newlp.x[indexRM], collapse = " + ") ) )
                regression(x) <- f
            }
      
            newlp[[1]]$link <- NULL # link is not part of the lp object. It has been constructed by the lp() extractor so it should be removed for updating using lp<-
            newlp[[1]]$con <- newlp[[1]]$con[-indexRM]
            newlp[[1]]$x <- newlp[[1]]$x[-indexRM]
      
            lp(x, lp = name.lp[iterLP]) <- newlp      
        }
    }

    ## clean object
    if(clean){
        x <- clean(x, ...)
    }    
  
    return(x)
}
# }}}

