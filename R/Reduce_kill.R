#' @title Remove variables from linear predictors
#' @description Remove one or several variable from the linear predictors
#' @name killLP
#' 
#' @param x \code{lvm}-object
#' @param var the names of the variables that should be removed
#' @param value the names of the variables that should be removed
#' @param expar should the external parameters be also removed
#' @param restaure should the link be kept while removing the linear predictor
#' @param ... argument passed to \code{clean}
#' 
#' @examples
#'
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y")
#' kill(m) <- "x7"
#' 
#' rm <- reduce(m, rm.exo = FALSE)
#'
#' kill(rm) <- ~x6
#' kill(rm, value = "x1")
#' kill(rm, value = c("x2","x3"))
#' kill(rm, value = ~x4+x5)
#' kill(rm, value = ~LPy)
#'
#'
#' m <- lvm.reduced()
#' regression(m) <- y1 ~ x1+x2+x3
#' regression(m) <- y2 ~ z
#' covariance(m) <- y1 ~ y2
#' cancel(m) <- y1~x1+x2+x3
#' m
#' cancel(m) <- "y1 ~~ y2"
#' coef(m)
#' 
#' # see test/testthat/test-cancel.R

# {{{ reduce.remove.hook
#' @rdname killLP
#' @export
lvmReduce.remove.hook  <- function(x, var, expar = TRUE, restaure = FALSE,
                                ...){

      ## lvm object with linear predictor
    if("lvm.reduced" %in% class(x)){
        ## remove linear predictors
        if(any(var %in% lp(x))){
            lp.rm <- var[var %in% lp(x)]
            x <- cancel(x, value = lp(x, lp = lp.rm, type = "link"), expar = expar, restaure = restaure, clean = FALSE)    
        }
           
        ## remove variables in the linear predictor
        if(any(var %in% lp(x, type = "x"))){ 
            var.rm <- var[var %in% lp(x, type = "x")]
            index.x <- lapply(lp(x, type = "x", format = "list2"), function(var){which(var$x %in% var.rm)})
            name.lp <- lp(x, type = "name")

            for(iterLP in 1:length(index.x)){
                if(length(index.x)==0){next}
                x <- cancel(x, value = lp(x, lp = name.lp[iterLP], type = "link")[index.x[[iterLP]]],
                            expar = expar, restaure = restaure, clean = FALSE)
            }
    
        }
      
    }
  
  return(x)
}

# }}}

# {{{ cancel.lvm.reduce
#' @name killLP
#' @export
cancel.lvm.reduced  <- function(x, value, expar = TRUE, restaure = FALSE,
                                ...){

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
  
    return(x)
}
# }}}

