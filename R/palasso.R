
#--- Workhorse function --------------------------------------------------------

#' @title
#' Paired lasso
#' 
#' @export
#' 
#' @description
#' The function \code{palasso} fits the paired lasso. Use this regression
#' technique if the covariates are numerous and occur in pairs.
#' 
#' @param y
#' response\strong{:}
#' vector of length \eqn{n},
#' 
#' @param X
#' covariates\strong{:}
#' list of matrices,
#' each with \eqn{n} rows (samples)
#' and \eqn{p} columns (variables)
#' 
#' @param ...
#' further arguments for \code{\link[glmnet]{cv.glmnet}} or
#' \code{\link[glmnet]{glmnet}}
#' 
#' @details
#' Let \code{x} denote one entry of the list \code{X}. See \link[glmnet]{glmnet}
#' for alternative specifications of \code{y} and \code{x}. Among the further
#' arguments, \code{family} must equal \code{"gaussian"}, \code{"binomial"},
#' \code{"poisson"}, or \code{"cox"},
#' and \code{penalty.factor} must not be used.
#' \emph{This method (hopefully) improves prediction if
#' \code{alpha}\eqn{=1} (lasso)
#' and \code{pmax}\eqn{<<p} (non-zero coefficients).}
#' 
#' @return
#' This function returns an object of class \code{palasso}.
#' Available methods include
#' \code{\link[=predict.palasso]{predict}},
#' \code{\link[=coef.palasso]{coef}},
#' \code{\link[=weights.palasso]{weights}},
#' \code{\link[=fitted.palasso]{fitted}},
#' \code{\link[=residuals.palasso]{residuals}},
#' \code{\link[=deviance.palasso]{deviance}},
#' \code{\link[=logLik.palasso]{logLik}},
#' and \code{\link[=summary.palasso]{summary}}.
#' Click \code{\link[=extra]{here}} for hidden functions.
#' 
#' @references
#' A Rauschenberger, RX Menezes, MA Jonker, and MA van de Wiel (2018).
#' "Penalised regression with paired covariates."
#' \emph{Manuscript in preparation.} \email{a.rauschenberger@vumc.nl}
#' 
#' @examples
#' set.seed(1)
#' n <- 40; p <- 1000
#' y <- rbinom(n=n,size=1,prob=0.5)
#' X <- lapply(1:2,function(x) matrix(rnorm(n*p),nrow=n,ncol=p))
#' object <- palasso(y=y,X=X,family="binomial",pmax=10)
#' 
palasso <- function(y,X,...){

    # checks
    base <- list(...)
    funs <- list(glmnet::glmnet,glmnet::cv.glmnet)
    formals <- unlist(lapply(funs,function(x) formals(x)),recursive=FALSE) # changed!
    if(any(!names(base) %in% names(formals))){stop("Invalid argument.")}
 
    # arguments
    base$y <- y
    base$x <- do.call(what="cbind",args=X)
    default <- list(family="gaussian",alpha=1,nfolds=10,type.measure="deviance")
    base <- c(base,default[!names(default) %in% names(base)])
    if(!base$family %in% c("gaussian","binomial","poisson","cox")){
        stop("Invalid argument \"family\".")
    }
    if(!is.null(base$pmax) & base$alpha==0){
        stop("Unexpected argument \"pmax\" as \"alpha=0\".")
    }
    
    # dimensionality
    k <- ifelse(is.list(X),length(X),1)
    if(k==1){
        stop("Invalid argument \"X\".")
    }
    n <- nrow(X[[1]])
    p <- ncol(X[[1]])
    
    # penalty factor (trial)
    if(!is.null(base$penalty.factor)){
        stop("Unexpected argument \"penalty.factor\".")
    }
    
    # distribution (trial)
    if(survival::is.Surv(y)){
        guess <- "cox"
    } else {
        guess <- "gaussian"
        guess[all(base$y%%1==0 & base$y>=0)] <- "poisson"
        guess[!is.vector(base$y) | length(unique(base$y))==2] <- "binomial"
    }
    if(guess!=base$family){
        warning(paste0("Consider family \"",guess,"\"."))
    }
    
    # fold identifier
    cond <- logical()
    cond[1] <- is.null(base$foldid)
    cond[2] <- base$family=="binomial"
    cond[3] <- is.vector(y)
    if(all(cond)){
        base$foldid <- rep(NA,times=n)
        base$foldid[y==0] <- sample(rep(seq_len(base$nfolds),
            length.out=sum(y==0)))
        base$foldid[y==1] <- sample(rep(seq_len(base$nfolds),
            length.out=sum(y==1)))
    }
    
    # marginal effects (correlation)
    mar <- list()
    for(i in seq_len(k)){
        if(base$family=="cox"){
            mar[[i]] <- apply(X[[i]],2,function(x) abs(2*survival::survConcordance(y~x)$concordance-1))
        } else {
            mar[[i]] <- suppressWarnings(as.vector(abs(stats::cor(X[[i]],y))))
        }
        mar[[i]][is.na(mar[[i]])] <- 0
    }
    
    # # marginal effects (univariate regression)
    # family <- eval(parse(text=base$family))()
    # mar <- list()
    # for(i in seq_len(k)){
    #     mar[[i]] <- abs(apply(X[[i]],2,function(x) stats::glm.fit(y=y,x=cbind(1,x),family=family)$coefficients[2]))
    #     mar[[i]][is.na(mar[[i]])] <- 0
    # }
    ## for object x of class "summary.glm":
    ## x$coefficients[,"Estimate"]/x$coefficients[,"Std. Error"]
    ## for object x of class "coxph":
    ## x$coefficients/sqrt(x$var)
    
    weight <- model <- list()
    
    # standard lasso
    for(i in seq_len(k)){
        weight[[i]] <- rep(1*(seq_len(k)==i),each=p)
    }
    weight[[k+1]] <- rep(1/k,times=k*p)
    
    # weighted lasso
    weight[[k+2]] <- rep(1/k*sapply(mar,mean)/mean(unlist(mar)),each=p)
    weight[[k+3]] <- unlist(mar)/rowSums(do.call(cbind,mar))
    weight[[k+3]][is.na(weight[[k+3]])] <- 0
    
    # adaptive lasso
    for(i in seq_len(k)){
        weight[[k+3+i]] <- weight[[i]]*mar[[i]]
    }
    weight[[2*k+4]] <- unlist(mar)
    
    # cross-validation
    args <- base
    for(i in seq_along(weight)){
        args$penalty.factor <- 1/weight[[i]]
        ### start original ###
        #model[[i]] <- tryCatch(expr=do.call(what=glmnet::cv.glmnet,args=args),
        #                       error=function(x) NA)
        #if(class(model[[i]])!="cv.glmnet"){
        #    # intercept-only model
        #    temp <- base
        #    temp$lambda <- c(99e99,99e98)
        #    model[[i]] <- do.call(what=glmnet::cv.glmnet,args=temp)
        #}
        ### end original ###
        model[[i]] <- palasso:::.cv.glmnet(args) ### trial
        if(i > 1){
            # free memory
            model[[i]]$glmnet.fit$call$x <- NULL
        }
    }
    
    # names
    if(k==2){
        names = c("x","z") 
    } else {
        names = letters[seq_len(k)]
    }  
    all = paste0(names,collapse="")
    names(model) = c(paste0("standard_",c(names,all)),
                     paste0(c("between_","within_"),all),
                     paste0("adaptive_",c(names,all)))
    
    # output
    call <- sapply(list(...),function(x) deparse(x))
    attributes(model)$info <- list(n=n,k=k,p=p,names=names,call=call)
    class(model) <- "palasso"
    return(model)
}

.error <- function(x,args){
    pattern <- c("Error in predmat\\[which, seq\\(nlami\\)\\] <- preds",
                 "replacement has length zero")
    cond <- sapply(X=pattern,FUN=function(p) grepl(pattern=p,x=x))
    if(all(cond)){
        warning("Fitting intercept-only model.")
        args$lambda <- c(99e99,99e98)
        do.call(what=glmnet::cv.glmnet,args=args)
    } else {
        stop(x)
    }
}

.warning <- function(x){
    pattern <- c("from glmnet Fortran code \\(error code",
                 "Number of nonzero coefficients along the path exceeds pmax=",
                 "lambda value; solutions for larger lambdas returned")
    cond <- sapply(X=pattern,FUN=function(p) grepl(pattern=p,x=x))
    if(all(cond)){
        invokeRestart("muffleWarning")
    }
}

.cv.glmnet <- function(args){
    withCallingHandlers(expr=tryCatch(expr=do.call(what=glmnet::cv.glmnet,args=args),
                        error=function(x) palasso:::.error(x,args)),
                        warning=function(x) palasso:::.warning(x))
}
