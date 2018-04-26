
#--- Workhorse function --------------------------------------------------------

#' @title
#' Paired lasso
#' 
#' @export
#' 
#' @description
#' The function \code{palasso} fits the paired lasso. Use this function if you
#' have \emph{paired covariates} and want a \emph{sparse model}.
#' 
#' @param y
#' response\strong{:}
#' vector of length \eqn{n}
#' 
#' @param X
#' covariates\strong{:}
#' list of matrices,
#' each with \eqn{n} rows (samples)
#' and \eqn{p} columns (variables)
#' 
#' @param max
#' maximum number of non-zero coefficients\strong{:}
#' positive numeric, or \code{NULL} \eqn{(no sparsity constraint)}
#' 
#' @param ...
#' further arguments for \code{\link[glmnet]{cv.glmnet}} or
#' \code{\link[glmnet]{glmnet}}
#' 
#' @details
#' Let \code{x} denote one entry of the list \code{X}. See \link[glmnet]{glmnet}
#' for alternative specifications of \code{y} and \code{x}. Among the further
#' arguments, \code{family} must equal \code{"gaussian"}, \code{"binomial"},
#' \code{"poisson"}, or \code{"cox"}, and \code{penalty.factor} must not be
#' used. Fit additional lasso models by setting the hidden argument
#' \code{standard} to \code{TRUE}.
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
#' 
#' @references
#' A Rauschenberger, RX Menezes, MA Jonker, and MA van de Wiel (2018).
#' "Sparse regression with paired covariates."
#' \emph{Manuscript in preparation.} \email{a.rauschenberger@vumc.nl}
#' 
#' @examples
#' set.seed(1)
#' n <- 50; p <- 20
#' y <- rbinom(n=n,size=1,prob=0.5)
#' X <- lapply(1:2,function(x) matrix(rnorm(n*p),nrow=n,ncol=p))
#' object <- palasso(y=y,X=X,family="binomial",max=10,standard=TRUE)
#' 
palasso <- function(y,X,max=NULL,...){
    
    # extract
    base <- list(...)
    standard <- !(is.null(base$standard)||!base$standard)
    #adaptive <- !(is.null(base$adaptive)||!base$adaptive)
    base$standard <- base$adaptive <- NULL
    
    # checks
    funs <- list(glmnet::glmnet,glmnet::cv.glmnet)
    formals <- unlist(lapply(funs,function(x) formals(x)),recursive=FALSE)
    if(any(!names(base) %in% names(formals))){stop("Unexpected argument.",call.=FALSE)}
 
    # arguments
    base$y <- y
    base$x <- do.call(what="cbind",args=X)
    default <- list(family="gaussian",alpha=1,nfolds=10,type.measure="deviance")
    base <- c(base,default[!names(default) %in% names(base)])
    if(!base$family %in% c("gaussian","binomial","poisson","cox")){
        stop("Invalid argument \"family\".",call.=FALSE)
    }
    if(!is.null(base$dfmax) & base$alpha==0){
        stop("Unexpected argument \"dfmax\" as \"alpha=0\".",call.=FALSE)
    }
    
    # dimensionality
    if(length(unique(lapply(X,dim)))!=1){
        stop("Invalid argument \"X\".",call.=FALSE)
    }
    k <- ifelse(is.list(X),length(X),1)
    if(k==1){
        stop("Invalid argument \"X\".",call.=FALSE)
    }
    n <- nrow(X[[1]])
    p <- ncol(X[[1]])
    
    # penalty factor
    if(!is.null(base$penalty.factor)){
        stop("Unexpected argument \"penalty.factor\".",call.=FALSE)
    }
    
    # distribution
    if(survival::is.Surv(y)){
        guess <- "cox"
    } else {
        guess <- "gaussian"
        guess[all(base$y%%1==0 & base$y>=0)] <- "poisson"
        guess[!is.vector(base$y) | length(unique(base$y))==2] <- "binomial"
    }
    if(guess!=base$family){
        warning(paste0("Consider family \"",guess,"\"."),call.=FALSE)
    }
    
    # fold identifier
    #cond <- logical()
    c1 <- is.null(base$foldid)
    c2 <- base$family=="binomial" & is.vector(y)
    c3 <- base$family=="cox"
    if(c1 & (c2 | c3)){
        base$foldid <- palasso:::.folds(y=y,nfolds=base$nfolds) # changed!
    }
    
    # names
    #names <- names(X)
    #if(is.null(names)||anyDuplicated(names)){
        if(k==2){
            names <- c("x","z") 
        } else {
            names <- letters[seq_len(k)]
        }
    #}
    
    # Pearson correlation
    cor <- list()
    for(i in seq_len(k)){
        if(base$family=="cox"){
            cor[[i]] <- apply(X[[i]],2,function(x) abs(2*survival::survConcordance(y~x)$concordance-1))
        } else {
            cor[[i]] <- suppressWarnings(as.vector(abs(stats::cor(X[[i]],y))))
        }
        cor[[i]][is.na(cor[[i]])] <- 0
    }

    weight <- list()
    
    # standard lasso
    if(standard){ 
        temp <- list()
        for(i in seq_len(k)){
            temp[[i]] <- rep(1*(seq_len(k)==i),each=p)
        }
        temp[[k+1]] <- rep(1/k,times=k*p)
        names(temp) <- paste0("standard_",c(names,paste(names,collapse="")))
        weight <- c(weight,temp)
    }
    
    # adaptive lasso
    #if(adaptive){
        # ### via ridge regression ### (check!)
        # temp <- list()
        # for(i in seq_len(k)){
        #     args <- base
        #     args$alpha <- 0
        #     args$penalty.factor <- 1/rep(1*(seq_len(k)==i),each=p)
        #     model <- palasso:::.cv.glmnet(args)
        #     temp[[i]] <- abs(glmnet::coef.cv.glmnet(model,s="lambda.min")[-1])
        # }
        # names(temp) <- paste0("adaptive_",names)
        # args$penalty.factor <- rep(1,times=k*p)
        # model <- palasso:::.cv.glmnet(args)
        # temp[[k+1]] <- abs(glmnet::coef.cv.glmnet(model,s="lambda.min")[-1])
        # names(temp)[k+1] <- paste0("adaptive_",paste(names,collapse=""))
        # weight <- c(weight,temp)
        # ### via univariate regression ### (works!)
        # if(base$family=="cox"){
        #     mar <- abs(apply(base$x,2,function(x) survival::coxph(y~x)$coefficients))
        # } else {
        #     family <- eval(parse(text=base$family))()
        #     mar <- abs(apply(base$x,2,function(x) stats::glm.fit(y=y,x=cbind(1,x),family=family)$coefficients[2]))
        # }
        # mar[is.na(mar)] <- 0
        # temp <- list()
        # for(i in seq_len(k)){
        #     temp[[i]] <- rep(1*(seq_len(k)==i),each=p)*mar
        # }
        # temp[[k+1]] <- mar
        # names(temp) <- paste0("adaptive_",c(names,paste(names,collapse="")))
        # weight <- c(weight,temp)
    #}
    
    # adaptive lasso
    temp <- list()
    for(i in seq_len(k)){
        temp[[i]] <- rep(1*(seq_len(k)==i),each=p)*cor[[i]] 
    }
    temp[[k+1]] <- unlist(cor)
    names(temp) <- paste0("adaptive_",c(names,paste(names,collapse="")))
    weight <- c(weight,temp)
    
    # weighted lasso
    temp <- list() 
    temp[[1]] <- rep(1/k*vapply(cor,mean,numeric(1))/mean(unlist(cor)),each=p)
    temp[[2]] <- unlist(cor)/rowSums(do.call(cbind,cor))
    temp[[2]][is.na(temp[[2]])] <- 0
    names(temp) <- paste0(c("between_","within_"),paste(names,collapse=""))
    weight <- c(weight,temp)
    
    # cross-validation
    model <- list()
    args <- base
    for(i in seq_along(weight)){
        args$penalty.factor <- 1/weight[[i]]
        model[[i]] <- palasso:::.cv.glmnet(args)
        if(i > 1){ # free memory
            model[[i]]$glmnet.fit$call$x <- NULL 
        }
    }
    names(model) <- names(weight)
    
    # output
    call <- lapply(list(...),function(x) unlist(x))
    attributes(model)$info <- list(n=n,k=k,p=p,names=names,call=call,max=max)
    class(model) <- "palasso"
    return(model)
}


.folds <- function(y,nfolds){
    if(survival::is.Surv(y)){y <- y[,"status"]}
    foldid <- rep(NA,times=length(y))
    foldid[y==0] <- sample(rep(seq_len(nfolds),length.out=sum(y==0)))
    foldid[y==1] <- sample(rep(seq_len(nfolds),length.out=sum(y==1)))
    return(foldid)
}

.error <- function(x,args){
    pattern <- c("Error in predmat\\[which, seq\\(nlami\\)\\] <- preds",
                 "replacement has length zero")
    cond <- vapply(X=pattern,FUN=function(p) grepl(pattern=p,x=x),
                   FUN.VALUE=logical(1))
    if(all(cond)){
        warning("Fitting intercept-only model.",call.=FALSE)
        args$lambda <- c(99e99,99e98)
        do.call(what=glmnet::cv.glmnet,args=args)
    } else {
        stop(x,call.=FALSE)
    }
}

.warning <- function(x){
    #pattern <- c("from glmnet Fortran code \\(error code",
    #             "Number of nonzero coefficients along the path exceeds pmax=",
    #             "lambda value; solutions for larger lambdas returned")
    pattern <- c("from glmnet Fortran code \\(error code",
                 "Convergence for",
                 "lambda value not reached after maxit=",
                 "iterations; solutions for larger lambdas returned")
    cond <- vapply(X=pattern,FUN=function(p) grepl(pattern=p,x=x),
                   FUN.VALUE=logical(1))
    if(all(cond)){
        invokeRestart("muffleWarning")
    }
}

.cv.glmnet <- function(args){
    withCallingHandlers(expr=tryCatch(expr=do.call(what=glmnet::cv.glmnet,args=args),
                                      error=function(x) palasso:::.error(x,args)),
                        warning=function(x) palasso:::.warning(x))
}

#' Toydata
#' 
#' This dataset allows you to reproduce
#' the examples shown in the vignette.
#' 
#' @docType data
#' @keywords internal
#' @name toydata
#' @usage data(toydata)
#' @return All entries are numeric.
#' @format A list of numeric vectors and numeric matrices.
NULL

# weighted lasso
#if(FALSE){
### trial start ### (remove this?)
## multiply between and within weights
#extra <- list()
#extra[[1]] <- temp[[1]]*temp[[2]]
#names(extra) <- "combine_xz"
#weight <- c(weight,temp1)
### trial end ### (remove this?)
### trial start ###
## obtain weights from xz correlations
#cc <- sapply(seq_len(p),function(i) abs(cor(X[[1]][,i],X[[2]][,i])))
#cc[is.na(cc)] <- 0
#extra <- list()
#extra[[1]] <- rep(x=(2-cc)/2,times=2)
#names(extra) <- "combine_xz"
#weight <- c(weight,extra)
### trial end ###
#}

# # naive lasso
# if(FALSE){
#     # standard deviation
#     sd <- list()
#     for(i in seq_len(k)){
#         sd[[i]] <- apply(X[[i]],2,stats::sd)
#     }
#     # scaling
#     cm1 <- lapply(X,function(x) apply(x,2,mean))
#     cm2 <- lapply(X,function(x) apply(x,2,stats::var))
#     cond1 <- sapply(cm1,function(x) all(x>-0.01 & x<+0.01))
#     cond2 <- sapply(cm2,function(x) all(x>+0.99 & x<+1.01))
#     if(any(cond1)|any(cond2)){ # too strict?
#         stop("Provide unstandardised covariates!")
#     }
#     # weighting
#     temp <- list()
#     groups <- rep(1/k*sapply(cor,mean)/mean(unlist(cor)),each=p)
#     pairs <- unlist(sd)/rowSums(do.call(cbind,sd))
#     pairs[is.na(pairs)] <- 0
#     temp[[1]] <- groups*pairs
#     names(temp) <- "naive_xz"
#     weight <- c(weight,temp)
# }

## Original function
# palasso <- function(y,X,trial=FALSE,...){
#     
#     # checks
#     base <- list(...)
#     funs <- list(glmnet::glmnet,glmnet::cv.glmnet)
#     formals <- unlist(lapply(funs,function(x) formals(x)),recursive=FALSE) # changed!
#     if(any(!names(base) %in% names(formals))){stop("Invalid argument.")}
#     
#     # arguments
#     base$y <- y
#     base$x <- do.call(what="cbind",args=X)
#     default <- list(family="gaussian",alpha=1,nfolds=10,type.measure="deviance")
#     base <- c(base,default[!names(default) %in% names(base)])
#     if(!base$family %in% c("gaussian","binomial","poisson","cox")){
#         stop("Invalid argument \"family\".")
#     }
#     if(!is.null(base$pmax) & base$alpha==0){
#         stop("Unexpected argument \"pmax\" as \"alpha=0\".")
#     }
#     
#     # dimensionality
#     k <- ifelse(is.list(X),length(X),1)
#     if(k==1){
#         stop("Invalid argument \"X\".")
#     }
#     n <- nrow(X[[1]])
#     p <- ncol(X[[1]])
#     
#     # penalty factor (trial)
#     if(!is.null(base$penalty.factor)){
#         stop("Unexpected argument \"penalty.factor\".")
#     }
#     
#     # distribution (trial)
#     if(survival::is.Surv(y)){
#         guess <- "cox"
#     } else {
#         guess <- "gaussian"
#         guess[all(base$y%%1==0 & base$y>=0)] <- "poisson"
#         guess[!is.vector(base$y) | length(unique(base$y))==2] <- "binomial"
#     }
#     if(guess!=base$family){
#         warning(paste0("Consider family \"",guess,"\"."))
#     }
#     
#     # fold identifier
#     cond <- logical()
#     cond[1] <- is.null(base$foldid)
#     cond[2] <- base$family=="binomial"
#     cond[3] <- is.vector(y)
#     if(all(cond)){
#         base$foldid <- rep(NA,times=n)
#         base$foldid[y==0] <- sample(rep(seq_len(base$nfolds),
#                                         length.out=sum(y==0)))
#         base$foldid[y==1] <- sample(rep(seq_len(base$nfolds),
#                                         length.out=sum(y==1)))
#     }
#     
#     # marginal effects (correlation)
#     mar <- list()
#     for(i in seq_len(k)){
#         if(base$family=="cox"){
#             mar[[i]] <- apply(X[[i]],2,function(x) abs(2*survival::survConcordance(y~x)$concordance-1))
#         } else {
#             mar[[i]] <- suppressWarnings(as.vector(abs(stats::cor(X[[i]],y))))
#         }
#         mar[[i]][is.na(mar[[i]])] <- 0
#     }
#     
#     # # marginal effects (univariate regression)
#     # family <- eval(parse(text=base$family))()
#     # mar <- list()
#     # for(i in seq_len(k)){
#     #     mar[[i]] <- abs(apply(X[[i]],2,function(x) stats::glm.fit(y=y,x=cbind(1,x),family=family)$coefficients[2]))
#     #     mar[[i]][is.na(mar[[i]])] <- 0
#     # }
#     
#     # # trial: marginal effects (univariate regression, s.e. correction)
#     #mar <- list()
#     #for(i in seq_len(k)){
#     #    reg <- apply(X,2,function(x) summary(glm(y~x,family=family)))
#     #    beta <- sapply(reg,function(x) x$coefficients["x","Estimate"])
#     #    sd <- sapply(reg,function(x) x$coefficients["x","Std. Error"])
#     #    mar[[i]] <- abs(beta)/sd
#     #}
#     ## for object x of class "summary.glm":
#     ## x$coefficients[,"Estimate"]/x$coefficients[,"Std. Error"]
#     ## for object x of class "coxph":
#     ## x$coefficients/sqrt(x$var)
#     
#     weight <- model <- list()
#     
#     # standard lasso
#     for(i in seq_len(k)){
#         weight[[i]] <- rep(1*(seq_len(k)==i),each=p)
#     }
#     weight[[k+1]] <- rep(1/k,times=k*p)
#     
#     # weighted lasso
#     weight[[k+2]] <- rep(1/k*sapply(mar,mean)/mean(unlist(mar)),each=p)
#     weight[[k+3]] <- unlist(mar)/rowSums(do.call(cbind,mar))
#     weight[[k+3]][is.na(weight[[k+3]])] <- 0
#     
#     # adaptive lasso
#     for(i in seq_len(k)){
#         weight[[k+3+i]] <- weight[[i]]*mar[[i]]
#     }
#     weight[[2*k+4]] <- unlist(mar)
#     
#     if(trial){
#         c1 = apply(X[[1]],2,stats::sd)
#         c2 = apply(X[[2]],2,stats::sd)
#         temp = c(c1/(c1+c2),c2/c(c1+c2))*weight[[k+2]]
#         temp[is.na(temp)] = 0
#         weight[[2*k+5]] = temp
#     }
#     
#     # cross-validation
#     args <- base
#     for(i in seq_along(weight)){
#         args$penalty.factor <- 1/weight[[i]]
#         ### start original ###
#         #model[[i]] <- tryCatch(expr=do.call(what=glmnet::cv.glmnet,args=args),
#         #                       error=function(x) NA)
#         #if(class(model[[i]])!="cv.glmnet"){
#         #    # intercept-only model
#         #    temp <- base
#         #    temp$lambda <- c(99e99,99e98)
#         #    model[[i]] <- do.call(what=glmnet::cv.glmnet,args=temp)
#         #}
#         ### end original ###
#         model[[i]] <- palasso:::.cv.glmnet(args) ### trial
#         if(i > 1){
#             # free memory
#             model[[i]]$glmnet.fit$call$x <- NULL
#         }
#     }
#     
#     # names
#     if(k==2){
#         names = c("x","z") 
#     } else {
#         names = letters[seq_len(k)]
#     }  
#     all = paste0(names,collapse="")
#     if(trial){
#         names(model) = c(paste0("standard_",c(names,all)),
#                          paste0(c("between_","within_"),all),
#                          paste0("adaptive_",c(names,all)),"trial")
#     } else {
#         names(model) = c(paste0("standard_",c(names,all)),
#                          paste0(c("between_","within_"),all),
#                          paste0("adaptive_",c(names,all)))
#     }
#     
#     # output
#     call <- sapply(list(...),function(x) deparse(x))
#     attributes(model)$info <- list(n=n,k=k,p=p,names=names,call=call)
#     class(model) <- "palasso"
#     return(model)
# }
