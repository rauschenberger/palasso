
#--- Workhorse function --------------------------------------------------------

#' @title
#' Paired lasso
#' 
#' @export
#' @aliases palasso-package
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
#' positive numeric, or \code{NULL} (no sparsity constraint)
#' 
#' @param shrink
#' adaptive shrinkage of the initial coefficients\strong{:}
#' logical
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
#' used. Deactivate adaptive lasso by setting \code{adaptive} to \code{FALSE},
#' and activate standard lasso by setting \code{standard} to \code{TRUE}.
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
#' object <- palasso(y=y,X=X,family="binomial")
#' 
palasso <- function(y,X,max=10,shrink=TRUE,...){
    
    # extract
    base <- list(...)
    adaptive <- is.null(base$adaptive)||base$adaptive
    standard <- !(is.null(base$standard)||!base$standard)
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
    if(base$alpha==0 && !is.null(c(max,base$dfmax,base$pmax))){
        stop("Unexpected argument \"max\", \"dfmax\" or \"pmax\" (\"alpha=0\").",call.=FALSE)
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
    c1 <- is.null(base$foldid)
    c2 <- base$family=="binomial" & is.vector(y)
    c3 <- base$family=="cox"
    if(c1 & (c2 | c3)){
        base$foldid <- .folds(y=y,nfolds=base$nfolds)
    }
    
    # names
    if(k==2){
        names <- c("x","z") 
    } else {
        names <- letters[seq_len(k)]
    }
    
    ## Pearson correlation (original) # 18 June 2018
    cor <- list()
    for(i in seq_len(k)){
        if(base$family=="cox"){
            cor[[i]] <- apply(X[[i]],2,function(x) abs(2*survival::survConcordance(y~x)$concordance-1))
        } else {
            cor[[i]] <- suppressWarnings(as.vector(abs(stats::cor(X[[i]],y))))
        }
        cor[[i]][is.na(cor[[i]])] <- 0
    }
    
    ### trial start ### 18 June 2018
    if(shrink){
        for(i in 1:2){
            cor[[i]] <- CorShrink::CorShrinkVector(corvec=cor[[i]],nsamp_vec=n)
            if(all(cor[[i]]==0)){cor[[i]] <- rep(0.001,times=p)} # warning?
        }
    }
    ### trial end ###
    
    # # shrinkage (trial) #  18 June 2018
    #cor <- list()
    #for(i in seq_len(k)){
    #     cor[[i]] <- abs(.mar(y=y,X=X[[i]],family=base$family,shrink=shrink)$beta_eb)
    #}
    
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
        #
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
    if(adaptive){
        temp <- list()
        for(i in seq_len(k)){
            temp[[i]] <- rep(1*(seq_len(k)==i),each=p)*cor[[i]] 
        }
        temp[[k+1]] <- unlist(cor)
        names(temp) <- paste0("adaptive_",c(names,paste(names,collapse="")))
        weight <- c(weight,temp)
    }
    
    # weighted lasso # replaced on 18 June 2018
    #temp <- list() 
    #temp[[1]] <- rep(1/k*vapply(cor,mean,numeric(1))/mean(unlist(cor)),each=p)
    ##temp[[2]] <- unlist(cor)/rowSums(do.call(cbind,cor)) # standard (old)
    #temp[[2]] <- unlist(cor)^2/rowSums(do.call(cbind,cor)) # adaptive (new)
    #temp[[2]][is.na(temp[[2]])] <- 0
    #names(temp) <- paste0(c("between_","within_"),paste(names,collapse=""))
    #weight <- c(weight,temp)
    if(standard){
        temp <- list()
        #temp[[1]] <- unlist(cor)/rowSums(do.call(cbind,cor))
        temp[[1]] <- rep(1/k*vapply(cor,mean,numeric(1))/mean(unlist(cor)),each=p)
        names(temp) <- paste0("between_",paste(names,collapse=""))
        weight <- c(weight,temp)
    }
    if(adaptive){
        temp <- list()
        temp[[1]] <- unlist(cor)^2/rowSums(do.call(cbind,cor))
        names(temp) <- paste0("within_",paste(names,collapse=""))
        weight <- c(weight,temp)
    }
    
    # cross-validation
    model <- list()
    args <- base
    for(i in seq_along(weight)){
        args$penalty.factor <- 1/weight[[i]]
        model[[i]] <- .cv.glmnet(args)
        if(i > 1){ # free memory
            model[[i]]$glmnet.fit$call$x <- NULL 
        }
    }
    names(model) <- names(weight)
    
    # output
    call <- lapply(list(...),function(x) unlist(x))
    attributes(model)$info <- list(n=n,k=k,p=p,names=names,call=call,max=max,
                                   standard=standard,adaptive=adaptive)
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

# .error <- function(x,args){
#     pattern <- c("Error in predmat\\[which, seq\\(nlami\\)\\] <- preds",
#                  "replacement has length zero")
#     cond <- vapply(X=pattern,FUN=function(p) grepl(pattern=p,x=x),
#                    FUN.VALUE=logical(1))
#     if(all(cond)){
#         warning("Fitting intercept-only model.",call.=FALSE)
#         args$lambda <- c(99e99,99e98)
#         do.call(what=glmnet::cv.glmnet,args=args)
#     } else {
#         stop(x,call.=FALSE)
#     }
# }

# .warning <- function(x){
#     pattern <- c("from glmnet Fortran code \\(error code",
#                  "Convergence for",
#                  "lambda value not reached after maxit=",
#                  "iterations; solutions for larger lambdas returned")
#     cond <- vapply(X=pattern,FUN=function(p) grepl(pattern=p,x=x),
#                    FUN.VALUE=logical(1))
#     if(all(cond)){
#         invokeRestart("muffleWarning")
#     }
# }

.error <- function(x,args){
    pattern1 <- c("Error in predmat\\[which, seq\\(nlami\\)\\] <- preds",
                 "replacement has length zero")
    pattern2 <- c("Matrices must have same number of columns in",
                 "rbind2\\(.Call\\(dense_to_Csparse, x\\), y\\)")
    cond1 <- vapply(X=pattern1,FUN=function(p) grepl(pattern=p,x=x),
                   FUN.VALUE=logical(1))
    cond2 <- vapply(X=pattern2,FUN=function(p) grepl(pattern=p,x=x),
                    FUN.VALUE=logical(1))
    if(all(cond1)|all(cond2)){
        warning("Modified lambda sequence!",call.=FALSE)
        args$lambda <- exp(seq(from=log(99e99),to=log(99),length.out=100))
        initial <- do.call(what=glmnet::cv.glmnet,args=args)
        lambda.max <- min(initial$lambda[initial$nzero==0])
        if(is.null(args$lambda.min.ratio)){args$lambda.min.ratio <- 0.01}
        if(is.null(args$nlambda)){args$nlambda <- 100}
        sequence <- exp(seq(from=log(lambda.max),
                            to=log(lambda.max*args$lambda.min.ratio),
                            length.out=args$nlambda))
        do.call(what=glmnet::cv.glmnet,args=args) 
    } else {
        stop(x,call.=FALSE)
    }
}


.warning <- function(x){
    pattern1 <- c("from glmnet Fortran code \\(error code",
                 "Convergence for",
                 "lambda value not reached after maxit=",
                 "iterations; solutions for larger lambdas returned")
    pattern2 <- c("In getcoef\\(fit, nvars, nx, vnames\\) :",
                  "an empty model has been returned; probably a convergence issue")
    cond1 <- vapply(X=pattern1,FUN=function(p) grepl(pattern=p,x=x),
                   FUN.VALUE=logical(1))
    cond2 <- vapply(X=pattern2,FUN=function(p) grepl(pattern=p,x=x),
                    FUN.VALUE=logical(1))
    if(all(cond1)|all(cond2)){
        invokeRestart("muffleWarning")
    }
}

.cv.glmnet <- function(args){
    withCallingHandlers(expr=tryCatch(expr=do.call(what=glmnet::cv.glmnet,args=args),
                                      error=function(x) .error(x,args)),
                        warning=function(x) .warning(x))
}


## ashr shrinkage
# .mar <- function(y,X,family,shrink=FALSE){
#     
#     beta <- se <- rep(NA,times=ncol(X))
#     if(family=="cox"){
#         #cox <- apply(X,2,function(x) summary(survival::coxph(y~x))$coefficients)
#         #beta <- sapply(cox,function(x) x["x","coef"])
#         #se <- sapply(cox,function(x) x["x","se(coef)"])
#         for(i in seq_len(ncol(X))){
#             x <- X[,i]
#             cox <- summary(survival::coxph(y~x))$coefficients
#             if(any(dim(cox)!=c(1,5))){
#                 beta[i] <- se[i] <- NA
#             } else {
#                 beta[i] <- cox["x","coef"]
#                 se[i] <- cox["x","se(coef)"]
#             }
#         }
#     } else {
#         for(i in seq_len(ncol(X))){
#             x <- X[,i]
#             glm <- summary(stats::glm(y~x,family=family))$coefficients
#             if(any(dim(glm)!=c(2,4))){
#                 beta[i] <- se[i] <- NA
#             } else {
#                 beta[i] <- glm["x","Estimate"]
#                 se[i] <- glm["x","Std. Error"]
#             }
#         }
#     }
#     beta_hat <- beta/se
#     
#     if(shrink){
#         se <- rep(x=1,times=length(beta_hat))
#         shrinkage <- ashr::ash(betahat=beta_hat,sebetahat=se,
#                                mixcompdist="normal",pointmass=TRUE,
#                                optmethod="mixEM")
#         beta_eb <- shrinkage$result[,"PosteriorMean"]
#     }
#     
#     if(!shrink||all(beta_eb==0)){
#         beta_eb <- beta_hat
#         beta_eb[is.na(beta_eb)] <- 0
#     }
#     
#     #beta_eb <- abs(beta_eb)
# 
#     return(list(beta_hat=beta_hat,beta_eb=beta_eb))
# }

#list <- .mar(y=y,X=X[[1]],family="binomial")
#plot(x=list$beta_hat,y=list$beta_eb)


## own shrinkage
# .mar <- function(y,X,family){
#     beta <- se <- trim <- rep(NA,times=ncol(X))
#     if(family=="cox"){
#         cox <- abs(apply(X,2,function(x) summary(survival::coxph(y~x))))
#         beta <- sapply(cox,function(x) x$coefficients["x","coef"])
#         se <- sapply(cox,function(x) x$coefficients["x","se(coef)"])
#         beta[is.na(beta)] <- 0
#         se[is.na(se)] <- 0
#     } else {
#         for(i in seq_len(ncol(X))){
#             x <- X[,i]
#             glm <- stats::glm(y~x,family=family)
#             trim[i] = any(glm$fitted.values>1-1e-14 | glm$fitted.values<1e-14)
#             temp = summary(glm)$coefficients
#             if(nrow(temp)==1){
#                 beta[i] = se[i] = 0
#             } else {
#                 beta[i] <- temp["x","Estimate"]
#                 se[i] <- temp["x","Std. Error"]
#             }
#         }
#     }
# 
#     # trimming
#     #if(family=="binomial"){
#     #    #trim <- sapply(glm,function(x) any(x$fitted.values>1-1e-14 | x$fitted.values<1e-14))
#     #    cutoff <- max(abs(beta[!trim]))
#     #    beta[beta < -cutoff] <- -cutoff
#     #    beta[beta > cutoff] <- cutoff
#     #    cutoff <- max(se[!trim])
#     #    se[se > cutoff] <- cutoff
#     #}
# 
#     # shrinkage
#     tausq <- pmax(0,var(beta)-mean(se^2)) # originally "mean"
#     weight <- (tausq)/(tausq+se^2)
#     eb <- mean(beta) + weight*(beta-mean(beta))
#     # eb <- weight*beta # trial
# 
#     # absolute value
#     cor <- as.numeric(cor(y,X))
#     
#     return(list(beta=beta,se=se,eb=eb,trim=trim,cor=cor))
# }


#' Toydata
#' 
#' Dataset for reproducing the vignette.
#' 
#' @docType data
#' @keywords internal
#' @name toydata
#' @usage data(toydata)
#' @return All entries are numeric.
#' @format A list of numeric vectors and numeric matrices.
NULL
