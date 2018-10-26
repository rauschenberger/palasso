
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
#' activate standard lasso by setting \code{standard} to \code{TRUE},
#' and deactivate shrinkage by setting \code{shrink} to \code{FALSE}.
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
#' A Rauschenberger, I Ciocanea-Teodorescu, RX Menezes, MA Jonker,
#' and MA van de Wiel (2018). "Sparse regression with paired covariates."
#' \emph{Manuscript in preparation.} \email{a.rauschenberger@vumc.nl}
#' 
#' @examples
#' set.seed(1)
#' n <- 50; p <- 20
#' y <- rbinom(n=n,size=1,prob=0.5)
#' X <- lapply(1:2,function(x) matrix(rnorm(n*p),nrow=n,ncol=p))
#' object <- palasso(y=y,X=X,family="binomial")
#' 
palasso <- function(y=y,X=X,max=10,...){
    
    ### start function ###
    ## INCLUDE: test whether y and family are compatible
    ## INCLUDE: test whether any arguments are redundant
    args <- args.trial(...)
    args <- c(args,dims.trial(y=y,X=X,args=args))
    x <- do.call(what="cbind",args=X) # fuse covariates
    ### end function ###

    # fold identifiers
    foldid <- folds.trial(y=y,nfolds=args$nfolds,foldid=args$foldid)
    
    # model fitting
    fit.full <- fit.trial(y=y,X=x,args=args)
    
    # lambda sequence
    if(is.null(args$lambda)){
        lambda <- lapply(fit.full,function(x) x$lambda)
    } else {
        lambda <- lapply(seq_len(args$num),function(x) args$lambda)
    }
    
    # internal cross-validation
    fit <- cv.trial(y=y,x=x,foldid=foldid,lambda=lambda,args=args)
    
    # loss sequence
    cvm <- loss.trial(y=y,fit=fit,family=args$family,type.measure=args$type.measure,foldid=foldid)
    
    # optimisation
    model <- extract.trial(fit=fit.full,lambda=lambda,cvm=cvm,type.measure=args$type.measure)
    
    # output
    call <- lapply(list(...),function(x) unlist(x))
    attributes(model)$info <- list(n=args$n,k=args$k,p=args$p,
                                   names=args$names,call=call,max=max,
                                   standard=args$standard,adaptive=args$adaptive)
    class(model) <- "palasso"
    return(model)
}

#' @title to do
#' @description to do
#' @param family to do
#' @param n sample size
#' @param p number of covariates
#' @return to do
#' @examples
#' NA
sim.trial <- function(family,n=200,p=150){
    X <- matrix(data=stats::rnorm(n=n*p),nrow=n,ncol=p)
    Z <- matrix(data=stats::rnorm(n=n*p),nrow=n,ncol=p)
    if(family=="gaussian"){
        y <- stats::rnorm(n=n)
    } else if(family=="binomial"){
        y <- stats::rbinom(n=n,size=1,prob=0.2)
    } else if(family=="poisson"){
        y <- stats::rpois(n=n,lambda=4)
    } else if(family=="cox"){
        event <- stats::rbinom(n=n,size=1,prob=0.3)
        time <- stats::rpois(n=n,lambda=4)+1
        y <- survival::Surv(time=time,event=event)
    } else {
        stop("Invalid family.")
    }
    return(list(y=y,X=X,Z=Z))
}

#' @title
#' dimensionality
#' 
#' @description
#' this function ...
#' 
#' @param y
#' to do
#' @param X
#' to do
#' @param args
#' to do
#' 
#' @return
#' to do
#' 
#' @examples
#' NA
dims.trial <- function(y,X,args=NULL){
    
    # dimensionality
    if(length(unique(lapply(X=X,FUN=dim)))!=1){
        stop("Invalid argument \"X\".",call.=FALSE)
    }
    k <- ifelse(is.list(X),length(X),1)
    if(k==1){
        stop("Invalid argument \"X\".",call.=FALSE)
    }
    
    # dimensionality
    n <- nrow(X[[1]])
    p <- ncol(X[[1]])
    k <- length(X)
    
    # names
    if(k==2){
        names <- c("x","z") 
    } else {
        names <- letters[seq_len(k)]
    }
    
    # number of models
    adaptive <- is.null(args$adaptive)||args$adaptive
    standard <- !(is.null(args$standard)||!args$standard)
    num <- (standard+adaptive)*(k+2)
    
    list <- list(n=n,p=p,k=k,num=num,names=names)
}

#' @title to do
#' @description to do
#' @param ... to do
#' @return to do
#' @examples
#' NA
args.trial <- function(...){
    
    args <- list(...)
    
    # model
    if(is.null(args$family)){args$family <- "gaussian"}
    if(is.null(args$alpha)){args$alpha <- 1}
    if(is.null(args$type.measure)){args$type.measure <- "deviance"}
    if(is.null(args$grouped)){args$grouped <- TRUE}
    
    # lambda sequence
    if(is.null(args$lambda)){
        if(is.null(args$nlambda)){args$nlambda <- 100}
    } else {
        args$nlambda <- length(args$lambda)
    }
    
    # fold identifier
    if(is.null(args$foldid)){
        if(is.null(args$nfolds)){args$nfolds <- 10}
    } else {
        args$nfolds <- length(unique(args$foldid))
    }
    
    # weighting schemes
    if(is.null(args$shrink)){args$shrink <- TRUE}
    if(is.null(args$standard)){args$standard <- FALSE}
    if(is.null(args$adaptive)){args$adaptive <- TRUE}
    
    # # model bag
    # adaptive <- is.null(base$adaptive)||base$adaptive
    # standard <- !(is.null(base$standard)||!base$standard)
    # 
    # # default cv.glmnet
    # names0 <- names(formals(glmnet::cv.glmnet))
    # args0 <- base[names(base) %in% names0]
    # def0 <- list(type.measure="deviance",nfolds=10)
    # base0 <- c(args0,def0[!names(def0) %in% names(base)])
    # 
    # # default glmnet
    # names1 <- names(formals(glmnet::glmnet))
    # args1 <- base[names(base) %in% names1]
    # def1 <- list(family="gaussian",alpha=1,nlambda=100)
    # base1 <- c(args1,def1[!names(def1) %in% names(base)])
    # 
    # base0$adaptive <- adaptive
    # base0$standard <- standard
    # base0$glmnet <- base1
    
    # check arguments
    #if(any(!names(args) %in% c(names0,names1,"adaptive","standard"))){
    #    stop("Unexpected argument.",call.=FALSE)
    #}
    if(!args$family %in% c("gaussian","binomial","poisson","cox")){
        stop("Invalid argument \"family\".",call.=FALSE)
    }
    #if(args$alpha==0 && !is.null(c(max,args$dfmax,args$pmax))){
    #    # missing argument "max"!!!
    #    stop("Unexpected argument \"max\", \"dfmax\" or \"pmax\" (\"alpha=0\").",call.=FALSE)
    #}
    if(!is.null(args$penalty.factor)){
        stop("Unexpected argument \"penalty.factor\".",call.=FALSE)
    }
    
    # distribution
    #if(survival::is.Surv(y)){
    #    guess <- "cox"
    #} else {
    #    guess <- "gaussian"
    #    guess[all(base$y%%1==0 & base$y>=0)] <- "poisson"
    #    guess[!is.vector(base$y) | length(unique(base$y))==2] <- "binomial"
    #}
    #if(guess!=base$family){
    #    warning(paste0("Consider family \"",guess,"\"."),call.=FALSE)
    #}
    
    return(args)
}

#' @title to do
#' @description to do
#' @param y response
#' @param nfolds ...
#' @param foldid ...
#' @return to do
#' @examples
#' NA
folds.trial <- function(y,nfolds,foldid){
    if(!is.null(foldid)){return(foldid)}
    if(survival::is.Surv(y)){y <- y[,"status"]}
    if(all(y %in% c(0,1))){
        foldid <- rep(x=NA,times=length(y))
        foldid[y==0] <- sample(x=rep(x=seq_len(nfolds),length.out=sum(y==0)))
        foldid[y==1] <- sample(x=rep(x=seq_len(nfolds),length.out=sum(y==1)))
    } else {
        foldid <- sample(x=rep(x=seq_len(nfolds),length.out=length(y)))
    }
    return(foldid)
}

#' @title fit models
#' @description this function ...
#' @param y response
#' @param X covariates
#' @param args see ...
#' @return list of glmnet-like objects
#' @examples
#' NA
fit.trial <- function(y,X,args){
    cor <- cor.trial(y=y,x=X,args=args)
    weight <- weight.trial(cor=cor,args=args)
    
    net <- list()
    for(j in seq_along(weight)){
        args <- args[c("alpha","family","nlambda")]
        args$y <- y
        args$x <- X
        args$lambda <- NULL # important!
        args$penalty.factor <- 1/weight[[j]]
        net[[j]] <- do.call(what=glmnet::glmnet,args=args)
        
        iter <- 0
        while((min(net[[j]]$df)>3)|(length(net[[j]]$lambda)==1)){
            warning("Modifying lambda sequence.",call.=FALSE)
            iter <- iter + 1
            if(iter>10){
                #browser()
                stop("Modifying lambda sequence failed.",call.=FALSE)
            }
            args$lambda <- exp(seq(from=log(99e99),to=log(0.01),length.out=100))
            initial <- do.call(what=glmnet::glmnet,args=args)
            lambda.max <- min(initial$lambda[initial$df==0])
            args$lambda <- exp(seq(
                from=log(lambda.max),
                to=log(0.01*lambda.max),
                length.out=pmax(args$nlambda,100)))
            net[[j]] <- do.call(what=glmnet::glmnet,args=args)
            
        }
        
        # #lambda[[j]] <- net[[j]]$lambda
        # 
        # # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - 
        # # this whole chunk has become redundant;
        # # call glmnet::glmnet again if there is a problem,
        # # maybe with a modified lambda sequence
        # 
        # if(flexible){
        #     if(min(net[[j]]$df)>3 | length(net[[j]]$lambda)==1){
        #         message("Modified all lambdas!")
        #         args$lambda <- exp(seq(from=log(99e99),to=log(0.01),length.out=100))
        #         initial <- do.call(what=glmnet::glmnet,args=args)
        #         lambda.max <- min(initial$lambda[initial$df==0])
        #         args$lambda <- sequence <- exp(seq(from=log(lambda.max),
        #                                 to=log(0.01*lambda.max),
        #                                 length.out=args$nlambda))
        #         #test <- lambda.max*exp((log(0.01*lambda.max)-log(lambda.max))/100)^seq_len(100)
        #         net[[j]] <- do.call(what=glmnet::glmnet,args=args)
        #     }
        # }
        # 
        # ### start trial ###
        # while(length(net[[j]]$lambda)==1){
        #     message("Modify first lambda!")
        #     args$lambda[1] <- 1.5*args$lambda[1]
        #     net[[j]] <- do.call(what=glmnet::glmnet,args=args)
        # }
        # ### end trial ###
        # 
        # #wait <- TRUE # trial
        # #while(wait){ # trial
        # #    net[[j]] <- do.call(what=glmnet::glmnet,args=args)
        # #    wait <- net[[j]]$lambda[1]==Inf # trial
        # #    if(wait){ # trial
        # #        #args$lambda[1] <- max(c(2,2*lambda[[j]][1])) # trial c(2,...)
        # #        args$lambda[1] <- 1.1*args$lambda[1] # trial
        # #        message("Modified first lambda!")
        # #    } # trial
        # #}
        # # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - 
        # if(length(net[[j]]$lambda)==1){browser(); stop("lambda length one")}
        
        if(j > 1){ # free memory
            net[[j]]$call$x <- NULL 
        }
    }
    names(net) <- names(weight)
    return(net)
}

#' @title correlation
#' @description this function ...
#' @param y response
#' @param x covariates
#' @param args see ...
#' @return list of vectors
#' @examples
#' NA
cor.trial <- function(y,x,args){
    if(args$family=="cox"){
        cor <- apply(X=x,MARGIN=2,FUN=function(x) abs(2*survival::survConcordance(y~x)$concordance-1))
    } else {
        cor <- suppressWarnings(as.vector(abs(stats::cor(x,y))))
    }
    cor[is.na(cor)] <- 0
    cor <- lapply(seq_len(args$k),function(x) cor[((x-1)*args$p+1):(x*args$p)]) # check this
    for(i in seq_len(args$k)){
        if(args$shrink){
            cor[[i]] <- CorShrink::CorShrinkVector(corvec=cor[[i]],nsamp_vec=nrow(x))
            if(all(cor[[i]]==0)){cor[[i]] <- rep(0.001,times=args$p)} # warning?
        }
    }
    return(cor)
}

#' @title compute weights
#' @description this function ...
#' @param cor see ...
#' @param args see ...
#' @return list of vectors
#' @examples
#' NA
weight.trial <- function(cor,args){
    k <- length(cor)
    p <- length(cor[[1]])
    weight <- list()
    if(args$standard){ 
        temp <- list()
        for(i in seq_len(k)){
            temp[[i]] <- rep(1*(seq_len(k)==i),each=p)
        }
        temp[[k+1]] <- rep(1/k,times=k*p)
        names(temp) <- paste0("standard_",c(args$names,paste(args$names,collapse="")))
        weight <- c(weight,temp)
    }
    if(args$adaptive){
        temp <- list()
        for(i in seq_len(k)){
            temp[[i]] <- rep(1*(seq_len(k)==i),each=p)*cor[[i]] 
        }
        temp[[k+1]] <- unlist(cor)
        names(temp) <- paste0("adaptive_",c(args$names,paste(args$names,collapse="")))
        weight <- c(weight,temp)
    }
    if(args$standard){
        temp <- list()
        temp[[1]] <- unlist(cor)/rowSums(do.call(cbind,cor))
        temp[[1]][is.na(temp[[1]])] <- 0
        names(temp) <- paste0("between_",paste(args$names,collapse=""))
        weight <- c(weight,temp)
    }
    if(args$adaptive){
        temp <- list()
        temp[[1]] <- unlist(cor)^2/rowSums(do.call(cbind,cor))
        temp[[1]][is.na(temp[[1]])] <- 0
        names(temp) <- paste0("within_",paste(args$names,collapse=""))
        weight <- c(weight,temp)
    }
    return(weight)
}

#' @title cross-validation
#' @description this function
#' @param y response
#' @param x covariates
#' @param foldid vector
#' @param lambda vector
#' @param args see ...
#' @return to do
#' @examples
#' NA
cv.trial <- function(y,x,foldid,lambda,args){
    fit <- lapply(seq_len(args$num),function(i) matrix(nrow=ifelse(args$family=="cox",args$nfolds,args$n),ncol=length(lambda[[i]])))
    for(i in seq_len(args$nfolds)){
        y0 <- y[foldid!=i]
        y1 <- y[foldid==i]
        X0 <- x[foldid!=i,,drop=FALSE]
        X1 <- x[foldid==i,,drop=FALSE]
        fit.sub <- fit.trial(y=y0,X=X0,args=args)
        for(j in seq_len(args$num)){
            if(args$family=="cox"){
                beta <- glmnet::predict.coxnet(object=fit.sub[[j]],newx=X1,type="coeff",s=lambda[[j]])
                # check whether these are the betas !!!
                if(args$grouped){
                    plfull <- glmnet::coxnet.deviance(x=x,y=y,beta=beta)
                    plminusk <- glmnet::coxnet.deviance(x=X0,y=y0,beta=beta)
                    temp <- plfull - plminusk
                } else {
                    temp <- glmnet::coxnet.deviance(x=X1,y=y1,beta=beta)
                }
                fit[[j]][i,seq_along(temp)] <- temp
            } else {
                #temp <- glmnet::predict.glmnet(object=fit.sub[[j]],newx=X1,type="response",s=lambda[[j]])
                temp <- stats::predict(object=fit.sub[[j]],newx=X1,type="response",s=lambda[[j]])
                # check whether 0 < temp < 1 in binomial
                if(any(is.na(temp))|(ncol(temp)==1)){
                    stop("Prediction problem.",call.=FALSE)
                }
                fit[[j]][foldid==i,seq_len(ncol(temp))] <- temp
            }
        }
    }
    return(fit) 
}

#' @title to do
#' @description to do
#' @param y response
#' @param fit see ..
#' @param family ...
#' @param type.measure ...
#' @param foldid to do
#' @return to do
#' @examples
#' NA
loss.trial <- function(y,fit,family,type.measure,foldid=NULL){
    if(!is.list(fit)){fit <- list(fit)}
    loss <- list()
    for(i in seq_along(fit)){
        if(is.vector(fit[[i]])){fit[[i]] <- as.matrix(fit[[i]])}
        if(is.null(foldid)&(family=="cox"|type.measure=="auc")){
            stop("Missing foldid.",call.=FALSE)
        }
        if(family=="gaussian"){
            if(type.measure %in% c("deviance","mse")){
                loss[[i]] <- apply(X=fit[[i]],MARGIN=2,FUN=function(x) mean((x-y)^2))
            } else if(type.measure=="mae"){
                loss[[i]] <- apply(X=fit[[i]],MARGIN=2,FUN=function(x) mean(abs(x-y)))
            } else {
                stop("Invalid type measure.",call.=FALSE)
            }
        } else if(family=="binomial"){
            if(type.measure=="deviance"){
                # EXAMINE BOUNDARIES
                limit <- 1e-05
                fit[[i]][fit[[i]]<limit] <- limit
                fit[[i]][fit[[i]]>1-limit] <- 1-limit
                loss[[i]] <- apply(X=fit[[i]],MARGIN=2,FUN=function(x) mean(-2*(y*log(x)+(1-y)*log(1-x))))
            } else if(type.measure=="mse"){
                loss[[i]] <- apply(X=fit[[i]],MARGIN=2,FUN=function(x) 2*mean((x-y)^2))
            } else if (type.measure=="mae"){
                loss[[i]] <- apply(X=fit[[i]],MARGIN=2,FUN=function(x) 2*mean(abs(x-y)))
            } else if(type.measure=="class"){
                loss[[i]] <- apply(X=fit[[i]],MARGIN=2,FUN=function(x) mean(abs(round(x)-y)))
            } else if(type.measure=="auc"){
                weights <- table(foldid)
                cvraw <- matrix(data=NA,nrow=length(weights),ncol=ncol(fit))
                for(k in seq_along(weights)){
                    cvraw[k,] <- apply(X=fit[[i]],MARGIN=2,FUN=function(x) glmnet::auc(y=y[foldid==k],prob=x[foldid==k]))
                }
                loss[[i]] <- apply(X=cvraw,MARGIN=2,FUN=function(x) stats::weighted.mean(x=x,w=weights,na.rm=TRUE))
            } else {
                stop("Invalid type measure.",call.=FALSE)
            }
        } else if(family=="poisson"){
            if(type.measure=="deviance"){
                loss[[i]] <- apply(X=fit[[i]],MARGIN=2,FUN=function(x) mean(2*(ifelse(y==0,0,y*log(y/x))-y+x),na.rm=TRUE))
            } else if(type.measure=="mse"){
                loss[[i]] <- apply(X=fit[[i]],MARGIN=2,FUN=function(x) mean((x-y)^2))
            } else if(type.measure=="mae"){
                loss[[i]] <- apply(X=fit[[i]],MARGIN=2,FUN=function(x) mean(abs(x-y)))
            } else {
                stop("Invalid type measure.",call.=FALSE)
            }
        } else if(family=="cox"){
            if(type.measure=="deviance"){
                weights <- tapply(X=y[,"status"],INDEX=foldid,FUN=sum)
                loss[[i]] <- apply(X=fit[[i]],MARGIN=2,FUN=function(x) stats::weighted.mean(x=x/weights,w=weights,na.rm=TRUE))
            } else {
                stop("Invalid type measure.",call.=FALSE)
            }
        } else {
            stop("Invalid family.",call.=FALSE)
        }
        if(sum(diff(is.na(loss[[i]])))==1){ # trial
            loss[[i]] <- loss[[i]][!is.na(loss[[i]])] # trial
        } # trial
    }
    return(loss)
}

#' @title extract
#' @description this function ...
#' @param fit ...
#' @param lambda ...
#' @param cvm ...
#' @param type.measure ...
#' @details to do
#' @return this ...
#' @examples
#' NA
extract.trial <- function(fit,lambda,cvm,type.measure){
    model <- sapply(X=fit,FUN=function(x) list())
    for(i in seq_along(fit)){
        if(type.measure=="auc"){
            sel <- which.max(cvm[[i]])
        } else {
            sel <- which.min(cvm[[i]])
        }
        model[[i]]$lambda <- lambda[[i]]
        model[[i]]$cvm <- cvm[[i]]
        model[[i]]$name <- type.measure # strange
        model[[i]]$glmnet.fit <- fit[[i]]
        model[[i]]$nzero <- Matrix::colSums(fit[[i]]$beta!=0)
        model[[i]]$lambda.min <- lambda[[i]][sel]
        #model[[i]]$lambda.1se <-  # maybe solve this ...
    }
    return(model)
}

.folds <- function(y,nfolds){
  if(survival::is.Surv(y)){y <- y[,"status"]}
  foldid <- rep(NA,times=length(y))
  foldid[y==0] <- sample(rep(seq_len(nfolds),length.out=sum(y==0)))
  foldid[y==1] <- sample(rep(seq_len(nfolds),length.out=sum(y==1)))
  return(foldid)
}

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

