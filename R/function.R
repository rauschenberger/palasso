
#' @export
#' @title
#' Paired lasso
#' 
#' @inheritParams arguments
#' 
#' @aliases palasso-package
#' 
#' @description
#' The function \code{palasso} fits the paired lasso.
#' Use this function if you have \emph{paired covariates}
#' and want a \emph{sparse model}.
#' 
#' @details
#' 
#' Let \code{x} denote one entry of the list \code{X}. See \link[glmnet]{glmnet}
#' for alternative specifications of \code{y} and \code{x}. Among the further
#' arguments, \code{family} must equal \code{"gaussian"}, \code{"binomial"},
#' \code{"poisson"}, or \code{"cox"}, and \code{penalty.factor} must not be
#' used.
#' 
#' Hidden arguments:
#' Deactivate adaptive lasso by setting \code{adaptive} to \code{FALSE},
#' activate standard lasso by setting \code{standard} to \code{TRUE},
#' and activate shrinkage by setting \code{shrink} to \code{TRUE}.
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
#' and MA van de Wiel (2020). "Sparse classification with paired covariates."
#' \emph{Advances in Data Analysis and Classification}. 14:571-588.
#' \doi{10.1007/s11634-019-00375-6},
#' \href{https://link.springer.com/content/pdf/10.1007/s11634-019-00375-6.pdf}{pdf},
#' \email{armin.rauschenberger@uni.lu}
#' 
#' @examples
#' set.seed(1)
#' n <- 50; p <- 20
#' y <- rbinom(n=n,size=1,prob=0.5)
#' X <- lapply(1:2,function(x) matrix(rnorm(n*p),nrow=n,ncol=p))
#' object <- palasso(y=y,X=X,family="binomial") # adaptive=TRUE,standard=FALSE
#' names(object)
#' 
palasso <- function(y=y,X=X,max=10,...){
    
    # arguments
    args <- .args(...)
    args <- c(args,.dims(y=y,X=X,args=args))
    x <- do.call(what="cbind",args=X) # fuse covariates

    # fold identifiers
    foldid <- .folds(y=y,nfolds=args$nfolds,foldid=args$foldid)
    
    # model fitting
    fit.full <- .fit(y=y,x=x,args=args)
    
    # lambda sequence
    if(is.null(args$lambda)){
        lambda <- lapply(fit.full,function(x) x$lambda)
    } else {
        lambda <- lapply(seq_len(args$num),function(x) args$lambda)
    }
    
    # internal cross-validation
    fit <- .cv(y=y,x=x,foldid=foldid,lambda=lambda,args=args)
    
    # loss sequence
    cvm <- .loss(y=y,fit=fit,family=args$family,type.measure=args$type.measure,foldid=foldid)
    
    # optimisation
    model <- .extract(fit=fit.full,lambda=lambda,cvm=cvm,type.measure=args$type.measure)
    
    # output
    call <- lapply(list(...),function(x) unlist(x))
    attributes(model)$info <- list(n=args$n,k=args$k,p=args$p,
                                   names=args$names,call=call,max=max,
                                   standard=args$standard,adaptive=args$adaptive)
    class(model) <- "palasso"
    return(model)
}


#' @title
#' Arguments for "palasso"
#' 
#' @name arguments
#' 
#' @description
#' This page lists the arguments for the (internal) "palasso" function(s).
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
#' @param x
#' covariates\strong{:}
#' matrix with \eqn{n} rows (samples)
#' and \eqn{k * p} columns (variables)
#'
#' @param args
#' options for paired lasso\strong{:}
#' list of arguments
#' (output from \link{.dims} and \link{.args})
#' 
#' @param nfolds
#' number of folds\strong{:}
#' positive integer
#' (\eqn{>= 10} recommended)
#' 
#' @param foldid
#' fold identifiers\strong{:}
#' vector of length \eqn{n},
#' with entries from \eqn{1} to \code{nfolds}
#' 
#' @param cor
#' correlation coefficients\strong{:}
#' list of \eqn{k} vectors of length \eqn{p}
#' (one vector for each covariate set with
#' one entry for each covariate)
#' 
#' @param lambda
#' lambda sequence\strong{:}
#' vector of decreasing positive values
#' 
#' @param family
#' model family\strong{:}
#' character "gaussian", "binomial", "poisson", or "cox"
#' 
#' @param type.measure ...
#' loss function\strong{:}
#' character "deviance", "mse", "mae", "class", or "auc"
#' 
#' @param fit
#' matrix with one row for each sample
#' ("gaussian", "binomial" and "poisson"),
#' or one row for each fold (only "cox"),
#' and one column for each \code{lambda}
#' (output from \link{.fit})
#'
#' @param cvm
#' mean cross-validated loss\strong{:}
#' vector of same length as \code{lambda}
#' (output from \link{.loss})
#' 
NULL

#' @title Arguments
#' 
#' @description
#' Checks the validity of the provided arguments.
#' 
#' @param ...
#' Arguments supplied to \code{\link{palasso}},
#' other than \code{y}, \code{X} and \code{max}.
#'  
#' @return
#' Returns the arguments as a list, including default values
#' for missing arguments.
#' 
#' @examples
#' NA
#' 
.args <- function(...){
    
    args <- list(...)
    
    # check arguments (trial)
    names0 <- names(formals(glmnet::cv.glmnet))
    names1 <- names(formals(glmnet::glmnet))
    if(any(!names(args) %in% c(names0,names1,"adaptive","standard","shrink","elastic"))){
      stop("Unexpected argument.",call.=FALSE)
    }
    
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
    if(is.null(args$shrink)){args$shrink <- FALSE} # was TRUE
    
    # conditional usage of CorShrink
    if(args$shrink & (!"CorShrink" %in% .packages(all.available=TRUE))){
      warning("Install \"CorShrink\", or shrink=FALSE.",call.=FALSE)
      args$shrink <- FALSE
    }
    
    if(is.null(args$standard)){args$standard <- FALSE}
    if(is.null(args$adaptive)){args$adaptive <- TRUE}
    if(is.null(args$elastic)){args$elastic <- FALSE}
    
    # check consistency
    if(!args$family %in% c("gaussian","binomial","poisson","cox")){
        stop("Invalid argument \"family\".",call.=FALSE)
    }
    #if(args$alpha==0 && !is.null(c(max,args$dfmax,args$pmax))){
    #    # argument "max" here unavailable
    #    stop("Unexpected argument \"max\", \"dfmax\" or \"pmax\" (\"alpha=0\").",call.=FALSE)
    #}
    if(!is.null(args$penalty.factor)){
        stop("Unexpected argument \"penalty.factor\".",call.=FALSE)
    }
    
    return(args)
}

#' @title
#' Dimensionality
#' 
#' @description
#' This function extracts the dimensions.
#' 
#' @inheritParams arguments
#' 
#' @return
#' The function \code{.dims} extracts the dimensionality.
#' It returns the numbers of samples,
#' covariate pairs and covariate sets.
#' It also returns the number of weighting schemes,
#' and the names of these weighting schemes.
#' 
#' @examples 
#' NA
#' 
.dims <- function(y,X,args=NULL){
  
  # distribution (trial)
  if(survival::is.Surv(y)){
      guess <- "cox"
  } else {
      guess <- "gaussian"
      guess[all(y%%1==0 & y>=0)] <- "poisson"
      guess[!is.vector(y) | length(unique(y))==2] <- "binomial"
  }
  if(guess!=args$family){
      warning(paste0("Consider family \"",guess,"\"."),call.=FALSE)
  }
  
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
  elastic <- !(is.null(args$elastic)||!args$elastic)
  num <- (standard+adaptive)*(k+2) + 1*elastic # was 4*elastic
  
  list <- list(n=n,p=p,k=k,num=num,names=names)
}

#' @title Cross-validation folds
#' 
#' @description
#' Assigns samples to cross-validation folds,
#' balancing the folds in the case of a binary or survival response.
#' 
#' @inheritParams arguments
#' 
#' @return
#' Returns the fold identifiers.
#' 
#' @examples
#' NA
#' 
.folds <- function(y,nfolds,foldid=NULL){
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

#' @title Model bag
#' 
#' @description Fits all models from the chosen bag.
#' 
#' @inheritParams arguments
#' 
#' @return list of \code{glmnet}-like objects
#' 
#' @examples
#' NA
.fit <- function(y,x,args){
    cor <- .cor(y=y,x=x,args=args)
    weight <- .weight(cor=cor,args=args)
    
    net <- list()
    for(j in seq_along(weight)){
        net[[j]] <- .fit.int(y=y,x=x,weight=weight[[j]],args=args)
        if(j > 1){ # free memory
            net[[j]]$call$x <- NULL 
        }
    }
    names(net) <- names(weight)
  
    if(args$elastic){
      alpha <- 0.95 # was c(0.25,0.50,0.75,1.00)
      #args$lambda.min.ratio <- 0.1 # not good for unconstrained version
      args$nlambda <- 2*args$nlambda # higher resolution
      for(k in seq_along(alpha)){
        args$alpha <- alpha[k]
        net[[j+k]] <- .fit.int(y=y,x=x,weight=weight$standard_xz,args=args)
        names(net)[j+k] <- paste0("elastic",100*args$alpha)
      }
      rm(args)
    }
    
    return(net)
}

.fit.int <- function(y,x,weight,args){
  args <- args[c("alpha","family","nlambda")]
  args$y <- y
  args$x <- x
  args$lambda <- NULL # important!
  args$penalty.factor <- 1/weight
  net <- do.call(what=glmnet::glmnet,args=args)

  ### start original ###
  iter <- 0
  while((min(net$df)>3)|(length(net$lambda)==1)){
   warning("Modifying lambda sequence.",call.=FALSE)
   iter <- iter + 1
   if(iter>10){
     stop("Modifying lambda sequence failed.",call.=FALSE)
   }
   args$lambda <- exp(seq(from=log(99e99),to=log(0.01),length.out=100))
   initial <- do.call(what=glmnet::glmnet,args=args)
   if(all(is.na(initial$lambda[initial$df==0]))){next}
   lambda.max <- min(initial$lambda[initial$df==0])
   args$lambda <- exp(seq(
     from=log(lambda.max),
     to=log(0.01*lambda.max),
     length.out=pmax(args$nlambda,100)))
   net <- do.call(what=glmnet::glmnet,args=args)
  }
  ### end original ###
  
  # # ### start trial (12 April 2019) ###
  # warning("TRIAL FUNCTION ACTIVE!",call.=FALSE)
  # max <- 10 # provide this as an argument!
  # iter <- 0
  # while((min(net$df)>3)|(length(net$lambda)==1)|(max(net$df)<=max)|(!max %in% net$df)){
  #   #warning("Modifying lambda sequence.",call.=FALSE)
  #   iter <- iter + 1
  #   if(iter>10){
  #     browser()
  #     stop("Modifying lambda sequence failed.",call.=FALSE)
  #   }
  #   # largest lambda not sparse
  #   if((min(net$df)>3)|(length(net$lambda)==1)){
  #     warning("largest lambda",call.=FALSE)
  #     args$lambda <- exp(seq(from=log(99e99),to=log(0.01),length.out=100))
  #     initial <- do.call(what=glmnet::glmnet,args=args)
  #     if(all(is.na(initial$lambda[initial$df==0]))){next}
  #     lambda.max <- 2*min(initial$lambda[initial$df==0]) # trial (2*)
  #     args$lambda <- exp(seq(
  #       from=log(lambda.max),
  #       to=log(0.01*lambda.max),
  #       length.out=pmax(args$nlambda,100)))
  #     net <- do.call(what=glmnet::glmnet,args=args)
  #     next
  #   }
  #   # smallest lambda too sparse
  #   #if(max(net$df)<=max){
  #   #  warning("smallest lambda",call.=FALSE)
  #   #  lambda.max <- max(net$lambda)
  #   #  lambda.min <- 0.0001*max(net$lambda)
  #   #  args$lambda <- exp(seq(from=log(lambda.max),
  #   #                         to=log(lambda.min),
  #   #                         length.out=pmax(args$nlambda,100)))
  #   #  net <- do.call(what=glmnet::glmnet,args=args)
  #   #  next
  #   #}
  #   # constraint excluded
  #   if(!max %in% net$df){
  #     warning("center",call.=FALSE)
  #     above <- net$lambda[net$df<(max)] # trial -2
  #     below <- net$lambda[net$df>(max)] # trial +2
  #     from <- log(min(above)); to <- log(max(below))
  #     if(is.na(to)){
  #       lambda.min <- min(net$lambda)
  #       extra <- exp(seq(from=log(lambda.min),
  #                        to=log(0.01*lambda.min),
  #                        length.out=pmax(args$nlambda,100)))
  #       args$lambda <- sort(c(above,extra),decreasing=TRUE)
  #     } else {
  #       focus <- exp(seq(from=from,to=to,length.out=50))
  #       args$lambda <- sort(c(above,focus,below),decreasing=TRUE)
  #     }
  #     net <- do.call(what=glmnet::glmnet,args=args)
  #   }
  #   break
  # }
  # # ### end trial (12 April 2019) ###
  
  ### start trial (on 11 April 2019) ###
  # max <- 10
  # if(!max %in% net$df){ # !max %in% net$df ### REACTIVATE THIS!
  #     warning("Refining lambda sequence!",call.=FALSE)
  #     above <- net$lambda[net$df<(max)] # trial -2
  #     below <- net$lambda[net$df>(max)] # trial +2
  #     from <- log(min(above))
  #     to <- log(max(below))
  #     if(is.na(to)){
  #       lambda.max <- max(net$lambda)
  #       args$lambda <- exp(seq(from=log(lambda.max),
  #                              to=log(0.0001*lambda.max),
  #                              length.out=pmax(args$nlambda,100)))
  #     } else {
  #       focus <- exp(seq(from=from,to=to,length.out=50))
  #       args$lambda <- sort(c(above,focus,below),decreasing=TRUE)
  #     }
  #     net <- do.call(what=glmnet::glmnet,args=args)
  # }
  # later: Replace "max <- 10" by argument "max"!
  # later: Split nlambda in two parts, e.g. 50 + 50? (No good idea!)
  # now: Put this inside security loop (above)!
  # now: Execute this conditionally (if max not in net$df)! (Done.)
  # Think about lambda.min.ratio! (Rather not!)
  ### end trial (on 11 April 2019) ###
  
  return(net)
}


#' @title Correlation
#' 
#' @description
#' Calculates the correlation between the response and the covariates.
#' Shrinks the correlation coefficients for each covariate set separately.
#' 
#' @inheritParams arguments
#' 
#' @param y
#' vector of length \eqn{n}
#' 
#' @param x
#' matrix with \eqn{n} rows and \eqn{p} columns
#' 
#' @return list of vectors
#' 
#' @examples
#' NA
#' 
.cor <- function(y,x,args){
    if(args$family=="cox"){
        # cor <- apply(X=x,MARGIN=2,FUN=function(x) abs(2*survival::survConcordance(y~x)$concordance-1)) # will depreciate
        cor <- apply(X=x,MARGIN=2,FUN=function(x) abs(2*survival::concordance(y~x)$concordance-1))
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

#' @title Weighting schemes
#' @description
#' Calculates the weighting schemes.
#' 
#' @inheritParams arguments
#' 
#' @return list of named vectors (one for each weighting scheme)
#' 
#' @examples
#' NA
#' 
.weight <- function(cor,args){
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

#' @title Cross-validation
#' 
#' @description
#' Repeatedly leaves out samples, and predicts their response.
#' 
#' @inheritParams arguments
#' 
#' @return Returns matrix of predicted values (except "cox")
#' 
#' @examples
#' NA
#' 
.cv <- function(y,x,foldid,lambda,args){
    fit <- lapply(seq_len(args$num),function(i) matrix(nrow=ifelse(args$family=="cox",args$nfolds,args$n),ncol=length(lambda[[i]])))
    for(i in seq_len(args$nfolds)){
        y0 <- y[foldid!=i]
        y1 <- y[foldid==i]
        X0 <- x[foldid!=i,,drop=FALSE]
        X1 <- x[foldid==i,,drop=FALSE]
        fit.sub <- .fit(y=y0,x=X0,args=args)
        for(j in seq_len(args$num)){
            if(args$family=="cox"){
                beta <- stats::predict(object=fit.sub[[j]],newx=X1,type="coeff",s=lambda[[j]])
                # check whether these are the betas!
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
                # check whether 0 < temp < 1 in binomial!
                if(any(is.na(temp))|(ncol(temp)==1)){
                    stop("Prediction problem.",call.=FALSE)
                }
                fit[[j]][foldid==i,seq_len(ncol(temp))] <- temp
            }
        }
    }
    return(fit) 
}

#' @title Cross-validation loss
#' 
#' @description
#' Calculates mean cross-validated loss
#' 
#' @inheritParams arguments
#' 
#' @return
#' Returns list of vectors, one for each model.
#' 
#' @examples
#' NA
#' 
.loss <- function(y,fit,family,type.measure,foldid=NULL){
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
                limit <- 1e-05 # examine boundaries?
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
                cvraw <- matrix(data=NA,nrow=length(weights),
                                ncol=ncol(fit[[i]])) # new
                                # ncol=ncol(fit)) # old
                for(k in seq_along(weights)){
                    cvraw[k,] <- apply(X=fit[[i]],MARGIN=2,FUN=function(x) glmnet.auc(y=y[foldid==k],prob=x[foldid==k]))
                }
                loss[[i]] <- apply(X=cvraw,MARGIN=2,FUN=function(x) stats::weighted.mean(x=x,w=weights,na.rm=TRUE))
                names(loss[[i]]) <- colnames(fit[[i]]) # new
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
        if(sum(diff(is.na(loss[[i]])))==1){
            loss[[i]] <- loss[[i]][!is.na(loss[[i]])]
        }
    }
    return(loss)
}

#' @title Extraction
#' 
#' @description
#' Extracts \code{cv.glmnet}-like object.
#' 
#' @inheritParams arguments
#' 
#' @examples
#' NA
#' 
.extract <- function(fit,lambda,cvm,type.measure){
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
        #model[[i]]$lambda.1se <-  # add?
    }
    return(model)
}

##' Toydata
##' 
##' Dataset for reproducing the vignette.
##' 
##' @docType data
##' @keywords internal
##' @name toydata
##' @usage data(toydata)
##' @return All entries are numeric.
##' @format A list of numeric vectors and numeric matrices.
#NULL

