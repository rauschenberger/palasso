
#' @export
#' @title
#' Paired lasso regression
#' 
#' @aliases palasso-package
#' 
#' @description
#' trial
#' 
#' @param y
#' trial
#' 
#' @param X
#' trial
#' 
#' @param max
#' trial
#' 
#' @param ...
#' trial
#' 
#' @details
#' trial
#' 
#' @return
#' trial
#' 
#' @references
#' trial
#' 
#' @examples
#' if(FALSE){
#'     start <- Sys.time()
#'     for(i in 1:100){
#'         cat("i =",i,"\n")
#'         set.seed(i)
#'         n <- 50; p <- 20
#'         y <- rbinom(n=n,size=1,prob=0.5)
#'         X <- lapply(1:2,function(x) matrix(rnorm(n*p),nrow=n,ncol=p))
#'         object <- palasso(y=y,X=X,family="binomial")
#'         a <- predict(object,newdata=X)
#'         b <- coef(object)
#'         c <- weights(object)
#'         d <- fitted(object)
#'         e <- residuals(object)
#'         f <- deviance(object)
#'         g <- logLik(object)
#'         #h <- summary(object)
#'     }
#'     end <- Sys.time()
#'     end-start
#' }
#' 
palasso <- function(y=y,X=X,max=10,...){
    
    ### start function ###
    ## INCLUDE: test whether y and family are compatible
    ## INCLUDE: test whether any arguments are redundant
    ## ALSO: fuse functions args.trial and dim.trial
    args <- args.trial(...) # (...)
    args <- c(args,dim.trial(y=y,X=X,args=args))
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


#' @export
#' @title to do
#' @description to do
#' @param family
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

# param y
# param X
# param args
#' @export
#' @title
#' @description
#' @param ...
#' @details
#' @return
#' @examples
#' NA
dim.trial <- function(y,X,args=NULL){
    
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

#' @export
#' @title to do
#' @description to do
#' @param ...
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


#' @export
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


#' @export
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


#' @export
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

#' @export
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

#' @export
#' @title cross-validation
#' @description this function
#' @param y response
#' @param x covariates
#' @param foldid vector
#' @param lambda vector
#' @param args see ...
#' @return
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
                if(args$grouped){
                    plfull <- glmnet::coxnet.deviance(x=x,y=y,beta=beta)
                    plminusk <- glmnet::coxnet.deviance(x=X0,y=y0,beta=beta)
                    temp <- plfull - plminusk
                } else {
                    temp <- glmnet::coxnet.deviance(x=X1,y=y1,beta=beta)
                }
                fit[[j]][i,seq_along(temp)] <- temp
            } else {
                temp <- glmnet::predict.glmnet(object=fit.sub[[j]],newx=X1,type="response",s=lambda[[j]])
                if(any(is.na(temp))|(ncol(temp)==1)){
                    stop("Prediction problem.",call.=FALSE)
                }
                fit[[j]][foldid==i,seq_len(ncol(temp))] <- temp
            }
        }
    }
    return(fit) 
}

#' @export
#' @title to do
#' @description to do
#' @param y response
#' @param fit see ..
#' @param family ...
#' @param type.measure ...
#' @param foldid
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

#' @export
#' @title extract
#' @description this function ...
#' @param fit ...
#' @param lambda ...
#' @param cvm ...
#' @param type.measure ...
#' @details
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








# ------------------------------------------------------------------------------

if(FALSE){
    
    set.seed(2)
    old <- palasso::palasso(y=y,X=x,family="binomial",standard=TRUE,nfolds=5)
    sapply(old,function(x) signif(x$lambda.min,2))
    sapply(old,function(x) x$lambda.min==rev(x$lambda)[1])
    sapply(old,function(x) signif(x$cvm[x$lambda==x$lambda.min],2))
    
    set.seed(2)
    max <- max(sapply(test,function(x) max(x$lambda)))
    min <- min(sapply(test,function(x) min(x$lambda)))
    prop <- exp((log(min)-log(max))/100)
    lambda <- max*prop^seq_len(100)
    #lambda <- old$standard_x$lambda
    new <- palasso::palasso(y=y,X=x,family="binomial",lambda=lambda,standard=TRUE,nfolds=5)
    sapply(new,function(x) x$lambda.min==rev(x$lambda)[1])
    sapply(new,function(x) signif(x$lambda.min,2))
    sapply(new,function(x) signif(x$cvm[x$lambda==x$lambda.min],2))
    
}



########################
### multiple lambdas ### -------------------------------------------------------
########################
# 
# rm(list=ls())
# 
# sel <- c(c(1,0,0,1,1),c(1,1,1,1,1),c(1,0,0,1,1),c(1,0,0,0,0))
# family <- rep(c("gaussian","binomial","poisson","cox"),each=5)[sel==1]
# type.measure <- rep(c("deviance","class","auc","mse","mae"),times=4)[sel==1]
# all <- one <- list()
# lambda <- c(1,0.1,0.001,0.001) # adapt this
# grouped <- TRUE
#     
# for(i in seq_along(family)){
#     
#     list <- sim.trial(family=family[i],n=200,p=300)
#     #foldid <- sample(seq_len(10),replace=TRUE,size=200)
#     foldid <- rep(seq_len(10),length.out=200)
#     y <- list$y; X <- list(X=list$X,Z=list$Z)
#     
#     set.seed(1)
#     all[[i]] <- palasso.trial(y=y,X=X,
#                             #lambda=NULL,
#                             lambda=lambda, # remove
#                             foldid=foldid,
#                             family=family[i],
#                             type.measure=type.measure[i],
#                             grouped=grouped)
#     
#     set.seed(1)
#     one[[i]] <- glmnet::cv.glmnet(y=y,x=X$X,
#                                   lambda=lambda,
#                                   foldid=foldid,
#                                   family=family[i],
#                                   type.measure=type.measure[i],
#                                   grouped=grouped)
# 
#         #set.seed(1)
#         #old[[i]] <- glmnet::cv.glmnet(y=y,x=X,
#         #                          #lambda=new[[i]]$lambda,
#         #                          lambda=lambda, # remove
#         #                          foldid=foldid,
#         #                          family=family[i],
#         #                          type.measure=type.measure[i],
#         #                          grouped=grouped)
# 
# }
# 
# i <- 5
# all[[i]][[1]]$cvm
# one[[i]]$cvm
# 
# #sapply(new,function(x) x$cvm[2])
# #sapply(old,function(x) x$cvm[2])
# 
# # sapply(a,function(x) x[[1]])==sapply(b,function(x) x[[1]])
# # cond <- logical()
# # for(i in 1:12){
# #     cond[[i]] <- mean(abs(a[[i]]-b[[i]]))
# # }
# # data.frame(family=families,type=type.measures,a,b,same=1*(abs(a-b)<1e-06))



# y: response
# x: covariates
# family:
# lambda: vector
# shrink: logical
if(FALSE){
    list <- sim.trial(family="cox")
    y <- list$y; X <- list(X=list$X,Z=list$Z)
    max <- 10
    shrink <- TRUE
    pal <- palasso(y=y,X=X,family="cox")
}
#sapply(pal,function(x) x$cvm[x$lambda==x$lambda.min])

# authorised extra arguments:
# family, alpha
# lambda, nlambda, foldid, nfolds, shrink, standard, adaptive





# ### trial glmnet ###
# rm(list=ls())
# x=matrix(rnorm(100*20),100,20)
# y=rnorm(100)
# offset <- NULL
# lambda <- NULL
# type.measure <- "deviance"
# grouped <- TRUE
# keep <- FALSE
# parallel <- FALSE
# nfolds <- 10
# glmnet <- glmnet::glmnet
# cvtype <- glmnet:::cvtype
# cv.elnet <- glmnet:::cv.elnet
# 
# type.measure = match.arg(type.measure)
# if (!is.null(lambda) && length(lambda) < 2)
#     stop("Need more than one value of lambda for cv.glmnet")
# N = nrow(x)
# weights = rep(1, N)
# y = drop(y)
# glmnet.call = match.call(expand.dots = TRUE)
# which = match(c("type.measure", "nfolds", "foldid", "grouped",
#                 "keep"), names(glmnet.call), F)
# if (any(which))
#     glmnet.call = glmnet.call[-which]
# glmnet.call[[1]] = as.name("glmnet")
# glmnet.object = glmnet(x, y, weights = weights, offset = offset,
#                        lambda = lambda)
# glmnet.object$call = glmnet.call
# subclass=class(glmnet.object)[[1]]
# type.measure=cvtype(type.measure,subclass)
# is.offset = glmnet.object$offset
# ###Next line is commented out so each call generates its own lambda sequence
# # lambda=glmnet.object$lambda
# if (inherits(glmnet.object, "multnet") && !glmnet.object$grouped) {
#     nz = predict(glmnet.object, type = "nonzero")
#     nz = sapply(nz, function(x) sapply(x, length))
#     nz = ceiling(apply(nz, 1, median))
# } else nz = sapply(predict(glmnet.object, type = "nonzero"),
#                  length)
# foldid = sample(rep(seq(nfolds), length = N))
# if (nfolds < 3)
#     stop("nfolds must be bigger than 3; nfolds=10 recommended")
# outlist = as.list(seq(nfolds))
# if (parallel) {
#     #  if (parallel && require(foreach)) {
#     outlist = foreach(i = seq(nfolds), .packages = c("glmnet")) %dopar%
#     {
#         which = foldid == i
#         #      if (is.matrix(y))
#         if (length(dim(y))>1)
#             y_sub = y[!which, ]
#         else y_sub = y[!which]
#         if (is.offset)
#             offset_sub = as.matrix(offset)[!which, ]
#         else offset_sub = NULL
#         glmnet(x[!which, , drop = FALSE], y_sub, lambda = lambda,
#                offset = offset_sub, weights = weights[!which],
#                ...)
#     }
# } else {
#     for (i in seq(nfolds)) {
#         which = foldid == i
#         if (is.matrix(y))
#             y_sub = y[!which, ]
#         else y_sub = y[!which]
#         if (is.offset)
#             offset_sub = as.matrix(offset)[!which, ]
#         else offset_sub = NULL
#         outlist[[i]] = glmnet(x[!which, , drop = FALSE],
#                               y_sub, lambda = lambda, offset = offset_sub,
#                               weights = weights[!which])
#     }
# }
# fun = paste("cv", subclass, sep = ".")
# lambda = glmnet.object$lambda
# cvstuff = do.call(fun, list(outlist, lambda, x, y, weights,
#                             offset, foldid, type.measure, grouped, keep))
# cvm = cvstuff$cvm
# cvsd = cvstuff$cvsd
# nas=is.na(cvsd)
# if(any(nas)){
#     lambda=lambda[!nas]
#     cvm=cvm[!nas]
#     cvsd=cvsd[!nas]
#     nz=nz[!nas]
# }
# cvname = names(cvstuff$type.measure)
# names(cvname)=cvstuff$type.measure# to be compatible with earlier version; silly, I know
# out = list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm +
#                cvsd, cvlo = cvm - cvsd, nzero = nz, name = cvname, glmnet.fit = glmnet.object)
# if (keep)
#     out = c(out, list(fit.preval = cvstuff$fit.preval, foldid = foldid))
# lamin=if(cvname=="AUC")getmin(lambda,-cvm,cvsd)
# else getmin(lambda, cvm, cvsd)
# obj = c(out, as.list(lamin))
# class(obj) = "cv.glmnet"
# obj

# check cv.elnet !!!
