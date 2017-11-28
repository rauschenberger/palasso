
#' @title
#' Pairwise-adaptive lasso
#' 
#' @export
#' @keywords methods
#' 
#' @description
#' The function \code{palasso} cross-validates
#' the pairwise-adaptive lasso.
#' Use this regression technique if
#' the covariates are numerous and occur in pairs.
#' 
#' @param y
#' response\strong{:}
#' vector of length \eqn{n}
#' 
#' @param X
#' covariates\strong{:}
#' list of matrices with \eqn{n} rows (samples)
#' and \eqn{p} columns (variables)
#' 
#' @param oser
#' one-standard-error rule\strong{:}
#' logical (temporary argument)
#' 
#' @param trial
#' development mode\strong{:}
#' logical (temporary argument)
#' 
#' @param ext
#' external weights\strong{:}
#' vector of length \eqn{p}
#' 
#' @param ...
#' Further arguments for \link[glmnet]{cv.glmnet}
#' or \link[glmnet]{glmnet}.
#' 
#' @details
#' Let \code{x} denote one entry of the list \code{X}.
#' See \link[glmnet]{glmnet} for alternative
#' specifications of \code{y} and \code{x}.
#' 
#' @return
#' This function returns an object
#' of type \link[glmnet]{cv.glmnet},
#' with the additional slot \code{weights}.
#' 
#' @examples
#' n <- 100; p <- 200
#' y <- rbinom(n=n,size=1,prob=0.5)
#' X <- lapply(1:2,function(x) matrix(rnorm(n*p),nrow=n,ncol=p))
#' reg <- palasso(y=y,X=X,family="binomial")
#' 
palasso <- function(y,X,oser=FALSE,trial=FALSE,ext=NULL,...){
    
    # checks
    base <- list(...)
    #funs <- list(glmnet::glmnet,glmnet::cv.glmnet)
    #formals <- unlist(lapply(funs,function(x) formals(x)))
    #formals <- unlist(lapply(funs,function(x) names(formals(x))))
    #if(any(!names(base) %in% names(formals))){stop("Invalid argument.")}
    
    # arguments
    base$y <- y
    base$x <- do.call(what="cbind",args=X)
    default <- list(family="gaussian",alpha=1,nfolds=10,type.measure="deviance")
    base <- c(base,default[!names(default) %in% names(base)])

    # dimensionality
    n <- length(y)
    k <- ifelse(is.list(X),length(X),1)
    p <- ncol(X[[1]])
    
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
    
    # empty model
    args <- base
    args$alpha <- 1
    args$lambda <- c(2,1)
    args$dfmax <- 0
    args$pmax <- 0
    null <- do.call(what=glmnet::cv.glmnet,args=args)$cvm[1]
    
    if(is.null(ext)){
    
        # simple models
        args <- base
        args$alpha <- 0
        weights <- model <- coef <- list()
        for(i in seq_len(k)){
            weights[[i]] <- rep(1*(seq_len(k)==i),each=p)
            args$penalty.factor <- 1/weights[[i]]
            model[[i]] <- do.call(glmnet::cv.glmnet,args=args)
            if(i==1){args$lambda <- model[[1]]$lambda}
        }
        cvm <- lapply(model,function(x) x$cvm)
        num <- min(sapply(cvm,length))
        cvm <- rowSums(sapply(cvm,function(x) x[seq_len(num)]))
        min <- which.min(cvm)
        for(i in seq_len(k)){
            coef[[i]] <- abs(stats::coef(model[[i]]$glmnet.fit))[seq(from=(i-1)*p+2,to=i*p+1,by=1),min]
        }
        
        # equal weights
        weights[[k+1]] <- rep(1/k,times=k*p)
        
        # between-group weights
        alt <- sapply(model,function(x) x$cvm[x$lambda==x$lambda.min])
        prop <- (null-alt)/sum(null-alt)
        weights[[k+2]] <- rep(prop,each=p)
        
        # within-pair weights
        a <- do.call(what="c",args=coef)
        b <- rowSums(do.call(what="cbind",args=coef))
        weights[[k+3]] <- a / b
        weights[[k+3]][b==0] <- 0
    
    } else {
        # external weights
        weights <- model <- list()
        for(i in seq_len(k)){
            weights[[i]] <- rep(1*(seq_len(k)==i),each=p)
        }
        weights[[k+1]] <- rep(1/k,times=k*p)
        if(!is.list(ext)){ext <- list(ext)}
        for(i in seq_len(ext)){
            weights[[k+1+i]] <- ext[[i]]
        }
    }

    # cross-validation
    args <- base
    for(i in seq_along(weights)){
        if(base$alpha==0 & i <= length(model)){
             if(!is.null(model[[i]])){next}
        }
        args$penalty.factor <- 1/weights[[i]]
        model[[i]] <- do.call(what=glmnet::cv.glmnet,args=args)
    }
    
    # optimal choice
    cvm <- sapply(model,function(x) x$cvm[x$lambda==ifelse(oser,x$lambda.1se,x$lambda.min)])
    if(base$type.measure=="auc"){
        i <- which.max(cvm)
    } else {
        i <- which.min(cvm)
    }
    
    # output
    if(k==2){
        palasso::scales(x=weights[[i]][1:p],
                y=weights[[i]][(p+1):(2*p)],
                main=paste0("i = ",i))
    }
    cat("\n i =",i,"\n")
    model[[i]]$weights <- weights[[i]]
    
    return(model[[i]])
}

#' @title
#' Weight scales
#' 
#' @export
#' @keywords internal
#' 
#' @description
#' The function \code{scales} plots weights.
#' 
#' @param x
#' weights \strong{:}
#' vector of length \eqn{n}
#' 
#' @param y
#' weights \strong{:}
#' vector of length \eqn{n},
#' or \code{NULL} \eqn{(y=1-x)}
#' 
#' @param ...
#' Graphical arguments\strong{:}
#' line width \code{lwd},
#' title \code{main},
#' colours \code{col}
#' 
#' @return
#' This function returns a plot.
#' 
#' @examples
#' x <- runif(10)
#' scales(x)
#' 
scales <- function(x,y=NULL,...){
    
    # arguments
    args <- list(...)
    if(is.null(y)){y <- 1-x}
    if(is.null(args$lwd)){args$lwd <- 50/length(x)}
    if(is.null(args$main)){args$main <- ""}
    if(is.null(args$col)){args$col <- c("#0000CD","#CD0000")}
    
    # initialise
    graphics::plot.new()
    graphics::plot.window(xlim=c(-1,1),ylim=c(1,length(x)))
    
    # segments
    graphics::segments(x0=-x,x1=0,y0=seq_along(x),
                       col=args$col[1],lwd=args$lwd)
    graphics::segments(x0=+y,x1=0,y0=seq_along(x),
                       col=args$col[2],lwd=args$lwd)
    graphics::segments(x0=0,x1=0,y0=seq_along(x),
                       lwd=args$lwd,col="grey")
    
    # mean
    at <- seq(from=-mean(x),to=mean(y),length.out=1000)
    graphics::axis(side=1,at=at,labels=FALSE,col.ticks="grey")
    
    # terminate
    graphics::abline(v=0,lty=1)
    at <- seq(from=-1,to=1,by=0.5)
    graphics::axis(side=1,at=at,labels=abs(at))
    graphics::mtext(text="weight",side=1,line=2,at=0)
    graphics::mtext(text="X",side=1,line=2,at=-1,col=args$col[1])
    graphics::mtext(text="Z",side=1,line=2,at=+1,col=args$col[2])
    graphics::mtext(text=args$main,side=3,line=0,at=0)
}

