
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
palasso <- function(y,X,trial=FALSE,ext=NULL,...){

    # checks
    base <- list(...)
    funs <- list(glmnet::glmnet,glmnet::cv.glmnet)
    formals <- unlist(lapply(funs,function(x) formals(x)))
    if(any(!names(base) %in% names(formals))){stop("Invalid argument.")}
    
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
    
    weights <- model <- list()
    
    # weights
    if(is.null(ext)){
        cor <- list()
        for(i in seq_len(k)){
            cor[[i]] <- as.vector(abs(stats::cor(X[[i]],y)))
            cor[[i]][is.na(cor[[i]])] <- 0
        }
    } else {
        cor <- ext
    }
        
    # standard lasso
    for(i in seq_len(k)){
        weights[[i]] <- rep(1*(seq_len(k)==i),each=p)
    }
    weights[[k+1]] <- rep(1/k,times=k*p)
    
    # paired lasso
    weights[[k+2]] <- rep(0.5*sapply(cor,mean)/mean(unlist(cor)),each=p)
    weights[[k+3]] <- unlist(cor)/rowSums(do.call(cbind,cor))
    weights[[k+3]][is.na(weights[[k+3]])] <- 0
    
    # adaptive lasso
    for(i in seq_len(k)){
        weights[[k+3+i]] <- weights[[i]]*cor[[i]]
    }
    weights[[2*k+4]] <- unlist(cor)
    
    # cross-validation
    args <- base
    for(i in seq_along(weights)){
        args$penalty.factor <- 1/weights[[i]]
        model[[i]] <- do.call(what=glmnet::cv.glmnet,args=args)
    }
    
    if(!trial){
        # optimal choice
        cvm <- sapply(model,function(x) x$cvm[x$lambda==x$lambda.min])
        if(base$type.measure=="auc"){
            i <- which.max(cvm)
        } else {
            i <- which.min(cvm)
        }
        # output
        if(k==2){
            tryCatch(palasso::scales(x=weights[[i]][1:p],
                                     y=weights[[i]][(p+1):(2*p)],
                                     main=paste0("i = ",i)),error=function(x) NULL)
        }
        cat("\n i =",i)
        model[[i]]$weights <- weights[[i]]
        return(model[[i]])
    }
    
    if(trial){
        if(length(model)!=8){stop("Not implemented!")}
        names(model) <- c("standard_x",
                          "standard_z",
                          "standard_xz",
                          "between_xz",
                          "within_xz",
                          "adaptive_x",
                          "adaptive_z",
                          "adaptive_xz")
        return(model)
    }
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

#' @title
#' Set comparison
#' 
#' @export
#' @keywords internal
#' 
#' @description
#' The function \code{scatter} plots the difference
#' between two sets of measurements.
#' 
#' @param x
#' vector of length \eqn{n}
#' 
#' @param y
#' vector of length \eqn{n}
#' 
#' @param p
#' confidence interval\strong{:}
#' numeric between \eqn{0} and \eqn{1}
#' 
#' @param ...
#' Graphical arguments
#' 
#' @return
#' This function returns a plot.
#' 
#' @examples
#' x <- runif(100)
#' y <- runif(100)
#' scatter(x,y)
#' 
scatter <- function(x,y,p=0.95,...){
    
    # difference
    diff <- x - y
    cutoff <- stats::quantile(abs(diff),p=p,na.rm=TRUE)
    index <- sign(diff)*(abs(diff) > cutoff) + 2
    pch <- c(16,1,16)[index]
    col <- c("blue","grey","red")[index]
    
    # visualisation
    graphics::par(mar=c(4,4,1,1))
    ylim <- 1.2*range(c(-diff,diff),na.rm=TRUE)
    graphics::plot(x=diff,ylim=ylim,
                   ylab="difference",col=col,pch=pch)
    graphics::abline(h=c(-1,0,1)*cutoff,lty=2)
    
    # legend: distribution
    mean <- signif(mean(diff,na.rm=TRUE),digits=3)
    median <- signif(median(diff,na.rm=TRUE),digits=3)
    sd <- signif(sd(diff,na.rm=TRUE),digits=3)
    m <- sum(!is.na(diff))
    l1 <- as.expression(bquote(bar(X)==.(mean)))
    l2 <- as.expression(bquote(tilde(X)==.(median)))
    l3 <- as.expression(bquote(hat(sigma)==.(sd)))
    l4 <- as.expression(bquote(ring(m)==.(m)))
    graphics::legend(x="topright",legend=c(l1,l2,l3,l4),bty="n")
    
    # legend: outliers
    nblue <- sum(col=="blue",na.rm=TRUE)
    nred <- sum(col=="red",na.rm=TRUE)
    graphics::legend(x="bottomleft",
                     legend=c(nblue," +",nred," =",nblue+nred),
                     text.col=c("blue","black","red","black","black"),
                     bty="n",
                     horiz=TRUE,
                     x.intersp=0.5,
                     text.width=1)
    
    # legend: p-values
    wilcox <- suppressWarnings(signif(stats::wilcox.test(diff,alternative="less",na.rm=TRUE)$p.value,digits=2))
    adhoc <- signif(1-stats::pbinom(q=nblue-1,size=nblue+nred,prob=0.5),digits=2)
    l1 <- paste0("wilcox = ",wilcox)
    l2 <- paste0("binom = ",adhoc)
    graphics::legend(x="bottomright",legend=c(l1,l2),bty="n")
}

