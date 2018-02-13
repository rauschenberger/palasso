
#--- Workhorse function --------------------------------------------------------

#' @title
#' Paired lasso
#' 
#' @export
#' 
#' @description
#' The function \code{palasso} cross-validates
#' the paired lasso.
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
#' @param ...
#' further arguments for \code{\link[glmnet]{cv.glmnet}}
#' or \code{\link[glmnet]{glmnet}}
#' 
#' @details
#' Let \code{x} denote one entry of the list \code{X}.
#' See \link[glmnet]{glmnet} for alternative
#' specifications of \code{y} and \code{x}.
#' 
#' @return
#' This function returns an object of class \code{palasso}.
#' 
#' @seealso
#' Available \code{\link[palasso]{methods}} for class \code{palasso} are
#' \code{\link[=predict.palasso]{predict}},
#' \code{\link[=coef.palasso]{coef}},
#' \code{\link[=fitted.palasso]{fitted}},
#' \code{\link[=weights.palasso]{weights}},
#' \code{\link[=deviance.palasso]{deviance}},
#' \code{\link[=summary.palasso]{summary}},
#' and \code{\link[=subset.palasso]{subset}}.
#' 
#' This package also includes hidden functions for
#' \code{\link[=extra]{comparing methods}} and
#' \code{\link[=plots]{plotting results}}.
#' 
#' @examples
#' set.seed(1)
#' n <- 100; p <- 200
#' y <- rbinom(n=n,size=1,prob=0.5)
#' X <- lapply(1:2,function(x) matrix(rnorm(n*p),nrow=n,ncol=p))
#' fit <- palasso(y=y,X=X,family="binomial",trial=TRUE)
#' 
palasso <- function(y,X,trial=FALSE,...){

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
    if(!base$family %in% c("gaussian","binomial","poisson")){stop("Invalid family.")}
    
    # dimensionality
    # n <- length(y) # Watch out! Object y must be a vector!
    k <- ifelse(is.list(X),length(X),1)
    n <- nrow(X[[1]])
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
    
    weight <- model <- list()
    
    # weight (correlation)
    mar <- list()
    for(i in seq_len(k)){
        mar[[i]] <- as.vector(abs(stats::cor(X[[i]],y)))
        mar[[i]][is.na(mar[[i]])] <- 0
    }
    
    # # marginal effects (univariate regression)
    # family <- eval(parse(text=base$family))()
    # mar <- list()
    # for(i in seq_len(k)){
    #     mar[[i]] <- abs(apply(X[[i]],2,function(x) stats::glm.fit(y=y,x=cbind(1,x),family=family)$coefficients[2]))
    #     mar[[i]][is.na(mar[[i]])] <- 0
    # }
    
    # # marginal effects (external)
    # mar <- ext
    
    # standard lasso
    for(i in seq_len(k)){
        weight[[i]] <- rep(1*(seq_len(k)==i),each=p)
    }
    weight[[k+1]] <- rep(1/k,times=k*p)
    
    # paired lasso
    weight[[k+2]] <- rep(0.5*sapply(mar,mean)/mean(unlist(mar)),each=p)
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
        model[[i]] <- tryCatch(do.call(what=glmnet::cv.glmnet,args=args),
                               error=function(x) NA)
        if(class(model[[i]])!="cv.glmnet"){
            # intercept-only model (verify this)
            temp <- base
            temp$lambda <- c(99e99,99e98)
            model[[i]] <- do.call(what=glmnet::cv.glmnet,args=temp)
        }
        ### trial start ###
        if(i > 1){
            model[[i]]$glmnet.fit$call$x <- NULL # save memory!
        }
        ### trial end ###
    }
    
    if(!trial){
        # optimal choice
        cvm <- sapply(model,function(x) x$cvm[x$lambda==x$lambda.min])
        # try: cvm <- sapply(model,function(x) min(x$cvm)) # not good!
        if(base$type.measure=="auc"){
            i <- which.max(cvm)
        } else {
            i <- which.min(cvm)
        }
        # output
        if(k==2){
            tryCatch(palasso:::plot_pairs(x=weight[[i]][1:p],
                                     y=weight[[i]][(p+1):(2*p)],
                                     main=paste0("i = ",i)),error=function(x) NULL)
        }
        cat("\n i =",i)
        model[[i]]$weight <- weight[[i]]
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
        class(model) <- "palasso"
        return(model)
    }
}

## trial start ##
#if(FALSE){
#print(object.size(fit),units="Mb")
#print(object.size(fit[[1]]$glmnet.fit$call$x)*8,units="Mb")
#x <- lapply(fit,function(x) x$glmnet.fit$call$x)
#all(x[[1]]==x[[2]])
#}
## trial end ##



#--- Generic functions ---------------------------------------------------------

#' @name methods
#' 
#' @title
#' Methods for class "palasso"
#' 
#' @description
#' generic functions
#' 
#' @param object,x
#' \link[palasso]{palasso} object
#' 
#' @param newdata
#' covariates\strong{:}
#' list of matrices with \eqn{n} rows (samples)
#' and \eqn{p} columns (variables)
#' 
#' @param s penalty parameter\strong{:}
#' character "lambda.min" or "lambda.1se",
#' or positive numeric
#' 
#' @param model
#' character "paired"
#' 
#' @param ... further arguments for \link[glmnet]{predict.cv.glmnet},
#' \link[glmnet]{coef.cv.glmnet}, or \link[glmnet]{deviance.glmnet}
#' # model character "paired", "standard", "adaptive"
#' # groups character "X" or "Z", or NULL (all groups)
#' 
#' @details
#' to do
#' 
#' @return
#' to do
#' 
#' @seealso
#' Use \link[palasso]{palasso} to fit the pairwise-adaptive lasso.
#' 
NULL

#' @rdname methods
#' @export
#' 
subset.palasso <- function(x,model="paired",...){
    object <- x
    if(length(list(...))!=0){warning("Ignoring argument.")}
    
    if(!inherits(x=object,what="palasso")){
        warning("Fake palasso object?")
    }
    
    name <- unique(sapply(X=object,FUN=function(x) x$name))
    if(length(name)!=1){
        stop("Different loss functions!")
    }
    
    if(model=="paired"){
        pattern <- "adaptive|between|within"
        cond <- grepl(pattern=pattern,x=names(object))
    } else {
        cond <- names(object)==model
    }
    
    temp <- object[[1]]$glmnet.fit$call$x # trial
    object <- object[cond]
    if(name=="AUC"){
        loss <- sapply(object,function(x) max(x$cvm))
        object <- object[[which.max(loss)]]
    } else {
        loss <- sapply(object,function(x) min(x$cvm))
        object <- object[[which.min(loss)]]
    }
    object$glmnet.fit$call$x <- temp # trial
    
    return(object)
}

#' @rdname methods
#' @export
#' 
predict.palasso <- function(object,newdata,s="lambda.min",model="paired",...){
    if(missing(newdata)||is.null(newdata)) {
        warning("Returning fitted values.")
        return(palasso:::fitted.palasso(object=object,s=s,model=model,...))
    }
    object <- palasso:::subset.palasso(x=object,model=model)
    newx <- do.call(what="cbind",args=newdata)
    glmnet::predict.cv.glmnet(object=object,newx=newx,s=s,...)
}

#' @rdname methods
#' @export
#' 
coef.palasso <- function(object,s="lambda.min",model="paired",...){
    object <- palasso:::subset.palasso(x=object,model=model)
    glmnet::coef.cv.glmnet(object=object,s=s,...)
}

#' @rdname methods
#' @export
#' 
fitted.palasso <- function(object,s="lambda.min",model="paired",...){
    object <- palasso:::subset.palasso(x=object,model=model)
    newx <- object$glmnet.fit$call$x
    glmnet::predict.cv.glmnet(object=object,newx=newx,s=s,...)
}

#' @rdname methods
#' @export
#' 
deviance.palasso <- function(object,model="paired",...){
    object <- palasso:::subset.palasso(x=object,model=model)
    glmnet::deviance.glmnet(object$glmnet.fit,...)
}

#' @rdname methods
#' @export
#' @importFrom stats weights
#' 
weights.palasso <- function(object,model="paired",...){
    object <- palasso:::subset.palasso(x=object,model=model)
    # Better split weights into X and Z group.
    # For this, adapt function <<subset.palasso>>.
    1/object$glmnet.fit$call$penalty.factor
}

#' @rdname methods
#' @export
#' 
summary.palasso <- function(object,model="paired",...){
    if(length(list(...))!=0){warning("Ignoring argument.")}
    
    object <- palasso:::subset.palasso(x=object,model=model)
    
    # head
    title <- paste(object$glmnet.fit$call$family,"palasso")
    line <- paste(rep("-",times=nchar(title)),collapse="")
    cat("",line,"\n",title,"\n",line,"\n\n")

    # dimensionality
    n <- object$glmnet.fit$dim[2]
    p <- object$glmnet.fit$dim[1]
    weight <- 1/object$glmnet.fit$call$penalty.factor
    p_in <- sum(weight!=0)
    p_out <- sum(weight==0)
    cat(paste0("dimensionality: n = ",n,", p = ",p))
    cat("\n \t       ",paste0("(",p_in," in, ",p_out," out)"),"\n\n")

    # cross-validation
    id <- list()
    id$min <- which(object$lambda==object$lambda.min)
    id$ose <- which(object$lambda==object$lambda.1se)
    frame <- matrix(NA,nrow=2,ncol=3)
    frame[,1] <- sapply(id,function(x) object$lambda[x])
    frame[,2] <- sapply(id,function(x) object$nzero[x])
    frame[,3] <- sapply(id,function(x) object$cvm[x])
    rownames(frame) <- c("prediction","explanation")
    colnames(frame) <- c("lambda","nzero",names(object$name))
    base::print(round(frame,digits=2))
    
    # weights
    cat("\n")
    M <- matrix(NA,nrow=1,ncol=3)
    M[1,1] <- stats::median(weight[weight!=0])
    M[1,2] <- mean(weight[weight!=0])
    M[1,3] <- max(weight)
    rownames(M) <- "weights (>0)"
    colnames(M) <- c("median","mean","max")
    base::print(round(M,digits=2))
    
    return(invisible("palasso"))
}

#--- Visualisation -------------------------------------------------------------

#' @name plots
#' 
#' @title
#' plots
#' 
#' @description
#' generic functions
#' 
#' @param X
#' matrix with \eqn{n} rows and \eqn{p} columns
#' 
#' @param x,y
#' vectors of equal length
#' 
#' @param choice
#' numeric between \eqn{1} and \eqn{p}
#' 
#' @param b
#' between-group correlation\strong{:}
#' vector of length \eqn{p}
#' 
#' @param w
#' within-group correlation\strong{:}
#' matrix with \eqn{p} rows and \eqn{p} columns
#' 
#' @param group
#' vector of length \eqn{p}
#' 
#' @param prob
#' confidence interval\strong{:}
#' numeric between \eqn{0} and \eqn{1}
#' 
#' @param cutoff
#' numeric between \eqn{0} and \eqn{1}
#' 
#' @param margin
#' \eqn{0} (none), \eqn{1} (rows), or \eqn{2} (columns)
#' 
#' @param ...
#' to do
#' 
#' @details
#' to do
#' 
#' @return
#' to do
#' 
#' @seealso
#' Use \link[palasso]{palasso} to fit the paired lasso.
NULL

# Graphical arguments\strong{:}
# line width \code{lwd},
# title \code{main},
# colours \code{col}
# 
# The function \code{scatter} plots the difference
# between two sets of measurements. The function \code{scales} plots weights.
# This function returns a plot.


#' @rdname plots
#' @keywords internal
#' @examples
#' ### score ###
#' 
#' n <- 10; p <- 4
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' palasso:::plot_score(X)
#' 
plot_score <- function(X,choice=NULL){
    
    # input
    n <- nrow(X); p <- ncol(X)
    if(is.null(rownames(X))){rownames(X) <- seq_len(n)}
    if(is.null(colnames(X))){colnames(X) <- seq_len(p)}
    if(is.null(choice)){choice <- p}
    if(is.character(choice)){choice <- which(colnames(X)==choice)}

    # score
    y <- list()
    temp <- X[,choice]
    y$gain <- apply(X,2,function(x) sum(temp<x))
    y$equal <- apply(X,2,function(x) sum(temp==x))
    y$loss <- apply(X,2,function(x) sum(temp>x))
    y <- lapply(y,function(x) x[-choice])
    
    # frame
    graphics::plot.new()
    graphics::plot.window(xlim=c(0.5,p-0.5),ylim=c(0,n))
    graphics::box()
    graphics::abline(h=n/2,lwd=2,col="grey")
    graphics::axis(side=2)
    graphics::title(ylab="count",line=2.5)
    palasso:::mtext(text=colnames(X)[-choice],side=1)
    
    # bars
    for(i in seq_len(p-1)){
        graphics::polygon(x=c(i-0.25,i-0.25,i+0.25,i+0.25),
                          y=c(0,y$gain[i],y$gain[i],0),
                          col="#0000CD")
        graphics::polygon(x=c(i-0.25,i-0.25,i+0.25,i+0.25),
                          y=c(y$gain[i],
                              y$gain[i]+y$equal[i],
                              y$gain[i]+y$equal[i],
                              y$gain[i]),
                          col="white")
        graphics::polygon(x=c(i-0.25,i-0.25,i+0.25,i+0.25),
                          y=c(y$gain[i]+y$equal[i],
                              y$gain[i]+y$equal[i]+y$loss[i],
                              y$gain[i]+y$equal[i]+y$loss[i],
                              y$gain[i]+y$equal[i]),
                          col="#CD0000")
    }
}

#' @rdname plots
#' @keywords internal
#' @examples
#' ### table ###
#' 
#' n <- 5; p <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' palasso:::plot_table(X,margin=1)
#' 
plot_table <- function(X,margin=2){
    
    n <- nrow(X); p <- ncol(X)
    if(is.null(rownames(X))){rownames(X) <- seq_len(n)}
    if(is.null(colnames(X))){colnames(X) <- seq_len(p)}
    
    v <- 0.5/(n-1)
    h <- 0.5/(p-1)
    
    graphics::plot.new()
    graphics::plot.window(xlim=c(-h,1+h),ylim=c(-v,1+v))
    graphics::par(usr=c(-h,1+h,-v,1+v))
    palasso:::mtext(text=rev(rownames(X)),unit=TRUE,side=2)
    palasso:::mtext(text=colnames(X),unit=TRUE,side=3)
    
    if(margin==0){
        temp <- matrix(rank(X),nrow=n,ncol=p) # overall rank
    }
    if(margin==1){
        temp <- t(apply(X,1,rank)) # rank per row
    }
    if(margin==2){
        temp <- apply(X,2,rank) # rank per column
    }
    
    image <- t(temp)[,seq(from=n,to=1,by=-1)]
    col <- grDevices::colorRampPalette(colors=c("darkblue","red"))(n*p)
    graphics::image(x=image,col=col,add=TRUE)

    if(margin==1){
        graphics::segments(x0=-h,x1=1+h,
                           y0=seq(from=-v,to=1+v,by=2*v),
                           col="white",lwd=3)
    }
    if(margin==2){
        graphics::segments(x0=seq(from=-h,to=1+h,by=2*h),
                           y0=1+v,y1=0-v,
                           col="white",lwd=3)
    }
    
    labels <- round(as.numeric(X),digits=2)
    xs <- rep(seq_len(p),each=n)
    ys <- rep(seq_len(n),times=p)
    graphics::text(x=(xs-1)/(p-1),y=(n-ys)/(n-1),labels=labels,col="white")
}

#' @rdname plots
#' @keywords internal
#' @examples
#' ### circle ###
#' 
#' n <- 50; p <- 100
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Z <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' 
#' b <- sapply(seq_len(p),function(i) abs(cor(X[,i],Z[,i])))
#' w <- pmax(abs(cor(X)),abs(cor(Z)),na.rm=TRUE)
#' 
#' palasso:::plot_circle(b,w)
#' 
plot_circle <- function(b,w,cutoff=NULL,group=NULL){
    
    # checks
    if(any(dim(w)!=length(b))){stop("Invalid dimensions!")}
    if(any(w!=t(w))){stop("Matrix X is asymmetric!")}
    w[row(w)>=col(w)] <- NA
    p <- length(b)
    if(is.null(cutoff)){
        cutoff <- numeric()
        cutoff[1] <- sort(b,decreasing=TRUE)[0.05*p]
        cutoff[2] <- sort(w,decreasing=TRUE)[0.05*p]
    }
    if(length(cutoff)==1){
        cutoff <- c(cutoff,cutoff)
    }
    
    # cutoff
    id <- list()
    id$b <- which(b>cutoff[1])
    id$w <- as.data.frame(which(w>cutoff[2],arr.ind=TRUE))
    
    # coordinates
    degree <- seq(from=0,to=2*pi,length.out=p+1)
    y <- sin(degree); x <- cos(degree)
    
    # initialisation
    graphics::plot.new()
    graphics::plot.window(xlim=1.2*c(-1,1),ylim=1.2*c(-1,1))
    graphics::par(usr=1.2*c(-1,1,-1,1),mar=c(1,1,1,1))
    graphics::lines(x=x,y=y)
    graphics::lines(x=0.8*x,y=0.8*y)
    
    # between
   #  col <- "#0000CD"
    col <- grDevices::rgb(red=0,green=0,blue=205,maxColorValue=255,alpha=50)
    graphics::segments(x0=x[id$b],y0=y[id$b],
                       x1=0.8*x[id$b],y1=0.8*y[id$b],col=col)
    graphics::segments(x0=0.8*x[id$w[,1]],y0=0.8*y[id$w[,1]],
                       x1=0.8*x[id$w[,2]],y1=0.8*y[id$w[,2]],col=col)
    
    # annotation
    if(!is.null(group)){
        names <- unique(group)
        size <- sapply(names,function(x) sum(x==group))
        cumsum <- cumsum(size)
        border <- cumsum+0.5
        graphics::segments(x0=x[border],y0=y[border],x1=1.1*x[border],y1=1.1*y[border],lwd=2)
        centre <- 0.5*(c(0,cumsum[-length(cumsum)])+cumsum)
        graphics::text(x=1.15*x[centre],y=1.15*y[centre],labels=names(cumsum))
    }
    
}

#' @rdname plots
#' @keywords internal
#' @examples
#' ### box ###
#' n <- 10; p <- 5
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' palasso:::plot_box(X,choice=5)
#' 
plot_box <- function(X,choice=NULL){
    
    # input
    n <- nrow(X); p <- ncol(X)
    if(is.null(rownames(X))){rownames(X) <- seq_len(n)}
    if(is.null(colnames(X))){colnames(X) <- seq_len(p)}
    if(is.null(choice)){choice <- p}
    if(is.character(choice)){choice <- which(colnames(X)==choice)}
    
    col <- rep(x="#CD0000",times=p)
    col[choice] <- "#0000CD"

    graphics::plot.new()
    graphics::plot.window(xlim=c(0.5,p+0.5),ylim=range(X))
    graphics::box()
    graphics::axis(side=2)
    graphics::title(ylab="",line=2.5)
    palasso:::mtext(text=colnames(X),side=1)
    
    for(i in seq_len(p)){
        # vioplot::vioplot(X[,i],at=i,add=TRUE,col="white")
        graphics::boxplot(x=X[,i],at=i,add=TRUE,col=col[i],boxwex=1)
        graphics::points(y=mean(X[,i]),x=i,col="white",pch=16)
    }
    
}

#' @rdname plots
#' @keywords internal
#' @examples
#' ### pairs ###
#' n <- 10
#' x <- runif(n)
#' y <- runif(n)
#' palasso:::plot_pairs(x,y)
#' 
plot_pairs <- function(x,y=NULL,...){
    
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

#' @rdname plots
#' @keywords internal
#' @examples
#' ### diff ###
#' n <- 100
#' x <- runif(n)
#' y <- runif(n)
#' palasso:::plot_diff(x,y)
#' 
plot_diff <- function(x,y,prob=0.95,...){
    
    # difference
    diff <- x - y
    cutoff <- stats::quantile(abs(diff),p=prob,na.rm=TRUE)
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






##--- Internal functions ----

#' @title
#' Split
#' 
#' @keywords internal
#' 
#' @description
#' generic functions
#' 
#' @param text character vector
#' 
#' @param unit logical
#' 
#' @param side to do
#' 
#' @details
#' to do
#' 
#' @return
#' to do
#' 
mtext <- function(text,unit=FALSE,side=1){
    p <- length(text)
    # separator
    pattern <- c("_",".","|","+","-",":","*","^","$"," ")
    number <- sapply(pattern,function(z) sum(grepl(x=text,pattern=paste0("\\",z))))
    split <- pattern[which.max(number)]
    # separation
    strsplit <- strsplit(x=text,split=split)
    groups <- sapply(strsplit,function(x) x[1])
    names <- sapply(strsplit,function(x) paste0(x[-1],collapse=split))
    names[names==""] <- groups[names==""]
    groups[groups==names] <- ""
    # location
    border <- which(groups[-1]!=groups[-p])+0.5
    temp <- c(0.5,border,p+0.5)
    centre <- 0.5*(temp[-1]+temp[-length(temp)])
    # checks
    cond <- logical()
    cond[1] <- length(unique(groups))>1
    cond[2] <- length(unique(groups))==length(border)+1
    cond[3] <- length(unique(groups))<p
    at <- function(x){
        if(unit){
            (x-1)/(p-1)
        } else {
            x
        }
    }
    if(all(cond)){
        graphics::mtext(text=names,side=side,at=at(seq_len(p)))
        graphics::mtext(text="|",side=side,at=at(border),font=2)
        graphics::mtext(text=groups[centre],side=side,at=at(centre),line=1)
    } else {
        graphics::mtext(text=text,side=side,at=at(seq_len(p)))
    }
}

#--- Application ---------------------------------------------------------------

#' @title
#' extra
#' 
#' @name extra
#' 
#' @description
#' generic functions
#' 
#' @param X
#' covariates\strong{:}
#' matrix with \eqn{n} rows and \eqn{p} columns
#' 
#' @param x
#' covariates\strong{:}
#' list of length \eqn{k},
#' including matrices with \eqn{n} rows and \eqn{p} columns
#' 
#' @param y
#' response\strong{:}
#' vector of length \eqn{n}
#' 
#' @param effects
#' number of causal covariates\strong{:}
#' vector of length \eqn{k}
#' 
#' @param index
#' indices of causal covariates\strong{:}
#' list of length \eqn{k},
#' including vectors
#' 
#' @param nfolds.ext
#' number of external folds
#' 
#' @param ...
#' arguments for \link[palasso]{palasso}
#' 
#' @details
#' By default, the function \code{predict} returns
#' the linear predictor (\code{type="link"}).
#' Consider predicting the response (\code{type="response"}).
#' 
#' @return
#' to do
#' 
#' @seealso
#' Use \link[palasso]{palasso} to fit the paired lasso.
#' 
#' @examples 
#' set.seed(1)
#' n <- 20; p <- 50
#' X <- matrix(rpois(n*p,lambda=4),nrow=n,ncol=p)
#' x <- .prepare(X)
#' y <- .simulate(x,effects=c(1,2))
#' .predict(y,x)
#' .select(y,x,attributes(y))
NULL

#' @rdname extra
#' @keywords internal
#' @export
#' @examples
#' 1+1
.prepare <- function(X,cutoff=NULL){
    
    # checks
    if(nrow(X)>=ncol(X)){
        stop("Low-dimensional data!")
    }
    if(any(X)<0 | any(X!=round(X))){
        stop("Raw counts required!")
    }
    
    # Remove features with low abundance.
    lib.size <- Matrix::rowSums(X)
    abundance <- Matrix::colSums(X)
    if(is.null(cutoff)){
        cutoff <- 0.05*nrow(X) 
    }
    X <- X[,abundance>=cutoff]
    X <- Matrix::as.matrix(X)
    
    # Adjust for different library sizes.
    norm.factors <- edgeR::calcNormFactors(object=t(X),lib.size=lib.size)
    gamma <- norm.factors*lib.size/mean(lib.size)
    gamma <- matrix(gamma,nrow=nrow(X),ncol=ncol(X))
    X <- X / gamma
    
    # transform Z
    Z <- matrix(integer(),nrow=nrow(X),ncol=ncol(X))
    Z[,] <- X > 0 # zero-indicator
    # alternative: Z[,] <- X > stats::quantile(X,p=0.8) # quantile
    prop <- mean(Z==0)
    
    # transform X
    X <- 2*sqrt(X+3/8) # Anscombe transform
    # alternative: X <- sqrt(X) # square root
    
    # scaling X
    X <- scale(X)
    cx <- apply(X,2,function(x) all(is.na(x)))
    X[,cx] <- 0
    
    # scaling Z
    Z <- scale(Z)
    cz <- apply(Z,2,function(z) all(is.na(z)))
    Z[,cz] <- 0
    
    # return
    x <- list(X=X,Z=Z)
    attributes(x)$info <- data.frame(n=nrow(X),p=ncol(X),prop=prop)
    return(x)
}

#' @rdname extra
#' @keywords internal
#' @export
#' @examples
#' 1+1
.simulate <- function(x,effects){
    
    # covariates
    if(length(x)!=length(effects)){stop("Invalid.")}
    if(ncol(unique(sapply(x,dim),MARGIN=2))!=1){stop("Invalid.")}
    k <- length(x)
    n <- nrow(x[[1]])
    p <- ncol(x[[1]])
    if(n>=p){warning("Low-dimensional data!")}
    
    # coefficients
    coef <- lapply(seq_len(k),function(i) 
        sample(rep(x=c(0,1),times=c(p-effects[i],effects[i]))))
    indices <- lapply(coef,function(x) which(x!=0))
    names(indices) <- names(x)
    
    # response
    eta <- rowSums(sapply(seq_len(k),function(i) x[[i]] %*% coef[[i]]))
    y <- stats::rbinom(n=n,size=1,prob=1/(1+exp(-eta)))
    # y <- stats::rnorm(n=n,mean=eta,sd=1)
    # y <- stats::rpois(n=n,lambda=exp(eta))
    
    attributes(y) <- indices
    return(y)
}

#' @rdname extra
#' @keywords internal
#' @export
#' @examples
#' 1+1
.predict <- function(y,X,pmax=NULL,nfolds.ext=5,nfolds.int=5){
    
    start <- Sys.time()
    
    # dimensionality
    p <- unique(sapply(X,ncol))
    k <- length(X)
    if(is.null(pmax)){pmax <- k*p}
    if(is.na(pmax)){pmax <- k*p}
    
    # external folds
    fold.ext <- rep(NA,times=length(y))
    fold.ext[y==0] <- sample(rep(seq_len(nfolds.ext),
                                 length.out=sum(y==0)))
    fold.ext[y==1] <- sample(rep(seq_len(nfolds.ext),
                                 length.out=sum(y==1)))
    
    # models
    names <- c(paste0("standard_",c("x","z","xz")),
            paste0("adaptive_",c("x","z","xz")),"paired")
    
    # predictions
    pred <- matrix(NA,nrow=length(y),ncol=length(names))
    deviance <- auc <- rep(NA,times=length(names))
    colnames(pred) <- names(deviance) <- names(auc) <- names

    # cross-validation
    for(i in seq_len(nfolds.ext)){
        
        y0 <- y[fold.ext!=i]
        X0 <- lapply(X,function(x) x[fold.ext!=i,,drop=FALSE])
        X1 <- lapply(X,function(x) x[fold.ext==i,,drop=FALSE])
        
        # internal folds
        fold.int <- rep(NA,times=length(y0))
        fold.int[y0==0] <- sample(rep(seq_len(nfolds.int),
                                      length.out=sum(y0==0)))
        fold.int[y0==1] <- sample(rep(seq_len(nfolds.int),
                                      length.out=sum(y0==1)))
        
        object <- palasso::palasso(y=y0,X=X0,foldid=fold.int,
                                  family="binomial",type.measure="deviance",
                                  pmax=pmax,trial=TRUE)
        
        pred[fold.ext==i,] <- sapply(names,function(x)
            palasso:::predict.palasso(object=object,newdata=X1,model=x,
                                      type="response"))
    }
    
    for(i in seq_along(names)){
        y_hat <- pmax(1e-05,pmin(pred[,i],1-1e-05))
        deviance[i] <- mean(-2*(y*log(y_hat)+(1-y)*log(1-y_hat)))
        auc[i] <- pROC::roc(response=y,predictor=y_hat)$auc
    }
    
    end <- Sys.time()
    
    info <- data.frame(nfolds.ext=nfolds.ext,nfolds.int=nfolds.int,
                       time=format(end-start))
    list <- list(info=info,deviance=deviance,auc=auc)
    
    return(list)
}

#' @rdname extra
#' @keywords internal
#' @export
#' @examples
#' 1+1
.select <- function(y,X,index,pmax=10,nfolds=5){
    
    p <- ncol(X[[1]])
    
    fit <- palasso::palasso(y=y,X=X,pmax=pmax,nfolds=nfolds,trial=TRUE)
    
    names <- c(paste0("standard_",c("x","z","xz")),
               paste0("adaptive_",c("x","z","xz")),"paired")
    
    coef <- sapply(names,function(i) palasso:::coef.palasso(object=fit,model=i)[-1])
    
    # coef <- split(x=coef,f=names) # trial A
    # select <- lapply(seq_along(names),function(i) which(coef[,i]!=0)) # trial B
    select <- apply(coef,2,function(x) which(x!=0))
    if(is.matrix(select)){select <- as.list(as.data.frame(select))}
    
    select <- lapply(select,function(x) x - p*(x %/% p))
    
    shots <- sapply(select,length)
    hits <- sapply(select,function(x) sum(unique(x) %in% unlist(index)))
    
    list <- list(shots=shots,hits=hits)
    return(list)
    
}

