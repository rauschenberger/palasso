
#--- Workhorse function --------------------------------------------------------

#' @title
#' Pairwise-adaptive lasso
#' 
#' @export
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
#' Available \code{\link[palasso]{methods}} for class \code{palasso}:
#' \code{\link[=predict.palasso]{predict}},
#' \code{\link[=coef.palasso]{coef}},
#' \code{\link[=fitted.palasso]{fitted}},
#' \code{\link[=weights.palasso]{weights}},
#' \code{\link[=deviance.palasso]{deviance}},
#' \code{\link[=summary.palasso]{summary}},
#' \code{\link[=subset.palasso]{subset}}
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
            tryCatch(palasso::scales(x=weight[[i]][1:p],
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

#--- Generic functions ---------------------------------------------------------

#' @name methods
#' 
#' @title
#' Methods for class "palasso"
#' 
#' @description
#' generic functions
#' 
#' @param object,x \link[palasso]{palasso} object
#' 
#' @param newdata covariates\strong{:}
#' list of matrices with \eqn{n} rows (samples)
#' and \eqn{p} columns (variables)
#' 
#' @param s penalty parameter\strong{:}
#' character "lambda.min" or "lambda.1se",
#' or positive numeric
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
subset.palasso <- function(x,...){
    object <- x
    if(length(list(...))!=0){warning("Ignoring argument.")}
    
    if(!inherits(x=object,what="palasso")){
        warning("Fake palasso object?")
    }
    
    name <- unique(sapply(X=object,FUN=function(x) x$name))
    if(length(name)!=1){
        stop("Different loss functions!")
    }
    
    cond <- grepl(pattern="adaptive|between|within",x=names(object))
    object <- object[cond]
    
    if(name=="AUC"){
        loss <- sapply(object,function(x) max(x$cvm))
        object <- object[[which.max(loss)]]
    } else {
        loss <- sapply(object,function(x) min(x$cvm))
        object <- object[[which.min(loss)]]
    }
    
    return(object)
}

#' @rdname methods
#' @export
#' 
predict.palasso <- function(object,newdata,s="lambda.min",...){
    if(missing(newdata)||is.null(newdata)) {
        return(palasso:::fitted.palasso(object=object,s=s,...))
    }
    object <- palasso:::subset.palasso(object)
    newx <- do.call(what="cbind",args=newdata)
    glmnet::predict.cv.glmnet(object=object,newx=newx,s=s,...)
}

#' @rdname methods
#' @export
#' 
coef.palasso <- function(object,s="lambda.min",...){
    object <- palasso:::subset.palasso(object)
    glmnet::coef.cv.glmnet(object=object,s=s,...)
}

#' @rdname methods
#' @export
#' 
fitted.palasso <- function(object,s="lambda.min",...){
    object <- palasso:::subset.palasso(object)
    newx <- object$glmnet.fit$call$x
    glmnet::predict.cv.glmnet(object=object,newx=newx,s=s,...)
}

#' @rdname methods
#' @export
#' 
deviance.palasso <- function(object,...){
    object <- palasso:::subset.palasso(object)
    glmnet::deviance.glmnet(object$glmnet.fit,...)
}

#' @rdname methods
#' @export
#' @importFrom stats weights
#' 
weights.palasso <- function(object,...){
    object <- palasso:::subset.palasso(object)
    # Better split weights into X and Z group.
    # For this, adapt function <<subset.palasso>>.
    1/object$glmnet.fit$call$penalty.factor
}

#' @rdname methods
#' @export
#' 
summary.palasso <- function(object,...){
    if(length(list(...))!=0){warning("Ignoring argument.")}
    
    object <- palasso:::subset.palasso(object)
    
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
#' Plots for class "palasso"
#' 
#' @description
#' generic functions
#' 
#' @param X,Y matrix with \eqn{n} rows and \eqn{p} columns
#' 
#' @param x,y vector
#' 
#' @param names vector of length \eqn{p}
#' 
#' @param groups vector of length \eqn{p}
#' 
#' @param choice numeric between \eqn{1} and \eqn{p}
#' 
#' @param ... to do
#' 
#' @details
#' to do
#' 
#' @return
#' to do
#' 
#' @seealso
#' Use \link[palasso]{palasso} to fit the pairwise-adaptive lasso.
NULL

#' @rdname plots
#' @export
#' 
plot_score <- function(X,choice){
    
    # input
    n <- nrow(X); p <- ncol(X)-1
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
    graphics::plot.window(xlim=c(0.5,p+0.5),ylim=c(0,n))
    graphics::box()
    graphics::abline(h=n/2,lwd=2,col="grey")
    graphics::axis(side=2)
    graphics::title(ylab="count",line=2.5)
    palasso:::mtext(text=colnames(X)[-choice])
    
    # bars
    for(i in seq_len(p)){
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
#' @export
#' 
plot_table <- function(X,margin=2,text=TRUE){
    
    n <- nrow(X)
    p <- ncol(X)
    
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
    
    if(text){
        labels <- round(as.numeric(X),digits=2)
        xs <- rep(seq_len(p),each=n)
        ys <- rep(seq_len(n),times=p)
        graphics::text(x=(xs-1)/(p-1),y=(n-ys)/(n-1),labels=labels,col="white")
    }
    
}




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



