
#--- Visualisation -------------------------------------------------------------

#' @name plots
#' 
#' @title
#' Plot functions (manuscript)
#' 
#' @description
#' Functions for the \code{palasso} manuscript.
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
#' The function \code{plot_score} compares a selected column to each of the
#' other columns. It counts the number of rows where the entry in the selected
#' column is smaller (blue), equal (white), or larger (red).
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
    palasso:::.mtext(text=colnames(X)[-choice],side=1)
    
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
#' n <- 5; p <- 3
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' palasso:::plot_table(X,margin=2)
#' 
plot_table <- function(X,margin=2,labels=TRUE,las=1){
    #par <- graphics::par(no.readonly=TRUE)
    
    n <- nrow(X); p <- ncol(X)
    if(is.null(rownames(X))){rownames(X) <- seq_len(n)}
    if(is.null(colnames(X))){colnames(X) <- seq_len(p)}
    
    v <- 0.5/(n-1)
    h <- 0.5/(p-1)
    
    graphics::plot.new()
    graphics::plot.window(xlim=c(-h,1+h),ylim=c(-v,1+v))
    par_usr <- graphics::par()$usr
    graphics::par(usr=c(-h,1+h,-v,1+v))
    palasso:::.mtext(text=rev(rownames(X)),unit=TRUE,side=2,las=las)
    palasso:::.mtext(text=colnames(X),unit=TRUE,side=3,las=las)
    
    if(margin==0){
        temp <- matrix(rank(X),nrow=n,ncol=p) # overall rank
    }
    if(margin==1){
        temp <- t(apply(X,1,rank)) # rank per row
    }
    if(margin==2){
        temp <- apply(X,2,rank) # rank per column
    }
    
    temp[is.na(X)] <- NA # temporary
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
    
    if(labels){
        labels <- round(as.numeric(X),digits=2)
        xs <- rep(seq_len(p),each=n)
        ys <- rep(seq_len(n),times=p)
        graphics::text(x=(xs-1)/(p-1),y=(n-ys)/(n-1),labels=labels,col="white")
    }
    
    graphics::par(usr=par_usr)
}

#' @rdname plots
#' @keywords internal
#' @examples
#' ### circle ###
#' n <- 50; p <- 25
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' Z <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' b <- sapply(seq_len(p),function(i) abs(cor(X[,i],Z[,i])))
#' w <- pmax(abs(cor(X)),abs(cor(Z)),na.rm=TRUE)
#' palasso:::plot_circle(b,w,cutoff=0)
#' 
plot_circle <- function(b,w,cutoff=NULL,group=NULL){
    # par <- graphics::par(no.readonly=TRUE)
    
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
    par_usr <- graphics::par()$usr
    par_mar <- graphics::par()$mar
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
    
    graphics::par(usr=par_usr,mar=par_mar)
}

#' @rdname plots
#' @keywords internal
#' @examples
#' ### box ###
#' n <- 10; p <- 5
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' palasso:::plot_box(X,choice=5)
#' 
plot_box <- function(X,choice=NULL,ylab="",ylim=NULL){
    
    # input
    n <- nrow(X); p <- ncol(X)
    if(is.null(rownames(X))){rownames(X) <- seq_len(n)}
    if(is.null(colnames(X))){colnames(X) <- seq_len(p)}
    if(is.null(choice)){choice <- p}
    if(is.character(choice)){choice <- which(colnames(X)==choice)}
    
    col <- rep(x="#CD0000",times=p)
    col[choice] <- "#0000CD"
    
    graphics::plot.new()
    if(is.null(ylim)){ylim <- range(X)}
    graphics::plot.window(xlim=c(0.5,p+0.5),ylim=ylim)
    graphics::box()
    graphics::axis(side=2)
    graphics::title(ylab=ylab,line=2.5)
    palasso:::.mtext(text=colnames(X),side=1)
    
    if(ylab %in% c("AUC","auc")){
        graphics::abline(h=max(apply(X,2,mean)),col="grey",lty=2)
        graphics::abline(h=max(apply(X,2,stats::median)),col="grey")
    } else if(ylab %in% c("MSE","deviance","mse")){
        graphics::abline(h=min(apply(X,2,mean)),col="grey",lty=2)
        graphics::abline(h=min(apply(X,2,stats::median)),col="grey")
    }

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
    graphics::mtext(text="weight",side=1,line=2.2,at=0)
    graphics::mtext(text="X",side=1,line=2.2,at=-1,col=args$col[1])
    graphics::mtext(text="Z",side=1,line=2.2,at=+1,col=args$col[2])
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
    # par <- graphics::par(no.readonly=TRUE)
    
    # difference
    diff <- x - y
    cutoff <- stats::quantile(abs(diff),p=prob,na.rm=TRUE)
    index <- sign(diff)*(abs(diff) > cutoff) + 2
    pch <- c(16,1,16)[index]
    col <- c("blue","grey","red")[index]
    
    # visualisation
    par_mar <- graphics::par()$mar
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
    
    graphics::par(mar=par_mar)
}

.mtext <- function(text,unit=FALSE,side=1,las=1){
    p <- length(text)
    # separator
    pattern <- c("_",".","|","+","-",":","*","^","$"," ")
    number <- sapply(pattern,function(z) sum(grepl(x=text,pattern=paste0("\\",z))))
    split <- pattern[which.max(number)]
    # separation
    strsplit <- strsplit(x=text,split=split)
    groups <- sapply(strsplit,function(x) x[1])
    names <- sapply(strsplit,function(x) paste0(x[-1],collapse=split))
    # names[names==""] <- groups[names==""]
    # groups[groups==names] <- ""
    # location
    border <- which(groups[-1]!=groups[-p])+0.5
    temp <- c(0.5,border,p+0.5)
    centre <- 0.5*(temp[-1]+temp[-length(temp)])
    # checks
    cond <- logical()
    cond[1] <- length(unique(groups))>=1 # trial! was >1 instead of >=1
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
        graphics::mtext(text=names,side=side,at=at(seq_len(p)),las=las)
        if(length(unique(groups))!=1){
            graphics::mtext(text="|",side=side,at=at(border),font=2)
        }
        graphics::mtext(text=groups[centre],side=side,at=at(centre),line=1)
    } else {
        graphics::mtext(text=text,side=side,at=at(seq_len(p)),las=las)
    }
}

#--- Application ---------------------------------------------------------------

#' @title
#' Other functions (manuscript)
#' 
#' @name other
#' 
#' @description
#' Functions for the \code{palasso} manuscript.
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
#' @param trial development option
#' 
#' @param ...
#' arguments for \link[palasso]{palasso}
#' 
#' @details
#' 
#' The function \code{.prepare} pre-processes sequencing data. It removes
#' features with a low total abundance, adjusts for different library sizes,
#' binarises the counts with cutoff zero or takes the Anscombe transform, and
#' scales all covariates to mean zero and unit variance.
#' (If the cutoff differed from zero, we would have to first normalise,
#' then binarise, and finally filter the raw counts.)
#' 
#' The function \code{.simulate} simulates a response. Exploiting an
#' experimental covariate matrix, it allows for different numbers of non-zero
#' coefficients for X and Z.
#' 
#' The function \code{.predict} estimates the predictive performance of
#' different lasso models (standard X and/or Z, adaptive X and/or Z, paired
#' X and Z). Minimising the deviance, it returns the deviance and the area
#' under the curve. Currently, only the binomial family is implemented.
#' 
#' The function \code{.select} estimates the selective performance of different
#' lasso models (standard X and/or Z, adaptive X and/or Z, paired X and Z).
#' Limiting the number of covariates to \eqn{10}, it returns the number of
#' selected covariates, and the number of correctly selected covariates.
#' 
#' @return
#' to do
#' 
#' @seealso
#' Use \link[palasso]{palasso} to fit the paired lasso.
#' 
#' @examples 
#' set.seed(1)
#' n <- 50; p <- 100
#' X <- matrix(rpois(n*p,lambda=4),nrow=n,ncol=p)
#' x <- palasso:::.prepare(X)
#' y <- palasso:::.simulate(x,effects=c(1,2))
#' predict <- palasso:::.predict(y,x)
#' select <- palasso:::.select(y,x,attributes(y))
NULL

#' @rdname other
#' @keywords internal
#' 
.prepare <- function(X,cutoff=NULL,quantile=NULL,scale=TRUE){
    
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
        cutoff <- nrow(X) # original: 0.05*nrow(X)
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
    if(is.null(quantile)){
        Z[,] <- X > 0 # zero-indicator 
    } else {
        # Z[,] <- X > stats::quantile(X,probs=quantile) # quantile
        Z[,] <- t(apply(X,1,function(x) x > stats::quantile(x=x,probs=quantile)))
    }
    
    # transform X
    X <- 2*sqrt(X+3/8) # Anscombe transform
    
    # properties
    prop <- mean(Z==0)
    sdx <- apply(X,2,stats::sd)
    sdz <- apply(Z,2,stats::sd)
    
    if(scale){
        # scaling X
        X <- scale(X)
        cx <- apply(X,2,function(x) all(is.na(x)))
        X[,cx] <- 0
        # scaling Z
        Z <- scale(Z)
        cz <- apply(Z,2,function(z) all(is.na(z)))
        Z[,cz] <- 0
    }

    # return
    x <- list(X=X,Z=Z)
    attributes(x)$info <- data.frame(n=nrow(X),p=ncol(X),prop=prop)
    return(x)
}

#' @rdname other
#' @keywords internal
#' @examples
#' 
.simulate <- function(x,effects){
    
    ### start trial ###
    for(i in seq_along(x)){
        x[[i]] <- scale(x[[i]])
        x[[i]][is.na(x[[i]])] <- 0
    }
    ### end trial ###
    
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

#' @rdname other
#' @keywords internal
#' @examples
#' 
.predict <- function(y,X,nfolds.ext=5,nfolds.int=5,standard=TRUE,adaptive=TRUE,...){
    
    start <- Sys.time()
    
    # dimensionality
    n <- unique(sapply(X,nrow))
    p <- unique(sapply(X,ncol))
    k <- length(X)
    
    # external folds
    fold.ext <- rep(NA,times=length(y))
    fold.ext[y==0] <- sample(rep(seq_len(nfolds.ext),
                                 length.out=sum(y==0)))
    fold.ext[y==1] <- sample(rep(seq_len(nfolds.ext),
                                 length.out=sum(y==1)))
    
    model <- character()
    if(standard){model <- c(model,paste0("standard_",c("x","z","xz")))}
    if(adaptive){model <- c(model,paste0("adaptive_",c("x","z","xz")))}
    model <- c(model,paste0("among_",c("x","z","xz")),
               "between_xz","within_xz","paired","trial1","trial2")
    nzero <- c(3,4,5,10,15,20,25,50,Inf)
    
    # predictions
    pred <- matrix(list(rep(NA,times=n)),nrow=length(nzero),
                 ncol=length(model),dimnames=list(nzero,model))
    #pred <- matrix(NA,nrow=n,ncol=length(model))
    #temp <- rep(NA,times=length(model))
    #names(temp) <- names
    deviance <- auc <- mse <- mae <- class <- matrix(NA,nrow=length(nzero),
            ncol=length(model),dimnames=list(nzero,model))
    
    # cross-validation
    for(k in seq_len(nfolds.ext)){
        
        y0 <- y[fold.ext!=k]
        X0 <- lapply(X,function(x) x[fold.ext!=k,,drop=FALSE])
        X1 <- lapply(X,function(x) x[fold.ext==k,,drop=FALSE])
        
        # internal folds
        fold.int <- rep(NA,times=length(y0))
        fold.int[y0==0] <- sample(rep(seq_len(nfolds.int),
                                      length.out=sum(y0==0)))
        fold.int[y0==1] <- sample(rep(seq_len(nfolds.int),
                                      length.out=sum(y0==1)))
        
        object <- palasso::palasso(y=y0,X=X0,foldid=fold.int,family="binomial",
                                   standard=standard,adaptive=adaptive,...)
       
        for(i in seq_along(nzero)){
            for(j in seq_along(model)){
                temp <- palasso:::predict.palasso(object=object,
                    newdata=X1,model=model[j],type="response",max=nzero[i])
                pred[i,j][[1]][fold.ext==k] <- temp
            }
        }
    }
    
    for(i in seq_along(nzero)){
        for(j in seq_along(model)){
            y_hat <- pred[i,j][[1]]
            mse[i,j] <- mean((y_hat-y)^2)
            mae[i,j] <- mean(abs(y_hat-y))
            class[i,j] <- mean(abs(round(y_hat)-y))
            y_hat <- pmax(1e-05,pmin(y_hat,1-1e-05))
            deviance[i,j] <- mean(-2*(y*log(y_hat)+(1-y)*log(1-y_hat)))
            auc[i,j] <- pROC::roc(response=y,predictor=y_hat)$auc
        }
    }
    
    end <- Sys.time()
    
    info <- data.frame(nfolds.ext=nfolds.ext,nfolds.int=nfolds.int,
                       time=format(end-start))
    list <- list(info=info,deviance=deviance,auc=auc,mse=mse,mae=mae,class=class)
    return(list)
}

#' @rdname other
#' @keywords internal
#' @examples
#' 
.select <- function(y,X,index,nfolds=5,standard=TRUE,adaptive=TRUE,...){
    
    fit <- palasso::palasso(y=y,X=X,family="binomial",nfolds=nfolds,
                            standard=standard,adaptive=adaptive,...)
    
    names <- unique(c(names(fit),"paired","trial1","trial2"))
    nzero <- c(3,4,5,10,15,20,25,50,Inf)
    
    shots <- hits1 <- hits2 <- matrix(integer(),
                            nrow=length(nzero),
                            ncol=length(names),
                            dimnames=list(nzero,names))
    
    for(i in seq_along(nzero)){
        for(j in seq_along(names)){
            coef <- palasso:::coef.palasso(fit,model=names[j],max=nzero[i])
            temp <- unlist(lapply(coef,function(x) Matrix::which(x!=0)))
            shots[i,j] <- length(temp)
            hits1[i,j] <- sum(unique(temp) %in% unlist(index)) # original
            hits2[i,j] <- sum(temp %in% unlist(index)) # new
        }
    }

    list <- list(shots=shots,hits1=hits1,hits2=hits2)
    return(list)
}

.logo <- function(...){
cat("
     __   __       __   __   __   __ 
    |__| |__| |   |__| |__  |__  |  |
    |    |  | |__ |  |  __|  __| |__|
    ")
}

