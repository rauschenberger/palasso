
glmnet.auc <- get("auc",envir=asNamespace("glmnet"))

#--- Visualisation -------------------------------------------------------------

#' @name plots
#' 
#' @title
#' Plot functions for manuscript
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
#' additional arguments
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
plot_score <- function(X,choice=NULL,ylab="count"){
    
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
    graphics::title(ylab=ylab,line=2.5)
    .mtext(text=colnames(X)[-choice],side=1)
    
    # bars
    for(i in seq_len(p-1)){
        graphics::polygon(x=c(i-0.25,i-0.25,i+0.25,i+0.25),
                          y=c(0,y$gain[i],y$gain[i],0),
                          col="#00007F") # was "#0000CD"
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
                          col="#FF3535") # was "CD0000"
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
plot_table <- function(X,margin=2,labels=TRUE,colour=TRUE,las=1,cex=1,cutoff=NA){
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
    .mtext(text=rev(rownames(X)),unit=TRUE,side=2,las=las,cex=cex)
    .mtext(text=colnames(X),unit=TRUE,side=3,las=las,cex=cex)
    
    if(margin==-1){
        breaks <- c(0,0.9,0.95,0.96,0.97,0.98,0.99,1)
        temp <- apply(X=X,MARGIN=1,FUN=function(x) cut(x=x,breaks=breaks,labels=seq_len(length(breaks)-1)))
        temp <- apply(X=temp,MARGIN=1,FUN=function(x) as.numeric(x))
    }
    if(margin==0){
        temp <- matrix(rank(X),nrow=n,ncol=p) # overall rank
    }
    if(margin==1){
        temp <- t(apply(X,1,rank)) # rank per row
    }
    if(margin==2){
        temp <- apply(X,2,rank) # rank per column
    }
    if(!is.na(cutoff)){
        temp <- X > cutoff
    }
    
    if(colour){
        temp[is.na(X)] <- NA # temporary
        image <- t(temp)[,seq(from=n,to=1,by=-1)]
        if(margin==-1){
            col <- grDevices::colorRampPalette(colors=c("navyblue","white"))(6)
        } else {
            col <- grDevices::colorRampPalette(colors=c("darkblue","red"))(n*p) 
        }
        graphics::image(x=image,col=col,add=TRUE)
    }
    
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
        labels <- format(labels,digits=2)
        xs <- rep(seq_len(p),each=n)
        ys <- rep(seq_len(n),times=p)
        graphics::text(x=(xs-1)/(p-1),y=(n-ys)/(n-1),labels=labels,
                       col=ifelse(colour,"white","black"))
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
    if(any(dim(w)!=length(b))){stop("Invalid dimensions!",call.=FALSE)}
    if(any(w!=t(w))){stop("Matrix X is asymmetric!",call.=FALSE)}
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
plot_box <- function(X,choice=NULL,ylab="",ylim=NULL,zero=FALSE,invert=FALSE){
    
    # input
    n <- nrow(X); p <- ncol(X)
    if(is.null(rownames(X))){rownames(X) <- seq_len(n)}
    if(is.null(colnames(X))){colnames(X) <- seq_len(p)}
    if(is.null(choice)){choice <- p}
    if(is.character(choice)){choice <- which(colnames(X)==choice)}
    
    graphics::plot.new()
    if(is.null(ylim)){ylim <- range(X)}
    graphics::plot.window(xlim=c(0.5,p+0.5),ylim=ylim)
    graphics::box()
    graphics::axis(side=2)
    graphics::title(ylab=ylab,line=2.5)
    .mtext(text=colnames(X),side=1)
    
    if(ylab %in% c("AUC","auc")){
        graphics::abline(h=max(apply(X,2,mean)),col="grey",lty=2)
        graphics::abline(h=max(apply(X,2,stats::median)),col="grey")
    } else if(ylab %in% c("MSE","deviance","mse")){
        graphics::abline(h=min(apply(X,2,mean)),col="grey",lty=2)
        graphics::abline(h=min(apply(X,2,stats::median)),col="grey")
    }
    
    if(zero){
        graphics::abline(h=0,col="grey",lty=1)
    }

    for(i in seq_len(p)){
        #vioplot::vioplot(X[,i],at=i,add=TRUE,col="white")
        #graphics::boxplot(x=X[,i],at=i,add=TRUE,col=col[i],boxwex=1)
        #graphics::points(y=mean(X[,i]),x=i,col="white",pch=16)
        .boxplot(x=X[,i],at=i,invert=invert)
    }
    
}


.boxplot <- function(x,at=1,wex=0.25,invert=FALSE){
    q <- stats::quantile(x,p=c(0.05,0.25,0.75,0.95))
    
    col <- c("#00007F","#FF3535")
    if(invert){col <- rev(col)}


    # outliers
    cond <- (x < q[1] | x > q[4]) & x > 0
    graphics::points(x=rep(at,sum(cond)),y=x[cond],col=col[2])
    cond <- (x < q[1] | x > q[4]) & x < 0
    graphics::points(x=rep(at,sum(cond)),y=x[cond],col=col[1])
    
    # box
    top <- max(0,q[3]); bot <- max(0,q[2])
    graphics::polygon(x=c(at-wex,at-wex,at+wex,at+wex),
            y=c(bot,top,top,bot),col=col[2],border=NA)
    top <- min(0,q[3]); bot <- min(0,q[2])
    graphics::polygon(x=c(at-wex,at-wex,at+wex,at+wex),
            y=c(bot,top,top,bot),col=col[1],border=NA)
    
    # median
    m <- stats::median(x)
    graphics::segments(x0=at-wex,x1=at+wex,y0=m,lwd=3,col="white")
    
    # box
    graphics::polygon(x=c(at-wex,at-wex,at+wex,at+wex),
            y=c(q[2],q[3],q[3],q[2]))
    
    # whiskers
    graphics::segments(x0=at-wex/2,x1=at+wex/2,y0=q[1],col="black",lwd=2)
    graphics::segments(x0=at,y0=q[1],y1=q[2],col="black",lwd=2)
    graphics::segments(x0=at,y0=q[3],y1=q[4],col="black",lwd=2)
    graphics::segments(x0=at-wex/2,x1=at+wex/2,y0=q[4],col="black",lwd=2)
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
    if(is.null(args$col)){args$col <- c("#00007F","#FF3535")}
    
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
plot_diff <- function(x,y,prob=0.95,ylab="",xlab="",...){
    # par <- graphics::par(no.readonly=TRUE)
    
    # difference
    diff <- x - y
    cutoff <- stats::quantile(abs(diff),p=prob,na.rm=TRUE)
    index <- sign(diff)*(abs(diff) > cutoff) + 2
    pch <- c(16,1,16)[index]
    col <- c("red","grey","blue")[index]
    
    # visualisation
    par_mar <- graphics::par()$mar
    #graphics::par(mar=c(4,4,1,1))
    ylim <- 1.2*range(c(-diff,diff),na.rm=TRUE)
    graphics::plot(x=diff,ylim=ylim,xlab=xlab,
                   ylab=ylab,col=col,pch=pch)
    graphics::abline(h=c(-1,1)*cutoff,lty=2)
    graphics::abline(h=0,lty=1)
    
    # legend: distribution
    mean <- signif(mean(diff,na.rm=TRUE),digits=3)
    median <- signif(median(diff,na.rm=TRUE),digits=3)
    sd <- signif(sd(diff,na.rm=TRUE),digits=3)
    m <- sum(!is.na(diff))
    l1 <- as.expression(bquote(bar(X)==.(mean)))
    l2 <- as.expression(bquote(tilde(X)==.(median)))
    l3 <- as.expression(bquote(hat(sigma)==.(sd)))
    l4 <- as.expression(bquote(ring(m)==.(m)))
    #graphics::legend(x="topright",legend=c(l1,l2,l3,l4),bty="n")
    
    # legend: outliers
    nblue <- sum(col=="blue",na.rm=TRUE)
    nred <- sum(col=="red",na.rm=TRUE)
    graphics::legend(x="bottomleft",
                     legend=c(nred,":",nblue), #," =",nblue+nred),
                     text.col=c("red","black","blue"), # "black","black"),
                     text.font=c(1,2,1),
                     bty="n",
                     horiz=TRUE,
                     x.intersp=0.5,
                     text.width=1)
    
    # legend: p-values
    wilcox <- suppressWarnings(signif(stats::wilcox.test(diff,alternative="less",na.rm=TRUE)$p.value,digits=2))
    adhoc <- signif(1-stats::pbinom(q=nblue-1,size=nblue+nred,prob=0.5),digits=2)
    l1 <- paste0("wilcox = ",wilcox)
    l2 <- paste0("binom = ",adhoc)
    #graphics::legend(x="bottomright",legend=c(l1,l2),bty="n")
    
    graphics::par(mar=par_mar)
}

.mtext <- function(text,unit=FALSE,side=1,las=1,cex=1){
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
    cond[1] <- length(unique(groups))>=1
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
        graphics::mtext(text=names,side=side,at=at(seq_len(p)),las=las,cex=cex)
        if(length(unique(groups))!=1){
            graphics::mtext(text="|",side=side,at=at(border),font=2,cex=cex)
        }
        graphics::mtext(text=groups[centre],side=side,at=at(centre),line=1,cex=cex)
    } else {
        graphics::mtext(text=text,side=side,at=at(seq_len(p)),las=las,cex=cex)
    }
}

#--- Application ---------------------------------------------------------------

#' @title
#' Analysis functions for manuscript
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
#' @param filter
#' numeric, multiplying the sample size
#' 
#' @param cutoff
#' character "zero", "knee", or "half"
#' 
#' @param scale
#' logical
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
#' 
#' \code{.prepare}\strong{:}
#' pre-processes sequencing data by
#' removing features with a low total abundance,
#' and adjusting for different library sizes;
#' obtains two transformations of the same data
#' by (1) binarising the counts with some cutoff
#' and (2) taking the Anscombe transform;
#' scales all covariates to mean zero and unit variance.
#' 
#' \code{.simulate}\strong{:}
#' simulates the response by
#' exploiting two experimental covariate matrices;
#' allows for different numbers of non-zero coefficients for X and Z.
#' 
#' \code{.predict}\strong{:}
#' estimates the predictive performance of different lasso models
#' (standard X and/or Z, adaptive X and/or Z, paired X and Z);
#' minimises the loss function "deviance", but also returns other loss functions;
#' supports logistic and Cox regression.
#' 
#' \code{.select}\strong{:}
#' estimates the selective performance of different lasso models
#' (standard X and/or Z, adaptive X and/or Z, paired X and Z);
#' limits the number of covariates to \eqn{10};
#' returns the number of selected covariates,
#' and the number of correctly selected covariates.
#' 
#' @seealso
#' Use \link[palasso]{palasso} to fit the paired lasso.
#' 
#' @examples
#' \dontrun{set.seed(1)
#' n <- 30; p <- 40
#' X <- matrix(rpois(n*p,lambda=3),nrow=n,ncol=p)
#' x <- palasso:::.prepare(X)
#' y <- palasso:::.simulate(x,effects=c(1,2))
#' predict <- palasso:::.predict(y,x)
#' select <- palasso:::.select(y,x,attributes(y))}
NULL

#' @rdname other
#' @keywords internal
#' 
.prepare <- function(X,filter=1,cutoff="zero",scale=TRUE){
    
    # checks
    if(nrow(X)>=ncol(X)){
        stop("Low-dimensional data!",call.=FALSE)
    }
    if(any(X)<0 | any(X!=round(X))){
        stop("Raw counts required!",call.=FALSE)
    }
    
    # Remove features with low abundance.
    lib.size <- Matrix::rowSums(X)
    abundance <- Matrix::colSums(X)
    X <- X[,abundance>=filter*nrow(X)]
    X <- Matrix::as.matrix(X)
    
    # Adjust for different library sizes.
    norm.factors <- edgeR::calcNormFactors(object=t(X),lib.size=lib.size)
    gamma <- norm.factors*lib.size/mean(lib.size)
    gamma <- matrix(gamma,nrow=nrow(X),ncol=ncol(X))
    X <- X / gamma
    
    # transform X
    temp <- X
    X <- 2*sqrt(X+3/8) # Anscombe transform
    
    # transform Z
    Z <- matrix(integer(),nrow=nrow(X),ncol=ncol(X))
    if(cutoff=="zero"){
        Z[,] <- temp > 0 # equivalent to X > 2*sqrt(3/8)
    } else if(cutoff=="knee"){
        Z[,] <- .knee(X)
    } else if(cutoff=="half"){
        Z[,] <- t(apply(X,1,function(x) x > stats::quantile(x=x,probs=0.5)))
    }

    # properties
    prop <- mean(Z==0)
    
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

.knee <- function(X){
    n <- nrow(X)
    p <- ncol(X)
    Z <- matrix(integer(),nrow=n,ncol=p)
    #prop <- seq(from=0,to=1,length.out=n-1) # old
    #weight <- log(prop)*log(1-prop) # old
    prop <- seq(from=1,to=n-1,by=1)/n
    weight <- -prop*log(prop,base=2)-(1-prop)*log(1-prop,base=2)
    for(j in seq_len(p)){
        step <- sort(X[,j])
        index <- which.max(diff(step)*weight)
        Z[,j] <- X[,j] > step[index]
        #plot(step); abline(v=index+0.5,col="blue",lty=2)
    }
    return(Z)
}

#' @rdname other
#' @keywords internal
#' 
.simulate <- function(x,effects){
    
    # constant covariates
    for(i in seq_along(x)){
        x[[i]] <- scale(x[[i]])
        x[[i]][is.na(x[[i]])] <- 0
    }
    
    # covariates
    if(length(x)!=length(effects)){stop("Invalid.",call.=FALSE)}
    if(ncol(unique(sapply(x,dim),MARGIN=2))!=1){stop("Invalid.",call.=FALSE)}
    k <- length(x)
    n <- nrow(x[[1]])
    p <- ncol(x[[1]])
    if(n>=p){warning("Low-dimensional data!",call.=FALSE)}
    
    # coefficients
    coef <- lapply(seq_len(k),function(i) 
        sample(rep(x=c(0,1),times=c(p-effects[i],effects[i]))))
    indices <- lapply(coef,function(x) which(x!=0))
    names(indices) <- names(x)
    
    # response
    eta <- rowSums(sapply(seq_len(k),function(i) x[[i]] %*% coef[[i]]))
    y <- stats::rbinom(n=n,size=1,prob=1/(1+exp(-eta))) # binomial
    # y <- stats::rnorm(n=n,mean=eta,sd=1) # gaussian
    # y <- stats::rpois(n=n,lambda=exp(eta)) # poisson
    
    attributes(y) <- indices
    return(y)
}

# # simulation all families (without covariate effects)
# simulation <- function(family,n=200,p=150){
#   X <- matrix(data=stats::rnorm(n=n*p),nrow=n,ncol=p)
#   Z <- matrix(data=stats::rnorm(n=n*p),nrow=n,ncol=p)
#   if(family=="gaussian"){
#     y <- stats::rnorm(n=n)
#   } else if(family=="binomial"){
#     y <- stats::rbinom(n=n,size=1,prob=0.2)
#   } else if(family=="poisson"){
#     y <- stats::rpois(n=n,lambda=4)
#   } else if(family=="cox"){
#     event <- stats::rbinom(n=n,size=1,prob=0.3)
#     time <- stats::rpois(n=n,lambda=4)+1
#     y <- survival::Surv(time=time,event=event)
#   } else {
#     stop("Invalid family.")
#   }
#   return(list(y=y,X=X,Z=Z))
# }

#' @rdname other
#' @keywords internal
#' 
.predict <- function(y,X,nfolds.ext=5,nfolds.int=5,adaptive=TRUE,
                     standard=TRUE,elastic=TRUE,shrink=TRUE,family="binomial",...){
    
    if(survival::is.Surv(y)!=(family=="cox")){stop("Survival?")}
    
    start <- Sys.time()
    
    # dimensionality
    n <- unique(sapply(X,nrow))
    p <- unique(sapply(X,ncol))
    k <- length(X)
    
    # external folds
    fold.ext <- .folds(y=y,nfolds=nfolds.ext)
    
    model <- character()
    if(adaptive){model <- c(model,paste0("adaptive_",c("x","z","xz")),
                            "within_xz","paired.adaptive")}
    if(standard){model <- c(model,paste0("standard_",c("x","z","xz")),
                            "between_xz","paired.standard")}
    if(adaptive&standard){model <- c(model,"paired.combined")}
    if(elastic){model <- c(model,"elastic","elastic95")}
    #if(elastic){model <- c(model,"elastic")}
    
    nzero <- c(3,4,5,10,15,20,25,50,Inf)
    
    # predictions
    if(family=="binomial"){
        pred <- matrix(list(rep(NA,times=n)),nrow=length(nzero),
                 ncol=length(model),dimnames=list(nzero,model))
        deviance <- auc <- mse <- mae <- class <- matrix(NA,nrow=length(nzero),
            ncol=length(model),dimnames=list(nzero,model))
    }
    if(family=="cox"){
        cvraw <- matrix(list(rep(NA,times=nfolds.ext)),nrow=length(nzero),
                       ncol=length(model),dimnames=list(nzero,model))
        loss <- matrix(NA,nrow=length(nzero),ncol=length(model),
                       dimnames=list(nzero,model))
    }
    
    # cross-validation
    for(k in seq_len(nfolds.ext)){
        
        y0 <- y[fold.ext!=k]
        X0 <- lapply(X,function(x) x[fold.ext!=k,,drop=FALSE])
        X1 <- lapply(X,function(x) x[fold.ext==k,,drop=FALSE])
        
        # internal folds
        fold.int <- .folds(y=y0,nfolds=nfolds.int)
        
        #if(elastic){
        #  x0 <- do.call(what="cbind",args=X0)
        #  x1 <- do.call(what="cbind",args=X1)
        #  elastic50 <- glmnet::cv.glmnet(alpha=0.50,y=y0,x=x0,foldid=fold.int,family=family,...)
        #  elastic95 <- glmnet::cv.glmnet(alpha=0.95,y=y0,x=x0,foldid=fold.int,family=family,...)
        #}
        
        #if(elastic){
        #  x0 <- do.call(what="cbind",args=X0)
        #  x1 <- do.call(what="cbind",args=X1)
        #  #enet <- palasso:::enet(y=y0,x=x0,alpha=c(0.25,0.5,0.75,1),foldid=fold.int,family=family,
        #  #                       dfmax=10,lambda.min.ratio=0.1,...)
        #  enet <- palasso:::enet(y=y0,x=x0,alpha=1,family=family,dfmax=10)
        #  # Set alpha to 0.5 or 0.95 !!!
        #}
        
        object <- palasso::palasso(y=y0,X=X0,foldid=fold.int,family=family,
                                   standard=standard,elastic=elastic,shrink=shrink,...)
        
        ### start trial ###
        max <- signif(sapply(object,function(x) max(x$lambda)),1)
        sel <- signif(sapply(object,function(x) x$lambda.min),1)
        min <- signif(sapply(object,function(x) min(x$lambda)),1)
        one <- paste0("[",max,",(",sel,"),",min,"]")
        min <- sapply(object,function(x) min(x$nzero))
        sel <- sapply(object,function(x) x$nzero[x$lambda==x$lambda.min])
        max <- sapply(object,function(x) max(x$nzero))
        two <- paste0("[",min,",(",sel,"),",max,"]")
        three <- sapply(object,function(x) length(x$lambda))
        print(data.frame(lambda=one,nzero=two,length=three))
        ### end trial ###
       
        for(i in seq_along(nzero)){
            for(j in seq_along(model)){
                #if(model[j]=="elastic" & nzero[i]!=10){next} # trial
                if(family=="binomial"){
                    #if(model[j]=="elastic50"){
                    #  temp <- glmnet:::predict.cv.glmnet(object=elastic50,newx=x1,type="response",max=nzero[i])
                    #  # BUG: max is not taken into account!
                    #} else if(model[j]=="elastic95"){
                    #  temp <- glmnet:::predict.cv.glmnet(object=elastic95,newx=x1,type="response",max=nzero[i])
                    #  # BUG: max is not taken into account!
                    #if(model[j]=="elastic"){
                    #  #temp <- palasso:::predict.enet(object=enet,newdata=x1)
                    #  temp <- palasso:::predict_enet(object=enet,newdata=x1,type="response")
                    #} else {
                      temp <- predict.palasso(object=object,
                        newdata=X1,model=model[j],type="response",max=nzero[i])
                    #}
                    pred[i,j][[1]][fold.ext==k] <- temp
                }
                if(family=="cox"){
                    #if(model[j] %in% c("elastic50","elastic95")){
                    #  stop("Elastic net not implemented for Cox model!")
                    #} else {
                        beta <- predict.palasso(object=object,newdata=X1,
                          model=model[j],type="coeff",max=nzero[i])
                    #}
                    newX <- do.call(what="cbind",args=X)
                    plfull <- glmnet::coxnet.deviance(x=newX,y=y,beta=beta)
                    newX0 <- do.call(what="cbind",args=X0)
                    plminusk <- glmnet::coxnet.deviance(x=newX0,y=y0,beta=beta)
                    cvraw[i,j][[1]][k] <- plfull - plminusk
                }
            }
        }
        
        #if(elastic){
        #  if(family=="binomial"){
        #          temp <- palasso:::predict.enet(object=enet,newdata=x1)
        #          pred["10","elastic"][[1]][fold.ext==k] <- temp
        #  }
        # }
        
    }
    
    for(i in seq_along(nzero)){
        for(j in seq_along(model)){
            #if(model[j]=="elastic" & nzero[i]!=10){next} # trial
            if(family=="binomial"){
                y_hat <- pred[i,j][[1]]
                mse[i,j] <- mean((y_hat-y)^2)
                mae[i,j] <- mean(abs(y_hat-y))
                class[i,j] <- mean(abs(round(y_hat)-y))
                y_hat <- pmax(1e-05,pmin(y_hat,1-1e-05))
                deviance[i,j] <- mean(-2*(y*log(y_hat)+(1-y)*log(1-y_hat)))
                auc[i,j] <- pROC::roc(response=y,predictor=y_hat)$auc
            }
            if(family=="cox"){
                weights <- tapply(X=y[,"status"],INDEX=fold.ext,FUN=sum)
                loss[i,j] <- mean(cvraw[i,j][[1]]/mean(weights))
            }
        }
    }
    
    end <- Sys.time()
    
    info <- data.frame(nfolds.ext=nfolds.ext,nfolds.int=nfolds.int,
                       time=format(end-start))
    if(family=="binomial"){
       list <- list(info=info,deviance=deviance,auc=auc,mse=mse,mae=mae,class=class) 
    }
    if(family=="cox"){
        list <- list(info=info,partial.likelihood=loss)
    }
    
    return(list)
}

#' @rdname other
#' @keywords internal
#' 
.select <- function(y,X,index,nfolds=5,standard=TRUE,adaptive=TRUE,...){
    
    fit <- palasso::palasso(y=y,X=X,family="binomial",nfolds=nfolds,
                            standard=standard,...)
    
    names <- unique(c(names(fit),"paired.adaptive","paired.standard","paired.combined"))
    nzero <- c(3,4,5,10,15,20,25,50,Inf)
    
    shots <- hits1 <- hits2 <- matrix(integer(),
                            nrow=length(nzero),
                            ncol=length(names),
                            dimnames=list(nzero,names))
    
    for(i in seq_along(nzero)){
        for(j in seq_along(names)){
            coef <- coef.palasso(fit,model=names[j],max=nzero[i])
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


.design <- function(x){
  
  if(length(x)==1){
    n <- x
    names <- seq_len(x)
  } else {
    n <- length(x)
    names <- x
  }
  
  # check input
  if(n!=round(n)){stop("Provide integer n.")}
  if(n<2){stop("Provide greater n.")}
  if(length(n)!=1){stop("Provide single n.")}
  
  # compute
  if(n%%2==1){
    data <- c(0,rep(seq_len(n),times=n),seq_len(n-1))
    X <- matrix(data=data,nrow=n+1,ncol=n)[-(n+1),]
  }
  if(n%%2==0){
    data <- c(0,rep(1:(n-1),times=n+1))
    X <- matrix(data=data,nrow=n,ncol=n)
    temp <- 2*seq_len(n-1)
    X[-1,n] <- ifelse(temp<n,temp,temp-n+1)
  }
  X[row(X)>=col(X)] <- 0
  rownames(X) <- names
  colnames(X) <- names
  
  # check output
  con <- list()
  con[[1]] <- dim(X)==c(n,n)
  table <- table(X)
  if(n%%2==1){
    con[[2]] <- names(table)==c(0,seq_len(n))
    con[[3]] <- table[-1]==choose(n,2)/n 
  }
  if(n%%2==0){
    con[[2]] <- names(table)==c(0,seq_len(n-1))
    con[[3]] <- table[-1]==choose(n,2)/(n-1)
  }
  con[[4]] <- table(X)[1]==choose(n,2)+n
  con[[5]] <- sapply(seq_len(n),function(x) rowSums(X==x)+colSums(X==x))<=1
  if(any(!unlist(con))){stop("Invalid output!")}
  
  # image(t(X)[,n:1])
  return(X)
}

#' @title
#' Combining p-values
#'
#' @description
#' This function combines local \eqn{p}-values to a global \eqn{p}-value.
#' 
#' @export
#' @keywords methods
#' 
#' @param x local \eqn{p}-values\strong{:}
#' numeric vector of length \eqn{k}
#' 
#' @param method
#' character \eqn{"fisher"}, \eqn{"tippet"}, \eqn{"sidak"}, or \eqn{"simes"}
#' 
#' @return
#' These functions return a numeric vector of length \eqn{p}
#' (main effects),
#' or a numeric matrix with \eqn{p} rows and \eqn{p} columns
#' (interaction effects).
#' 
#' @references
#' Westfall, P. H. (2005). "Combining p-values".
#' Encyclopedia of Biostatistics
#' \doi{10.1002/0470011815.b2a15181}
#'
#' @examples
#' # independence
#' p <- runif(10)
#' palasso:::.combine(p)
#' 
#' ## dependence 
#' #runif <- function(n,cor=0){
#' #    Sigma <- matrix(cor,nrow=n,ncol=n)
#' #     diag(Sigma) <- 1
#' #     mu <- rep(0,times=n)
#' #     q <- MASS::mvrnorm(n=1,mu=mu,Sigma=Sigma)
#' #     stats::pnorm(q=q)
#' #}
#' #p <- runif(n=10,cor=0.8)
#' #combine(p)
#' 
.combine <- function(x,method="simes"){
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  if(any(x>1)){stop("Invalid p-value.")}
  k <- length(x)
  if(method=="fisher"){
    1 - stats::pchisq(q=sum(-2*log(x)),df=2*k)
  } else if(method=="tippet"){
    min(1,min(x)*k)
  } else if(method=="sidak"){
    1 - (1-min(x))^k
  } else if(method=="simes"){
    min(k*sort(x)/(1:k))
  }
}


