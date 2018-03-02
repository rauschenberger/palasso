
#--- Generic functions ---------------------------------------------------------

#' @name methods
#' 
#' @title
#' Methods for class "palasso"
#' 
#' @description
#' This page lists the main methods for class "palasso".
#' 
#' @param object
#' \link[palasso]{palasso} object
#' 
#' @param newdata
#' covariates\strong{:}
#' list of matrices, each with \eqn{n} rows (samples)
#' and \eqn{p} columns (variables)
#' 
#' @param s penalty parameter\strong{:}
#' character \code{"lambda.min"} or \code{"lambda.1se"},
#' positive numeric,
#' or \code{NULL} (entire sequence)
#' 
#' @param model
#' character \code{"paired"},
#' or an entry of \code{names(object)}
#' 
#' @param ...
#' further arguments for
#' \code{\link[glmnet]{predict.cv.glmnet}},
#' \code{\link[glmnet]{coef.cv.glmnet}},
#' or \code{\link[glmnet]{deviance.glmnet}}
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
NULL

#' @rdname methods
#' @export
#' 
predict.palasso <- function(object,newdata,model="paired",s="lambda.min",...){
    if(missing(newdata)||is.null(newdata)) {
        stop("Fitted values?")
    }
    x <- palasso:::subset.palasso(x=object,model=model)
    newx <- do.call(what="cbind",args=newdata)
    if(is.null(s)){s <- x$glmnet.fit$lambda}
    glmnet::predict.cv.glmnet(object=x,newx=newx,s=s,...)
}

#' @rdname methods
#' @export
#' 
coef.palasso <- function(object,model="paired",s="lambda.min",...){
    x <- palasso:::subset.palasso(x=object,model=model)
    if(is.null(s)){s <- x$glmnet.fit$lambda}
    coef <- glmnet::coef.cv.glmnet(object=x,s=s,...)
    if(rownames(coef)[1]=="(Intercept)"){
        # intercept <- coef[1,]
        coef <- coef[-1,,drop=FALSE]
    }
    palasso:::.split(x=coef,info=x$palasso)
}
# CONTINUE HERE: Maybe it is better to return a matrix of coefficients
# (including intercept). Then rows are covariates and columns are lambdas.
# Use x and z for rownames. However, indices are less meaningful.

#' @rdname methods
#' @export
#' @importFrom stats weights
#' 
weights.palasso <- function(object,model="paired",...){
    if(length(list(...))!=0){warning("Ignoring argument.")}
    x <- palasso:::subset.palasso(x=object,model=model)
    weights <- 1/x$glmnet.fit$call$penalty.factor
    palasso:::.split(x=weights,info=x$palasso)
}

#' @rdname methods
#' @export
#' 
fitted.palasso <- function(object,model="paired",s="lambda.min",...){
    x <- palasso:::subset.palasso(x=object,model=model)
    if(x$glmnet.fit$call$family=="cox"){stop("Use \"predict\" for Cox regression.")}
    newx <- x$glmnet.fit$call$x
    if(is.null(s)){s <- x$glmnet.fit$lambda}
    glmnet::predict.cv.glmnet(object=x,newx=newx,s=s,type="response",...)
}

#' @rdname methods
#' @export
#' 
residuals.palasso <- function(object,model="paired",s="lambda.min",...){
    x <- palasso:::subset.palasso(x=object,model=model)
    if(x$glmnet.fit$call$family=="cox"){stop("Use \"predict\" for Cox regression.")}
    newx <- x$glmnet.fit$call$x
    if(is.null(s)){s <- x$glmnet.fit$lambda}
    y <- x$glmnet.fit$call$y
    y_hat <- glmnet::predict.cv.glmnet(object=x,newx=newx,s=s,type="response",...)
    y - y_hat
}

#' @rdname methods
#' @export
#' 
deviance.palasso <- function(object,model="paired",...){
    x <- palasso:::subset.palasso(x=object,model=model)
    glmnet::deviance.glmnet(x$glmnet.fit,...)
}

#' @rdname methods
#' @export
#' 
logLik.palasso <- function(object,model="paired",...){
    if(length(list(...))!=0){warning("Ignoring argument.")}
    x <- palasso:::subset.palasso(x=object,model=model)$glmnet.fit
    if(x$call$family=="cox"){
        cox <- survival::coxph(x$call$y~1,weights=x$call$weights)
        ll0 <- cox$loglik # survival:::logLik.coxph.null(cox)
    } else {
        glm <- stats::glm(x$call$y~1,weights=x$call$weights,family=x$call$family)
        ll0 <- stats::logLik(glm)
    }
    ll1 <- x$nulldev/2 + ll0 - glmnet::deviance.glmnet(x)/2
    attributes(ll1)$df <- palasso:::df.residual.glmnet(x)
    attributes(ll1)$nobs <- x$nobs
    class(ll1) <- c("logLik.palasso","logLik")
    return(ll1)
}

# Consider modifying function "logLik.palasso"
# by adding the logical argument "glmnet". 
# If TRUE then calculate logLik from deviance,
# if FALSE from fitted values.
# logLik.palasso <- function(object,s=NULL,model="paired",...){
# 
#     x <- palasso:::subset.palasso(x=object,model=model)
#     if(is.null(s)){s <- x$glmnet.fit$lambda}
#     y <- x$glmnet.fit$call$y
#     newx <- x$glmnet.fit$call$x
#     family <- x$glmnet.fit$call$family
# 
#     eta <- glmnet::predict.cv.glmnet(object=x,newx=newx,s=s,type="link",...)
#     if(is.vector(eta)){eta <- as.matrix(eta)}
# 
#     ll <- rep(x=NA,times=length(s))
#     for(i in seq_along(ll)){
#         if(family=="gaussian"){
#             mu <- eta[,i]
#             sd <- sqrt(sum((y-mu)^2)/length(y))
#             ll[i] <- sum(log(1/sqrt(2*pi*sd^2)*exp(-(y-mu)^2/(2*sd^2))))
#         } else if(family=="binomial"){
#             p <- exp(eta[,i])/(1+exp(eta[,i]))
#             ll[i] <- sum(log(p^y*(1-p)^(1-y)))
#         } else if(family=="poisson"){
#             lambda <- exp(eta[,i])
#             ll[i] <- sum(log(lambda^y*exp(-lambda)/factorial(y)))
#         }
#     }
#     return(ll)
# }

#' @rdname methods
#' @export
#' 
summary.palasso <- function(object,model="paired",...){
    if(length(list(...))!=0){warning("Ignoring argument.")}
    
    # header
    title <- paste(object[[1]]$glmnet.fit$call$family,"palasso")
    line <- paste(rep("-",times=nchar(title)),collapse="")
    cat("",line,"\n",title,"\n",line,"\n\n")
    
    # dimensions
    palasso:::print.palasso(object)
    cat("\n")
    
    # non-zero weights
    weights <- palasso:::weights.palasso(object)
    name <- colnames(weights)
    number <- colSums(weights!=0)
    cat("non-zero weights:",paste(number,name,collapse=", "),"\n\n")
    
    # cross-validation
    x <- palasso:::subset.palasso(x=object,model=model)
    id <- list()
    id$min <- which(x$lambda==x$lambda.min)
    id$ose <- which(x$lambda==x$lambda.1se)
    frame <- matrix(NA,nrow=2,ncol=3)
    frame[,1] <- sapply(id,function(i) x$lambda[i])
    frame[,2] <- sapply(id,function(i) x$nzero[i])
    frame[,3] <- sapply(id,function(i) x$cvm[i])
    rownames(frame) <- c("min","1se")
    colnames(frame) <- c("lambda","nzero",names(x$name))
    base::print(round(frame,digits=2))
    return(invisible(NULL))
}

#' @export
#' 
print.palasso <- function(x,...){
    if(length(list(...))!=0){warning("Ignoring argument.")}
    info <- attributes(x)$info
    cat("palasso object: ")
    cat(info$n,"samples, ")
    cat(paste0(info$k,"*",info$p),"covariates\n")
    if(length(info$call)>0){
        cat("(",paste(names(info$call),info$call,sep="=",collapse=", "),")\n",sep="")
    }
    return(invisible(NULL))
}

#' @export
#' 
subset.palasso <- function(x,model="paired",...){
    if(length(list(...))!=0){warning("Ignoring argument.")}
    
    if(!inherits(x=x,what="palasso")){
        warning("Fake palasso object?")
    }
    
    name <- unique(sapply(X=x,FUN=function(x) x$name))
    if(length(name)!=1){
        stop("Different loss functions!")
    }
    
    if(model=="paired"){
        pattern <- "adaptive|between|within"
        cond <- grepl(pattern=pattern,x=names(x))
    } else {
        cond <- names(x)==model
    }
    
    object <- x[cond]
    if(name=="AUC"){
        loss <- sapply(object,function(x) max(x$cvm))
        select <- which.max(loss)
    } else {
        loss <- sapply(object,function(x) min(x$cvm))
        select <- which.min(loss)
    }
    object <- object[[select]]
    object$glmnet.fit$call$x <- x[[1]]$glmnet.fit$call$x
    object$palasso <- attributes(x)$info
    object$palasso$select <- names(select)
    
    return(object)
}

#' @export
#' @importFrom stats df.residual
#' 
df.residual.glmnet <- function(object,...){
    if(length(list(...))!=0){warning("Ignoring argument.")}
    if(object$call$alpha==1){
        # df <- Matrix::colSums(object$beta!=0)
        df <- Matrix::colSums(glmnet::coef.glmnet(object=object)!=0)
    } else {
        d <- svd(object$call$x)$d^2
        df <- sum(d^2/(d^2+object$lambda))
    }
    return(df)
}

#' @export
#' 
print.logLik.palasso <- function(x,...){
    if(length(list(...))!=0){warning("Ignoring argument.")}
    X <- rbind(x,attributes(x)$df)
    rownames(X) <- c("log Lik.","eff. df")
    print(X)
}

.split <- function(x,info){
    k <- info$k; p <- info$p
    if(is.vector(x)){x <- as.matrix(x)}
    split <- sapply(seq_len(k),function(i) x[seq(from=(i-1)*p+1,to=i*p,by=1),,drop=FALSE])
    if(is.list(split)){names(split) <- info$names}
    if(is.matrix(split)){colnames(split) <- info$names}
    # as.data.frame(split)
    split
}
