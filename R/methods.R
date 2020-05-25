
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
#' @param max
#' maximum number of non-zero coefficients,
#' positive integer,
#' or \code{NULL}
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
#' @seealso
#' Use \link[palasso]{palasso} to fit the paired lasso.
#' 
NULL

#' @rdname methods
#' @export
#' 
predict.palasso <- function(object,newdata,model="paired",s="lambda.min",max=NULL,...){
    x <- subset.palasso(x=object,model=model,max=max)
    newx <- do.call(what="cbind",args=newdata)
    if(is.null(s)){s <- x$glmnet.fit$lambda}
    if(s=="lambda.min"){s <- x$lambda.min}
    stats::predict(object=x$glmnet.fit,newx=newx,s=s,...)
}

#' @rdname methods
#' @export
#' 
coef.palasso <- function(object,model="paired",s="lambda.min",max=NULL,...){
    x <- subset.palasso(x=object,model=model,max=max)
    if(is.null(s)){s <- x$glmnet.fit$lambda}
    if(s=="lambda.min"){s <- x$lambda.min}
    coef <- stats::coef(object=x$glmnet.fit,s=s,...)
    if(rownames(coef)[1]=="(Intercept)"){
        # intercept <- coef[1,]
        coef <- coef[-1,,drop=FALSE]
    }
    .split(x=coef,info=x$palasso)
}

#' @rdname methods
#' @export
#' @importFrom stats weights
#' 
weights.palasso <- function(object,model="paired",max=NULL,...){
    if(length(list(...))!=0){warning("Ignoring argument.",call.=FALSE)}
    x <- subset.palasso(x=object,model=model,max=max)
    weights <- 1/x$glmnet.fit$call$penalty.factor
    .split(x=weights,info=x$palasso)
}

#' @rdname methods
#' @export
#' 
fitted.palasso <- function(object,model="paired",s="lambda.min",max=NULL,...){
    x <- subset.palasso(x=object,model=model,max=max)
    if(x$glmnet.fit$call$family=="cox"){stop("Use \"predict\" for Cox regression.",call.=FALSE)}
    newx <- x$glmnet.fit$call$x
    if(is.null(s)){s <- x$glmnet.fit$lambda}
    if(s=="lambda.min"){s <- x$lambda.min}
    stats::predict(object=x$glmnet.fit,newx=newx,s=s,type="response",...)
}

#' @rdname methods
#' @export
#' 
residuals.palasso <- function(object,model="paired",s="lambda.min",max=NULL,...){
    x <- subset.palasso(x=object,model=model,max=max)
    if(x$glmnet.fit$call$family=="cox"){stop("Use \"predict\" for Cox regression.",call.=FALSE)}
    newx <- x$glmnet.fit$call$x
    if(is.null(s)){s <- x$glmnet.fit$lambda}
    if(s=="lambda.min"){s <- x$lambda.min}
    y <- x$glmnet.fit$call$y
    y_hat <- stats::predict(object=x$glmnet.fit,newx=newx,s=s,type="response",...)
    y - y_hat
}

#' @rdname methods
#' @export
#' 
deviance.palasso <- function(object,model="paired",max=NULL,...){
    x <- subset.palasso(x=object,model=model,max=max)
    stats::deviance(x$glmnet.fit,...)
}

#' @rdname methods
#' @export
#' 
logLik.palasso <- function(object,model="paired",max=NULL,...){
    if(length(list(...))!=0){warning("Ignoring argument.",call.=FALSE)}
    x <- subset.palasso(x=object,model=model,max=max)$glmnet.fit
    if(x$call$family=="cox"){
        cox <- survival::coxph(x$call$y~1,weights=x$call$weights)
        ll0 <- cox$loglik # survival:::logLik.coxph.null(cox)
    } else {
        glm <- stats::glm(x$call$y~1,weights=x$call$weights,family=x$call$family)
        ll0 <- stats::logLik(glm)
    }
    ll1 <- x$nulldev/2 + ll0 - stats::deviance(x)/2
    attributes(ll1)$df <- df.residual.glmnet(x)
    attributes(ll1)$nobs <- x$nobs
    class(ll1) <- c("logLik.palasso","logLik")
    return(ll1)
}

#' @rdname methods
#' @export
#' 
summary.palasso <- function(object,model="paired",...){
    if(length(list(...))!=0){warning("Ignoring argument.",call.=FALSE)}
    
    # header
    title <- paste(object[[1]]$glmnet.fit$call$family,"palasso")
    line <- paste(rep("-",times=nchar(title)),collapse="")
    cat("",line,"\n",title,"\n",line,"\n\n")
    
    # dimensions
    print.palasso(object)
    cat("\n")
    
    # non-zero weights
    weights <- weights.palasso(object,model=model)
    name <- colnames(weights)
    number <- colSums(weights!=0)
    cat("non-zero weights:",paste(number,name,collapse=", "),"\n\n")
    
    # cross-validation
    x <- subset.palasso(x=object,model=model)
    id <- list()
    id$min <- which(x$lambda==x$lambda.min)
    #id$ose <- which(x$lambda==x$lambda.1se)
    frame <- matrix(NA,nrow=length(id),ncol=3)
    frame[,1] <- vapply(id,function(i) x$lambda[i],numeric(1))
    frame[,2] <- vapply(id,function(i) x$nzero[i],integer(1))
    frame[,3] <- vapply(id,function(i) x$cvm[i],numeric(1))
    if(length(id)>1){rownames(frame) <- names(id)} #c("min","1se")
    colnames(frame) <- c("lambda","nzero",x$name) # was names(x$name)
    base::print(round(frame,digits=2))
    return(invisible(NULL))
}

#' @export
#' 
print.palasso <- function(x,...){
    if(length(list(...))!=0){warning("Ignoring argument.",call.=FALSE)}
    info <- attributes(x)$info
    cat("palasso object: ")
    cat(info$n,"samples, ")
    cat(paste0(info$k,"*",info$p),"covariates\n")
    if(length(info$call)>0){
        cond <- vapply(X=info$call,FUN=function(x) length(x)>1,FUN.VALUE=logical(1))
        info$call[cond] <- "..."
        call <- vapply(X=info$call,FUN=deparse,FUN.VALUE=character(1))
        call <- paste0(names(call),"=",call)
        call <- paste(call,collapse=", ")
        call <- paste0("(",call,")")
        call <- gsub(x=call,pattern="\"\\.\\.\\.\"",replacement="...")
        cat(call)
    }
    return(invisible(NULL))
}

#' @export
#' 
subset.palasso <- function(x,model="paired",max=NULL,...){
    
    if(length(list(...))!=0){warning("Ignoring argument.",call.=FALSE)}
    
    if(!inherits(x=x,what="palasso")){
        warning("Fake palasso object?",call.=FALSE)
    }
    
    if(!model %in% c(names(x),"elastic",paste0("paired",c("",".adaptive",".standard",".combined",".adaptive1",".standard1")))){
        stop("Invalid argument \"model\".",call.=FALSE)
    }
    
    name <- unique(vapply(X=x,FUN=function(x) x$name,FUN.VALUE=character(1)))
    if(length(name)!=1){
        stop("Different loss functions!",call.=FALSE)
    }
    
    if(is.null(max)){
        max <- attributes(x)$info$max
    }
    
    if(!is.null(max)){
        for(i in seq_along(x)){
            cond <- x[[i]]$nzero<=max
            if(length(cond)==0){stop("Adapt lambda sequence!",call.=FALSE)}
            for(j in c("lambda","cvm","cvsd","cvup","cvlo","nzero")){
                x[[i]][[j]] <- x[[i]][[j]][cond] 
            }
            if(name=="AUC"){
                x[[i]]$lambda.min <- x[[i]]$lambda[which.max(x[[i]]$cvm)]
                #x[[i]]$lambda.1se <- max(x[[i]]$lambda[x[[i]]$cvm>=max(x[[i]]$cvlo[which.max(x[[i]]$cvm)])])
            } else {
                x[[i]]$lambda.min <- x[[i]]$lambda[which.min(x[[i]]$cvm)]
                #x[[i]]$lambda.1se <- max(x[[i]]$lambda[x[[i]]$cvm<=min(x[[i]]$cvup[which.min(x[[i]]$cvm)])])
            }
            cond <- x[[i]]$glmnet.fit$df<=max
            for(j in c("a0","df","lambda","dev.ratio")){
                x[[i]]$glmnet.fit[[j]] <- x[[i]]$glmnet.fit[[j]][cond]
            }
            x[[i]]$glmnet.fit$beta <- x[[i]]$glmnet.fit$beta[,cond,drop=FALSE]
        }
    }
    
    if(model=="paired"){
        adaptive <- attributes(x)$info$adaptive
        standard <- attributes(x)$info$standard
        if(adaptive & !standard){
            model <- "paired.adaptive"
        } else if(!adaptive & standard){
            model <- "paired.standard"
        } else if(adaptive & standard){
            model <- "paired.combined" # original
            warning("Consider model=\"paired.adaptive\" or model=\"paired.standard\".",call.=FALSE)
        }
    }
    
    if(model=="paired.adaptive"){
        pattern <- "adaptive|within"
        cond <- grepl(pattern=pattern,x=names(x))
        if(sum(cond)!=attributes(x)$info$k+2){stop("Mismatch.")}
    } else if(model=="paired.standard"){
        pattern <- "standard|between"
        cond <- grepl(pattern=pattern,x=names(x))
        if(sum(cond)!=attributes(x)$info$k+2){stop("Mismatch.")}
    } else if(model=="paired.combined"){
        pattern <- "standard|adaptive|between|within"
        cond <- grepl(pattern=pattern,x=names(x))
        if(sum(cond)!=2*attributes(x)$info$k+4){stop("Mismatch.")}
    } else if(model=="elastic"){
        pattern <- "elastic"
        cond <- grepl(pattern=pattern,x=names(x))
        if(sum(cond)!=1){stop("Mismatch.")} # was 4
    } else {
        cond <- names(x)==model # important
    }
   
    object <- x[cond]
    if(name=="AUC"){
        loss <- vapply(X=object,FUN=function(x) max(x$cvm),FUN.VALUE=numeric(1)) # trial: na.rm=TRUE
        select <- which.max(loss)
    } else {
        loss <- vapply(X=object,FUN=function(x) min(x$cvm),FUN.VALUE=numeric(1)) # trial: na.rm=TRUE
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
    if(length(list(...))!=0){warning("Ignoring argument.",call.=FALSE)}
    if(object$call$alpha==1){
        # df <- Matrix::colSums(object$beta!=0)
        df <- Matrix::colSums(stats::coef(object=object)!=0)
    } else {
        d <- svd(object$call$x)$d^2
        df <- sum(d^2/(d^2+object$lambda))
    }
    return(df)
}

#' @export
#' 
print.logLik.palasso <- function(x,...){
    if(length(list(...))!=0){warning("Ignoring argument.",call.=FALSE)}
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
    split
}
