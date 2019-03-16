
if(FALSE){

enet <- function(y,x,alpha=NULL,nalpha=21,foldid=NULL,nfolds=10,...){
  if(is.null(alpha)){
    alpha <- seq(from=0,to=1,length.out=nalpha)
  }
  if(is.null(foldid)){
    foldid <- palasso:::.folds(y=y,nfolds=nfolds)
  }
  model <- list()
  for(i in seq_along(alpha)){
    model[[i]] <- glmnet::cv.glmnet(y=y,x=x,foldid=foldid,alpha=alpha[i],...)
    model[[i]]$alpha <- alpha[i]
  }
  
  names(model) <- paste0("alpha",alpha)
  class(model) <- "enet"
  return(model)
}

predict.enet <- function(object,newdata,model="elastic",nalpha=NULL,alpha=NULL){
  if(!is.null(nalpha)&!is.null(alpha)){stop("Invalid.",call.=FALSE)}
  if(!is.null(nalpha)|!is.null(alpha)){
    given <- names(object)
    if(!is.null(nalpha)){
      select <- paste0("alpha",seq(from=0,to=1,length.out=nalpha))
    }
    if(!is.null(alpha)){
      select <- paste0("alpha",alpha)
    }
    if(any(!select %in% given)){stop("Missing alpha.",call.=FALSE)}
    cond <- which(given %in% select)
    object <- object[cond]
  }
  
  if(model=="elastic"){
    cvm <- sapply(object,function(x) x$cvm[x$lambda==x$lambda.min])
    id <- which.min(cvm)
  } else if(model=="ridge"){
    id <- which(names(object)=="alpha0")
  } else if(model=="lasso"){
    id <- which(names(object)=="alpha1")
  } else {
    stop("Invalid model")
  }
  
  #cat(names(object)[id])
  #graphics::plot(x=alpha,y=cvm,main="elastic net"); rm(cvm)
  #graphics::axis(side=1,at=c(0,1),label=c("ridge","lasso"),line=1.5,tick=FALSE)
  
  out <- glmnet::predict.cv.glmnet(object=object[[id]],
                                   newx=newdata,type="response",s="lambda.min")
  return(out)
}


}

##simulate data
if(FALSE){
  n <- 100; p<- 500;
  x <- matrix(rnorm(n*p),n,p)
  x <- t((t(x) - apply(t(x),1,mean))/apply(t(x),1,sd))
  Beta <- c(rep(1,10),rep(-1,10), rep(0,p-20))
  y <- x%*%Beta
  y <- y>median(y)
  net <- enet(y=y,x=x,family="binomial")
  max(net$df[net$lambda>=net$lambda.sel])

  #sum(stats::predict(object=net,newx=x,type="coef",s=net$lambda.sel)!=0)
  
  #sum(predict(net,newdata=x,type="coef")!=0)
  #cbind(net$lambda,net$df)
}


#setwd("C:/Users/arra/Desktop/MATHS/palasso_desktop")

enet <- function(y,x,alpha=0.5,dfmax=10,family="gaussian"){
  if(is.list(x)){x <- cbind(x[[1]],x[[2]])}

  net <- glmnet::glmnet(x=x,y=y,alpha=alpha,family=family)
  #nsel <- length(which(net$df<=dfmax))
  #if(nsel==dfmax){
  #  lambda <- net$lambda[nsel]
  #} else {
  
  nzero <- Matrix::colSums(net$beta!=0)
  if(any(nzero!=net$df)){warning("Invalid.")}
   
  lower <- min(net$lambda[net$df<dfmax])
  upper <- max(net$lambda[net$df>dfmax])
    
  lambda <- seq(from=upper,to=lower,length.out=100)
  net <- glmnet::glmnet(x=x,y=y,alpha=alpha,family=family,lambda=lambda)
    
    #fsel <- function(lambda,dfmax=10,alpha=0.5){
    #  if(lambda==0){
    #    Inf
    #  } else {
    #    net <- glmnet::glmnet(x=x,y=y,nlambda=1,lambda=lambda,alpha=alpha,family=family)
    #    return(net$df-dfmax)
    #  }
    #}
    #lambda <- suppressMessages(stats::uniroot(fsel,interval=c(lower,upper),maxiter=50,
    #                  alpha=alpha,dfmax=dfmax)$root)
  #}

  #verify; few cases it selects a model with e.g. 1 covariate less
  #net <- glmnet::glmnet(x=x,y=y,nlambda=1,lambda=lambda,alpha=alpha,family=family)
  
  nzero <- Matrix::colSums(net$beta!=0)
  if(any(nzero!=net$df)){warning("Unequal.")}
  
  net$lambda.sel <- net$lambda[net$df<=dfmax]
  
  if(any(net$df[net$lambda %in% net$lambda.sel]>dfmax)){warning("Too many.")}
  #cond1 <- 
  #cond2 <- nzero[net$lambda %in% net$lambda.sel]
  
  net$lambda.sel <- stats::median(net$lambda[net$df==dfmax]) # == and median?
  if(is.na(net$lambda.sel)){
    net$lambda.sel <- min(net$lambda[net$df<=dfmax])
  }
  
  #class(net) <- "enetsel"
  
  return(net)
}

predict_enet <- function(object,newdata,type="response"){
  #glmnet::predict.glmnet(object=object,newx=newdata,type=type,s=object$lambda.sel)
  stats::predict(object=object,newx=newdata,type=type,s=object$lambda.sel)
}

