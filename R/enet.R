

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








# #simulate data
# n<- 100; p<- 500;
# X <- matrix(rnorm(n*p),n,p)
# X <- t((t(X) - apply(t(X),1,mean))/apply(t(X),1,sd))
# Beta <- c(rep(1,10),rep(-1,10), rep(0,p-20))
# Y <- X%*%Beta
# 
# net <- enet(Y,X)
# 
# setwd("C:/Users/arra/Desktop/MATHS/palasso_desktop")
# 
# enet <- function(Y,X,alpha0=0.5,maxsel0=10,fam="gaussian"){
#   if(is.list(X)){X <- cbind(x[[1]],x[[2]])}
#   
#   #code
#   library(glmnet)
#   nfeat <- ncol(X)
#   penselEN0 <- glmnet(X,Y,alpha=alpha0,family=fam,standardize=FALSE)
#   nsel = length(which(penselEN0$df<=maxsel0))
#   if(nsel==maxsel0){
#     lam <- penselEN0$lambda[nsel]
#   } else {
#     whlamlow <- penselEN0$lambda[length(which(penselEN0$df <= maxsel0))+1]
#     whlamup <- penselEN0$lambda[1]
#     
#     fsel <- function(lam,maxsel=10,alpha=0.5){
#       if(lam==0) return(nfeat) else {
#         penselEN <- glmnet(X,Y,nlambda=1,lambda=lam,alpha=alpha,family=fam,standardize=FALSE)
#         coef <- penselEN$beta
#         return(length(coef[coef!=0])-maxsel)
#       }
#     }
#     lam <- uniroot(fsel,interval=c(whlamlow,whlamup),maxiter=50,alpha=alpha0,maxsel=maxsel0)$root
#   }
#   
#   #verify; few cases it selects a model with e.g. 1 covariate less
#   penselEN2 <- glmnet(X,Y,nlambda=1,lambda=lam,alpha=alpha0,family=fam,standardize=FALSE)
#   return(penselEN2)
# }
# 
# penselEN2$df
