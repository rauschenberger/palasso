

# An extra difficulty is to account for the maximum number of non-zero
# coefficients. If multiple alphas are used, proceed with selection
# as with palasso.
# 
# 
# for(j in seq_along(alpha)){
#   model[[j]] <- glmnet::cv.glmnet(y=y0,x=X0,family="binomial",alpha=alpha[j],
#                                   foldid=foldid.int) # type.measure
# }
# rm(j)
# cvm <- sapply(model,function(x) x$cvm[x$lambda==x$lambda.min])
# id <- which.min(cvm)
# graphics::plot(x=alpha,y=cvm,main="elastic net")
# graphics::axis(side=1,at=c(0,1),label=c("ridge","lasso"),line=1.5,tick=FALSE)
# 
# pred[foldid==i,"lasso"] <- glmnet::predict.cv.glmnet(object=model[[which(alpha==1)]],
#                                                      newx=X1,type="response",s="lambda.min")