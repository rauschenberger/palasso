
#--- initialisation ------------------------------------------------------------

set.seed(1)

n <- 100
p <- 200
k <- 2
pmax <- 10
family <- "gaussian"

if(family=="gaussian"){
    y <- stats::rnorm(n=n)
    mu <- mean(y)
    sd <- sqrt(sum((y-mu)^2)/n)
    y <- (y-mu)/sd
}
if(family=="binomial"){
    y <- 1*(stats::rbinom(n=n,size=1,prob=0.5)>=0.5)
}
if(family=="poisson"){
    y <- stats::rpois(n=n,lambda=5)
}
X <- lapply(seq_len(k),function(x) matrix(stats::rnorm(n*p),nrow=n,ncol=p))
X <- lapply(X,function(x) scale(x))

fit <- palasso::palasso(y=y,X=X,family=family,pmax=pmax)

names <- c(names(fit),"paired")
weights <- lapply(X=names,FUN=function(x) weights(object=fit,model=x))
coef <- lapply(X=names,FUN=function(x) coef(object=fit,model=x))
deviance <- lapply(X=names,FUN=function(x) deviance(object=fit,model=x))
logLik <- lapply(X=names,FUN=function(x) logLik(object=fit,model=x))
fitted <- sapply(X=names,FUN=function(x) fitted(object=fit,model=x))
predict <- sapply(X=names,FUN=function(x) predict(object=fit,model=x,newdata=X,type="response"))
residuals <- sapply(X=names,FUN=function(x) residuals(object=fit,model=x))

#--- unit tests ----------------------------------------------------------------

testthat::test_that("testthat works",{
    testthat::expect_true(TRUE)
})

testthat::test_that("weights are large",{
    x <- all(sapply(weights,function(x) all(x>=0)))
    testthat::expect_true(x)
})

testthat::test_that("weights are small",{
    x <- all(sapply(weights,function(x) all(x<=1)))
    testthat::expect_true(x)
})

testthat::test_that("pmax is effective",{
    x <- all(sapply(coef,function(x) sum(x$x!=0)+sum(x$z!=0))<=pmax)
    testthat::expect_true(x)
})

testthat::test_that("deviance decreases",{
    x <- all(sapply(deviance,function(x) all(diff(x)<0)))
    testthat::expect_true(x)
})

testthat::test_that("logLik increaes",{
    x <- all(sapply(logLik,function(x) all(diff(x)>0)))
    testthat::expect_true(x)
})

testthat::test_that("deviance and logLik are perfectly correlated",{
    diff <- 1+sapply(seq_along(names),function(i) cor(deviance[[i]],logLik[[i]],method="spearman"))
    x <- all(abs(diff)<1e-06)
    testthat::expect_true(x)
})

testthat::test_that("fitted equals predict",{
    testthat::expect_identical(object=fitted,expected=predict)
})

testthat::test_that("fitted plus residuals equals observed",{
    diff <- (fitted+residuals)-y
    x <- all(abs(diff)<1e-06)
    testthat::expect_true(x)
})

testthat::test_that("weights sum to one",{
    cond <- grepl(x=names,pattern="standard|between|within")
    diff <- 1-sapply(weights[cond],rowSums)
    x <- all(abs(diff)<1e-06)
    testthat::expect_true(x)
})

# low dimensionality
X <- lapply(X,function(x) x[,seq_len(n/(5*k))]) 
fit <- palasso::palasso(y=y,X=X,lambda=c(99e99,0),family=family)
glm0 <- stats::glm(y~1,family=family)
glm1 <- stats::glm(y~do.call(what="cbind",args=X),family=family)

testthat::test_that("deviance stats",{
    diff <- deviance(fit,model="adaptive_xz")-c(deviance(glm0),deviance(glm1))
    x <- all(abs(diff)<1e-06)
    testthat::expect_true(x)
})

testthat::test_that("logLik stats",{
    diff <- logLik(fit,model="adaptive_xz")-c(logLik(glm0),logLik(glm1))
    if(family=="gaussian"){
        x <- abs(diff[1])<1e-06 & abs(diff[2])<2
    } else {
        x <- all(abs(diff)<1e-06)
    }
    testthat::expect_true(x)
})
