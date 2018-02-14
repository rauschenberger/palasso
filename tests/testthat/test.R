
#--- initialisation ------------------------------------------------------------

set.seed(1)
n <- 100
p <- 200
k <- 2
y <- stats::rbinom(n=n,size=1,prob=0.5)
X <- lapply(seq_len(2),function(x) matrix(stats::rnorm(n*p),nrow=n,ncol=p))
pmax <- 10
fit <- palasso(y=y,X=X,family="binomial",pmax=pmax)

names <- c(names(fit),"paired")

weights <- sapply(X=names,FUN=function(x) weights(object=fit,model=x))
coef <- sapply(X=names,FUN=function(x) coef(object=fit,model=x)[-1])
deviance <- sapply(X=names,FUN=function(x) deviance(object=fit,model=x))
fitted <- sapply(X=names,FUN=function(x) fitted(object=fit,model=x))
predict <- sapply(X=names,FUN=function(x) predict(object=fit,model=x,newdata=X,type="response"))

#--- unit tests ----------------------------------------------------------------

testthat::test_that("testthat works",{
    testthat::expect_true(TRUE)
})

testthat::test_that("weights are large",{
    x <- all(weights >= 0)
    testthat::expect_true(x)
})

testthat::test_that("weights are small",{
    x <- all(weights <= 1)
    testthat::expect_true(x)
})

testthat::test_that("pmax is effective",{
    x <- all(colSums(coef!=0)<=pmax)
    testthat::expect_true(x)
})

testthat::test_that("deviance decreases",{
    x <- all(sapply(deviance,function(x) all(diff(x)<0)))
    testthat::expect_true(x)
})

testthat::test_that("fitted equals predict",{
    testthat::expect_identical(object=fitted,expected=predict)
})

testthat::test_that("fitted equals predict",{
    cond <- grepl(x=names,pattern="standard|between|within")
    group <- rep(seq_len(k),each=p)
    pair <- rep(seq_len(p),times=k)
    sum <- sapply(seq_len(p),function(x) colSums(weights[pair==x,cond]))
    x <- all(sum>1-1e-06 & sum<1+1e-06)
    testthat::expect_true(x)
})
