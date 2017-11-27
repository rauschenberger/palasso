
testthat::test_that("testthat works",{
    testthat::expect_identical(object=1,expected=1)
})


testthat::test_that("weights are in unit interval",{
    set.seed(4)
    n <- 100; p <- 200
    y <- rbinom(n=n,size=1,prob=0.5)
    X <- lapply(1:2,function(x) matrix(rnorm(n*p),nrow=n,ncol=p))
    reg <- palasso(y=y,X=X,family="binomial")
    testthat::expect_identical(object=all(reg$weights >= 0 & reg$weights<=1),expected=TRUE)
})
