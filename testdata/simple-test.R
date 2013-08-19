######
# testing with mean-centered multivariate gaussians
#  X = A Z  ~ N(0,I)
#  cov(X)[i,j] = sum_k A_ik A_jk = A A^T
#  Y = X - V X = (I-V) A Z
#  cov(Y)[i,j] = (I-V) A A^T (I-V^T)

make.sqrtcov <- function (par) {
    matrix( c(par[1:3],
            0,par[4:5],
            0,0,par[6]), nrow=3 )
}

make.fullmat <- function (par) {
    sqrtcovmat <- make.sqrtcov(par)
    return( tcrossprod(sqrtcovmat) )
}

initpar <- c(1,.5,.25,1,.5,1)
sqrtcovmat <- make.sqrtcov(initpar)
full.covmat <- make.fullmat(initpar)

xx <- t( sqrtcovmat %*% matrix(rnorm(300000),nrow=3) )
var(xx) - full.covmat

yy <- (xx - rowMeans(xx))
sub.covmat <- tcrossprod( (diag(3)-1/3) %*% sqrtcovmat )
var(yy) - sub.covmat

fchol <- chol(sub.covmat,pivot=TRUE)
sub.covmat - crossprod( fchol[,order(attr(fchol,"pivot"))] )


## test mle
xx <- t( sqrtcovmat %*% matrix(rnorm(300000),nrow=3) )
datavec <- t(xx)
llfun <- function (par) {
    fullmat <- make.fullmat( par )
    fchol <- chol(fullmat)
    return( sum( backsolve( fchol, datavec )^2 )/2 + ncol(datavec) * sum(log(diag(fchol))) ) 
}
llfun(initpar)

mlestim1 <- optim( par=initpar, fn=llfun, method="Nelder-Mead", control=list( maxit=500, fnscale=llfun(initpar)/10 ) )
rbind(initpar,mlestim1$par)

## and on projected data
qr.pmat <- qr( (diag(3)-1/3) )
stopifnot( all.equal( t(yy), tcrossprod((diag(3)-1/3), xx) ) )
pmat <- qr.Q(qr.pmat)[,1:2]
stopifnot( all.equal( xx %*% pmat, yy %*% pmat ) )

xx <- t( sqrtcovmat %*% matrix(rnorm(300000),nrow=3) )
yy <- (xx - rowMeans(xx))
datavec <- t(yy)
llfun <- function (par) {
    fullmat <- make.fullmat( par )
    submat <- ( ( crossprod( pmat, fullmat) %*% pmat ) )
    fchol <- chol(submat)
    return( sum( backsolve( fchol, datavec )^2 )/2 + ncol(datavec) * sum(log(diag(fchol))) ) 
}
llfun(initpar)

mlestim1 <- optim( par=initpar, fn=llfun, method="Nelder-Mead", control=list( maxit=500, fnscale=llfun(initpar)/10 ) )
rbind(initpar,mlestim1$par)


###
# ok, really simple
# 1-D:
xx <- rnorm(100000)
llfun <- function(sigma) {
    return( sum( xx^2 ) / (2*sigma^2)  + length(xx)*sum(log(sigma)) )
}
optimize( f=llfun, interval=c(0,10) )

# 2-D, diagonal:
xx <- c(1,2)*matrix(rnorm(200000),nrow=2)
llfun <- function(sigma) {
    return( sum( xx^2/sigma^2 ) / 2  + ncol(xx)*sum(log(sigma)) )
}
optim( par=c(1,2), fn=llfun, method="Nelder-Mead", control=list( maxit=500, fnscale=llfun(c(1,2))/10 ) )

# 2-D, general, variance-covariance parametrized

true.vars <- c(1,2)
true.cors <- c(.2)
true.sigma <- c(1,true.cors,true.cors,1) * outer(sqrt(true.vars),sqrt(true.vars),"*")
true.sqrt.sigma <- t( chol(true.sigma) )  # NOTE TRANSPOSE
stopifnot( all.equal( tcrossprod(true.sqrt.sigma) , true.sigma ) )

xx <- true.sqrt.sigma %*% matrix(rnorm(200000),nrow=2)
llfun <- function(par) {
    vars <- par[1:2]
    cors <- par[3]
    sigma <- c(1,cors,cors,1) * outer(sqrt(vars),sqrt(vars),"*")
    fchol <- chol(sigma)
    return( sum( backsolve(fchol,xx,transpose=TRUE)^2 ) / 2  + ncol(xx)*sum(log(diag(fchol))) )
}
initpar <- c(true.vars,true.cors)
optim( par=initpar, fn=llfun, method="L-BFGS-B", lower=c(0,0,-1)+1e-6, upper=c(Inf,Inf,1-1e-6), control=list( maxit=500, fnscale=llfun(initpar)/10 ) )

# 2-D, general, variance-covariance parametrized
true.sqrt.sigma <- matrix( c(1,.5,.5,sqrt(2)), nrow=2 )
true.sigma <- crossprod(true.sqrt.sigma)

xx <- t(true.sqrt.sigma) %*% matrix(rnorm(200000),nrow=2)
llfun <- function(par) {
    sqrt.sigma <- matrix( par[c(1,2,2,3)], nrow=2 )
    sigma <- crossprod(sqrt.sigma)
    fchol <- chol(sigma)
    return( sum( backsolve(fchol,xx,transpose=TRUE)^2 ) / 2  + ncol(xx)*sum(log(diag(fchol))) )
}
initpar <- true.sqrt.sigma[c(1,3,4)]
optim( par=initpar, fn=llfun, method="L-BFGS-B", control=list( maxit=500, fnscale=llfun(initpar)/10 ) )

###
# check the math
require(MASS)
true.fchol <- chol(true.sigma)
true.siginv <- ginv(true.sigma)
true.fchol.inv <- backsolve(true.fchol,diag(2))
zz <- matrix( rnorm(10),nrow=2 )
all.equal( crossprod(true.fchol), true.sigma )
all.equal( t(zz) %*% t(true.fchol) %*% true.fchol %*% zz, t(zz) %*% true.sigma %*% zz )
all.equal( sum( diag( t(zz) %*% t(true.fchol) %*% true.fchol %*% zz ) ), sum( (true.fchol %*% zz)^2 ) )
all.equal( sum( (true.fchol%*%zz)^2 ), sum( diag(crossprod(zz,true.sigma%*%zz)) ) )
# 
all.equal( tcrossprod(backsolve(true.fchol,diag(2))), true.siginv )
all.equal( t(backsolve(true.fchol,zz,transpose=TRUE)) %*% backsolve(true.fchol,zz,transpose=TRUE), t(zz) %*% true.siginv %*% zz )
all.equal( sum( backsolve(true.fchol,zz,transpose=TRUE)^2 ), sum( diag(crossprod(zz,true.siginv%*%zz)) ) )
#
all.equal( log(sqrt(det(true.sigma))), sum(log(diag(true.fchol))) )
#
all.equal( sum( apply( zz, 2, function (z) sum( backsolve(true.fchol,z,transpose=TRUE)^2 ) / 2  + sum(log(diag(true.fchol))) ) ),
    sum( backsolve(true.fchol,zz,transpose=TRUE)^2 )/2 + ncol(zz) *  sum(log(diag(true.fchol))) )
