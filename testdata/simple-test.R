######
# testing with mean-centered multivariate gaussians
#  X = A Z  ~ N(0,I)
#  cov(X)[i,j] = sum_k A_ik A_jk = A A^T
#  Y = X - V X = (I-V) A Z
#  cov(Y)[i,j] = (I-V) A A^T (I-V^T)

sqrtcovmat <- matrix( c(1,.5,.25,
                    0,1,.5,
                    0,0,1), nrow=3 )

full.covmat <- tcrossprod(sqrtcovmat)

xx <- t( sqrtcovmat %*% matrix(rnorm(300000),nrow=3) )
var(xx) - full.covmat

yy <- (xx - rowMeans(xx))
sub.covmat <- tcrossprod( (diag(3)-1/3) %*% sqrtcovmat )
var(yy) - sub.covmat

fchol <- chol(sub.covmat,pivot=TRUE)
sub.covmat - crossprod( fchol[,order(attr(fchol,"pivot"))] )
