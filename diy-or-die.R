###PETER: just change this path to wherever you put 42_rib.FEMALES.combined.residuals and consensusTree_9species.txt
#everything should run fine after that.
path="./"
data_file="42_rib.FEMALES.combined.residuals"
tree_file="consensusTree_9species.txt"

require(MASS)
library(ape)
library(nlme)

mydata<-read.table(file=paste(path, data_file, sep=""), header=TRUE)
mytree<-read.nexus(file=paste(path, tree_file, sep=""))

color_choices <- rainbow(nlevels(mydata$species))
with(mydata, plot(residual_testis, residual_rib_volume, col=color_choices[species], main="RIB, FEMALES, COMBINED" ) ) 
legend("topleft", levels(mydata$species), pch=1, col=color_choices, cex=0.5, bty="n")

within_length <- .01

newdata <- mydata
rownames(newdata) <- make.names(newdata$species,unique=TRUE)

newtree <- mytree
for (sp in levels(newdata$species)[table(newdata$species)>1]) {
    theseones <- (newdata$species==sp)
    nsamps <- sum(theseones)
    tmptree <- list(
        edge = cbind( rep(nsamps+1,nsamps), 1:nsamps ),
        edge.length = rep(within_length,nsamps),
        Nnode = 1,
        tip.label = rownames(newdata)[theseones],
        root.edge = within_length/10
        )
    class(tmptree) <- "phylo"
    newtree <- bind.tree(x=newtree,y=tmptree,where=match(sp,newtree$tip.label))
}

plot(newtree)

stopifnot( setequal( rownames(newdata), newtree$tip.label ) )


# fit model where y ~ a x + b + eps
#  where eps is N(0,Sigma)
#  and Sigma is specified between species by treecor times sigma2
downweight <- TRUE
pendant.edges <- (newtree$edge[,2] <= nrow(newdata)) & (newtree$edge.length == within_length)
d <- nrow(mydata)
sqrt_nsamples <- sqrt(table(mydata$species)[mydata$species])
loglik <- function (params) {
    a <- params[1]  # slope 
    b <- params[2]  # intercept
    sigma <- params[3]  # variance
    newtree$edge.length[pendant.edges] <- params[4]  # pendant edge length
    new_corstr <- corBrownian(phy=newtree)
    cor_matrix <- corMatrix(Initialize(new_corstr,newdata))
    if (downweight) { cor_matrix <- cor_matrix / outer(sqrt_nsamples,sqrt_nsamples,"*") }
    cov_chol <- chol(cor_matrix*sigma,pivot=TRUE)
    x <- mydata$residual_rib_volume - b - (a * mydata$residual_testis)
    z <- sum( (solve(cov_chol) %*% x )^2 )/ (2*sigma)
    logdet <- log(prod(diag(cov_chol)^2))
    ans <- z + (d/2)*logdet
    # if (!is.numeric(ans) | is.na(ans) | !is.finite(ans) ) { browser() }
    return( ans )
}

initvals <- c(0,-0.1,var(mydata$residual_rib_volume),within_length)
ans <- optim( par=initvals, fn=loglik, lower=c(-Inf,-Inf,0,within_length/10), method="L-BFGS-B", control=list(parscale=c(.01,.1,.01,.01)) )

# ok, plot

# with(newdata, points( tapply(residual_testis,species,mean), tapply(residual_rib_volume,species,mean), pch=20, col=color_choices, cex=2 ) )
abline(ans$par[2],ans$par[1])
mtext(paste("(V)=", round(ans$par[1], digits=8), "*(L)+", round(ans$par[2], digits=8), sep=""))

