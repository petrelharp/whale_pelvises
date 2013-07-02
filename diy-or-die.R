#PETER - there is no more variable called path.  Write in full paths to data_file and tree_file below.  
rm(list=ls())
#this is just for testing... 
#males, pelvics: 
argument_list<-list(
    path="./",
    outpath="./",
    data_file="30_pelvic_volume.FEMALES.Rin",
    tree_file="consensusTree.pelvic.FEMALES.tree",
    x_variable="body_length",
    y_variable="pelvic_volume",
    logify="TRUE",
    main_title="MALES,pelvic",
    outfile="47.males"
)
argument_list[c("data_file","tree_file")] <- paste(argument_list$path,argument_list[c("data_file","tree_file")],sep='')
argument_list['outfile'] <- paste(argument_list$outpath,argument_list['outfile'],sep='')

#then parse them... 
attach(argument_list)

require(MASS)
library(ape)
library(nlme)

mydata<-read.table(file=data_file, header=TRUE)
mytree<-read.nexus(file=tree_file)

################################################################
#for grabbing (and plotting) some initial values... 
if(logify){
	simple.lm <- with(mydata, lm( log10(mydata[,y_variable]) ~ log10(mydata[,x_variable]) ) )
	} else {
	simple.lm <- with(mydata, lm( (mydata[,y_variable]) ~ (mydata[,x_variable]) ) )
	}
color_choices <- rainbow(nlevels(mydata$species))
if(logify){
	with(mydata, plot(mydata[,x_variable], mydata[,y_variable], log='xy', col=color_choices[species], main=main_title ) ) 
	} else {
	with(mydata, plot(mydata[,x_variable], mydata[,y_variable], col=color_choices[species], main=main_title ) ) 
	}
legend("topleft", levels(mydata$species), pch=1, col=color_choices, cex=0.5, bty="n")
abline(coef(simple.lm)[1],coef(simple.lm)[2],lty=2,col='blue')

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
x <- if (logify) { log10(newdata[,x_variable]) } else { newdata[,x_variable] }
y <- if (logify) { log10(newdata[,y_variable]) } else { newdata[,y_variable] }

sigma_guess <- sqrt(var(resid(simple.lm))) 
pendant_guess <- within_length

downweight <- TRUE
pendant_edges <- (newtree$edge[,2] <= nrow(newdata)) 
orig_lengths <- newtree$edge.length[pendant_edges]
orig_lengths[orig_lengths==within_length] <- 0
d <- nrow(newdata)
sqrt_nsamples <- sqrt(table(newdata$species)[newdata$species])
loglik <- function (params) {
    a <- params[1]  # slope 
    b <- params[2]  # intercept
    sigma <- params[3]  # SD
    newtree$edge.length[pendant_edges] <- orig_lengths + params[4]  # pendant edge length
    new_corstr <- corBrownian(phy=newtree)
    cor_matrix <- corMatrix(Initialize(new_corstr,newdata))
    if (downweight) { cor_matrix <- cor_matrix / outer(sqrt_nsamples,sqrt_nsamples,"*") }
    cov_chol <- chol(cor_matrix,pivot=TRUE)*sigma
    resids <- (y - b - (a * x))
    z <- sum( (solve(cov_chol) %*% resids )^2 )/ 2
    logdet <- sum(log(diag(cov_chol)^2))
    ans <- z + (d/2)*logdet
    # if (!is.numeric(ans) | is.na(ans) | !is.finite(ans) ) { browser() }  ## Uncomment this for debugging.
    return( ans )
}

# do grid search to pick good starting points for variance parameters
sigma_fac <- .9 * sigma_guess
pendant_fac <- .9 * pendant_guess
for (k in 1:3) {
    vargrid <- expand.grid( sigma=sigma_guess+seq(-sigma_fac,3*sigma_fac,length.out=20), pendant=pendant_guess+seq(-pendant_fac,20*pendant_fac,length.out=20) )
    vargrid$loglik <- apply( vargrid,1,function (x) loglik( c(rev(coef(simple.lm)),x) ) )
    if (interactive()) { with( vargrid, plot( sigma, pendant, cex=3/(1+loglik-min(loglik)) ) ) }
    best_guess <- which.min( vargrid$loglik )
    sigma_guess <- vargrid$sigma[best_guess]
    pendant_guess <- vargrid$pendant[best_guess]
    if ( sigma_guess < max(vargrid$sigma) & sigma_guess > min(vargrid$sigma) ) { 
        sigma_fac <- min( sigma_fac/5, .9*sigma_guess ) 
    } else {
        sigma_fac <- .9*sigma_guess 
    }
    if ( pendant_guess < max(vargrid$pendant) & pendant_guess > min(vargrid$pendant) ) { 
        pendant_fac <- min( pendant_fac/5, .9*pendant_guess ) 
    } else {
        pendant_fac <- .9*pendant_guess 
    }
}

initvals <- c( coef(simple.lm)[2], coef(simple.lm)[1], sigma_guess, pendant_guess )
parscale <- abs(initvals)
ans <- optim( par=initvals, fn=loglik, lower=c(-Inf,-Inf,initvals[3],initvals[4]), method="L-BFGS-B", 
    control=list(parscale=abs(initvals)/10,fnscale=max(1,abs(loglik(initvals))),trace=1,maxit=1000) )
if (ans$convergence==52) {
    ans <- optim( par=initvals, fn=loglik, lower=c(-Inf,-Inf,initvals[3],initvals[4]), method="L-BFGS-B", 
        control=list(parscale=abs(initvals)/100,fnscale=max(1,abs(loglik(initvals))),trace=1,maxit=1000) )
}
ans #output the output

stopifnot(ans$convergence==0)

# ok, plot
anstree <- newtree
anstree$edge.length[pendant_edges] <- orig_lengths + ans$par[4]
plot(anstree)

# with(newdata, points( tapply(body_length,species,mean), tapply(rib_volume,species,mean), pch=20, col=color_choices, cex=2 ) )
color_choices <- rainbow(nlevels(mydata$species))

pdf(file=paste(outfile, ".", y_variable, ".pdf", sep=""))
with(mydata, plot(mydata[,x_variable], mydata[,y_variable], xlab=x_variable, ylab=y_variable, log=if (logify){'xy'}else{''}, col=color_choices[species], main=main_title ) )
mtext(paste(if(logify){"log10 "}else{''}, y_variable, "=", if(logify){"log10 "}else{''}, round(ans$par[1], digits=8), "*", x_variable, "+", round(ans$par[2], digits=8), sep=""))
legend("topleft", levels(mydata$species), pch=1, col=color_choices, cex=0.5, bty="n")
abline(ans$par[2],ans$par[1])
abline(coef(simple.lm)[1],coef(simple.lm)[2],col='blue',lty=2)
dev.off()

if(logify){
	#calculate residuals by hand... 
	expected_values<-ans$par[1]*log10(mydata[,x_variable])+ans$par[2]
	expected_values<-10**expected_values
	#points(mydata$body_length, expected_values, col='black', pch=19)
	myresiduals<-log10(mydata[,y_variable]-log10(expected_values))
} else {
	expected_values<-ans$par[1]*(mydata[,x_variable])+ans$par[2]
	myresiduals<-(mydata[,y_variable]-(expected_values))
}

#add residuals into data and output as a table... 
mydata[,paste(y_variable, "_residuals", sep="")]<-myresiduals
write.table(mydata, file=paste(outfile, ".", y_variable, ".residuals", sep=""), quote=FALSE, row.names=FALSE, sep="\t", append=FALSE)

