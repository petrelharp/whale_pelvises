###PETER: just change this path to wherever you put 42_rib.FEMALES.combined.residuals and consensusTree_9species.txt
#everything should run fine after that.
path="./"
data_file="42_rib.FEMALES.combined.residuals"
tree_file="consensusTree_9species.txt"

library(ape)
library(nlme)

coef.corIntraspecific <- function (object, unconstrained = TRUE, ...) 
{
	if (!("corIntraspecific" %in% class(object))) 
		stop("object is not of class \"corIntraspecific\"")
	if (unconstrained) {
		if (attr(object, "fixed")) {
			return(numeric(0))
		}
		else {
			return(as.vector(object))
		}
	}
	aux <- as.vector(object)
	names(aux) <- "intraSpecificCorrelation"
	aux
}

corMatrix.corIntraspecific <- function (object, covariate = getCovariate(object), 
	...) 
{
	if (!("corIntraspecific" %in% class(object))) 
		stop("object is not of class \"corIntraspecific\"")
	if (!any(attr(object, "index"))) 
		stop("object has not been initialized.")
	tree <- attr(object, "tree")
	index <- attr(object, "index")
	baseCorMatrix <- cov2cor(vcv(tree))
	expandedCorMatrix <- baseCorMatrix[index,index]
	intraspEntries <- outer(index,index,"==")
	diag(intraspEntries) <- FALSE
	nsamples <- sapply(index,function(x)sum(index==x))
	expandedCorMatrix[intraspEntries] <- object[1]
	if (!all(is.finite(expandedCorMatrix))) stop("something is wrong.")
	return(expandedCorMatrix)
}

corIntraspecific <- function (value, phy, form = ~1, fixed = FALSE) 
{
	if (length(value) > 1) 
		stop("only one parameter is allowed")
	if (abs(value) >  1) 
		stop("the parameter alpha must be no greater than one")
	if (!inherits(phy, "phylo")) 
		stop("object \"phy\" is not of class \"phylo\"")
	attr(value, "formula") <- form
	attr(value, "fixed") <- fixed
	attr(value, "tree") <- phy
	class(value) <- c("corIntraspecific", "corPhyl", "corStruct")
	value
}


Initialize.corIntraspecific <- function (object, data, species, ...) 
{
	if (missing(species)) { species <- rownames(data) }
   form <- formula(object)
	if (!is.null(getGroupsFormula(form))) {
		attr(object, "groups") <- getGroups(object, form, data = data)
		attr(object, "Dim") <- Dim(object, attr(object, "groups"))
	}
	else {
		attr(object, "Dim") <- Dim(object, as.factor(rep(1, nrow(data))))
	}
	attr(object, "covariate") <- getCovariate(object, data = data)
	phy <- attr(object, "tree")
	index <- match(species, phy$tip.label)
	if (any(is.na(index))) {
		stop("species in data frame do not match tree tip names")
	}
	else {
		attr(object, "index") <- index
	}
	object
}

gls_with_Intraspecific <- function (model, data = sys.frame(sys.parent()), correlation = NULL, 
	weights = NULL, subset, method = c("REML", "ML"), na.action = na.fail, 
	control = list(), verbose = FALSE, ...) 
{
	Call <- match.call()
	controlvals <- glsControl()
	if (!missing(control)) {
		if (!is.null(control$nlmStepMax) && control$nlmStepMax < 
			0) {
			warning("negative control$nlmStepMax - using default value")
			control$nlmStepMax <- NULL
		}
		controlvals[names(control)] <- control
	}
	if (!inherits(model, "formula") || length(model) != 3) {
		stop("\nmodel must be a formula of the form \"resp ~ pred\"")
	}
	method <- match.arg(method)
	REML <- method == "REML"
	if (!is.null(correlation)) {
		groups <- getGroupsFormula(correlation)
	}
	else groups <- NULL
	glsSt <- glsStruct(corStruct = correlation, varStruct = varFunc(weights))
	model <- terms(model, data = data)
	mfArgs <- list(formula = asOneFormula(formula(glsSt), model, 
		groups), data = data, na.action = na.action)
	if (!missing(subset)) {
		mfArgs[["subset"]] <- asOneSidedFormula(Call[["subset"]])[[2]]
	}
	mfArgs$drop.unused.levels <- TRUE
	dataMod <- do.call("model.frame", mfArgs)
	origOrder <- row.names(dataMod)
	if (!is.null(groups)) {
		groups <- eval(parse(text = paste("~1", deparse(groups[[2]]), 
			sep = "|")))
		grps <- getGroups(dataMod, groups, level = length(getGroupsFormula(groups, 
			asList = TRUE)))
		ord <- order(grps)
		grps <- grps[ord]
		dataMod <- dataMod[ord, , drop = FALSE]
		revOrder <- match(origOrder, row.names(dataMod))
	}
	else grps <- NULL
	X <- model.frame(model, dataMod)
	contr <- lapply(X, function(el) if (inherits(el, "factor")) 
		contrasts(el))
	contr <- contr[!unlist(lapply(contr, is.null))]
	X <- model.matrix(model, X)
	if (ncol(X) == 0) 
		stop("no coefficients to fit")
	y <- eval(model[[2]], dataMod)
	N <- nrow(X)
	p <- ncol(X)
	parAssign <- attr(X, "assign")
	fTerms <- terms(as.formula(model), data = data)
	namTerms <- attr(fTerms, "term.labels")
	if (attr(fTerms, "intercept") > 0) {
		namTerms <- c("(Intercept)", namTerms)
	}
	namTerms <- factor(parAssign, labels = namTerms)
	parAssign <- split(order(parAssign), namTerms)
	attr(glsSt, "conLin") <- list(Xy = array(c(X, y), c(N, ncol(X) + 
		1), list(row.names(dataMod), c(colnames(X), deparse(model[[2]])))), 
		dims = list(N = N, p = p, REML = as.integer(REML)), logLik = 0)
	glsEstControl <- controlvals[c("singular.ok", "qrTol")]
	glsSt <- Initialize(glsSt, dataMod, glsEstControl, ...)
	parMap <- attr(glsSt, "pmap")
	numIter <- numIter0 <- 0
	repeat {
		oldPars <- c(attr(glsSt, "glsFit")[["beta"]], coef(glsSt))
		if (length(coef(glsSt))) {
			optRes <- if (controlvals$opt == "nlminb") {
				nlminb(c(coef(glsSt)), function(glsPars) -logLik(glsSt, 
				  glsPars), control = list(trace = controlvals$msVerbose, 
				  iter.max = controlvals$msMaxIter), lower=-1, upper=1)
			}
			else {
				optim(c(coef(glsSt)), function(glsPars) -logLik(glsSt, 
				  glsPars), method = controlvals$optimMethod, 
				  control = list(trace = controlvals$msVerbose, 
					maxit = controlvals$msMaxIter, reltol = if (numIter == 
					  0) controlvals$msTol else 100 * .Machine$double.eps), lower=-1, upper=1)
			}
			coef(glsSt) <- optRes$par
		}
		else {
			optRes <- list(convergence = 0)
		}
		attr(glsSt, "glsFit") <- nlme:::glsEstimate(glsSt, control = glsEstControl)
		if (!needUpdate(glsSt)) {
			if (optRes$convergence) 
				stop(optRes$message)
			break
		}
		numIter <- numIter + 1
		glsSt <- update(glsSt, dataMod)
		aConv <- c(attr(glsSt, "glsFit")[["beta"]], coef(glsSt))
		conv <- abs((oldPars - aConv)/ifelse(aConv == 0, 1, aConv))
		aConv <- c(beta = max(conv[1:p]))
		conv <- conv[-(1:p)]
		for (i in names(glsSt)) {
			if (any(parMap[, i])) {
				aConv <- c(aConv, max(conv[parMap[, i]]))
				names(aConv)[length(aConv)] <- i
			}
		}
		if (verbose) {
			cat("\nIteration:", numIter)
			cat("\nObjective:", format(optRes$value), "\n")
			print(glsSt)
			cat("\nConvergence:\n")
			print(aConv)
		}
		if (max(aConv) <= controlvals$tolerance) {
			break
		}
		if (numIter > controlvals$maxIter) {
			stop("maximum number of iterations reached without convergence")
		}
	}
	glsFit <- attr(glsSt, "glsFit")
	namBeta <- names(glsFit$beta)
	attr(parAssign, "varBetaFact") <- varBeta <- glsFit$sigma * 
		glsFit$varBeta * sqrt((N - REML * p)/(N - p))
	varBeta <- crossprod(varBeta)
	dimnames(varBeta) <- list(namBeta, namBeta)
	Fitted <- fitted(glsSt)
	if (!is.null(grps)) {
		grps <- grps[revOrder]
		Fitted <- Fitted[revOrder]
		Resid <- y[revOrder] - Fitted
		attr(Resid, "std") <- glsFit$sigma/(varWeights(glsSt)[revOrder])
	}
	else {
		Resid <- y - Fitted
		attr(Resid, "std") <- glsFit$sigma/(varWeights(glsSt))
	}
	names(Resid) <- names(Fitted) <- origOrder
	if (controlvals$apVar) {
		apVar <- nlme:::glsApVar(glsSt, glsFit$sigma, .relStep = controlvals[[".relStep"]], 
			minAbsPar = controlvals[["minAbsParApVar"]], natural = controlvals[["natural"]])
	}
	else {
		apVar <- "Approximate variance-covariance matrix not available"
	}
	dims <- attr(glsSt, "conLin")[["dims"]]
	dims[["p"]] <- p
	attr(glsSt, "conLin") <- NULL
	attr(glsSt, "glsFit") <- NULL
	estOut <- list(modelStruct = glsSt, dims = dims, contrasts = contr, 
		coefficients = glsFit[["beta"]], varBeta = varBeta, sigma = glsFit$sigma, 
		apVar = apVar, logLik = glsFit$logLik, numIter = if (needUpdate(glsSt)) numIter else numIter0, 
		groups = grps, call = Call, method = method, fitted = Fitted, 
		residuals = Resid, parAssign = parAssign, na.action = attr(dataMod, 
			"na.action"))
	if (inherits(data, "groupedData")) {
		attr(estOut, "units") <- attr(data, "units")
		attr(estOut, "labels") <- attr(data, "labels")
	}
	attr(estOut, "namBetaFull") <- colnames(X)
	class(estOut) <- "gls"
	estOut
}

Initialize.glsStruct <- function (object, data, control = list(singular.ok = FALSE, qrTol = .Machine$single.eps), 
	...) 
{
	if (length(object)) {
		object[] <- lapply(object, Initialize, data, ...)
		theta <- lapply(object, coef)
		len <- unlist(lapply(theta, length))
		num <- seq_along(len)
		if (sum(len) > 0) {
			pmap <- outer(rep(num, len), num, "==")
		}
		else {
			pmap <- array(FALSE, c(1, length(len)))
		}
		dimnames(pmap) <- list(NULL, names(object))
		attr(object, "pmap") <- pmap
		attr(object, "glsFit") <- nlme:::glsEstimate(object, control = control)
		if (needUpdate(object)) {
			object <- update(object, data)
		}
	}
	object
}

#print("Data: /Users/mattdean/Working/papers/pelves/latest_experiments_pelvics/REPLICATES_EXPERIMENT/42_rib.FEMALES.combined.residuals")
#print("Tree: /Users/mattdean/Working/papers/pelves/latest_experiments_pelvics/REPLICATES_EXPERIMENT/CONSENSUS_TREES_OLD/consensusTree_9species.txt")

mydata<-read.table(file=paste(path, data_file, sep=""), header=TRUE)
mytree<-read.nexus(file=paste(path, tree_file, sep=""))

corStr_Intraspecific<-corIntraspecific(1, phy=mytree, fixed=TRUE)
num_species <- as.vector( table(mydata$species)[match(mydata$species,levels(mydata$species))] )
weights_Intraspecific <- varFixed( ~ num_species )
# corStr_Intraspecific<-corIntraspecific(0.5, phy=mytree, fixed=TRUE)
#note that you are not "logging" the variables here since they are already the residuals on a log scale.  
gls_fit<-gls_with_Intraspecific((residual_rib_volume)~(residual_testis), mydata, corStr_Intraspecific, weights=weights_Intraspecific, species=mydata$species, method="ML")
summary(gls_fit)
gls_fit

####################################
#everything below here is plotting
color_choices <- rainbow(nlevels(mydata$species))
with(mydata, plot(residual_testis, residual_rib_volume, col=color_choices[species], main="RIB, FEMALES, COMBINED", ylim=range(residual_rib_volume,fitted(gls_fit)) ) ) #took out the "log='xy'" for 42...py
with(mydata, points( residual_testis, fitted(gls_fit), col=color_choices[species], pch=20 ) )
abline(gls_fit$coefficients["(Intercept)"], gls_fit$coefficients["residual_testis"], untf=TRUE)
mtext(paste("(V)=", round(gls_fit$coefficients["residual_testis"], digits=8), "*(L)+", round(gls_fit$coefficients["(Intercept)"], digits=8), sep=""))
legend("topleft", levels(mydata$species), pch=1, col=color_choices, cex=0.5, bty="n")

###
# Huh, maybe the problem is non-pos-def'ness?
# ok, fine.  Make a new tree with samples in as edges.
#  THIS DOES NOT fix the problem.
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

stopifnot( setequal( rownames(newdata), newtree$tip.label ) )

new_corstr <- corBrownian(.01, phy=newtree)
num_species <- as.vector( table(newdata$species)[match(newdata$species,levels(newdata$species))] )
weights_Intraspecific <- varFixed( ~ num_species )
gls_fit<-gls( residual_rib_volume~residual_testis, newdata, new_corstr, weights=weights_Intraspecific, method="ML")
summary(gls_fit)
gls_fit

color_choices <- rainbow(nlevels(newdata$species))
with(newdata, plot(residual_testis, residual_rib_volume, col=color_choices[species], main="RIB, FEMALES, COMBINED", ylim=range(residual_rib_volume,fitted(gls_fit)) ) ) #took out the "log='xy'" for 42...py
with(newdata, points( residual_testis, fitted(gls_fit), col=color_choices[species], pch=20 ) )
abline(gls_fit$coefficients["(Intercept)"], gls_fit$coefficients["residual_testis"], untf=TRUE)
mtext(paste("(V)=", round(gls_fit$coefficients["residual_testis"], digits=8), "*(L)+", round(gls_fit$coefficients["(Intercept)"], digits=8), sep=""))
legend("topleft", levels(newdata$species), pch=1, col=color_choices, cex=0.5, bty="n")


#####
# ok, DIY or die
require(MASS)
cor_matrix <- corMatrix(Initialize(new_corstr,newdata))

stopifnot(all(eigen(cor_matrix)$values>=-1e-16))
cov_matrix <- cor_matrix * var(mydata$residual_rib_volume)
fixed_gls_fit <- lm.gls( (residual_rib_volume)~(residual_testis), data=mydata, W=cov_matrix, inverse=TRUE )

with(newdata, points( tapply(residual_testis,species,mean), tapply(residual_rib_volume,species,mean), pch=20, col=color_choices, cex=2 ) )
with(newdata, points( residual_testis, fitted(fixed_gls_fit), col=color_choices[species], pch=3 ) )
abline(fixed_gls_fit$coefficients["(Intercept)"], fixed_gls_fit$coefficients["residual_testis"], untf=TRUE)

group_vars <- with( mydata, tapply( residual_rib_volume, species, var ) )
withingroup_var <- mean( group_vars, na.rm=TRUE )
withingroup_cor <- 1-withingroup_var/var(mydata$residual_rib_volume,)


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

abline(ans$par[2],ans$par[1],col='purple')
