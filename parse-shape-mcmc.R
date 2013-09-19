#/usr/bin/R --vanilla
# parse the mcmc output to infer parameters for shape evolution

source("correlated-traits-fns.R")
require(ape)
require(Matrix)
require(colorspace)

mcmcfiles <- list( pelvic='pelvic-shape-brlens/5678-mcmc-run.RData',
            rib="rib-shape-brlens/0806-mcmc-run.RData",
            'sub-pelvic'="sub.pelvic-shape-brlens/7401-mcmc-run.RData" )
mcruns <- lapply( mcmcfiles, function (x) { tmpenv <- environment(); load(x,envir=tmpenv); as.list(tmpenv) } )
sapply(lapply(mcruns, "[[", 'mcrun'), "[[", "accept")

layout(matrix(1:6,nrow=3))
lapply( mcruns, function (x) matplot(x$mcrun$batch[floor(seq(1,nrow(x$mcrun$batch),length.out=500)),],type='l' ) )
lapply( mcruns, function (x) matplot(x$mcrun$batch[floor(seq(1,nrow(x$mcrun$batch),length.out=500)),],type='l', log='y') )

load("shape-brlens-stuff.RData")
load("spmapping.RData")
load("shared-num-bones.RData")
load("sampled-edge-testes.RData")

# which branches are irrelevant for rib analysis?
missing.taxa <- lapply( list(pelvic=pelvic.speciesdiffsq, rib=rib.speciesdiffsq), function (x) names( which( apply(is.na(x),1,all) ) ) )
stopifnot( all( missing.taxa[[1]] %in% missing.taxa[[2]] ) )
spdesc <- get.descendants(sptree)
sp.edge.indices <- sptree$edge[,2] # associate each edge with the downstream node
spedgedesc <- spdesc[sp.edge.indices,]
missing.tipindices <- lapply( missing.taxa, function (x) match(x,sptree$tip.label) )
missing.edges <- lapply( missing.tipindices, function (x) { rowSums( spedgedesc[,setdiff(1:Ntip(sptree),x)] ) == 0 } )
missing.edges['sub-pelvic'] <- missing.edges['rib']

testes.tree <- sptree
testes.tree$edge.length <- abs(sp.edge.testes) * sptree$edge.length

branch.pars <- 4:ncol(mcruns[[1]]$mcrun$batch)
burnin <- min(sapply(mcruns,function(x)nrow(x$mcrun$batch)))/2
posterior.means <- lapply( seq_along(mcruns), function (k) {
            x <- mcruns[[k]]$mcrun$batch
            z <- colMeans(x[burnin:nrow(x),]) 
            # z[branch.pars][ missing.edges[[ names(mcruns)[k] ]] ] <- 0
            return(z)
        } )
names(posterior.means) <- names(mcruns)
mean.trees <- lapply( posterior.means, function (x) { sptree$edge.length <- x[branch.pars]; sptree } )

edge.widths <- lapply( seq_along(mean.trees), function (k) {
            medges <- if (names(mean.trees)[k]=='pelvic') { missing.edges[['pelvic']] } else { missing.edges[['rib']] }
            ifelse( medges, 0, 3 )
        } )
edge.cols <- diverge_hcl(32,c=100,l=c(20,50))[cut(sp.edge.testes,32)]


layout(matrix(1:6,nrow=3))
plot(sptree, edge.color=edge.cols, edge.width=3 ) 
mtext('phylogeny',3)
plot( testes.tree, edge.color=edge.cols, edge.width=3 )
mtext('testes-weighted', 3)
legend("bottomright",fill=edge.cols[c(which.max(sp.edge.testes),which.min(sp.edge.testes))], legend=c(max(sp.edge.testes),min(sp.edge.testes)), title='value' )
for (k in seq_along(mean.trees)) {
    plot(mean.trees[[k]], edge.color=edge.cols, edge.width=edge.widths[[k]] )
    mtext(names(mean.trees)[k],3)
}

layout(matrix(1:4,nrow=2))
whichedges <- lapply( seq_along(mean.trees), function (k) {
        pmeans <- posterior.means[[k]][branch.pars]
        usethese <- !missing.edges[[ names(mean.trees)[k] ]]
        plot( sp.edge.testes[usethese], (pmeans/sptree$edge.length)[usethese], main=names(posterior.means)[k], xlab='relative testes size', ylab='speed of shape evolution' )
        which(usethese)[identify( sp.edge.testes[usethese], (pmeans/sptree$edge.length)[usethese], labels=which(usethese) )]
    } )
usethese <- !missing.edges[[ 'rib' ]]
plot( posterior.means[['rib']][branch.pars][usethese], posterior.means[['sub-pelvic']][branch.pars][usethese], xlab='speed of rib shape evolution', ylab='speed of pelvic shape evolution' )
whichedges <- c( whichedges, list( which(usethese)[identify( posterior.means[['rib']][branch.pars][usethese], posterior.means[['sub-pelvic']][branch.pars][usethese], labels=which(usethese) )] ) )
plot(sptree)
edgelabels(edge=do.call(c,whichedges))

###
# back to the ol' lm
edgedata <- as.data.frame( 
