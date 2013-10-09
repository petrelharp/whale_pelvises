# pairwise testes difference
load("sampled-edge-testes.RData")
layout(1:2)
lapply( brlen.fits, function (x) {
            rlens <- x[-(1:2)] / sptree$edge.length
            plot( rlens, sp.edge.testes, xlab='speed of shape change', xlim=c(0,5) )
            tmp.lm <- lm( sp.edge.testes ~ rlens, subset=(rlens<=5) )
            abline(coef(tmp.lm))
            summary(tmp.lm)
        } )


with( as.data.frame(fullmean), pairs( data.frame( bodylength, pelvic=(left.pelvic+right.pelvic)/2, rib=(left.rib+right.rib)/2, testes=actual_testes_mass_max, scaled.testes ) ) )

with( as.data.frame(fullmean), plot( actual_testes_mass_max ~ bodylength ) )

time.dist <- treedist(adjtree)
length.dist <- with( as.data.frame(thedata), outer( bodylength, bodylength, "-" ) )
testes.dist <- with( as.data.frame(thedata), outer( actual_testes_mass_max, actual_testes_mass_max, "-" ) )

write.csv(time.dist,file="time-dist.csv")
write.csv(length.dist,file="length-dist.csv")
write.csv(testes.dist,file="testes-dist.csv")

tmp1 <- with( subset(shapediff,bone1==bone2 & bone1=='pelvic'), tapply( shape_difference, list(specimen1,specimen2), length ) )
tmp2 <- with( subset(shapediff,bone1==bone2 & bone1=='pelvic'), tapply( shape_difference, list(specimen2,specimen1), length ) )
denoms <- ifelse( is.na(tmp1), 0, tmp1 ) + ifelse( is.na(tmp2), 0, tmp2 )

tmp1 <- with( subset(shapediff,bone1==bone2 & bone1=='pelvic'), tapply( shape_difference, list(specimen1,specimen2), sum, na.rm=TRUE ) )
tmp2 <- with( subset(shapediff,bone1==bone2 & bone1=='pelvic'), tapply( shape_difference, list(specimen2,specimen1), sum, na.rm=TRUE ) )
pelvic.dist  <- ( ifelse( is.na(tmp1), 0, tmp1 ) + ifelse( is.na(tmp2), 0, tmp2 ) ) / denoms

tmp1 <- with( subset(shapediff,bone1==bone2 & bone1=='rib'), tapply( shape_difference, list(specimen1,specimen2), length ) )
tmp2 <- with( subset(shapediff,bone1==bone2 & bone1=='rib'), tapply( shape_difference, list(specimen2,specimen1), length ) )
denoms <- ifelse( is.na(tmp1), 0, tmp1 ) + ifelse( is.na(tmp2), 0, tmp2 )

tmp1 <- with( subset(shapediff,bone1==bone2 & bone1=='rib'), tapply( shape_difference, list(specimen1,specimen2), sum, na.rm=TRUE ) )
tmp2 <- with( subset(shapediff,bone1==bone2 & bone1=='rib'), tapply( shape_difference, list(specimen2,specimen1), sum, na.rm=TRUE ) )
rib.dist  <- ( ifelse( is.na(tmp1), 0, tmp1 ) + ifelse( is.na(tmp2), 0, tmp2 ) ) / denoms

thedata.species <- whales$species[ match(rownames(thedata),whales$specimen) ]
thedata.species[ rownames(thedata) %in% levels(whales$species) ] <- rownames(thedata)[ rownames(thedata) %in% levels(whales$species) ]
sprow <- thedata.species[row(time.dist)]
spcol <- thedata.species[col(time.dist)]
thedata.species <- whales$sex[ match(rownames(thedata),whales$specimen) ]
thedata.species[ rownames(thedata) %in% levels(whales$species) ] <- "M"
sexrow <- thedata.species[row(time.dist)]
sexcol <- thedata.species[col(time.dist)]

sp.time.dist <- tapply( time.dist, list( sprow,spcol,sexrow,sexcol), mean, na.rm=TRUE )
sp.length.dist <- tapply( length.dist, list( sprow,spcol,sexrow,sexcol), mean, na.rm=TRUE )
sp.testes.dist <- tapply( testes.dist, list( sprow,spcol,sexrow,sexcol), mean, na.rm=TRUE )
sp.pelvic.dist <- tapply( pelvic.dist, list( sprow,spcol,sexrow,sexcol), mean, na.rm=TRUE )
sp.rib.dist <- tapply( rib.dist, list( sprow,spcol,sexrow,sexcol), mean, na.rm=TRUE )

tmp.alldists <- lapply( list( 'F-F'=c(1,1), 'F-M'=c(1,2), 'M-F'=c(2,1), 'M-M'=c(2,2) ), function (kk) {
        do.call( cbind, lapply( list( bodylength=sp.length.dist[,,kk[1],kk[2]], testes=sp.testes.dist[,,kk[1],kk[2]], time=sp.time.dist[,,kk[1],kk[2]], pelvic=sp.pelvic.dist[,,kk[1],kk[2]], rib=sp.rib.dist[,,kk[1],kk[2]] ), as.vector ) )
    } )
alldists <- do.call( rbind, lapply( 1:4, function (k) { 
                    x <- as.data.frame(tmp.alldists[[k]]); 
                    x$sex <- factor(names(tmp.alldists)[k],levels=names(tmp.alldists)); x 
                } ) )

pairs( alldists, 
        lower.panel=function(x,y,...) points( x[alldists$sex=="M-M"], y[alldists$sex=="M-M"], ... ),
        upper.panel=function(x,y,...) { if (any(!is.na(x[alldists$sex=="F-F"]) & !is.na(y[alldists$sex=="F-F"]))) { points( x[alldists$sex=="F-F"], y[alldists$sex=="F-F"], ... ) } else { NULL } },
        col=ifelse(pair.cols,'red','black'), pch=ifelse(pair.cols,20,1), cex=ifelse(pair.cols,2,1)
        )
dev.set( dev.next() )
pairs( alldists, 
        lower.panel=function(x,y,...) { if (any(!is.na(x[alldists$sex=="F-M"]) & !is.na(y[alldists$sex=="M-F"]))) { points( x[alldists$sex=="F-M"], y[alldists$sex=="M-F"], ... ) } else { NULL } },
        upper.panel=function(x,y,...) { if (any(!is.na(x[alldists$sex=="M-F"]) & !is.na(y[alldists$sex=="F-M"]))) { points( x[alldists$sex=="M-F"], y[alldists$sex=="F-M"], ... ) } else { NULL } },
        col=ifelse(pair.cols,'red','black'), pch=ifelse(pair.cols,20,1), cex=ifelse(pair.cols,2,1)
        )


sp.pairs <- list(
        c("INIA_GEOFFRENSIS", "PONTOPORIA_BLAINVILLEI"), 
        c("PHOCOENOIDES_DALLI", "PHOCOENA_PHOCOENA"), 
        c("DELPHINUS_DELPHIS", "DELPHINUS_CAPENSIS"), 
        c("TURSIOPS_ADUNCUS", "STENELLA_FRONTALIS"), 
        c("TURSIOPS_TRUNCATUS", "STENELLA_COERULEOALBA"), 
        c("STENELLA_ATTENUATA", "STENELLA_LONGIROSTRIS"), 
        c("STENO_BREDANENSIS", "LAGENORHYNCHUS_ACUTUS"), 
        c("FERESA_ATTENUATA", "GRAMPUS_GRISEUS"), 
        c("LISSODELPHIS_BOREALIS", "LAGENORHYNCHUS_OBLIQUIDENS"), 
        c("ESCHRICHTIUS_ROBUSTUS", "BALAENOPTERA_MUSCULUS"), 
        c("BALAENOPTERA_ACUTOROSTRATA", "EUBALAENA_GLACIALIS")
        )
tmp1 <- match( sapply( sp.pairs, "[", 1 ), rownames(sp.time.dist) )
tmp2 <- match( sapply( sp.pairs, "[", 2 ), rownames(sp.time.dist) )

ut.pair.cols <- matrix( FALSE, nrow=nrow(sp.time.dist), ncol=ncol(sp.time.dist) )
ut.pair.cols[ rbind(cbind(tmp1,tmp2)) ] <- TRUE
pair.cols <- ( ut.pair.cols | t(ut.pair.cols) )

layout(1:2)
with( subset(alldists, sex=="M-M"), plot( pelvic ~ abs(testes), col=ifelse(pair.cols,'red',adjustcolor('black',.15)), pch=ifelse(pair.cols,20,1) ) )
with( subset(alldists, sex=="M-M"), identify( abs(testes), pelvic, labels=outer(dimnames(sp.time.dist)[[1]],dimnames(sp.time.dist)[[2]],paste) ) )
abline(coef(with( subset(alldists, sex=="M-M"), lm( pelvic ~ abs(testes)) )))
with( subset(alldists, sex=="M-M"), plot( rib ~ abs(testes), col=ifelse(pair.cols,'red',adjustcolor('black',.25)), pch=ifelse(pair.cols,'o','.') ) )
abline(coef(with( subset(alldists, sex=="M-M"), lm( rib ~ abs(testes)) )))

summary(with( subset(alldists, sex=="M-M"), lm( rib ~ abs(testes)), subset=ut.pair.cols ))
summary(with( subset(alldists, sex=="M-M"), lm( pelvic ~ abs(testes)), subset=ut.pair.cols ))

summary(with( subset(alldists, sex=="M-M"), lm( rib ~ abs(testes) + time, subset=ut.pair.cols ) ))
summary(with( subset(alldists, sex=="M-M"), lm( pelvic ~ abs(testes) + time, subset=ut.pair.cols ) ))

layout( matrix(1:4,nrow=2,byrow=TRUE) )
for (bone in c("rib","pelvic")) {
    y <- subset(alldists, sex=="M-M")[,bone]
    with( subset(alldists, sex=="M-M"), plot( y ~ abs(testes), col=ifelse(pair.cols,'red',adjustcolor('black',.15)), pch=ifelse(pair.cols,20,1), ylab=bone ) )
    with( subset(alldists, sex=="M-M"), plot( y ~ time, col=ifelse(pair.cols,'red',adjustcolor('black',.15)), pch=ifelse(pair.cols,20,1), ylab=bone ) )
}


