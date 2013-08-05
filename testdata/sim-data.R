#!/usr/bin/R

require(ape)
if (!exists("scriptdir")) {scriptdir <- "."}
source(paste(scriptdir,"correlated-traits-fns.R",sep="/"))

tree_file <- "consensusTree_ALL_CETACEA.tree"
species_tree<-read.nexus(file=tree_file)

