ALL : correlated-traits-analysis.R thedata-and-covmatrices.Rdata morphology_table_2013_June_27-plr.txt

morphology_table_2013_June_27-plr.txt : morphology_table_2013_June_27.txt
	cat morphology_table_2013_June_27.txt | cut -f 1-7 -d '    ' > morphology_table_2013_June_27-plr.txt

thedata-and-covmatrices.Rdata : parse-correlated-trait-data.R
	Rscript parse-correlated-trait-data.R

tree-stuff.RData : parse-correlated-trait-data.R
	Rscript $<

all-sample-tree.RData: parse-correlated-trait-data.R
	Rscript $<

correlated-traits-analysis.R : thedata-and-covmatrices.Rdata analysis-results.RData all-sample-tree.RData

mcmc-setup.RData : correlated-traits-analysis.R
	Rscript $<

analysis-results.RData : correlated-traits-analysis.R
	Rscript $<

parse-mcmc.R : mcmc-setup.RData thedata-and-covmatrices.Rdata

results.RData : parse-mcmc.R
	Rscript $<

shape-analysis.R : tree-stuff.RData thedata-and-covmatrices.Rdata mcmc-setup.RData results.RData
