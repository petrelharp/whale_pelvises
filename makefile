ALL : correlated-traits-analysis.R thedata-and-covmatrices.Rdata morphology_table_2013_June_27-plr.txt

morphology_table_2013_June_27-plr.txt : morphology_table_2013_June_27.txt
	cat morphology_table_2013_June_27.txt | cut -f 1-7 -d '    ' > morphology_table_2013_June_27-plr.txt

thedata-and-covmatrices.Rdata : parse-correlated-trait-data.R
	Rscript parse-correlated-trait-data.R

tree-stuff.RData : parse-correlated-trait-data.R
	Rscript $<

all-sample-tree.RData: parse-correlated-trait-data.R
	Rscript $<

mcmc-setup.RData : correlated-traits-analysis.R thedata-and-covmatrices.Rdata analysis-results.RData all-sample-tree.RData
	Rscript $<

analysis-results.RData : correlated-traits-analysis.R thedata-and-covmatrices.Rdata all-sample-tree.RData
	Rscript $<

results.RData : parse-mcmc.R mcmc-setup.RData thedata-and-covmatrices.Rdata
	Rscript $<

shape-stuff.RData : shape-analysis-setup.R tree-stuff.RData thedata-and-covmatrices.Rdata mcmc-setup.RData results.RData
	Rscript $<

