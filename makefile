ALL : females/thedata-and-covmatrices.Rdata males/thedata-and-covmatrices.Rdata morphology_table_2013_June_27-plr.txt

morphology_table_2013_June_27-plr.txt : morphology_table_2013_June_27.txt
	cat morphology_table_2013_June_27.txt | cut -f 1-7 -d '    ' > morphology_table_2013_June_27-plr.txt

# males
males/thedata-and-covmatrices.Rdata : parse-correlated-trait-data.R
	cd males
	Rscript parse-correlated-trait-data.R

males/tree-stuff.RData : parse-correlated-trait-data.R
	cd males
	Rscript $<

males/all-sample-tree.RData: parse-correlated-trait-data.R
	cd males
	Rscript $<

males/mcmc-setup.RData : correlated-traits-analysis.R males/thedata-and-covmatrices.Rdata males/all-sample-tree.RData
	cd males
	Rscript $<

males/analysis-results.RData : correlated-traits-analysis.R males/thedata-and-covmatrices.Rdata males/all-sample-tree.RData
	cd males
	Rscript $<

males/results.RData : parse-mcmc.R males/mcmc-setup.RData males/thedata-and-covmatrices.Rdata
	cd males
	Rscript $<

# females
females/thedata-and-covmatrices.Rdata : parse-correlated-trait-data.R
	cd females
	Rscript parse-correlated-trait-data.R

females/tree-stuff.RData : parse-correlated-trait-data.R
	cd females
	Rscript $<

females/all-sample-tree.RData: parse-correlated-trait-data.R
	cd females
	Rscript $<

females/mcmc-setup.RData : correlated-traits-analysis.R females/thedata-and-covmatrices.Rdata females/all-sample-tree.RData
	cd females
	Rscript $<

females/analysis-results.RData : correlated-traits-analysis.R females/thedata-and-covmatrices.Rdata females/all-sample-tree.RData
	cd females
	Rscript $<

females/results.RData : parse-mcmc.R females/mcmc-setup.RData females/thedata-and-covmatrices.Rdata
	cd females
	Rscript $<

# # shape stuff
# shape-stuff.RData : shape-analysis-setup.R tree-stuff.RData thedata-and-covmatrices.Rdata mcmc-setup.RData results.RData
# 	Rscript $<

