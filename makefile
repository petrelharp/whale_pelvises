ALL : correlated-traits-analysis.R thedata-and-covmatrices.Rdata morphology_table_2013_June_27-plr.txt

morphology_table_2013_June_27-plr.txt : morphology_table_2013_June_27.txt
	cat morphology_table_2013_June_27.txt | cut -f 1-7 -d '    ' > morphology_table_2013_June_27-plr.txt

thedata-and-covmatrices.Rdata : parse-correlated-trait-data.R
	Rscript parse-correlated-trait-data.R


