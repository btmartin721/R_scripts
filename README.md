# R_scripts
Scripts written in R

## parseDtests.R 

Takes a directory of D-test files from Steve Mussmann's [Comp-D_MPI](https://github.com/smussmann82/Comp-D_MPI) or [Comp-D](https://github.com/smussmann82/Comp-D) software.  

It then writes a summary table for all the files to "summary.tsv"    
The script assumes that there are two files for each test with extensions .txt and .z  

The .txt file should have 13 columns and contain one row for each permutation of individuals in the O, P3, P2, and P1 populations. Column 12 should have the Z-scores and column 8 the D-statistics.  

The .z file should have one row with three columns containing the test type (must be 4-taxon D-test) and population average Z-score and P-value.  

```
Example *.txt file:  

O              P3           P2            P1            ABBA     BABA  nloci   D      STDEV    X^2    X^2_pval  Z-score Z pval
ONIL_BXON25   EANC_BX318    TTLA_BX422    MXMX_BX1196   25.8125 27.0625 2146  -0.0236 0.0883  0.0296  0.8635    0.2678  1.0000
ONIL_BXON25   EANC_BX318    TTLA_BX422    MXMX_BX1195   25.3125 24.0625 2125  0.0253  0.0960  0.0316  0.8588   -0.2638  1.0000
...
```

```
Example *.z file:  

Statistic   Z-score   P-val
D           -1.2633   0.317311
```

```

Example output (summary.tsv):

Test	              Mean.D (STDEV.D)	  Z-Range	          Chi-Square.Sig	Z.Significant	Bonferroni.Alpha	Bonferroni.Significant	Population-Z	Pop_P-value
ON+EAAL+GULA+MX	    0.22677 (0.06088)	  [-4.5141,-0.1234]	36/90	                  50/90	         0.00056	                 15/90	             -3.74605	     0.00018
ON+EAAL+GULA+TTAR	  0.13508 (0.09668)	 [-3.1064,0.6354  14/90	                 15/90	        0.00056	                          0/90              -1.4051       0.15999
...
```

**NOTE**: You must specify a working directory in R prior to running the script.


## get_dnaex_data.R

Takes a directory of .xls files in a specific format (Douglas Lab DNAex forms) and outputs data from several columns.  
must setwd() in R prior to using.  
Must be in .xls, and NOT .xlsx  
Also, Cannot have hidden .xls files in directory.  
If the format of your DNAex files differs from mine, you might have to modify the code.  

## get_sequenced.R

Takes a directory of Douglas lab ligation files as input, and outputs sequenced samples to .txt file.  
Must setwd() in R prior to using.  All files must be .xlsx, NOT .xls  
If the formatting for your ligation sheets differs from mine, you might have to modify the code.  

## hybridDetective / parallelnewhybrid scripts

Takes genepop input file and runs through the hybriddetective pipeline (https://github.com/bwringe/hybriddetective) conjunction with parallelnewhybrid (https://github.com/bwringe/parallelnewhybrid).  

There are four R scripts associated with this pipeline: Steps 1-2, step 3, steps 4-7, and step 8. I did it that way because I was running parallelnewhybrid (steps 3 and 8) on a different computer, and then I'd collect the output and run the next steps.  

Each R script requires that you change the settings and set the working directory in the SETTINGS sections at the top of the scripts. It is also preferable that you run each script one function at a time, so you can make sure everything works OK.  

Take note of where the forward slashes are when specifying filenames and directories.  


## genepop2hzar.R  

This script takes a genepop file as input and generates the four files needed to run my HZAR_runSingle.R script (see below).  

You need to modify the settings and working directory at the top of the script prior to running it.  
If you want to change the names of the four output files, you can specify them towards the bottom of the script. This was a quick and dirty script, so they are hard-coded in.  

The script will subset diagnostic (clinal) loci based on allele frequencies >= 0.8 at one end of the cline and <=0.2 at the other end.  


## HZAR_runSingle.R  

R script to run HZAR - Hybrid Zone Analysis using R.  

The script only runs one locus at a time, so you have to make sure all the settings are correct and run it from the command-line.  

When running from the command-line, it takes one argument: The index of the locus from the loci.txt file

Example usage:

`Rscript --vanilla HZAR_runSingle.R 2`

The above runs the locus at index 2 (from the locinames.csv file) through HZAR. If you want to run multiple loci sequentially, say all 12 loci, do as follows:  

```for i in `seq 1 12`; do Rscript --vanilla HZAR_runSingle.R $i; done >> log.txt 2>&1```   

If you want to run multiple loci (e.g., all 12 of them) in parallel, you can use GNU parallel as follows:  

```parallel "Rscript --vanilla HZAR_runSingle.R {} >loc{}.logfile.txt 2>&1" :::: <(seq 1 12);```

This requires GNU parallel to be installed. For an example bash script that I used to run HZAR in parallel (on a TORQUE PBS Job Scheduler), see my [pbs_scripts repository](https://github.com/btmartin721/pbs_scripts)  

The script needs four csv files to run:   
1. a locinames file with the locus names for all loci placed on one comma-separated row.  
2. a one-column distances file with the distances for each population, starting at 0 and ending at the furthest distance.  
3. a refAlleles file with the allele frequencies from the populations.  
4. an nsamples file with the sample sizes for each allele from each population.  

I generated the input files using the genepop2hzar.R script in this repository.  
 
## mle2bfd.Rmd  

This is an Rmarkdown script to parse the output from Bayes Factor Delimitation - BFD - runs. It outputs two tables: One supplementary table containing a lot of info about the analysis, and another more final table to summarize the BFD results and rank the Bayes Factors.  

The script requires the presence of a comma-separated file containing specific column names in the following order:  

```
BFD_run
model
K
posterior_ESS
likelihood_ESS
MLE_ESS_min
MLE_ESS_max
MLE_ESS_avg
MLE_ESS_median
MLE	MLE_SD
chainLength
burnin
MLE_burnin%
MLE_CV_intervals
MLE_CV_reps	
Extended?
BF	
rank	
support
```

I have included a template CSV file that you can use: BFDdata.csv   

All the column names have to be present, but the values in some of the rows can be NA. For example, if you have runs that aren't finished running yet. If the script detects any NA values in a row, it doesn't include that run in the tables.    

The values for all the columns were obtained either by copying the results from the PathSamplingAnalysis in BEAST2's app launcher (MODEL_SELECTION package) to excel. I then used Excel functions to calculate e.g., min, max, mean, and median ESS values from the PathSampling results.  The likelihood_ESS and posterior_ESS were obtained by assessing the likelihood traces in Tracer.  

BFD_run is just the number of runs, sequentially 1 to Nruns. The model values can be any string identifying which model you are testing (e.g., East vs. West). K is the number of tips (species) in the BFD analysis.  

The Extended?, rank, BF, and support column headers must be present, but the values are left empty. The Rmarkdown script will fill them in.  
