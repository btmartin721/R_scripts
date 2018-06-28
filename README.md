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

Test                D-Range             Z-Range            No.Significant  Bonferroni      Population-Z    Pop_P-value
ON+EAAL+GULA+MX     [0.0176,0.3735]     [-4.5141,-0.1234]  50/90           15/90          -3.74605         0.000179640801360612
ON+EAAL+GULA+TTAR   [-0.0817,0.3577]    [-3.1064,0.6354]   15/90           0/90           -1.4051          0.159991544126064
...
```

**NOTE**: You must specify a working directory in R prior to running the script.

