# Analyzing Raw Genomic Sequences Using RStudio
Analyzing the raw fastq files prior to trimming and after trimming. The same code format will be used with minor adjustments to where the files are pulled from.

## FASTQC 
The installation only needs to be completed once.

```
install.packages("fastqcr")
```

I ran another install just incase and following the information at the following link:

https://cran.r-project.org/web/packages/fastqcr/readme/README.html

```
fastqc_install()
```

Below is the code that was run the FastQC on each individual .fastq file. All of the .fastq files were placed in a singular folder and that folder was indicated as the working directory

```
library(fastqcr)

fastqc(fq.dir = "#where your files are stored",
       qc.dir = "#where the reports are going to go"
       threads = 4 #how much CPU (cores) is dedicated to this function, default is 4
)
```

The FastQC reports that were generated are put into a folder indicated in the qc.dir portion of the function. The files generated are in both .html and .zip files.

To interepret the results of the FastQC reports reference the following link:

All results:
https://rtsf.natsci.msu.edu/genomics/technical-documents/fastqc-tutorial-and-faq.aspx

Per Tile Sequence Quality:
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/12%20Per%20Tile%20Sequence%20Quality.html

Sequence Length Distribtution:
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/7%20Sequence%20Length%20Distribution.html

## MultiQC
This is done to compile all of the FastQC report that were generated during the FastQC. This is known as aggregating the files. For additional information on this function: 

https://rpkgs.datanovia.com/fastqcr/qc-reports/fastqcr-multi-qc-report.html

```
library(fastqcr)
library(dplyr)

qc <- qc_aggregate(qc.dir, progressbar = FALSE) #qc.dir will be the same was what was used above to store the reports
qc #prints the summary table inline

summary(qc)

qc_stats(qc)
```
### Interpretation of the summary() Results
* _Modules_: fastqc modules
* _nb_samples_: number of samples tested
* _nb_pass,nb_fail, nb_warn_: the number of samples that passed. failed, and warned
* _failed, warned_: name of samples that failed and warned

### General Statistics with qc_stats()
* _pct.dup_: percentage of duplicate reads
* _pct.gc_: percentage of GC content
* _tot.seq_: total sequences or the number of reads
* _seq.length_: sequence length or the length of reads
