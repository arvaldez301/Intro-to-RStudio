# What is Bioconductor?
## Introduction to the Bioconductor Project
### Finding the Bioconductor Version
```
# Load the BiocManager package
library(BiocManager)

# Explicitly check the Bioconductor version
version()
```
### BiocManager to Install
```BSgenome``` is a Bioconductor data package that contains representations of several genomes. This package has already been installed for you, as installing the dependencies usually takes some time, using the following code:
```
# Load the BSgenome package
library(BSgenome)

# Check the versions of the packages loaded in the session
sessionInfo()
```
## The role of S4 in Bioconductor
### S4 class definition
We will use the class BSgenome, which is already loaded for you.

Let's check the formal definition of this class by using the function ```showClass("className")```. Check the ```BSgenome``` class results and find its parent classes (Extends) and the classes that inherit from it (Subclasses).
Answer: Annotated and MaskedBSgenome
## Interaction with classes
Let's say we have an object called ```a_genome``` from class ```BSgenome```. With ```a_genome```, you can ask questions like these:
```
# What is a_genome's main class?
class(a_genome)  # "BSgenome"

# What is a_genome's other classes?
is(a_genome)  # "BSgenome", "Annotated"

# Is a_genome an S4 representation?
isS4(a_genome)  # TRUE
If you want to find out more about the a_genome object, you can use the accessor show(a_genome) or use other specific accessors from the list of .S4methods(class = "BSgenome").
```
```
# Investigate the a_genome using show()
show(a_genome)

# Investigate some other accesors
organism(a_genome)
provider(a_genome)
seqinfo(a_genome)
```
## Introducing biology of genomic datasets
### Discovering the yeast genome
You'll begin to explore the yeast genome for yourself using the package ```BSgenome.Scerevisiae.UCSC.sacCer3```, which is already installed for you.

As with other data in R, you can use ```head()``` and ```tail()``` to explore the yeastGenome object. You can also subset the genome by chromosome by using $ syntax as follows: ```object_name$chromosome_name```.

The names of the chromosomes can be obtained using the ```names()``` function, and ```nchar()``` can be used to count the number of characters in a sequence.
```
# Load the yeast genome
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Assign data to the yeastGenome object
yeastGenome <- BSgenome.Scerevisiae.UCSC.sacCer3

# Get the head of seqnames and tail of seqlengths for yeastGenome
head(seqnames(yeastGenome))
tail(seqlengths(yeastGenome))

# Print chromosome M, alias chrM
print(yeastGenome$chrM)

# Count characters of the chrM sequence
nchar(yeastGenome$chrM)
```
### Partitioning the yeast genome
Genomes are often big, but interest usually lies in specific regions of them. Therefore, we need to subset a genome by extracting parts of it. To pick a sequence interval, use ```getSeq()``` and specify the name of the chromosome and the start and end of the sequence interval.

The following example will select the bases of "chrI" from 100 to 150.

```getSeq(yeastGenome, names = "chrI", start = 100, end = 150)```
Note: names is optional; if not specified, it will return all chromosomes. The parameters start and ```end``` are also optional and, if not specified, will take the default values ```1``` and the length of the sequence, respectively.
```
# Load the yeast genome
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Assign data to the yeastGenome object
yeastGenome <- BSgenome.Scerevisiae.UCSC.sacCer3

# Get the first 30 bases of chrM
getSeq(yeastGenome, name = "chrM", end = 30)
```
# Biostrings and When to Use Them?
## Introduction to Biostrings
### Exploring the Zika Virus
It's your turn to explore the Zika virus genome, which has been loaded in your workspace as ```zikaVirus```. The genome was downloaded from NCBI and you can apply ```Biostrings``` functions to learn more about it.

Start by checking the alphabet of the sequence.
```
alphabet() # Shows the letters included in the sequence
alphabetFrequency() # Shows the counts per letter
```
Remember from the video that each alphabet corresponds to a specific biostring container, and each alphabet usually has extra code letters and symbols.
```
# Load packages
library(Biostrings)

# Check the alphabet of the zikaVirus
alphabet(zikaVirus)

# Check the alphabetFrequency of the zikaVirus
alphabetFrequency(zikaVirus)

# Check alphabet of the zikaVirus using baseOnly = TRUE
alphabet(zikaVirus, baseOnly = TRUE)
```
### Manipulating BioStrings
Using a short sequence (```dna_seq```) from the ```zikaVirus``` object, it is your turn to have a go with the two biological processes of transcription and translation.

In the first two parts of this exercise, you will apply these processes separately. In the last part, you'll apply them in one step.

You'll be using a very small sequence in this exercise, but remember that the power of ```Biostrings``` comes to light when manipulating much larger sequences.

The ```Biostrings``` package has already been loaded for you. Using ```zikaVirus```, you will translate the first 21 characters into an ```AAString```.
```
# Unlist the set, select the first 21 letters, and assign to dna_seq
dna_seq <- subseq(unlist(zikaVirus), end = 21)
dna_seq

# Transcribe dna_seq into an RNAString object and print it
rna_seq <- RNAString(dna_seq) 
rna_seq

# Translate rna_seq into an AAString object and print it
aa_seq <- translate(rna_seq)
print(aa_seq)

# Transcribe and translate in one step the dna_seq into an AAString object and print it
aa_seq <- (translate(DNAString(dna_seq)))
print(aa_seq)
```
## Sequence Handling
### From a Set to a Single Sequence
From the video, you know that sequences can be loaded into R using the function ```readDNAStringSet()```. The ```zikaVirus``` has been read into your environment using this function.

It is your turn to convert this set into a single sequence, explore the new sequence, and subset it.
```
# Create zikv with one collated sequence using zikaVirus
zikv <- unlist(zikaVirus)

# Check the length of zikaVirus and zikv
length(zikaVirus)
length(zikv)

# Check the width of zikaVirus
width(zikaVirus)

# Subset zikv to only the first 30 bases
subZikv <- subseq(zikv, end = 30)
subZikv
```
### Common sequence manipulation functions
So far, you've learned the the most common sequence manipulation functions:
```
reverse()
complement()
reverseComplement()
translate()
```
In real life, you can manipulate really large sequences using these functions.

However, to see the value and the results of these functions in this exercise, you will use a small subset. The object ```zikv```, which you have previously subsetted from the Zika genome, has only 30 bases.
```
# Reverse the zikv sequence
reverse(zikv)

# Complement the zikv sequence
complement(zikv)

# Reverse complement the zikv sequence
reverseComplement(zikv)

# Translate the zikv sequence
translate(zikv)
```
## Why are we interested in patterns?
### Finding Palindromes
It is your turn to find palindromic sequences using the ```zikv``` sequence. Remember, ```findPalindromes()``` can only search a single sequence, and does not work with a set.
```
# Find palindromes in zikv
findPalindromes(zikv)
```
### Finding a conserved region within six frames
Now you will be able to look for the NS5 protein sequence in the Zika virus sequence. The NS5 is a very conserved virus protein. It was downloaded and loaded for you from Uniprot.

The Zika virus DNA sequence has been transcribed into an RNAStringSet, called   ```rnaframesZikaSet```. The set has six reading frames (one per sequence) for you to translate into amino acids. When doing the search, you will set the ```max.mismatch``` argument in your call of ```vcountPattern()``` to add flexibility to your search.
```
# Print rnaframesZikaSet
rnaframesZikaSet

# Translate rnaframesZikaSet
AAzika6F <- translate(rnaframesZikaSet)
AAzika6F

# Count NS5 protein matches in AAzika6F, allowing 15 mismatches
vcountPattern(pattern = NS5, subject = AAzika6F, max.mismatch = 15)

# Subset the frame that contains the match from AAzika6F
selectedSet <- AAzika6F[3] 
  
# Convert selectedSet into a single sequence
selectedSeq <- unlist(selectedSet)
```
### Looking for a match
From the previous exercise, you have two objects: selectedSet(a set) and selectedSeq (a single sequence).

You have recently discovered that pattern ns5 is on frame 3 of the AAzika6F. Now, you will discover what the match looks like using selectedSet and selectedSeq.
```
# Use vmatchPattern() with the set
vmatchPattern(pattern = ns5, subject = selectedSet, max.mismatch = 15)

# Use matchPattern() with the single sequence
matchPattern(pattern = ns5, selectedSeq, max.mismatch = 15)

# Take your time to see the similarities/differences in the result.
```
# IRanges and GenomicRanges
## IRanges and Genomic Structures
### Constructing IRanges
In the video, some ```IRanges``` constructor examples were provided. This is your turn to practice creating sequence ranges with different arguments and see how these arguments are reused or complemented.

Using the ```IRanges()``` function, you can specify parameters such as ```start```, ```end```, or ```width```. These parameter inputs can fall into one of two categories:
```start```, ```end```, and ```width``` are numeric vectors.
The ```start``` parameter is a logical vector.
Missing arguments will be resolved using the equation ```width = end - start + 1```.

The ```IRanges()``` constructor indicates that all of the parameters are optional with default NULL:

```IRanges(start = NULL, end = NULL, width = NULL, names = NULL)```
```
# Load IRanges package
library(IRanges)

# IRnum1: start - vector 1 through 5, end - 100 
IRnum1 <- IRanges(start = 1:5, end = 100)

# IRnum2: end - 100, width - 89 and 10
IRnum2 <- IRanges(end = 100, width = c(89, 10))

# IRlog1: start = Rle(c(F, T, T, T, F, T, T, T))
IRlog1 <- IRanges(start = Rle(c(F, T, T, T, F, T, T, T)))

# Print objects in a list
print(list(IRnum1 = IRnum1, IRnum2 = IRnum2, IRlog1 = IRlog1))
```
### Interacting with IRanges
You can use the ```IRanges()``` function to create a single sequence. You can also provide vectors to ```IRanges()``` to create multiple sequence ranges at the same time. This is both fascinating and useful! The creation of objects ```seq_1``` and ```seq_2``` are examples of this.

For this exercise, check the width of each of the sequences provided here, using ```width()``` and ```lengths()```. Notice the difference between the two types of outputs.

Remember that ```width = end - start + 1```.
```
# Create the first sequence seq_1
seq_1 <- IRanges(start = 10, end = 37)

# Create the second sequence seq_2
seq_2 <- IRanges(start = c(5, 35, 50),
                 end = c(12, 39, 61),
                 names = LETTERS[1:3])

# Check the width of seq_1 and seq_2
width(seq_1)
width(seq_2)

# Check the width of seq_1 and seq_2
lengths(seq_1)
lengths(seq_2)
```
## Gene of Interest
### From tabular data to Genomic Ranges
In the video, you learned ways to create ```GRanges``` objects. You can define a GRange with a range's name, start, and end positions (```seqnames```, ```start```, and ```end```). If you have data in table format, you can also transform it into a ```GRanges``` object. Let's use a data frame, called ```seq_intervals```, as this is most likely where you have stored your sequence intervals. Note: you can also use a ```tibble``` if you are more familiar with them.

You will use the predefined ```seq_intervals``` data frame to transform into a ```GRanges``` object using the ```as()``` function. The ```as()``` function was introduced in the last video - it takes in an object and the name of the class to convert the object to.
```
# Load GenomicRanges package
library(GenomicRanges)

# Print seq_intervals
print(seq_intervals)

# Create myGR
myGR <- as(seq_intervals, "GRanges")

# Print myGR
print(myGR)
```
### GenomicRanges accessors
In the previous exercise, you created a ```GRanges``` object from a data frame containing the basic information. You will discover that ```GRanges``` can store much more!

Use the accessor method to explore the ```GRanges``` object, ```myGR```. You can extract characteristics from a ```GRanges``` object such as chromosome names, the number of sequences, the names of each sequence, information about strand, score, length, and more.

The following are basic accessors for ```GRanges```:
```
seqnames(gr)
ranges(gr)
mcols(gr)
genome(gr)
seqinfo(gr)
```
For a complete list of accessors, you can check ```methods(class = "GRanges")```.
```
# Load GenomicRanges
library(GenomicRanges)

# Print the seqinfo of myGR
print(seqinfo(myGR))

# Check the metadata
mcols(myGR)
```
### Human genome chromosome X
It is your turn to use the ```TxDb.Hsapiens.UCSC.hg38.knownGene``` package and extract information from it. Like in the video, you will subset all genes in chromosome X, using the function ```genes()``` with the argument ```filter``` set to select instances where ```tx_chrom = "chrX"```. Then, you will explore this subset of genes.

Remember that ```filter``` receives a ```list()``` of filter conditions to select specific genome intervals.

If you would like to test other filters, valid names for this list are: ```"gene_id", "tx_id", "tx_name", "tx_chrom", "tx_strand", "exon_id", "exon_name", "exon_chrom", "exon_strand", "cds_id", "cds_name", "cds_chrom", "cds_strand", and "exon_rank".```
```
# Load human reference genome hg38
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Assign hg38 to hg, then print it
hg <- TxDb.Hsapiens.UCSC.hg38.knownGene
hg

# Extract all the genes in chromosome X as hg_chrXg, then print it
hg_chrXg <- genes(hg, filter = list(tx_chrom = c("chrX")))
hg_chrXg

# Extract all positive stranded genes in chromosome X, assign to hg_chrXgp, then sort it
hg_chrXgp <- genes(hg, filter = list(tx_chrom = c("chrX"), tx_strand = "+"))
sort(hg_chrXgp)
```
## Manipulating collections of GRanges
### More about ABCD1
Now that you know that there is an overlap between chromosome X and the gene ```ABCD1```, you can find its gene id and location, also called locus.
```
# Store the overlapping range in rangefound
rangefound <- subsetByOverlaps(hg_chrX, ABCD1)

# Print names of rangefound
names(rangefound)

# Print the gene of interest 
print(ABCD1)

# Print rangefound
print(rangefound)
```
### How many transcripts?
Remember in the video how we discovered the length of the exons in one of the transcripts of our gene of interest? It is your turn to find out how many transcripts the gene ```ABCD1``` has. You can find out by using:

```transcriptsBy(x, by = "gene")```

Once you get all the transcripts by gene, you can then select any gene by its id. The gene id of ```ABCD1``` is ```215```. A little tip to select the information on a specific gene is to use back-tick marks around the gene id, for example ```transcripts$`215````.
```
# Load the human transcripts DB to hg
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
hg <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Prefilter chromosome X "chrX" using seqlevels()
seqlevels(hg) <- c("chrX")

# Get all transcripts by gene and print it
hg_chrXt <- transcriptsBy(hg, by = "gene")
hg_chrXt

# Select gene `215` from the hg_chrXt
hg_chrXt$'215'
```
### From GRangesList object into a GRanges object
The ```unlist()``` operation is fast and serves to partition a ```GRangesList```.

You can unlist the ```hg_ChrX``` and then check how the lengths differ between the ```GRangesList``` and the ```GRanges``` objects.
```
# Unlist hg_ChrX and save result as myGR
myGR <- unlist(hg_ChrX)

# Compare classes of hg_ChrX and myGR
class(hg_ChrX)
class(myGR)

# Compare length of hg_ChrX and myGR
length(hg_ChrX)
length(myGR)
```
# Introducing ShortRead
## Sequence files
### Exploring a fastq file
Fastq files usually contain thousands or millions of reads, and can become very large in size! For this exercise, you will use a small ```fastq``` sub sample of 500 reads, which fits easily into memory and can be read entirely using the function ```readFastq()```.

The original sequence file comes from Arabidopsis thaliana, provided by the UC Davis Genome Center. The accession number is SRR1971253 and was downloaded from the Sequence Read Archive (SRA). It contains DNA from leaf tissues, pooled and sequenced on Illumina HiSeq 2000. These sequences are single-read sequences with 50 base pairs (bp) length.

```fqsample``` is a``` ShortReadQ``` object and contains information about reads, quality scores, and ids. It's your turn to explore it!
```
# Load ShortRead
library(ShortRead)

# Print fqsample
fqsample

# Check class of fqsample
class(fqsample)

# Check class sread fqsample
class(sread(fqsample))

# Check ids of fqsample
id(fqsample)
```
### Extract a sample from a fastq file
It is your turn to draw a sample piece from a sequence of many reads.

You will use the same file you've read in the previous exercise. This file has 500 reads, each of 50 bp. The file path is stored in an object called ```f```.

Using ```FastqSampler(con = file_path, n = length)```, ```set.seed()```, and ```yield()``` you can subset 100 reads from your sequence file.
```
# Load ShortRead
library(ShortRead)

# Set a seed for sampling
set.seed(1234)

# Use FastqSampler with f and select 100 reads
fs <- FastqSampler(con = f, n = 100)

# Generate new sample yield
my_sample <- yield(fs)

# Print my_sample
print(my_sample)
```
## Seqeucne quality
### Exploring sequence quality
It is your turn to perform a quality control check on the ```fqsample```. This is an important step before starting further analyses to quickly identify data problems.

To check the encoding values for each letter in ```quality()```, use ```encoding()```:

encoding(quality(fqsample))
For a quality assessment (QA) summary use qa():
```
qaSummary <- qa(fqsample, type = "fastq", lane = 1)
```
This ```qaSummary``` has already been created for you. QA elements can be accessed with ```qaSummary[["nameElement"]]```, where ```nameElement``` is the name of the element you wish to inspect.
```
# load ShortRead
library(ShortRead)

# Check quality
quality(fqsample)

# Check encoding of quality
encoding(quality(fqsample))

# Check baseQuality
qaSummary[["baseQuality"]]
```
### Try your own nucleotide frequency plot
Now it's time to take a closer look at the frequency of nucleotides per cycle. The best way to do this is by making a visualization. Usually, the first cycles are a bit random, and then the frequency of nucleotides should stabilize with the coming cycles.

This exercise uses the complete ```fastq``` file SRR1971253 with some pre-processing done for you:
```
library(ShortRead)
fqsample <- readFastq(dirPath = "data", 
                      pattern = "SRR1971253.fastq")
# extract reads                      
abc <- alphabetByCycle(sread(fqsample))

# Transpose nucleotides A, C, G, T per column
nucByCycle <- t(abc[1:4,]) 

# Tidy dataset
nucByCycle <- nucByCycle %>% 
  as_tibble() %>% # convert to tibble
  mutate(cycle = 1:50) # add cycle numbers
```
Your task is to make a Nucleotide Frequency by Cycle plot using tidyverse functions!
```
# Glimpse nucByCycle
glimpse(nucByCycle)

# Create a line plot of cycle vs. count
nucByCycle %>% 
  # Gather the nucleotide letters in alphabet and get a new count column
  pivot_longer(-cycle, names_to = "alphabet", values_to = "count") %>% 
  ggplot(aes(x = cycle, y =  count, color = alphabet)) +
  geom_line(size = 0.5 ) +
  labs(y = "Frequency") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank())
```
## Match and Filter
### filtering reads on the go!
What if, from all of the reads in a file, you are only interested in some of those reads? You can use a filter!

Let's say that you are interested only in those reads that start with the pattern ```"ATGCA"```. A tiny filtering function can do the job, making use of the ```srFilter()``` function:
```
myStartFilter <- srFilter(function(x) substr(sread(x), 1, 5) == "ATGCA")
```
This function, which has been created for you, takes a ```ShortRead``` derived object as an input and outputs the reads starting with the pattern ```"ATGCA"```. Let's put this function to use!
```
# Load package ShortRead
library(ShortRead)

# Check class of fqsample
class(fqsample)

# Filter reads into selectedReads using myStartFilter
selectedReads <- fqsample[myStartFilter(fqsample)]

# Check class of selectedReads
class(selectedReads)

# Check detail of selectedReads
detail(selectedReads)
```
### More filtering
Awesome! Now that you've had some practice with filtering reads, let's use the function ```polynFilter()```. This function selects reads that contain less than a given number of duplicate nucleotides. For example, ```polynFilter(threshold = 20, nuc = c("A"))``` will select all reads that contain less than 20 A's. The parameter ```nuc``` is a character vector containing IUPAC symbols for nucleotides or the value "other" for all non-nucleotide symbols.

The ```fqsample``` object is available in your workspace.
```
# Check reads of fqsample
sread(fqsample)

# Create myFil using polynFilter
myFil <- polynFilter(threshold = 3, nuc = c("A"))

# Apply your filter to fqsample
filterCondition <- myFil(fqsample)

# Use myFil with fqsample
filteredSequences <- fqsample[filterCondition]

# Check reads of filteredSequences
sread(filteredSequences)
```
## Multiple assessment
### Plotting cycle average quality
```
# Load package Rqc
library(Rqc)

# Average per cycle quality plot
rqcCycleAverageQualityPlot(qa)

# Average per cycle quality plot with white background
rqcCycleAverageQualityPlot(qa) + theme_minimal()

# Read quality plot with white background
rqcReadQualityPlot(qa) + theme_minimal()
```
