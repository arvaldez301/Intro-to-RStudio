# Multi-Sequence-Alignments

## Step1: Load the necessary packages for sequence analysis:
```
library(msa)
```
If your R verison is not allowing for the download of the msa package please use the following:
```
BiocManager::install("msa")
```
## Step 2: Set Path to the Fasta File
We are using the ```msa``` package to pull example files. To view the files that the ```msa``` package has to over you can use the following code:
```
package_dir <- system.file("extdata", package = "DECIPHER")
files <- list.files(package_dir)
print(files)
```
Once selecting a file from the practice data you can load it using the code below:

```
seq <- system.file("examples", "exampleAA.fasta", package = "msa")
```
If you wish to up load your own file path use the following:
```
seq <- "<<<path to FASTA file>>>"
```
You need to replace ```"<<path to FASTA file>>"``` with the actual path to your FASTA file on your computer if you are using your own

## Step 3: Read the Sequences and Perfrom the Alignment

```
sequences <- readDNAStringSet(seq)
```
Note that you can modify this line depending on the type of sequence you are working with. For example, you can use ```readAAStringSet``` for amino acid sequences or ```readRNAStringSet``` for RNA sequences.

After reading the sequences, it performs multiple sequence alignment (MSA) using the `msa` function:

```
alignment <- msa(sequences)
```

This line performs a default alignment using the `msa` function without specifying a specific alignment method. Alternatively, you can specify the alignment method by adding the `method` argument, such as `method = "ClustalOmega"` or `method = "Muscle"`. This determines the algorithm used for alignment.

To print the alignment, you can use the `print` function:

```
print(alignment)
```

### Trimming sequence if you need to!

#### Using the Decipher Package
```
library(msa)
library(DECIPHER)
```

Set the filepath to the FASTA file you wish to use. For this we are using an example data set from the DECIPHER package

```
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package = "DECIPHER")
```

This line provides an example file path from the ```DECIPHER``` package. You need to replace it with the path to your own FASTA file.

Then, it reads the DNA sequences from the FASTA file using the ```readDNAStringSet``` function from the ```DECIPHER``` package:

```
dna <- readDNAStringSet(fas)
```

The code proceeds to trim the DNA sequences using the ```TrimDNA``` function from the ```DECIPHER``` package:

```
trimmed_dna <- TrimDNA(dna,
                        leftPatterns = "",
                        rightPatterns = "",
                        type = "sequences",
                        quality = NULL)
```

In this example, the ```TrimDNA``` function is used with empty strings for ```leftPatterns``` and ```rightPatterns```. You can provide specific patterns to trim from the sequences if desired. The `type` argument is set to "sequences" to specify that the trimming should be applied to the entire sequences. The ```quality``` argument is set to ```NULL``` assuming that no quality scores are available for the sequences.

Finally, the trimmed DNA sequences are printed:

```
print(trimmed_dna)
```

#### Using the Biostrings Package

```
library(Biostrings)
```

Then, it sets the path to a FASTQ file containing DNA sequences:

```
fpath <- system.file("extdata", "s_1_sequence.txt", package = "Biostrings")
```

This line provides an example file path from the ```Biostrings``` package. Again, you need to replace it with the path to your own FASTQ file.

The code reads the quality-scaled DNA sequences from the FASTQ file using the ```readQualityScaledDNAStringSet``` function from the ```Biostrings``` package
