# BIOCONDUCTOR
Bioconductor (https://www.bioconductor.org/) is an open source software for bioinformatics that is full of a variety programs that can help with the analysis of biological data. For a full list of the packages available please visit the following: https://www.bioconductor.org/packages/release/BiocViews.html#___Software
Further information on the other packages that we will be using can also be found at the link.

# Building your Tree
Before we get started we need to install a few packages. Bioconductor requires an additional software to help with the install process, called BiocManager. You will only need to do this install once and then can load each of the librarys as normal.
```
install.packages(c("rentrez", "BiocManager"))
BiocManager::install(c("DECIPHER", "Biostrings", "phangorn", "treeio", "ggtree"))
```

The following libraries will be used to build our tree. Each libary serves a specific function:
* ```rentrez``` Used directly to retrieve nucleotide sequenes from NCBI GenBank via accession numbers. This automates sequence aquisition without downloading FASTA files.
* ```DECIPHER``` Part of the Bioconductor Software. Provides functions for multiple sequence alignment.
* ```Biostrings``` Part of the Bioconductor Software. Provides efficient data structures and functions for handling large DNA and protein sequence datasets. It ensures the sequences are stored and manipulated properly for downstream analyses.
* ```phangorn``` Used for phylogenetic reconstruction and statistical support estimation.
* ```treeio``` Facilitates the manipulation and storage of phylogenetic trees with associated data, ensuring compatibility between phylogenetic objects and visualization tools.
* ```ggtree``` Provides advanced visualization of phylogenetic trees within the ```ggplot2``` framework. This package was used to produce publication-quality trees, add bootstrap values to nodes, and adjust label formatting for readability.
* ```ggplot2``` The general plotting system in R that underlies ```ggtree```. It provides fine control over graphical elements, ensuring that the tree layout and labels are clear and customizable.

```
#Load each of these libaries prior to beginning your code
library(rentrez)
library(DECIPHER)
library(Biostrings)
library(phangorn)
library(treeio)
library(ggtree)
library(ggplot2)
```

## Step 1: Define your accession numbers
Now that all of our libraries are loaded, we can begin determining our accesion numbers and logging those into our code. Within NCBI (https://www.ncbi.nlm.nih.gov/), you can do a search of the organsisms that you are looking. Each organism logged has a unique number that identifies them. This number will be used in our code and will let R know what sequences to pull. For this code, you can include as many accession numbers as you would like, just make sure that each number is in quotations and seperated by a comma.

```
accessions <- c("HQ013317.1", "KX233874", "DQ297506.1")
```

## STEP 2: Function to fetch sequences and scientific names
The next step is to pull all of that information out of NCBI that we want. The first part of the code below will pull out the FASTA File, then will pull out the scientific name that is registered with the accession number, and finally will rename the pulled information to also include the accesstion number. THis cam be modified however you would like, so if you dont want to include the accession number in your tree you would remove ```(", acc, ")``` from the code.
```
fetch_seqs_and_names <- function(accessions) {
  seqs_list <- list()
  name_map <- character()
  
  for (acc in accessions) {
    # Fetch sequence in FASTA format
    fasta <- entrez_fetch(db = "nucleotide", id = acc, rettype = "fasta")
    writeLines(fasta, paste0(acc, ".fasta"))
    seq <- readDNAStringSet(paste0(acc, ".fasta"))
    
    # Get scientific name
    summary <- entrez_summary(db = "nucleotide", id = acc)
    sci_name <- summary$organism
    
    # Rename sequence with scientific name + accession to ensure uniqueness
    new_name <- paste0(sci_name, " (", acc, ")")
    names(seq) <- new_name
    name_map[acc] <- new_name
    
    seqs_list[[acc]] <- seq
  }
  
  list(sequences = Reduce(c, seqs_list), name_map = name_map)
}
```

## STEP 3: Fetch sequences and names
```
data <- fetch_seqs_and_names(accessions)
seqs <- data$sequences
name_map <- data$name_map
```

## STEP 4: Align sequences
From the library ```DECIPHER```, this function is used to generate a high-quality alignment of COI gene sequences prior to tree construction.
```
aligned <- AlignSeqs(seqs)
```

## STEP 5: Create phylogenetic data
```
phydat <- phyDat(as.matrix(aligned), type = "DNA")
```

## STEP 6: Distance matrix and initial NJ tree
```
dm <- dist.ml(phydat)
treeNJ <- NJ(dm)
```

## STEP 7: Bootstrap support (500 replicates)
```
bs <- bootstrap.phyDat(phydat, FUN = function(x) NJ(dist.ml(x)), bs = 500)
```

## STEP 8: Create tree with bootstrap values, without plotting
```
tree_with_bs <- plotBS(treeNJ, bs, p = 50, type = "phylogram", plot = FALSE)
```

## Step 9: OutGroup
```
out_tips <- grep("^Haliotis", tree_with_bs$tip.label, value = TRUE)
rtree <- root(tree_with_bs, outgroup = out_tips, resolve.root = TRUE)
```

## STEP 10: Create and plot the tree
```
gtree <- ggtree(tree_with_bs, layout = "rectangular", branch.length = "branch.length")

gtree$data$y <- gtree$data$y * 30

gtree$data$label_y <- gtree$data$y + runif(nrow(gtree$data), -0.5, 0.5)  # jitter y pos

p <- gtree +
  geom_text2(aes(x = x, y = label_y, label = label), size = 3, hjust = -0.3) +  # bootstrap labels above nodes
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +  # vertical padding at top/bottom
  xlim(0, max(gtree$data$x, na.rm = TRUE) + 0.5) +  # leave space for long labels
   theme_tree2(base_size = 12, base_family = "Times New Roman") +
  theme(
    text = element_text(family = "Times New Roman", size = 12),
    plot.margin = margin(2, 6, 2, 1, "cm")
  )

print(p)
```
