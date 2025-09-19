# Building a Phylogenetic Tree

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
accessions <- c("____", "____", "______")
```
We will be using the following accession numbers so please copy and paste the following code into your markdown file.

```
accessions <- c("DQ673917.1", "MH253817.1", "KF384926.1", "MH260579.1", "JAATJU010022200.1", "BI431008.1", "NC_027278.1", "U39315.1", "HM102293.1", "HM102307.1", "MK773649.1", "JF459178.1")
```

## STEP 2: Function to fetch sequences and scientific names
The next step is to pull all of that information out of NCBI that we want. The first part of the code below will pull out the FASTA File, then will pull out the scientific name that is registered with the accession number, and finally will rename the pulled information to also include the accesstion number. THis cam be modified however you would like, so if you dont want to include the accession number in your tree you would could you the alternate code. Take a look to see how the codes differ!

**Scientific Name and Accession Number**

* ```entrez_fetch``` pulls the raw sequence data from NCBI
* ```entrez_summary``` retrieves metadata, such as the organism name.
* ```DNAStringSet``` allows efficient handling of sequences in R.
* ```Reduce(c, seqs_list)``` combines multiple DNAStringSet objects into one.

The function ensures each sequence has a unique, descriptive name by combining the scientific name with its accession.

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

**Scientific Name Only**
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
    
    # Use only the scientific name as the label
    new_name <- sci_name
    names(seq) <- new_name
    name_map[acc] <- new_name
    
    seqs_list[[acc]] <- seq
  }
  
  list(sequences = Reduce(c, seqs_list), name_map = name_map)
}
```

## STEP 3: Fetch sequences and names
Here we use the custom function fetch_seqs_and_names() (built in Step 2) to retrieve DNA sequences and their corresponding scientific names from GenBank.

```data$sequences``` contains the DNA sequences in DNAStringSet format
```data$name_map``` stores the mapping between accession numbers and species names

```
data <- fetch_seqs_and_names(accessions)
seqs <- data$sequences 
name_map <- data$name_map
```

## STEP 4: Align sequences
From the library ```DECIPHER```, this function is used to generate a high-quality alignment of COI gene sequences prior to tree construction. Multi Sequence Alignment ensures that homologous nucleotides are compared across species.

```
aligned <- AlignSeqs(seqs)
```

## STEP 5: Create phylogenetic data
Convert the aligned sequences into a ```phyDat object``` (from the ```phangorn``` package). This is a special format required for phylogenetic analysis.

```
phydat <- phyDat(as.matrix(aligned), type = "DNA")
```

## STEP 6: Distance matrix and initial NJ tree
Moving on to conducting a maximum liklihood and Neighbor Joining analysis. Each serves a specific purpose and is important to the construction of our tree. Below I have provided short explanatins regarding what each of these methods do.

**Maximum Likelihood**
* Statistical apporach
    * Finds the tree that most likely explains the observed sequence data given a model of evolution
* Model Based
  * Incorporates subsitution models to account for different rates of evolution among sites
* Branch Support
  * Combined with bootstrapping to assess confidence in tree branches
* Accurate but computationally intensive
  * Provides high-resolution trees but can be slow for large datasets
 
**Neighbor-Joining (NJ)**
* Distance-based method
    * Builds a tree using pairwise distances between sequences
* Fast and efficient
  * Good for large datasets where computational speed is important
* Produces unrooted trees
  * Typically requires an outgroup to root the tree. We will be doing this later on!
* Less model-dependent
  * Simpler than ML by may be less accurate for complex evolutionary scenarios

```dist.ml()``` → calculates a maximum-likelihood (ML) based distance matrix for all sequences.

```NJ()``` → builds an initial Neighbor-Joining (NJ) tree based on those distances.

```
dm <- dist.ml(phydat)
treeNJ <- NJ(dm)
```

## STEP 7: Bootstrap support (500 replicates)
Bootstrapping is a statistical method used to assess the reliability of branches in a phylogenetic tree. It involves repeatedly resampling the sequence alignment with replacement to create many pseudo-replicate datasets. For each replicate, a tree is reconstructed, and the frequency with which a particular branch appears across all replicates is calculated as the bootstrap support value. High bootstrap values (typically ≥70%) indicate strong support for that branch, while low values suggest uncertainty. Bootstrapping provides a way to quantify confidence in the inferred relationships and helps researchers distinguish well-supported clades from more tentative ones.

In short, we will be doing a bootstrap with 500 replicates and remember, higher bootstrap values indicate a stronger level of confidence in that clade
  
```
bs <- bootstrap.phyDat(phydat, FUN = function(x) NJ(dist.ml(x)), bs = 500)
```

## STEP 8: Create tree with bootstrap values, without plotting
Attach the bootstrap values to the tree topology.

* ```plotBS()``` adds bootstrap values to nodes, but with ```plot = FALSE``` we keep the tree object without drawing it.

```
tree_with_bs <- plotBS(treeNJ, bs, p = 50, type = "phylogram", plot = FALSE)
```

## Step 9: OutGroup
Root the tree using an outgroup (a distantly related taxon). This outgroup provides a reference point to determine the direction of evolutionary change, allowing the tree to be properly rooted. Some important points that an outgroup can help you determine are:
* Helps distingues ancestral traits from derived (apomorphic) trains in the ingroup
* Clarifys relationships by providing a comparison outside the group of interest. Outgroups help resolve branching order among ingroup taxa
* They help to avoid misinterpretation by preventing misidentification of trait evolution patterns due to lack of context
* It also serves as a benchmark for divergence. It serves as a baseline to estimate relative evolutionary distances among ingroups
  
For this example all sequences beginning with "Haliotis" are used as outgroup taxa. ```resolve.root = TRUE``` ensures the tree is properly rooted.
  
```
out_tips <- grep("^Haliotis", tree_with_bs$tip.label, value = TRUE)
rtree <- root(tree_with_bs, outgroup = out_tips, resolve.root = TRUE)
```

## STEP 10: Create and plot the tree
Now lets visualize our tree! For this we will be utilizing ```ggtree```. Everything is fully customizable and can be altered to fit your preferences. Some aspects of the ```plot.margin``` function will need to be adjusted depending on the amount of seqeunces you are including as the tree can become clutered. Playing around with the numbers and replotting the tree is a great way to learn what each function does!

* ```layout = "rectangular"``` → standard phylogram layout.
* ```geom_text2()``` → adds bootstrap labels above branches.
* Axis scaling and spacing adjustments ensure legibility (especially for large trees).
* The tree is formatted in Times New Roman, size 12, for publication-ready figures.
  
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
