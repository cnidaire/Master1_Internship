

## DADA2

https://www.bioconductor.org/packages/release/bioc/vignettes/dada2/inst/doc/dada2-intro.html#overview-of-the-dada2-pipeline

```R
library(dada2); packageVersion("dada2")

fnF1 <- system.file("extdata", "sam1F.fastq.gz", package="dada2")
fnR1 <- system.file("extdata", "sam1R.fastq.gz", package="dada2")
filtF1 <- tempfile(fileext=".fastq.gz")
filtR1 <- tempfile(fileext=".fastq.gz")
```

### 1 Filter and Trim

check quality:

```R
plotQualityProfile(fnF1) # Forward
plotQualityProfile(fnR1) # Reverse
```

```R
filterAndTrim(fwd=fnF1, filt=filtF1, rev=fnR1, filt.rev=filtR1,
                  trimLeft=10, truncLen=c(240, 200), 
                  maxN=0, maxEE=2,
                  compress=TRUE, verbose=TRUE)
```

you don't need to the forward and the reverse reads at the same position, resulting in reads of different length.

The FilterAndTrim() function filter forward and reverse reads jointly

- *trimLeft=10* unable to remove the first 10 nucleotides of each reads(as there usually is a drop for the first nucleotides)
- *truncLen=c(240, 200)* trunk the forward at nucleotide 240 and the backward at nucleotide 200
- *maxN=0* we filter out all the reads with more than 0 ambiguous nucleotides

If using a pair-end sequencing data, must have at the very least 20 nucleotides that overlap after the trimming. 



### 2 Dereplicate

```R
derepF1 <- derepFastq(filtF1, verbose=TRUE)
derepR1 <- derepFastq(filtR1, verbose=TRUE)
```

Condense the data by collapsing together all reads that encode the same  sequence, which significantly reduces later computation times

derepFastq maintain a summary of the quality information for each dereplicated sequence in $quals



### 3 Learn the errors rates

```R
errF <- learnErrors(derepF1, multithread=FALSE) # multithreading is available on many functions
errR <- learnErrors(derepR1, multithread=FALSE)
```



### 4 Infer sample composition

```R
dadaF1 <- dada(derepF1, err=errF, multithread=FALSE)
dadaR1 <- dada(derepR1, err=errR, multithread=FALSE)
print(dadaF1)
```

 

### 5 Merge forward/reverse reads 

We’ve inferred the sample sequences in the forward and reverse reads  independently. Now it’s time to merge those inferred sequences together, throwing out those pairs of reads that don’t match.

It will return a data frame corresponding to each successfully merged sequences.

```R
merger1 <- mergePairs(dadaF1, derepF1, dadaR1, derepR1, verbose=TRUE)
```



### 6 Remove chimeras

when sequence incompletely amplified, a chimera is formed and will result in the next generation in half of this amplicon and half of an other. Hence we have to remove those sequences.

```R
merger1.nochim <- removeBimeraDenovo(merger1, multithread=FALSE, verbose=TRUE)
```



### 7 A second sample

repeat the same steps for an other sample, pretty straight forward.



### 8 Create a sequence table

If we want to combine all the inferred samples into one unified table: give a matrix.

Each row is a processed sample and each column is a non-chimeric inferred sample sequence.

```R
seqtab <- makeSequenceTable(list(merger1, merger2))
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
dim(seqtab.nochim)
```



### 9 further analysis

They recommend the "*phyloseq*" package: a tool to import, store, analyze, and graphically display complex  phylogenetic sequencing data that has already been clustered into  Operational Taxonomic Units.

I'm not sure if it is really useful in our case.



# Notes from Mar

## Goal of metabarcoding

Characterization of the taxonomic diversity inhabiting various ecosystems using direct environmental DNA

We can use the full genome, hence, we have to use 16S region (~200-400bp): a marker for bacterial and Archeal identification. But there is PCR and sequencing (Illumina) errors in the sequences.

<img src="./images/overview.png" alt="10" style="zoom:50%;" /> 

Pair-end sequencing is usually better, as max read of illumina is 300bp, we obtain 2 x 300bp.

We pool all the sample together(multiplexing and add an index/barcode to the sequence so that we can then separate the samples: demultiplexing). It cost less, is faster and we can avoid as much batch effect as possible.

==Demultiplex:== using the barcodes

==Filter==: remove bad quality sequencing reads

==Cluster(OTUs):== one cluster = sequence identity > 97%. but the can be some ambiguous OTU assignement when a sequence can match in more than one OTU. And there is a centroid sequence taht is the sequence that minimize the sum of the distances to the other sequences of the cluster (commonly the most abundant one). One of the problem is taht it tends to overestimate the bidiversity and is efficient if there is not much errors, which is not the case for PCR and Illumina sequencing because we have chimera for exemple. A cut-off of 98.5% or even a bit high (98.7%) would be better and is used in bacterial identification.



 ==Taxonomic assignation:== 





# [Youtube video from Brown uni](https://www.youtube.com/watch?v=wV5_z7rR6yw)

## before putting data in dada2

remove barcodes and adaptators (cutadapt or Trimmomatic)

samples should bedemultiplex into individual fastQfiles (Idemp)

put reverse and forward reads in the same order

### Check the overlap

if nearly fully overlap, you can agressivly trim the read to avoid bad sequences

## Dada2 workflow

<img src="./images/workflow.png" style="zoom:40%;" />

### Load Packages

```R
library(dada2)
library(DECIPHER) # sequence alignement that will be needed to construct the phylogenic tree
library(phangorn) # to construct the phylogenic tree
library(ggplot2) # data visualization and analysis (at the end)
library(phyloseq) # data visualization and analysis (at the end)
```



### Raw FASTQ Files

```R
# Specify path where FASTQ files are located
path <- "~/DADA2_Tutorial" # change to the diretory containing the FASTQ files after unzipping
list.files(path) # list all the reads. Usually, R1 files contains the forward reads and R2 the corresponding resverse reads.

# Sort the files to ensure reverse and forward reads are in the same order
fnFs <- sort(list.files(path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq"))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1) # here I don't get the end of he code

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
```



### Trim and Filter Reads

Plot the quality score and determine the trimming cut-offs (usually 30 is the lowest accepted quality score but I will makr that we only keep data with really high quality)

We don't have to check for every file (just a couple), there is usually not much variation from sample to sample.

Usually, the quality of the reverse sequence is rather bad compared to the forward, that's why we will have to trimm it more. However, we have to make sure that the Forward and reverse sequences are overlaping (at the very least 20 bp for DADA2 to successfully merge sequences).

```R
# Create a new file path to store filtered and trimmed reads
filt_path <- file.path(path, "filtered") # place the filtered files in a "filtered" subdirectory

# Rename filtered files
filtFs <- files.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filsRs <- files.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Quality Filtering and trimming
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), # fnFs: where F reads/filtFs: where put filtered
                    maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, # MaxN: max ambiguous base allowed/ maxErrorEstimator
                    compress=TRUE, multithread=TRUE) # keep compress if the files should be zipped

head(out)
```

==maxEE== can be lowered if the quality of the reads is good or increased if the quality is bad but the samples are preciousand you can't discard the bad ones.

==truncQ== get rid of the read containing a quality score of 2 (~63% chance of a base call being incorrect), just filter out the worst samples.

### Error Rate estimation

To start, makes the assumption you have the maximum error rate in your sample (takes your most abundant sequence and assumes that's the only true sequence and that all the others are caused by errors).

```R
# Estimate the error model for DADA2 algorithm using reverse reads
errR <- learnErrors(filtRs, multithread=TRUE)
# do the same on the reverse read and it will take more cycles because the reverse reads are worse

# Plot error rates for all possible bases transitions
plotErrors(errF, nominalQ=TRUE) # black:observed error rate, red:expected
```

In general, for the error plot, the frequency of the errors rate decrease as the quality score increase.



### Dereplicate Reads

So that DADA2 doesn't have to work on every single read we have to speed up and simplify the computation

```R
# Dereplicate FASTQ files to speed up computation
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```

### Sample Inference with DADA2 Algorithm

testing the null-hypothesis that the sequence is too abundant in the sample to be solely explained by errors in the data set.

-> low p-value sequence can be considered as real sequences that are not caused by random errors.

If the sequence has a high p-value, it won't be kept 

```R
# Apply core sequence-variant inference algorithm to reverse reads
dadaFs <- dada(derepFs, err=errR, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

==Now all the Forward and Reverse reads have been denoised== we can finally merge all the forward and reverse sequences.

### Merge Reads

```R
# Merge denoised reads

```



### Chimera Checking and Removal

### Assign Taxonomy

### Construct Phylogenic tree

### Import into Phylliseq

### Analysis and visualization





## Setting the environment

### Store the paths

```R
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names=TRUE))
# List the files names located at path containing the pattern and gives the full name (full path)
# and all these names are stored in fnFs for R1 and fnFr for R2
```

### Retrieve sample names

```R
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
# basename retrieve the name before the last separator "/"
# strsplit split by removing the "_"
```

### Create a folder "filtered" with the file names

```R
filtFs <- file.path(path, "Filtered", basename(fnFs))
# will create a new folder "Filtered" withe all the file names inside (the files are empty).
```

### Set working directory

##  Phyloseq objects

### Abundance table

```R
otu_table()
# OTU abundance in the different samples
```

### Metadata

```R
sample_data()
# contains group, treatments, localisations, sample names, env data
```

### Taxonomy classification

```R
tax_table()
# devide each taxa in taxonomic ranks: Kingdom, Phylum, Class, Order, Family (and the sample name of course)
```

### Phylogenic tree

```R
phy_tree()
# number of nodes, number of tips
```

### DNA sequence of ASV/OTU

```R
ref_seq
# gives the reference sequences I think
```





## Table with abundance and taxonomic assignement

we have to fuse two tables: otu_table() and tax_table()

```R
tableau=cbind(t(otu_table(Myselection2)),tax_table(Myselection2))
# t() is to transpose the table
# but you don't always need to transpose when you are working on microbes data apparently
```



```R

```



```R

```



```R

```



```R

```



```R

```



```R

```



```R

```



```R

```



```R

```



```R

```



```R

```



```R

```



```R

```



```R

```



```R

```



```R

```



```R

```



```R

```



```R

```



```R

```



```R

```



```R

```



```R

```



```R

```

