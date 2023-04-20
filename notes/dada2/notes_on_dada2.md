

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