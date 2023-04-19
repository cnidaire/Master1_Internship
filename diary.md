# Day 1

Had a tour of the building and joined the different group chats, read bibliography and the is preparing for the meting at 4:30 where Elena presented the weekly paper. It was interesting.

I received the keys from Christan and I will have access to the computer.  I will have a meting with Christan and Molly on 28 so it would be nice if I knew a bit more about what they talk because they expect me to participate and I know nothing for the moment.



# Day 2

continue to work on understanding the Dada2 pipeline

It might be a good idea to pout the different parts in different files to makes modules, instead of a huge code.

I filtered the un-annotated sequences, It doesn't seem like there is any major mistake.

I tried to convert the csv file into a fasta file by making my own code because I didn't succeed to run Christian one and I think it's better to avoid for loops. But I two main problems, first I get some spaces between > and the number and the more sequences I add, the more spaces I get and the second is I have a space before the 2 and I wonder why and if I should remove it manually.



# Day 3

Continue debugging the csv to fasta code. And stopped using the basic R and used Tidyverse, it's faster (less or around a second) and is way cleaner to me.

I'm now trying to download the data from ENA and NCBI and to really understand Dada2



## Evaluation of Primers Targeting the Diazotroph Functional Gene and Development of NifMAP – A Bioinformatics Pipeline for Analyzing nifH Amplicon Data

Marker gene for the Nitrogen fixation : *dinitrogenase NifH* is genetically well conserved.

Here: presentation of a pipeline for processing NifH amplicon datasets: NifMAP(“NifH MiSeq Illumina Amplicon Analysis Pipeline”).

Then choose which primer is the most adapted to our study. (for NifMAP, we used 4 forward primers and 4 reverse primers)

The fixation of nitrogen among the bacteria is quite old and hence it is found among different branches. There is 20 functionnal and regulatory genes organized in several operons termed together Nif regulons. Among theme, the most conserved is NifH coding for the dinitrogenase reductase.

NifH form four clusters:

1. canonical MoFe nitrogenase
2. nearly all sequences of the alternative vnfH nitrogenase + sequences of the second alternative nitrogenase, anfH
3. anaerobic bacteria and archaea such as methanogens, spirochetes, sulfate reducers, non-sulfur purple bacteria, green sulfur bacteria, acetogens and Clostridia
4.  most of the known cluster IV nitrogenases do not encode for a protein involved in N2 fixation

PCR-based surveys continue to serve as an important tool for studying the diazotroph diversity and, if the nifH transcripts are targeted (via cDNA), can also reveal transcription patterns in the environment.

### Construction of a Hidden Markov Chain Model

We will filter out the non-NifH genes and align the NifH genes.

### NifMAP – An Automated Pipeline for Analyzing nifH Amplicon Reads

raw reads ->  merge reads and quality filter -> filter the seq using a NifH HMM -> generate OTUS -> translate and correct seq using a framebot (translation to AA and correct eventual frameshift) -> filter homologous genes using a HMM (only OTU representative that scored the highest with the NifH model) -> classify OTUs using blast, or other





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