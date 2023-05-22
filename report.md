# Format

Abstract, Introduction, Methodology/Results, Discussion, References



# Why nifH



Nitrogen fixation is catalyzed by the enzyle nitrogenase. Nitrogenase  is composed of two sub-units metabolloproteins: molybdenum iron protein and the iron protein. Likely an ancient protein ,since it's distributed widely through bacteria and archea. Both of the nitrgenase prot are highly conserved but the Fe protein encoded by nifH is the most highl conserved.



nif genes are a collection of genes involved in nitrogen fixation

N2 fixation: role in global biogeo chemical cycling of N, interlinked cycling of C

provide primary natural source for nitrogen ecosystems

Nitrogenase: enzyme that catalyze N2 fixation

Nitrogen availability fuels the biological carbon pump in otherwise nitrogen-limited systems such as subtropical gyres (Karl et al. 2012), and hence boost the ability of the oceans to cope with excess CO2

reduction of atmospheric N2 to ammonia is important to maintain fertility of ocean/ to support primary organic matter production (CO2 fixation for example)

understanding the balance of N cycle has implication to understand past current and future foodwebs and role for marine N2 fixation in sequestration of atmospheric CO2 and the production and consumption of other greenhouse gases such as nitrous oxide.

Nitrogen fixation is important because many living organisms are unable  to metabolize directly the atmospheric nitrogen and would require the  nitrogen fixation capability of certain bacteria in order to produce a  form of nitrogen (e.g. ammonia) that can be readily utilized.



## How is the data generally collected (metabarcoding)?

explain the primers





# Methodology

- explain metabarcoding?
- Present DADA2 Pipiline
- Present NCBI
  - The American National Center for Biotechnology Information science and health by providing access to and genomic information.
  - Well known for GenBank, Pubmed and BLAST (sequence alignment program)
  - there is different databases available, the ones we will be interested in are:
    - gene:
      - Gene integrates information from a wide range of species. A record may  include nomenclature, Reference Sequences (RefSeqs), maps, pathways,  variations, phenotypes, and links to genome-, phenotype-, and  locus-specific resources worldwide.
    - proteins:
      - The Protein database is a collection of sequences from several sources,  including translations from annotated coding regions in GenBank, RefSeq  and TPA, as well as records from SwissProt, PIR, PRF, and PDB. Protein  sequences are the fundamental determinants of biological structure and  function.
      - It is really interesting in our case because the nifH gene is coding for a protein and hence we can access to the taxonomy and many other interesting informations on different databases. The downside is that we obtain the Amino Acid sequences and not the DNA one and we have duplicates
    - Identical Protein Groups (IPG):
      - take the best of both worlds: we obtain all the sequences clustured when these are identical, we obtain the DNA sequences and it's gathering the information from different databases such as Swissprot, GenBank, PIR, etc.
      - downside: Entrez is not really suitable for it
- Present entrez (Global Query Cross-Database Search System)
  - entrez is a searching tool provided by NCBI enabling to the NCBI databases
  - indexing and retrieving system having data from various sources. Is mostly used to integrate information from different sources databases and  formats into a uniform information model and retrieval system which can efficiently retrieve that relevant references, sequences, and structures.
  - can access information with a query







# DADA2 pipeline

## Workflow

**Raw FASTQ file -> Trim and filter the reads -> Error rate estimation -> Sample inference with DADA2 algorithm -> Merge reads -> Chimera Checking and removal -> Assign taxonomy**



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

Usually, the quality of the reverse sequence is rather bad compared to the forward, that's why we will have to trim it more. However, we have to make sure that the Forward and reverse sequences are overlapping (at the very least 20 bp for DADA2 to successfully merge sequences).

We usually remove the first 10 nucleotides with *trimleft=10* because there is a drop of quality and because it contains the primers that will be annoying for later. 

```R
# Create a new file path to store filtered and trimmed reads
filt_path <- file.path(path, "filtered") # place the filtered files in a "filtered" subdirectory

# Rename filtered files
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

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

Condense the data by collapsing together all reads that encode the same sequence so that DADA2 doesn't have to work on every single read we have to speed up and simplify the computation

```R
# Dereplicate FASTQ files to speed up computation
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```

### Sample Inference with DADA2 Algorithm

Testing the null-hypothesis that the sequence is too abundant in the sample to be solely explained by errors in the data set.

-> low p-value sequence can be considered as real sequences that are not caused by random errors.

If the sequence has a high p-value, it won't be kept 

```R
# Apply core sequence-variant inference algorithm to reverse reads
dadaFs <- dada(derepFs, err=errR, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

==Now all the Forward and Reverse reads have been denoised== we can finally merge all the forward and reverse sequences.

### Merge Reads

We’ve inferred the sample sequences in the forward and reverse reads independently. Now it’s
time to merge those inferred sequences together, throwing out those pairs of reads that don’t match.
It will return a data frame corresponding to each successfully merged sequences.

```R
# Merge denoised reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE) # will only merge perfectly overlapped seq so can addmaxMismatch 1 or 2 if lots of reads are not merging

# Inspect the merger data.frame from the first sample
head(mergers[[1]])



# Tabulate Denoised and Merged data
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# View the total length of al total RSVs (Ribosomal Sequence Variants) (close to an OTU table)
table(nchar(getSequences(seqtab))) # gives how sequences there is for the length of the sequence
```



### Chimera Checking and Removal

Chimera: fusion of two or more parents sequences

Perform multiple sequence alignment from the least abundant read and for all the more abundant read, it will do sequence alignment with all possible combinations. When a chimera is detected, it is removed from the sequence table.

```R
# Perform de novo chimera sequence detection and removal
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Calculate the proportion of the non-chimeric RSVs (reads)
sum(seqtab.nochim)/sum(seqtab)
```

More than 90% of the unique sequences identified are bimeras and most of the total reads shouldn't be identified as chimera (less than 70% is really bad).

If too much chimera, check if trim the 20 first base pairs of the reads because contain primers that can artificially increase the number of chimera. If not, try trimming more the low quality bp.

### Assign Taxonomy

```R
# Assign taxonomy using RDP database (greengenes and Silva also available)
# This is performed in two steps: this first one assigns phylum to genus
taxa <- assignTaxonomy(seqtab.cochim,"~DADA2_Tutorial/Taxonomy/rdp_train_set_16.fa.gz", multithread=TRUE) # database we want to use
unname(head(taxa))

# Assign species (when possible)
system.time({taxa.plus <- addSpecies(taxa, "~DADA2_Tutorial/Taxonomy/rdp_species_assignement_16.fa.gz", verbose=TRUE)
colnames(taxa.plus) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
unname(taxa.plus)}) # Usually hard to go down to species level assignment with V4 region (hypervariable region)
```

==using a reference database==

## Reference database



The current database was first created in 2014 and updated until 2017 by the Zehr Lab (Heller 2014). It was then normalized and kept updated manually by Molly Moynihan (put publication). My goal would be to improve the database in order to identify more sequences coming from metabarcoding. The main problems I see in the curent data base are:

1. created manually: possibility of typo
2. not updated automatically: if some sequences attribution change, will not be updated
3. the taxonomy changed: they had to modify the 18,700 sequences manually: lot of work and possible mistakes
4. a lot of identical sequences and hence we hope that the sequences have the exact same taxonomy associated otherwise (normally, there is an algorithm to check it and automatically annotate it in NCBI implemented in 2001 but as it's created manually, there is a possibility of typo) , randomly chosen by DADA2
5. hard to know which sequences are not in the reference yet

That' s why I think building a new database from scratch would be a better idea than improving and building on top of something we are not fully confident in.



How did I proceed?

1. try to understand how the DADA2 pipeline works by re-running it on some already obtained data
2. tried to quantify how much sequences are non annotated at the different levels
3. try to annotate the sequences by using blast on the un-annotated ones to obtain some aligned sequences that I can then put into the reference database bu three problems:
   1. To obtain a full database, you will have to re-run it every time on the un annotated sequences and then incorporate it to the database but it's not what they asked: there will be a lot of different databases annotated by different people
   2. Blast and DADA2 alignment method are not identical, hence it's possible that the sequence was already inside of the reference and that every time we run it, we add a duplicate.
   3. Molly used a database and improved it but the link to the original database is dead so it's not even possible to know how it was made in the first place
4. try to create a database from the DNA sequences by making a query request and then obtaining the right portion of the DNA sequence but soon realize that there is (once we get rid of all the un-annotated ones) 7576 sequences in the nucleotide between 100 and 350 base pairs (to obtain the partial CDS too).
   1. don't exactly remember why I stoped, maybe I should start again
5. Use something called Identical protein groups and which is pretty cool and perfectly suited for what I want but there is few people using it and obtaining the Amino Acids sequence is easy by the DNA sequence is way harder because there is no package to do it directly so I ad to tinker:
   1. get all the ID of the IPG reports
   2. Download the IPG reports
   3. extract from these IPG reports the nucleotide accession number (often linking to he whole genome), the stand, the beginning of the sequence and the end of the sequence
   4. then download so custom GeneBank report containing the taxonomy, other info but only the part of the genome I'm interested in
   5. extract the taxonomy and the sequence for each file to put it into the reference file in the fasta format
   6. **realize that even though I limited the length of the sequences, I have some that are way longer than what they are supposed to and that for some report, inside of these reports, all the sequences don't have the same length even if they are supposed to be identical.**









# Structure

1. Why is nifH stduy interesting for (explain metabarcoding)
   1. gene size, variable like the ribosome, genetically conserved between the species
   2. use intro seen previously
2. How do we annotate the data we obtain (DADA2 pipeline)
3. How was the database created and what are the problems with this db
   1. created manually: possibility of typo
   2. not updated automatically: if some sequences attribution change, will not be updated
   3. the taxonomy changed: they had to modify the 18,700 sequences manually: lot of work and possible mistakes
   4. a lot of identical sequences and hence we hope that the sequences have the exact same taxonomy associated otherwise, randomly chosen by DADA2
   5. hard to know which sequences are not in the reference yet
4. How did I proceed?
   1. try to understand how the DADA2 pipeline works by re-running it on some already obtained data
   2. tried to quantify how much sequences are non annotated at the different levels
   3. try to annotate the sequences by using blast on the un-annotated ones to obtain some aligned sequences that I can then put into the reference database bu three problems:
      1. To obtain a full database, you will have to re-run it every time on the un annotated sequences and then incorporate it to the database but it's not what they asked: there will be a lot of different databases annotated by different people
      2. Blast and DADA2 alignment method are not identical, hence it's possible that the sequence was already inside of the reference and that every time we run it, we add a duplicate.
      3. Molly used a database and improved it but the link to the original database is dead so it's not even possible to know how it was made in the first place
   4. try to create a database from the DNA sequences by making a query request and then obtaining the right portion of the DNA sequence but soon realize that there is (once we get rid of all the un-annotated ones) 7576 sequences in the nucleotide between 100 and 350 base pairs (to obtain the partial CDS too).
      1. don't exactly remember why I stoped, maybe I should start again
   5. Use something called Identical protein groups and which is pretty cool and perfectly suited for what I want but there is few people using it and obtaining the Amino Acids sequence is easy by the DNA sequence is way harder because there is no package to do it directly so I ad to tinker:
      1. get all the ID of the IPG reports
      2. Download the IPG reports
      3. extract from these IPG reports the nucleotide accession number (often linking to he whole genome), the stand, the beginning of the sequence and the end of the sequence
      4. then download so custom GeneBank report containing the taxonomy, other info but only the part of the genome I'm interested in
      5. extract the taxonomy and the sequence for each file to put it into the reference file in the fasta format
      6. **realize that even though I limited the length of the sequences, I have some that are way longer than what they are supposed to and that for some report, inside of these reports, all the sequences don't have the same length even if they are supposed to be identical.**



**Maybe I should try to redo it with nucleotide to obtain a result**

but when I looked at the DNA sequences, there is some letters that are not DNA (R,Y,M,K, etc: IUPAC notation) and I'm not sure of how DADA2 will handle it

If DADA2 can use the IUPAC notation, then It can be really awesome but I kinda doubt of it and in the documentation, it look like this is only the case for the primers





























### Citations of articles:

Marine nitrogen (N 2) fixation is important in the global biogeo chemical cycling of N, and the interlinked cycling of carbon (C). Understanding and predicting global marine N 2 fixation requires information on diazotroph species specific rates of growth and N2 fixation, their biogeography, and the physiology of nutrient limitation in diazotrophs and nondiazotrophs (for the basis of competition in models). Because many diazotrophs are not easily visualized or quantified and many are yet uncultivated, the application of polymerase chain reaction (PCR) methods to amplify the nifH gene, which encodes nitrogenase, the enzyme that catalyzes N2 fixation, have led to multiple discoveries including new microorganisms (Zehr et al. 1998; Zehr and Capone 2020). Moreover, use of quantitative PCR (qPCR) has revealed unexpected biogeography of key marine diazotrophs regionally and globally (Church et al. 2005; Bentzon-Tilia et al. 2015; Langlois et al. 2015; Messer et al. 2015; Shiozaki et al. 2017; Harding et al. 2018; Mulholland et al. 2019). Biogeochemical modeling approaches have also been invaluable for providing hypothetical dynamic biogeography of the biomass associated with size classes of diazotrophs based on a number of assumptions (e.g., growth, nutrient uptake characteristics, and mortality; Dutkiewicz et al. 2014). It is, however, difficult to validate these models since there are few comprehensive datasets of N2-fixing organism biogeography.



Nitrogen gas (N2 ) fixing microorganisms (diazotrophs) are one of the most ecologically important functional guilds on Earth, providing the primary natural source for nitrogen to ecosystems through biological N2 fixation (BNF; Fowler et al., 2013). Isotopic evidence suggests that BNF has emerged as early as ca. 3.2 Gyr ago (Stüeken et al., 2015). It is thought to have evolved in an anaerobic archaeon and was later transferred to an aerobic bacterium (Boyd et al., 2015)



We investigated dinitrogen (­N2) fixation activity and diazotroph community composition across the Cape Verde Frontal Zone



The global nitrogen inventory is set by the balance between ﬁxed nitrogen gains in the form of biological dinitrogen (N2) ﬁxation and losses due to denitriﬁcation and annamox, but current measurements of both processes involve large uncertainties (Landolﬁ et al. 2018). Nitrogen availability fuels the biological carbon pump in otherwise nitrogen-limited systems such as subtropical gyres (Karl et al. 2012), and hence boost the ability of the oceans to cope with excess CO2 (Hutchins and Fu 2017). Thus, an accurate estimate of N2 ﬁxation is key to assessing the global ocean’s nitrogen inventory and its role in climate regulation. N2 ﬁxation represents a central reactive nitrogen supply in warm oligotrophic ocean regions, such as the subtropical gyres (Mahaffey et al. 2005), where dissolved inorganic nitrogen concentrations are typically close to the detection limit and surface temperatures exceed 20 C year round (Luo et al. 2012). Nutrient-rich coastal shelf seas and temperate estuaries were, however, long thought to lack signiﬁcant diazotrophic activity due to the prevalence of nitrogen-rich inputs from rivers, watershed runoff, and anthropogenic nutrient loading



Biological dinitrogen (N2) fixation, the reduction of atmospheric N2 to ammonia, is important for maintaining the fertility of the oceans by providing biologically useful nitrogen to support primary organic matter production (i.e., carbon dioxide fixation). N2 fixation offsets the removal of combined nitrogen by microbial denitrification and anaerobic ammonium oxidation (anammox) and export to the deep sea

Understanding the balance of the N cycle in the sea has wide-ranging implications for past, current, and future foodwebs, as well as for the role of marine N2 fixation in the sequestration of atmospheric CO2 and the production and consumption of other greenhouse gases such as nitrous oxide.

Heterotrophic marine diazotrophs are diverse but have not yet been definitively shown to fix N2 or to contribute substantially to water column N2 fixation.







# References

- Angel R, Nepel M, Panhölzl C,Schmidt H, Herbold CW, Eichorst SA and Woebken D (2018) Evaluation of Primers Targeting the Diazotroph Functional Gene and Development of NifMAP – A Bioinformatics Pipeline for Analyzing nifH Amplicon Data. Front. Microbiol. 9:703. doi: 10.3389/fmicb.2018.00703
- Hallstrøm, S., Benavides, M., Salamon, E.R. *et al.* Activity and distribution of diazotrophic communities across the Cape Verde Frontal Zone in the Northeast Atlantic Ocean. *Biogeochemistry* **160**, 49–67 (2022). https://doi.org/10.1007/s10533-022-00940-w
- Hallstrom S., Benavides Mar, Salamon E. R., Evans C. W., Potts L. J., Granger J., Tobias C. R.,  Moisander P. H., Riemann L.  (2022).     Pelagic N-2 fixation dominated  by sediment diazotrophic communities in a shallow temperate estuary.          *Limnology and Oceanography*,    67 (2),    364-378.      ISSN 0024-3590.
- Zehr JP, Capone DG. Changing perspectives in  marine nitrogen fixation. Science. 2020 May 15;368(6492):eaay9514. doi:10.1126/science.aay9514. PMID: 32409447.
- **Quantification of gene copy numbers is valuable in marine microbial ecology: A comment to Meiler et al. (2022)**
- https://www.jzehrlab.com/nifh