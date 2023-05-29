---
title: "The Document Title"
author: [Example Author, Another Author]
date: "2017-02-20"
keywords: [Markdown, Example]
abstract: |
	Short version of the introduction: to do later
titlepage: true
---

<!-- Voir -->
<!-- https://github.com/Wandmalfarbe/pandoc-latex-template/tree/v2.1.0 -->
<!-- pour les options de mise en page-->

# Introduction

## Why study dinitrogen fixating bacterias?

Many living organisms are unable to metabolize directly the atmospheric nitrogen. Thus, the main source of nitrogen in the ocean is coming from Dinitrogen (N2) fixating Bacteria (Cyanobacteria and non-Cyanobacteria) that are transforming it into other forms that can be used by the other organisms.

This fixation has a role in the global biogeochemical cycling of Nitrogen, and is interlinked withe the cycling of Carbon. Therefore, understanding the Nitrogen cycle would help to understand the current, past, and future food webs. As well as the storage dynamics of atmospheric CO2, and other greenhouse gases like nitrous oxide (that is a greenhouse gas 300 times more effective than carbon dioxide).

## Why nifH

Nitrogen fixation is catalyzed by the enzyme Nitrogenase. The *nif* genes are a collection of genes involved in nitrogen fixation. The Nitrogenase protein is composed of two sub-units metalloproteins: molybdenum iron protein and iron protein. It is likely to be an ancient protein since it's distributed widely through Bacteria and Archaea. Both of the Nitrogenase proteins are highly conserved but the Fe protein encoded by nifH is the most highly conserved.

Thus NifH is the ideal gene to study nitrogen-fixating bacteria as it is highly conserved but not identical between the Bacteria's species and hence for environmental samples, we can identify the specie the protein belong to.



## How to collect environmental data: Metabarcoding

Metabarcoding is a tool to characterize the taxonomic diversity of an ecosystem using environmental DNA. Usually, the region of the gene we focus on for the identification of the Bacterias and Archaeas is the 16S region due to various properties such as the size (~200-400bd) and the good conservation over the different species. However, it is possible to use some other parts of the genome if it is suitable for identification, like the nifH gene.

The metabarcoding itself is divided into a few steps:

First of all, we sample the middle we are working with, in this case, it's ocean or seawater. Then we extract the DNA and amplify using Polymerization Chain Reaction (PCR) of the region we are interested in using the appropriate primers. The next step is to sequence the data using high-throughput sequencing (like Illumina sequencing). We usually use pair hand sequencing and pool all the samples together to avoid as much batch effect as possible. Finally, the result is processed and demultiplexed through an adequate pipeline to obtain the annotation of the sequences.



## Goals of the internship

To annotate the sequences coming from the nifH environmental data, the DADA2 pipeline is generally used and a reference file is needed. This reference file contains the sequences and the taxonomy associated. The current database was first created in 2014 and kept updated until 2017 by the Zehr Lab (Heller 2014). It was then normalized and kept updated manually by Molly Moynihan (put publication).

However, this database annotates quite poorly the data obtained. Resulting in a lot of uncertainty and the goal would be to either improve the current database or build a new one and find a way to automatize it.

There are many problems I noticed with the current database:

1. It has been created manually and hence there is a possibility of typos
2. It is not updated automatically: if some sequence's taxonomy attribution changes it will not be updated
3. Recently the taxonomy changed and they had to modify the 18,700 sequences manually, resulting in a lot of work and possible mistakes
4. Lack of control over what is inside the database: there are a lot of identical sequences and hence we hope that the sequences have the same taxonomy associated otherwise (normally, there is an algorithm to check it and automatically annotate it in NCBI implemented in 2001 but as it's created manually, there is a possibility of typos once again), randomly chosen by DADA2
5. It is hard to know which sequences are not in the reference yet

For all these reasons, I think building a new fully automatized database from scratch would be more reasonable than improving and building on top of something we are not fully confident in.



# Methods

**Put the code in the annexes**

## Dada2 Workflow

For the our study, we are interested in the workflow until "Assign Taxonomy".

![](notes/dada2/images/workflow.png)

### Raw FASTQ Files

We first have to import the Raw FASTQ files. In most of the cases, we have pair hand sequencing so we have to be careful to associate the reads pairwise.

### Trim and Filter Reads

Then we plot the quality score and determine the trimming cut-offs (traditionally we take 30 as the lowest accepted quality score). Usually, the quality of the reverse sequence is not as good as the forward, that's why we have to trim it more. However, we have to make sure that the Forward and reverse sequences are overlapping because DADA2 need at the very least 20 base pair to successfully merge the two sequences.

This steps helps to remove the sequencing errors and artifacts, in order to keep only the high quality sequences.

**Annotate + explain figure**

![Workflow](img/quality_score.jpg)

### Error Rate estimation

DADA2 constructs an error model by learning the sequencing errors from  the quality-filtered reads. The error model is then used to correct  errors in the reads, improving the sequence accuracy.

In general, for the error plot, the frequency of the errors rate decrease as the quality score increase.

**Annotate figure if kept**

![Workflow](img/error_rate.jpg)

### Dereplicate Reads

Condense the data by collapsing together all reads encoding for the same sequence to reduce the redundancy. It will speed up and simplify the computation.



### Sample Inference with DADA2 Algorithm

Testing the null-hypothesis that the sequence is too abundant in the sample to be solely explained by errors in the data set. Hence, a low p-value sequence can be considered as real sequences that are not caused by random errors and in contrary, if the sequence has a high p-value, it won't be kept 

**Now all the Forward and Reverse reads have been denoised** we can finally merge all the forward and reverse sequences.


​	

### Merge Reads

We have inferred the sample sequences in the forward and reverse reads independently. Now it’s time to merge those inferred sequences together, throwing out those pairs of reads that don’t match. It will return a data frame corresponding to each successfully merged sequences.



### Chimera Checking and Removal

A chimera is a fusion of two or more parents sequences and can be created during the PCR amplification process. In order to spot the chimera sequences, we perform multiple sequence alignment from the least abundant read and for all the more abundant read, it will do sequence alignment with all possible combinations. When a chimera is detected, it is removed from the sequence table.

More than 90% of the unique sequences identified are bimeras and most of the total reads shouldn't be identified as chimera.

If there is too much chimera, we have to check if the 20 first base pairs of the reads were really trimmed because contain primers that can artificially increase the number of chimera. But if it is not the case, we should try to trim more of the low quality bp.



### Assign Taxonomy

This is the part I'm interested in. We should give as an input a reference database at the FASTA format containing as the first line the Taxonomy of the sequence at the format: Domain / Pylum / Class / Order / Family / Genus; and as a second line the corresponding sequence of the nifH gene.

Then DADA2 will infer the taxonomy based on the sequence similarity.



## NCBI databases

NCBI: National Center of Biotechnology Information, is an American institute developing software to analyze the genome data. NCBI is not a database but it develops many databases such as GenBank, RefSeq, and PubMed which are the most famous ones. Here I'm only going to present the main ones that might be interesting for the study.

The **gene** database integrates information from a wide range of species. A record may include nomenclature, Reference Sequences (RefSeqs), maps, pathways, variations, phenotypes, and links to genome-, phenotype-, and locus-specific resources worldwide.

The **proteins** database is a collection of sequences from several sources, including translations from annotated coding regions in GenBank, RefSeq, and TPA, as well as records from SwissProt, PIR, PRF, and PDB. Protein sequences are the fundamental determinants of biological structure and function.

It is really interesting in our case because the nifH gene is coding for a protein and hence we can access the taxonomy and many other interesting information on different databases. The downside is that we obtain the Amino Acid sequences and not the DNA ones and we have duplicates.

The **Identical Protein Groups (IPG)** database takes in my opinion the best of both worlds: we obtain all the sequences clustered when these are identical, we obtain the DNA sequences and it's gathering the information from different databases such as Swissprot, GenBank, PIR, etc.

The downside is that Entrez is not really suitable for what I want because you have to pass by gene db to obtain the DNA sequence and there is not downloadable directly from the IPG db on NCBI so I have to send a request to consult the gene database.



## Entrez to access NCBI library from python

Entrez, Global Query Cross-Database Search System, is a searching tool provided by NCBI enabling one to browse through the information of over 20 databases such as Swiss-Prot, PRF, and PIR-International. It is an indexing and retrieving system gathering data from various sources and handing it back in a uniform format such as FlatFile, Fasta, or XML. Thanks to it, it is possible to retrieve information without going through the NCBI website and it makes it easier to automatize it when working on large datasets like the one we want to make.



## Implementation with python

I want to make a program that is gathering the sequences and the taxonomy from NCBI and maybe later from other databases for all the record available of the DNA sequences of the nifH gene. Here, I didn't worry about the computational time at all because this program is not ought be re-run on a daily basis but more every year or so if it work to create a new version of the database that is uniform and used by all the community.

I first choose the gene proteins but there were not a lot of sequences and then when I looked at the nucleotide database, I realized that it is at the IUPAC format. The IUPAC DNA format is often used because instead of only containing A, T, C, G and U, it contains R,Y, S, W, etc. It gives the information when two or more nucleotides can be used. However, DADA2 is not able to use well such a FASTA reference file.

In the end, I think the best idea is to use the Identical Protein Report. It will give us individual reports for each identical Amino Acids sequences. Even if it doesn't give access right away to the DNA sequence, inside of each report, there is the database it is present in, the reference of the DNA file, the Start and Stop because most of the time this DNA sequence is part of a bigger DNA sequence and the strand. With this information, we can download the corresponding DNA sequence for each references in the report. 

However, I first thought that only one sequence of each report was enough because I thought that the DNA sequences were identical, I mistake that I would have been able to spot on my own if I took the time to take a step back. When I obtain my reference file, I obtained some inconsistencies inside. That's why I contacted NCBI and only recently received a reply that answered most of my questions. but I haven't been able to implement it yet.



Therefore, here is how I proceeded with my python code so far.



- Gather all the ID of the IPG reports from a Query request

- Download the corresponding IPG report

- Extract the first line of the IPG report present in NCBI and put it in the table
  - normally all the sequences are identical and the taxonomy should be the same because NCBI check it and automatically annotate so It doesn't matter which sequence I select. As I'm currently dealing mostly with NCBI, it's easier to retrieve the sequences from NCBI (there is only 2 IPG reports that don't have any sequence in NCBI so I either I will have to do it manually as it's only two or find an other way because it's an other tool than Entrez)
- Fetch with entrez the XML report for each sequence in the previous table (**what is rettype='gb' then?**)
  - There is different formats that we can obtain but the most practical one is XML because it contains the taxonomy and the sequence in the same file. Moreover, it has a structure making it really easy to access to specific data inside of the file.
  - As most of the DNA sequences the IPG report is redirecting me too are genome or group of genes, I have to ask specifically the DNA sequence starting at a given point and stopping at an other (the information is in the IPG report). And we have to be careful of the DNA strand.
- Make a FASTA file by extracting for each file the taxa and the sequence.
  - the fasta file start with a ">" followed by a description of the of the sequence on the first line (in our case the taxonomy) and on the second line the sequence itself. And so on for all the sequences in the file
  - Sadly, there is some sequences or taxonomy that are not found and I didn't have the time to find a way to fix it.



# Results

As a point of comparison, I use the data coming from DUPE exploration. 90% of the DUPE data has a kingdom associated after the annotation, and in this case we only have Bacteria which is given because we are working on Bacteria. However, at the phylum level we already have 57% that are not annotated and this is increasing until the family level where we have nearly 80% of the sequences that are unannotated. And this is a major problem when we want to characterize a population.

I haven't been able to finish the database. In the current database I created, there is some inconsistencies and some of the sequences I fetch are too big so I asked questions to NCBI and now I think I know how to continue to improve the database and solve my problems.

For each DNA sequence associated with each IPG report, I have to split in different panda dataframe containing the id and location of the DNA file. As having the same protein sequence doesn't mean having the same DNA sequence, I will have to take all the report. I will still use IPG even if the reason I used it in the beginning doesn't stand anymore, because it enable me to obtain information for all different main databases in one go. Then I will have to download the sequences for the other databases but it will mean to learn an other too.



| Taxonomy | With old data base (% of unannotated) | With new database (% of unannotated) |
| -------- | ------------------------------------- | ------------------------------------ |
| Kingdom  | 10.39 %                               |                                      |
| Phylum   | 57.07 %                               |                                      |
| Class    | 65.78 %                               |                                      |
| Order    | 75.94 %                               |                                      |
| Family   | 78.59 %                               |                                      |
| Genus    | 80.17 %                               |                                      |



# Conclusion and Discussion

In conclusion, I am not so far from being able to create the database containing the sequences from NCBI. However, to have a complete database, I still need to add the DNA sequences of the nifH gene coming from the other databases such as UKProt and SwissProt. It will take time because the syntax and the packages to automatically retrieve the information will be different and will need more time than the 7 week I was given.



<!---

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



# What I still have to do



**Maybe I should try to redo it with nucleotide to obtain a result**

but when I looked at the DNA sequences, there is some letters that are not DNA (R,Y,M,K, etc: IUPAC notation) and I'm not sure of how DADA2 will handle it

If DADA2 can use the IUPAC notation, then It can be really awesome but I kinda doubt of it and in the documentation, it look like this is only the case for the primers



From the answer I received from NCBi, I think all the sequences are important and It would be easier to just take all the NCBI sequences in the IPG report and make other tables for each databases and then as a project for after the internship, retrieve fetch the information for the other database.

However, I have to be careful to limit the length of the DNA sequences I obtain because it still doesn't explain why there is some sequences that a way bigger than some others.

-->



# References

- Angel R, Nepel M, Panhölzl C,Schmidt H, Herbold CW, Eichorst SA and Woebken D (2018) Evaluation of Primers Targeting the Diazotroph Functional Gene and Development of NifMAP – A Bioinformatics Pipeline for Analyzing nifH Amplicon Data. Front. Microbiol. 9:703. doi: 10.3389/fmicb.2018.00703
- Hallstrøm, S., Benavides, M., Salamon, E.R. *et al.* Activity and distribution of diazotrophic communities across the Cape Verde Frontal Zone in the Northeast Atlantic Ocean. *Biogeochemistry* **160**, 49–67 (2022). https://doi.org/10.1007/s10533-022-00940-w
- Hallstrom S., Benavides Mar, Salamon E. R., Evans C. W., Potts L. J., Granger J., Tobias C. R.,  Moisander P. H., Riemann L.  (2022).     Pelagic N-2 fixation dominated  by sediment diazotrophic communities in a shallow temperate estuary.          *Limnology and Oceanography*,    67 (2),    364-378.      ISSN 0024-3590.
- Zehr JP, Capone DG. Changing perspectives in  marine nitrogen fixation. Science. 2020 May 15;368(6492):eaay9514. doi:10.1126/science.aay9514. PMID: 32409447.
- **Quantification of gene copy numbers is valuable in marine microbial ecology: A comment to Meiler et al. (2022)**
- https://www.jzehrlab.com/nifh
- M. A. Moynihan. 2020. nifHdada2 GitHub repository. Zenodo. http://doi.org/10.5281/zenodo.3958370



# appendix

Put the Python code

```{python}
print(2+2)
```

