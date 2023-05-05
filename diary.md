# Day 1: 17/04/2023

Had a tour of the building and joined the different group chats, read bibliography and the is preparing for the meting at 4:30 where Elena presented the weekly paper. It was interesting.

I received the keys from Christan and I will have access to the computer.  I will have a meting with Christan and Molly on 28 so it would be nice if I knew a bit more about what they talk because they expect me to participate and I know nothing for the moment.



# Day 2: 18/04/2023

continue to work on understanding the Dada2 pipeline

It might be a good idea to pout the different parts in different files to makes modules, instead of a huge code.

I filtered the un-annotated sequences, It doesn't seem like there is any major mistake.

I tried to convert the csv file into a fasta file by making my own code because I didn't succeed to run Christian one and I think it's better to avoid for loops. But I two main problems, first I get some spaces between > and the number and the more sequences I add, the more spaces I get and the second is I have a space before the 2 and I wonder why and if I should remove it manually.



# Day 3: 19/04/2023

Continue debugging the csv to fasta code. And stopped using the basic R and used Tidyverse, it's faster (less or around a second) and is way cleaner to me.

I'm now trying to download the data from ENA and NCBI and to really understand Dada2

So for Dada2, am I supposed to 

**NTP** NCBI for accessing via terminal to NCBI



# Day 4: 20/04/2023

Working on Dada2 with the materials provided by Mar

And on watching some youtube videos and other to fix my knowledge in DNA sequencing and analysis of Omics data : https://www.youtube.com/watch?v=wV5_z7rR6yw

I think I'm should now be able to use DADA2



# Day 5: 21/04/2023

Trying to use DADA2

Setting up the other computer with the data, updating R, struggling to install DADA2

Then I switched to my computer but there is problems with the files:

- Dupe: 551_R1_001.fastq is in data_002 and there is a negative sample that is weird
- Swings: single end sequenced? There is only one file for each sample
- Tonga: a lot of non complementary sequences: where are the others?



To progress a bit, I added manually the 551_R1_001.fastq



# Day 6: 24/04/2023

I'm trying to install docker  to lab's computer but I'm struggling a lot, I wonder why nothing works as I want on this computer.



so there is forward and reverse primer for the TONGA too so I don't get why I miss so much sample files.

I tried to look at Subhadeep paper to find how it was sequences because to me it seems that it's not pair end sequencing but I'm not sure if it's he right one.

I'm trying to learn more about docker:

- to pull an image : docker pull image_name:version_number
- to go in a docker image: docker run -t -i image
- to quit this image: exit

Re run dada2 with Molly's reference file https://github.com/moyn413/nifHdada2

We obtain 80% of the sequences that are not identified at the genus level

So I have to gather more sequences to construct a bigger reference file contains 18 711 individuals

## group meeting

It was interesting, we talked about what we did and Mar presented a future project for 2025 and exposed the founding problems, etc



# Day 7: 25/04/2023

goals: either understand docker and start harvesting the sequences

## Python

```python
variable = Entrez.efetch(db='database',id='id_of_seq',rettype='fasta') # need Entrz and SeqIO in biopython
print(variable.read()) # shows at fasta format but not the right one for ref db
# would work if I have a database and an id for each sequences
# database can be 'nucleotides' so probably ''
variable = Entrez.efetch(db='nucleotide',id='id_of_seq',rettype='fasta', retmode='text')
record = SeqIO.read(variable, 'fasta')
record.id # gives the sequence id
record.name
record.seq # then I should be able to make a custom reference Fasta dataset

# Download multiple genomes
variable = Entrez.efetch(db='database',id='id_of_seq1,id_of_seq2,id_of_seq3',rettype='fasta')
record = SeqIO.read(variable, 'fasta')
record = [i for i in record] # to count the number of sequences then len(record)
record[1] # to access to the fasta file of the first sequence

# save in an output file (for fasta and multifasta)
SeqIO.write(record,'output_name','fasta') # writes a number correspoding to the number of sequences inside the file
```

otherwise, if I have the url, I can download it with wget in shell

## Docker

[cheat sheet docker](https://learninglab.gitlabpages.inria.fr/mooc-rr/mooc-rr2-ressources/module2/seq2-package_mgmt/unit2-references.html)

```shell
docker run -it image_name:version # run interactive
docker run -it -d image_name:version # run in background
docker stop id # stop the docker with this id in the background
docker ps # gives the active dockers
docker system prune # will clean up everything
docker run -p image_name:version # can define the "ports" but idk what it does
```



### Docker file

```shell
FROM debian:9 # FROM is only usable once

RUN apt-get update -yq \ # to execute something in the container
&& apt-get install curl gnupg -yq \
&& curl -sL https://deb.nodesource.com/setup_10.x | bash \
&& apt-get install nodejs -yq \
&& apt-get clean -y

ADD . /app/ # to copy or download files in the image
WORKDIR /app # change the working directory, equivalent of cd
RUN npm install

EXPOSE 2368
VOLUME /app/logs

CMD npm run start
```



install.package("name") to install a package in R



# Day 8: 26/04/2023

I extracted the sequences of the dupe database that weren't annotated up to the genus level and then I tried to use blastn on the NCBI websit because I'm not able to do the blast on my computer with all the references datasets of ncbi. I specified that I exclude exclude Uncultured/environmental sample sequences and Models(XM/XP) to avoid as much un-annotated data as possible. Now I have a CSV file withe the ID, E-value, %query, etc.

Now I have to choose which of the results are plausible.



I plan to use blastn to align he sequences. Depending on the output, I want to gather the first sequences, check the E-value, how much is aligned, and for the first ones if it has the same phylogeny and if It is something different than transposons or uncultured bacterium. Be careful to put a threshold to the E-values like I think at least e-90/e-100

exclude Uncultured/environmental sample sequences and Models(XM/XP)



If the first sequences are as likely to be the real ones and have different phylogeny, then we keep only the parts that are identical.

By trying on a few sequences manually, I checked that I works really well for some sequences



I think I finally understood why there is so much non identified sequences, it's because the primers didn't fix at the right position and amplified some other sequences that are not nifh related. But It's still possible to identify it with blastn because we align on the whole genome.



```shell
blastn -query input_file_name -db reference_db -evalue 10 -out output_file_name # with outfmt, youcan choose the format of the ouput and then you can form exemple obtain the sequences id to then put it in the reference file.
```

I will run blastn manually on the NCBI website and then extract the data manually. I think the csv format even if I hate csv migth be the best because it's easily be opened in R



# Day 9: 27/04/2023

try to understand the accession number to know if it's a gene in order to filter out the full genome, etc.



I tried to specify in "Entrez Query" nifH[TW] so that it's only aligned with the annotated sequences containing nifh in the text (title and abstract, and MeSH terms, subheadings, chemical substance names, personal name as subject, and MEDLINE Secondary Source (SI) field).

nifH[TW] NOT unidentified[TI]

do I have to add CDS too in order to only have the part of the genome I'm interested in?

All the taxonomic names ending in ota for phylum names (like Cyanobacteriota) are the new taxonomy. Hopefuly, the phyluum names have been updated in NCBI on January 2023



for exemple, https://www.ncbi.nlm.nih.gov/nuccore/AY786992.1 contain the nifH but not only. Then it's annoying because if I put it in the nifH reference database, DADA2 can align it with nifH but not only.



**In the sequences I have selected, there is some with an % identity that is around 72% it's of course way too low but It would be a good idea to include it because the sequence is still relevant (It's a sequence between 100 and 1000 bp, containing nifH and that is not unidentified)**



## Goal

Find all of these sequences that are relevant and aligned with the non identified sequences with a blast

- put the un-identified sequences in a FASTA file
- use blast n on these seq (nifH[TW] AND 100:1000[SLEN])
- take all the unique sequences id (I wonder if I only take the most relevant one as I'm constructing a db and not assigning the taxonomy)
- get the sequences and the taxonomy (from GenBank) and make a FASTA file out of it 

Then fuse the two databases by omitting the sequences that were already inside of the reference (if there is one).

- check if the sequence of the new file is already in the old file
- if not, add in the new file

Re-run DADA2 on the sequences with the new reference DataBase

- the pipeline is already ready



# Day 10: 28/04/2023

The references sequences I obtained yesterday seems pretty good, at least, it is sequences that are all relevant.

Maybe I can do what I did before with biopython to even automatize that. It seems more adapted to what I want than Biomartr

It would be nice in the next database to put the accession number to have a reference from where this sequence is coming from. It will help for example if the sequence is updated, so that we don't keep the old version in the reference.





I looked at what are the questions/ problems I have with he nifH reference database.

## Problems of the database

there is a lot of sequences to whom there is only the up to the kingdom or Phylum level (mostly cyanobacteria, proteobacteria or archea)

- there is a lot of identical sequences for the same taxonomy... where was it found? Wouldn't it be better if there was only one verified sequence instead of a lot with a few bp difference? (example: line 5054 to 5058 where the sequences are exactly the same)
- The accession number correspond to a protein but the sequences stored in the database are nucleotides
- I'm sure that the duplicates are completely useless, but I don't know for the sequences that differ from a few databases and that are associated to the same phylogeny... as this is a variable sequence, it should be significantly different from other species (there is only 875 species with the same taxonomy)
- there is three non unique accession numbers



add taxonomy as an output by incorporating it in the taxonomic level

## Questions

- how was the database constructed (which criteria to select the sequences)?
  - 
- there is sequences duplicates: why? (there is 8 thousand duplicate sequences)
  - 
- sometime, same sequence but one is a bit longer, wouldn't it be better to keep only the longest one as it's assigning to the same taxonomy anyway?
  - 
- might seem a bit drastic but there is only 569 genus, can't we only have 569 sequences or do we keep it in order to maintain more variability in the sequences to have a more chance to attribute a sequence with DADA2?
  - 
- I didn't find the database you used as a starting point: from the Zehr Lab
  - 
- Good point is that everything is annotated to the class level at least.
  - 



# Day 11: 02/05/2023

I think building a new database from scratch is a better option because there is a lot of redundant sequences and I feel like it isn't optimal to have such a non controlled database (for example, some sequences are really old and my have a different annotation now but there is no control over it).

So I will download the protein package on blast then make a query search to gather all the data with the right filter corresponding to nifH and make a FASTA file out of it.

After, I will have to check all the identical sequences to check if the same classification is assigned, if not I make a new taxonomy with the common part (among the most annotated because there can be sequences with non complete annotation and in this case, we would loose information). And if it is, then I remove arbitrarily all the useless sequences to keep only one sequence.

Then, if I have done everything, I retry with DADA2 and if I obtain to see if it works.

Then I can add the sequences from DNA too in an other file and only add the ones hat are not already inside. 

## Attempt

There is a database called Identical protein group that is giving me all the unique sequences coming from different databases sources: refseq, INSDC, PDB, UniProtKB/Swiss-Prot, PIR, PRF. It is really amazing!!! I obtain 6783 sequences.

Then with the file, I obtain the acession number but it means that I will have to deal with different web interfaces, that's the main down side. But most of the sequences are coming from INSDC and RefSeq (except 59 sequences)



My problem so far is that I obtain the whole genome, it is really slow and I don't even have the taxonomy



# Day 12: 03/05/2023

If I can obtain the  CDS Region in Nucleotide from identical protein group, then I can obtain directly the sequence if I find where I have to search.

But It seems that NCBI didn't implement any search tool to obtain the right sequence directly from an ID like "NZ_FQZT01000024.1 24762-25643 (-)".

I will have to obtain the sequence of the whole genome and then only keep the parts I'm interested in. furthermore, I have to be carfull if it's the 6 or + strand. I think using identical protein group is the good option. The easiest way might be to retrieve the whole sequence and then to make a fasta file out of the part I'm interested in. The problem is that take the coding sequences of the whole genome would take around 4Mo and I have a bit less than 10k sequences, I would make downloading 40 go instead of 1Ko per sequences making 10 Mo to download which is reasonable. I can't use the first option, I would put NCBI servers down, restrict the other users from doing important stuff or just being stopped in the middle.

I think it would be better to retrieve the id directly in NCBI because it's not possible to use identical protein groups eiter in python or in R but after, maybe using python might be more suitable to extract and filter the data.



For the moment, I have to search in nucleotides with the Nucleotide.accession but it gives me the whole genome so now I'll have to see if possible to only obtain the pat I'm interested in.



## Bio-python

**to recap, I want two thing, giving a Nucleotid.accession, a start + a stop(if I only want to take the whole genome and only the CDS) of just the protein name (if I only take the CDS, then I have to select the correct sub-unit), the taxonomy and the sequence of the gene of interest**



```python
from Bio.Seq import Seq
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
my_seq[4:12] # if only want bp from 4 t 12, gives Seq('GATGGGCC')

```

### SeqRecord object

offers the following information as attributes:

- **.seq – The sequence itself, typically a Seq object.**
- **.id – The primary ID used to identify the sequence – a string. In most cases this is something like an accession number.**
- **.name – A “common” name/id for the sequence – a string. In some cases this will be the same as the accession number, but it could also be a clone name. I think of this as being analogous to the LOCUS id in a GenBank record.**
- .description – A human readable description or expressive name for the sequence – a string.
- .letter annotations – Holds per-letter-annotations using a (restricted) dictionary of additional information about the letters in the sequence. The keys are the name of the information, and the information is contained in the value as a Python sequence (i.e. a list, tuple or string) with the same length as the sequence itself. This is often used for quality scores (e.g. Section 20.1.6) or secondary structure information (e.g. from Stockholm/PFAM alignment files).
- .annotations – A dictionary of additional information about the sequence. The keys are the name of the information, and the information is contained in the value. This allows the addition of more “unstructured” information to the sequence.
- .features – A list of SeqFeature objects with more structured information about the features on a sequence (e.g. position of genes on a genome, or domains on a protein sequence). The structure of sequence features is described below in Section 4.3.
- .dbxrefs - A list of database cross-references as strings.

```python
from Bio import SeqIO
record = SeqIO.read("NC_005816.gb", "genbank") # if we have the database localy on the computer
record
```

### SeqFeature objects

the Biopython SeqFeature class attempts to encapsulate as much of the information about the sequence as possible: describe a region on a parent sequence, typically a SeqRecord object:

- **.type (CDS or gene)**
- **.location**
- .ref
- .ref_db
- **.strand**

### Slicing a SeqRecord (page 42)

```python
from Bio import SeqIO
record = SeqIO.read("NC_005816.gb", "genbank")
record
len(record)
len(record.features)
print(record.features[20]) # access to the feature 21
print(record.features[21])

```



# Day13 : 05/05/2023

Continue the tutorial of biopython.

I think I will need the genbank format because it contains the taxonomy and the sequence and then make a Fasta file out of it.

Apparently, It might be better to first save the files on my computer and then work on it.

I just realized that the accession number is usually giving me the whole genome and in  Amino Acids but contain the taxonomy. And on the other hand, the id contains only the sequences (DNA) I'm interested in but the taxonomy is  unidentified. So I want the sequences coming from the id and the taxonomy of the Accession number.



The thing is that I have a lot of similar sequences (CDS coming form the same genome). Hence I have to find a satisfying file format containing the whole genome (because when I only have the CDS, It's the CDS for the + strand and hence I can't obtain the complementfor the minus strand: with full sequences, I can obtain the complement with Biopython). **THE FASTA FORMAT, OF COURSE!!!** Then I have a way to proceed:

-  collect the full fasta files for the unique sequences: only around 2K for NCBI and store all of these files individually locally but I don't know yet how to obtain it from the nucleotide accession number
- After It, find a way to

For the moment, I'm against the wall because I wanna use the Identical Protein Group but there is a huge lack of information on it making that I don't really know what I'm doing.

https://astrobiomike.github.io/unix/ncbi_eutils might save me





### Parsing sequences from the net (page 54)

```python
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "A.N.Other@example.com"
with Entrez.efetch(
	db="nucleotide", rettype="fasta", retmode="text", id="6273291") as handle:
	seq_record = SeqIO.read(handle, "fasta")
print("%s with %i features" % (seq_record.id, len(seq_record.features)))

# for multiple records

from Bio import Entrez
from Bio import SeqIO
Entrez.email = "remi.legrand38@gmail.com.com"
with Entrez.efetch(
	db="nucleotide", rettype="gb", retmode="text", id="6273291,6273290,6273289") as handle:
	for seq_record in SeqIO.parse(handle, "gb"):
		print("%s %s..." % (seq_record.id, seq_record.description[:50]))
		print(
			"Sequence length %i, %i features, from: %s"
			% (
				len(seq_record),
				len(seq_record.features),
				seq_record.annotations["source"],
			)
		)
```

```python
# for SwissProt sequences

from Bio import ExPASy
from Bio import SeqIO
with ExPASy.get_sprot_raw("O23729") as handle:
seq_record = SeqIO.read(handle, "swiss")
print(seq_record.id)
print(seq_record.name)
print(seq_record.description)
print(repr(seq_record.seq))
print("Length %i" % len(seq_record))
print(seq_record.annotations["keywords"])
```



## Interest of Identical Protein Group

In 2014 NCBI introduced the ‘Identical Protein Report’ to the Protein  database to clarify the relationships between WP sequences and the set  of individual Nucleotide CDS sequences they represent ([8](javascript:;)). Now these reports have been improved and collected in a new resource called Identical Protein Groups ([www.ncbi.nlm.nih.gov/ipg/](http://www.ncbi.nlm.nih.gov/ipg/)). This IPG resource includes all NCBI protein sequences, including  records from INSDC, RefSeq, Swiss-Prot, and PDB, with links to  nucleotide coding sequences from GenBank and RefSeq. The title of each  record is derived from the ‘best’ sequence in each group, where the  hierarchy for determining the best sequence is RefSeq > Swiss-Prot  > PIR, PDB > GenBank > patent ([www.ncbi.nlm.nih.gov/ipg/docs/faq/](http://www.ncbi.nlm.nih.gov/ipg/docs/faq/)). Searches in this database can be filtered by database source, taxonomy, and the number of sequences in the group. These reports continue to be  available through the E-utility EFetch with *&db = protein&rettype = ipg* ([eutils.ncbi.nlm.nih.gov](http://eutils.ncbi.nlm.nih.gov)).



# Day 14: 05/05/2023

I gave up on using IPG, I want to obtain some exploitable results first.



## Gather ids of sequence of interest in the protein database

```python
from Bio import Entrez
# function to search data on NCBI
def search_ipg(query):
    handle = Entrez.esearch(db="protein", term=query, retmax=20000) # number of sequences
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]
# create a list containing the ID
query = "((nifH[Gene Name]) AND 100:350[Sequence Length]) NOT uncultured"
id_list = search_ipg(query)
print(result)
```

## Obtain the CDS (DNA) for a single sequence

around 1Ko

```python
from Bio import SeqIO
Entrez.email = "remi.legrand38@gmail.com"
with Entrez.efetch(
	db="protein", rettype="fasta_cds_na", retmode="text", id="2497592678") as handle:
	seq_record = SeqIO.read(handle, "fasta")
print("%s with %i features" % (seq_record.id, len(seq_record.features)))

seq_record
# in a function to make it cleaner and call it in a loop after:
def fetch_cds_fasta(protein_id):
    Entrez.email = "remi.legrand38@gmail.com"
    handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta_cds_na", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return record

# Usage
protein_id = "2497592678"  # Replace with the protein ID you want to retrieve
result = fetch_cds_fasta(protein_id)
print(result)
result.seq # to obtain only the sequence
```

# Obtain the taxonomy for a single sequence

around 2Ko for a genpept file

```python
from Bio import Entrez
from Bio import SeqIO

def fetch_cds_fasta(protein_id):
    Entrez.email = "remi.legrand38@gmail.com"
    handle = Entrez.efetch(db="protein", id=protein_id, rettype="gp", retmode="text")
    # genpept_data = handle.read()
    record = SeqIO.read(handle, "gb")
    handle.close()
    return record

# Usage
protein_id = "2497592678"  # Replace with the protein ID you want to retrieve
result = fetch_cds_fasta(protein_id)
# print(result)
result.annotations["taxonomy"]
```



# Pipeline

- [ ] Collect marine? nifH amplicon sequences from NCBI and ENA
- [ ] Run DADA2 with the latest nifH DB of Molly
- [ ] Generate fasta file, taxonomic and count table
- [ ] Upload into R package Ampvis
- [ ] Identify and extract unidentified sequences (e.g., sequence only identified at phylum level)
- [ ] Diamond blast and add hits to the database (remember to note down which was added).
- [ ] Manually, add the taxonomic annotation (i.e. kingdom, phylum, order and so on). Not sure
  have this can be done automatically. 