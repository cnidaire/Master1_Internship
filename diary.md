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



## Questions

- how was the database constructed (which criteria to select the sequences)?
- there is sequences duplicates: why? (there is 8 thousand duplicate sequences)
- sometime, same sequence but one is a bit longer, wouldn't it be better to keep only the longest one as it's assigning to the same taxonomy anyway?
- might seem a bit drastic but there is only 569 genus, can't we only have 569 sequences or do we keep it in order to maintain more variability in the sequences to have a more chance to attribute a sequence with DADA2?
- I didn't find the database you used as a starting point: from the Zehr Lab
- Good point is that everything is annotated to the class level at least.



# Pipeline

- [ ] Collect marine? nifH amplicon sequences from NCBI and ENA
- [ ] Run DADA2 with the latest nifH DB of Molly
- [ ] Generate fasta file, taxonomic and count table
- [ ] Upload into R package Ampvis
- [ ] Identify and extract unidentified sequences (e.g., sequence only identified at phylum level)
- [ ] Diamond blast and add hits to the database (remember to note down which was added).
- [ ] Manually, add the taxonomic annotation (i.e. kingdom, phylum, order and so on). Not sure
  have this can be done automatically. 