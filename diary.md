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

So I have to gather more sequences to construct a bigger reference file

# Pipeline

- [ ] Collect marine? nifH amplicon sequences from NCBI and ENA
- [ ] Run DADA2 with the latest nifH DB of Molly
- [ ] Generate fasta file, taxonomic and count table
- [ ] Upload into R package Ampvis
- [ ] Identify and extract unidentified sequences (e.g., sequence only identified at phylum level)
- [ ] Diamond blast and add hits to the database (remember to note down which was added).
- [ ] Manually, add the taxonomic annotation (i.e. kingdom, phylum, order and so on). Not sure
  have this can be done automatically. 