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

So for Dada2, am I supposed to 

**NTP** NCBI





# Pipeline

- [ ] Collect marine? nifH amplicon sequences from NCBI and ENA
- [ ] Run DADA2 with the latest nifH DB of Molly
- [ ] Generate fasta file, taxonomic and count table
- [ ] Upload into R package Ampvis
- [ ] Identify and extract unidentified sequences (e.g., sequence only identified at phylum level)
- [ ] Diamond blast and add hits to the database (remember to note down which was added).
- [ ] Manually, add the taxonomic annotation (i.e. kingdom, phylum, order and so on). Not sure
  have this can be done automatically.