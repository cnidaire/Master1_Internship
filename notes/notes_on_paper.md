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