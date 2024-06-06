# TileCreator

TileCreator is a python script generated to automatically design primers for overlapping tiles sweeping a given transcript. It forces the location of the first forward primer at the beginning of the transcript and the last reverse primer at the end of it, extending them until achieving the desired Tm. Then, it slices the transcript according to the desired tile size, and adds reverse primers 35 nt downstream of that site and forward primers 35 upstream, also extending them until achieving the desired Tm. 

It does not check for primer specificity nor primer dimers. 

It is based on the original olaygo.py script found in seismic-rna by Matthew F. Allan (please check the Rouskin Lab repository). 

0) Before you start: you will need to have installed python and biopython. 

1) Input data: a fasta file with the nucleotide sequence of interest (does not recognize "U", so please use the DNA sequence with "T").

2) Output data: the script will output a fasta file with all the pairs of primers. 

3) Usage:

`python TileCreator.py input.fasta output.fasta tile_size tm`

Where tile_size is the desired size for the tiles, so a number of nucleotides; and tm is the desired melting temperature (number). 

Example: 

`python TileCreator.py mRNA.fa primers.fa 400 60`

4) Posterior analysis: seismic-rna can be used to stitch together the tiles. 

# TileCreatorApp

TileCreatorApp is just the addition of a user interface to TileCreator. 

0) Before you start: you will need to have installed python, biopython, ttkbootstrap and tkinter. 

1) Usage:

`python TileCreatorApp.py`

Then just work your way through the user interface. 

# Standalone TileCreatorApp

Alternatively, an standalone app can be used that requires no interaction with the command line. Only available for MacOS for the moment. Please bear in mind that it takes a long time to open, but then is works quickly. 

Download from here: [https://sourceforge.net/projects/tidecreator/files/latest/download](https://sourceforge.net/projects/tilecreator/)


