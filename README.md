# Pairwise Global Alignment

This project implements a pairwise global alignment using a linear gap penalty and an Affine gap penalty.

## Testing

I've provided a file called cullpdb_pc30_res3.0_R1.0_d191017_chains18877.gz that contains 1.5k proteins and their chains. Their protein sequences will be fetched from rcsb.org.

The project also takes another input which is a protein sequence, this sequence will be globally aligned with all the proteins in the file and the best alignment will be returned.

## Results

Given cullpdb_pc30_res3.0_R1.0_d191017_chains18877.gz and a protein sequence PAD with a basic scoring scheme and a linear gap penalty, an output tempelate should look something like this:

Best Global Alignment:

XDPXGGGGGNGDFEEIPEYL
  P        D        
--P-------AD--------

Score = -16
The time it took to perform the pairwise global sequence alignment was 548.3037974834442 seconds.