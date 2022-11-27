# biopipes
Bioinformatics pipelines

## CNV

Idea - get coverage on a per position basis and feed it into a NN.
* hg38 positions: 3.08829e+09
    - Keep an array for each chromosome's position. Storing the value as each datatype below allows for increased max coverage at the cost of increased space.
    - uint8:  Max value = 255, space = ~3 GB.
    - uint16: Max value = ~65K, space = ~6 GB.
    - uint32: MAx value = ~4.3B, space = ~12 GB.
