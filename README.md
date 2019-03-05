# motif-mark
Tool for visualizing motifs in nucleotide sequences.

## Inputs
 * A FASTA file of nucleotide sequences. If lowercase is used for intron and uppercase for exons, this will be reflected in the resulting SVG.
 * A text file of motifs to look for, 1 per line with no other delimiter.
 * (optional) A text file of colors to use for motifs in hex and/or RGB (255) format, 1 per line with no other delimiter.

## Output
 * A single SVG (scalable vector graphic) image with motif locations drawn on diagrams of the sequences.

## Parameters
  -f/--fasta  path to FASTA file of sequences, optionally with introns lowercase and exons uppercase
  
  -m/--motif  path to file of motifs, 1 per line
  
  -c/--colors (optional) path to file of colors in hex and/or RGB (255) format, 1 per line
  
  -o/--output (optional) name of output file, default is motifs.svg
  
  -a/--allowOverlaps (optional) whether to allow overlapping copie sof a single motif, default is True
  
## Requirements
 * Python3
 * PyCairo
 
