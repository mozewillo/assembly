# OLC reads assembler
Assembling genetics sequences(reads) into contigs based on Overlap Layout Consensus approach.


**Usage**   
assembly <input_reads> <output_contigs> [optional arguments]

```
Help
assembly -h

positional arguments:
  input_reads                   input reads to assemble (fasta format)
  output_contigs                output file name (fasta format)
  
optional arguments:
  -cov, --coverage              sequencing coverage
  -minlap, --minmial_overlap    minimal overlap for the sequences
```
