# Indel Inference using dynamic programming

Implementation of indel inference based on Phylogenetic Profiling of Insertions and Deletions in Vertebrate Genomes [[1]](#1).

### Running the method
```
python spindel.py -n phylogenetic tree -a alignment file -o output folder

arguments:
  -a  fasta format alignment file
  -n  ancestors annotated phylogenetic tree in newick format
  -o  folder location where the output files will be stored

example:
python spindel.py -a /media/WorkingSpace/Share/indelmip/data/CYP2U_165/CYP2U_165.aln -n  /media/WorkingSpace/Share/indelmip/data/CYP2U_165/CYP2U_annotated_165.nwk -o /media/WorkingSpace/Share/indelmip/data/CYP2U_165/

```

## References
<a id="1">[1]</a> 
Snir et. al. and Patcher et. al. (2006). 
Phylogenetic Profiling of Insertions and Deletions in Vertebrate Genomes. 
https://link.springer.com/chapter/10.1007/11732990_23
