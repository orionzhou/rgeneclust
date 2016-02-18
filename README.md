# rgeneclust
A pipeline to characterize and cluster plant NBS-LRR (NB-LRR) genes, phylogenetic analysis and visualization

Usage:
  python rosar.py cfg.csv --out test

Config file (cfg.csv):
  A text file with species identifier followed by the absolute path of 
  CDS fasta in each line, for example:
    
    Fv,/home/zhoup/test/fvesca_v1.0_genemark_hybrid.fna
    Md,/home/zhoup/test/Malus_x_domestica.v1.0.consensus_CDS.fa
    Pp,/home/zhoup/test/Prunus_persica_v1.0_CDS.fa
  
  with paths of the grape, apple and pear CDS sequences.

Output:
  test/
