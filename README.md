# rgeneclust
A pipeline to characterize and cluster plant NBS-LRR (NB-LRR) genes, phylogenetic analysis and visualization

## Dependencies
This package is dependent on a list of open source packages:
 * Required packages
   * Perl5 https://www.perl.org
   * Python 2.7 https://www.python.org/downloads/
   * HMMER v3.1 http://hmmer.janelia.org
   * ClustalO v1.1.0 http://www.clustal.org/omega/
   * GNU Parallel http://www.gnu.org/software/parallel/
   * MUSCLE https://www.drive5.com/muscle/
   * PhyML http://www.atgc-montpellier.fr/phyml/
   * Most of these packages are available on MSI and can be checked and loaded by:
```bash
module show perl/python2/hmmer/clustalo/parallel/muscle/phyml
module load perl/python2/hmmer/clustalo/parallel/muscle/phyml
```

 * Required perl modules
   * Bioperl
   * Data::Table
   * List::MoreUtils
   * Time::HiRes
   * Data::Dumper

The BioPerl module has been installed on MSI and can be loaded by:
```bash
module load bioperl
```
Need to set the PERL5LIB environmental variable to include the source folder:
```bash
export export PERL5LIB=$PATH_TO_rgeneclust:$PERL5LIB
```
where *$PATH_TO_rgeneclust* is the absolute path of the source directory

 * Required python packages:
   * numpy http://www.numpy.org
   * pyfasta 0.5 https://pypi.python.org/pypi/pyfasta/

Both python packages can be installed by [pip](https://pip.pypa.io/en/stable/user_guide/#installing-packages):
```bash
pip install pyfasta
pip install numpy
```

By default *pip* will install packages into */usr/local*, but you will probably need to [install into user directories](https://pip.pypa.io/en/stable/user_guide/#user-installs):
```bash
pip install --user pyfasta
```

#Usage
```bash
usage: rosar.py [-h] [--cpu NCPU] cfgfile outdir

Identify, cluster and characterize plant NBS-LRR genes

positional arguments:
  cfgfile     config file (a text file with species identifier followed by the
              absolute path of CDS fasta in each line)
  outdir      output directory

optional arguments:
  -h, --help  show this help message and exit
  --cpu NCPU  number processors to use (default: all/24)
```
Config file (test.csv):
  A text file with species identifier followed by the absolute path of 
  CDS fasta in each line, for example:
    
    Fv,/home/zhoup/test/fvesca_v1.0_genemark_hybrid.fna
    Md,/home/zhoup/test/Malus_x_domestica.v1.0.consensus_CDS.fa
    Pp,/home/zhoup/test/Prunus_persica_v1.0_CDS.fa
  
  with paths of the grape, apple and pear CDS sequences.

