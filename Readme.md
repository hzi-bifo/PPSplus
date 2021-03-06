### Method description

https://peerj.com/articles/1603/

### Installation instructions, Tutorial

https://github.com/hzi-bifo/PPSplus/tree/master/doc 
(https://github.com/algbioi/ppsp)

### PPS+ Docker image with instructions and latest reference database (updated 2019)

https://hub.docker.com/r/cami/ppsp

https://osf.io/heksw/

### Benchmark Datasets

https://github.com/algbioi/datasets

### Repositories (source code)

PPS+ python scripts:
https://github.com/hzi-bifo/PPSplus/tree/master/ppsplus
(https://github.com/algbioi/ppsplus)

The k-mer counting algorithm:
https://github.com/hzi-bifo/PPSplus/tree/master/kmer_counting
(https://github.com/algbioi/kmer_counting)

### (PPS+ VM - alternative installation instructions)

(https://research.bifo.helmholtz-hzi.de/software/ppsp/readme_vm_release.txt)

### This section describes, how to use additional marker genes in the analysis:

To add a marker gene, you need to provide three files:
* Model HMM: a file generated by the _HMMER 3.0_ software from a particular protein seed alignment that contains sequences characteristic for a particular marker gene.
* Template FASTA: a FASTA file containing DNA reference sequences.
* Taxonomy file: a tab separated mapping file containing a taxonomic identifier for each sequence from (Template FASTA), first column ~ sequence name, second column ~ NCBI taxonomic identifier.

Rename the files, for a marker gene _X_, use the same prefix, e.g.:
* _marker_gene_X.hmm_
* _marker_gene_X.fna_
* _marker_gene_X.tax_

Add the files in the following steps:

1. Locate the current reference that you are using for _PPS+_, e.g. _reference_NCBI201502_. 
2. Go to subdirectory starting with "_mg_" (e.g. _mg5_). 
3. Put the _marker_gene_X.hmm_ to directory that starts with "_hmm_" (e.g. _hmm_). 
4. Put the _marker_gene_X.fna_ and _marker_gene_X.tax_ files to directory starting with "_db_" (e.g. _db_). 
5. Edit tab separated file "_content.csv_", add, analogously, four lines, e.g.:

* _marker_gene_X     taxonomyDNA     db/marker_gene_X.tax_
* _marker_gene_X     templateDNA     db/marker_gene_X.fna_
* _marker_gene_X     hmmPROTPrim     hmm/marker_gene_X.hmm_
* _marker_gene_X     hmmPROTSec      hmm/marker_gene_X.hmm_

Now, when you run _PPS+_ with option "_-g_", it will also consider the newly added marker genes.

Note that it is also possible to remove marker genes from the analysis by removing the corresponding files from the reference directories and by editing "_content.csv_").

### News
* **The data was moved to a different server, please use this link to download the VM:** https://research.bifo.helmholtz-hzi.de/software/ppsp/1_4/ppsp_1_4_vm_64bit.ova. Before running the "update" command, run "wget https://raw.githubusercontent.com/hzi-bifo/PPSplus/master/ppsplus/algbioi/ref/update.py?token=AHgcgMzH02drqhFnu-TsSVLNlgcuL9fgks5Z0RLcwA%3D%3D -O /apps/pps/sys/update.py" - this will tell the VM where the source and reference files are located.
* PPS+ version 1.5c in a Docker container is in preparation, see: https://hub.docker.com/r/ivan737/ppsp/
* There is an updated reference database available. Run "update -r NCBI201502" to get the new reference. The paths in the configuration file has to be set accordingly.

### Note
The source of this repository was:
https://github.com/algbioi/ppsp/wiki
