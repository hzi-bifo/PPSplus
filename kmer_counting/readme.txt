K-mer counting algorithm

COMPILE
To compile the program from the source code on a Linux system, run:
sh compile_kmers_linux.sh

(note that this must be run from the directory where the .sh script and directory "kmers" are located)

The result of this step is binary file: fasta2kmers2

Note that this compile command was tested on Debian and Ubuntu, you may need to modify the command for a different platform. 


RUN
To get all the options, with which you can run the software, run:
./fasta2kmers2 -h


PPS+ INTEGRATION
This binary is used in the PPS+ distribution as a subprogram and located in:
/apps/pps/tools/PhyloPythiaS/1_3/bin/fasta2kmers2


DOCUMENTATION
Gregor et al. (2016) PhyloPythiaS+: a self-training method for the rapid reconstruction of low-ranking taxonomic bins from metagenomes. PeerJ. 
(https://peerj.com/articles/1603/)
(See: section "Simultaneous counting of multiple short k-mers", also please see the minor correction)
