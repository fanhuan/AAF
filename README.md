AAF (Alignment and Assembly Free)
===

This is a package for constructing phylogeny without doing alignment or assembly.

###Bootstrap
The most feedback I received about AAF are around bootsrap. It is very computationally intensive to do the two-step nonparametric bootstrap. In case you have a higher coverage (>8X), we assume that the incomplete coverage problem is minor. To reduce the computational load, you can choose only to carry out the seconde step of the bootstrap (nonparametric\_bootstrap\_s2only.py): sample the kmer table with replacement 1/k of the number of the rows of the table. To further reduce the computaiton, here is a version to sample from the shared kmer table (nonparametric\_bootstrap\_s2only_skt.py). Singletons from each sample (i.e. kmers that only appear in one sample) are calculated from the difference between the total diversity file and the shared kmer table. Then those singletons are added back during the the calculation of pariwise distance, following a poisson distribution with a mean of 1/k of each singletone number.  
####nonparametric\_bootstrap\_s2only.py
####nonparametric\_bootstrap\_s2only_skt.py


