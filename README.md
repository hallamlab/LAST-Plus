# LAST+: Optimized threading  for fast annotation

Dongjae Kim, Aria S. Hahn. Niels W. Hanson, Kishori M. Konwar, and Steven J. Hallam

Comparative genomic research relies heavily on protein alignment to infer the metabolic potential of an organism or community. Local Alignment Search Tool LAST, uses adaptive seed lengths to gain approximately 50 times the speed of standard seed-and-extend algorithm BLAST. However, LAST is currently capable of utilizing a single CPU. Here, we present LAST+, a multi-threaded, IO optimized implementation of LAST that uses a model for algorithmic thread synchronization to allow efficient use of multiple CPU processes while requiring only 4-8GB of memory. We demonstrate that LAST+ is approximately 6-times faster than LAST using real world environmental data with > 50,000$ sequences, and show that LAST+ is over 2-times faster than DIAMOND, the fastest available aligner prior to LAST+. Finally, we implement new features such as the calculation of BLAST like e-values, and tabular format outputs compatible with multiple downstream analyses tools. 

More information can be found on the [Wiki](https://github.com/hallamlab/LAST-Plus/wiki).

 Test datasets used to produce the claimed results can be found on Dropbox https://www.dropbox.com/sh/phgpjur66fhvej3/AAAgpIfERMdGXux5yXhPQL-Ta?dl=0.
