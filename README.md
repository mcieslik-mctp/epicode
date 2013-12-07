EpiCODE
=======

epicode.py - Discovers "epigenetic codes" within ChIP-seq datasets.

```
$ epicode.py absolute -bed [BED6+ file] -bams [BAM files] [options]
$ epicode.py differential -bed [BED6+ file] -abams [BAM files] -bbams [BAM files] [options]
$ epicode.py discriminatory -beds [BED6+ files] -bams [BAM files] [options]
```

To get help specific to these three methods see:

```
$ epicode.py {absolute, differential, discriminatory} --help
```

The goal of epicode is to discover patterns of histone modifications from
aligned sequence data. Epicode looks for combinations (subsets) of marks
that tend to occur in (at least) sub-portions of the data. Alternatively
it identifies combinations of marks that change coordinately i.e. are
"gained" or "lost" frequently at the same time.

The algorithm provides three modes ```absolute```, ```discriminatory```, and
```differential```. The first two modes identify co-occurring marks within
one or many sets of genomic losi, respectively. The ```differential``` mode
attempts to find patterns of coordinated mark changes. In ```discriminatory```
mode two (or more) genomic loci are differentiated based their associated
patterns.


Walk-through example
--------------------

To see how epicode could be used with ENCODE data see the full example in.


EpiCODE modes
-------------

Each of the provided modes corresponds to a specific subcommand of epicode.

* ```absolute``` for experiments with multiple histone modifications or 
  epigenetics marks mapped in a single condition. Epicode finds patterns
  (epigenetic codes) of frequently co-occurring marks.

* ```differential``` for experiments with the same marks mapped in two conditions.
  Epicode finds patterns of coordinated mark changes i.e. subsets of marks
  that are often "gained" or "lost" at the same time.

* ```discriminatory``` for experiments where one is interested in the epigenetic
  patterns that distinguish different sets of (preferably non-overlapping) genomic
  loci. Multiple histone modifications are mapped in a single condition and
  quantified for two two or more (experimental) sets of loci.

As input epicode expects at least one BED6+ file of reference genomic regions
(-bed or -beds) and at least one set of aligned sequence reads in coordinate
sorted BAM files. Epicode is not filtering duplicate reads, please run
```samtools dedup``` to create deduplicated input files if this is desired.


EpiCODE tasks
-------------

The three provided high-level modes are wrappers around tasks with more fine-grained options.
A list of all the available tasks can be seen by:

```
$ epicode.py --help
...
| extract_absolute      Processes multiple bam files in "absolute" mode.
|
| discriminatory        Discriminatory "epigenetic codes" that emphasize
|                       epigenetic differences between two (or more
|                       [experimental]) sets of sites.
|
| scale_features        Scales features of any input array (loci x features)
|                       using any of the supported algorithms.
|
| extract_diff          Processes multiple bam files in "differential" mode.
|
| differential          Differential "epigenetic codes" from "gain-loss"
|                       changes in levels of epigenetic marks from two
|                       experimental conditions in single set of sites.
|
| scale_pairs           Scales observed paired columns of read-overlap counts.
|
| code_sklearn          Non-negative matrix factorization using scikits-learn.
|
| multi_code_sklearn    Multi-array Non-negative matrix factorization using
|                       scikits-learn.
|
| scale_diff            Calculates differential features from paired counts
|                       array (paired loci x features).
|
| absolute              Absolute "epigenetic codes" from levels of epigenetic
|                       marks in a single experimental condition and in a
|                       single set of sites.
...
``` 


### High-level interface

#### Absolute mode: ```absolute```

Designed to work on epigenetic marks mapped in one condition and quantified within one type of 
loci. The sites should be provided as a BED6+ file e.g. a promoter file or an enhancer file. The 
sequencing data should be provided as coordinate sorted bam files (e.g. using samtools or novosort). 
The algorithm results (files) are saved into the ```odn``` directory (default ```absolute_out``` and
prefixed with ```runid```, which defaults to an automatically generated and likely unique integer. Column
names within all output files that are data matrices are generated from input BAM filenames and are optionally 
shortened (```--shorten``` option) by removing redundant substrings. 

This is a wrapper for the following chain of tasks, each task saves the generated data as intermediate 
files in a sigle ```run``` directory:

  1. extract_absolute - Extracts counts of reads overlapping genomic regions and 
     normalizes by region length.
  2. scale_features - Scales features (columns) within each ```lvl.arr``` array. 
  3. code_sklearn - Learn discriminatory epigenetic codes.

The procedure creates two files for the two matrices ```{parameters}.arr``` (optionally) and ```{parameters}.epi```.

Example command line:

```
epicode.py absolute -c 6 -par 4 -bed <<repo dir>>/data/h1esc_prom_5000.bed -bams <<bam dir>>/*.bam
```

Here ```<<repo dir>>``` is where you checked out the git reposity and ```<<bam dir>>``` is a directory with bam files. The parameters ```-c 6``` and ```-par 4``` mean that six codes in "absolute" mode will be learned and the bam file processing will happen using four cores.


#### Differential mode: ```differential```

Designed to work on epigenetic marks mapped in two conditions (A and B) quantified in one type of locus.
The sites should be provided as a BED6+ file e.g. a promoter file or an enhancer  file. The sequencing
data should be provided as coordinate sorted bam files (e.g. using samtools or novosort). The algorithm
results (files) are saved into the ```odn``` directory (default ```differential_out``` and prefixed with
```runid```, which defaults to an automatically generated and likely unique integer. Column names within all
output files that are data matrices are generated from input BAM filenames and are optionally shortened
(```--shorten``` option) by removing redundant substrings. 

The input to the ```differential``` mode are two sets of BAM files. Each set should contain several BAM 
files, typically more than 4. Within each set the bam files should have identical names i.e. they should
be the same marks mapped in two conditions.

This is a wrapper for the following chain of tasks, each task saves the generated data as intermediate 
files in a sigle ```run``` directory:

  1. extract_differential - Extracts paired counts within genomic regions in ```step``` resultion.
  2. scale_pairs - Normalizes counts paired samples for sequencing depth.
  3. scale_differential - Converts scaled absolute counts to ```gain-loss``` levels.
  4. scale_features - Scales features (columns) within each ```{parameters}lvl.arr``` array.
  5. code_sklearn - Learn absolute epigenetic codes.
    
The procedure creates two files for the two matrices ```{parameters}.arr``` (optionally) and
```{parameters}.epi```.

Example command line:

```
epicode.py differential -c 8 -par 4 <<repo dir>>data/hsmm_prom_1000.bed -abams <<A bams dir>>/*.bam -bbams <<B bams dir>>/*.bam
```
Here ```<<repo dir>>``` is where you checked out the git reposity and ```<<A bams dir>>``` and  ```<<B bams dir>>``` are directories with BAM files. Epicode assumes that both directories contain BAM files with identical names, such that listing their contents returns identical filenames (basenames). The parameters ```-c 6``` and ```-par 4``` mean that six codes in "differential" mode will be learned and the BAM file processing will happen using four cores.

#### Discriminatory mode: ```discriminatory```

Designed to work on epigenetic marks mapped in a single condition, quantified within two 
(or more) types of loci. The types of sites are provided as BED6+ files e.g. a promoter file and an enhancer 
file. The sequencing data should be provided as coordinate sorted bam files (e.g. using samtools or novosort). 
The algorithm results (files) are saved into the ```odn``` directory and prefixed with ```runid```, which defaults
to ```discriminatory```. Column names within all output files that are data matrices are generated from input BAM 
filenames and are optionally shortened (```--shorten``` option) by removing redundant substrings. 

The algorithm takes three important parameters ```c``` the number of expected histone codes and also rank of the 
factored matrices, ```colsca``` the algorithm used to scale the levels (columns) of the final input matrices
to the NMF algorithm, and ```init``` the algorithm used to initialize matrices. The two latter parameters are
best left as defaults. The ```c``` parameter has no default as it depends both on the number of input bam files, 
redundancy (correlation) of the assayed epigenetic marks and the biological complexity of the genomic regions 
(bed files). Typically a value between 4-10 gives interpretable results, but please see our publication for 
some recommendations and properties of NMF applied to epigenomic data.

This is a wrapper for the following chain of tasks, each task saves the generated data as intermediate files 
in a sigle ```runid``` directory:
    
  1. extract_absolute - Extracts mark lavels ```*lvl.arr``` for all marks (BAM files) and regions (BED files)
     combinations.
  2. scale_features - Scales features (columns) within each ```*lvl.arr``` array.
  3. multi_code_sklearn - Learn discriminatory epigenetic codes.
  
The procedure creates two files for the two matrices ```{parameters}.arr``` (optionally) and
```{parameters}.epi```. The files can be used as input to the logistic regression classifier task
``logistic_classifier``.

```
epicode.py discriminatory -par 4 -beds <<repo dir>>data/pol2_1000_e25_*.bed -c 7 -bams <<bam dir>>/*.bam 
```

Here ```<<repo dir>>``` is where you checked out the git reposity and ```<<bam dir>>``` is a directory with bam files. The parameters ```-c 7``` and ```-par 4``` mean that seven codes in "differential" mode will be learned and the bam file processing will happen using four cores.


### Low-level interface


#### ```code_sklearn```

NMF factorization of asingle array. Creates two new files with ```.epi``` and ```.arr``` suffixes (only with
```--transform```) that include algorithm parameters in their names in the same directory as the input.
Method and remaining arguments are currently ignored from the command line. (see: ``sklearn.decomposition.NMF``).

#### ```multi_code_sklearn```

NMF factorization of multiple arrays. Creates two new files with ```.epi``` and ```.arr``` suffixes (only with
```--transform```) that include algorithm parameters in their names in the same directory as the input.
Method and remaining arguments are currently ignored from the command line. (see: ``sklearn.decomposition.NMF``).

#### ```scale_features```

Scales features of any input array (loci x features) using any of the supported algorithms. Creates new
array file with scaled columns scaling algorithm name is included as suffix. See: ``scarr``function
documentation for details.

#### ```scale_pairs```

Scales observed paired columns of read-overlap counts (paired samples). Currently only the ```deseq``` algorithm is
implemented.

#### ```scale_diff```

Calculates differential features from paired counts array (paired loci x features). A proper input for
this function is obtained by the ``extract_differential`` sub-command. Matches pairs by an ```:a``` and ```:b```
suffix. Creates a new file with ```_lvl.arr``` suffix in the same directory as the input. The gain and
loss columns are suffixed ```:g``` and ```:l```, respectievely. 

#### ```extract_absolute```

It estimates enrichment levels for each mark at each input genomic regions from the BED6+ file. Read count
extraction can be done in parallel (per-chromosome parallelism). The output is saved into the ```odn``` directory
with ```runid``` prefix (optional) and has an extension of ```{runid}_lvl.arr```. If shorten is enabled \
```--shorten```the algorithm will try to shorten the file names when producing column names by removing common
substrings. In the case of erros see the log messages. Program aborts if output files are present.

Example command-line:

```epicode.py extract_absolute -bed test.bed -bams directory_with_bams/*bam -odn test_extract -par 8```


#### ```extract_differential```

For each bam file it counts the number reads overlapping each genomic region from the BED6+ file.
Read count extraction can be done in parallel (per-chromosome parallelism). The output is saved into the
```odn``` directory with ```runid``` prefix (optional) and has an extension of ```{runid}_cnt.arr```. If shorten
is enabled ```--shorten``` the algorithm will try to shorten the file names when producing column names by
removing common substrings. Column names receive a ```:a``` or ```:b``` suffix when they are from tha A-list or
B-list, respectively. In the case of erros see the log messages. Program aborts if output files are present.
Output array is proper input for the ```scale_features``` task.



Plotting
--------

Two scripts are provided to facilitate plotting. They can be found in the ```epicode/scripts``` directory. Each
script take an epigenetic code file ```.epi``` and the name of an output graphics file supported by ggplot e.g.
```png```, ```pdf``` etc.

Example usage:

```
$ Rscript <<repo dir>>/scripts/absolute_plot.r epicode/data/absolute_codes.epi absolute_codes.png
```

```
$ Rscript <<repo dir>>/scripts/gainloss_plot.r epicode/data/differential_codes.epi differential_codes.png
```

The above scripts require the following packages: ```ggplot2```, ```reshape2```, ```stringr```.
All of these can be installed from CRAN by the following command:


```R
install.packages("<<package name>>")
```

For example:

```R
install.packages("ggplot2")
```


Logging Configuration 
---------------------

Epicode.py can be configured to use an external configuration file and to log to an arbitrary file stream

See ```moke``` help for details:

```
$ epicode.py --help
```



Installation
------------

### Automatic Installation

In the simple case installing Epicode requires only:

```bash
$ easy_install-2.7 epicode
```

If the above command is not found first try:

```
$ easy_install epicode
```

If the installation procedure runs, but fails at some point follow the manual installation guide.
If ```easy_install``` is still not found please try installing ```setuptools``` first (see below). 


### Manual Installation

Since the automatic installation procedure failed we have to make sure that we are running the correct 
version of ```Python``` and ```easy_install``` (```setuptools```):

#### Python and Setuptools

```
$ python2 --version
| Python 2.7.5
```

Verify that ```easy_install``` can be found:

```
$ which easy_install
| .../bin/easy_install
```

If you cannot find ```easy_install``` you have to install ```setuptools``` this is best done system-wide,
using the specific mechanisms. We will use ```easy_install``` to install additional dependencies
later on.

Arch:

```
$ pacman -S python2-setuptools
```

Fedora: 

```
$ yum install python-setuptools python-setuptools-devel
```

Ubuntu:

```
$ sudo apt-get install python-setuptools python-dev
```

Alternatively one can try to follow the instructions at:

https://pypi.python.org/pypi/setuptools/0.6c11

#### Additional Dependencies

Epicode has a small number of dependencies. On many systems they will be successfully installed
using the ```setuptools/easy_install``` mechanism. However, we recommend to install ```numpy``` and
```scipy``` using the sytem-wide mechanism. As their compilation is particularly involved. Epicode
was tested and developed for ```numpy-1.7.1``` and ```scipy-0.12.0```, but should also work on
other relatively recent releases.

For Arch linux:

```
$ pacman -S extra/python2-numpy community/python2-scipy 
```

Fedora: 

```
$ yum install numpy scipy
```

Ubuntu:

```
$ sudo apt-get install python-numpy python-scipy
```

For other operating systems follow the instructions at:

http://www.scipy.org/install.html 

Next, we will install ```pysam```, ```scikit-learn```, and ```moke``` from PyPI:

```
$ easy_install pysam==0.7.5 scikit-learn==0.14.1 moke==1.1.5
```

If all the commands returned correctly you should be able to start ```python```:

```
$ python2
```

And issue the following statements:

```python
>>> import numpy
>>> import scipy
>>> import sklearn
>>> import moke
>>> import pysam
```

Now you are ready to install epicode

```
$ easy_install epicode
```
