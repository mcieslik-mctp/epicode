EpiCODE
=======

epicode.py - discover epigenetic "codes" from ChIP-seq data.

The goal of epicode is to discover patterns of histone modifications from
aligned sequence data. Epicode looks for combinations (subsets) of marks
that tend to occur in (at least) sub-portions of the data. Alternatively
it identifies combinations of marks that change coordinately i.e. are
"gained" or "lost" frequently at the same time.

The algorithm provides three modes "absolute", "discriminatory", and
"differential". The first two modes identify co-occurring marks within
one or many sets of genomic losi, respectively. The "differential" mode
attempts to find patterns of coordinated mark changes. In "discriminatory"
mode two (or more) genomic loci are differentiated based their associated
patterns.

EpiCODE modes
-------------

Each of the following modes corresponds to a specific subcommand of epicode.
To get help specific to these three methods see:

```bash
$ epicode.py {absolute, differential, discriminatory} --help
```

* "absolute" for experiments with multiple histone modifications or 
  epigenetics marks mapped in a single condition. Epicode finds patterns
  (epigenetic codes) of frequently co-occurring marks.

* "differential" for experiments with the same marks mapped in two conditions.
  Epicode finds patterns of coordinated mark changes i.e. subsets of marks
  that are often "gained" or "lost" at the same time.

* "discriminatory" for experiments where one is interested in the epigenetic
  patterns that distinguish different sets of (preferably non-overlapping) genomic
  loci. Multiple histone modifications are mapped in a single condition and
  quantified for two two or more (experimental) sets of loci.

As input epicode expects at least one BED6+ file of reference genomic regions
(-bed or -beds) and at least one set of aligned sequence reads in coordinate
sorted BAM files. Epicode is not filtering duplicate reads, please run
``samtools dedup`` to create deduplicated input files.




```bash
epicode.py absolute -bed [BED6+ file] -bams [BAM files] [options]
epicode.py differential -bed [BED6+ file] -abams [BAM files] -bbams [BAM files] [options]
epicode.py discriminatory -beds [BED6+ files] -bams [BAM files] [options]
```



Installation
------------


### Automatic Installation

In the simple case installing Epicode requires only:

```bash
$ easy_install-2.7 epicode
```

If the above command is not found first try:

```bash
$ easy_install epicode
```

If the installation procedure runs, but fails at some point follow the manual installation guide.
If ```easy_install``` is still not found please try installing ```setuptools``` first (see below). 


### Manual Installation

Since the automatic installation procedure failed we have to make sure that we are running the correct 
version of ```Python``` and ```easy_install``` (```setuptools```):

#### Python and setuptools

```bash
$ python2 --version
| Python 2.7.5
```

Verify that ```easy_install``` can be found:

```bash
$ which easy_install
| .../bin/easy_install
```

If you cannot find ```easy_install``` you have to install ```setuptools``` this is best done system-wide,
using the specific mechanisms. We will use ```easy_install``` to install additional dependencies
later on.

Arch:

```bash
$ pacman -S python2-setuptools
```

Fedora: 

```bash
$ yum install python-setuptools python-setuptools-devel
```

Ubuntu:

```bash
$ sudo apt-get install python-setuptools python-dev
```

Alternatively one can try to follow the instructions at:

https://pypi.python.org/pypi/setuptools/0.6c11

#### Installation of Dependencies

Epicode has a small number of dependencies. On many systems they will be successfully installed
using the ```setuptools/easy_install``` mechanism. However, we recommend to install ```numpy``` and
```scipy``` using the sytem-wide mechanism. As their compilation is particularly involved. Epicode
was tested and developed for ```numpy-1.7.1``` and ```scipy-0.12.0```, but should also work on
other relatively recent releases.

For Arch linux:

```bash
$ pacman -S extra/python2-numpy community/python2-scipy 
```

Fedora: 

```bash
$ yum install numpy scipy
```

Ubuntu:

```bash
$ sudo apt-get install python-numpy python-scipy
```

For other operating systems follow the instructions at:

http://www.scipy.org/install.html 

Next, we will install ```pysam```, ```scikit-learn```, and ```moke``` from PyPI:

```bash
$ easy_install pysam==0.7.5 scikit-learn==0.14.1 moke==1.1.5
```

If all the commands returned correctly you should be able to start ```python```:

```bash
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

```bash
$ easy_install epicode
```
