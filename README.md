EpiCODE
=======

epicode.py - discover epigenetic "codes" from ChIP-seq data.

The goal of epicode is to discover patterns of histone modifications.
We are looking for subsets of marks that tend to occur in sub-portions 
of the data ["absolute" and "discriminatory" modes] or coordinately 
change ("gain" or "loss" at the same time) ["differential" mode]. 
    
The algorithm finds frequently co-occurring or coordinately changing marks. 
In addition it is possible to differentiate genomic loci based their 
associated patterns.
    
Epicode provides three modes of operation:

* "absolute" for experiments with multiple histone modifications or 
epigenetics marks mapped in a single condition. Epicode finds "codes" 
of frequently co-occurring marks. 
* "differential" for experiments with the same marks mapped in two conditions.
Epicode finds patterns of coordinated marke changes i.e. subsets of marks
that are often gained or lost together.
* "discriminatory" for experiments where one is interested in the features
that distinguish two sets of genomic loci. Multiple histone modifications 
are mapped in a single condition and quantified for two sets of loci.

As input it expects a BED6+ files of reference genomic regions (-bed or -beds)
and one ("absolute", "discriminatory") or two "differential" sets of aligned 
sequence reads in sorted BAM files.

```bash
epicode.py absolute -bed <<BED6+ file>> -bams <<BAM files>> [options]

epicode.py differential -bed <<BED6+ file>> -abams <<BAM files>> -abams <<BAM files>> [options]

epicode.py discriminatory -beds <<BED6+ files>> -bams <<BAM files>> [options]
```


To get help specific to the two methods see:

```bash
epicode.py {absolute, differential, discriminatory} --help
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
If "easy_install" is still not found please try installing "setuptools" first (see below). 


### Manual Installation

Since the automatic installation procedure failed we have to make sure that we are running the correct 
version of "Python" and easy_install (setuptools):

```bash
$ python2 --version
| Python 2.7.5
```

Verify that "easy_install" can be found:

```bash
$ which easy_install
| .../bin/easy_install
```


Epicode has a number of dependencies. The correct way of installing them is distribution specific.
At minimum one should install "numpy" and "scipy" using the sytem-wide mechanism. Epicode was 
tested and developed for "numpy-1.7.1" and "scipy-0.12.0", but should also work on other recent
releases.

For some common linux distributions follow the instructions at: http://www.scipy.org/install.html 

For Arch linux:

```bash
$ pacman -S extra/python2-numpy community/python2-scipy 
````

To install additional dependencies we will use easy_install. This toos is part of the "setuptools" package, 
which is provided on most systems, but if not it should be installed system wide:

Fedora: 

```bash
$ yum install python-setuptools python-setuptools-devel
```

Ubuntu:

```bash
$ sudo apt-get install python-setuptools python-dev
```

Arch:

```bash
$ pacman -S python2-setuptools
```


We will install "pysam", "scikit-learn", and "moke" from PyPI:

```bash
$ easy_install pysam==0.7.5 scikit-learn==0.14.1 moke==1.1.5
```



If all the commands returned correctly you should be able to start python:

```python
$ python2
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
