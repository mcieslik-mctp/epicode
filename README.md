EpiCODE
=======

XXX

Installation
------------

First, please make sure you are running the correct version of Python and easy_install:

$ python2 --version
| Python 2.7.5

Epicode has a number of dependencies. The correct way of installing them is distribution specific.
At minimum one should install "numpy" and "scipy" using the sytem-wide mechanism. Epicode was 
tested and developed for "numpy-1.7.1" and "scipy-0.12.0", but should work on others.

For some common linux distributions follow the instructions at: http://www.scipy.org/install.html 

For Arch linux:

$ pacman -S extra/python2-numpy community/python2-scipy 

To install additional dependencies we will use easy_install. This toos is part of the "setuptools" package, 
which is provided on most systems, but if not it should be installed system wide:

Fedora: 
$ yum install python-setuptools python-setuptools-devel

Ubuntu:
$ sudo apt-get install python-setuptools python-dev

Arch:
$ pacman -S python2-setuptools

Verify that "easy_install" can be found:

$ which easy_install
| .../bin/easy_install

We will install "pysam", "scikit-learn", and "moke" from PyPI:

$ easy_install pysam==0.7.5 scikit-learn==0.14.1 moke==1.1.5

If all the commands returned correctly you should be able to start python:

$ python2
>>> import numpy
>>> import scipy
>>> import sklearn
>>> import moke
>>> import pysam

Now you are ready to install epicode

$ easy_install epicode

