Step 1
------

Install epicode, either from PyPI or from the github repository (this step requires a that your OS has python2.7, with setuptools (or PIP) pre-installed). Other dependencies will be downloaded and compiled if needed, but it is recommended that numpy, scipy, and Cython are also pre-installed by the OS-specific package manager.

```bash
$ git clone https://github.com/mcieslik-mctp/epicode.git
$ cd epicode
$ python2 setpu.py install
```

This step should finish without error and should install any missing dependencies (moke, pysam, scikit-learn, etc).

To verify the installation:

```bash
$ cd ..
$ epicode.py
```

The last command tests if epicode.py can be found in your ```$PATH``` and all dependecies can be imported. If no command is found or epicode.py returns an error please read the detailed installation instructions in ```epicode/README.md``` and if it does not help file an issue at (https://github.com/mcieslik-mctp/epicode)

Step 2
------

Download the necessary files (requires ```xargs``` and ```wget``` both pre-installed on most Linux distributions)

```bash
$ mkdir example
$ cd example
$ cp epicode/data/a549_rep1_etoh.txt .
$ cat a549_rep1_etoh.txt | xargs -n1 -P6 -I '{}' wget -c '{}'
```

This process can take several hours depending on the speed of your internet connection as the total size od the downloaded data is 8G. Six files (``-P6``) are downloaded at the same time. You can resume downloading by running the command again.

Finally, we need a bed file with the regions of interest (here ~19k 1000bp promoters)

```bash
$ cp epicode/data/hsmm_prom_1000.bed . 
```

Step 3
------

Now we can run epicode.py in ```absolute``` mode.

```bash
$ epicode.py absolute -c 6 -bed hsmm_prom_1000.bed -bams *bam -runid encode
```

This step can take some time depending on the speed of your disk. By default 4 processes are used to read the data. This often a good choice, but might slow down your computer considerably.

To visualize the result

```bash
$ Rscript epicode/scripts/absolute_plot.r absolute_out/*epi encode_absolute.png
```

The output should match the one included in the repository at ```epicode/data/encode_absolute.png```.

