#!/usr/bin/env python2
from setuptools import setup, find_packages
from glob import glob

def main():
    setup(name = "epicode",
          url="https://github.com/i000/epicode",
          version = "1.0.0",
          author="Marcin Cieslik",
          author_email="marcin.cieslik@gmail.com",
          packages = find_packages(),
          scripts = glob("bin/*.py"),
          install_requires = ["scikit-learn==0.14.1", "pysam==0.7.5", "moke==1.1.5"]
          )

if __name__ == "__main__":
    main()
