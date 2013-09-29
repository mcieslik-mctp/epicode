#!/usr/bin/env python2.7
from moke import *

import numpy as np
from numpy.testing import assert_approx_equal as equal
from numpy.testing import assert_array_almost_equal as aequal


@task
def test_scale_pairs():
    mkdir(path("test_scale_pairs"))
    assert sh("cp ../data/test_ab_cnt.arr test_scale_pairs") == 0
    assert sh("epicode.py scale_pairs test_scale_pairs/test_ab_cnt.arr") == 0
    with open("test_scale_pairs/test_ab_cnt_deseq.arr") as fh:
        head = fh.readline()
        firs = fh.readline()
        equal(float(firs.split("\t")[1]), 1.632993161855451847e+00)
    assert sh("rm -rf test_scale_pairs") == 0

@task
def test_scale_diff():
    mkdir(path("test_scale_diff"))
    assert sh("cp ../data/test_ab_cnt_deseq.arr test_scale_diff") == 0
    assert sh("epicode.py scale_diff test_scale_diff/test_ab_cnt_deseq.arr") == 0
    with open("test_scale_diff/test_ab_cnt_deseq_lvl.arr") as fh:
        head = fh.readline()
        firs = fh.readline()
        equal(float(firs.split("\t")[1]), 3.635961336943780431e-01)
    assert sh("rm -rf test_scale_diff") == 0

@task
def test_scale_features():
    mkdir(path("test_scale_features"))
    assert sh("cp ../data/test_ab_cnt_deseq_lvl.arr test_scale_features") == 0
    assert sh("epicode.py scale_features -scalgo sig95 test_scale_features/test_ab_cnt_deseq_lvl.arr") == 0
    with open("test_scale_features/test_ab_cnt_deseq_lvl_sig95.arr") as fh:
        head = fh.readline()
        firs = fh.readline()
        equal(float(firs.split("\t")[1]), 1.211047172546386719e-02)
    assert sh("rm -rf test_scale_features") == 0

@task
def test_code_sklearn():
    mkdir(path("test_code_sklearn"))
    assert sh("cp ../data/a549_start_lvl_sig95.arr test_code_sklearn") == 0
    assert sh("epicode.py code_sklearn -c 6 test_code_sklearn/a549_start_lvl_sig95.arr") == 0
    with open("test_code_sklearn/a549_start_lvl_sig95_pgnmf-c#6-i#None-p#.epi") as fh:
         head = fh.readline()
         firs = fh.readline()
         equal(float(firs.split("\t")[0]), 4.0693399501)
    assert sh("rm -rf test_code_sklearn") == 0

@task
def test_multi_code_sklearn():
    mkdir(path("test_multicode_sklearn"))
    assert sh("cp ../data/{0,1}*sig95.arr test_multicode_sklearn") == 0
    assert sh("epicode.py multi_code_sklearn -base test_multicode_sklearn/tve -c 6 test_multicode_sklearn/{0,1}*sig95.arr") == 0
    with open("test_multicode_sklearn/tve_pgnmf-c#6-i#None-p#.epi") as fh:
         head = fh.readline()
         firs = fh.readline()
         equal(float(firs.split("\t")[1]), 2.75326630104)
    with open("test_multicode_sklearn/tve_pgnmf-c#6-i#None-p#.arr") as fh:
         head = fh.readline()
         firs = fh.readline()
         equal(float(firs.split("\t")[0]), 1.250747926106863528e-01)
    assert sh("rm -rf test_multicode_sklearn") == 0

@task
def test_recode_sklearn():
    mkdir(path("test_recode_sklearn"))
    assert sh("cp ../data/0*sig95.arr test_recode_sklearn") == 0
    assert sh("cp ../data/tss_vs_enh_pgnmf-c#6-i#None-p#.epi test_recode_sklearn") == 0
    assert sh("epicode.py recode_sklearn -arr test_recode_sklearn/0_tss_vs_enh_lvl_sig95.arr -epi test_recode_sklearn/tss_vs_enh_pgnmf-c#6-i#None-p#.epi -base recode -odn test_recode_sklearn") == 0
    with open("test_recode_sklearn/recode.arr") as fh:
        head = fh.readline()
        firs = fh.readline()
        equal(float(firs.split("\t")[0]), 1.250747926106863528e-01)
    assert sh("rm -rf test_recode_sklearn") == 0


if __name__ == "__main__":
    task()
