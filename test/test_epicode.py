#!/usr/bin/env python2.7
import os
import sys
import unittest

# path to epicode binary
test_script = os.path.realpath(__file__)
pkg_dir = test_script.rsplit("/", 2)[0]
bin_dir = os.path.join(pkg_dir, "bin")
sys.path.insert(0, bin_dir)

import numpy as np
from numpy.testing import assert_approx_equal as equal
from numpy.testing import assert_array_almost_equal as aequal

from epicode import * # this is bin/epicode.py not epicode/__init__.py

class EpicodeTest(unittest.TestCase):
    
    def test_loadarr(self):
        k,a = load_arr(os.path.join(pkg_dir, "data", "test_ab_cnt.arr"))
        assert k[0] == "dnase:a"
        aequal(a[0], [0.,2.,1.,2.,0.,0.,0.,12.,0.,2.,0.,2.,0.,0.,0.,0.,0.,0.,0.,1.,1.,1.,1.,0.])

    def test_loadarr(self):
        k,a = load_epi(os.path.join(pkg_dir, "data", "absolute_codes.epi"))
        assert k[0] == "h2az"
        equal(a[0][3], 1.30686752)

    def test_sparsevec(self):
        equal(sparsevec(np.array([0,1,1,1,1,2])), 0.22640339557688832)
        equal(sparsevec(np.array([1,1,1,1,1,1])), 4.08e-10, 3)

    def test_sparsemat(self):
        res = sparsemat(np.array([[0,1,1,1,1,2],
                                  [1,1,1,1,1,1]]))
        equal(res, 0.113201697993, 6)

    def test_dsig(self):
        inp = range(100)
        res50 = dsig(inp, 0, 0, 50)
        assert min(res50) >= 0.5
        assert max(res50) < 1.0
        res99 = dsig(inp, 0, 0, 99)
        assert res99[10] - res50[10] < 0
        equal(res99[5], 0.52523106, 6)

    def test_scarr(self):
        inp = np.array([np.arange(100),
                        np.arange(100, 200)]).reshape(100,2)
        inp95 = scarr(inp, "sig95")
        inpwt = scarr(inp, "whiten")
        aequal(inp95[:,0], inp95[:,1])
        equal(inp95[1,1], 0.0105715990067)
        assert (inpwt[:,0] - inpwt[:,1] < 0).all()

    def test_scapair(self):
        inp = np.array([np.arange(100),
                        np.arange(100, 200)]).reshape(100,2)
        meds = np.floor(np.median(scapair(inp, "deseq"), axis=0))
        assert (meds == 99.).all()
        
    def test_common_sub(data):
        data = ["aslkajdlasjd_dontlikedoingthat_budsdijaslkdsa",
                "lkjsalkjdlaksjdlaj_dontlikedoingthat_kjdhsakjdhaskjd"]
        assert common_sub(data) == "_dontlikedoingthat_"

    def test_parse_params(self):
        res = parse_params("rocky:2,web:2.5,epicode:1.0,pynchon:vice", {})
        assert res["pynchon"] == "vice"
        assert res["web"] == 2.5
        assert res["epicode"] == 1.0
        assert res["rocky"] == 2
        assert type(res["epicode"]) == float
        assert type(res["rocky"]) == int

if __name__ == '__main__':
    unittest.main()
