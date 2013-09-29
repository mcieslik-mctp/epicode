#!/usr/bin/env python2
from moke import *
from itertools import izip, chain
from multiprocessing import Pool

import pickle
import numpy as np
import scipy.stats as ss
from sklearn import decomposition, cross_validation, grid_search, linear_model, metrics
from sklearn.decomposition.nmf import nnls
from pysam import Samfile


def load_epi(epi):
    """(internal) load epi file
    """
    chk_exit(*inp_file(path(epi)))
    with open(epi) as fh:
        marks = fh.readline().strip().split("\t")
        h = np.loadtxt(fh, delimiter="\t")
    return (marks, h)

def load_arr(arr):
    """(internal) load arr file
    """
    chk_exit(*inp_file(path(arr)))
    with open(arr) as fh:
        marks = fh.readline().strip().split("\t")
        x = np.loadtxt(fh, delimiter="\t")
    return (marks, x)

def load_arrs(arrs):
    """(internal) load multiple arr files, assumes same marks (columns)
    """
    xs = []
    for arr in arrs:
        marks, x = load_arr(arr)
        xs.append(x)
    return (marks, xs)

## statistical functions

def sparsevec(x):
    """(internal) Calculates the sparsity of a vector.
    """
    eps = np.finfo(x.dtype).eps if 'int' not in str(x.dtype) else 1e-9
    n = x.shape[0]
    x1 = np.sqrt(n) - (np.abs(x).sum() + eps) / (np.sqrt(np.multiply(x, x).sum()) + eps)
    x2 = np.sqrt(n) - 1
    return x1 / x2 

def sparsemat(X):
    """(internal) Calculates the average sparsity of a matrx.
    """
    return np.mean([sparsevec(x) for x in X])

def dsig(a, lq, loc, uq):
    """(internal) Double sigmoid function to normalize features (columns).
    see:

    Jain, A., Nandakumar, K. and Ross, A. Score normalization in multimodal biometric systems. 
    Pattern Recognition 38, 2270-2285 (2005).

    """
    a = np.asanyarray(a, dtype="f")
    alpha_l = loc - lq
    alpha_r = uq - loc
    a = a - loc
    lsel = (a < 0.)
    rsel = (a >= 0.)    
    if alpha_l:
        a[lsel] = np.divide(a[lsel], -0.5 * alpha_l)
    if alpha_r:
        a[rsel] = np.divide(a[rsel], -0.5 * alpha_r)
    np.exp(a, a)
    np.add(a, 1, a)
    np.power(a, -1, a)
    return a

def scarr(arr, method):
    """(internal) Normalizes features of array (samples x features) through sigmoid scaling or whitening. 
    """
    if method.startswith("sig"):
        hi = float(method.split("sig")[1]) / 100
        data = np.array(arr, dtype=float).T
        qs = ss.mstats.mquantiles(data, (0.0, 0.0, hi), axis=1).T
        for row, lq, mu, uq in izip(data, qs[0], qs[1], qs[2]):
            row[:] = (dsig(row, lq, mu, uq) - 0.5) * 2.
    elif method == "whiten":
        data = np.array(arr, dtype=float).T
        dev = np.std(data, axis=1, ddof=1)[np.newaxis].T
        dev[dev == 0.] = np.nan
        data /= dev
    else:
        raise ValueError("unknown method")
    return data.T


def scapair(raw, method):
    """(internal) Normalizes paired featues (columns) of array (sample x features) through (currently) 
    the DESeq method. It computes size factors by adjusting medians. see: 
    
    Anders, S. & Huber, W. Differential expression analysis for sequence count data. 
    Genome Biology 11, R106 (2010).

    """

    def size_factors(counts):
        counts = counts[np.alltrue(counts, axis=1)]
        logcounts = np.log(counts)
        loggeommeans = np.mean(logcounts, axis=1).reshape(len(logcounts), 1)
        sf = np.exp(np.median(logcounts - loggeommeans, axis=0)) 
        return sf

    if method == "deseq":
        sel = np.alltrue(raw != -1, axis=1)
        scaled = np.array(raw, dtype=float)
        scaled[raw == -1] = np.nan
        tmp = raw[sel]
        for col1 in xrange(0, tmp.shape[1], 2):
            pair = tmp[:, col1:col1+2]
            sf = size_factors(pair)
            scaled[sel, col1:col1+2] = pair / sf
    else:
        raise ValueError("unknown method")

    return scaled


## helper functions

def run_par(fn, args, par=4):
    """(internal) Applies function onto argument list in parallel (multiprocessing).
    """
    pool = Pool(par)
    results = pool.map_async(fn, args).get()
    return results

def common_sub(data):
    """(internal) Finds longest common substring for a list of strings.
    """
    substr = ''
    if len(data) > 1 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0])-i+1):
                if j > len(substr) and all(data[0][i:i+j] in x for x in data):
                    substr = data[0][i:i+j]
    return substr

def parse_bed(fn):
    """(internal) Parses a BED6+ file.
    """
    regions = []
    with open(fn) as fh:
        for line in fh:
            fields = line.strip().split("\t")
            fields[1:3] = map(int, fields[1:3])
            bed6 = fields[:6]
            if fields[1] < 0:
                # samtools does not deal well with negative indices
                continue
            regions.append(bed6)
    return regions

def parse_params(params, kwargs):
    """(internal) Tries to guess the data type (int>float>string) of values in 
    a parameter dict.
    """
    if params:
        for param in params.split(","):
            k, v = param.split(":")
            try:
                v = int(v)
            except ValueError:
                try:
                    v = float(v)
                except ValueError:
                    pass
            kwargs[k] = v
    return kwargs

def write_codes(fn, nmf, bam_names):
    """(internal) Writes codes to file.
    """
    with open(fn, "wb") as wh:
        wh.write("\t".join(bam_names) + "\n")
        for row in nmf:
            wh.write("\t".join(map(str, list(row))) + "\n")

def write_values(fn, values, c, header=None):
    """(internal) Writes values (or any numeric array with optional header) to file.
    """
    with open(fn, "wb") as wh:
        header = "\t".join(header or ["c%s" % (cc + 1) for cc in range(c)])
        wh.write(header + "\n")
        np.savetxt(wh, values, delimiter="\t")


## absolute mode functions

def parse_bam_absolute(fn, regs):
    """(internal) Parses bam file in absolute mode. Proceeds by counting reads mapping 
    onto a segment (chr, start, end) and normalizes the count by the segment's length.
    """
    bam = Samfile(str(fn), "rb")
    count = []
    for reg in regs:
        chr, start, end = reg[:3]
        n = bam.count(chr, start, end)
        count.append(float(n) / (end - start))
    return count

def parse_bam_absolute_star(fn_regs):
    """(internal) unpack and evaluate
    """
    fn, regs = fn_regs
    return parse_bam_absolute(fn, regs)

def process_bam_absolute(bams, regs, shorten, par):
    """(internal) Processes multiple bam files in parallel using ``parse_bam_absolute``.
    """
    names = [bam_file.basename().splitext()[0] for bam_file in bams]
    fx = common_sub(names)
    if shorten and len(fx) > 6:
        names = [bam.replace(fx, "") for bam in names]
    args = [(bam, regs) for bam in bams]
    tmp = run_par(parse_bam_absolute_star, args, par=par)
    counts = np.column_stack(tmp)
    return (names, counts)

@task
def extract_absolute(bed=None, bams=None, odn=None, runid=None, shorten=False, par=None):
    """Processes multiple bam files in "absolute" mode. 

     - bed(``path``) input genomic regions in the BED6+ file format 
     - bams(``path+``) input sequencing data in sorted BAM files requiers BAI index files
     - odn(``path``) output directory name
     - runid(``str``) run id to prefix all output files
     - shorten(``bool``) truncate BAM file names to unambigous strings
     - par(``int``) number of parallel processes for bam extraction

    """
    # checks input
    chk_exit(bed is None, "a BED6+ file is required")
    chk_exit(bams is None, "a set of BAM files is required")
    chks([inp_file(bam) for bam in bams])
    mkdir(odn)
    
    # run id
    if not runid:
        runid = str(hash((bed,tuple(bams))))
    log("runid: %s" % runid)

    # process bed
    bed_regions = parse_bed(bed)
    log("number of query regions: %s" % len(bed_regions))

    # process bams
    chk_exit(bool(odn.listdir("%s_lvl.arr" % runid)), "error: %s exists in %s" % (runid, odn))
    names, counts = process_bam_absolute(bams, bed_regions, shorten, par)
    log("bam number: %s" % len(bams))
    log("bam names: %s" % ", ".join(names))
    fn = odn / (runid + "_%s.arr" % "lvl")
    with open(fn, "wb") as wh:
        wh.write("\t".join(names) + "\n")
        np.savetxt(wh, counts, delimiter="\t")
    log("saved: %s" % fn)
    return fn


## differential mode functions

def parse_bam_differential(afn, bfn, regs, step):
    """(internal) Parses bam file in absolute mode. Proceeds by counting reads mapping 
    onto a segment (chr, start, end). No normalization is done at this step.
    """
    abam = Samfile(str(afn), "rb")
    bbam = Samfile(str(bfn), "rb")
    acount = []
    bcount = []
    oldchr = "chr1"
    for reg in regs:
        chr, start, end = reg[:3]
        if chr != oldchr:
            log("files: %s - %s : %s counted" % (afn, bfn, oldchr))
            oldchr = chr
        # this could be improved
        for s in xrange(start, end, step):
            e = s + step
            an = abam.count(chr, s, e)
            bn = bbam.count(chr, s, e)
            acount.append(an)
            bcount.append(bn)
        acount.append(-1)
        bcount.append(-1)
    log("files: %s - %s : %s counted (finished)" % (afn, bfn, oldchr))
    return acount, bcount

def parse_bam_differential_star(afn_bfn_regs_step):
    """(internal) unpack and evaluate
    """
    afn, bfn, regs, step = afn_bfn_regs_step
    return parse_bam_differential(afn, bfn, regs, step)

def process_bam_differential(abams, bbams, regs, shorten, par, step):
    """(internal) Processes multiple paired bam files in parallel using ``parse_bam_differential``. 
    Bam files are expected to have the same ``basenames``.
    """
    anames = [bam_file.basename().splitext()[0] for bam_file in abams]
    bnames = [bam_file.basename().splitext()[0] for bam_file in bbams]
    anames = [a.split("_", 1)[0] for a in anames]
    bnames = [b.split("_", 1)[0] for b in bnames]
    assert (anames == bnames)
    if shorten:
        fx = common_sub(anames + bnames)
        if len(fx) > 6:
            anames = [a.replace(fx, "") for a in anames]
            bnames = [b.replace(fx, "") for b in bnames]
    anames = [a + ":a" for a in anames]
    bnames = [b + ":b" for b in bnames]
    args = [(abam, bbam,regs, step) for abam,bbam in zip(abams, bbams)]
    tmp = run_par(parse_bam_differential_star, args, par=par)
    tmp = list(chain(*tmp))
    counts = np.column_stack(tmp)
    names = list(chain(*zip(anames, bnames)))
    return (names, counts)


@task
def recode_sklearn(arr=None, epi=None, odn=path("."), base=None):
    """(internal) projects arr onto codes 

     - arr(``path``) 
     - epi(``path``)

    """
    arr_marks, X = load_arr(arr)
    epi_marks, H = load_epi(epi)
    assert arr_marks == epi_marks
    W = np.zeros((X.shape[0], len(H)))
    for j in range(0, X.shape[0]):
        W[j, :], _ = nnls(H.T, X[j, :])
    base = base or arr.basename().splitext()[0] + "_" + epi.basename().splitext()[0]
    ofn = odn / (base + ".arr")
    # write 
    write_values(ofn, W, W.shape[1])


@task
def extract_diff(bed=None, abams=None, bbams=None, odn=None, runid=None, shorten=False, step=None, par=None):
    """Processes multiple bam files in "differential" mode. 

     - bed(``path``) input genomic regions in the BED6+ file format 
     - abams(``path+``) sample A sequencing data in sorted BAM files requiers BAI index files
     - bbams(``path+``) sample B sequencing data in sorted BAM files requiers BAI index files
     - odn(``path``) output directory name
     - runid(``str``) run id to prefix all output files
     - shorten(``bool``) truncate BAM file names to unambigous strings
     - par(``int``) number of parallel processes for bam extraction

    """
    # checks input
    chk_exit(*inp_file(bed))
    chks([inp_file(bam) for bam in abams + bbams])
    abams = tuple(sorted(abams))
    bbams = tuple(sorted(bbams))
    mkdir(odn)
    if not runid:
        runid = str(hash((bed, abams, bbams)))
    log("runid: %s" % runid)

    # process bed
    bed_regions = parse_bed(bed)
    log("number of query regions: %s" % len(bed_regions))

    # process bams
    chk_exit(bool(odn.listdir("%s_cnt.arr" % runid)), "error: %s exists in %s" % (runid, odn))
    names, counts = process_bam_differential(abams, bbams, bed_regions, shorten, par, step)
    log("bam pair number: %s" % len(abams))
    log("bam names: %s" % ", ".join(names))
    fn = odn / (runid + "_%s.arr" % "cnt")
    with open(fn, "wb") as wh:
        wh.write("\t".join(names) + "\n")
        np.savetxt(wh, counts, fmt="%d", delimiter="\t")
    log("saved: %s" % fn)
    return fn


@task
def scale_pairs(arr, scalgo="deseq"):
    """Scales observed paired columns of read-overlap counts.

     - arr(``path``) input array regions x (markX in A, markX in B, markY in A, markY in B ...) 
     - scalgo(``str``) scaling algorithm

    """
    chk_exit(*inp_file(path(arr)))
    with open(arr) as fh:
        names = fh.readline().strip().split("\t")
        raw = np.loadtxt(fh, delimiter="\t")

    scaled = scapair(raw, scalgo)

    ofn = arr.replace(".arr", "_%s.arr" % (scalgo,))
    with open(ofn, "wb") as wh:
        wh.write("\t".join(names) + "\n")
        np.savetxt(wh, scaled, delimiter="\t")
    log("saved: %s" % ofn)
    return ofn

@task
def scale_diff(arr):
    """Calculates differential features from paired counts array (paired loci x features). 

     - arr(``path``) input array sites x marks

    """
    chk_exit(*inp_file(path(arr)))
    with open(arr) as fh:
        names = fh.readline().strip().split("\t")
        scaled = np.loadtxt(fh, delimiter="\t")

    # python does not have rreplace like rsplit and a: or b: might be in the name
    glnames = [n[::-1].replace("a:", "g:", 1).replace("b:", "l:", 1)[::-1] for n in names] 

    acols = scaled[:,0::2]
    bcols = scaled[:,1::2]

    gl = [[] for _ in xrange(scaled.shape[1])] # gain loss columns
    dcols = bcols - acols
    i = 0
    while i < dcols.shape[0]:
        j = 0
        for col in gl:
            col.append(0.0)
        while True:
            row = dcols[i]
            i += 1
            j += 1
            if np.isnan(row).any():
                break
            for c, v in enumerate(row):
                if v > 0:
                    gl[(c*2) + 0][-1] += v
                if v < 0:
                    gl[(c*2) + 1][-1] -= v
        for col in gl:
            col[-1] /= float(j)
        
    gla = np.column_stack(gl)
    ofn = arr.replace(".arr", "_%s.arr" % ("lvl",))
    with open(ofn, "wb") as wh:
        wh.write("\t".join(glnames) + "\n")
        np.savetxt(wh, gla, delimiter="\t")
    log("saved: %s" % ofn)
    return ofn


@task
def scale_features(arr, scalgo=None):
    """Scales features of any input array (loci x features) using any of the supported algorithms. 
    
     - arr(``path``) input array regions x features (mark levels or mark gain loss)
     - scalgo(``str``) scaling algorithm: sig95, whiten

    """
    chk_exit(*inp_file(path(arr)))
    with open(arr) as fh:
        names = fh.readline().strip().split("\t")
        raw = np.loadtxt(fh, delimiter="\t")
    scaled = scarr(raw, scalgo)
    ofn = arr.replace(".arr", "_%s.arr" % (scalgo,))
    with open(ofn, "wb") as wh:
        wh.write("\t".join(names) + "\n")
        np.savetxt(wh, scaled, delimiter="\t")
    log("saved: %s" % ofn)
    return ofn


@task
def code_sklearn(arr, method=None, init=None, c=None, params=None, transform=True):
    """Non-negative matrix factorization using scikits-learn. 

     - arr(``path``) input array (loci x scaled features) see: ``scale_features``.
     - c(``int``) number of expected histone codes (factorization rank).
     - init(``str``) matrix factorization initialization method.

    """
    chk_exit(*inp_file(path(arr)))
    with open(arr) as fh:
        bam_names = fh.readline().strip().split("\t")
        bam_scaled = np.loadtxt(fh, delimiter="\t")
    kwargs = parse_params(params, {"max_iter":1000})
    nmf = decomposition.NMF(n_components=c, init=init, sparseness='components', **kwargs)
    nmf.fit(bam_scaled)
    ofn_epi = arr.replace(".arr", "_%s-c#%s-i#%s-p#%s.epi" % ("pgnmf", c, init, (params or "")))
    ofn_arr = arr.replace(".arr", "_%s-c#%s-i#%s-p#%s.arr" % ("pgnmf", c, init, (params or "")))
    write_codes(ofn_epi, nmf.components_, bam_names)
    if transform:
        bam_transformed = nmf.transform(bam_scaled)
        write_values(ofn_arr, bam_transformed, c)
    return ofn_epi, ofn_arr

@task
def multi_code_sklearn(arrs, base=None, method=None, init=None, c=None, params=None):
    """Multi-array Non-negative matrix factorization using scikits-learn. 

     - arrs(``path+``) input arrays (loci x scaled features) see: ``scale_features``.
     - base(``str``) common basename for the output.
     - c(``int``) number of expected histone codes (factorization rank).
     - init(``str``) matrix factorization initialization method.

    """
    kwargs = parse_params(params, {"max_iter":1000})
    marks, xs = load_arrs(arrs)
    hs = []
    for x in xs:
        nmf = decomposition.NMF(n_components=c, init=init, sparseness='components', **kwargs)
        nmf.fit(x)
        hs.append(nmf.components_)
    H = np.vstack(hs)
    X = np.vstack(xs)
    W = np.zeros((X.shape[0], len(H)))
    for j in range(0, X.shape[0]):
        W[j, :], _ = nnls(H.T, X[j, :])
    # write codes
    ofnc = base + ("_%s-c#%s-i#%s-p#%s.epi" % ("pgnmf", c, init, (params or "")))
    write_codes(ofnc, H, marks)
    # write 
    ofna = base + ("_%s-c#%s-i#%s-p#%s.arr" % ("pgnmf", c, init, (params or "")))
    write_values(ofna, W, len(arrs)*c)
    return ofnc, ofna

@task
def absolute(bed=None, bams=None, odn=path("absolute_out"), runid=None, shorten=False, par=4, 
             colsca="sig95", method="pgnmf", init="nndsvd", c=None, params=None):
    """Absolute "epigenetic codes" from levels of epigenetic marks in a single experimental 
    condition and in a single set of sites.
 
      - bed(``path``) Genomic regions in the BED6+ file format.
      - bams(``path+``) Sequencing data in coordinated sorted BAM files (requiers BAI index files).
      - odn(``path``) Output directory name.
      - runid(``str``) Run id to prefix all output files.
      - shorten(``bool``) Truncate BAM file names to unambigous strings.
      - par(``int``) Number of parallel processes for BAM extraction.
      - colsca(``str``) Column rescaling method one of: sig95, whiten.
      - method(``str``) currently only pgnmf from scikit-learn is supprted.
      - init(``str``) NMF initialization method (see: scikit-learn documentation for alternatives).
      - c(``int``) number of expected histone codes (factorization rank).
      - params(``str``) Specific parameters for the sklearn PGNMF algorithm (see: scikit-learn for options).

    """
    chk_exit(c is None, "error: c (number of codes) not specified")
    abslvl = extract_absolute(bed, bams, odn, runid, shorten, par)
    abssca = scale_features(abslvl, colsca)
    codes = code_sklearn(abssca, method, init, c, params)
    return codes

@task
def differential(bed=None, abams=None, bbams=None, odn=path("differential_out"), runid=None, shorten=False, step=100, 
                 par=4, pairsca="deseq", colsca="sig95", method="pgnmf", init="nndsvd", c=None, params=None):
    """Differential "epigenetic codes" from "gain-loss" changes in levels of epigenetic marks from two 
    experimental conditions in single set of sites. 
   
      - bed(``path``) Genomic regions in the BED6+ file format.
      - abams(``path+``) Sample A sequencing data in sorted BAM files (requiers BAI index files).
      - bbams(``path+``) Sample B sequencing data in sorted BAM files (requiers BAI index files).
      - step(``int``) Step size (in bp) for coverage calculation with regions.
      - odn(``path``) Output directory name.
      - runid(``str``) Run id to prefix all output files.
      - shorten(``bool``) Truncate BAM file names to unambigous strings.
      - par(``int``) Number of parallel processes for BAM extraction.
      - pairsca(``str``) paired samples scaling method, one of: deseq
      - colsca(``str``) Column rescaling method, one of: sig95, whiten.
      - method(``str``) currently only pgnmf from scikit-learn is supprted.
      - init(``str``) NMF initialization method (see: scikit-learn documentation for alternatives).
      - c(``int``) Number of expected histone codes (factorization rank).
      - params(``str``) Specific parameters for the sklearn PGNMF algorithm (see: scikit-learn for options).

    """
    chk_exit(c is None, "error: c (number of codes) not specified")
    abcnt = extract_diff(bed, abams, bbams, odn, runid, shorten, step, par)
    ablvl = scale_pairs(abcnt, pairsca) # adjust for readdepth
    gllvl = scale_diff(ablvl) # from two sample to gain loss
    glsca = scale_features(gllvl, colsca)
    codes = code_sklearn(glsca, method, init, c, params)
    return codes

@task
def discriminatory(beds=None, bams=None, odn=path("discriminatory_out"), runid=None, shorten=False, par=4, 
             colsca="sig95", init="nndsvd", c=None, params=None):
    """Discriminatory "epigenetic codes" that emphasize epigenetic differences between two (or more [experimental]) 
    sets of sites.

     - beds(``path+``) Two BED6+ files of different sets of genomic sites.
     - bams(``path+``) Sequencing data in sorted BAM files (requiers BAI index files).
     - odn(``path``) Output directory name.
     - runid(``str``) Run id to prefix all output files.
     - shorten(``bool``) truncate BAM file names to unambigous strings.
     - par(``int``) Number of parallel processes for BAM extraction.
     - colsca(``str``) Column rescaling method currently one of: sig95, whiten.
     - init(``str``) NMF initialization method (see: scikit-learn for options).
     - c(``int``) Number of expected histone codes (factorization rank).
     - params(``str``) Specific parameters for the sklearn PGNMF algorithm (see: scikit-learn for options).

    """
    chk_exit(c is None, "error: c (number of codes) not specified")
    arrs = []
    for i, bed in enumerate(beds):
        abslvl = extract_absolute(bed, bams, odn, "%s_%s" % (i, runid), shorten, par)
        abssca = scale_features(abslvl, colsca)
        arrs.append(abssca)
    base = odn / (runid or "discriminatory")
    multepi, multarr = multi_code_sklearn(arrs, base=base, method="pgnmf", init=init, c=c, params=params)
    return multepi, multarr


if __name__ == "__main__":
    DOC = \
"""epicode.py - Discovers "epigenetic codes" within ChIP-seq datasets. The 
goal of epicode is to discover patterns of histone modifications from
aligned sequence data. Epicode looks for combinations (subsets) of marks
that tend to occur in (at least) sub-portions of the data. Alternatively
it identifies combinations of marks that change coordinately i.e. are
"gained" or "lost" frequently at the same time. The algorithm provides three 
modes "absolute", "discriminatory", and "differential". The first two modes 
identify co-occurring marks within one or many sets of genomic losi, 
respectively. The "differential" mode attempts to find patterns of 
coordinated mark changes. In "discriminatory" mode two (or more) genomic 
loci are differentiated based their associated patterns.
"""
    task(DOC)

# @task
# def code_pymf(arr, method=None, init=None, c=None, params=None, transform=True):
#     """(internal) non-negative matrix factorization using scikits-learn
#      - arr(``path``)
#      - c(``int``) number of archetype rows.
#      - c(``int``) number of histone codes.
#      - init(``str``) matrix initialization method.
#      - params(``str``) parameter string [max_iter]
#     """
#     from pymf.aa import AA
#     from pymf.cnmf import CNMF
#     from pymf.chnmf import CHNMF
#     chk_exit(*inp_file(path(arr)))
#     with open(arr) as fh:
#         names = fh.readline().strip().split("\t")
#         scaled = np.loadtxt(fh, delimiter="\t")
#     kwargs = parse_params(params, {"max_iter":1000})
#     data = scaled
#     if method == "aa":
#         model = AA(data, num_bases=c)
#     elif method == "cnmf":
#         model = CNMF(data, num_bases=c)
#     elif method == "chnmf":
#         model = CHNMF(data, num_bases=c)
#     else:
#         raise ValueError("unknow method")
#     model.factorize(niter=kwargs["max_iter"])
#     ofn = arr.replace(".arr", "_%s-c#%s-p#%s.epi" % (method, c, params or ""))
#     write_codes(ofn, model.H, names)
#     if transform:
#         ofn = arr.replace(".arr", "_%s-c#%s-p#%s.arr" % (method, c, params or ""))
#         write_values(ofn, model.W, c)

# @task
# def code_nimfa(arr, method=None, init=None, c=None, params=None):
#     """(internal) non-negative matrix factorization using nimfa
#      - arr(``path``)
#      - method(``str``) NMF factorization method
#      - c(``int``) number of histone codes.
#      - init(``str``) matrix initialization method.
#      - params(``str``) parameter string
#     """
#     from nimfa import mf, mf_run
#     chk_exit(*inp_file(arr))
#     with open(arr) as fh:
#         bam_names = fh.readline().strip().split("\t")
#         bam_scaled = np.loadtxt(fh, delimiter="\t")
#     kwargs = parse_params(params, {"max_iter":1000})
#     decomp = mf(bam_scaled.T, 
#               rank = c, 
#               seed = init, 
#               method = method, 
#               initialize_only = True,
#               **kwargs  
#               )
#     decomp.run()
#     basis = decomp.basis()
#     try:
#         basis = basis.todense()
#     except:
#         pass
#     codes = basis.T.tolist()
#     ofn = arr.replace(".arr", "_%s-c#%s-i#%s-p#%s.epi" % (method, c, init, params or ""))
#     write_codes(ofn, codes, bam_names)
#     if transform:
#         bam_transformed = decomp.fitted()
#         ofn = arr.replace(".arr", "_%s-c#%s-i#%s-p#%s.arr" % (method, c, init, params or ""))
#         write_values(ofn, bam_transformed, c)
#     # elif method in ("archetype",):
#     #     codes = code_archetype(abssca, method, init, c, params)
#     # else:
#     #     codes = code_nimfa(abssca, method, init, c, params)

