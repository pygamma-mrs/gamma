from __future__ import division

import pygamma as pg

infile = "gsh_test.sys"
outfile = "gsh_fid_exchange_pytest.txt"

s1 = "FID/Exchange  Test"
s2 = "using input sys file: " + infile

runname = "test_lines"

header = (s1, s2)

sys = pg.sys_dynamic()

sys.read(infile)

specfreq = sys.Omega()

mx = pg.TTable1D()

H = pg.Ho(sys)

detect = pg.Fm(sys)

sigma0 = pg.sigma_eq(sys)

sigmap = pg.Iypuls(sys, sigma0, 90.)

L = pg.Hsuper(H)
L *= pg.complex(0,1)

L += pg.Kex(sys, H.get_basis());

ACQ1 = pg.acquire1D(pg.gen_op(detect), L)

mx = ACQ1.table(sigmap);

#mx.dbwrite_old(outfile, "test_lines", -10, 10, specfreq, .1, 0, header)
mx.dbwrite(outfile, runname, specfreq, sys.spins(), 0, header)

