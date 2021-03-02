from __future__ import division

import pygamma as pg 

infile = 'gsh_test.sys' 
outfile = "gsh_fid_pytest.txt" 
h1 = "FID Simulation Test"
h2 = "using input sys file: " + infile
outname = "test_lines"
header = (h1, h2)
sys = pg.spin_system() 
sys.read(infile) 
specfreq = sys.Omega() 
H = pg.Hcs(sys) + pg.HJ(sys) 
D = pg.Fm(sys) 
ac = pg.acquire1D(pg.gen_op(D), H, 0.001) 
ACQ = ac 
sigma = pg.sigma_eq(sys) 
sigma0 = pg.Ixpuls(sys, sigma, 90.0) 
mx = ACQ.table(sigma0) 
mx.dbwrite(outfile, outname, specfreq, sys.spins(), 0, header);   # Print Table
