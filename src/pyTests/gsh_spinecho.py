from __future__ import division

import pygamma as pg 

dir = ""

sysfile = dir + "gsh_test.sys"
outfile = dir + "gsh_spinecho_pytest.txt"

runname = "test_lines"

sysfilestr = "using input sys file: " + sysfile

header = ("Spin Echo Simulation Test", sysfilestr)


sys = pg.spin_system()
sys.read(sysfile)

specfreq = sys.Omega()

H = pg.Hcs(sys) + pg.HJ(sys)
D = pg.Fm(sys)

t1 = 0.00833333
t2 = 0.00833333

Udelay1 = pg.prop(H, t1)
Udelay2 = pg.prop(H, t2)

ac = pg.acquire1D(pg.gen_op(D), H, 0.001)

ACQ = ac

sigma0 = pg.sigma_eq(sys)
sigma1 = pg.Ixpuls(sys, sigma0, 90.0)
sigma0 = pg.evolve(sigma1, Udelay1)
sigma1 = pg.Ixpuls(sys, sigma0, 180)
sigma0 = pg.evolve(sigma1, Udelay2)

mx = ACQ.table(sigma0)

mx.dbwrite(outfile, runname, specfreq, sys.spins(), 0, header)
