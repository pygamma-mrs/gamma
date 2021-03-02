from __future__ import division

import pygamma as pg

dir = ""

insysfile = dir + "gsh_test.sys"
inpulse180file = dir + "bjs180_1.txt"

outfile = dir + "gsh_spinecho_realpulses_pytest.txt";

out_name = "test_lines";


h1 = "Spin Echo Simulation Test With Real 180 Pulse"
h2 = "using input sys file: " + insysfile
h3 = "and input 180 pulse file: " + inpulse180file

header = (h1, h2, h3)

sys = pg.spin_system()

t1 = 0.025 
t2 = 0.025
pulsestep = 0.00001

sys.read(insysfile)
specfreq = sys.Omega()

# read_pulse replaces ReadPulse() and is a static member of row_vector
pulse = pg.row_vector.read_pulse(inpulse180file, pg.row_vector.ASCII_MT_DEG)

ptime = pg.row_vector(pulse.size())

#need to use pg.complex() so it can find correct function to call.
for j in range(pulse.size()):
    ptime.put(pg.complex(pulsestep, 0), j)

pwf = pg.PulWaveform(pulse, ptime, "TestPulse")

pulc = pg.PulComposite(pwf, sys, "1H")

H = pg.Hcs(sys) + pg.HJ(sys);
D = pg.Fm(sys);

Udelay1 = pg.prop(H, t1);
Udelay2 = pg.prop(H, t2);

# Neet to effectively typecast D as a gen_op.
ac = pg.acquire1D(pg.gen_op(D), H, 0.001)

ACQ = ac;

sigma0 = pg.sigma_eq(sys)

sigma1 = pg.Iypuls(sys, sigma0, 90.0)   #Apply a 90y pulse

sigma0 = pg.evolve(sigma1, Udelay1)     #Evolve through T1

Ureal180  = pulc.GetUsum(-1)            #Get the propagator for steps of 180

sigma1 = Ureal180.evolve(sigma0)        #Evolve through pulse

sigma0 = pg.evolve(sigma1, Udelay2)     #Evolve through T2

mx = ACQ.table(sigma0)                  #Transitions table (no lb)

mx.dbwrite(outfile, out_name, specfreq, sys.spins(), 0, header)
