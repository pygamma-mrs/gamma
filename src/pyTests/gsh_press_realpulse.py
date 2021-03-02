from __future__ import division

import pygamma as pg
import os

dir = ""

outfile = os.path.join(dir, "gsh_press_realpulses_pytest.txt")
sysfile = os.path.join(dir, "gsh_test.sys")
pulse180file = os.path.join(dir, "bjs180_1.txt")

runname = "test_lines"

h1 = "Press Simulation Test With Real 90,180 Pulses"
h2 = "using input sys file: " + sysfile
h3 = "and input 180 pulse file: " + pulse180file

header = (h1, h2, h3) 


sys = pg.spin_system()

tinit = 0.005           # evolution after 90 before first 180
TE = 0.025              # TE in sec
TE2 = TE/2.0            # TE/2
pulsestep = 0.00001     # 1 msec pulse steps


sys.read(sysfile)

specfreq = sys.Omega()

pulse = pg.row_vector.read_pulse(pulse180file, pg.row_vector.ASCII_MT_DEG)

ptime = pg.row_vector(pulse.size())

for j in range(pulse.size()):
    ptime.put(pg.complex(pulsestep, 0), j)      # pulse steps


pwf = pg.PulWaveform(pulse, ptime, "TestPulse")

pulc = pg.PulComposite(pwf, sys, "1H")

H = pg.Hcs(sys) + pg.HJ(sys)

D = pg.Fm(sys)

Udelay1 = pg.prop(H, tinit)
Udelay2 = pg.prop(H, TE2)
Udelay3 = pg.prop(H, TE2)

ac = pg.acquire1D(pg.gen_op(D), H, 0.001)       # Set up acquisition
ACQ = ac

sigma0 = pg.sigma_eq(sys)                       #Equilibrium density matrix
sigma1 = pg.Iypuls(sys, sigma0, 90.0)           #Apply a 90y pulse
sigma0 = pg.evolve(sigma1, Udelay1)             #Evolve through TINIT

Ureal180  = pulc.GetUsum(-1)                    #Get the propagator for steps of 180

sigma1 = Ureal180.evolve(sigma0)                #Evolve through pulse
sigma0 = pg.evolve(sigma1, Udelay2)             #Evolve through TE/2
sigma1 = Ureal180.evolve(sigma0)                #Evolve through pulse
sigma0 = pg.evolve(sigma1, Udelay3)             #Evolve through TE/2

mx = ACQ.table(sigma0)                          # Transitions table (no lb)

mx.dbwrite(outfile, runname, specfreq, sys.spins(), 0, header)      # Print Table

