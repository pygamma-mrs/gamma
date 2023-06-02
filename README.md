# Welcome to the PyGAMMA/GAMMA Repo

**GAMMA** is a C++ library for the simulation of magnetic resonance spectroscopic (NMR, MRS) experiments. GAMMA is an acronym for a ** *G*eneral *A*pproach to *M*agnetic resonance *M*athematical *A*nalysis **. It provides a simple and intuitive means to construct simulation programs to suit researchers' individual needs. 

**PyGamma** is a Python wrapper around GAMMA that makes almost all of GAMMA's API available via [Python](http://www.python.org/). PyGamma users can skip the C++ compile and link steps and can even run GAMMA commands interactively line-by-line.


## What is GAMMA

With GAMMA it is possible to model a large variety of nuclear magnetic resonance phenomena: This includes the ability to easily represent and simulate spin systems, density operators, Hamiltonians, detection operators, and the effect of various rf pulses - both ideal pulses and more realistic examples characterized by complex waveforms - by using (and extending) the C++ classes provided by GAMMA. Additional built in components of GAMMA are represented as C++ classes, include propagators, superoperators, Liouvillians, transition tables, Bloch equations, spin exchange, quadrupolar, dipolar and Hyperfine interactions, as well as many others.

## What is PyGAMMA

PyGAMMA is a library that exposes the GAMMA public intefaces within Python. It provides the functionality of GAMMA within the interactive and scriptable Python environment. 

Using PyGAMMA, users may access GAMMA commands interactively from a Python command prompt, or create complex Python programs. The Object Oriented design of GAMMA is not sacrificed in PyGAMMA. 

## Citation

If you publish material that makes use of GAMMA, please cite the following publication:

Smith SA, Levante TO, Meier BH, and Ernst RR. Computer Simulations in Magnetic Resonance. An Object Oriented Programming Approach. J Magn. Res. 1994; 106a:75-105.

If you publish material that makes use of PyGAMMA, please cite the following publication:

Soher B, Semanchuk P, Todd D, Ji X, Deelchand D, Joers J, Oz G and Young K. Vespa: Integrated applications for RF pulse design, spectral simulation and MRS data analysis. Magn Reson Med. 2023;1-16. epub doi: 10.1002/mrm.29686

GAMMA was originally written by Scott A. Smith and Tilo Levante under the guideance of B.H. Meier and R.R. Ernst at the ETH in Zürich. It has since been modifed and updated with contributions from a number of individuals. The expansion of PyGAMMA using SWIG was accomplished as part of the [Vespa](https://github.com/vespa-mrs/vespa) project.	 

## Licensing

BSD, specifically a "three-clause" BSD license

## Platforms

GAMMA is officially supported on Windows, Linux, OSX, and Solaris.

PyGAMMA is officially supported on Windows, Linux, and OSX.

GAMMA and PyGAMMA are part of and supported by the Vespa project (Versatile Simulation, Pulses and Analysis). For information about Vespa, please go to this site: https://github.com/vespa-mrs/vespa.io

For the most up to date information on how to download, install, use, and contribute to GAMMA or PyGamma, go to this site: 
https://pygamma-mrs.github.io/gamma.io
