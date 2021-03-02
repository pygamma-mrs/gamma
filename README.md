# Welcome to the PyGAMMA/GAMMA Repo

**GAMMA** is a C++ library for the simulation of magnetic resonance spectroscopic (NMR, MRS) experiments. GAMMA is an acronym for a ** *G*eneral *A*pproach to *M*agnetic resonance *M*athematical *A*nalysis **. It provides a simple and intuitive means to construct simulation programs to suit researchers' individual needs. 

**PyGamma** is a Python wrapper around GAMMA that makes almost all of GAMMA's API available via [Python](http://www.python.org/). PyGamma users can skip the C++ compile and link steps and can even run GAMMA commands interactively line-by-line.


## What is GAMMA

With GAMMA it is possible to model a large variety of nuclear magnetic resonance phenomena: This includes the ability to easily represent and simulate spin systems, density operators, Hamiltonians, detection operators, and the effect of various rf pulses - both ideal pulses and more realistic examples characterized by complex waveforms - by using (and extending) the C++ classes provided by GAMMA. Additional built in components of GAMMA are represented as C++ classes, include propagators, superoperators, Liouvillians, transition tables, Bloch equations, spin exchange, quadrupolar, dipolar and Hyperfine interactions, as well as many others.

GAMMA was originally written by Scott A. Smith and Tilo Levante under the guideance of B.H. Meier and R.R. Ernst at the ETH in Zürich. The reference to the original paper is given here:

    "Computer Simulations in Magnetic Resonance. An Object Oriented Programming Approach", S.A. Smith, T.O. Levante, B.H. Meier, and R.R. Ernst, Journal of Magnetic Resonance, 106a, 75-105, (1994). 

It has since been modifed and updated with contributions from a number of individuals.	 

## What is PyGAMMA

PyGAMMA is a library that exposes the GAMMA public intefaces within Python. It provides the functionality of GAMMA within the interactive and scriptable Python environment. 

Using PyGAMMA, users may access GAMMA commands interactively from a Python command prompt, or create complex Python programs. The Object Oriented design of GAMMA is not sacrificed in PyGAMMA. 


For licensing information see the LICENSE file in this directory.

## Platforms

GAMMA is supported on Windows, Linux, OSX, and Solaris.

PyGAMMA is officially supported on Windows, Linux, and OSX.

GAMMA and PyGAMMA are part of the Vespa project (Versatile Simulation, Pulses and Analysis). For information about Vespa, please go to this site: https://githum.com/vespa-mrs/vespa.io

For the most up to date information on how to download, install, use, and contribute to GAMMA or PyGamma, go to this site: 
https://githum.com/gamma-mrs/gamma.io
