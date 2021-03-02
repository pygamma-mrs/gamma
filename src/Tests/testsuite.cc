// This file contains a set of gamma tests
#include "testsuite.h"
#include <vector>
#include <string> 

using namespace std;

int GammaTest::iso_test()
{
  const char* outfile = "iso_test.txt";
  ofstream fp;

  fp.open(outfile);        // Open output file

  Isotope I1("13C");  
  fp << "\nCarbon 13 Has I Value of " << I1.qn() << "\n";

  Isotope I2("1H");
  fp << "\nHydrogen Has I Value of " << I2.qn() << "\n";

  return -1;
  
  fp.close();             // Close output file  

  return 0;
}


int GammaTest::fid_test(string & sysfile)
{
  const char* outfile = "fid_test.txt";
  string outname = "test_lines";

  vector<string> header;

  header.push_back("FID Simulation Test");
  header.push_back("using input sys file: " + sysfile);
  
  spin_system	sys;
  gen_op		H;
  gen_op		D;
  TTable1D		mx;
  gen_op		sigma0;
  acquire1D		ACQ;
  
  sys.read(sysfile);
  double specfreq = sys.Omega();

  H = Hcs(sys) + HJ(sys);
  D = Fm(sys);

  acquire1D ac(D, H, 0.001);     // Set up acquisition
  ACQ = ac;				  

  sigma0 = sigma_eq(sys);	       // Equilibrium density matrix
  sigma0 = Ixpuls(sys, sigma0, 90.0);    // Apply a 90x pulse

  mx = ACQ.table(sigma0);	       // Transitions table (no lb)
  //mx.dbwrite_old(outfile, outname, -10, 10, specfreq, 0.1, 0, header);   // Print Table
  mx.dbwrite(outfile, outname, specfreq, sys.spins(), 0, header);   // Print Table

  return -1;
}


int GammaTest::spinecho_test(string & sysfile)
{
  const char* outfile = "spinecho_test.txt";
  string outname = "test_lines";

  vector<string> header;

  header.push_back("Spin Echo Simulation Test");
  header.push_back("using input sys file: " + sysfile);

  spin_system	sys;
  gen_op		H;
  gen_op		D;
  TTable1D		mx;
  gen_op		sigma0;
  gen_op        sigma1;
  gen_op        Udelay1;
  gen_op        Udelay2;
  double        t1 = 0.00833333;
  double        t2 = 0.00833333;
  acquire1D		ACQ;
  
  sys.read(sysfile);
  double specfreq = sys.Omega();

  H = Hcs(sys) + HJ(sys);
  D = Fm(sys);

  Udelay1 = prop(H,t1);
  Udelay2 = prop(H,t2);

  acquire1D ac(D, H, 0.001);     // Set up acquisition
  ACQ = ac;				  

  sigma0 = sigma_eq(sys);	       // Equilibrium density matrix
  sigma1 = Ixpuls(sys,sigma0,90.0);    // Apply a 90x pulse
  sigma0 = evolve(sigma1,Udelay1);
  sigma1 = Ixpuls(sys,sigma0,180);
  sigma0 = evolve(sigma1,Udelay2);

  mx = ACQ.table(sigma0);	       // Transitions table (no lb)
  //mx.dbwrite_old(outfile, outname, -10, 10, specfreq, 0.1, 0, header);   // Print Table
  mx.dbwrite(outfile, outname, specfreq, sys.spins(), 0, header);   // Print Table

  return -1;
}


int GammaTest::spinecho_realpulse_test(string & sysfile, string & pulse180file)
{
  const string outfile = "spinecho_realpulses_test.txt";
  string outname = "test_lines";

  vector<string> header;

  header.push_back("Spin Echo Simulation Test With Real 180 Pulse");
  header.push_back("using input sys file: " + sysfile);
  header.push_back("and input 180 pulse file: " + pulse180file);

  spin_system	sys;
  gen_op		H;
  gen_op		D;
  TTable1D		mx;
  gen_op		sigma0;
  gen_op        sigma1;
  gen_op        Udelay1;
  gen_op        Udelay2;
  HSprop        Ureal180;
  double        t1 = 0.025; 
  double        t2 = 0.025; 
  double        pulsestep = 0.00001; // 1 msec pulse steps
  acquire1D		ACQ;
  
  sys.read(sysfile);
  double specfreq = sys.Omega();

  // Read Pulse Files, Initialize Waveform --------------------------

  row_vector pulse = row_vector::read_pulse(pulse180file, row_vector::ASCII_MT_DEG);

  row_vector ptime(pulse.size());

  for(int j=0; j<pulse.size(); j++) 
  {
	ptime.put(pulsestep, j); // pulse steps
  }

  PulWaveform	pwf(pulse, ptime, "TestPulse");

  PulComposite	pulc = PulComposite(pwf, sys, "1H");

  H = Hcs(sys) + HJ(sys);
  D = Fm(sys);

  Udelay1 = prop(H,t1);
  Udelay2 = prop(H,t2);

  acquire1D ac(D, H, 0.001);     // Set up acquisition
  ACQ = ac;				  

  sigma0 = sigma_eq(sys);	       // Equilibrium density matrix
  sigma1 = Iypuls(sys,sigma0,90.0);    // Apply a 90y pulse
  sigma0 = evolve(sigma1,Udelay1);     // Evolve through T1

  Ureal180  = pulc.GetUsum(-1);       // Get the propagator for steps of 180
  sigma1 = Ureal180.evolve(sigma0);    // Evolve through pulse
  sigma0 = evolve(sigma1,Udelay2);     // Evolve through T2

  mx = ACQ.table(sigma0);	       // Transitions table (no lb)
  //mx.dbwrite_old(outfile, outname, -10, 10, specfreq, 0.1, 0, header);   // Print Table
  mx.dbwrite(outfile, outname, specfreq, sys.spins(), 0, header);   // Print Table  

  return -1;
}


int GammaTest::press_realpulses_test(string & sysfile, string & pulse180file)
{
  const string outfile = "press_realpulses_test.txt";
  string outname = "test_lines";

  vector<string> header;

  header.push_back("Press Simulation Test With Real 90,180 Pulses");
  header.push_back("using input sys file: " + sysfile);
  header.push_back("and input 180 pulse file: " + pulse180file);
  
  spin_system	sys;
  gen_op		H;
  gen_op		D;
  TTable1D		mx;
  gen_op		sigma0;
  gen_op        sigma1;
  gen_op        Udelay1;
  gen_op        Udelay2;
  gen_op        Udelay3;
  HSprop        Ureal180;
  double        tinit = 0.005; // evolution after 90 before first 180
  double        TE = 0.025; // TE in sec
  double        TE2 = TE/2.0; // TE/2
  double        pulsestep = 0.00001; // 1 msec pulse steps
  acquire1D		ACQ;

  sys.read(sysfile);
  double specfreq = sys.Omega();
  
  // Read Pulse Files, Initialize Waveform --------------------------

  row_vector	pulse = row_vector::read_pulse(pulse180file, row_vector::ASCII_MT_DEG);
  row_vector	ptime(pulse.size());

  for(int j=0; j<pulse.size(); j++) 
  {
    ptime.put(pulsestep, j); // pulse steps
  }

  PulWaveform	pwf(pulse, ptime, "TestPulse");
  PulComposite	pulc = PulComposite(pwf, sys, "1H");

  H = Hcs(sys) + HJ(sys);
  D = Fm(sys);

  Udelay1 = prop(H,tinit);
  Udelay2 = prop(H,TE2);
  Udelay3 = prop(H,TE2);
  
  acquire1D ac(D, H, 0.001);     // Set up acquisition
  ACQ = ac;				  

  sigma0 = sigma_eq(sys);	       // Equilibrium density matrix
  sigma1 = Iypuls(sys,sigma0,90.0);    // Apply a 90y pulse
  sigma0 = evolve(sigma1,Udelay1);     // Evolve through TINIT

  Ureal180  = pulc.GetUsum(-1);        // Get the propagator for steps of 180

  sigma1 = Ureal180.evolve(sigma0);    // Evolve through pulse
  sigma0 = evolve(sigma1,Udelay2);     // Evolve through TE/2
  sigma1 = Ureal180.evolve(sigma0);    // Evolve through pulse
  sigma0 = evolve(sigma1,Udelay3);     // Evolve through TE/2

  mx = ACQ.table(sigma0);	       // Transitions table (no lb)
  //mx.dbwrite_old(outfile, outname, -10, 10, specfreq, 0.1, 0, header);   // Print Table
  mx.dbwrite(outfile, outname, specfreq, sys.spins(), 0, header);   // Print Table

  return -1;
}


int GammaTest::fid_exchange_test(string & sysfile)
{
  const char* outfile = "fid_exchange_test.txt";
  string outname = "test_lines";

  vector<string> header;

  header.push_back("FID/Exchange  Test");
  header.push_back("using input sys file: " + sysfile);  

  // spin_system sys;
  sys_dynamic sys;				// Declare a system
  sys.read(sysfile);
  double specfreq = sys.Omega();

  TTable1D mx;                          // transition table for output

  gen_op H = Ho(sys);			// Set isotropic Hamiltonian 
  gen_op detect = Fm(sys);		// Set detection operator to F-
  gen_op sigma0 = sigma_eq(sys);		// Set density mx equilibrium
  gen_op sigmap = Iypuls(sys,sigma0, 90.);	// This is 1st 90 pulse

  super_op L = complexi*Hsuper(H);
  L += Kex(sys, H.get_basis());                // Add in exchange

  acquire1D ACQ1(detect,L); 			// Set up the acquisition

  mx = ACQ1.table(sigmap);
  //mx.dbwrite_old(outfile, outname, -10, 10, specfreq, 0.1, 0, header);   // Print Table
  mx.dbwrite(outfile, outname, specfreq, sys.spins(), 0, header);   // Print Table
  
  return -1;
}
