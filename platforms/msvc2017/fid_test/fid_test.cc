// This file contains a set of gamma tests
#include "gamma.h"
#include <vector>
#include <string> 

using namespace std;



int main ()
{
  string sysfile = "../src/Tests/gsh_test.sys";

  const char* outfile = "fid_msvc_test.txt";
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

  return 0;
}

