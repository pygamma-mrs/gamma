/*

decouple5R.cc

decoupling with a continous +-phi phase modulated sequence
modulation frequency must be n times the spinning speed
static simulation
additional spin diffusion simulated by Liouville space
relaxation matrix
use MAS small step integration

arbitrary orientation of the CAS tensor

*/

#include "gamma.h"
//#include <sys/time.h>
//#include <sys/resource.h>

using namespace std;

#define NPROP 100

//spin 0 is detected, spin 1 is irradiated

int main(int argc, char *argv[])

{
  spin_system ax(2);
  gen_op temp, Ham, H0, H[5], detect, sigma, sigma1;
  super_op UT, Iop, L, L0, U[NPROP];
  space_T Adip, Acsa, Acsa_R, Adip_R;
  double D, delta_CSA, etha_CSA;
  double k, gamB1;
  int i,j,Fnp,count,qu,steps;
  string name, names;
  const double thetam=54.73561032;
  double mas_freq, time_gamB1, mod_freq;
  double ltime, time, time_prop, scale;
  int index, nprop, ntimes;
  double alpha_CSA,beta_CSA,gamma_CSA;
  double alpha,beta,gamma;
  double phi;
  //struct rusage me;
  //char hostname[200];
  matrix exchange(16,16,0);


  int value1[] = {2, 50, 100, 144, 200, 300, 538, 1154};
  int value2[] = {1,  7,  27,  11,  29,  37,  55,  107};
  int value3[] = {1, 11,  41,  53,  79,  61, 229,  271};

  count = 1;
  query_parameter(argc,argv,count++,"Spin System Name            ? ", names);
  query_parameter(argc,argv,count++,"Dipolar Coupling Constant   ? ", D);
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// this is -delta/4Pi = omegaD/2Pi = +mu/4Pi gamma1*gamma2*hbar/2Pi /r^3
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  query_parameter(argc,argv,count++,"CSA tensor (delta CSA)      ? ", delta_CSA);
  query_parameter(argc,argv,count++,"CSA tensor (etha CSA)       ? ", etha_CSA);
  query_parameter(argc,argv,count++,"Euler angle alpha           ? ", alpha_CSA);
  query_parameter(argc,argv,count++,"Euler angle beta            ? ", beta_CSA);
  query_parameter(argc,argv,count++,"Euler angle gamma           ? ", gamma_CSA);
  query_parameter(argc,argv,count++,"Field Strength gamB1        ? ", gamB1);
  query_parameter(argc,argv,count++,"Phase angle +-phi           ? ", phi);
  query_parameter(argc,argv,count++,"modulation frequency        ? ", ntimes);
  query_parameter(argc,argv,count++,"Rate constant k (1/s)       ? ", k);
  query_parameter(argc,argv,count++,"Powder Quality (cheng)      ? ", qu);
  query_parameter(argc,argv,count++,"Number of time steps        ? ", steps);
  query_parameter(argc,argv,count++,"spinning speed              ? ", mas_freq);
  query_parameter(argc,argv,count++,"Number of sampling per rotor? ", nprop);
  query_parameter(argc,argv,count++,"Number of Data Points       ? ", Fnp);
  query_parameter(argc,argv,count++,"Output Filename             ? ", name);

  time_prop  = (1.0/mas_freq)/nprop;
  time       = (1.0/mas_freq)/steps;
  time_gamB1 = 1.0/gamB1;
  mod_freq = ntimes*mas_freq;
  //(void) gethostname(hostname,199);

  cout << "\n\nSimulation of isotropic chemical shift by dipolar coupling\n";
  cout << "==========================================================\n\n";
  //  cout << "Program version: " << __FILE__ << " compiled at " << __DATE__ ", "
  //       << __TIME__ << "\n";
  //  cout << "running on machine " << hostname << "\n\n";
  cout << "Parameters:\n";
  cout << "rotation angle thetam    : " << thetam << " Degree\n";
  cout << "dipolar coupling constant: " << D << " Hz\n";
  cout << "CSA tensor (delta)       : " << delta_CSA << " Hz\n";
  cout << "CSA tensor (etha)        : " << etha_CSA << "\n";
  cout << "relativ orientation of CSA tensor: (" << alpha_CSA << "," <<
          beta_CSA << "," << gamma_CSA << ")\n";
  cout << "rf-Field Strength:         " << gamB1 << " Hz\n";
  cout << "Rate constant (k):         " << k << " 1/s\n";
  cout << "# of sampling points:      " << nprop << "\n";
  cout << "Powder Quality Number:     " << qu << " (" << value1[qu] << " points)\n";
  cout << "Number of data points:     " << Fnp << " points\n";
  cout << "basic rf pulse frequency:  " << gamB1 << " Hz\n";

	char tprop[100];
#ifdef _MSC_VER
	sprintf_s ( tprop, "%.10f", time_prop);
#else
	sprintf ( tprop, "%.10f", time_prop);
#endif

  cout << "length of one prop. (dw):  " << tprop << " s\n";
  cout << "phase angle (phi):         " << phi << " degree\n";
  cout << "Modulation frequency:      " << mod_freq << " Hz\n";
  cout << "MAS frequency:             " << mas_freq << " Hz\n";

	char timex[100];
#ifdef _MSC_VER
	sprintf_s ( timex, "%.12f", time);
#else
	sprintf ( timex, "%.12f", time);
#endif

  cout << "time increments:           " << timex << "s\n";
  cout << "Output filename:           " << name << "\n";
  cout << "\n";


  block_1D data(Fnp);
  for(i=0;i<Fnp;++i)
    data(i) = 0;


//setup for the spin system
  ax.read(names);

//these are the constant parts of H and L
  for(i=0;i<16;++i)
    exchange.put_h(1,i,i);
  Iop = super_op(exchange);

  temp = Iz(ax,1);
  L0  = k*d_commutator(temp);
  temp = Ix(ax,1);
  L0 += k*d_commutator(temp);
  temp = Iy(ax,1);
  L0 += k*d_commutator(temp);

  detect = Im(ax,0);
  sigma  = Ix(ax,0);
  H0     = ax.shift(0) * Iz(ax,0) + ax.shift(1) * Iz(ax,1);
  phi = phi/180.0*PI;

//setup for the space tensor
  matrix help(3,3,0);
  help.put_h(-1.0,0,0);
  help.put_h(-1.0,1,1);
  help.put_h( 2.0,2,2);
  help = - (complex) D * help;
  Adip = A2(help);
  help.put_h(-1.0/2.0*(1.0+etha_CSA),0,0);
  help.put_h(-1.0/2.0*(1.0-etha_CSA),1,1);
  help.put_h( 1.0,2,2);
  help = (complex) delta_CSA * help;
  Acsa = A2(help);
//rotate CSA tensor into the PA system of the dipolar tensor
  Acsa = Acsa.rotate(alpha_CSA,beta_CSA,gamma_CSA);

  string name1 = name+".mat";
  string name2 = name;

//here starts the powder loop
//reference JCP 59 (8), 3992 (1973).

  for(count=1; count<value1[qu]; ++count)
  { beta  = 180.0 * count/value1[qu];
    alpha = 360.0 * ((value2[qu]*count) % value1[qu])/value1[qu];
    gamma = 360.0 * ((value3[qu]*count) % value1[qu])/value1[qu];
    if(count % 10 == 1)
    { //getrusage(0, & me);
      cout << count << "\tbeta = " << beta << "\talpha = "
           << alpha << "\tgamma = " << gamma 
           << "\n";
      cout.flush();
    }
    scale = sin(beta/180.0*PI);

//now we rotate the space tensor
    Adip_R = Adip.rotate(alpha,beta,gamma);
    Acsa_R = Acsa.rotate(alpha,beta,gamma);

//now we can combine the hamiltonians for the different side diagonals

//zero all components
    for(i=0;i<5;++i)
      H[i] = gen_op();

//this is the dipolar part
    for(j=-2;j<=2;++j)
    { H[j+2] += Adip_R.component(2,j) * d2(j,0,thetam) * 1/sqrt(6.0)*2*Iz(ax,0)*Iz(ax,1);
      H[j+2] += Acsa_R.component(2,j) * d2(j,0,thetam) * 1/sqrt(6.0)*2*Iz(ax,1);
    }

    for(i=0;i<nprop;++i)
      U[i] = Iop;

//now we calculate the propagator over one cycle of the MAS
    for(ltime=0.5*time;ltime<=steps*time;ltime += time)
    { Ham = H0;
      Ham += gamB1 * (Ix(ax,1)*cos(phi*sin(mod_freq*2*PI*ltime)) +
                      Iy(ax,1)*sin(phi*sin(mod_freq*2*PI*ltime)));
      for(i=-2;i<=2;++i)
        Ham += exp(complex(0,i*2.0*PI*ltime*mas_freq)) * H[i+2];
      index = int(ltime/time_prop);
      Ham.set_DBR();
      L  = complex(0,2.0*PI)*commutator(Ham);
      L  = L+L0;
      UT = exp(L,-time);
      UT.set_HBR();
      U[index] = UT*U[index];
    }
    for(i=0;i<nprop;++i)
    { U[i].set_HBR();
    }
    detect.set_DBR();
    sigma1=sigma;
    sigma1.set_DBR();
    for(i=0;i<Fnp;++i)
    { data(i) += trace(detect,sigma1)*scale;
      sigma1 = U[i%nprop]*sigma1;
    }
  } // end of powder loop
  MATLAB(name1,name2,data,1);

  return 0;
}
