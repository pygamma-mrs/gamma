/*
homo_stat.cc

Simulation of homonuclear dipolar coupled spin system.
Can run up to 10 spins with the CSA tensors. 
*/

#include "gamma.h"
#include <string>
//#include <time.h>
//#include <resource.h>

#define NPROP 100
#define MAXSPINS 12

using namespace std;

int main(int argc, char *argv[])
{
  spin_system ax(2);
  gen_op Ham, det[MAXSPINS], sigma;
  spin_T  Hdip[MAXSPINS][MAXSPINS];
  space_T Adip[MAXSPINS][MAXSPINS], Adip_R[MAXSPINS][MAXSPINS];
  space_T Acsa[MAXSPINS], Acsa_R[MAXSPINS];
  double D[MAXSPINS][MAXSPINS];
  double iso_CSA[MAXSPINS], delta_CSA[MAXSPINS], eta_CSA[MAXSPINS];
  int i,j,count,qu,Fnp;
  string name, names;
  double dw;
  double scale;
  int nspins;
  double alpha,beta,gamma;
  double alpha_D[MAXSPINS][MAXSPINS],beta_D[MAXSPINS][MAXSPINS];
  double gamma_D[MAXSPINS][MAXSPINS];
  double alpha_CSA[MAXSPINS],beta_CSA[MAXSPINS],gamma_CSA[MAXSPINS];
  //struct rusage me;

  int value1[] = {1, 50, 100, 144, 200, 300, 538, 1154, 2000, 5000, 10000, 50000, 100000};
  int value2[] = {1,  7,  27,  11,  29,  37,  55,  107,  297, 1197,  3189, 14857,  38057};
  int value3[] = {1, 11,  41,  53,  79,  61, 229,  271,  479, 1715,  4713,  9027,  27205};

  count = 1;
  query_parameter(argc,argv,count++,"Spin System Name            ? ", names);
//setup for the spin system
  ax.read(names);
  nspins = ax.spins();
  if(nspins < 2 || nspins > 11)
  { cerr << "This program is written for a two to eleven spin system.\n";
    cerr << "Please change your spin system definition in \n";
    cerr << "the file " << names << ".\n";
    cerr << "Aborting ....\n\n";
    exit(-1);
  }
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// this is -delta/4Pi = omegaD/2Pi = +mu/4Pi gamma1*gamma2*hbar/2Pi /r^3
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  for(i=0;i<nspins-1;++i)
  { for(j=i+1;j<nspins;++j)
    { query_parameter(argc,argv,count++,"Dipolar Coupling Constant   ? ", D[i][j]);
      query_parameter(argc,argv,count++,"Euler angle alpha           ? ", alpha_D[i][j]);
      query_parameter(argc,argv,count++,"Euler angle beta            ? ", beta_D[i][j]);
      query_parameter(argc,argv,count++,"Euler angle gamma           ? ", gamma_D[i][j]);
    }
  }
  for(i=0;i<nspins;++i)
  { query_parameter(argc,argv,count++,"isotropic chemical shift    ? ", iso_CSA[i]);
    query_parameter(argc,argv,count++,"CSA tensor (delta CSA)      ? ", delta_CSA[i]);
    query_parameter(argc,argv,count++,"CSA tensor (eta CSA)        ? ", eta_CSA[i]);
    query_parameter(argc,argv,count++,"Euler angle alpha           ? ", alpha_CSA[i]);
    query_parameter(argc,argv,count++,"Euler angle beta            ? ", beta_CSA[i]);
    query_parameter(argc,argv,count++,"Euler angle gamma           ? ", gamma_CSA[i]);
  }
  query_parameter(argc,argv,count++,"Powder Quality (cheng)      ? ", qu);
  query_parameter(argc,argv,count++,"Number of Data Points       ? ", Fnp);
  query_parameter(argc,argv,count++,"Dwell time                  ? ", dw);
  query_parameter(argc,argv,count++,"Output Filename             ? ", name);


  cout << "\n\nSimulation of isotropic chemical shift by dipolar coupling\n";
  cout << "==========================================================\n\n";
  //cout << "Program version: " << __FILE__ << " compiled at " << __DATE__ ", "
  //     << __TIME__ << "\n\n";
  cout << "Parameters:\n";
  cout << "size of spin system            : " << nspins << " spins\n";
  for(i=0;i<nspins-1;++i)
  { for(j=i+1;j<nspins;++j)
    { cout << "dipolar coupling constant (" << i << "," << j << ") : " << 
              D[i][j] << " Hz\n";
      cout << "relativ orientation of D tensor: (" << alpha_D[i][j] << "," <<
               beta_D[i][j] << "," << gamma_D[i][j] << ")\n";
    }
  }
  for(i=0;i<nspins;++i)
  { cout << "isotropic chem shift(" << i << ") : " << 
            iso_CSA[i] << " Hz\n";
    cout << "delta of CSA tensor (" << i << ") : " << 
            delta_CSA[i] << " Hz\n";
    cout << "eta of CSA tensor (" << i << ") : " << 
            eta_CSA[i] << "\n";
    cout << "relativ orientation of CSA tensor: (" << alpha_CSA[i] << "," <<
             beta_CSA[i] << "," << gamma_CSA[i] << ")\n";
  }

	char dws[100];
#ifdef _MSC_VER
	sprintf_s ( dws, "%f", dw);
#else
	sprintf ( dws, "%f", dw);
#endif

  cout << "Powder Quality Number:         " << qu << "  (" << value1[qu] <<
          " orientations)\n";
  cout << "Number of data points:         " << Fnp << " points\n";
  cout << "dwell time:                    " << dws << "s\n";
  cout << "Output filename:               " << name << "\n";
  cout << "\n";

  block_1D datax(Fnp);
  block_2D data(nspins,Fnp);

  for(i=0;i<nspins;++i)
  { for(j=0;j<Fnp;++j)
    { data(i,j)=0;
    }
  }

//setup for the hamiltonian
  for(i=0;i<nspins-1;++i)
  { for(j=i+1;j<nspins;++j)
    { Hdip[i][j] = T_D(ax,i,j);
    }
  }

//setup for the space tensor
  matrix help(3,3,0);
  for(i=0;i<nspins-1;++i)
  { for(j=i+1;j<nspins;++j)
    { help.put_h(-1.0,0,0);
      help.put_h(-1.0,1,1);
      help.put_h( 2.0,2,2);
      help   = - (complex) D[i][j] * help;
      Adip[i][j] = A2(help);
      Adip[i][j] = Adip[i][j].rotate(alpha_D[i][j],beta_D[i][j],gamma_D[i][j]);
    }
  }
  for(i=0;i<nspins;++i)
  { help.put_h(-1.0/2.0*(1.0+eta_CSA[i]),0,0);
    help.put_h(-1.0/2.0*(1.0-eta_CSA[i]),1,1);
    help.put_h( 1.0,2,2);
    help = (complex) delta_CSA[i] * help;
    Acsa[i] = A2(help);
    Acsa[i] = Acsa[i].rotate(alpha_CSA[i],beta_CSA[i],gamma_CSA[i]);
  }

  string name1 = name + ".mat";
  string name2 = name;

//this is the detection operator
  for(i=0;i<nspins;++i)
    det[i]  = Ip(ax,i);

//here starts the powder loop
//reference JCP 59 (8), 3992 (1973).

  for(count=1; count<=value1[qu]; ++count)
  { 
		beta  = 180.0 * count/value1[qu];
    alpha = 360.0 * ((value2[qu]*count) % value1[qu])/value1[qu];
    gamma = 360.0 * ((value3[qu]*count) % value1[qu])/value1[qu];
    if(count % 1000 == 1)
    { 
			// DCT, 2009-11-26. 
			// getrusage() is a Unix/Linux only command, and does not work on windows.
			// Also, the time elapsed would be very machine specific so should
			// not be part of this testing library (other libraries, maybe!).
			// getrusage(0, & me);
      cout << count << "\tbeta = " << beta << "\talpha = "
           << alpha << "\tgamma = " << gamma << "\n";
      // << ",\ttime used: " << me.ru_utime.tv_sec << " seconds\n";
      cout.flush();
    }
    scale = sin(beta/180.0*PI);

    sigma = Fx(ax);
  
//now we rotate the space tensor
    for(i=0;i<nspins-1;++i)
    { for(j=i+1;j<nspins;++j)
      { Adip_R[i][j] = Adip[i][j].rotate(alpha,beta,gamma);
      }
    }
    for(i=0;i<nspins;++i)
    { Acsa_R[i] = Acsa[i].rotate(alpha,beta,gamma);
    }

    Ham=gen_op();
    for(i=0;i<nspins-1;++i)
    { for(j=i+1;j<nspins;++j)
      { Ham += Adip_R[i][j].component(2,0)*Hdip[i][j].component(2,0);
      }
    }
    for(i=0;i<nspins;++i)
    {  Ham += Acsa_R[i].component(2,0) * 1/sqrt(6.0)*2*Iz(ax,i)+iso_CSA[i]*Iz(ax,i);
    }
    for(i=0;i<nspins;++i)
    { FID(sigma,det[i],Ham,dw,Fnp,datax);
      for(j=0;j<Fnp;++j)
      { data(i,j) += datax(j)*scale;
      }
    }

  } // end of powder loop

  //cout << "\nEntering matlab with parameters:  " << name1 << ", " << name2 << ", " << "data" << ", 1\n";

  MATLAB(name1, name2, data, 1);
  exit(0);
}
