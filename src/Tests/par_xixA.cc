/*

cp_xixA.cc

Double XiX sequence for CP
Can run up to 10 spins without the CSA tensors. 
Brute force integration of the MAS rotation (one cycle)

*/

#include "gamma.h"
//#include <sys/time.h>
//#include <sys/resource.h>

#define NPROP 8000
#define MAXSPINS 10

/* lcm.c */
/* GCD (Gratest Common Devisor) and LCM (Least/Lowest Common Multiple) */

int lcm(int a, int b)
{ int c;

  a = abs(a);
  b = abs(b);
  c = a*b;
  while (a != b)
  { if (a < b)
      b -= a;
    else
      a -= b;
  }
  return(c/a) ;
}

using namespace std;

int main(int argc, char *argv[])

{
  spin_system ax;
  gen_op Hrf,Ham, U[NPROP], Uevolve[NPROP], H[5], sigma1, sigma, detect[MAXSPINS];
  spin_T  Hdip[MAXSPINS][MAXSPINS];
  space_T Adip[MAXSPINS][MAXSPINS], Adip_R[MAXSPINS][MAXSPINS];
  space_T Acsa[MAXSPINS], Acsa_R[MAXSPINS];
  double D[MAXSPINS][MAXSPINS];
  double J[MAXSPINS][MAXSPINS];
  double iso_CSA[MAXSPINS];
  double eta_CSA[MAXSPINS];
  double delta_CSA[MAXSPINS];
  int i,j,k,k1,k2,Fnp,count,qu,steps;
  int n_pwI, n_pwS, n_evolve, n_cycle, n_prop;
  string name, names;
  const double thetam=54.73561032;
  double gamB1I, gamB1S, mas_freq;
  double ltime, time, scale;
  double time_cycle, time_prop, time_mas, time_pwI, time_pwS;
  int nstart, nspins, index;
  double alpha,beta,gamma;
  double angleI, angleS;
  double alpha_CSA[MAXSPINS],beta_CSA[MAXSPINS];
  double gamma_CSA[MAXSPINS];
  double alpha_D[MAXSPINS][MAXSPINS],beta_D[MAXSPINS][MAXSPINS];
  double gamma_D[MAXSPINS][MAXSPINS];
  //struct rusage me;
  //char hostname[200];

  int value1[] = {1, 50, 100, 144, 200, 300, 538, 1154};
  int value2[] = {1,  7,  27,  11,  29,  37,  55,  107};
  int value3[] = {1, 11,  41,  53,  79,  61, 229,  271};

  count = 1;
  query_parameter(argc,argv,count++,"Spin System Name            ? ", names);
//setup for the spin system
  ax.read(names);
  nspins = ax.spins();
  if(nspins < 2 || nspins >= 10)
  { cerr << "This program is written for a two to nine spin system.\n";
    cerr << "Please change your spin system definition in \n";
    cerr << "the file " << names << ".\n";
    cerr << "Aborting ....\n\n";
    exit(-1);
  }
  for(i=0;i<nspins-1;++i)
  { for(j=i+1;j<nspins;++j)
    { query_parameter(argc,argv,count++,"Scalar Coupling Constant    ? ", J[i][j]);
      query_parameter(argc,argv,count++,"Dipolar Coupling Constant   ? ", D[i][j]);
      query_parameter(argc,argv,count++,"Euler angle alpha           ? ", alpha_D[i][j]);
      query_parameter(argc,argv,count++,"Euler angle beta            ? ", beta_D[i][j]);
      query_parameter(argc,argv,count++,"Euler angle gamma           ? ", gamma_D[i][j]);
    }
  }
  for(i=0;i<nspins;++i)
  { query_parameter(argc,argv,count++,"Isotropic chemical shift    ? ", iso_CSA[i]);
    query_parameter(argc,argv,count++,"anisotropy chemical shift   ? ", delta_CSA[i]);
    query_parameter(argc,argv,count++,"asymmetry chemical shift    ? ", eta_CSA[i]);
    query_parameter(argc,argv,count++,"Euler angle alpha           ? ", alpha_CSA[i]);
    query_parameter(argc,argv,count++,"Euler angle beta            ? ", beta_CSA[i]);
    query_parameter(argc,argv,count++,"Euler angle gamma           ? ", gamma_CSA[i]);
  }
  query_parameter(argc,argv,count++,"Powder Quality (cheng)      ? ", qu);
  query_parameter(argc,argv,count++,"Number of time steps        ? ", steps);
  query_parameter(argc,argv,count++,"rf field amplitude I spins  ? ", gamB1I);
  query_parameter(argc,argv,count++,"rf field amplitude S spins  ? ", gamB1S);
  query_parameter(argc,argv,count++,"MAS cycle time              ? ", time_mas);
  query_parameter(argc,argv,count++,"Pulse length XiX pulse I    ? ", time_pwI);
  query_parameter(argc,argv,count++,"Pulse length XiX pulse S    ? ", time_pwS);
  query_parameter(argc,argv,count++,"Number of rotor periods     ? ", n_cycle);
  query_parameter(argc,argv,count++,"Number of sampling per rotor? ", n_prop);
  query_parameter(argc,argv,count++,"Number of data points       ? ", Fnp);
  query_parameter(argc,argv,count++,"Initial density operator?   ? ", nstart);
  query_parameter(argc,argv,count++,"Output Filename             ? ", name);
 
  if(nstart >= nspins)
  { cerr << "The starting density operator has a spin number which is\n";
    cerr << "larger than the number of spins.\n";
    cerr << "Aborting ....\n\n";
    exit(-1);
  }

  time_prop  = time_mas/n_prop;
  time_cycle = n_cycle*time_mas;
  time       = time_mas/steps;
  mas_freq   = 1.0/time_mas;
  n_pwI      = int(time_pwI/time_prop+0.5);
  n_pwS      = int(time_pwS/time_prop+0.5);
  n_evolve   = lcm(lcm(2*n_pwI,n_prop),2*n_pwS)/n_prop;
  if(n_evolve>Fnp)
    n_evolve=Fnp;
  angleI     = n_pwI*time_prop*gamB1I*360;
  angleS     = n_pwS*time_prop*gamB1S*360;

  //(void) gethostname(hostname,199);

  cout << "\n\nSimulation of CP Under Double XiX Irradiation             \n";
  cout << "==========================================================\n\n";
  //cout << "Program version: " << __FILE__ << " compiled at " << __DATE__ ", "
  //     << __TIME__ << "\n\n";
  //cout << "running on machine " << hostname << "\n\n";
  cout << "Parameters:\n";
  cout << "rotation angle thetam          : " << thetam << " Degree\n";
  cout << "size of spin system            : " << nspins << " spins\n";
  cout << "initial density operator       : Iz(ax," << nstart << ")\n"; 
  for(i=0;i<nspins-1;++i)
  { for(j=i+1;j<nspins;++j)
    { cout << "scalar  coupling constant (" << i << "," << j << ") : " << 
              J[i][j] << " Hz\n";
      cout << "dipolar coupling constant (" << i << "," << j << ") : " << 
              D[i][j] << " Hz\n";
      cout << "relativ orientation of D tensor: (" << alpha_D[i][j] << "," <<
               beta_D[i][j] << "," << gamma_D[i][j] << ")\n";
    }
  }
  for(i=0;i<nspins;++i)
  { cout << "isotropic chemical shift (" << i << ") : " << iso_CSA[i] << " Hz\n";
    cout << "anisotropy chemical shift (" << i << "): " << delta_CSA[i] << " Hz\n";
    cout << "asymmetry chemical shift (" << i << ") : " << eta_CSA[i] << "\n";
    cout << "relativ orientation of CSA tensor: (" << alpha_CSA[i] << "," <<
             beta_CSA[i] << "," << gamma_CSA[i] << ")\n";
  }
  cout << "# of sampling points:          " << n_prop << "\n";
  cout << "# of rotor periods:            " << n_cycle << "\n";
  cout << "rf field amplitude on I spins: " << gamB1I << " Hz\n";
  cout << "rf field amplitude on S spins: " << gamB1S << " Hz\n";
  cout << "Powder Quality Number:         " << qu << "  (" << value1[qu] <<
          " orientations)\n";
  cout << "Number of data points:         " << Fnp << " points\n";
  cout << "MAS frequency:                 " << mas_freq << " Hz\n";
  
  char timex[100];
#ifdef _MSC_VER
	sprintf_s ( timex, "%.6f", time_pwI);
#else
	sprintf ( timex, "%.6f", time_pwI);
#endif
  
  cout << "pulse length (XiX I):          " << timex << " s\n";
  cout << "pulse length (XiX I):          " << n_pwI << " propagators\n";
  
#ifdef _MSC_VER
	sprintf_s ( timex, "%.6f", time_pwS);
#else
	sprintf ( timex, "%.6f", time_pwS);
#endif
  
  cout << "pulse length (XiX S):          " << timex << " s\n";
  cout << "pulse length (XiX S):          " << n_pwS << " propagators\n";
  cout << "# of propagators/rotor cycle:  " << n_prop << "\n";
  
#ifdef _MSC_VER
	sprintf_s ( timex, "%.8f", time_prop);
#else
	sprintf ( timex, "%.8f", time_prop);
#endif  
  
  cout << "length of one propagator:      " << timex << " s\n";
  cout << "# of rotor cycles:             " << n_cycle << "\n";
  cout << "# of evolution steps:          " << n_evolve << "\n";
  cout << "flip angle per pulse (I):      " << angleI <<" degree\n";
  cout << "flip angle per pulse (S):      " << angleS <<" degree\n";
  
#ifdef _MSC_VER
	sprintf_s ( timex, "%.8f", time);
#else
	sprintf ( timex, "%.8f", time);
#endif  

  cout << "time increments:               " << timex << "s\n";
  cout << "Output filename:               " << name << "\n";
  cout << "\n";
  cout.flush();

  block_2D data(nspins,Fnp);
  for(k1=0;k1<nspins;++k1)
  { for(k2=0;k2<Fnp;++k2)
    { data(k1,k2) = 0;
    }
  }

//setup for the hamiltonian
  for(i=0;i<nspins-1;++i)
  { for(j=i+1;j<nspins;++j)
    { Hdip[i][j] = T_D(ax,i,j);
    }
  }

  Hrf = gamB1I * Fx(ax,"1H") + gamB1S * Fx(ax,"13C");

//setup for the space tensor
  matrix help(3,3,0);
  for(i=0;i<nspins-1;++i)
  { for(j=i+1;j<nspins;++j)
    { help.put_h(-1.0/2.0,0,0);
      help.put_h(-1.0/2.0,1,1);
      help.put_h( 1.0,2,2);
      help   = - (complex) D[i][j] * help;
      Adip[i][j] = A2(help);
      Adip[i][j] = Adip[i][j].rotate(alpha_D[i][j],beta_D[i][j],gamma_D[i][j]);
    }
  }

  for(i=0;i<nspins;++i)
  { help.put(-1.0/2.0*(1.0+eta_CSA[i]),0,0);
    help.put(-1.0/2.0*(1.0-eta_CSA[i]),1,1);
    help.put( 1.0,2,2);
    help = (complex) delta_CSA[i] * help;
    Acsa[i] = A2(help);
    Acsa[i] = Acsa[i].rotate(alpha_CSA[i],beta_CSA[i],gamma_CSA[i]);
  }

  string name1 = name+".mat";
  string name2 = name;

//here starts the powder loop
//reference JCP 59 (8), 3992 (1973).

  for(count=1; count<=value1[qu]; ++count)
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

    sigma=  Ix(ax,nstart);
    for(k=0;k<nspins;++k)
    { detect[k] = Ix(ax,k);
    }
  
//now we rotate the space tensor
    for(i=0;i<nspins;++i)
    { for(j=i+1;j<nspins;++j)
      { Adip_R[i][j] = Adip[i][j].rotate(alpha,beta,gamma);
      }
    }
    for(i=0;i<nspins;++i)
    { Acsa_R[i] = Acsa[i].rotate(alpha,beta,gamma);
    }

//zero all components
    for(i=0;i<5;++i)
      H[i] = gen_op();

//this is the dipolar part
    for(i=0;i<nspins;++i)
    { H[2] += iso_CSA[i]*Iz(ax,i);
    }
    for(k=-2;k<=2;++k)
    { for(i=0;i<nspins-1;++i)
      { H[k+2] += Acsa_R[i].component(2,k) * d2(k,0,thetam)*2.0/sqrt(6.0)*Iz(ax,i);
        for(j=i+1;j<nspins;++j)
        { if(ax.isotope(i) != ax.isotope(j))
          { H[k+2] += Adip_R[i][j].component(2,k) * d2(k,0,thetam) * 1.0/sqrt(6.0)*2*Iz(ax,i)*Iz(ax,j);
	    if(k==0)
	      H[k+2] += J[i][j] * Iz(ax,i)*Iz(ax,j);
          }
          else
          { H[k+2] += Adip_R[i][j].component(2,k) * d2(k,0,thetam) * Hdip[i][j].component(2,0);
	    if(k==0)
	      H[k+2] += J[i][j] * (Iz(ax,i)*Iz(ax,j)+Ix(ax,i)*Ix(ax,j)+Iy(ax,i)*Iy(ax,j));
          }
        }
      }
    }

    for(i=0;i<n_prop;++i)
      U[i] = Ie(ax,0);

//now we calculate the propagator over one cycle of the MAS
    for(ltime=0.5*time;ltime<=steps*time;ltime += time)
    { Ham = Hrf;
      for(i=-2;i<=2;++i)
        Ham += exp(complex(0,i*2.0*PI*ltime*mas_freq)) * H[i+2];
      index = int(ltime/time_prop);
      U[index] &= prop(Ham,time);
    }

    for(i=0;i<n_evolve;++i)
      Uevolve[i] = Ie(ax,0);

    for(i=0;i<n_prop*n_evolve;++i)
      Uevolve[i/n_prop] &= evolve(U[i%n_prop],Rz(ax,"1H",180.0*((i/n_pwI)%2)));

    for(i=0;i<n_evolve;++i)
    { Uevolve[i].set_DBR();
    }

    for(i=0;i<nspins;++i)
      detect[i].set_DBR();
    sigma.set_DBR();
    sigma1=sigma;
    for(i=0;i<Fnp;++i)
    { for(j=0;j<nspins;++j)
      { data(j,i) += proj(sigma1,detect[j])*scale;
      }
      for(j=0;j<n_cycle;++j)
        sigma1.sim_trans_ip(Uevolve[(i*n_cycle+j)%(n_evolve)]);
    }
  } // end of powder loop
  MATLAB(name1,name2,data,1);
  exit(0);
}
