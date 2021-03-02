/*

sl_exchange5.cc

3 spin simulation with exchange in Hilbert space

*/

#include "gamma.h"
//#include <sys/time.h>
//#include <sys/resource.h>

#define NPROP 1000
#define MAXSPINS 3

using namespace std;

int main(int argc, char *argv[])

{
  sys_dynamic ax;
  gen_op Hrf,Ham, H[5], sigma, sigma1, detect[MAXSPINS];
  super_op Iop,K,L,U[NPROP],UT;
  spin_T  Hdip[MAXSPINS][MAXSPINS];
  space_T Adip[MAXSPINS][MAXSPINS], Adip_R[MAXSPINS][MAXSPINS];
  double D[MAXSPINS][MAXSPINS];
  double iso_CSA[MAXSPINS];
  int i,j,k,k1,k2,Fnp,count,qu,steps,nrotor;
  string name, names;
  const double thetam=54.73561032;
  double mas_freq,gamB1;
  double ltime, time, time_prop, scale;
  int nstart, nspins, index, nprop;
  double alpha,beta,gamma;
  double alpha_D[MAXSPINS][MAXSPINS],beta_D[MAXSPINS][MAXSPINS];
  double gamma_D[MAXSPINS][MAXSPINS];
  double tauc, kex;
  //struct rusage me;
  //struct timeval tp;


  int value1[] = {1, 50, 100, 144, 200, 300, 538, 1154};
  int value2[] = {1,  7,  27,  11,  29,  37,  55,  107};
  int value3[] = {1, 11,  41,  53,  79,  61, 229,  271};

/*initialize the random number generator*/
  //gettimeofday(&tp,NULL);
  //srand48(tp.tv_usec);
  enable_blockdiag();

  count = 1;
  query_parameter(argc,argv,count++,"Spin System Name            ? ", names);
//setup for the spin system
  ax.read(names);
  nspins = ax.spins();
  if(nspins != 3)
  { cerr << "This program is written for a MAXSPINS spin system.\n";
    cerr << "Please change your spin system definition in \n";
    cerr << "the file " << names << ".\n";
    cerr << "Aborting ....\n\n";
    exit(-1);
  }
  for(i=0;i<2;++i)
  { for(j=i+1;j<3;++j)
    { query_parameter(argc,argv,count++,"Dipolar Coupling Constant   ? ", D[i][j]);
      query_parameter(argc,argv,count++,"Euler angle alpha           ? ", alpha_D[i][j]);
      query_parameter(argc,argv,count++,"Euler angle beta            ? ", beta_D[i][j]);
      query_parameter(argc,argv,count++,"Euler angle gamma           ? ", gamma_D[i][j]);
    }
  }
  for(i=0;i<3;++i)
  { query_parameter(argc,argv,count++,"Isotropic chemical shift    ? ", iso_CSA[i]);
  }
  query_parameter(argc,argv,count++,"Powder Quality (cheng)      ? ", qu);
  query_parameter(argc,argv,count++,"Number of time steps        ? ", steps);
  query_parameter(argc,argv,count++,"spinning speed              ? ", mas_freq);
  query_parameter(argc,argv,count++,"Number of rotor periods     ? ", nrotor);
  query_parameter(argc,argv,count++,"Number of sampling per rotor? ", nprop);
  query_parameter(argc,argv,count++,"Initial density operator    ? ", nstart);
  query_parameter(argc,argv,count++,"jump correlation time       ? ", tauc);
  query_parameter(argc,argv,count++,"rf-field amplitude          ? ", gamB1);
  query_parameter(argc,argv,count++,"Output Filename             ? ", name);

  if(nstart >= nspins)
  { cerr << "The starting density operator has a spin number which is\n";
    cerr << "larger than the number of spins.\n";
    cerr << "Aborting ....\n\n";
    exit(-1);
  }

  Fnp = nprop*nrotor;

  time_prop  = (1.0/mas_freq)/nprop;
  time       = (1.0/mas_freq)/steps;
  kex        =1.0/(2.0*tauc);

  cout << "\n\nSimulation of isotropic chemical shift by dipolar coupling\n";
  cout << "==========================================================\n\n";
  cout << "Parameters:\n";
  cout << "rotation angle thetam          : " << thetam << " Degree\n";
  cout << "size of spin system            : " << nspins << " spins\n";
  cout << "initial density operator       : Ix(ax," << nstart << ")\n"; 
  for(i=0;i<2;++i)
  { for(j=i+1;j<3;++j)
    { cout << "dipolar coupling constant (" << i << "," << j << ") : " <<
              D[i][j] << " Hz\n";
      cout << "relativ orientation of D tensor: (" << alpha_D[i][j] << "," <<
               beta_D[i][j] << "," << gamma_D[i][j] << ")\n";
    }
  }
  for(i=0;i<3;++i)
  { cout << "isotropic chemical shift (" << i << ") : " << iso_CSA[i] << " Hz\n";
  }
  cout << "# of sampling points:          " << nprop << "\n";
  cout << "# of rotor periods:            " << nrotor << "\n";
  cout << "Powder Quality Number:         " << qu << "  (" << value1[qu] <<
          " orientations)\n";
  cout << "Number of data points:         " << Fnp << " points\n";
  cout << "MAS frequency:                 " << mas_freq << " Hz\n";

	char timex[100];
#ifdef _MSC_VER
	sprintf_s ( timex, "%.6f", time);
#else
	sprintf ( timex, "%.6f", time);
#endif

  cout << "time increments:               " << timex << "s\n";

	char tprop[100];
#ifdef _MSC_VER
	sprintf_s ( tprop, "%.5f", time_prop);
#else
	sprintf ( tprop, "%.5f", time_prop);
#endif

  cout << "dwell time:                    " << tprop << "s\n";
  cout << "jump rate constant             " << kex << " s-1\n";
  cout << "correlation time               " << tauc << " s\n";
  cout << "rf-field amplitude:            " << gamB1 << " Hz\n";
  cout << "Output filename:               " << name << "\n";
  cout << "\n";
  cout.flush();

  ax.Kex(kex,0);
  matrix mx(ax.HS()*ax.HS(), ax.HS()*ax.HS(), i_matrix_type); // ExpLOp is the identity superop
  Iop = super_op(mx);

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

//setup for the space tensor
  matrix help(3,3,0);
  for(i=0;i<2;++i)
  { for(j=i+1;j<3;++j)
    { help.put_h(-1.0/2.0,0,0);
      help.put_h(-1.0/2.0,1,1);
      help.put_h( 1.0,2,2);
      help   = - (complex) D[i][j] * help;
      Adip[i][j] = A2(help);
      Adip[i][j] = Adip[i][j].rotate(alpha_D[i][j],beta_D[i][j],gamma_D[i][j]);
    }
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
           << alpha << "\tgamma = " << gamma << "\n";
      cout.flush();
    }
    scale = sin(beta/180.0*PI);
    sigma = iso_CSA[nstart]/sqrt(iso_CSA[nstart]*iso_CSA[nstart]+gamB1*gamB1)*Iz(ax,nstart) +
            gamB1/sqrt(iso_CSA[nstart]*iso_CSA[nstart]+gamB1*gamB1)*Ix(ax,nstart);
    for(k=0;k<nspins;++k)
    { detect[k] = Ix(ax,k);
    }
  
//now we rotate the space tensor
    for(i=0;i<2;++i)
    { for(j=i+1;j<3;++j)
      { Adip_R[i][j] = Adip[i][j].rotate(alpha,beta,gamma);
      }
    }

//zero all components
    for(i=0;i<5;++i)
      H[i] = gen_op();

//this is the dipolar part
    for(i=0;i<nspins;++i)
      H[2] += iso_CSA[i]*Iz(ax,i);
    for(i=0;i<2;++i)
    { for(j=i+1;j<3;++j)
      { for(k=-2;k<=2;++k)
        { H[k+2] += Adip_R[i][j].component(2,k) * d2(k,0,thetam) * Hdip[i][j].component(2,0);
        }
      }
    }
    for(i=0;i<nprop;++i)
      U[i] = Iop;

//now we calculate the propagator over one cycle of the MAS
    for(ltime=0.5*time;ltime<=steps*time;ltime += time)
    { index = int(ltime/time_prop);
      Ham = gamB1*Fx(ax,"1H");
      for(i=-2;i<=2;++i)
        Ham += exp(complex(0,i*2.0*PI*ltime*mas_freq)) * H[i+2];
      L = complex(0,2.0*PI)*commutator(Ham);
//    cout << "\nHam = \n";
//    cout << Ham;
//    cout << "\nL = \n";
//    cout << L;
      K = Kex(ax, Ham.get_basis());
//    cout << "\nK = \n";
//    cout << K;
//    cout.flush();
      L += K;
      UT = exp(L,-time);
      UT.set_HBR();
      U[index] = UT*U[index];
    }
    for(i=0;i<nprop;++i)
      U[i].set_DBR();
    sigma.set_DBR();
    for(i=0;i<Fnp;++i)
    { for(j=0;j<nspins;++j)
      { data(j,i) += proj(sigma,detect[j])*scale;
      }
      sigma=U[i%nprop]*sigma;
    }
  } // end of powder loop
  MATLAB(name1,name2,data,1);
  exit(0);
}
