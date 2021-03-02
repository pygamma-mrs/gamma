#include <gamma.h>
#include <vector>

int main()
  {
//                  This Is A Typical Dipolar Spatial Cartesian Tensor (KHz)

  matrix Dmx(3,3,h_matrix_type);
  Dmx.put  ( 0.79,0,0); Dmx.put_h( 0.00,0,1); Dmx.put_h( 0.00,0,2); 
                        Dmx.put  (-1.13,1,1); Dmx.put_h(-0.93,1,2); 
                                              Dmx.put  ( 0.34,2,2);

//              First We Diagonalize It To Determine delzz and eta
//       We Must Sort The Diagonal Values So That |Dzz| >= |Dyy| >= |Dxx|

  matrix Aeval, Aevec;			// Eigvenvalues, Eigenvectors
  diag(Dmx, Aeval, Aevec);		// Diagonalize Dmx
  double Axx = Aeval.getRe(0,0);	// Get diagonal Axx
  double Ayy = Aeval.getRe(1,1);	// Get diagonal Ayy
  double Azz = Aeval.getRe(2,2);	// Get diagonal Azz
  double tmp;
  if(fabs(Axx) > fabs(Azz)) { tmp=Azz; Azz=Axx; Axx=tmp; }
  if(fabs(Ayy) > fabs(Azz)) { tmp=Azz; Azz=Ayy; Ayy=tmp; }
  if(fabs(Ayy) > fabs(Axx)) { tmp=Ayy; Ayy=Axx; Axx=tmp; }
  double Aeta = (Axx-Ayy)/Azz;

  cout << "\n\n\tCartesian Tensor\n"    << Dmx; 
  cout << "\n\n\tDiagonalized Tensor\n" << Aeval; 
  cout << "\n\n\tThe trace value is:            " << trace(Aeval);
  cout << "\n\n\tThe delxx value appears to be: " << Axx;
  cout << "\n\n\tThe delyy value appears to be: " << Ayy;
  cout << "\n\n\tThe delzz value appears to be: " << Azz;
  cout << "\n\n\tThe asymmetry appears to be:   " << Aeta;

//           Solve For Angle Beta Assuming No Asymmetry

/*
                      1             2
               A   =  - del   [ 3cos (beta) - 1 ]
                zz    2    zz
        
                 2         1  [                  ]
              cos (beta) = -  | 2A   / del   + 1 |
                           3  [   zz      zz     ]                  */
         
  double Csqbeta = 2.0*Dmx.getRe(2,2)/Azz + 1.0; 
  double Cbeta = sqrt(Csqbeta);
  double beta = acos(Cbeta);
  cout << "\n\tCosine Of Beta Squared Is " << Csqbeta;
  cout << "\n\tCosine Of Beta Is " << Cbeta;
  cout << "\n\tAngle Beta Determined As " << beta*RAD2DEG << " Degrees";

/*
                      1         
               A   =  - del   [ 3sin(alpha)*sin(2beta) ]
                yz    4    zz
        
                           1  [              ]   /
              sin(alpha) = -  | 4A   / del   |  / sin(2.0*beta)
                           3  [   yz      zz ] /                 */
          
  double SalphaC2beta3 = 4.0*Dmx.getRe(1,2)/Azz; 
  double Salpha      = (1.0/3.0)*SalphaC2beta3/sin(2.0*beta);
  double alpha = asin(Salpha);
  cout << "\n\tThree Sin Alpha Cos 2 Beta Is " << SalphaC2beta3;
  cout << "\n\tSine Of Alpha Is " << Salpha;
  cout << "\n\tAngle Alpha Determined As " << alpha*RAD2DEG << " Degrees";


  coord X(Axx, Ayy, Azz);
  IntRank2A IR2A(X);
  IR2A.printCartesian(cout,alpha*RAD2DEG,beta*RAD2DEG);

  }
