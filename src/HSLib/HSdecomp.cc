/* HSdecomp.cc **************************************************-*-c++-*-
**                                  					**
**                              G A M M A				** 
**									**
**      Base Decomposition			Implementation		**
**									**
**      Copyright (c) 1991						**
**      Beat H. Meier, Scott A. Smith					**
**      Eidgenoessische Technische Hochschule				**
**      Labor fuer physikalische Chemie					**
**      8092 Zurich / Switzerland					**
**									**
**      $Header: $
**                		              				**
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** The GAMMA Platform Provides Functions for Simulation of Magnetic     ** 
** Resonance Experiments and Other Associated Mathematical              **
** Capabilities.  The Set of Functions Herein Provides For The		**
** Decomposition of a General Operator Into its Components Proportional	**
** to an Orthogonal Basis Set.   		                        **
**									**
*************************************************************************/

#ifndef   HSdecomp_cc_ 			// Is file already included?
#  define HSdecomp_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)               	// Using the GNU compiler
#    pragma implementation 		// this is the implementation
#  endif

#include <HSLib/HSdecomp.h>		// Include the header file
#include <HSLib/SpinOp.h>		// Knowledge of spin operators
#include <HSLib/SpinOpCmp.h>		// Know composite spin operators
#include <Basics/StringCut.h>		// Know about Gdec function

#define BD_SMALL 1.0e-15		// Floats < SMALL are zero

// used for rounding off to the nearest integer.
const double EPSILON_HSDECOMP = 0.000001;


void int_to_xbase(int retvec[], int x, const int base)
  {
  int z, i = 0;
  z = x;
  while ( z!=0)
    {
    retvec[i] = z % base;
    z = z / base;
    i++;
    }
  }


void Prod_base_dec(const spin_sys &sys, const gen_op &Op, double thres)

	// Input	sys  : Spin_system
	//		Op   : Operator to be decomposed
        //		thres: coefficients below thres
	//		       (absolute value) will be set to 0
	// Output	     : Decomposition is written to stdout

  {
  int nspins = sys.spins();		// # of spins in the system
  int i=0;
  for(i=0; i<nspins; i++)		// Check system for only Spin 1/2 
    if(sys.qn(i) != 0.5)
      {
      std::cout << "Prod_base_dec works only for spin 1/2: aborted\n";
      return;
      } 
 
  int* basef_nu;
  basef_nu = new int[nspins];	// Array of Base function codes
  spin_op B;				// Base operator
  std::string name;				// Base operator name
  int q;				// Index for number of spins

// not sure about destructor now
//  spin_op Ixyz[nspins][3];		// Array of Ix,Iy,Iz operators
  spin_op** Ixyz;
  Ixyz = new spin_op*[nspins];
  for(i=0; i<nspins; i++)
    {
    Ixyz[i] = new spin_op[3];
    Ixyz[i][0] = Ix(sys,i);
    Ixyz[i][1] = Iy(sys,i);
    Ixyz[i][2] = Iz(sys,i);
    }
       

   std::cout << "---Product base decomposition---------------------------\n"; 

     int n = nspins;
	 // DCT: 11/20/09... made changes here.
	 // for windows compatibility, and accuracy of comparison.
	 int fourToN;
	 fourToN = static_cast<int>(pow(4.0, n) + EPSILON_HSDECOMP);
   for ( i=1; i<fourToN; i++)	// Loop over all 4**N base operators
      { 
      for (int ii=0; ii<nspins; ii++)	// Initialize all base function
	basef_nu[ii]= 0;		// codes to contain zeros
      int_to_xbase(basef_nu, i, 4); 	// get base function code
					// for each spin: 0=none, 1=x, 2=y, 3=z
					// first spin is varying the fastest

      B = Ie(sys,0);			// Initialize base operator to unity
      name = std::string("");			// Initialize base operator name
      q=0;				// Initailize # of operators in product
      for(int j=nspins-1; j>-1; j--)	// Assemble product operator from code
       {				// Note : these are not yet normalized
       switch (basef_nu[j])  
         {
         case 0:			// No component for this spin
		name = "     " + name;
      		break;
         case 1:			// Ix component for this spin
		B *= Ixyz[j][0];
                if(nspins == 2)
                  {
                  if(j==0) name = std::string("Ix") + name;
                  else     name = std::string("Sx") + name;
                  }
                else
		  name = "Ix(" + Gdec(j) + ")" + name;
	 	q++;
   		break;
         case 2:			// Iy component for this spin
		B *= Ixyz[j][1];
                if(nspins == 2)
                  {
                  if(j==0) name = std::string("Iy") + name;
                  else     name = std::string("Sy") + name;
                  }
                else
		  name = "Iy(" + (std::string)Gdec(j) + ")" + name;
	 	q++;
		break;
         case 3:			// Iz component for this spin
		B *= Ixyz[j][2];
                if(nspins == 2)
                  {
                  if(j==0) name = std::string("Iz") + name;
                  else     name = std::string("Sz") + name;
                  }
                else
		  name = "Iz(" + (std::string)Gdec(j) + ")" + name; 
	 	q++;
		break;      
         default:
		std::cout << "error in Prod_base_dec: basef_nu["
	             << j << "]=" << basef_nu[j] <<" is illegal \n";
		break;
          }
        }

	double scaling = pow(2.0,q-1);	// B normalization via EBW (2.1.87)
//gen_op BOp = B;			// Projection of Op on unnormalized B
//gen_op BOp(B.matrix());			// Projection of Op on unnormalized B
gen_op BOp(B);			// Projection of Op on unnormalized B
gen_op OpX(Op);
complex coeff = proj(OpX,BOp);	// Projection of Op on unnormalized B
//	complex coeff = proj(Op,B);	// Projection of Op on unnormalized B
	coeff /= scaling;		// Add normalization into the projection

//	if (abs(Im(coeff)) > BD_SMALL)	// If non-zero immaginary part, flag
//	  std::cout << "complex coeficient for " << name << "!\n";

//	double coef = Re(coeff); 	// If real, non-zero coefficient, output
	if (norm(coeff) > thres) 
          {
	  std::cout << coeff << "*\t" ;
	  if (pow(2.0,q-1) != 1)
            std::cout << Gdec(int(pow(2.0,q-1))) << "*";
          std::cout << "\t" << name << "\n";
          }	
    }
  std::cout << "--------------------------------------------------------\n";
  delete [] basef_nu;
  delete [] Ixyz;
  }

#endif							// HSdecomp.cc
