/* CartMx2A.cc **************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Cartesian Matrix To Rank 2 Spatial Tensor           Interface	**
**                                                                      **
**      Scott Smith                                                     **
**      Copyright (c) 2001                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/
 
/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** This file contains functions that allow a 3x3 matrix representing a	**
** Cartesian rank 2 tensor to be converted into a GAMMA irredicuble	**
** rank 3 spherical tensor (class IntRank2A), an isotropic rank 0 	**
** component (Aiso), and a set of Euler angles that relate the rank 2	**
** PAS to the axes that were used to represent the input tensor.	**
**                                                                      **
** In short these function are intended to do the following:		**
**                                                                      **
**      A Cartesian                       A Spherical			** 
**   -----------------        ---------------------------------------	**
**                                                                      **
**   [ A    A    A   ]        A     = Scalar, Rank 0, Isotropic Term 	**
**   |  xx   xy   xz |         iso					**
**   |               |                                                  **
**   | A    A    A   |  ===>  del   = Scaling factor			**
**   |  yx   yy   yz |           zz					**
**   |               |                                                  **
**   | A    A    A   |        eta   = Asymmetry (in Class IntRank2)	**
**   [  zz   zy   zy ]                                                  **
**			      EA    = Euler Angles {alpha,beta,gamma}	**
**                                                                      **
** The values on the right are an equivalent representation of the 	**
** Cartesian values on the left. GAMMA can easily regenerate the array 	**
** by rotating the PAS representation back to the input orientation 	**
** using the Euler angles, then properly scaling it and adding back the **
** isotropic component. Using A	(eta, EA) for the rotated irreducible	**
** rank 2 spatial tensor       2					**
**									**
**       A Cartesian = 
**									**
**									**
*************************************************************************/

#ifndef   CartMx2A_cc_			// Is file already included?
#  define CartMx2A_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation 		// This is the implementation
#  endif

#include <IntRank2/CartMx2A.h>		// Include our interface
#include <Basics/Gconstants.h>		// Include RAD2DEG conversion
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <IntRank2/IntRank2A.h>		// Include spatial tensors
#include <string>			// Include stdlibc++ strings
#include <Basics/StringCut.h>		// Include Gform and Gdec functions
#include <stdlib.h>
#include <cmath>			// Inlcude HUGE_VAL_VAL

#if defined(_MSC_VER) || defined(__SUNPRO_CC)				// If we are using MSVC++
 #pragma warning (disable : 4786)       //   Kill STL namelength warnings
#endif

using std::string;			// Using libstdc++ strings
using std::vector;			// Using libstdc++ STL vectors
using std::list;			// Using libstdc++ STL lists
using std::cout;			// Using libstdc++ standard output
using std::ostream;			// Using libstdc++ output streams

typedef list<double> DoubleList;	// Simplfy usings lists of doubles

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                    Class CartMx2A Error Handling
// ____________________________________________________________________________

/*       Input                C2A     : Cartesian To GAMMA A Conversion (this)
                             eidx    : Flag for error type
                             noret   : Flag for return (0=return)
                             pname   : String included in message
        Output               none    : Error message
                                       Program execution stopped if fatal   */

void CartMx2A::C2Aerror(int eidx, int noret) const
  {
  string hdr("CartMx2A");
  string si(Gdec(Iter));
  switch(eidx)
    {
    case  0: GAMMAerror(hdr,"Program Aborting.....",       noret);break; // (0)
    case  2: GAMMAerror(hdr,"Problems During Construction",noret);break; // (2)
    case 10: GAMMAerror(hdr,"Max. Iterations Must Be > 19",noret);break; //(10)
    case 11: GAMMAerror(hdr,"Max. Iterations # Excessive", noret);break; //(11)
    case 12: GAMMAerror(hdr,"Array A Is Not Symmetric",    noret);break; //(12)
    case 13: GAMMAerror(hdr,"Cant Treat Rank 1 Components",noret);break; //(13)
    case 14: GAMMAerror(hdr,"Anti-Symmetric Terms In Mx",  noret);break; //(14)
    case 15: GAMMAerror(hdr,"Array Not Of Dimension 3x3",  noret);break; //(15)
    case 16: GAMMAerror(hdr,"Mx Cannot Be Rank 2 Tensor",  noret);break; //(16)
    case 17: GAMMAerror(hdr,"Cartesian Array Is Zero",     noret);break; //(17)
    case 18: GAMMAerror(hdr,"Anisotropy Is Zero",          noret);break; //(18)
    case 19: GAMMAerror(hdr,"Array Is Not Diagonal",       noret);break; //(19)
    case 20: GAMMAerror(hdr,"Value Left Unchanged",        noret);break; //(20)
    case 21: GAMMAerror(hdr,"Array Is Diagonal",           noret);break; //(21)
    case 25: GAMMAerror(hdr,"Asymmetry Is Non-Zero",       noret);break; //(25)
    case 26: GAMMAerror(hdr,"Improper Function Call",      noret);break; //(26)
    case 27: GAMMAerror(hdr,"Asymmetry Is Zero",           noret);break; //(27)
    case 30: GAMMAerror(hdr,"Cannot Obtain Euler Angles",  noret);break; //(30)
    case 35: GAMMAerror(hdr,"Cartesian Conversion Failure",noret);break; //(35)
    case 36: GAMMAerror(hdr,"Minimization Iterations:\n",  noret);break; //(36)
    case 37: GAMMAerror(hdr,"Can't Get Spherical Tensor\n",noret);break; //(37)
    case 39: GAMMAerror(hdr,"Please Report This Trouble..",noret);break; //(39)
    case 40: GAMMAerror(hdr,"Failed To Obtain Jacobian",   noret);break; //(40)
    case 41: GAMMAerror(hdr,"Current Iteration No.: " + si,noret);break; //(41)
    case 50: GAMMAerror(hdr,"Bad Jacobian Row Selected",   noret);break; //(50)
    case 51: GAMMAerror(hdr,"Jacobian Row Range [0,2]",    noret);break; //(51)
    case 52: GAMMAerror(hdr,"Bad Jacobian Funct. Selected",noret);break; //(52)
    case 53: GAMMAerror(hdr,"Jacobian Funct. Range [0,5]", noret);break; //(53)
    case 58: GAMMAerror(hdr,"Cannot Set Jacobian Function",noret);break; //(58)
    case 59: GAMMAerror(hdr,"Repeated Jacobian Function",  noret);break; //(59)
    case 60: GAMMAerror(hdr,"Could Not Handle Special A",  noret);break; //(60)
    case 61: GAMMAerror(hdr,"Some Zero Off-Diagonals In A",noret);break; //(61)
    case 70: GAMMAerror(hdr,"Error Converting Special A",  noret);break; //(70)
    case 71: GAMMAerror(hdr,"Case: Axy=Axz=0, Eta=Ayz=#",  noret);break; //(71)
    case 72: GAMMAerror(hdr,"Case: Axy=Ayz=0, Eta=Axz=#",  noret);break; //(72)
    case 73: GAMMAerror(hdr,"Case: Axz=Ayz=0, Eta=Axy=#",  noret);break; //(73)
    case 74: GAMMAerror(hdr,"Case: Axy=Axz=#, Eta=#",      noret);break; //(74)
    case 75: GAMMAerror(hdr,"Case: Axy=Ayz=#, Eta=#",      noret);break; //(75)
    case 76: GAMMAerror(hdr,"Case: Axz=Ayz=#, Eta=#",      noret);break; //(76)
    case 77: GAMMAerror(hdr,"Case: Axy = -Axz = #, Eta=#", noret);break; //(77)
    case 78: GAMMAerror(hdr,"Case: Axy = -Ayz = #, Eta=#", noret);break; //(78)
    case 79: GAMMAerror(hdr,"Case: Axz = -Ayz = #, Eta=#", noret);break; //(79)
    }
  }

volatile void CartMx2A::C2Afatal(int eidx) const
  {
  C2Aerror(eidx, 1);                            // Output error message
  if(eidx) C2Aerror(0);                         // State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

void CartMx2A::C2Aerror(int eidx, const string& pname, int noret) const
  {
  string hdr("CartMx2A");
  string si(Gdec(Iter));
  switch(eidx)
    {
    case 20: GAMMAerror(hdr,"Value Unchanged At "  +pname, noret);break; //(20)
    case 26: GAMMAerror(hdr,"Bad Function Call To "+pname, noret);break; //(26)
    case 54: GAMMAerror(hdr,"Cant Set Row <"+pname+"|J>",  noret);break; //(54)
    case 55: GAMMAerror(hdr,"Cant Use Auv Funct. " +pname, noret);break; //(55)
    }
  }

volatile void CartMx2A::C2Afatal(int eidx, const string& pname) const
  {
  C2Aerror(eidx, pname, 1);			// Output error message
  if(eidx) C2Aerror(0);                         // State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii                   Class CartMx2A SetUp Functions
// ____________________________________________________________________________

void CartMx2A::SetDefaults()
  {
  A         = matrix(3,3,0.0);		// Cartesian Array
  AisoVal   = 0.0;			// Isotropic Component
  PASDelzz  = 0.0;			// PAS delzz value
  EtaVal    = 0.0;			// Asymmetry (default is zero)
  EA        = EAngles(0,0,0);		// Euler Angles (no orientation)
 
  MaxAngs   = 5; 			// # each angle for min search
  MaxIter   = 500;			// Maximum # of Iterations Allowed
  MaxSSteps = 5;			// No. Repeated Minimizations Allowed
  Step      = 0;			// Current minimizations done
  Angle     = 0;			// Current angle set count
  Iter      = 0;			// Current iteration count
  JStep     = 0;			// Current Jacobian Evaluation
  SStep     = 0;			// Current Repeat Minimization

  J         = A;			// Set Jacobian to Zero
  F         = col_vector(6, 0.);	// Set Function Vector to zero
  Av        = vector<double> (6,0.);	// Set Vector Of A Components
  EAo       = EAngles(PI/4.,PI/3., 0.);	// Intitial Euler Angles Guess
  X         = col_vector(3, 0.0);	// Initialize Scratch X Vector
  Y         = X;			// Initialize Scratcy Y Vector

  EtaCut    = 5.e-3;			// Cutoff for eta zero
  DifCut    = 1.e-15;			// Cutoff for differential zero
  NormCut   = 1.e-7;			// Cutoff for norm |Y> zero
  ZeroMxCut = 1.e-6;			// A consistency cutoff

  Plevel    = 0;			// Minimization output level
  JOK       = true;			// Flag for Jacobian OK
  ConvMeth  = 0;			// Conversion method (none)
  Flist     = vector<int> (3);		// List of functions in J
  Flist[0]  = 0;			// Set first function to Axx
  Flist[1]  = 1;			// Set second function to Axy
  Flist[2]  = 2;			// Set third function to Axz
  Fuse      = "";			// No function use in J
  }

bool CartMx2A::SetAMx(const matrix& AC, int warn)
  {
//           First We Insure Cartesian Array Is Symmetric

  if(!AC.is_symmetric(DifCut))
    {
    if(warn)
      {
      C2Aerror(12, 1);				// Array A is not symmetric
      C2Aerror(14, 1);				// Anti-symmetric terms exist
      if(warn>1) C2Afatal(13);			// Can't Treat Rank 1 Comps.
      else       C2Aerror(13,1);		// Can't Treat Rank 1 Comps.
      }
    return false;
    }

//               Next We Insure Cartesian Array Is 3x3

  if(AC.rows() != 3)
    {
    if(warn)
      {
      C2Aerror(15, 1);				// Array A is not 3x3
      if(warn>1) C2Afatal(16);			// Can't Be A Rank 2 Tensor
      else       C2Aerror(16,1);
      }
    return false;
    }

//     All Appears In Order, Store Matrix A & 5 Components In Vector Av

  A = AC;
  matrix Aunsc = A-matrix(3,3,trace(A));// Remove isotropic part
  Av = vector<double> (6);		// Vector of A ==> |A>
  Av[0] = Aunsc.getRe(0,0);		// These are the irreducible
  Av[1] = Aunsc.getRe(0,1);		// rank 2 components of A
  Av[2] = Aunsc.getRe(0,2);
  Av[3] = Aunsc.getRe(1,1);		// GAMMA scaling is not used here
  Av[4] = Aunsc.getRe(1,2); 		// as it is unneccesary
  Av[5] = Aunsc.getRe(2,2);
  return true;
  }

void CartMx2A::SetAngles()
  {
//                   Input Angles { Alpha, Beta, Gamma }

  double alpha = EA.alpha();
  double beta  = EA.beta();
  double gamma = EA.gamma();

//                      Variables Involving Only Alpha

  Sa     = sin(alpha);				// sin(2a)
  Ca     = cos(alpha);				// cos(2a)
  S2a    = sin(2.0*alpha);			// sin(2a)
  C2a    = cos(2.0*alpha);			// cos(2a)

//                       Variables Involving Only Beta

  Cb     = cos(beta);				// sin(b)
  Sb     = sin(beta);				// cos(b)
  C2b    = cos(2.0*beta);			// sin(2b)
  S2b    = sin(2.0*beta);			// cos(2b)
  SbCb   = Sb*Cb;				// sin(b)*cos(b)
  Ssqb   = Sb*Sb;				// sin(b)*sin(b)
  Csqb   = Cb*Cb;				// cos(b)*cos(b)
  Csqbp1 = Csqb + 1.0;				// cos(b)*cos(b)+1

//                       Variables Involving Only Gamma

  S2g    = sin(2.0*gamma);			// sin(2g)
  C2g    = cos(2.0*gamma);			// cos(2g)
  }


// ____________________________________________________________________________
// iii          Class CartMx2A Jacobian Matrix Singularity Test
// ____________________________________________________________________________

/* Given a possible row in a 3x3 Jacobian matrix, <i|a|0>, <i|a|1>, & <i|a|2>,
   this function does a very simple check to insure the row will not produce
   a singlular Jacobian matrix.  It does so by checking that 1.) row i is not
   zero and 2.) row j where j<i is not a constant multiple of row i. The
   latter check is done for all rows j where j<i.                            */

bool CartMx2A::CheckSing(double ai0, double ai1, double ai2, int i)
  {
  if(ai0==0 && J.getRe(i,0) != 0) return false;
  if(ai0!=0 && J.getRe(i,0) == 0) return false;
  double sf = ai0/J.getRe(i,0);
  if(fabs(ai1/J.getRe(i,1) - sf) > 1.e-10) return false;
  if(fabs(ai2/J.getRe(i,2) - sf) > 1.e-10) return false;
  return true;
  }

// ____________________________________________________________________________
// iv         Class CartMx2A Caretesian Matrix Conversion Routine
// ____________________________________________________________________________

/* This attempts to perform the conversion of the current Cartesian tensor into
   a GAMMA spherical tensor. It is separated from other routiens and private
   because it uses the current class settings which will dictate how the 
   conversion process will be attempted.                                     */

bool CartMx2A::Convert(int warn)
  {
  TrackConv(0);				// Track conversion process (if needed)
  bool T = true;
  if(A.is_zero()) return T;		// If no A, we are done....
  AisoDelzEta();			// Set Aiso, delzz, & eta
  if(!PASDelzz) return T;		// If no anisotropy, were done
  if(!EtaVal) 				// If symmetric tensor, simple 
    {					// function for getting angles
    if(A.is_diagonal()) 		//   If diagonal, either in PAS
      { DiagSymCartEA(); return T; }	//   or simple rotation & were done
    else				//   If not diagonal, still easy
      { SymCartEA();     return T; }	//   to get Euler angles & were done
    }
  else					// If asymmetric, maybe some
    {					// structure we can use to simplify
    if(A.is_diagonal())			//   Asymmetric but diagonal, maybe
     { if(DiagASymCartEA()) return T; }	//   a simple way to Euler angles
    else if(OffDiagonals())		//   Asymmetric not diagonal, maybe
      { 				//   some off-diagonal zeros?
      if(ASymCartODZEA(warn?1:0)) 	//     If OK just return, we have
         return T;			//     found the Euler angles
      if(warn)				//     If not OK issue warnings
        {				//     if desired
        C2Aerror(60,1);			//       Cannot handle array
        if(warn>2) C2Afatal(61);	//       A with some off diagonals
        else       C2Aerror(61,1);	//       (done in ASymCartODZEA)
        }
      return false;
      }
    else if(SymOffDiags())		//   Asymmetric not diagonal, maybe
      {					//   symmetry in off-diagonals?
      if(ASymCartODSEA(warn?1:0)) 	//     If OK just return, we have
        return T;
else
return false;
//cout << "\n\tTrying Generic Conversion..";
      }
    if(ASymCartEA()) {      return T; }	// Nope, try hard way by fitting
    if(warn)				// Failed! What should we do? If
      {					// warnings are desired then we
      C2Aerror(35,1);			//   Conversion Failure
      if(warn > 2) C2Afatal(39);	//     If desired, exit
      else         C2Aerror(39,1);	//     else just issue warning
      }
    }
  return false;
  }

// ____________________________________________________________________________
// v           Class CartMx2A Jacobian Function List Check
// ____________________________________________________________________________

/* This insures that the mapping of {Auv} into a Jacobian row is not improper.
   The value of findex is: {0:Axx, 1:Axy, 2:Axz, 3:Ayy, 4:Ayz, 5:Azz} and we
   make sure its range is [0,5]. There are only three rows of the Jacobian
   so we insure that the range of frow is [0,2].                             */

bool CartMx2A::CheckF(int findex, int frow, int warn) const
  {
  if(frow<0 || frow>2)
    {
    if(warn)
      {
      C2Aerror(50, 1);			// Bad Jacobian row selected
      C2Aerror(51, 1);			// Jacobian row range [0,2]
      if(warn>2) C2Afatal(54,Gdec(frow));
      else       C2Aerror(54,Gdec(frow),1);
      }
    return false;
    }
  if(findex<0 || findex>5)
    {
    if(warn)
      {
      C2Aerror(52, 1);			// Bad Jacobian function selected
      C2Aerror(53, 1);			// Jacobian function range [0,5]
      if(warn>2) C2Afatal(55,Gdec(findex));
      else       C2Aerror(55,Gdec(findex),1);
      }
    return false;
    }
  return true;
  }

bool CartMx2A::CheckNorms() const
  {
  DoubleList::const_iterator nrm1, nrm2;	// Pixs into norms list
  DoubleList::const_iterator nend	;	// Pixs into norms list
  nrm1 = Norms.begin();				// Set to 1st norm
  nrm2 = Norms.begin();				// Set to 2nd norm
  nrm2++;
  nend = Norms.end();
  while(nrm2!=nend && nrm1!=nend) 
    {
    if(*nrm2 > *nrm1)
      if((*nrm2) - (*nrm1) > 0.50*(*nrm1))
        return true;
//        return false;
    nrm1++;
    nrm2++;
    }
  return true;
  }

// ____________________________________________________________________________
// vi          Class CartMx2A Euler Angles For Minimizaiton Search
// ____________________________________________________________________________

double CartMx2A::NewAlpha(int i) const
  {
  if(i<=0)       return   1.0*DEG2RAD;
  if(i>=MaxAngs) return 359.0*DEG2RAD;
  return (1.0 + 358.0*double(i)/double(MaxAngs-1))*DEG2RAD;
  }

double CartMx2A::NewBeta(int i) const
  {
  if(i<=0)       return   1.0*DEG2RAD;
  if(i>=MaxAngs) return 179.0*DEG2RAD;
  return (1.0 + 178.0*double(i)/double(MaxAngs-1))*DEG2RAD;
  }
 
double CartMx2A::NewGamma(int i) const
  {
  if(i<=0)       return   1.0*DEG2RAD;
  if(i>=MaxAngs) return 359.0*DEG2RAD;
  return (1.0 + 358.0*double(i)/double(MaxAngs-1))*DEG2RAD;
  }
 
vector<EAngles> CartMx2A::AngSeeds() const
  {
  int N = MaxAngs*MaxAngs*MaxAngs;
  double a, b, g;
  vector<EAngles> EAs(N);
  int i,j,k,ijk=0;
  for(i=0; i<MaxAngs; i++)
    {
    a = NewAlpha(i);
    for(j=0; j<MaxAngs; j++)
      {
      b = NewBeta(j);
      for(k=0; k<MaxAngs; k++,ijk++)
        {
        g = NewGamma(k);
        EAs[ijk] = EAngles(a,b,g);
        }
      }
    }
  return EAs;
  }

vector<int> CartMx2A::FctSeeds() const
  {
  vector<int> FSeeds;
  int i,j,k;
  for(i=0; i<4; i++)
    for(j=i+1; j<5; j++)
      for(k=j+1; k<6; k++)
        FSeeds.push_back(i+10*j+100*k);
  return FSeeds;
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A           Class CartMx2A Construction, Assignment, Destruction
// ____________________________________________________________________________


CartMx2A::CartMx2A() { SetDefaults(); }
CartMx2A::CartMx2A(const matrix& AC, int warn)
  {
  SetDefaults();			// Set default parameters
  if(!SetAMx(AC,warn?1:0)) C2Afatal(2); // Set AC as our Cartesian array		
  Convert(warn);
  }


void CartMx2A::operator= (const CartMx2A& C2A)
  {
  A         = C2A.A;			// Copy Cartesian tensor
  AisoVal   = C2A.AisoVal;		// Copy Isotropic value 
  PASDelzz  = C2A.PASDelzz;		// Copy PAS delzz
  EtaVal    = C2A.EtaVal;		// Copy Asymmetry eta
  EA        = C2A.EA;			// Copy Euler Angles

  MaxAngs   = C2A.MaxAngs;		// Copy Max. # of Angle Sets Allowed
  MaxIter   = C2A.MaxIter;		// Copy Max. # of Iteraction Allowed
  MaxSSteps = C2A.MaxSSteps;		// Copy Max. # of Repeat Steps
  Step      = C2A.Step;			// Copy Current Minimization Step
  Angle     = C2A.Angle;		// Copy Current Angle Set Count
  Iter      = C2A.Iter;			// Copy Current Iteration Count
  JStep     = C2A.JStep;		// Copy Current Jacobian Evaluation
  SStep     = C2A.SStep;		// Copy Current Repeat Step Count

  J         = C2A.J; 			// Copy Jacobian Matrix J   (3x3)
  F         = C2A.F;			// Copy Function Vector |F> (3x1)
  Av        = C2A.Av;			// Copy Vector |A> of 6 Auv Values
  EAo       = C2A.EAo;			// Copy Intitial Euler Angle Guesses
  X         = C2A.X;			// Copy Current Euler Angle Vector |X>
  Y         = C2A.Y;			// Copy Euler Angle Corrections    |Y>  
  EtaCut    = C2A.EtaCut;		// Copy Eta Zero Cutoff Value
  DifCut    = C2A.DifCut;		// Copy Differential Zero Cutoff
  NormCut   = C2A.NormCut;		// Copy Minimization Ending Cutoff
  ZeroMxCut = C2A.ZeroMxCut;		// Copy A Consistency Cutoff
  Plevel    = C2A.Plevel;		// Copy Print Level
  Norms     = C2A.Norms;		// Copy Iteration Norms
  JOK       = C2A.JOK;			// Copy Jacobian OK Flag		
  ConvMeth  = C2A.ConvMeth;		// Copy Conversion Method
  Flist     = C2A.Flist;		// Copy Function List
  Fuse      = C2A.Fuse;			// Copy String of Functions Used
  }

CartMx2A::~CartMx2A () { }

 
// ____________________________________________________________________________
// B                  Class CartMx2A Access Functions
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                   Functions To Obtain Class Values
// ----------------------------------------------------------------------------

matrix     CartMx2A::ACart()         const { return A;        }
double     CartMx2A::Aiso()          const { return AisoVal;  }
double     CartMx2A::delzz()         const { return PASDelzz; }
double     CartMx2A::Eta()           const { return EtaVal;   }
EAngles    CartMx2A::EulerAngles()   const { return EA;       }
int        CartMx2A::MaxIterations() const { return MaxIter;  }
int        CartMx2A::Iteration()     const { return Iter;     }
matrix     CartMx2A::Jacobian()      const { return J;        }
row_vector CartMx2A::Functional()    const { return F;        }
EAngles    CartMx2A::StartAngles()   const { return EAo;      }
double     CartMx2A::EtaCutoff()     const { return EtaCut;   }
double     CartMx2A::DifCutoff()     const { return DifCut;   }
int        CartMx2A::PrintLevel()    const { return Plevel;   }
int        CartMx2A::JFunct(int i)   const { return Flist[i]; }

// ----------------------------------------------------------------------------
//                     Functions To Set Class Values
// ----------------------------------------------------------------------------

bool CartMx2A::ACart(const matrix& AC, int warn)
  {
  if(!SetAMx(AC, warn?1:0))		// Try & Set AC as our array
    { 					// If that fails, then we will
    if(warn)				// issue warnings & perhaps quit
      {
      if(warn>1) C2Afatal(35);		//   Bad conversion, quit program
      else       C2Aerror(35);		//   Bad conversion, only warning
      }
    return false;
    }
  if(!CartMx2A::Convert(warn?1:0))	// Array OK, try conversions
    { 					// If that fails, then we will
    if(warn)				// issue warnings & perhaps quit
      {
      if(warn>1) C2Afatal(35);		//   Bad conversion, quit program
      else       C2Aerror(35);		//   Bad conversion, only warning
      }
    return false;
    }
  return true;
  }

void CartMx2A::MaxIterations(int mi)
  {
  if(mi<20)   { C2Aerror(10, 1); C2Aerror(20, Gdec(MaxIter), 1); return; }
  if(mi>5000) { C2Aerror(11, 1); C2Aerror(20, Gdec(MaxIter), 1); return; }
  MaxIter = mi;
  }

void CartMx2A::StartAngles(const EAngles& EAin) { EAo = EAin; }

void CartMx2A::EtaCutoff(double ec) { EtaCut = fabs(ec); }
void CartMx2A::DifCutoff(double dc) { DifCut = fabs(dc); }
void CartMx2A::PrintLevel(int   pl) { Plevel = abs(pl);  }
void CartMx2A::JFuncts(int f1, int f2, int f3)
  { 
  if(CheckF(f1,0,1)			// Insure all three functions
  && CheckF(f2,1,1)			// have the proper range
  && CheckF(f3,2,1))
    {
    if(f1==f2 || f1==f3 || f2==f3)	//   Insure no repeat functions
      {
      C2Aerror(59,1);			//     Repeated Jacobian function
      C2Afatal(58);			//     Cannot set Jacobian function
      }
    int tmp;				//   Insure successive ordering
    if(f1>f2) {tmp=f1; f1=f2; f2=tmp; }	//   This just simplifies picking
    if(f1>f3) {tmp=f1; f1=f3; f3=tmp; }	//   a set of functions to fill up
    if(f2>f3) {tmp=f2; f2=f3; f3=tmp; }	//   the Jacobian in a minimization
    Flist[0] = f1;
    Flist[1] = f2;
    Flist[2] = f3;
    }
  else C2Afatal(58);			// Cannot set Jacobian function
  }

void CartMx2A::JFuncts(int f1f2f3)
  {
  int f3 = 0;
  int tmp = f1f2f3;
  while(tmp >= 100) { tmp-=100; f3++; }
  int f1f2 = f1f2f3 - f3*100; 
      tmp = f1f2;
  int f2 = 0;
  while(tmp >= 10)  { tmp-=10;  f2++; }
  int f1   = f1f2   - f2*10;
  JFuncts(f1, f2, f3);
  }
 
// ____________________________________________________________________________
// C            CARTESIAN MATRIX ISOTROPY, ANISOTROPY, ASYMMETRY
// ____________________________________________________________________________

/* This function will glean the three values { Aiso, delzz, eta } from the
   input Cartesian rank 2 spatial tensor. The tensor is represented by a 3x3
   symmetric matrix. The values are determined
                                                                   A     A
          1                     2 [      1 [           ]]           xx -  yy
   A    = - Trace(A)    del   = - | A  - - | A   - A   ||    eta = ---------
    iso   3                zz   3 [  zz  2 [  xx    yy ]]              A
                                                                        zz

   where we adhere to the convention |A  | >= |A  | >= |A  | or eta = [0,1]
                                        zz       yy       xx                 */

void CartMx2A::AisoDelzEta(const matrix& AC)
  { ACart(AC); AisoDelzEta(); }

void CartMx2A::AisoDelzEta()
  {
  AisoVal = (1.0/3.0)*Re(trace(A));		// Get isotropic component
  if(fabs(AisoVal) < DifCut) AisoVal = 0.0;	// In sure zero if below delcut 
  matrix APAS, Aevect;                          // For diagonal mx & eigenvects
  diag(A, APAS, Aevect);			// Diagonalize A (put in PAS)
  double APASxx = APAS.getRe(0,0);              // Get Axx (PAS)
  double APASyy = APAS.getRe(1,1);              // Get Ayy (PAS)
  double APASzz = APAS.getRe(2,2);              // Get Azz (PAS)
  double tmp = 0;
  if(fabs(APASxx) > fabs(APASzz))               // Make sure we abide by rule
    { tmp=APASzz; APASzz=APASxx; APASxx=tmp; }
  if(fabs(APASyy) > fabs(APASzz)) 		//  |Azz| >= |Ayy| >= |Axx|
    { tmp=APASzz; APASzz=APASyy; APASyy=tmp; }
  if(fabs(APASxx) > fabs(APASyy))
    { tmp=APASyy; APASyy=APASxx; APASxx=tmp; }

  EtaVal   = (APASxx-APASyy)/APASzz;		// Calculate the asymmetry
  if(fabs(EtaVal) < EtaCut)       EtaVal = 0.0;	// -> Small values set zero
  if(fabs(EtaVal - 1.0) < EtaCut) EtaVal = 1.0;	// -> Set exactly 1 if close
  PASDelzz = APASzz - 0.5*(APASxx+APASyy);	// Calculate delzz value
  PASDelzz *= (2./3.);				//  THIS IS NOT GAMMA SCALED
  if(fabs(PASDelzz) < DifCut) PASDelzz = 0.0; 	//  (Small values set zero)
  if(EtaVal == 1) PASDelzz = fabs(PASDelzz);	// Use + delzz if eta is 1
  }						//   (where Ayy = - Azz)
 
// ____________________________________________________________________________
// D                      CARTESIAN MATRIX EULER ANGLES
// ____________________________________________________________________________

// ---------------------------------------------------------------------------- 
//           Euler Angles For Diagonal Irreducible Cartesian Tensor
// ---------------------------------------------------------------------------- 

/* This function will assume that the irreducible Cartesian spatial tensor
   array A (symmetric) has no asymmetry and is diagonal. That is, eta = 0 for A
   and A contains no non-zero off-diagonal elements. It will return three Euler
   angles EA:{alpha,beta,gamma} that relate to the tensor PAS to the input 
   array axes. In this case, i.e. eta=0 & A diagonal, either the tensor is
   already in its PAS or only a simple rotation is needed to switch diagonal
   element ordering so it matches the PAS.  Here are the possible Euler angles:

         EA = {0,0,0}            EA = {0,90,0}            EA = {90,90,0}

     1      [ -1  0  0 ]      1      [  2  0  0 ]      1      [ -1  0  0 ]
     - del  |  0 -1  0 |      - del  |  0 -1  0 |      - del  |  0  2  0 |
     2   zz [  0  0  2 ]      2   zz [  0  0 -1 ]      2   zz [  0  0 -1 ]  */
  
void CartMx2A::DiagSymCartEA()
  {
  string fct("DiagSymCartEA");				// Our function name
  ConvMeth = 1;						// Our conversion type
  TrackConv(ConvMeth);					// Track conv. process
  if(EtaVal)
    {
    C2Aerror(25, 1);					// Eta is not zero
    C2Aerror(26, fct, 1);				// Bad function call
    C2Afatal(30);					// Cant get EA
    }
  if(A.is_zero())					// Zero array
    {
    C2Aerror(17, 1);					// Cartesian array 0
    C2Aerror(26, fct, 1);				// Bad function call
    C2Afatal(30);					// Cant get EA
    }
  if(!PASDelzz)
    {
    C2Aerror(18, 1);					// Anisotropy 0
    C2Aerror(26, fct, 1);				// Bad function call
    C2Afatal(30);					// Cant get EA
    }
  if(!A.is_diagonal())					// Array not diagonal
    {
    C2Aerror(19, 1);					// Non-diag. array
    C2Aerror(26, fct, 1);				// Bad function call
    C2Afatal(30);					// Cant get EA
    }
  double axx = A.getRe(0,0) - AisoVal;
  double ayy = A.getRe(1,1) - AisoVal;
  double azz = A.getRe(2,2) - AisoVal;
  if(fabs(azz - PASDelzz) < DifCut) { EA = EAngles(0,      0     , 0); return; }
  if(fabs(axx - PASDelzz) < DifCut) { EA = EAngles(0,      PI/2.0, 0); return; }
  if(fabs(ayy - PASDelzz) < DifCut) { EA = EAngles(PI/2.0, PI/2.0, 0); return; }
  TrackConv(16);					// Track conv. process
  }

// ---------------------------------------------------------------------------- 
//           Euler Angles For Symmetric Irreducible Cartesian Tensor
// ---------------------------------------------------------------------------- 

/* This function will assume that the irreducible Cartesian spatial tensor
   array A (symmetric) has no asymmetry. That is, eta = 0 for A.  It will 
   return three Euler angles EA:{alpha,beta,gamma} that relate to the tensor
   PAS to the input array axes. In this case, i.e. eta=0, the Euler angle
   gamma is ill-defined and set to zero. The two other angles alpha and beta
   may be generated from the following:
                                                                           1
                                                   [ [   [   A         ] ] - ]
      1       [     2         ]                    | | 1 |    zz       | | 2 |
 A  = - del   | 3cos (beta)-1 |        beta = acos | | - | 2 ----- + 1 | |   |
  zz  2    zz [               ]                    | | 3 |   del       | |   |
                                                   [ [   [      zz     ] ]   ]

                                  [           A          ] 
                                  | 4          xz        |
                     alpha = acos | - ------------------ |
                                  | 3 del  * sin(2*beta) |
                                  [      zz              ]                   */


void CartMx2A::SymCartEA()
  {
double TFCut = 1.e-14;
  string fct("SymCartEA");				// Our function name
  ConvMeth = 2;						// Our conversion type
  TrackConv(ConvMeth);					// Track conv. process

//              Insure We Are Dealing With Asymmetric, Non-Zero Tensor

  if(EtaVal)
    {
    C2Aerror(25, 1);					// Eta is not zero
    C2Aerror(26, fct, 1);				// Bad function call
    C2Afatal(30);					// Cant get EA
    }
  if(A.is_zero())
    {
    C2Aerror(17, 1);					// Cartesian array 0
    C2Aerror(26, fct, 1);				// Bad function call
    C2Afatal(30);					// Cant get EA
    }
  if(!PASDelzz)
    {
    C2Aerror(18, 1);					// Anisotropy is 0
    C2Aerror(26, fct, 1);				// Bad function call
    C2Afatal(30);					// Cant get EA
    }

  double Axx = A.getRe(0,0) - AisoVal;			// Oriented Axx value
  double Axy = A.getRe(0,1) - AisoVal;			// Oriented Axy value
  double Axz = A.getRe(0,2) - AisoVal;			// Oriented Axz value
  double Ayy = A.getRe(1,1) - AisoVal;			// Oriented Ayy value
  double Ayz = A.getRe(1,2) - AisoVal;			// Oriented Ayz value
  double Azz = A.getRe(2,2) - AisoVal;			// Oriented Azz value

//              Determine  The Possible Values Of Beta From Azz
//         Use Of cos(beta)*cos(beta) Makes The Result MultiValued
//              Use Of acos Produces A Result Between [-90,90]

  double beta;						// Here is angle beta
  double alpha;						// Here is angle alpha
  beta  = acos(sqrt((1.0/3.0)*(2.0*Azz/PASDelzz+1.0)));	// Beta based on Azz
  beta  = fabs(beta);					// [-90,90] ==> [0,PI] 
  alpha = acos((4.0/3.0)*Axz/(PASDelzz*sin(2.0*beta)));	// Alpha based on Axz
  Sb    = sin(beta);					// sin(beta)
  S2b   = sin(2.0*beta);				// sin(2.0*beta)
  Ssqb  = Sb*Sb;					// sin(beta)*sin(beta)
  Cb    = cos(beta);					// cos(beta)
  Csqb  = Cb*Cb;					// cos(beta)*cos(beta)
  Ca    = cos(alpha);					// cos(alpha)
  S2a   = sin(2.0*alpha); 				// sin(2.0*alpha)
  C2a   = cos(2.0*alpha); 				// sin(2.0*alpha)

  if((fabs(Axy - 0.75*PASDelzz*S2a*Ssqb)                 < TFCut)
  && (fabs(Axz - 0.75*PASDelzz*Ca*S2b)                   < TFCut)
  && (fabs(Axx - 0.25*PASDelzz*(3.*C2a*Ssqb-3.*Csqb+1.)) < TFCut)
  && (fabs(Ayy + 0.25*PASDelzz*(3.*C2a*Ssqb+3.*Csqb-1.)) < TFCut))
    { EA = EAngles(alpha, beta, 0.0); return; }


  if(beta < PI/2.0) beta = PI-beta;			// Try alternate beta
  else              beta = PI/2.0-beta;
  alpha = acos((4.0/3.0)*Axz/(PASDelzz*sin(2.0*beta)));	// Alpha based on Axz
  Sb    = sin(beta);					// sin(beta)
  S2b   = sin(2.0*beta);				// sin(2.0*beta)
  Ssqb  = Sb*Sb;					// sin(beta)*sin(beta)
  Cb    = cos(beta);					// cos(beta)
  Csqb  = Cb*Cb;					// cos(beta)*cos(beta)
  Ca    = cos(alpha); 					// cos(alpha)
  S2a   = sin(2.0*alpha); 				// sin(2.0*alpha)
  C2a   = cos(2.0*alpha); 				// sin(2.0*alpha)

  if((fabs(Axy - 0.75*PASDelzz*S2a*Ssqb)                 < TFCut)
  && (fabs(Axz - 0.75*PASDelzz*Ca*S2b)                   < TFCut)
  && (fabs(Axx - 0.25*PASDelzz*(3.*C2a*Ssqb-3.*Csqb+1.)) < TFCut)
  && (fabs(Ayy + 0.25*PASDelzz*(3.*C2a*Ssqb+3.*Csqb-1.)) < TFCut))
    { EA = EAngles(alpha, beta, 0.0); return; }

// sosiz - don't think (or know if) the rest of this is used anymore?
  double beta2;
  if(beta < PI/2.0) beta2 = PI-beta;
  else              beta2 = PI/2.0-beta;
  beta2 = PI - beta + PI/2.0;				// Here is other beta

  if(beta2 > PI) beta2 = beta - PI/2.0;		// Both are possible
  double b1alphaxz, b2alphaxz;
  b1alphaxz = acos((4.0/3.0)*Axz/(PASDelzz*sin(2.0*beta)));
  b2alphaxz = acos((4.0/3.0)*Axz/(PASDelzz*sin(2.0*beta2)));
  double b1alphayz, b2alphayz;
  b1alphayz = asin((4.0/3.0)*Ayz/(PASDelzz*sin(2.0*beta)));
  b2alphayz = asin((4.0/3.0)*Ayz/(PASDelzz*sin(2.0*beta2)));
  double diff1 = fabs(b1alphaxz - b1alphayz);
  double diff2 = fabs(b2alphaxz - b2alphayz);
  if(diff1 < diff2)
    { EA = EAngles(b1alphaxz, beta, 0.0); }
  else
    { EA = EAngles(b2alphaxz, beta2, 0.0); }
  TrackConv(16);					// Track conv. process
  }

// ---------------------------------------------------------------------------- 
//           Euler Angles For Diagonal Asymmetric Cartesian Tensor
// ---------------------------------------------------------------------------- 

/* This function will assume that the irreducible Cartesian spatial tensor
   array A (symmetric) has non-zero asymmetry but that it is diagonal. That is,
   eta != 0 but A contains no off-diagonal elements. It will attempt to 
   generate three Euler angles EA:{alpha,beta,gamma} that relate to the tensor
   PAS to the input array axes. In this case, eta=# & A diagonal, there are
   five angle combinations that will maintain a symmetric tensor under rotation
   from the PAS. The function will return false if we do not catch the proper
   Euler angles.

                I: EA = {0,0,0}                   II: EA = {0,90,0}           

         1      [ eta-1   0    0 ]          1      [  2    0    0   ] 
         - del  |   0  -eta-1  0 |          - del  |  0 -eta-1  0   |
         2   zz [   0     0    2 ]          2   zz [  0    0  eta-1 ]


             III: EA = {0,90,180}                   II: EA = {0,90,0}           

         1      [-eta-1   0    0 ]          1      [  2    0    0   ] 
         - del  |   0   eta-1  0 |          - del  |  0 -eta-1  0   |
         2   zz [   0     0    2 ]          2   zz [  0    0  eta-1 ]

                        
                               III: EA = {alpha,0,0}

             1      [ eta*cos(2*alpha)-1           0          0 ]         
             - del  |         0          -eta*cos(2*alpha)-1  0 |
             2   zz [         0                    0          2 ]         
*/
  
bool CartMx2A::DiagASymCartEA()			
  {
  string fct("DiagASymCartEA");				// Our function name
  ConvMeth = 3;						// Our conversion type
  TrackConv(ConvMeth);					// Track conv. process
  if(!EtaVal)
    {
    C2Aerror(27, 1);					// Eta is zero
    C2Aerror(26, fct, 1);				// Bad function call
    C2Afatal(30);					// Cant get EA
    }
  if(A.is_zero())					// Zero array
    {
    C2Aerror(17, 1);					// Cartesian array 0
    C2Aerror(26, fct, 1);				// Bad function call
    C2Afatal(30);					// Cant get EA
    }
  if(!PASDelzz)
    {
    C2Aerror(18, 1);					// Anisotropy 0
    C2Aerror(26, fct, 1);				// Bad function call
    C2Afatal(30);					// Cant get EA
    }
  if(!A.is_diagonal())					// Array not diagonal
    {
    C2Aerror(19, 1);					// Non-diag. array
    C2Aerror(26, fct, 1);				// Bad function call
    C2Afatal(30);					// Cant get EA
    }

  double axx = A.getRe(0,0) - AisoVal;
  double ayy = A.getRe(1,1) - AisoVal;
  double azz = A.getRe(2,2) - AisoVal;
  double alpha;
  double dx, dy, dz;
  if(fabs(azz - PASDelzz) < DifCut)			// if PASDelzz = azz
    {							// Look for true PAS
    dx = -0.5*PASDelzz*(1.0-EtaVal);			// i.e. I: EA = {0,0,0}
    dy = -0.5*PASDelzz*(1.0+EtaVal);
    if(fabs(axx-dx)<DifCut && fabs(ayy-dy)<DifCut)
      { EA = EAngles(0,0,0); return true; }
    alpha = .5*acos((1./EtaVal)*(2.*axx/PASDelzz + 1.));// Look for
    dx = -0.5*PASDelzz*(1.0-EtaVal*cos(2.0*alpha));
    dy = -0.5*PASDelzz*(1.0+EtaVal*cos(2.0*alpha));
    if(fabs(axx-dx)<DifCut && fabs(ayy-dy)<DifCut)	// iLook for
      { EA = EAngles(alpha,0,0); return true; }
    alpha = 0.5*acos((1.0/EtaVal)*(-2.0*axx/PASDelzz + 1.0));
    dx =- 0.5*PASDelzz*(1.0-EtaVal*cos(2.0*alpha));
    dy =- 0.5*PASDelzz*(1.0+EtaVal*cos(2.0*alpha));
    if(fabs(axx-dx)<DifCut && fabs(ayy-dy)<DifCut)
      { EA = EAngles(alpha,PI,0); return true; }
    dx = -0.5*PASDelzz*(1.0-EtaVal);			// i.e. I: EA = {0,0,0}
    dy = -0.5*PASDelzz*(1.0+EtaVal);
    if(fabs(axx-dy)<DifCut && fabs(ayy-dx)<DifCut)	// Look for 3rd case
      { EA = EAngles(0,PI,PI/2.); return true; }	// III: EA = {0,180,90}
cout << "\nTROUBLE CONVERTING ASYMMETRIC TENSOR WITH Azz = PAS delzz";
    }
  if(fabs(axx - PASDelzz) < DifCut)			// If PASDelzz = axx
    {							// Look for 2nd case
    dy = -0.5*PASDelzz*(1.0+EtaVal); 			// II: EA = {0,90,0}
    dz = -0.5*PASDelzz*(1.0-EtaVal);
    if(fabs(ayy-dy)<DifCut && fabs(azz-dz)<DifCut)
      { EA = EAngles(0,PI/2.0,0); return true; }
    if(fabs(azz - 0.5*PASDelzz*(EtaVal-1.0)) < DifCut)
      { EA = EAngles(0,PI/2.0,0); return true; }
    if(fabs(azz - 0.5*PASDelzz*(-EtaVal-1.0)) < DifCut)
      { EA = EAngles(0,PI/2.0,PI/2.0); return true; }

    double gamma = 0.5*acos((1.0/EtaVal)*(2.0*azz/PASDelzz + 1.0));
           dy = -0.25*PASDelzz*(3.0+EtaVal*cos(2.0*gamma) + cos(2.0*gamma) - 1);
           dz =  0.50*PASDelzz*(EtaVal*cos(2.0*gamma) - 1.0);
    if(fabs(ayy-dy)<DifCut && fabs(azz-dz)<DifCut)
      { EA = EAngles(0,PI/2.0,gamma); return true; }
cout << "\nTROUBLE CONVERTING ASYMMETRIC TENSOR WITH Axx = PAS delzz";
    }
  if(fabs(ayy - PASDelzz) < DifCut)			// If PASDelzz = ayy
    {
    dx = -0.5*PASDelzz*(1.0+EtaVal); 			// IV: EA = {90,90,0}
    dz = -0.5*PASDelzz*(1.0-EtaVal);
    if(fabs(axx-dx)<DifCut && fabs(azz-dz)<DifCut)
      { EA = EAngles(PI/2.0,PI/2.0,0.0); return true; }
    if(fabs(axx-dz)<DifCut && fabs(azz-dx)<DifCut)
      { EA = EAngles(PI/2.0,PI/2.0,PI/2.0); return true; }
cout << "\nTROUBLE CONVERTING ASYMMETRIC TENSOR WITH Ayy = PAS delzz";
    }
  TrackConv(16);					// Track conv. process
  return false;
  }
// sosi


// ---------------------------------------------------------------------------- 
// Euler Angles For Asymmetric Cartesian Tensor With Some Off-Diagonals At Zero
// ---------------------------------------------------------------------------- 

bool CartMx2A::ASymCartODZEA(int warn)
  {
  string fct("ASymCartODZEA");				// Our function name
  ConvMeth = 4;						// Our conversion type
  TrackConv(ConvMeth);					// Track conv. process
  if(!EtaVal)
    {
    C2Aerror(27, 1);					// Eta is zero
    C2Aerror(26, fct, 1);				// Bad function call
    C2Afatal(30);					// Cant get EA
    }
  if(A.is_zero())					// Zero array
    {
    C2Aerror(17, 1);					// Cartesian array 0
    C2Aerror(26, fct, 1);				// Bad function call
    C2Afatal(30);					// Cant get EA
    }
  if(!PASDelzz)
    {
    C2Aerror(18, 1);					// Anisotropy 0
    C2Aerror(26, fct, 1);				// Bad function call
    C2Afatal(30);					// Cant get EA
    }
  if(A.is_diagonal())					// Array is diagonal
    {
    C2Aerror(21, 1);					// Diagonal array
    C2Aerror(26, fct, 1);				// Bad function call
    C2Afatal(30);					// Cant get EA
    }

  double axy = A.getRe(0,1) - AisoVal;			// Get irreducible
  double axz = A.getRe(0,2) - AisoVal;			// rank 2 components
  double ayz = A.getRe(1,2) - AisoVal;			// to work with
  double axx = A.getRe(0,0) - AisoVal;
  double ayy = A.getRe(1,1) - AisoVal;
  double azz = A.getRe(2,2) - AisoVal;

//                    Section When Axy = Axz = 0 (Ayz Non-Zero)

/*          These Are The Euler Angle Sets Looked For In This Section

       {  0, 90, 45}, {  0, 90,135}, { 90, 45,  0}, { 90, 45, 90}
       { 90,135,  0}, { 90,135, 90}
                                                                             */
// ---------------------------------------------------------------------------- 

//{0,45,0}, {0,45,90}, {0, 135, 90}, {180,45,0}

  if(fabs(axy)<DifCut && fabs(axz)<DifCut)
    {
    TrackConv(7);					// Track conv. process
    if(fabs(axx - PASDelzz) < DifCut)			// This @ {0, 90, 45} 
      {							// and  @ {0, 90,135}
      if((fabs(ayy + 0.50*PASDelzz)         < DifCut)
      && (fabs(azz + 0.50*PASDelzz)         < DifCut)
      && (fabs(ayz + 0.50*PASDelzz*+EtaVal) < DifCut))
        { EA = EAngles(0,PI/2.,PI/4.); return true; }
      if((fabs(ayy + 0.50*PASDelzz)         < DifCut)
      && (fabs(azz + 0.50*PASDelzz)         < DifCut)
      && (fabs(ayz - 0.50*PASDelzz*+EtaVal) < DifCut))
        { EA = EAngles(0,PI/2.,3*PI/4.); return true; }
      }

    if(fabs(ayz-0.25*PASDelzz*(3.0-EtaVal)) < DifCut)	// This @ {90,45,0}
      {
      if((fabs(axx+0.50*PASDelzz*( 1.0+EtaVal)) < DifCut)
      && (fabs(ayy-0.25*PASDelzz*( 1.0+EtaVal)) < DifCut)
      && (fabs(azz-0.25*PASDelzz*( 1.0+EtaVal)) < DifCut))
        { EA = EAngles(PI/2.,PI/4.,0.0); return true; }
      }

    if(fabs(ayz-0.25*PASDelzz*(3.0+EtaVal)) < DifCut)	// This @ {90,45,90}
      {
      if((fabs(axx+0.50*PASDelzz*( 1.0-EtaVal)) < DifCut)
      && (fabs(ayy-0.25*PASDelzz*( 1.0-EtaVal)) < DifCut)
      && (fabs(azz-0.25*PASDelzz*( 1.0-EtaVal)) < DifCut))
        { EA = EAngles(PI/2.,PI/4.,PI/2.); return true; }
      }

    if(fabs(ayz-0.25*PASDelzz*(-3.0+EtaVal)) < DifCut)	// This @ {90,135,0}
      {
      if((fabs(axx+0.50*PASDelzz*( 1.0+EtaVal)) < DifCut)
      && (fabs(ayy-0.25*PASDelzz*( 1.0+EtaVal)) < DifCut)
      && (fabs(azz-0.25*PASDelzz*( 1.0+EtaVal)) < DifCut))
        { EA = EAngles(PI/2.,3.*PI/4.,0.0); return true; }
      }

    if(fabs(ayz+0.25*PASDelzz*(3.0+EtaVal)) < DifCut)	// This @ {90,135,90}
      {
      if((fabs(axx+0.50*PASDelzz*( 1.0-EtaVal)) < DifCut)
      && (fabs(ayy-0.25*PASDelzz*( 1.0-EtaVal)) < DifCut)
      && (fabs(azz-0.25*PASDelzz*( 1.0-EtaVal)) < DifCut))
        { EA = EAngles(PI/2.,3.*PI/4.,PI/2.); return true; }
      }

/*
    if((fabs(axx-0.50*PASDelzz*(-1.0-EtaVal*cos(2.0*gamma))) < DifCut)
    && (fabs(ayy-0.25*PASDelzz*( 1.0+EtaVal*cos(2.0*gamma))) < DifCut))
      {
      if(fabs(ayz - 0.25*PASDelzz*(3.0-EtaVal)) < DifCut)
        { EA = EAngles(PI/2.,PI/4.,0.0); return true; }
      else if(fabs(ayz - 0.25*PASDelzz*(-3.0+EtaVal)) < DifCut)
        { EA = EAngles(PI/2.,3.0*PI/4.,0.0); return true; }
      }
    double gamma = 0.0;					//   or {90, 135, 0}
    gamma = PI/2.0;					// Try {90,45,90}
    if((fabs(axx-0.50*PASDelzz*(-1.0-EtaVal*cos(2.0*gamma))) < DifCut)
    && (fabs(ayy-0.25*PASDelzz*( 1.0+EtaVal*cos(2.0*gamma))) < DifCut))
      {
      if(fabs(ayz - 0.25*PASDelzz*(3.0+EtaVal)) < DifCut)
        { EA = EAngles(PI/2.,PI/4., gamma); return true; }
      else if(fabs(ayz - 0.25*PASDelzz*(-3.0-EtaVal)) < DifCut)
        { EA = EAngles(PI/2.,3.0*PI/4.,gamma); return true; }
      }

    gamma = 0.5*acos((4.*azz/PASDelzz-1)/EtaVal);
    if((fabs(axx-0.50*PASDelzz*(-1.0-EtaVal*cos(2.0*gamma))) < DifCut)
    &&(fabs(ayy-0.25*PASDelzz*(1.0+EtaVal*cos(2.0*gamma))) < DifCut))
      { EA = EAngles(PI/2.,PI/4.,gamma); return true; }
    gamma = 0.5*acos((4.*azz/PASDelzz-1)/EtaVal);
*/
    TrackConv(16); 					// Track conv. process
    if(warn)						// Failed conversion
      {							// If errors desired
      C2Aerror(70,1);					//   Bad conversion
      if(warn>1) C2Afatal(71);				//   Axy=Axz=0
      else       C2Aerror(71,1);
      }
    return false; 					// Did not succeed
    }

//                  Section When Axy = Ayz = 0 (Axz Non-Zero)

/*          These Are The Euler Angle Sets Looked For In This Section

       {  0, 45,  0}, {  0, 45, 90}, {  0, 135,  0}, {0, 135, 90}
       { 90, 90, 45}, { 90, 90,135}, {  0,beta,  0}                          */
// ---------------------------------------------------------------------------- 


  if(fabs(axy)<DifCut && fabs(ayz)<DifCut)	
    { 							
    TrackConv(8);					// Track conv. process
    if(fabs(axz - 0.25*PASDelzz*(3.0-EtaVal)) < DifCut)	// This @ {0, 45, 0} 
      {
      if((fabs(axx - 0.25*PASDelzz*(1.0+EtaVal)) < DifCut)
      && (fabs(ayy + 0.50*PASDelzz*(1.0+EtaVal)) < DifCut)
      && (fabs(azz - 0.25*PASDelzz*(1.0+EtaVal)) < DifCut))
        { EA = EAngles(0,PI/4.,0); return true; }
      }
    if(fabs(axz - 0.25*PASDelzz*(3.0+EtaVal)) < DifCut)	// This @ {0, 45, 90} 
      {
      if((fabs(axx - 0.25*PASDelzz*(1.0-EtaVal)) < DifCut)
      && (fabs(ayy + 0.50*PASDelzz*(1.0-EtaVal)) < DifCut)
      && (fabs(azz - 0.25*PASDelzz*(1.0-EtaVal)) < DifCut))
        { EA = EAngles(0,PI/4.0,PI/2.0); return true; }
      }
    if(fabs(axz + 0.25*PASDelzz*(3.0-EtaVal)) < DifCut)	// This @ {0, 135, 0} 
      {
      if((fabs(axx - 0.25*PASDelzz*(1.0+EtaVal)) < DifCut)
      && (fabs(ayy + 0.50*PASDelzz*(1.0+EtaVal)) < DifCut)
      && (fabs(azz - 0.25*PASDelzz*(1.0+EtaVal)) < DifCut))
        { EA = EAngles(0,3.*PI/4.,0); return true; }
      }
    if(fabs(axz + 0.25*PASDelzz*(3.0+EtaVal)) < DifCut)	// This @ {0, 135, 90} 
      {
      if((fabs(axx - 0.25*PASDelzz*(1.0-EtaVal)) < DifCut)
      && (fabs(ayy + 0.50*PASDelzz*(1.0-EtaVal)) < DifCut)
      && (fabs(azz - 0.25*PASDelzz*(1.0-EtaVal)) < DifCut))
        { EA = EAngles(0,3.*PI/4.,PI/2.); return true; }
      }
    if(fabs(axz - 0.50*PASDelzz*EtaVal) < DifCut)	// This @ {90, 90, 45} 
      {	
      if((fabs(axx + 0.50*PASDelzz) < DifCut)
      && (fabs(ayy - PASDelzz)      < DifCut)
      && (fabs(azz + 0.50*PASDelzz) < DifCut))
        { EA = EAngles(PI/2.,PI/2.,PI/4.); return true; }
      }
    if(fabs(axz + 0.50*PASDelzz*EtaVal) < DifCut)	// This @ {90, 90, 135} 
      {	
      if((fabs(axx + 0.50*PASDelzz) < DifCut)
      && (fabs(ayy - PASDelzz)      < DifCut)
      && (fabs(azz + 0.50*PASDelzz) < DifCut))
        { EA = EAngles(PI/2.,PI/2.,3.*PI/4.); return true; }
      }
    if(fabs(ayy + 0.50*PASDelzz*(1.0+EtaVal)) < DifCut)	// This @ {0, beta, 0} 
      {	
      double em3   = EtaVal - 3.0;
      double betao = 0.5*asin(-4.0*axz/(PASDelzz*(em3)));
      double beta  = betao;
      if((fabs(axx - 0.50*PASDelzz*(2.0+cos(beta)*cos(beta)*em3)) < DifCut)
      && (fabs(azz - 0.50*PASDelzz*(2.0+sin(beta)*sin(beta)*em3)) < DifCut))
        { EA = EAngles(0,beta,0); return true; }
      beta = PI/2.0 - betao;
      if((fabs(axx - 0.50*PASDelzz*(2.0+cos(beta)*cos(beta)*em3)) < DifCut)
      && (fabs(azz - 0.50*PASDelzz*(2.0+sin(beta)*sin(beta)*em3)) < DifCut))
        { EA = EAngles(0,beta,0); return true; }
      }


/*
    if(fabs(ayy - PASDelzz) < DifCut)			// This @ {90, 90, g} 
      {
      double g = 0.5*asin(2.0*axz/(PASDelzz*EtaVal));
      if((fabs(axx + 0.50*PASDelzz*(1.0+EtaVal*cos(2.0*g))) < 1.e-8)
      && (fabs(azz + 0.50*PASDelzz*(1.0-EtaVal*cos(2.0*g))) < 1.e-8))
        { EA = EAngles(PI/2.,PI/2.,g); return true; }
      }

    if((fabs(axz - PASDelzz) < DifCut) 			//   or {180, 45, 0}
    && (fabs(EtaVal - 1) < DifCut))
      { EA = EAngles(0,PI/4.,PI/2.); return true; }
    if((fabs(axz + PASDelzz) < DifCut) 
    && (fabs(EtaVal - 1) < DifCut))
      { EA = EAngles(PI,PI/4.,0.0); return true; }
*/
    if(warn)						// Failed conversion
      {							// If errors desired
      C2Aerror(70,1);					//   Bad conversion
      if(warn>1) C2Afatal(72);				//   Axy=Axz=0
      else       C2Aerror(72,1);
      }
    TrackConv(16);					// Track conv. process
    return false; 					// Did not succeed
    }

//                    Section When Axz = Ayz = 0 (Axy Non-Zero)
/*          These Are The Euler Angle Sets Looked For In This Section

    {0,0,45}, {0,0,135}, 
              {  0,  0, gamma}, {  0, 90,  gamma}, {  0,180,gamma}
*/

  if(fabs(axz)<DifCut && fabs(ayz)<DifCut)	
    {
    TrackConv(9);					// Track conv. process
    if(fabs(axy - 0.5*PASDelzz*EtaVal) < DifCut)	// This @ {0, 0, 45} 
      {
      if((fabs(axx + 0.50*PASDelzz) < DifCut)
      && (fabs(ayy + 0.50*PASDelzz) < DifCut)
      && (fabs(azz - PASDelzz)      < DifCut))
        { EA = EAngles(0,0,PI/4.); return true; }
      }
    if(fabs(axy + 0.5*PASDelzz*EtaVal) < DifCut)	// This @ {0, 0, 135} 
      {
      if((fabs(axx + 0.50*PASDelzz) < DifCut)
      && (fabs(ayy + 0.50*PASDelzz) < DifCut)
      && (fabs(azz - PASDelzz)      < DifCut))
        { EA = EAngles(0,0,3.*PI/4.); return true; }
      }
    if(fabs(azz - PASDelzz) < DifCut)			// This @ {0, 0/PI, g} 
      {
      double gamma = 0.5*acos((-2.0*ayy/PASDelzz -1.0)/EtaVal);
      if((fabs(axx + 0.5*PASDelzz*(1.0-EtaVal*cos(2.0*gamma)))   < DifCut)
      && (fabs(axy - 0.5*PASDelzz*EtaVal*sin(2.0*gamma)) < DifCut))
        { EA = EAngles(0,0,gamma); return true; }
      if((fabs(axx + 0.5*PASDelzz*(1.0-EtaVal*cos(2.0*gamma)))   < DifCut)
      && (fabs(axy + 0.5*PASDelzz*EtaVal*sin(2.0*gamma)) < DifCut))
        { EA = EAngles(0,PI,gamma); return true; }
      }

    if(warn)						// Failed conversion
      {							// If errors desired
      C2Aerror(70,1);					//   Bad conversion
      if(warn>1) C2Afatal(73);				//   Axz=Ayz=0
      else       C2Aerror(73,1);
      }
    TrackConv(16); 					// Track conv. process
    return false; 					// Did not succeed
    }
  TrackConv(10);					// Track conv. process
  return false;						// Did not succeed
  }

// ---------------------------------------------------------------------------- 
//  Euler Angles For Asymmetric Cartesian Tensor With Symmetric Off-Diagonals 
// ---------------------------------------------------------------------------- 

bool CartMx2A::ASymCartODSEA(int warn)
  {
  string fct("ASymCartODSEA");				// Our function name
  ConvMeth = 5;						// Our conversion type
  TrackConv(ConvMeth);					// Track conv. process
  if(!EtaVal)
    {
    C2Aerror(27, 1);					// Eta is zero
    C2Aerror(26, fct, 1);				// Bad function call
    C2Afatal(30);					// Cant get EA
    }
  if(A.is_zero())					// Zero array
    {
    C2Aerror(17, 1);					// Cartesian array 0
    C2Aerror(26, fct, 1);				// Bad function call
    C2Afatal(30);					// Cant get EA
    }
  if(!PASDelzz)
    {
    C2Aerror(18, 1);					// Anisotropy 0
    C2Aerror(26, fct, 1);				// Bad function call
    C2Afatal(30);					// Cant get EA
    }
  if(A.is_diagonal())					// Array is diagonal
    {
    C2Aerror(21, 1);					// Diagonal array
    C2Aerror(26, fct, 1);				// Bad function call
    C2Afatal(30);					// Cant get EA
    }

  double axy = A.getRe(0,1) - AisoVal;			// Get irreducible
  double axz = A.getRe(0,2) - AisoVal;			// rank 2 components
  double ayz = A.getRe(1,2) - AisoVal;			// to work with
  double axx = A.getRe(0,0) - AisoVal;
  double ayy = A.getRe(1,1) - AisoVal;


//                         Section When Axy = Axz
//                               {90,135,gamma}

  if(fabs(axy-axz)<DifCut)
    {
    double gamma = 0.5*acos((4.0*ayy/PASDelzz-1.0)/EtaVal);
    double fact = 0.5/sqrt(2.0);
    if((fabs(axy-fact*EtaVal*PASDelzz*sin(2.0*gamma))       < DifCut)
    && (fabs(ayz+0.25*PASDelzz*(3.0-EtaVal*cos(2.0*gamma))) < DifCut)
    && (fabs(axx+0.50*PASDelzz*(1.0+EtaVal*cos(2.0*gamma))) < DifCut))
        { EA = EAngles(PI/2.0,3.0*PI/4.0,gamma); return true; }
    if(gamma < PI) gamma += PI/2;
    else           gamma -= PI/2;       
    if((fabs(axy-fact*EtaVal*PASDelzz*sin(2.0*gamma))       < DifCut)
    && (fabs(ayz+0.25*PASDelzz*(3.0-EtaVal*cos(2.0*gamma))) < DifCut)
    && (fabs(axx+0.50*PASDelzz*(1.0+EtaVal*cos(2.0*gamma))) < DifCut))
        { EA = EAngles(PI/2.0,3.0*PI/4.0,gamma); return true; }
    if(warn)						// Failed conversion
      {							// If errors desired
      C2Aerror(70,1);					//   Bad conversion
      if(warn>1) C2Afatal(74);				//   Axz=Ayz=#
      else       C2Aerror(74,1);
      }
    TrackConv(16); 					// Track conv. process
    return false; 					// Did not succeed
    }

//                         Section When Axy = Ayz
//                               {0,135,gamma}

  if(fabs(axy-ayz)<DifCut)
    {
    double gamma = 0.5*acos((4.0*axx/PASDelzz-1.0)/EtaVal);
    double fact = 0.5/sqrt(2.0);
    if((fabs(axy+fact*EtaVal*PASDelzz*sin(2.0*gamma))       < DifCut)
    && (fabs(axz+0.25*PASDelzz*(3.0-EtaVal*cos(2.0*gamma))) < DifCut)
    && (fabs(ayy+0.50*PASDelzz*(1.0+EtaVal*cos(2.0*gamma))) < DifCut))
        { EA = EAngles(0,3.0*PI/4.0,gamma); return true; }
    if(gamma < PI) gamma += PI/2;
    else           gamma -= PI/2;       
    if((fabs(axy+fact*EtaVal*PASDelzz*sin(2.0*gamma))       < DifCut)
    && (fabs(axz+0.25*PASDelzz*(3.0-EtaVal*cos(2.0*gamma))) < DifCut)
    && (fabs(ayy+0.50*PASDelzz*(1.0+EtaVal*cos(2.0*gamma))) < DifCut))
        { EA = EAngles(0,3.0*PI/4.0,gamma); return true; }
    if(warn)						// Failed conversion
      {							// If errors desired
      C2Aerror(70,1);					//   Bad conversion
      if(warn>1) C2Afatal(75);				//   Axz=Ayz=#
      else       C2Aerror(75,1);
      }
    TrackConv(16); 					// Track conv. process
    return false; 					// Did not succeed
    }

//                         Section When Axz = Ayz

  if(fabs(axz-ayz)<DifCut)
    {
    if(warn)						// Failed conversion
      {							// If errors desired
      C2Aerror(70,1);					//   Bad conversion
      if(warn>1) C2Afatal(76);				//   Axz=Ayz=#
      else       C2Aerror(76,1);
      }
    TrackConv(16); 					// Track conv. process
    return false; 					// Did not succeed
    }

//                         Section When Axy = -Axz
//                               {90,45,gamma}

  if(fabs(axy+axz)<DifCut)
    {
    double gamma = 0.5*acos((4.0*ayy/PASDelzz-1.0)/EtaVal);
    double fact = 0.5/sqrt(2.0);
    if((fabs(axy+fact*EtaVal*PASDelzz*sin(2.0*gamma))       < DifCut)
    && (fabs(ayz-0.25*PASDelzz*(3.0-EtaVal*cos(2.0*gamma))) < DifCut)
    && (fabs(axx+0.50*PASDelzz*(1.0+EtaVal*cos(2.0*gamma))) < DifCut))
        { EA = EAngles(PI/2.0,PI/4.0,gamma); return true; }
    if(gamma < PI) gamma += PI/2;
    else           gamma -= PI/2;       
    if((fabs(axy+fact*EtaVal*PASDelzz*sin(2.0*gamma))       < DifCut)
    && (fabs(ayz-0.25*PASDelzz*(3.0-EtaVal*cos(2.0*gamma))) < DifCut)
    && (fabs(axx+0.50*PASDelzz*(1.0+EtaVal*cos(2.0*gamma))) < DifCut))
        { EA = EAngles(PI/2.0,PI/4.0,gamma); return true; }
    if(warn)						// Failed conversion
      {							// If errors desired
      C2Aerror(70,1);					//   Bad conversion
      if(warn>1) C2Afatal(77);				//   Axz=-Ayz=#
      else       C2Aerror(77,1);
      }
    TrackConv(16); 					// Track conv. process
    return false; 					// Did not succeed
    }

//                         Section When Axy = -Ayz
//                               {0,45,gamma}

  if(fabs(axy+ayz)<DifCut)
    {
    double gamma = 0.5*acos((4.0*axx/PASDelzz-1.0)/EtaVal);
    double fact = 0.5/sqrt(2.0);
    if((fabs(axy-fact*EtaVal*PASDelzz*sin(2.0*gamma))       < DifCut)
    && (fabs(axz-0.25*PASDelzz*(3.0-EtaVal*cos(2.0*gamma))) < DifCut)
    && (fabs(ayy+0.50*PASDelzz*(1.0+EtaVal*cos(2.0*gamma))) < DifCut))
        { EA = EAngles(0,PI/4.0,gamma); return true; }
    if(gamma < PI) gamma += PI/2;
    else           gamma-=PI/2;       
    if((fabs(axy-fact*EtaVal*PASDelzz*sin(2.0*gamma))       < DifCut)
    && (fabs(axz-0.25*PASDelzz*(3.0-EtaVal*cos(2.0*gamma))) < DifCut)
    && (fabs(ayy+0.50*PASDelzz*(1.0+EtaVal*cos(2.0*gamma))) < DifCut))
        { EA = EAngles(0,PI/4.0,gamma); return true; }
    if(warn)						// Failed conversion
      {							// If errors desired
      C2Aerror(70,1);					//   Bad conversion
      if(warn>1) C2Afatal(78);				//   Axz=-Ayz=#
      else       C2Aerror(78,1);
      }
    TrackConv(16); 					// Track conv. process
    return false; 					// Did not succeed
    }

//                         Section When Axz = -Ayz

  if(fabs(axy-ayz)<DifCut)
    {
    if(warn)						// Failed conversion
      {							// If errors desired
      C2Aerror(70,1);					//   Bad conversion
      if(warn>1) C2Afatal(79);				//   Axz=-Ayz=#
      else       C2Aerror(79,1);
      }
    TrackConv(16); 					// Track conv. process
    return false; 					// Did not succeed
    }

// sosi
  return false;						// Did not succeed
  }

// ---------------------------------------------------------------------------- 
//                 Euler Angles For Generic Cartesian Tensor
// ---------------------------------------------------------------------------- 

/* This function will assume that the Cartesian spatial tensor array A 
   (symmetric) has some asymmetry. That is, eta != 0 for A. It will also assume
   that the three tensor values {Aiso, delzz, eta} have been previously found.
   It will then attempt to obtain three Euler angles EA:{alpha,beta,gamma} that
   relate the irreducible rank 2 part of tensor A PAS to the input array axes.

   For this general case, the angles must be found through an iterative
   minimization process. Currently this uses Newtons Method (see function
   Minimize) which explicitly calculates Jacobian matrices and functional
   vectors. The Euler angles in Eo are used for an initial guess and placed
   into the vector |X>. This vector is then given to the minimization function
   which will return a revised |X> that contains the best fit values of the
   Euler angles. The "best fit" angles are then set to the tensor Euler angles. 

   Note that, since the Cartesian tensor elements are composed of multi-valued
   trigonometric functions, false minima can easily occur in an iterative 
   process to determine the Euler angles. Thus, the resulting Euler angles are
   ALWAYS checked by regeneration of the Cartesian matrix using the determined
   spherical PAS components {Ais, delzz, eta} and the minimized Euler angles
   {alpha,beta,gamma}, and a subsequent comparison with the original array A.

   If it is found that the minimization process has failed, there are three
   possible avenues the function automatically procede down. These are 

               1.) Try to minimize with new starting angles.
               2.) Try to minimize with more iterations.
               3.) Try to minimize on different Cartesian elements.

   the possibilities:

          Failure           
   =====================  ====================================================
   Did Not Reach Minimum  1.) If the vector norm of |Y> is small, then perhaps
                              it is simply very slow to converge on the EA set
                              Use more iterations and hope that works.
                          2.) If the vector norm of |Y> is not small, or more
                              iterations fail to generate a better norm, then
                              a different starting angle set should be used.
   Reached Wrong Minimum  1.) This can only be caused by falling into a local
                              minimum.  Retry the minimization with another
                              initial guess at the staring angle set.        */


bool CartMx2A::ASymCartEA()
  {
  string fct("ASymCartEA");			// Our function name
  ConvMeth = 6;					// Our conversion type
  TrackConv(ConvMeth);				// Track conv. process

//         Insure Array A Has The Proper Structure For This Function
//       (Or Lack Thereof, Any Structured A Handled In Other Routines)

  if(!EtaVal)
    {
    C2Aerror(27, 1);				// Eta is zero
    C2Aerror(26, fct, 1);			// Bad function call
    C2Afatal(30);				// Cant get EA
    }
  if(A.is_zero())
    {
    C2Aerror(17, 1);				// Cartesian array 0
    C2Aerror(26, fct, 1);			// Bad function call
    C2Afatal(30);				// Cant get EA
    }
  if(!PASDelzz)
    {
    C2Aerror(18, 1);				// Anisotropy 0
    C2Aerror(26, fct, 1);			// Bad function call
    C2Afatal(30);				// Cant get EA
    }

/*                         Minimization Section

    We Loop Over A Series Of Starting Angles And Iterate At Each Starting Guess
    To Get The Proper Euler Angles. Before Trying A New Guess At The Angles, We
    Try Another Three Functions (Of The Six Possible) To Build The Jacobian
    And Functional Vector.                                                   */

  bool TF = false;				// Flag if Done
  bool TFM; 					// Flag if Minization OK
  bool TFA;					// Check A consistency
  Step  = 0;					// No. Minimization tries
  SStep = 0;					// No. Repeat tries
  JStep = 0;					// No. J evaluations

//           These Handle Changes To The Three Auv Functions

  vector<EAngles> ASeeds = AngSeeds();		// Get our angle seeds
  vector<int>     FSeeds = FctSeeds();		// Get our function seeds
  int NA = ASeeds.size();			// Number of angle seeds
  int NJ = FSeeds.size();			// Number of J function sets			
  int action = 0;				// Action to take on next step
  double lastnorm = HUGE_VAL;			// Set |Y> norm, last step
  EAngles EAoo = EAo;				// Store starting angles
  while(!TF && Step<=NA)			// Begin minimization attampts
    {
    TrackMin(9);				//   Track if desired
    TrackConv(11); 				//   Track conv. process
    TFM = Minimize();	 			//   Minimize for EAngles
    TFA = Check() < ZeroMxCut;			//   Check A consistency
    TF  = (TFM & TFA);				//   Overall flag if OK
    action = MinAction(TFM,TFA,NJ,JStep,lastnorm);	//   Get next action to take
    switch(action)				//   Set for next minmize. step
      {
      case 0: 					// 0: Minimization complete!
        Step  = 0;				//   No Minimization steps
        JStep = 0;				//   No Jacobian evaluations
        SStep = 0;				//   No repeated steps now
        EAo = EAoo;				//   Back to original start angles
        return true;
        break;
      case 1:					// 1: Try New Angle Set
        EAo = ASeeds[Step++];			//   Try these starting EA's
        TrackMin(21);				//   Track if desired
        SStep = 0;				//   No repeated steps now
        lastnorm = HUGE_VAL;			//   Set last norm as big
        break;
      case 2:					// 2: Just try more iterations
        EAo = EA;				//   Set new start angles
        SStep++;				//   Inc. repeat step  count
        lastnorm = Y.norm();			//   Store last step norm
        break;
      case 3:					// 3: Try new function set
        TrackMin(22);				//   Track if desired
        TrackMin(47);				//   Track if desired
        JStep++;				//   Increment J evalution
        JFuncts(FSeeds[JStep]);			//   Set new J functions
        TrackMin(48);				//   Track if desired
        SStep = 0;				//   No repeated steps now
        lastnorm = HUGE_VAL;			//   Set last norm as big
        break;

      default:
       cout << "\n\nNo Action Yet Specified...";
       exit(-1);
       break;
      }
    }
  TrackConv(16); 				// Track conv. process
  return TF;					// Return if we succeeded
  }

int CartMx2A::MinAction(bool TFM, bool TFA, int NJ, int JStep, double lastnorm) const
  {
  if(TFM && TFA)			// Good Minimization, Self-Consistent
    { 
    TrackConv(12);			//   Track if desired (Min. Complete!)
    return 0;  				//   Recommended Action: Exit
    }
  if(!JOK)				// Bad Jacobian
    {
    TrackMin(30);			//   Track if desired (Bad Jacobian)
    if(JStep < NJ-1)			//   If # J evaluations not exceed
      return 3; 		 	//   Recommended Action: New J 
    return 1;			 	//   Recommended Action: New Angles
    }
  if(TFM && !TFA)			// Good Minimiz., Not Self-Consistent
    {
    TrackMin(20, Check());		//   Track if desired (Not Consistent)
    if(JStep < NJ-1)			//   If # J evaluations not exceed
      return 3; 		 	//   Recommended Action: New J 
    return 1;			 	//   Recommended Action: New Angles
    }
  if(!TFM)				// Bad Minimization (|Y> Norm Not Met)
    {
    TrackMin(11, lastnorm);		//   Track if desired (Norm Not Met)
    if(lastnorm-Y.norm() > NormCut)	//   If norm is small & slope decent
      {					//
      if(SStep < MaxSSteps)
        {
        TrackMin(13);			//   Track if desired
        return 2;			//   Recommended Action: More Iterations
        }
      else if(JStep < NJ-1)		//   If # J evaluations not exceed
        return 3; 		 	//   Recommended Action: New J 
      else
        {
        TrackMin(14);			//   Track if desired
        return 1;		 	//   Recommended Action: New Angles
        }
      }
    else				///  Minimization not minimizing
      {
      TrackMin(10);			//   Track if desired
      if(JStep < NJ-1)			//   If # J evaluations not exceed
        return 3; 		 	//   Recommended Action: New J 
      return 1;			 	//   Recommended Action: New Angles
      }
    }
  return 1;
  }

// ____________________________________________________________________________
// E                  Class CartMx2A Minimization Functions
// ____________________________________________________________________________

/* This function performs a Newton Minimization on an irreducible rank 2 
   Cartesian tensor to obtain a set of Euler angles {alpha, beta, gamma} that
   relate the tensor orientation to its PAS. Prior to the call of this function
   the user should set the Cartesian tensor A (which may be reducible), and
   have obtained the three components {Aiso, delzz, eta} independently.

   The minimizaton begins with an initial guess as to the Euler angles, EAo.
   It then iterates to find the "true" Euler angles, EA, using a minimization
   proceedure involving a Jacobian matrix and Function vector. "True" angles 
   are when the minimization angle corrections have a norm less than the value
   set for NormCut. The number of iterations will not be allowed to exceed
   MaxIter.

   Note that, since the Jacobian matrix and Functional vector are based on
   multi-valued trigonometric functions, false minima can easily occur. The
   resulting Euler angles should ALWAYS be checked by regeneration of the
   Cartesian matrix using spherical PAS components and the determined angle
   and subsequent comparison with the original array that was to be converted.
   The minimization process should be repeated with different starting angles
   if a bad angle set is produced even though the minimization has succeeded
   in getting the correction norm below NormCut. The same can be said when
   the minimization fails because it cannot find a minimium.                 */
   
bool CartMx2A::Minimize()
  {
  matrix Jinv;					// Jacobian Inverse
  double Ynorm = 0;				// Norm of vector <Y|
  double a, b, g;				// Adjusted angles
  EAngles EAY;
  TrackMin(0);					// Track minimization as desired
  EA = EAo;					// Initialize Euler angles
  X.put(EAo.alpha(), 0);			// Initialize the first pick
  X.put(EAo.beta(),  1);			// for the Euler angles
  X.put(EAo.gamma(), 2);
  Y = col_vector(3,complex0);
  for(int k=0; k<MaxIter; k++, Iter++)		// Iterate with minimization
    {
    TrackMin(1);				//   Track iteration if wanted
    if(!JacobianF()) { JOK=false; return JOK; }	//   Get Jacobian & Functional
    JOK = true;					//   Flag Jacobian Is OK
    Jinv = inv(J);				//   Generate inverse of J
    Y    = -1.0*(Jinv*F);			//   Generate update vector |Y>
    TrackMin(-1);				//   Track iteration as desired
    Ynorm = Y.norm();				//   Get norm of correction |Y>
    if(Ynorm < NormCut) 			//   If no angle corrections
      { TrackMin(2); TrackMin(3); return true; }//   we are finished iterating
    else 					//   Update angles next iterate 
      {
 //     if(Norms.size() >= 20) Norms.pop_front();	//     Keep 20 stored norms max
 //     Norms.push_back(Ynorm);			//     Store calculate |Y> norm
 //     if(!CheckNorms())				//     Insure norms decreasing
 //       { TrackMin(31); return false; }
      X += Y;					//     Update angles in X
      a=X.getRe(0); b=X.getRe(1); g=X.getRe(2);	//     Individual angle updates
      EA=EAngles(a,b,g);			//     Set update to EA
      }
    TrackMin(2);				//   Track iteration if needed
    }						//   Go to next iteration
  return false;					// Failed! Can't find angles
  }						// within MaxIter iterations
 
// ---------------------------------------------------------------------------- 
//          Class CartMx2A Jacobian Matrix & Functional Vector Generators
// ---------------------------------------------------------------------------- 

/* This function will generate a 3x3 Jacobian matrix and a 3x1 Functional
   vector for use in a Newton minimization to determine the three Euler angles
   EA:{alpha, beta, gamma} associated with an irreducible rank 2 Cartesian
   spatial tensor A. The angles relate the tensor PAS to the orientation in 
   which it exists. The Jacobian array J & Funcitonal vector F produced are


          [  @f1      @f1     @f1   ]                 [                      ]
          | ------   -----   ------ |                 | f1(alpha,beta,gamma) |
          | @alpha   @beta   @gamma |                 |                      |
          |                         |                 |                      |
          |  @f2      @f2     @f2   |                 |                      |
  J(EA) = | ------   -----   ------ |       |F(EA)> = | f2(alpha,beta,gamma) | 
          | @alpha   @beta   @gamma |                 |                      |
          |                         |                 |                      |
          |  @f3      @f3     @f3   |                 |                      |
          | ------   -----   ------ |                 | f3(alpha,beta,gamma) |
          [ @alpha   @beta   @gamma ]                 [                      ]


   The function takes the EAngles EA containing a possible set of Euler angles.
   Were these the true angles of the problem the functional vector will be zero
   by design, i.e. |F> = 0. In addition, the function uses the vector |Av>
   which contains 6 of the original components of the Cartesian matrix, A, being
   treated 
                       T
                   |Av>  = [ Axx, Axy, Axz, Ayy, Ayz, Azz ]
  
   Also used are the values of delzz (PASDelzz) and eta (EtaVal), both of 
   which may be determined directly from the Cartesian array A. The function
   takes an optional print flag allows the user to view the process involved
   in generating J and |F>.

   There are 6 possible functions that may be used for generating the rows of
   J and |F>, one for each input Cartesian element of A (Auv in vector |Av>). 
   Hard-coded in are the formula for generating the elements of the Jacobian
   and functional vector for all six. However, the minimization to obtain
   the Euler angles demands we use only three of these. Hence, the routine 
   will select out three of the six that are deemed good functions.... namely
   ones that will not produce a singular Jacobian matrix. In this problem a
   singular Jacobian matrix could easily occur since the six functions are
   inter-related and may produce identical rows in J or may produce a zero
   row in J because on often finds zero angles in the set EA.

   For more information on the Jacobian matrix, Functional vector F, and 
   overall minimization scheme see the GAMMA documentation on class IntRank2A.


	   Input	*	Av	: Vector - { Axx,Axy,Axz,Ayy,Ayz,Azz }
			*	PASDelzz: PAS delzz value
	  		*	EtaVal  : Asymmetry value
	   		*	EA	: Euler Angles {alpha,beta,gamma} rad
	  		*	F       : Functional Vector - {f1,f2,f3}
	  		*	J	: Jacobian Matrix (3x3)
	  			pl      : Print level for output
	   Output		TF      : Given the values in EA & |Av> the
	  				  vector |F> and Jacobian matrix J
	  				  are generated.  The function will
	  				  return true if this procedure is
	  				  deemed OK, false if it appears that
	  				  one cannot build a proper Jacobian
	   Note				: * Input variables are class internal
	   Note				: USED FOR FITTING {alpha,beta,gamma}
	  				  WHEN WE KNOW {delzz,eta}
	   Note				: The values in Av must come from a
					  symmetric traceless matrix
	   Note				: If the function returns false, one
	  				  should try a different set of 
                                          Euler Angles.                      */

 
void CartMx2A::Axx(double& f, double& dfda, double& dfdb, double& dfdg) const
  {
  double term1, term2, term3;

//                    Generate Function Value f Based On Axx  

  term1 = C2a*(3.0*Ssqb + EtaVal*Csqbp1*C2g);	// cos(2a)*{3*sin(b)*sin(b) +
						// eta*[cos(b)*cos(b)+1]cos(2g)}
  term2 = -2.0*EtaVal*S2a*Cb*S2g;		// -2*eta*sin(2a)*cos(b)*sin*2g)
  term3 = -3.0*Csqb + 1.0 - EtaVal*Ssqb*C2g;	// -3*cos(b)*cos(b)+1 + 
						//    eta*sin(b)*sin(b)*cos(2g) 
  f = term1+term2+term3;			// This is unscaled fxx

//                 Generate Jacobian Value df/da Based On Axx  

  term1 = -2.0*S2a*(3.0*Ssqb+EtaVal*Csqbp1*C2g);// -2*sin(2a)*{3*sin(b)*sin(b)+
						// eta*[1+cos(b)*cos(b)]cos(2g)}
  term2 = -4.0*EtaVal*C2a*Cb*S2g;		// -4*eta*sin(2a)*cos(b)*sin(2g)
  dfda =  term1 + term2;			// This is unscaled @fxx/@alpha

//                 Generate Jacobian Value df/db Based On Axx  

  term1 =  C2a*(6.0*SbCb - 2.*EtaVal*SbCb*C2g);	// sin(2a)*{6*sin(b)*cos(b) -
						// 2*eta*sin(b)*cos(b)*cos(2g)}
  term2 =  2.0*EtaVal*S2a*Sb*S2g;		// 2*eta*cos(2a)*sin(b)*sin(2g)
  term3 =  6.0*SbCb - 2.0*EtaVal*SbCb*C2g;	// 6*cos(b)*sin(b) - 
						//  2*eta*sin(b)*cos(b)*cos(2g)
  dfdb =  term1 + term2 + term3;		// This is unscaled @fxx/@beta

//                 Generate Jacobian Value df/dg Based On Axx  
      
  term1 = -2.0*EtaVal*C2a*Csqbp1*S2g;		// -2.0*eta*sin(2a)*
						//    [cos(b)*cos(b)+1]*sin(2g)
  term2 = -4.0*EtaVal*S2a*Cb*C2g;		//  4*eta*cos(2a)*cos(b)*cos
  term3 =  2.0*EtaVal*Ssqb*S2g;			//  2*eta*sin(b)*sin(b)*sin(2g)
  dfdg = term1 + term2;				// This is unscaled @fxx/@gamma
  }
 

void CartMx2A::Axy(double& f, double& dfda, double& dfdb, double& dfdg) const
  {
  double term1, term2;

//                    Generate Function Value f Based On Axy  

  term1 = S2a*(3.0*Ssqb + EtaVal*Csqbp1*C2g);	// sin(2a)*{3.0*sin(b)*sin(b) +
						// eta*[cos(b)*cos(b)+1]cos(2g)}
  term2 = 2.0*EtaVal*C2a*Cb*S2g;		// 2*eta*cos(2a)*cos(b)*sin(2g)
  f = term1 + term2;				// This is unscaled function fxy

//                 Generate Jacobian Value df/da Based On Axy  

  term1 =  2.0*C2a*(3.0*Ssqb+EtaVal*Csqbp1*C2g);// 2.0*cos(2a)*{3.0*sin(b)*sin(b)
						// +eta*[1+cos(b)*cos(b)]cos(2g)}
  term2 = -4.0*EtaVal*S2a*Cb*S2g;		// -4.0*eta*sin(2a)*cos(b)*sin(2g)
  dfda =  term1 + term2;			// This is unscaled @fxy/@alpha

//                 Generate Jacobian Value df/db Based On Axy  

  term1 =  S2a*(6.0*SbCb - 2.0*EtaVal*SbCb*C2g);// sin(2a)*{6.0*sin(b)*cos(b) -
						// 2.0*eta*sin(b)*cos(b)*cos(2g)}
  term2 = -2.0*EtaVal*C2a*Sb*S2g;		// -2.0*eta*cos(2a)*sin(b)*sin(2g)
  dfdb =  term1 + term2;			// This is unscaled @fxy/@beta

//                 Generate Jacobian Value df/dg Based On Axy  
      
  term1 = -2.0*EtaVal*S2a*Csqbp1*S2g;		// -2.0*eta*sin(2a)
						// *[cos(b)*cos(b)+1]*sin(2g)
  term2 =  4.0*EtaVal*C2a*Cb*C2g;		//  4.0*eta*cos(2a)*cos(b)*cos
  dfdg = term1 + term2;				// This is unscaled @fxy/@gamma
  }

void CartMx2A::Axz(double& f, double& dfda, double& dfdb, double& dfdg) const
  {
  double term1, term2;

//                    Generate Function Value f Based On Axz  

  term1 = Ca*(3.0*S2b - 2.0*EtaVal*SbCb*C2g);	// cos(a)*{3.0*sin(2b) -
						// 2*eta*sin(b)*cos(b)*cos(2g)}
  term2 = 2.0*EtaVal*Sa*Sb*S2g;			// 2*eta*sin(a)*sin(b)*sin(2g)
  f = term1 + term2;				// This is unscaled function fxz
      
//                 Generate Jacobian Value df/da Based On Axz  

//sosi switched sign here before 2.0*EtaVal
  term1 = -1.0*Sa*(3.0*S2b-2.0*EtaVal*SbCb*C2g);// -1.0*sin(a)*{3.0*sin(2b) -
						// 2.0*eta*sin(b)*cos(b)cos(2g)}
  term2 =  2.0*EtaVal*Ca*Sb*S2g;		//  2.0*eta*cos(a)*sin(b)*sin(2g)
  dfda =  term1 + term2;			// This is unscaled @fxz/@alpha

//                 Generate Jacobian Value df/db Based On Axz  
      
  term1=Ca*(6.*C2b-2.*EtaVal*(1.-2.*Ssqb*C2g));	// cos(a)*{6.0*cos(2b) - 
						// 2*eta*(1-2*sin(b)*sin(b))*cos(2g)}
  term2 = 2.0*EtaVal*Sa*Cb*S2g;			// 2.0*eta*sin(a)*cos(b)*sin(2g)
  dfdb =  term1 + term2;			// This is unscaled @fxz/@beta

//                 Generate Jacobian Value df/dg Based On Ayz  
      
  term1 =  4.0*EtaVal*Ca*SbCb*S2g;		//  4*eta*cos(a)*sin(b)*cos(b)*sin(2g)
  term2 =  4.0*EtaVal*Sa*Sb*C2g;		//  4*eta*sin(a)*sin(b)*cos(2g)
  dfdg = term1 + term2;				// This is unscaled @fxz/@gamma
  }

void CartMx2A::Ayy(double& f, double& dfda, double& dfdb, double& dfdg) const
  {
  double term1, term2, term3;

//                    Generate Function Value f Based On Ayy  

  term1 = -C2a*(3.0*Ssqb + EtaVal*Csqbp1*C2g);	// -cos(2a)*{3.0*sin(b)*sin(b) +
						// eta*[cos(b)*cos(b)+1]cos(2g)}
  term2 =  2.0*EtaVal*S2a*Cb*S2g;		//  2*eta*sin(2a)*cos(b)*sin*2g)
  term3 = -3.0*Csqb + 1.0 - EtaVal*Ssqb*C2g;	// -3*cos(b)*cos(b) + 1
						// + eta*sin(b)*sin(b)*cos(2g) 
  f = term1+term2+term3;			// This is unscaled fyy

//                 Generate Jacobian Value df/da Based On Ayy  
      
  term1 =  2.*S2a*(3.*Ssqb + EtaVal*Csqbp1*C2g);// 2.*sin(2a)*{3.0*sin(b)*sin(b)
						// +eta*[1+cos(b)*cos(b)]cos(2g)}
  term2 =  4.0*EtaVal*C2a*Cb*S2g;		// 4*eta*sin(2a)*cos(b)*sin(2g)
  dfda =  term1 + term2;			// This is unscaled @fyy/@alpha

//                 Generate Jacobian Value df/db Based On Ayy  
      
  term1 = -C2a*(6.0*SbCb - 2.0*EtaVal*SbCb*C2g);// -cos(2a)*{6*sin(b)*cos(b) -
						// 2*eta*sin(b)*cos(b)*cos(2g)}
  term2 = -2.0*EtaVal*S2a*Sb*S2g;		// -2*eta*cos(2a)*sin(b)*sin(2g)
  term3 =  6.0*SbCb - 2.0*EtaVal*SbCb*C2g;	// 6*cos(b)*sin(b) - 
						// 2*eta*sin(b)*cos(b)*cos(2g)
  dfdb =  term1 + term2 + term3;		// This is unscaled @fyy/@beta

//                 Generate Jacobian Value df/dg Based On Ayy  
      
  term1 =  2.0*EtaVal*C2a*Csqbp1*S2g;		//  2*eta*sin(2a)*
						//  [cos(b)*cos(b)+1]*sin(2g)
  term2 =  4.0*EtaVal*S2a*Cb*C2g;		//  4*eta*cos(2a)*cos(b)*cos
  term3 =  2.0*EtaVal*Ssqb*S2g;			//  2*eta*sin(b)*sin(b)*sin(2g)
  dfdg = term1 + term2;				// This is unscaled @fyy/@gamma
  }

void CartMx2A::Ayz(double& f, double& dfda, double& dfdb, double& dfdg) const
  {
  double term1, term2, term3;

//                    Generate Function Value f Based On Ayz  

  term1 =  Sa*(3.0*S2b - 2.0*EtaVal*SbCb*C2g);	// sin(a)*[3.*sin(2b)-
						// 2.*eta*sin(b)*cos(b)*cos(2g)]
  term2 = -2.0*EtaVal*Ca*Sb*S2g;		// -2*eta*cos(a)*sin(b)*sin(2g)
  f = term1 + term2;				// This is unscaled function fyz

//                 Generate Jacobian Value df/da Based On Ayz  
      
  term1 =  Ca*(3.0*S2b - 2.0*EtaVal*SbCb*C2g);	// cos(a)*[3.*sin(2b)-
						// 2*eta*sin(b)*cos(b)*cos(2g)]
  term2 =  2.0*EtaVal*Sa*Sb*S2g;		// 2*eta*sin(a)*sin(b)*sin(2g)
  dfda = term1 + term2;				// This is is unscaled @fyz/@alpha

//                 Generate Jacobian Value df/db Based On Ayz  
      
  term1 = Sa*(-6.*C2b-2.*EtaVal*(1.-Ssqb)*C2g);	// sin(a)*{-6*cos(2b)-
						// 2*eta*[1-sin(b)*sin(b)]*cos(2g)}
  term2 = 2.0*EtaVal*Ca*Cb*S2g;			// 2*eta*cos(a)*cos(b)*cos(2g)
  dfdb = term1 + term2;				// This is unscaled @fyz/@beta

//                 Generate Jacobian Value df/dg Based On Ayz  
      
  term1 =  4.0*EtaVal*Sa*SbCb*S2g;		// 4*eta*sin(a)*sin(b)*cos(b)*sin(2g)	
  term3 = -4.0*EtaVal*Ca*Sb*C2g;		// -4*eta*cos(a)*sin(b)*cos(2g)
  dfdg = term1 + term2;				// This is unscaled @fyz/@gamma
  }

void CartMx2A::Azz(double& f, double& dfda, double& dfdb, double& dfdg) const
  {
  f    = 6.0*Csqb - 2.0 + 2.0*EtaVal*Ssqb*C2g;	// Unscaled function fzz
  dfda = 0.0;					// Unscaled @fzz/@alpha
  dfdb = SbCb*(4.0*EtaVal*C2g-12.0);		// Unscaled @fzz/@beta
  dfdg = -4.0*EtaVal*Ssqb*S2g; 			// Unscaled @fzz/@gamma
  }



// ---------------------------------------------------------------------------- 

bool CartMx2A::JacobianF()
  {

// ----------------------------------------------------------------------------
//              Set Up Angle Relatioships Common To F & J Components
// ----------------------------------------------------------------------------

  SetAngles();						// Set need angle fcts.

// ----------------------------------------------------------------------------
//               Set Up Jacobian Matrix and Functional Vector
// ----------------------------------------------------------------------------

  double f, dfda, dfdb, dfdg;
  J = matrix(3,3,complex0);				// Jacobian Matrix
  F = col_vector(3,complex0);				// Functional Vector
  Fuse = "";						// String of functions
  double  delzzo4 =  0.25*PASDelzz;			// Scaling Factor
  bool sing = false;					// Flag for singular J
  string suv;
  int i,j;
  for(i=0; i<3; i++)					// Loop Cartesian elems
    {
    switch(Flist[i])					// Fill Jacobian Row &
      {							// Functional for Auv
      default: 						// that is choosen
      case 0: suv="xx"; Axx(f,dfda,dfdb,dfdg); break;	// USE Axx
      case 1: suv="xy"; Axy(f,dfda,dfdb,dfdg); break;	// USE Axy
      case 2: suv="xz"; Axz(f,dfda,dfdb,dfdg); break;	// USE Axz
      case 3: suv="yy"; Ayy(f,dfda,dfdb,dfdg); break;	// USE Ayy
      case 4: suv="yz"; Ayz(f,dfda,dfdb,dfdg); break;	// USE Ayz
      case 5: suv="zz"; Azz(f,dfda,dfdb,dfdg); break;	// USE Azz
      }

    f = delzzo4*f - Av[Flist[i]];			// Form function --> 0
    dfda *= delzzo4;					// Form Jacobian <i|J|0>
    dfdb *= delzzo4;					// Form Jacobian <i|J|1>
    dfdg *= delzzo4;					// Form Jacobian <i|J|2>
    if(fabs(f)    < DifCut) f    = 0.0;			// Insure small f is 0
    if(fabs(dfda) < DifCut) dfda = 0.0;			// Insure <i|J|0> is 0
    if(fabs(dfdb) < DifCut) dfdb = 0.0;			// Insure <i|J|1> is 0
    if(fabs(dfdg) < DifCut) dfdg = 0.0;			// Insure <i|J|2> is 0

    TrackMin(49, Flist[i]);				// J row OK. Track
    if(Plevel>3)
      {
      cout << "\n\t\tAlpha Is:       " << EA.alpha()*RAD2DEG;
      cout << "\n\t\tBeta Is:        " << EA.beta()*RAD2DEG;
      cout << "\n\t\tGamma Is:       " << EA.gamma()*RAD2DEG;
      cout << "\n\t\tPAS delzz Is:   " << PASDelzz;
      cout << "\n\t\tAsymmetry Is:   " << EtaVal;
      cout << "\n\t\tA"  << suv << " Is:         " << Av[i];
      cout << "\n\t\tf"  << suv << " Is:         " << f;
      cout << "\n\t\t@f" << suv << "/@alpha Is: "  << dfda;
      cout << "\n\t\t@f" << suv << "/@beta Is:  "  << dfdb;
      cout << "\n\t\t@f" << suv << "/@gamma Is: "  << dfdg;
      }

    if(!dfda && !dfdb && !dfdg)			// Reject all zero row
      { TrackMin(50); return false; }		//   Track, return fail
    sing = false;				// Flag for singular J
    for(j=0; j<i && !sing; j++)			// Check for row making
      { 					// a singular matrix
      sing = CheckSing(dfda, dfdb, dfdg, j); 	//   If array J singular
      if(sing) { TrackMin(51); return false; }	//   Track, return fail
      }
    TrackMin(52, Flist[i]);			// J row OK. Track
    Fuse += Gdec(Flist[i]);			// Asjust function string
    F.put(f,i);					// Fill Function row i
    J.put(dfda,i,0); 				// Fill Jacobian row i
    J.put(dfdb,i,1);
    J.put(dfdg,i,2);
    }						// Next Jacobian/Fct row
  return true;					// True, found all 3 rows
  }


// sosi
//Plevel
void CartMx2A::TrackMin(int info, double val) const
  {
  string spc(1, ' ');
  if(!Plevel) return;
  if(Plevel)
    {
    switch(info)
      {
      case  9: 
        cout << "\n\t* Minimization Step " << Step+1;
        break;
      case 10:
        cout << "\n\t    Minimization Failed To Converge";
        cout << "\n\t    The Norm Of |Y> Is "          << Y.norm();
        cout << "\n\t    Greater Than Norm Cutoff Of " << NormCut;
        cout << "\n\t    Greater Than Last Norm Of "   << val;
        break;
      case 11:
        cout << "\n\t    Minimizing Slowly?";
        cout << "\n\t    The Norm Of |Y> Is "          << Y.norm();
        cout << "\n\t    Greater Than Norm Cutoff Of " << NormCut;
        cout << "\n\t    Less Than Last Norm Of "      << val;
        cout << "\n\t    Same Step Repeat Count Of "   << SStep;
        break;
      case 12:
        cout << "\n\t    ==> Minimization Not Converging";
        cout << "\n\t    Trying With Other Auv Functions";
        break;
      case 13:
        cout << "\n\t    Less Than Max Allowed Step Repeats Of "   << MaxSSteps;
        cout << "\n\t    ==> Minimization Converging Slowly?";
        cout << "\n\t    Trying More Iteration Steps";
        break;
      case 14:
        cout << "\n\t    Surpassed Max Allowed Step Repeats Of "   << MaxSSteps;
        cout << "\n\t    ==> Minimization Not Converging";
        break;
      case 20:
        cout << "\n\t    Minimization Norm Has Been Met OK";
        cout << "\n\t    The Norm Of |Y> Is " << Y.norm();
        cout << "\n\t    Less Than Norm Cutoff Of " << NormCut;
        cout << "\n\t    ==> Array A Not Consistent With Minimized A";
        cout << "\n\t    ==> Difference A-A' Deemed To Be " << val;
        cout << "\n\t    ==> Difference A-A' Greater Than " << ZeroMxCut;
        cout << "\n\t    ==> This Minimum Cannot Be The One We Seek";
        break;
      case 21:
        cout << "\n\t    Trying New Starting Angles {"
             << EAo.alpha()*RAD2DEG << ", " 
             << EAo.beta()*RAD2DEG  << ", " 
             << EAo.gamma()*RAD2DEG << "}"; 
        break;
      case 22:
        cout << "\n\t    Trying New Jacobian Function Set";
        break;
      case 30: cout << "\n\t    ==> Cannot Build Valid Jacobian";
        break;
      case 31: cout << "\n\t    ==> Minimization Not Converging";
        break;
      case 47:
        cout << "\n\t    Previous Functions " 
             << JFunctName(Flist[0]) << ", "
             << JFunctName(Flist[1]) << ", "
             << JFunctName(Flist[2]);
        break;
      case 48:
        cout << "\n\t    New Functions "
             << JFunctName(Flist[0]) << ", "
             << JFunctName(Flist[1]) << ", "
             << JFunctName(Flist[2]);
        break;
      }
    }
  if(Plevel > 1)
    {
    switch(info)
      {
      case 0: cout << "\n\n" << spc << "                   Cartesian Tensor Euler Angle Fitting\n";
              cout << "\n" << spc << "             Starting                 Correction                 Final          Fct |Y> Norm";
              cout << "\n" << spc << "Step    Alpha   Beta  Gamma      Alpha   Beta  Gamma      Alpha   Beta  Gamma   Use";
              cout << "\n" << spc << "====  =======================  =======================  ======================= === ========";
              break;
      case 1: cout << "\n" << spc << Gdec(Iter,4) << "  " 
                           << Gform("%7.3f", EA.alpha() * RAD2DEG) << " "
                           << Gform("%7.3f", EA.beta()  * RAD2DEG) << " "
                           << Gform("%7.3f", EA.gamma() * RAD2DEG) << " ";
              break;
      case 2: cout << Gform("%8.3f", fmod(Y.getRe(0) * RAD2DEG,360.0)) << ""
                   << Gform("%8.3f", fmod(Y.getRe(1) * RAD2DEG,360.0)) << ""
                   << Gform("%8.3f", fmod(Y.getRe(2) * RAD2DEG,360.0)) << "  "
                   << Gform("%7.3f", EA.alpha() * RAD2DEG) << " "
                   << Gform("%7.3f", EA.beta()  * RAD2DEG) << " "
                   << Gform("%7.3f", EA.gamma()  * RAD2DEG) << " " << Fuse << " ";
              if(Y.norm() > 1000) cout << "Real Big";
              else                cout << Y.norm();
              break;
      case 3: cout << "\n";
              break;
      default: break;
      }
    }
  if(Plevel > 2)
    {
    switch(info)
      {
      default:
        {
        matrix Jinv = inv(J);
        cout << "\n\n\tJacobian Matrix\n"   << J;
        cout << "\n\n\tFunctional Vector\n" << F;
        cout << "\n\n\tJ Inverse\n"         << Jinv;
        cout << "\n\n\tJ*inv(J)\n"          << J*Jinv;
        cout << "\n\n\tCorrection Y\n"      << Y;
        cout << "\n\n\tNext Update\n"       << X+Y; 
        cout << "\n\t\t\tNorm |Y> = "       << Y.norm();
        }
        break;
      case 49:
        cout << "\n\n\tWorking With "   << JFunctName(int(val)) << ":";
        break;
      case 50:
        cout << "\n\t\tJacobian Row Rejected, All Zero";
        break;
      case 51:
        cout << "\n\t\tJacobian Row Rejected, Singular";
        break;
      case 52:
        cout << "\n\t\t" << JFunctName((int)val)
             << " Values Accepted As Jacobian Row";
        break;
      }
    }
// sosi
  return;
  }

// ____________________________________________________________________________
// F                  Class CartMx2A Auxiliary Functions
// ____________________________________________________________________________

matrix CartMx2A::Regenerate() const
  {
  IntRank2A A(EtaVal);
  matrix mx(3,3,complex0,h_matrix_type);

/*
IR2ACart CC = A.CartCmp(EA);
mx.put(CC.Axx(),0,0); mx.put_h(CC.Axy(),0,1); mx.put_h(CC.Axz(),0,2);
                        mx.put(CC.Ayy(),  1,1); mx.put_h(CC.Ayz(),1,2);
                                                  mx.put(CC.Azz(),  2,2); 
*/
  mx.put(A.Axx(EA),0,0); mx.put_h(A.Axy(EA),0,1); mx.put_h(A.Axz(EA),0,2);
                         mx.put(A.Ayy(EA),  1,1); mx.put_h(A.Ayz(EA),1,2);
                                                  mx.put(A.Azz(EA),  2,2); 

  mx *= PASDelzz/A.delzz();
  mx += matrix(3,3,AisoVal,h_matrix_type);
  return mx;
  }

double CartMx2A::Check() const
  {
  matrix ACart = Regenerate();
  matrix Adiff = ACart - A;
  double nd = 0;
  if(!Adiff.is_zero())
    {
    int i, j;
    for(i=0; i<3; i++)
      for(j=0; j<3; j++)
        nd += fabs(Adiff.getRe(i,j));
    }
  return nd;
  }

bool CartMx2A::Check(double cutoff) const { return (Check() < cutoff); }

string CartMx2A::Type() const
  {
  if(A.is_zero())     return string("Zero");
  if(!PASDelzz)       return string("Isotropic");	
  if(!EtaVal)         return string("Symmetric");
                      return string("Asymmetric");
  if(A.is_diagonal()) return string("PAS");
  }

string CartMx2A::Method() const
  {
  string cv;
  switch(ConvMeth)
    {
    default:
    case 0: cv = string("None Required");                       break; 
    case 1: cv = string("Symmetric Diagonal");                  break; 
    case 2: cv = string("Symmetric NonDiagonal");               break; 
    case 3: cv = string("Asymmetric Diagonal");                 break; 
    case 4: cv = string("Asymmetric NonDiagonal With Zeros");   break; 
    case 5: cv = string("Asymmetric NonDiagonal Symmetric OD"); break; 
    case 6: cv = string("Asymmetric NonDiagonal");              break; 
    }
  return cv;
  }

string CartMx2A::JFunctName(int i) const
  {
  string fn;
  switch(i)
    {
    default: fn = string("No Idea"); break;			
    case 0:  fn = string("Axx");     break; 
    case 1:  fn = string("Axy");     break; 
    case 2:  fn = string("Axz");     break; 
    case 3:  fn = string("Ayy");     break; 
    case 4:  fn = string("Ayz");     break; 
    case 5:  fn = string("Azz");     break; 
    }
  return fn;
  }

void CartMx2A::TrackConv(int info) const
  {
  string spc(20, ' ');
  if(!Plevel) return;
  if(Plevel >= 1)
    {
    switch(info)
      {
      default:
      case  0: cout << "\n" << spc << "Tracking Cartesian Tensor Conversion";
               cout << "\n" << spc << "------------------------------------";
               cout << "\n" << spc << "------------------------------------\n";
               break;
      case  1: cout << "\n\t" << "* Treating Diagonal Symmetric A";
               break;
      case  2: cout << "\n\t" << "* Treating Non-Diagonal Symmetric A";
               break;
      case  3: cout << "\n\t" << "* Treating Diagonal Asymmetric A";
               break;
      case  4: cout << "\n\t" << "* Treating Non-Diagonal Asymmetric A With Zeros";
               break;
      case  5: cout << "\n\t" << "* Treating Non-Diagonal Asymmetric A With Special Off Diagonals";
               break;
      case  6: cout << "\n\t" << "* Treating Non-Diagonal Asymmetric A";
               break;
      case  7: cout << "\n\t" << "    Both Axy & Axz Are 0, Ayz Non-Zero!";
               break;
      case  8: cout << "\n\t" << "    Both Axy & Ayz Are 0, Axz Non-Zero!";
               break;
      case  9: cout << "\n\t" << "    Both Axz & Ayz Are 0, Axy Non-Zero!";
               break;
      case 10: cout << "\n\t" << "    No Non-Zero Off-Diagonal Elements Found";
               break;
      case 11: cout << "\n\t" << "    Preferred Functions: "
                              << JFunctName(Flist[0]) << ", "
                              << JFunctName(Flist[1]) << ", "
                              << JFunctName(Flist[2])
                    << "\n\t" << "    Starting Angles: "
                              << EAo.alpha()*RAD2DEG << ", "
                              << EAo.beta() *RAD2DEG << ", "
                              << EAo.gamma()*RAD2DEG;
               break;
      case 12: cout << "\n\t" << "    Minimized To Within Norm Of " << NormCut
                    << "\n\t" << "    A Is Self-Consistent To Within " << ZeroMxCut
		    << "\n\t" << "* Minimization Completed\n\n"; 
               break;
      case 16: cout << "\n\t" << "* Failed Using This Conversion Method!\n\n";
               break;
      }
    }
  }

bool CartMx2A::OffDiagonals() const
  {
  double axy = A.getRe(0,1) - AisoVal;			// Get irreducible
  double axz = A.getRe(0,2) - AisoVal;			// rank 2 components
  double ayz = A.getRe(1,2) - AisoVal;			// to work with
  if((fabs(axy)<DifCut) 
  || (fabs(axz)<DifCut)
  || (fabs(ayz)<DifCut))
    return true;
  return false;
  }

bool CartMx2A::SymOffDiags() const
  {
  double axy = fabs(A.getRe(0,1) - AisoVal);		// Get irreducible
  double axz = fabs(A.getRe(0,2) - AisoVal);		// rank 2 components
  double ayz = fabs(A.getRe(1,2) - AisoVal);		// to work with
  if(fabs(axy-axz) < DifCut) return true; 
  if(fabs(axy-ayz) < DifCut) return true; 
  if(fabs(axz-ayz) < DifCut) return true; 
  return false;
  }

// ____________________________________________________________________________
// Z               Class CartMx2A Formatted Output Functions
// ____________________________________________________________________________

/* These functions will output information concerning the conversion of the
   Cartesian spatial tensor into a GAMMA spatial tensor. The information will
   be output into any specified output stream.

                Input           C2A     : Cartesian 2 GAMMA A Conversion (this)
                                ostr    : Output stream
                                fflag   : Format flag
                                            0 - Sum rotation only
                                           !0 - Full composite rotation
           Output               none    : Conversion information placed
                                          into the output stream            */


ostream& CartMx2A::print(ostream& ostr, int fflag) const
  {
  string hdr("Rank 2 Cartesian Spatial Tensor Conversion");
  string spcr = string(2, ' ');
  string splt = string(5, ' ');
  string numl = string(8, ' ');
  string nl = string("\n") + spcr;
  ostr << "\n\n" << string((80-hdr.length())/2, ' ') << hdr << "\n";
  
  ostr << nl   << " GAMMA Components "
       << splt << string(18, ' ') << "Cartesian Array";
  ostr << nl   << "------------------"
       << splt << string(18, '-') << "---------------" << string(17, '-');
  ostr << nl << "Aiso:  " << Gform("%7.3f", AisoVal)            << "    "
       << splt << "[ A    A    A   ]   [ "
       << Gform("%8.3f", A.getRe(0,0)) << " "
       << Gform("%8.3f", A.getRe(0,1)) << " "
       << Gform("%8.3f", A.getRe(0,2))
       << " ]";
  ostr << nl << "delzz: " << Gform("%7.3f", PASDelzz)           << "    "
       << splt << "|  xx   xy   xz |   | "
       << numl                        << " "
       << numl                        << " "
       << numl
       << " |";
  ostr << nl << "eta:   " << Gform("%7.3f", EtaVal)             << "    "
       << splt << "| A    A    A   | = | "
       << Gform("%8.3f", A.getRe(1,0)) << " "
       << Gform("%8.3f", A.getRe(1,1)) << " "
       << Gform("%8.3f", A.getRe(1,2))
       << " |";
  ostr << nl << "alpha: " << Gform("%7.3f", EA.alpha()*RAD2DEG) << " deg"
       << splt << "|  yx   yy   yz |   | "
       << numl                        << " "
       << numl                        << " "
       << numl
       << " |";
  ostr << nl << "beta:  " << Gform("%7.3f", EA.beta()*RAD2DEG)  << " deg"
       << splt << "| A    A    A   |   | "
       << Gform("%8.3f", A.getRe(2,0)) << " "
       << Gform("%8.3f", A.getRe(2,1)) << " "
       << Gform("%8.3f", A.getRe(2,2))
       << " |";
  ostr << nl << "gamma: " << Gform("%7.3f", EA.gamma()*RAD2DEG) << " deg"
       << splt << "[  zx   zy   zz ]   [ "
       << numl                        << " "
       << numl                        << " "
       << numl
       << " ]";
  if(fflag)
    {
    hdr = string("Minimization Parameters");
    ostr << "\n\n" << string((80-hdr.length())/2, ' ') 
         << hdr;
    ostr << "\n"   << string((80-hdr.length())/2, ' ') 
         << string(hdr.length(), '-') << "\n";

    ostr << "\n\tMaximum No. Iterations: " << MaxIter;
    ostr << "\n\tCompleted Iterations:   " << Iter;
    ostr << "\n\tEta Zero Cutoff:        " << EtaCut;
    ostr << "\n\tElement Zero Cutoff:    " << DifCut;
    ostr << "\n\t|Y> Norm Zero Cutoff:   " << NormCut;
    ostr << "\n\tInital Alpha Angle:     " << EAo.alpha()*RAD2DEG;
    ostr << "\n\tInital Beta Angle:      " << EAo.beta()*RAD2DEG;
    ostr << "\n\tInital Gamma Angle:     " << EAo.gamma()*RAD2DEG;
    ostr << "\n\tConversion Method:      " << Method();
    }
  ostr << "\n\n";
  return ostr;
  }

ostream& operator<< (ostream& ostr, const CartMx2A& C2A)
  { return C2A.print(ostr); }

#endif							// CartMx2A.cc
