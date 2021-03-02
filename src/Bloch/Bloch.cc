/* Bloch.cc *****************************************************-*-c++-*-
**									**
**                                 G A M M A				**
**      Bloch Equation Auxiliaryi                      Implementation	**
**									**
**      Copyright (c) 1995, 2002					**
**      S.A. Smith 							**
**      National High Magnetic Field Laboratory				**
**      1800 E. Paul Dirac Drive					**
**      Tallahassee Florida, 32306-4005					**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
**  Description 							**
**									**
**  This module contains functions which aid in simulations involving 	**
**  the phenomenological Bloch equations. Most of the Bloch module	**
**  functionality resides in files assocaited with Bloch related 	**
**  classes. Thus, the functions herein are largely cosmetic and/or	**
**  deprecated.								**
**									**
*************************************************************************/

#ifndef   GBloch_cc_ 			// Is file already included?
#  define GBloch_cc_ 1      		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation   
#  endif

#include <GamGen.h>			// Inlcude system specifics
#include <Bloch/Bloch.h>		// Include our interface
#include <Matrix/matrix.h>		// Include matrices
#include <Basics/Gutils.h>		// Include query
#include <Basics/Gconstants.h>		// Include PI
#include <Basics/ParamSet.h>		// Include Parameter Sets
#include <Basics/StringCut.h>		// Include Gdec function

// ____________________________________________________________________________
// A                  Interactive I/O Auxiliary Functions
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                      Output Magnitude Functions
// ----------------------------------------------------------------------------

int DoubleMag(double x)
  {
  if(fabs(x) < 1.e-13) return 0;
  int mag=0;
  x = fabs(x);
  while(x >= 10) { mag++; x/=10; }
  while(x <  1 ) { mag--; x*=10; }
  return mag;
  }

// ----------------------------------------------------------------------------
//                        Output Unit Functions
// ----------------------------------------------------------------------------


std::string SecUnits(int mag, double& sf)
  {
  std::string timeu;
  switch(mag)
    {
    case  8:  case  9:  case 10:  sf = 1.e-6; return std::string("Msec"); break;
    case  5:  case  4:  case  3:  sf = 1.e-3; return std::string("Ksec"); break;
    case  2:  case  1:  case  0:  sf = 1.0;   return std::string("sec");  break;
    case -1:  case -2:  case -3:  sf = 1.e3;  return std::string("msec"); break;
    case -4:  case -5:  case -6:  sf = 1.e6;  return std::string("usec"); break;
    case -7:  case -8:  case -9:  sf = 1.e9;  return std::string("nsec"); break;
    case -10: case -11: case -12: sf = 1.e12; return std::string("psec"); break;
    case -13: case -14: case -15: sf = 1.e15; return std::string("fsec"); break;
    default:                      sf = 1.0;   return std::string("sec"); break;
    }
  return timeu;
  }

std::string HzUnits(int mag, double& sf)
  {
  std::string Hzu;
  switch(mag)
    {
    case  8:  case  9:  case 10:  sf = 1.e-6; return std::string("MHz"); break;
    case  5:  case  4:  case  3:  sf = 1.e-3; return std::string("KHz"); break;
    case  2:  case  1:  case  0:  sf = 1.0;   return std::string("Hz");  break;
    case -1:  case -2:  case -3:  sf = 1.e3;  return std::string("mHz"); break;
    case -4:  case -5:  case -6:  sf = 1.e6;  return std::string("uHz"); break;
    case -7:  case -8:  case -9:  sf = 1.e9;  return std::string("nHz"); break;
    case -10: case -11: case -12: sf = 1.e12; return std::string("pHz"); break;
    case -13: case -14: case -15: sf = 1.e15; return std::string("fHz"); break;
    default:                      sf = 1.0;   return std::string("Hz"); break;
    }
  return Hzu;
  }



// ----------------------------------------------------------------------------
//             Initial Magnetization Vector Functions
// ----------------------------------------------------------------------------




	// Input	    Mox	: Initial magnetization x-component
	//      	    Moy	: Initial magnetization y-component
	//      	    Moz	: Initial magnetization z-component
	// Output	    Mo	: 3x1 column vector (unitless),
	//			  the initial magnetization vector
	// Note			: The return is normalized.
	// Note			: The default is Meq

matrix Mo_vector(double Mox, double Moy, double Moz)
 {
 matrix Mo(3, 1, complex0);
 Mo.put(Mox, 0, 0);
 Mo.put(Moy, 1, 0);
 Mo.put(Moz, 2, 0);
 double Mmag = Mox*Mox+Moy*Moy+Moz*Moz;
 Mo *= complex(1/sqrt(Mmag));
 return Mo;
 }



// Set Up The Initial Magnetization Vector
matrix Mo_vector(int argc, char* argv[], matrix& Meq, int& qn)
 {
 matrix Mo(3, 1, complex0);
 Mo = Meq;					// Set Mo to equilibrium
 std::string yn;
 double Mcomp, Mmag;
 query_parameter(argc, argv, qn,
     "\n\t* Begin With a Non-equilibrium Magnetization Vector [y,n]? ", yn);
 qn++;
 if(yn == "y")
   {
   std::cout << "\n\n\t\t* X Component? ";		// Get the new X component
   std::cin >> Mcomp;
   Mo.put(Mcomp, 0, 0);
   Mmag = Mcomp*Mcomp;
   std::cout << "\n\n\t\t* Y Component? ";		// Get the new Y component
   std::cin >> Mcomp;
   Mo.put(Mcomp, 1, 0);
   Mmag += Mcomp*Mcomp;
   std::cout << "\n\n\t\t* Z Component? ";		// Get the new Z component
   std::cin >> Mcomp;
   Mo.put(Mcomp, 2, 0);
   Mmag += Mcomp*Mcomp;
   Mmag = sqrt(Mmag);
   Mo *= complex(1/Mmag);			// Scale to magnitude 1
   }
 return Mo;
 }


// ----------------------------------------------------------------------------
//                      Analysis Functions
// ----------------------------------------------------------------------------

 void analyze(double tinc, int& ntimes,
          int& do_ss, int& qn, double T1, double gamB1, double w)

// sosi - note do_ss and qn aren't yet used.  Is this unfinished?
//	 	 Set Up Time Parameters

 {
int itmp = do_ss;
itmp = qn;
 double time;
 if(T1)
   {
   time = 5*T1/(ntimes-1);
   std::cout <<"\n";
   std::cout << "\n\tFor The Specified T1 and Number of Trajectory Points";
   std::cout << "\n\tA Minimum Time Increment of " << time << " Seconds"; 
   std::cout << "\n\tIs Required to Span 5 T1 (to equilibrium or steady-state";
   }
 if(gamB1)
   {
   time =1/gamB1;
   std::cout <<"\n";
   std::cout << "\n\tTo Rotate 1 Time Around B1 Axis Takes " 
        << time << " Seconds"; 
   }
 time = 2*PI/w;
 std::cout << "\n\tTo Rotate 1 Time around B0 (z) Axis in the"
      << " Rotating Frame Takes " << time << " Seconds"; 
 if(gamB1)
   {
   std::cout << "\n\tThe Trajectory Should Circle Somewhere Between "
        << int(double(ntimes-1)*tinc/time) << " and "
        << int(double(ntimes-1)*tinc*gamB1) << " Times";
   }
 else
   {
   std::cout << "\n\tThe Trajectory Should Circle About "
        << int(double(ntimes-1)*tinc/time) << " Times";
   }
 return;
 }

// ----------------------------------------------------------------------------
//                    Bloch Related Parameters
// ----------------------------------------------------------------------------

void bloch_T1T2(const ParameterSet& pset, std::ostream& ostr,
                                                        double& T1, double& T2)
  {
  ParameterSet::const_iterator item;		// A pix into parameter list
  std::string pname, ssfile, pstate;		// Items in each pset entry
  pname = std::string("T1");			// Longitudinal relaxation time
  item = pset.seek(pname);			// Parameter for T1
  if(item != pset.end())			// If T1 is a valid parameter
    {
    (*item).parse(pname,T1,pstate);		// Read in the T1 value
    pname = std::string("T2");			// Transverse relaxation time
    item = pset.seek(pname);			// Parameter for T2
    if(item != pset.end())			// If T2 is a valid parameter
      (*item).parse(pname,T2,pstate);		// Read in the T2 value
     else
       {
       ostr << "\n\tCant Find T2 Time."
            << "\n\tSetting Both T1 and T2 To Zero.\n";
       T1 = 0;
       T2 = 0;
       }
     }
   else
     {
     ostr << "\n\tCant Find T1 Time."
          << "\n\tSetting Both T1 and T2 To Zero.\n";
     T1 = 0;
     T2 = 0;
     }
   }


 void bloch_Mo(const ParameterSet& pset, std::ostream& ostr, double& Mx, double& My, double& Mz)

  {
  ParameterSet::const_iterator item; // A pix into parameter list
  std::string pname, ssfile, pstate;		// Items in each pset entry
  pname = std::string("Mx");			// Magnetization x-component
  item=pset.seek(pname);		// Parameter for Mx
  if(item != pset.end())		// If Mx is a valid parameter
     (*item).parse(pname,Mx,pstate);	// Read in the Mx value
   else
    {
     ostr << "\n\tCant Find Mx Value,"
           << " Setting Mx to Zero.";
     Mx = 0;
     }
   pname = std::string("My");		// Magnetization y-component
   item = pset.seek(pname);		// Parameter for My
   if(item != pset.end())		// If My is a valid parameter
     (*item).parse(pname,My,pstate);	// Read in the My value
   else
     {
     ostr << "\n\tCant Find My Value,"
           << " Setting My to Zero.";
     My = 0;
     }
   pname = std::string("Mz");		// Magnetization z-component
   item = pset.seek(pname);		// Parameter for Mz
   if(item != pset.end())		// If Mz is a valid parameter
     (*item).parse(pname,Mz,pstate);	// Read in the Mz value
   else
     {
     ostr << "\n\tCant Find Mz Value,"
           << " Setting Mz to Zero.";
     Mz = 0;
     }
   if(Mx==0 && My==0 && Mz==0)
     {
     ostr << "\n\tCant Find Any Initial Magnetization Values,"
           << " Setting Mx=My=0, Mz=1.";
     Mz = 1.0;
     }
   }

 void bloch_B1(const ParameterSet& pset, std::ostream& ostr, double& gamB1, double& phi)

  {
  ParameterSet::const_iterator item;	// A pix into parameter list
  std::string pname, ssfile, pstate;		// Items in each pset entry
  pname = std::string("gamB1");		// RF-field strength
  item = pset.seek(pname);		// Parameter for gamB1
   if(item != pset.end())		// If gamB1 is a valid parameter
     {
     (*item).parse(pname,gamB1,pstate);	// Read in the gamB1 value
     pname = std::string("phi");		// RF-field phase
     item = pset.seek(pname);		// Parameter for phi
     if(item != pset.end())		// If phi is a valid parameter
       (*item).parse(pname,phi,pstate);	// Read in the phase value
     else
       {
       ostr << "\n\tCant Find RF-Field Phase."
            << "\n\tSetting The RF-Field Phase To Zero.\n";
       phi = 0;
       }
     }
   else
     {
     ostr << "\n\tCant Find RF-Field Strength."
          << "\n\tSetting Both gamB1 and phi To Zero.\n";
     gamB1 = 0;
     phi = 0;
     }
   }


 void bloch_Woff(const ParameterSet& pset, std::ostream& ostr, double& Woff)

  {
  ParameterSet::const_iterator item;         // A pix into parameter list
  SinglePar par;
  std::string pname, ssfile, pstate;		// Items in each pset entry
  pname = std::string("Woff");			// Offset from RF-field
  item = pset.seek(pname);			// Parameter for phi
  if(item != pset.end())			// If Woff is a valid parameter
    (*item).parse(pname,Woff,pstate);		// Read in the offset value
  else
    {
    ostr << "\n\tCant Find Offset From RF-Field."
         << "\n\tSetting Woff To Zero.\n";
    Woff = 0;
    }
  }



// ----------------------------------------------------------------------------
//                      Trajectory Timing Functions
// ----------------------------------------------------------------------------

void T1Blurb(std::ostream& ostr, double T1, int N)
  {
  if(T1)
    {
    double tinc = 5*T1/(N-1);
    ostr << "\n\tThe Longitudinal Relaxation Time Is  " << T1;
    ostr << "\n\tThe Number Of Points To Sample Is    " << N;
    ostr << "\n\tMinimum Time Increment For 5 T1's Is " << tinc << " Sec";
    }
  }

void B1Blurb(std::ostream& ostr, double gamB1)
  {
  if(gamB1)
    {
    double tinc =1/gamB1;
    ostr << "\n\tThe Applied Field Strength Is        " << gamB1 << " Hz";
    ostr << "\n\tOne Rotatation About B1 Axis Takes   " << tinc  << " Sec"; 
    }
   }

void WBlurb(std::ostream& ostr, double W)
  {
  if(W)
    {
    double tinc = PIx2/W;
    ostr << "\n\tThe Vector Offset Frequency Is       " << W     << " Hz";
    ostr << "\n\tOne Rotatation About B0 Axis Takes   " << tinc  << " Sec"; 
    }
  }

void TrajTiming(int argc, char* argv[], 
     double& tinc, int& ntimes, int& qn, double T1, double gamB1, double w)

  {
  std::string msg = std::string("\n\tNumber of Time Increments[") + Gdec(ntimes) + "]? ";
  ask_set(argc, argv, qn++, msg, ntimes);
  T1Blurb(std::cout, T1, ntimes);
  B1Blurb(std::cout, gamB1);
  WBlurb(std::cout,  w);
  msg = "\n\tTime Increment Between Magnetization Sampling (sec) ["
      + Gform("%8.3f", tinc) + "]? ";
  ask_set(argc, argv, qn, msg, tinc);
  } 

#endif								// Bloch.cc

