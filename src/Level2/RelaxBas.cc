/* RBasic.cc ****************************************************-*-c++-*-
**									**
**	                           G A M M A				**
**									**
**	Basic Relaxation                            Implementation	**
**									**
**	Scott A. Smith							**
**	Copyright (c) 1999						**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
** Description								**
**									**
** This module deals with a very primitive relaxation treatement. In 	**
** this case the users ascribes a single T1 and T2 value to each spin	**
** in a system.  The class RBasic provides the means of reading such	**
** values from an external file, generating a pseudo-relaxation super-	**
** operator, and evolving a density operator under simple relaxation.	**
**									**
** In this treatment, the off-diagonal density operator elements in the	**
** product basis (single and multiple quantum coherences) will decay	**
** with a rate of 1/T2 = R2 and the Fz components of each spin, Izi,    **
** will relax at a rate 1/T1 = R1.                                      **
**									**
*************************************************************************/

///Chapter Basic Relaxation
///Section Overview
///Body The ...
///Section Available Basic Relaxation Functions

#ifndef   RBasic_cc_			// Is this file already included?
#  define RBasic_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <stdlib.h>
#include <Level2/RelaxBas.h>		// Include interface
#include <vector>			// Include libstdc++ STL vectors
#include <Basics/Gutils.h>		// Include GAMMA error messages
#include <Basics/StringCut.h>		// Include GAMMA Gdec function
#include <Basics/Gconstants.h>		// Include GAMMA constants (PI)
#include <Basics/ParamSet.h>            // Include parameter sets
#include <Matrix/complex.h>		// Include GAMMA complex numbers
#include <Matrix/row_vector.h>		// Include GAMMA row_vectors
#include <HSLib/SpinOpCmp.h>		// Include single spin operators
#include <HSLib/SpinOp.h>		// Include spin operators
#include <HSLib/SpinSys.h>		// Include isotropic spin systems
#include <HSLib/HSprop.h>		// Include Hilbert space evolution
#include <HSLib/HSauxil.h>		// Include density operator stuff
#include <LSLib/SuperOp.h>		// Include superoperators
#include <LSLib/LSacquire.h>		// Include Liouville acquisitions
#include <Level2/acquire1D.h>

static double DefR12 = 0.0; 		// Default R1,R2 value
static double DefT12 = 1.e6;            // Default T1,T2 value
static double DefLW  = 1.e-4;           // Default LW    value

// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                        Relax Basic Error Handling
// ____________________________________________________________________________

/* These functions take care of any errors encountered when reading, writing,   
   and setting parameters in Bruker acquisition parameter files.
 
	Input		RBas    : Basic relaxation (this)
			eidx    : Error index
                        pname   : Additional error message
                        noret   : Flag for linefeed (0=linefeed)
        Output          void    : An error message is output                 */


/* The following error messages use the defaults set in the Gutils package
 
                Case                          Error Message
 
                (0)                     Program Aborting.....
                (3)                     Can't Construct From Parameter Set
                (4)                     Cannot Construct From Input File
                (5)                     Cannot Write To Parameter File
                (6)                     Cannot Write To Output FileStream
                (7)                     Cannot Produce Any Output
                default                 Unknown Error                        */
/* The following error messages use the defaults set in the Gutils package
 
                Case                          Error Message
 
                (1)                     Problems With File PNAME
                (2)                     Cannot Read Parameter PNAME
                default                 Unknown Error - PNAME                */  
 

void RBasic::RBasErr(int   eidx, int noret) const
  {
  std::string hdr("Basic Relaxation");
  std::string msg;
  switch (eidx)
    {
    case 40: GAMMAerror(hdr,"Negative Relaxation Rate!",        noret);
             break;                              //                        (40)
    case 41: GAMMAerror(hdr,"Non-Positive Relaxation Time!",    noret);
             break;                              //                        (41)
    case 42: GAMMAerror(hdr,"Negative Line Width Disallowed",   noret);
             break;                              //                        (42)
    case 50: GAMMAerror(hdr,"Accessed Spin Index Out of Range", noret);
             break;                              //                        (50)
    case 59: GAMMAerror(hdr,"Accessed Spin Out Of Range!",      noret);
             break;                              //                        (59)
    case 60: GAMMAerror(hdr,"Spin Size MisMatch!",              noret);
             break;                              //                        (60)
    case 61: GAMMAerror(hdr,"Incompatible Spin System",         noret);
             break;                              //                        (61)
    case 62: GAMMAerror(hdr,"Cannot Generate Iz Coefficients",  noret);
             break;                              //                        (62)
    case 63: GAMMAerror(hdr,"Cannot Set Up Iz Spin Operators",  noret);
             break;                              //                        (63)
    case 64: GAMMAerror(hdr,"Cannot Perform System Evolution",  noret);
             break;                              //                        (64)
    case 65: GAMMAerror(hdr,"Cannot Perform Acquisition",       noret);
             break;                              //                        (65)
    case 70: GAMMAerror(hdr,"Hilbert Space Dimension MisMatch!",noret);
             break;                              //                        (70)
    case 71: GAMMAerror(hdr,"Incompatible Density Operator",    noret);
             break;                              //                        (71)
    case 72: GAMMAerror(hdr,"Cannot Make Iz Inf. Coefficients", noret);
             break;                              //                        (72)
    case 80: GAMMAerror(hdr,"Problems Setting The Hamiltonian", noret);
             break;                              //                        (80)
    case 81: GAMMAerror(hdr,"Problems Setting Detection Op.",   noret);
             break;                              //                        (81)
    case 82: GAMMAerror(hdr,"Must First Specify A Spin System", noret);
             break;                              //                        (82)
    case 83: GAMMAerror(hdr,"Incompatible General Operator",    noret);
             break;                              //                        (83)
    case 85: GAMMAerror(hdr,"Problems Setting Sigma Infinte",   noret);
             break;                              //                        (85)
    case 90: GAMMAerror(hdr,"Problems Setting Spin System", noret);
             break;                              //                        (90)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  if(!noret) std::cout << "\n";
  }  

void RBasic::RBasErr(int   eidx, const std::string& pname, int noret) const
  {
  std::string hdr("Basic Relaxation");
  std::string msg;
  switch(eidx)
    {
    case 50: msg = std::string("Cannot Get ");
             GAMMAerror(hdr,msg+pname,noret);  break;   //                 (50)
    case 51: msg = std::string("Cannot Set ");
             GAMMAerror(hdr,msg+pname,noret);  break;   //                 (51)
    default: GAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error   (-1)
    }
  if(!noret) std::cout << "\n";
  }

void RBasic::RBasFatal(int eidx) const
  {
  RBasErr(eidx, 1);				// Normal non-fatal error
  if(eidx) RBasErr(0);				// Program aborting error
  std::cout << "\n";                                 // Keep screen nice
  exit(-1);                                     // Quit program
  }

void RBasic::RBasFatal(int eidx, const std::string& pname) const
  {
  RBasErr(eidx, pname, 1);			// Normal non-fatal error
  if(eidx) RBasErr(0);				// Program aborting error
  std::cout << "\n";                                 // Keep screen nice
  exit(-1);                                     // Quit program
  }
 
// ____________________________________________________________________________
// ii                        Relax Basic Setup Functions
// ____________________________________________________________________________

        // Input                RB      : Basic relaxation(this)
        //                      pset    : A parameter set
        //                      warn    : Warning level
        //                                      0 = no warning
        //                                      1 = warning, non-fatal
        //                                     >1 = fatal
        // Output               none    : RB number of parameters is set
        //                                from parameter in pset
        // Note                         : If # parameters read, returns TRUE
        //                                whether or not parameters exist

bool RBasic::GetNPoints(const ParameterSet& pset, int& ns, bool warn)
  {
  std::string pstate;
  ParameterSet::const_iterator item;         // Pix in parameter list
  std::string pname = std::string("NSpins");		// # spins par name
  item = pset.seek(pname);                      // Try & find NSpins par
  if(item != pset.end())                        // Retrieve the number spins
    (*item).parse(pname,ns,pstate);
 
// Can't Find Parameter NSpins. This Will Evoke Non-Fatal Warnings If Desired.
//    Then We'll Try & Just Count The Basic Relaxation Parameters Directly
 
  return false;
  }

// ____________________________________________________________________________
// iii                  Basic Relaxation Checking Functions
// ____________________________________________________________________________

/* These checking functions insure that a specficed relaxation rate,
   relaxation time, or linewidth is a reasonable value.                      */

void RBasic::CheckR1(double r1) { if(r1 < 0)  RBasFatal(40); }
void RBasic::CheckR2(double r2) { if(r2 < 0)  RBasFatal(40); }
void RBasic::CheckT1(double t1) { if(t1<= 0)  RBasFatal(41); }
void RBasic::CheckT2(double t2) { if(t2<= 0)  RBasFatal(41); }
void RBasic::CheckLW(double lw) { if(lw < 0)  RBasFatal(42); }

bool RBasic::CheckSpin(double i, bool warn) const
  {
  if(i<0 || i>=R1rates.size())			// If spin out of range
    {						//   Issue warning if desired
    if(warn) RBasErr(59, 1);			//   and return we failed
    return false;
    }
  return true;					// Spin count matches OK
  }

bool RBasic::CheckSpins(double N, bool warn) const
  {
  if(R1rates.size() != N)			// If spin count mismatch
    {						//   Issue warning if desired
    if(warn) RBasErr(60, 1);			//   and return we failed
    return false;
    }
  return true;					// Spin count matches OK
  }

bool RBasic::CheckHS(double hs, bool warn) const
  {
  if(HS() != hs)				// If Hilbert space mismatch
    {						//   Issue warning if desired
    if(warn) RBasErr(70, 1);			//   and return we failed
    return false;
    }
  return true;					// Hilbet space matches OK
  }

//bool RBasic::CheckAcq(bool warn)
//  {
// sosik
//  }

// ____________________________________________________________________________
// iv                      Longitudinal Setup Functions
// ____________________________________________________________________________

/*    Function                            Result
      ========  ===============================================================
       SetIzs   Makes/Stores {Izi} For A Specified System. This MUST Be Done If
                A New System Is Specified. Note That, As These Will Be In The
                Product Basis They Are Diagonal Operators & Use Little Storage.
      SetCinfs  Generates the Izi Coefficients At Infinite Time. The Spin 
                System MUST Have Been Previously Specified And All The {Izi}
                Operators Generated.
      SetCsigs  Generates the Izi Coefficients For Input Sigma. The Spin 
                System MUST Have Been Previously Specified And All The {Izi}
                Operators Generated.
      TestLong  Determines If Longitudinal Relaxation Will Be Detected. The
                Spin System MUST Have Been Previously Specified And All The
                {Izi} Operators Generated.                                   */

bool RBasic::SetIzs(const spin_sys& sys, bool warn)
  {
  int NS = sys.spins();				// Number of system spins
  if(!CheckSpins(NS, warn))			// Insure spin count matches
    { 						//   If mismatch, issue warnings
    if(warn) { RBasErr(61,1); RBasErr(63,1); }	//   if desired then return
    return false;				//   failed
    }
  Izis.clear();					// Remove any exitsting Izs
  for(int i=0; i<NS; i++)			// Loop over all spins
    Izis.push_back(Iz(sys,i));			// Store each Izi
  return true;
  }

bool RBasic::SetCinfs(const gen_op& sigma, bool warn)
  {
  if(!CheckHS(sigma.dim(), warn))		// Insure Hilbert space matches
    { 						//   If mismatch, issue warnings
    if(warn) { RBasErr(71,1); RBasErr(72,1); }	//   if desired then return
    return false;				//   failed
    }
  int NS = spins();				// Number of system spins
  Cinf = std::vector<double> (NS);			// Vector for coeffecients
  for(int i=0; i<NS; i++)			// Loop over all spins
    Cinf[i] = Re(proj(sigma,Izis[i]));		// Component Izi in sigma 
  return true;
  }

bool RBasic::SetCsigs(const gen_op& sigma, bool warn)
  {
  if(!CheckHS(sigma.dim(), warn))		// Insure Hilbert space matches
    { 						//   If mismatch, issue warnings
    if(warn) { RBasErr(71,1); RBasErr(62,1); }	//   if desired then return
    return false;				//   failed
    }
  int NS = spins();				// Number of system spins
  Csig = std::vector<double>(NS);			// Vector for coeffecients
  for(int i=0; i<NS; i++)			// Loop over all spins
    Csig[i] = Re(proj(sigma,Izis[i]));		// Component Izi in sigma 
  return true;
  }

bool RBasic::TestLong(double cutoff)
  {
  int LongRlx = 0;				// Assume no longitudinal rlx.
  Det.set_DBR();				// Insure Det in default basis
  int i,a;					// Spin, basis indices
  complex z;
  int hs = Det.dim();				// Get spin Hilbert space
  int ns = spins();				// Get number of spins
  for(i=0; i<ns; i++)				// Loop over all spins and
    {						// see if any detected relax.
    (Izis[i]).set_DBR();			//   Insure Izi comp. in DBR
    for(a=0; a<hs && !LongRlx; a++)		//   Loop diagonal elements
      {						//   that go into the trace 
      z = Det.get(a,a)*(Izis[i]).get(a,a);	//   and flag if anything is
      if(norm(z) > cutoff) return true;		//   detected at all
      }
    }
  return false;
  }


// ____________________________________________________________________________
// v                      Transverse Setup Functions
// ____________________________________________________________________________

/* This function sets up an array containing R2 times for each coherence in
   the system Hilbert space. If the associated spin system changes (i.e. its
   Hilbert space) or any relaxation rates are altered then this must be
   regnerated!                                                               */

void RBasic::SetR2Mx()
  {
  int hs = HS();				// The spin Hilbert space
  int NS = spins();				// Spins in the system
  R2mx = matrix(hs,hs,complex0,h_matrix_type);	// Our system R2 array
  double r2ave, delmtot, delm, msumsp;
  int a,b,i;					// Basis function indices
  for(a=0; a<hs-1; a++)				// Loop over U.T. density
    for(b=a+1; b<hs; b++) 			// operator elements 
      {
      r2ave = 0;				// Begin with no relaxation
      delmtot = Fzs.getRe(a) - Fzs.getRe(b);	// Total del Fz between a&b
      msumsp = 0;				// Begin
      for(i=0; i<NS; i++)			// Now loop over the spins
        { 					// indices & see which flip
        delm=FZxNS.getRe(a,i)-FZxNS.getRe(b,i); //   delmz of spin i <a| & |b>
        if(delm)				//   If spin flip, add in
          { 					//   relaxation contrib.
          r2ave += R2rates[i]/fabs(delm);	//   Add scaled to R2
          msumsp++;				//   Total spins flipped
          }
        }
      R2mx.put_h(r2ave,a,b);			// Store averaged R2 value
      }
  }

void RBasic::ZeroR2Mx()
  { R2mx = matrix(); }

// ----------------------------------------------------------------------------
// ------------------------------ PUBLIC FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// A                  RBasic File Constructors, Destructor
// ____________________________________________________________________________

/* These are the constructors of the class handling Bruker Relax Basic 1D
   acquisition binary files.  The default name such files is "fid". The
   output byte order is automatically set for the computer architecture.

        //                      name    : External filename 
	//			vx	: Data vector
	//			TD	: Total points (real + imag)
	//			byteord : Input byte order                   */

// ----------------------------------------------------------------------------
//                            Simple Constructors
// ----------------------------------------------------------------------------
 
RBasic::RBasic() { }

RBasic::RBasic(const row_vector& vx)
  { 
  for(int i=0; i<vx.size(); i++)
    {
    R1rates.push_back(vx.getRe(i));
    R2rates.push_back(vx.getIm(i));
    }
  }

RBasic::RBasic(const RBasic& RB)
  {
  R1rates = RB.R1rates;				// Copy R1 relax. rates (T1)
  R2rates = RB.R2rates;				// Copy R2 relax. rates (T2)
  FZxNS   = RB.FZxNS;				// Copy spin Fz values  (T2)
  Fzs     = RB.Fzs;				// Copy b.f. total Fzs  (T2)
  Izis    = RB.Izis;				// Copy the Iz ops.     (T1)
  SigInf  = RB.SigInf;				// Copy inf. time dnsop (T1)
  Cinf    = RB.Cinf;				// Copy Iz coeffs. tinf (T1)
  R2mx    = RB.R2mx;				// Copy R2 relax. mx    (T2)
  H0      = RB.H0;				// Copy static Ham.     (Ev/Aq)
  Det     = RB.Det;				// Copy detect op.	(Aq)
  Csig    = RB.Csig;				// Copy Iz coeffs. sigma
  }

// ----------------------------------------------------------------------------
//                          Destruction & Assignment
// ----------------------------------------------------------------------------

RBasic::~RBasic() { }

RBasic & RBasic::operator= (const RBasic& RB) 
  {
  R1rates = RB.R1rates;				// Copy R1 relax. rates (T1)
  R2rates = RB.R2rates;				// Copy R2 relax. rates (T2)
  FZxNS   = RB.FZxNS;				// Copy spin Fz values  (T2)
  Fzs     = RB.Fzs;				// Copy b.f. total Fzs  (T2)
  Izis    = RB.Izis;				// Copy the Iz ops.     (T1)
  SigInf  = RB.SigInf;				// Copy inf. time dnsop (T1)
  Cinf    = RB.Cinf;				// Copy Iz coeffs. tinf (T1)
  R2mx    = RB.R2mx;				// Copy R2 relax. mx    (T2)
  H0      = RB.H0;				// Copy static Ham.     (Ev/Aq)
  Det     = RB.Det;				// Copy detect op.	(Aq)
  return *this;
  }

// ____________________________________________________________________________
// B                    Basic Relaxation Access Functions
// ____________________________________________________________________________
 
// ----------------------------------------------------------------------------
//                     Access To Rates, Times, Linewidths
// ----------------------------------------------------------------------------

/* Note that the array R2mx contains elements related to the array R2rates
   based on the associated spin system. If either the array R2rates or the
   spin system are changed then R2mx must be regenerated!                    */

int    RBasic::HS() 	     const { return Izis.size()?Izis[0].dim():0; }
int    RBasic::spins() 	     const { return R1rates.size(); }	// No. spins
double RBasic::T1(int i)     const { return RB(i,0);  } 	// Single T1
double RBasic::T2(int i)     const { return RB(i,1);  }		// Single T2
double RBasic::R1(int i)     const { return RB(i,2);  }		// Single R1
double RBasic::R2(int i)     const { return RB(i,3);  }		// Single R2
double RBasic::LW(int i)     const { return RB(i,4);  }		// Single LW
double RBasic::RB(int i, int type) const
  {
  if(!CheckSpin(i))					// If spin out of range
    {
    std::string msg;
    switch(type)
      {
      case 0:  msg = "Longitudinal Relaxation Time"; break;
      case 1:  msg = "Transverse Relaxation Time";   break;
      case 2:  msg = "Longitudinal Relaxation Rate"; break;
      case 3:  msg = "Transverse Relaxation Rate";   break;
      case 4:  msg = "Half-Height Linewith";         break;
      }
    RBasFatal(50);	
    }
  switch(type)
    {
    case 0:  return 1./R1rates[i]; break;
    case 1:  return 1./R2rates[i]; break;
    case 2:  return R1rates[i];    break;
    case 3:  return R2rates[i];    break;
    case 4:  return R2rates[i]/PI; break;
    }
  return 0.0;
  }

void RBasic::T1(double val, int i) { RB(val,i,0); }	// Set T1
void RBasic::T2(double val, int i) { RB(val,i,1); }	// Set T2
void RBasic::R1(double val, int i) { RB(val,i,2); }	// Set R1
void RBasic::R2(double val, int i) { RB(val,i,3); }	// Set R2
void RBasic::LW(double val, int i) { RB(val,i,4); }	// Set LW
void RBasic::RB(double val, int i, int type)
  {
  if(!CheckSpin(i))					// If spin out of range
    {
    std::string msg;
    switch(type)
      {
      case 0:  msg = "Longitudinal Relaxation Time"; break;
      case 1:  msg = "Transverse Relaxation Time";   break;
      case 2:  msg = "Longitudinal Relaxation Rate"; break;
      case 3:  msg = "Transverse Relaxation Rate";   break;
      case 4:  msg = "Half-Height Linewith";         break;
      }
    RBasFatal(51);	
    }
  switch(type)
    {
    case 0:  R1rates[i] = 1.0/val;             break;
    case 1:  R2rates[i] = 1.0/val; ZeroR2Mx(); break;
    case 2:  R1rates[i] = val;                 break;
    case 3:  R2rates[i] = val;     ZeroR2Mx(); break;
    case 4:  R2rates[i] = PI*val;  ZeroR2Mx(); break;
    }
  }

std::vector<double> RBasic::T1s() const { return RBRates(0); }	// Set of T1s
std::vector<double> RBasic::T2s() const { return RBRates(1); }	// Set of T2s
std::vector<double> RBasic::R1s() const { return R1rates;    }	// Set of R1s
std::vector<double> RBasic::R2s() const { return R2rates;    }	// Set of R2s
std::vector<double> RBasic::LWs() const { return RBRates(2); }	// Set of LWs
std::vector<double> RBasic::RBRates(int type) const
  {
  int nspins = R1rates.size();
  std::vector<double> vals;
  int i=0;
  switch(type)
    {
    default:
    case 0:
      for(; i<nspins; i++)
        vals.push_back(T1(i));
      break;
    case 1:
      for(; i<nspins; i++)
        vals.push_back(T2(i));
      break;
    case 2:
      for(; i<nspins; i++)
        vals.push_back(PI*R2(i));
    }
  return vals;
  }
 
// ----------------------------------------------------------------------------
//                     Access To System & System Operators
// ----------------------------------------------------------------------------
  
/* Note that the spin system is NOT stored in class RBasic. Class RBasic is
   associated with a systems spin Hilbert space via the set of spin operators
   {Izi} as well as the arrays that track individual spin mz and total Fz 
   values for the basis functions, FZxNS and Fzs respectively. When the system
   is "set", the arrays FZxNS, Fzs, and {Izi} must be generated generated and
   at that point RBasic becomes associated with a spin Hilbert space. Prior to
   the system being set this is not true and without a Hilbert space we do not
   allow the evolution Hamiltonian nor the detection operator to be specfied.
   Indeed, if the spin system is set it will immediately remove any existing
   Hamiltonian (H0), density operator (D), and infinite time density operator
   (in Cinf).                                                                */

bool RBasic::SetSystem(const spin_sys& sys, int warn)
  {
//                        Set Up For T1 Relaxation
//                  System ==> {Izi}, SigmaInf == Czis       

  if(!SetIzs(sys, warn?true:false))		// Try to make all {Izi}  (T1)
    {						// If this fails then the (T1)
    if(warn)					// system has a different (T1)
      {						// spins count in it!     (T1)
      if(warn >= 2) RBasFatal(90);		  
      else          RBasErr(90);		  
      }						  
    return false;				  
    }						  
  if(!SetSigInf(sigma_eq(sys), warn?true:false))// Default SigInf [eq.]   (T1)
    { RBasFatal(64); return false; }

//                        Set Up For T2 Relaxation
//                      System ==> FZxNs, Fzs, R2Mx

  FZxNS = sys.qStates();			// Array of basis Fz's    (T2)
  Fzs   = sys.qnStates();			// Vector of total Fz's   (T2)
  SetR2Mx();					// Transverse R2's array  (T2)

//             Remove All Other System Dependent Parameters
//                          FZxNs, Fzs, R2Mx

  H0     = gen_op();				// Remove existing H
  Det    = gen_op();				// Remove existing D
  Csig.clear();					// Remove existing Iz coeffs.
  return true;
  }

bool RBasic::SetH0(const gen_op& H, int warn)
  {
  if(!CheckHS(H.dim(), warn?true:false))
    {
    if(warn)
      {
      if(!HS())     RBasErr(82, 1);
      else          RBasErr(83, 1);
      if(warn >= 2) RBasFatal(80);
      else          RBasErr(80);
      }
    return false;
    }
  H0 = H;
  return true;
  }

bool RBasic::SetDet(const gen_op& D, int warn)
  {
  if(!CheckHS(D.dim(), warn?true:false))
    {
    if(warn)
      {
      if(!HS())     RBasErr(82, 1);
      else          RBasErr(83, 1);
      if(warn >= 2) RBasFatal(81);
      else          RBasErr(81);
      }
    return false;
    }
  Det = D; 
  return true;
  }

// ----------------------------------------------------------------------------
//                    The Infinite Time Density Operator
// ----------------------------------------------------------------------------

/* The infinite time density operator represents the state of the system after
   an infinite evolution time. In simple cases this will be equilibrium, in
   other cases this will be some steady state matrix. The infinite time 
   density operator is used to generate the coefficients Cinfs.		     */

bool RBasic::SetSigInf(const gen_op& Sinf, int warn)
  {
  if(!CheckHS(Sinf.dim(), warn?true:false))	// Insure Hilbert space exists
    {						// and matches our own HS
    if(warn)
      {
      if(!HS())     RBasErr(82, 1);
      else          RBasErr(83, 1);
      if(warn >= 2) RBasFatal(85);
      else          RBasErr(85);
      }
    return false;
    }
  if(!SetCinfs(Sinf, warn?true:false))		// Set Izi inf coefficients
    {
    if(warn >= 2) RBasFatal(85);
    else          RBasErr(85);
    return false;
    }
  return true;
  }

// ----------------------------------------------------------------------------
//                    Transverse Relaxation System Values
// ----------------------------------------------------------------------------

/* The transverse relaxation matrix & superoperator depend upon both the spin
   relaxation rates (R2rates) and the spin system on which they are applied.
   Whenever any transverse rate changes, or if the spin system changes, these
   values must be recalculated. As such, the array R2mx MUST be zeroed if 
   such changes occur. The functions below will regenerate the values if the
   array is empty. If the array is not empty it has been previously made for
   the current R2rates and spin system. The matrix R2Mx will be HSxHS 
   whereas the R2 superoperator will be LSxLS.                               */

matrix RBasic::R2Mx()
  {
  if(!R2mx.rows())				// Rebuild R2 matrix if needed
    if(spins() && HS()) SetR2Mx(); 		//   Rebuild if Hilbert space
  return R2mx;					//   and set # of spins
  }

matrix RBasic::R2LOp()
  {
  int hs = HS();				// Our spin Hilbert space
  if(!R2mx.rows())				// Rebuild R2 matrix if needed
    if(spins() && HS()) SetR2Mx(); 		//   Rebuild if Hilbert space
  int ls = hs*hs;				// Our spin Liouville space
  int i,j,IJ;					// Basis function indices 
  matrix R2L(ls,ls,d_matrix_type);		// Set empty R2 LOp matrix
  for(i=0,IJ=0; i<hs; i++)			// Loop over Hilbert space
    for(j=0; j<hs; j++,IJ++) 			// basis functions
      R2L.put(R2mx.get(i,j), IJ, IJ);
  return R2L;
  }

// ----------------------------------------------------------------------------
//                    Longitudinal Relaxation System Values
// ----------------------------------------------------------------------------

/* The longitudinal relaxation rates are strictly specified in the vector
   R1rates. These change only when that vector is altered. The R1 superoperator
   simply contains the specified R1 values ordered by spin. It will be of
   dimension NSxNS and diagonal.                                             */

matrix RBasic::R1LOp()
  {
  int ns = spins();
  matrix R1L(ns,ns,d_matrix_type);
  for(int i=0; i<ns; i++)
    R1L.put(R1rates[i], i, i);
  return R1L;
  }

// ____________________________________________________________________________
// C                  Basic Relaxation Parameter Set Functions
// ____________________________________________________________________________

/*                 Single Spin Relaxation Parameter Functions

     Function                          Result
     ---------        ------------------------------------------
      ReadT2          T2 Value From pset Parameter [idx]T2(spin)
      ReadT1          T1 Value From pset Parameter [idx]T1(spin)
      ReadLW          LW Value From pset Parameter [idx]LW(spin)
      ReadR2          R2 Value From pset Parameter [idx]R2(spin)
      ReadR1          R1 Value From pset Parameter [idx]R1(spin)             */

double RBasic::ReadT2(const ParameterSet& pset,  int spin, int idx, int pf)
   { return ReadPar(pset, spin, 0, idx, pf); }
double RBasic::ReadT1(const ParameterSet& pset,  int spin, int idx, int pf)
   { return ReadPar(pset, spin, 1, idx, pf); }
double RBasic::ReadLW(const ParameterSet& pset,  int spin, int idx, int pf)
   { return ReadPar(pset, spin, 2, idx, pf); }
double RBasic::ReadR2(const ParameterSet& pset,  int spin, int idx, int pf)
   { return ReadPar(pset, spin, 3, idx, pf); }
double RBasic::ReadR1(const ParameterSet& pset,  int spin, int idx, int pf)
   { return ReadPar(pset, spin, 4, idx, pf); }
double RBasic::ReadPar(const ParameterSet& pset,int spin,
                                                 int type, int idx, int pf)
  {
  ParameterSet::const_iterator item;         // A pix into parameter list
  std::string pname, ssfile, pstate;                 // Items in each pset entry
  std::string Pre; 					// Name adjustment if indexed
  if(idx != -1) Pre = std::string("[") + Gdec(idx)	// to value but -1
                    + std::string("]");		// (this is a prefix)
  std::string Suff = std::string("(") + Gdec(spin)	// Spin adjustment for name
              + std::string(")");			// (this is a suffix)
  std::string Base;
  double val;
  switch(type)
    {
    case 0: 					// Shift anisotropy T2
      Base = std::string("T2"); val = DefT12; break;
    case 1: 					// Shift anisotropy T1
      Base = std::string("T1"); val = DefT12; break;
    case 2: 					// Shift anisotropy LWhh
      Base = std::string("LW"); val = DefLW; break;
    case 3: 					// Shift anisotropy R2
      Base = std::string("R2"); val = DefR12; break;
    case 4: 					// Shift anisotropy R1
      Base = std::string("R1"); val = DefR12; break;
    }
  pname = Pre + Base + Suff;			// Adjust name as desired
  item = pset.seek(pname);                      // Pix in parameter list
  if(item != pset.end())                        // Get node parameter
    (*item).parse(pname,val,pstate);
  else
    if(pf) std::cout << "\n\tCant Find " << pname
                << " in Parameter Set."
                << "\n\tSetting Value To "
                << val;
  return val;
  }

/*               Multiple Spin Relaxation Parameter Functions
 
            Function                      Result
            ---------   -------------------------------------------
            ReadT2s     T2 Values From pset Parameter [idx]T2(spin)
            ReadT1s     T1 Values From pset Parameter [idx]T1(spin)
            ReadLWs     LW Values From pset Parameter [idx]LW(spin)
            ReadR2s     R2 Values From pset Parameter [idx]R2(spin)
            ReadR1s     R1 Values From pset Parameter [idx]R1(spin)          */

std::vector<double> RBasic::ReadT2s(const ParameterSet& pset,int N,int idx,int pf)
   { return ReadPars(pset, N, 0, idx, pf); }
std::vector<double> RBasic::ReadT1s(const ParameterSet& pset,int N,int idx,int pf)
   { return ReadPars(pset, N, 1, idx, pf); }
std::vector<double> RBasic::ReadLWs(const ParameterSet& pset,int N,int idx,int pf)
   { return ReadPars(pset, N, 2, idx, pf); }
std::vector<double> RBasic::ReadR2s(const ParameterSet& pset,int N,int idx,int pf)
   { return ReadPars(pset, N, 3, idx, pf); }
std::vector<double> RBasic::ReadR1s(const ParameterSet& pset,int N,int idx,int pf)
   { return ReadPars(pset, N, 4, idx, pf); }
std::vector<double> RBasic::ReadPars(const ParameterSet& P,  int N,
                                                      int type,int idx,int pf)
  {
  std::vector<double> vals;
  for(int i=0; i<N; i++)
    vals.push_back(ReadPar(P,i,type,idx,pf));
  return vals;
  }

// ____________________________________________________________________________
// D                    Basic Relaxation Input Functions
// ____________________________________________________________________________


int RBasic::read(const std::string& filein, int idx, int warn)
  {
  ParameterSet pset;
  pset.read(filein);
  return read(pset, idx, warn?0:1);
  }

int RBasic::read(const ParameterSet& pset,int idx, int warn)
  {
  int ns = 0;
  int TF = GetNPoints(pset, ns, warn?0:1);	// Get number of spins

//            Try & Read In All Possible Parameters First

  std::vector<double> R1s, R2s;			// Arrays we want
  DefR12 = -1.0;				// Set default to odd value
  R1s = ReadR1s(pset,  ns, idx, 0);		// Try & read R1 values
  R2s = ReadR2s(pset,  ns, idx, 0);		// Try & read R2 values
  DefT12 = -1.0;				// Set default to odd value
  std::vector<double> T1s, T2s;			// Arrays we want
  T1s = ReadT1s(pset,  ns, idx, 0);		// Try & read T1 values
  T2s = ReadT2s(pset,  ns, idx, 0);		// Try & read T2 values
  DefLW = -1.0;					// Set default to odd value
  std::vector<double> LWs;				// Array we want
  LWs = ReadLWs(pset,  ns, idx, 0);		// Try & read LW values

//    Now Set The Array "rates" Using These Parameters (Hierarchically)

  double R1val, R2val;
  for(int i=0; i<ns; i++)			// Loop over # parameters
    {
    if(R1s[i] != DefR12)	{ CheckR2(R1s[i]); R1val = R1s[i];     }
    else if(T1s[i] != DefT12)   { CheckR2(T1s[i]); R1val = 1.0/T1s[i]; }
    else                                           R1val = 0.0;  
    if(R2s[i]!= DefR12)         { CheckR2(R2s[i]); R2val = R2s[i];     }
    else if(T2s[i] != DefT12)   { CheckT2(T2s[i]); R2val = 1.0/T2s[i]; }
    else if(LWs[i] != DefLW)    { CheckLW(LWs[i]); R2val = PI*LWs[i];  }
    else                                           R2val = 0.0;
    R1rates.push_back(R1val);
    R2rates.push_back(R2val);
    }
  return TF;
  }  

// ____________________________________________________________________________
//      Functions For Relaxation/Evolution Using Longitudinal Coefficients
// ____________________________________________________________________________

gen_op RBasic::SigmaT1(const gen_op& sigma)
  {
  gen_op sigT1;					// Empty operator
  int ns = spins();				// Number of spins
  SetCsigs(sigma);				// Genereate coefficients
  for(int i=0; i<ns; i++)			// Loop over spins
    sigT1 += Csig[i]*Izis[i];			// and sum sigT1i's
  return sigT1;
  }

gen_op RBasic::SigmaT2(const gen_op& sigma)
  { return sigma - SigmaT1(sigma); }

col_vector RBasic::SigmaC(const gen_op& sigma)
  {
  int ns = spins();				// Number of spins
  int hs = HS();				// Spin Hilbert space 
  int dim = spins() + hs*hs;			// Liouville dimension
  col_vector sigmaC(dim);			// Make blank vector
  SetCsigs(sigma);				// Genereate coefficients
  int i,j,k;					// Spin, basis indices
  for(i=0; i<ns; i++)				// Loop over the spins
    sigmaC.put(Csig[i], i);			// & set coefficients
  gen_op sigT2 = SigmaT2(sigma);		// Get sigma T2 operator
  sigT2.set_DBR();				// Do this in default basis
  for(i=0, k=ns; i<hs; i++)			// Loop over basis functions
    for(j=0; j<hs; j++, k++)			// Loop over basis functions
      sigmaC.put(sigT2(i,j), k);		// & set coefficients
  return sigmaC;
  }

col_vector RBasic::SigmaCEq(const gen_op& sigmaeq)
  {
  SetSigInf(sigmaeq);
  return SigmaCEq();
  }

col_vector RBasic::SigmaCEq()
  {
  int ns = spins();				// Number of spins
  int hs = HS();				// Spin Hilbert space 
  int dim = spins() + hs*hs;			// Liouville dimension
  col_vector sigmaCeq(dim, 0);			// Make blank vector
  for(int i=0; i<ns; i++)			// Loop over the spins
    sigmaCeq.put(Cinf[i], i);			// & set coefficients
  return sigmaCeq;
  }

matrix RBasic::RC()
  {
  int ns = spins();				// Number of spins
  int hs = HS();				// Spin Hilbert space 
  int dim = spins() + hs*hs;			// Liouville dimension
  matrix Rc(dim,dim,d_matrix_type);		// Zero relaxation matrix
  int i,j,k;					// Spin, basis indices
  for(i=0; i<ns; i++)				// Loop over the spins
    Rc.put(R1rates[i], i, i);			// & set longitudinal values
  matrix R2sm = R2Mx();				// Insure R2Mx exists
  for(i=0, k=ns; i<hs; i++)			// Loop over basis functions
    for(j=0; j<hs; j++, k++)			// Loop over basis functions
      Rc.put(R2sm.get(i,j), k, k);		// & set longitudinal values
  return Rc;
  }

matrix RBasic::HC(const gen_op& H)
  {
  SetH0(H);					// Set this as the Hamiltonian
  return HC();					// User overloaded function
  }

matrix RBasic::HC()
  {
  int ns = spins();				// Number of spins
  int hs = HS();				// Spin Hilbert space 
  int dim = spins() + hs*hs;			// Liouville dimension
  matrix Hc(dim,dim,d_matrix_type);		// Zero relaxation matrix
  H0.set_DBR();					// Put Hamiltonian in DBR
  matrix H2 = Hsuper(H0).get_mx();		// Convert to Liouville space
  Hc.put_block(ns, ns, H2);			// Set transverse block only
  return Hc;
  }

gen_op RBasic::Sigma(const col_vector& sigmaC)
  {
  int ns = spins();				// Number of spins
  int hs = HS();				// Spin Hilbert space 
  gen_op sigma;					// Make blank array
  int i,j,k;					// Spin, basis indices
  for(i=0, k=ns; i<hs; i++)			// Loop over basis functions
    for(j=0; j<hs; j++, k++)			// Loop over basis functions
      sigma.put(sigmaC.get(k), i, j);		// & set transverse elements
  for(i=0; i<ns; i++)				// Add in the longitudinal
    sigma += sigmaC.get(i)*Izis[i];		// part via coefficients
  return gen_op(sigma);
  }

// sosik

// ____________________________________________________________________________
// E               Basic Relaxation Time Evolution Functions
// ____________________________________________________________________________

/* These functions will perform an ad-hoc evolution of the density operator
   under a constant Hamiltonian and basic relaxation. The initial density 
   operator sig0 will evolve for the specified time under the effects of this
   simple relaxation scheme as well as under the effects of any input 
   Hamiltonian. It is tacitly assumed that evolution under the Hamiltonian is
   not coupled to evolution under relaxation. Hence, the evolution will first
   occur under the static Hamiltonian HO if one has been set, then under the
   phenomenological relaxation scheme of this class.                         */


gen_op RBasic::Evolve(const spin_sys& sys, const gen_op& sigmap, double t)
  {
//            Begin With Evolution Under The Active Hamiltonian

  gen_op sigmat = sigmap;			// The sigma to be evolved
  gen_op siginf = sigma_eq(sys);		// The equilibrium operator
  if(H0.dim()) evolve_ip(sigmat, H0, t);	// Evolve under H0 if set
  sigmat.set_DBR();				// Put into default basis

//                    First We Perform Ad Hoc T1 Evolution  

  if(!SetIzs(sys))      RBasFatal(64); 		// Set up all Izi operators
  if(!SetCinfs(siginf)) RBasFatal(64);		// Set siginf Izi coefficients
  return Evolve(sigmap, t);			// Use overload

/*
  if(!SetCsigs(sigmap)) RBasFatal(64);		// Set sigmap Izi coefficients
  int ns = spins();
  double KIzi, eR1t;
  for(int i=0; i<ns; i++)			// Loop over all spins
    {
    eR1t = exp(-R1rates[i]*t);			// The exponential factor
    KIzi = (Csig[i]-Cinf[i])*eR1t + Cinf[i];	// Relaxed Izi coefficient
    sigmat += (KIzi-Csig[i])*Izis[i];		// Set relaxed Izi component
    }

//                    Second We Perform Ad Hoc T2 Evolution  

  matrix R2val = R2Mx();			// Get array of R2 values
  int hs = sigmap.dim();			// Spin Hilbert space
  int a,b;					// Basis function indices
  complex z;					// Density operator element
  for(a=0; a<hs-1; a++)				// Loop over U.T. density
    for(b=a+1; b<hs; b++) 			// operator elements 
      {
      z = sigmat.get(a,b) * exp(-R2val(a,b)*t);	// Element to relax
      sigmat.put(z,a,b);			// Set relaxed <a|sigma(t)|b>
      sigmat.put(conj(z),b,a);			// Set relaxed <b|sigma(t)|a>
      }
  return sigmat;
*/
  }

gen_op RBasic::Evolve(const gen_op& sigmap, double t)
  {
//         First Evolve Under The Static Hamiltonian (If Existing)

  gen_op sigmat = sigmap;			// The sigma to be evolved
  if(H0.dim()) evolve_ip(sigmat, H0, t);	// Evolve under H0 if set

//              Next We Perform Phenomenological T1 Evolution  
//      Only Affects T1 Components, These Are Izi On The Diagonal

  sigmat.set_DBR();				// Put into default basis
  if(!SetCsigs(sigmap)) RBasFatal(64);		// Set sigmap Izi coefficients
  int ns = spins();				// Number of spins in system
  double KIzi, eR1t;
  for(int i=0; i<ns; i++)			// Loop over all spins
    {
    eR1t = exp(-R1rates[i]*t);			// The exponential factor
    KIzi = (Csig[i]-Cinf[i])*eR1t + Cinf[i];	// Relaxed Izi coefficient
    sigmat += (KIzi-Csig[i])*Izis[i];		// Set relaxed Izi component
    }

//              Last We Perform Phenomenological T2 Evolution  
//        This Should Leave The T1 Components Un-Touched (Diagonals)

  matrix R2val = R2Mx();			// Get array of R2 values
  int hs = sigmap.dim();			// Spin Hilbert space
  int a,b;					// Basis function indices
  complex z;					// Density operator element
  for(a=0; a<hs-1; a++)				// Loop over U.T. density
    {
    for(b=a+1; b<hs; b++) 			// operator elements 
      {
      z = sigmat.get(a,b) * exp(-R2val(a,b)*t);	// Element to relax
      sigmat.put(z,a,b);			// Set relaxed <a|sigma(t)|b>
      sigmat.put(conj(z),b,a);			// Set relaxed <b|sigma(t)|a>
      }
    }

  return sigmat;				// Return evolved density op.
  }

// ____________________________________________________________________________
// F                        Basic Acquisition Functions
// ____________________________________________________________________________

/* These functions will perform an ad-hoc evolution of the density operator
   under basic relaxation during an acquisition.  The density operator sig0
   will evolve for the time specified.  The evolution will first occur under
   the static Hamiltonian HO if one has be set, then under the pseudo-
   relaxation scheme of this class.

		Note			Assumes that the static Hamiltonain
				  	and the detection operator are set!
   									     */
  
void RBasic::FID(const gen_op& sigmap, double td, row_vector& fid, int N)
  {
  gen_op sigmat = sigmap;			// The evolved sigma
  bool LongRlx  = TestLong();			// Flag if longitudinal relax.
//bool RBasic::CheckAcq(bool warn)
//  if(!LongRlx)					// If T1 processes not needed
//    {						// we just work with T2 only
// sosik
//cout << "\n\tWorking With T2 Only";
//    }

//  gen_op siginf = sigma_eq(sys);		// The equilibrium operator
  if(LongRlx)					// If we include long. relax
    {						// then we'll need Izi coeffs.
//    if(!SetCinfs(siginf)) RBasFatal(65);	//   Set siginf Izi coeffs.
    if(!SetCsigs(sigmap)) RBasFatal(65);	//   Set sigmap Izi coeffs.
    }
  int hs = sigmap.dim();			// Spin Hilbert space
  int i;					// Spin index
  double KIzi, eR1t;				// For longitudinal relax.
  int a,b;					// Basis function indices
  complex z, eR2t;				// For transverse relax.

//              Begin Looping Over the Acquisition Points
  double time = 0;
  int k;
  int ns = spins();

  if(LongRlx)
  {
  for(k=0; k<N; k++)
    {
    fid(k) = trace(Det, sigmat);		// Get point k   
    if(k<N-1)					// Evolve to next point
      {
      time += td;				//   Time at point k
      if(H0.dim())				//   Evolve under H0 
        sigmat = evolve(sigmap, H0, time);	//   if Hamiltonian is set
      sigmat.set_DBR();				// Put into default basis
      for(i=0; i<ns; i++)			// Loop over all spins and
        {					// relax longitidunally
        eR1t = exp(-R1rates[i]*time);		//   The exponential factor
        KIzi = (Csig[i]-Cinf[i])*eR1t + Cinf[i];//   Relaxed Izi coefficient
        sigmat += (KIzi-Csig[i])*Izis[i];	//   Set relaxed Izi component
        }
      for(a=0; a<hs-1; a++)			// Loop over U.T. density
        for(b=a+1; b<hs; b++) 			// operator elements and relax
          {					// in transverse plane
          eR2t = exp(-R2mx(a,b)*time); 		//   The exponential factor
          z = eR2t * sigmat.get(a,b);		//   Relax the element
          sigmat.put(z,a,b);			//   Set relaxed <a|sigma(t)|b>
          sigmat.put(conj(z),b,a);		//   Set relaxed <b|sigma(t)|a>
          }
      }
    }
    }

//			   Only Using T2 Block For Relaxation

  if(!LongRlx)
    {
    super_op L = R2LOp();			// Liouvillian with R2 relax.
    if(H0.dim())				// If a Hamiltonian is present
      L += complexi*Hsuper(H0);			// add it to the Liouvillian
    acquire1D ACQ(Det,L);			// Set up for acquisition
    fid = ACQ.T(sigmap, N, td);			// Generate FID
    }

  }

void RBasic::FID(const spin_sys& sys, const gen_op& sigmap,
                                            double td, row_vector& fid, int N)
  {
//                      First Set Up Working Variables

  gen_op sigmat = sigmap;			// The evolved sigma
  gen_op siginf = sigma_eq(sys);		// The equilibrium operator
  if(!SetIzs(sys)) RBasFatal(65); 		// Set up all Izi operators
  int LongRlx = TestLong();			// Check if long. relax. needed
  if(LongRlx)					// If we include long. relax
    {						// then we'll need Izi coeffs.
    if(!SetCinfs(siginf)) RBasFatal(65);	//   Set siginf Izi coeffs.
    if(!SetCsigs(sigmap)) RBasFatal(65);	//   Set sigmap Izi coeffs.
    }
  matrix R2val = R2Mx();			// Get array of R2 values
  int hs = sigmap.dim();			// Spin Hilbert space
  int i;					// Spin index
  double KIzi, eR1t;				// For longitudinal relax.
  int a,b;					// Basis function indices
  complex z, eR2t;				// For transverse relax.
  int ns = spins();


//              Begin Looping Over the Acquisition Points

  double time = 0;
  for(int k=0; k<N; k++)
    {
    fid(k) = trace(Det, sigmat);		// Get point k   
    if(k<N-1)					// Evolve to next point
      {
      time += td;				//   Time at point k
      if(H0.dim())				//   Evolve under H0 
        sigmat = evolve(sigmap, H0, time);	//   if Hamiltonian is set
      sigmat.set_DBR();				// Put into default basis
      for(i=0; i<ns && LongRlx; i++)		// Loop over all spins and
        {					// relax longitidunally
        eR1t = exp(-R1rates[i]*time);		//   The exponential factor
        KIzi = (Csig[i]-Cinf[i])*eR1t + Cinf[i];//   Relaxed Izi coefficient
        sigmat += (KIzi-Csig[i])*Izis[i];	//   Set relaxed Izi component
        }
      for(a=0; a<hs-1; a++)			// Loop over U.T. density
        for(b=a+1; b<hs; b++) 			// operator elements and relax
          {					// in transverse plane
          eR2t = exp(-R2val(a,b)*time); 	//   The exponential factor
          z = eR2t * sigmat.get(a,b);		//   Relax the element
          sigmat.put(z,a,b);			//   Set relaxed <a|sigma(t)|b>
          sigmat.put(conj(z),b,a);		//   Set relaxed <b|sigma(t)|a>
          }
      }
    }
  }


row_vector RBasic::FID(const spin_sys& sys, const gen_op& sigmap,
                                                           double td, int N)
  {
  row_vector data(N);
  FID(sys, sigmap, td, data, N);
  return data;
  }
  
row_vector RBasic::FID(const gen_op& sigmap, double td, int N)
  {
  row_vector data(N);
  FID(sigmap, td, data, N);
  return data;
  }
  



// ____________________________________________________________________________
// G                    Basic Relaxation Output Functions
// ____________________________________________________________________________


std::ostream& RBasic::print(std::ostream& ostr, bool hdr) const
  {
  double R1, R2, T1, T2, LW;				// For printed values
  int ns = spins();					// Number of spins
  int i = 0;						// Just an index
  std::string Title("Basic Relaxation Parameters");		// Title to Output
  int ncols = 5;					// Number of columns
  std::string hdrS[5] = { "R1 Rate", "R2 Rate", 		// Column headers
                     "T1 Time", "T2 Time",
                                "Linewidth" };
  int hdrlens[5]; 					// Column header widths
  for(; i<ncols; i++) hdrlens[i] = hdrS[i].length();	// Store header widths

  int sl = 7;						// Spacer width
  std::string space(sl, ' ');				// String for spacer

  std::string Header, HeadUL;				// Column header line
  for(i=0; i<ncols; i++)				// Construct strings
    {							// for column header
    Header   += hdrS[i] + space;			// and header underline
    HeadUL += std::string(hdrlens[i], '-') + space;
    }

  int ml = 40-Header.length()/2;			// Margin width
  if(hdr)
  ostr << std::string(40-Title.length()/2, ' ') 		// Output the title
       << Title << "\n\n";				// centered
  if(!ns)
    {
    ostr << std::string(30, ' ')
         << "No Parameters Defined";
    return ostr;
    }
  ostr << std::string(ml, ' ') << Header << "\n";		// Output column hdrs
  ostr << std::string(ml, ' ') << HeadUL; 			// Output col. hdr uls

  for(i=0; i<ns; i++)				// Loop over the spins
    {
    R1 = R1rates[i];					//   R1 for this spin
    R2 = R2rates[i];					//   R2 for this spin
    T1 = 1/R1;						//   R1 for this spin
    T2 = 1/R2;						//   R2 for this spin
    LW = R2/PI;						//   LW for this spin
    

/*                   <--marg-->R1 Rate<--sl-->
                     <--marg-->-------<--sl-->
                            xxxxxx.xxx/uu                                    */

    ostr << "\n" << std::string(ml-3, ' ');
    if(fabs(R1) > 1.e9)      ostr << Gform("%10.3f/ns", R1*1.e-9);
    else if(fabs(R1) > 1.e6) ostr << Gform("%10.3f/us", R1*1.e-6);
    else if(fabs(R1) > 1.e3) ostr << Gform("%10.3f/ms", R1*1.e-3);
    else                     ostr << Gform("%10.3f/s ", R1);

/*                   <--sl-->R2 Rate<-sl>
                     <--sl-->-------<-sl>
                     /uu  xxxxxx.xxx/uu                                    */

    ostr << std::string(sl-6, ' ');
    if(fabs(R2) > 1.e9)      ostr << Gform("%10.3f/ns", R2*1.e-9);
    else if(fabs(R2) > 1.e6) ostr << Gform("%10.3f/us", R2*1.e-6);
    else if(fabs(R2) > 1.e3) ostr << Gform("%10.3f/ms", R2*1.e-3);
    else                     ostr << Gform("%10.3f/s ", R2);

/*                   <--sl-->T1 Time<-sl>
                     <--sl-->-------<-sl>
                     /uu  xxxxxx.xxx uu                                    */

    ostr << std::string(sl-6, ' ');
    if(fabs(T1) < 1.e-9)      ostr << Gform("%10.3f ns", T1*1.e9);
    else if(fabs(T1) < 1.e-6) ostr << Gform("%10.3f us", T1*1.e6);
    else if(fabs(T1) < 1.e-3) ostr << Gform("%10.3f ms", T1*1.e3);
    else                      ostr << Gform("%10.3f s ", T1);

/*                   <--sl-->T2 Time<-sl>
                     <--sl-->-------<-sl>
                      uu  xxxxxx.xxx uu                                    */

    ostr << std::string(sl-6, ' ');
    if(fabs(T2) < 1.e-9)      ostr << Gform("%10.3f ns", T2*1.e9);
    else if(fabs(T2) < 1.e-6) ostr << Gform("%10.3f us", T2*1.e6);
    else if(fabs(T2) < 1.e-3) ostr << Gform("%10.3f ms", T2*1.e3);
    else                      ostr << Gform("%10.3f s ", T2);

/*                   <--sl-->Linewidth
                     <--sl-->---------
                      uu    xxxxxx.xxx uuu                                   */


//                     Output the Linewidth (Width 14)

    ostr << std::string(sl-4, ' ');
    if(fabs(LW) > 1.e11)     ostr << Gform("%10.3f GHz", LW*1.e-9);
    else if(fabs(LW) > 1.e8) ostr << Gform("%10.3f MHz", LW*1.e-6);
    else if(fabs(LW) > 1.e5) ostr << Gform("%10.3f KHz", LW*1.e-3);
    else                     ostr << Gform("%10.3f Hz ", LW);
    }
  ostr << "\n";
  return ostr;
  }

std::ostream& operator << (std::ostream& ostr, const RBasic& RB)
  { return RB.print(ostr); }

// ____________________________________________________________________________
// H                   Basic Relaxation Auxiliary Functions
// ____________________________________________________________________________

// sosik
std::vector<double> RBasic::FzCoeffs(const spin_sys& sys, const gen_op& sigma)
  {
  int NS = sys.spins();				// Number of system spins
  if(!CheckSpins(NS))				// Insure spin count matches
    { RBasErr(61, 1); RBasFatal(62); }		// If not, quit after errors
  std::vector<double> Coeffs(NS);			// Vector for coeffecients
  gen_op IZi;					// For Izi operators
  for(int i=0; i<NS; i++)			// Loop over all spins
    {
    IZi = Iz(sys,i);				// Get Izi
    Coeffs[i] = Re(proj(sigma,IZi));		// Component Izi in sig0 
    }
  return Coeffs;
  }

#endif 						// RelaxBas.cc
