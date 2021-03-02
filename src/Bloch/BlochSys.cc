/* BlochSys.cc **************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Bloch Spin System                           Implementation 	**
**                                                                      **
**      Copyright (c) 2001                                              **
**      Dr. Scott A. Smith						**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** The class BlochSys defines a spin system used as a basis for         **
** simulations based on the phenomenological Bloch equations. The       **
** Bloch equations account for spin magnetization evolution under       **
** a static Bo field, an applied rf-field (B1), simplistic relaxation   **
** (T1 & T2), and possibly exchange between magnetizaiton vectors.      **
**                                                                      **
** Bloch spin system tracks any number of spins and their associated    **
** magnetization vectors. To each spin, and all of its associated       **
** magnetizationv vectors, it assigns single longitudinal & transverse  **
** relaxation rates. Between any two spins (and their magnetization     **
** vectors) it assigns an exchange rate.                                **
**                                                                      **
** Based on this Bloch spin system, functions can easily built that     **
** generate magnetization vectors, relaxation matrices, exchange        **
** matrices, and evolution matrices with/without an rf-field present.   **
**                                                                      **
*************************************************************************/

#ifndef   BlochSys_cc_			// Is file already included?
#  define BlochSys_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Bloch/BlochSys.h>		// Includes the interface 
#include <Basics/Gutils.h>              // Include GAMMA errors/queries
#include <Basics/StringCut.h>		// Include Gdec function
#include <HSLib/HSLibIF.h>		// Include Hilbert space module
#include <Level2/acquire1D.h>		// Inlcude 1D rapid acquisitions
//#include <Matrix/matrix.h>		// Include GAMMA matrices
//#include <string>			// Include libstdc++ strings

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                     Bloch System Error Handling
// ____________________________________________________________________________


        // Input                sys     : Spin system (this)
        //                      eidx    : Error index
        //                      noret   : Flag for linefeed (0=linefeed)
	// Output		void	: An error message is output

/* The following error messages use the defaults set in the Gutils package

                Case                          Error Message

		(0) 			Program Aborting.....
    		(3) 			Can't Construct From Parameter Set
		(4)			Cannot Construct From Input File
    		(5)			Cannot Write To Parameter File
    		(6)			Cannot Write To Output FileStream
    		(7)			Cannot Produce Any Output
    		default			Unknown Error                        */

        // Input                sys     : Spin system (this)
        //                      eidx    : Error index
	//			pname	: Additional error message
        //                      noret   : Flag for linefeed (0=linefeed)
        // Output               none    : Error message output

/* The following error messages use the defaults set in the Gutils package

                Case                          Error Message

		(1) 			Problems With File PNAME
		(2)			Cannot Read Parameter PNAME
    		default			Unknown Error - PNAME                */
        // Input                sys     : Spin system (this)
        //                      eidx    : Error index
	//			pname	: Additional error message
        // Output               none    : Error message output
        //                                Program execution stopped


void BlochSys::BSerror(int eidx, int noret) const
  {
  std::string hdr("Bloch Spin System");
  std::string msg;
  switch (eidx)
    {
    case  9: msg = std::string("Spin Count Disagrees With System");	// (9)
             GAMMAerror(hdr,msg,noret);  break;
    case 10: msg = std::string("Dimensioning Mismatch"); 		// (10)
             GAMMAerror(hdr,msg,noret);  break;
    case 11: msg = std::string("Improper Number Of R1 Values"); 	// (11)
             GAMMAerror(hdr,msg,noret);  break;
    case 12: msg = std::string("Improper Number Of R2 Values"); 	// (12)
             GAMMAerror(hdr,msg,noret);  break;
    case 13: msg = std::string("Improper Number Of Isotope Values");	// (13)
             GAMMAerror(hdr,msg,noret);  break;
    case 14: msg = std::string("Improper Number Of Exchange Rates");	// (14)
             GAMMAerror(hdr,msg,noret);  break;
    case 15: msg = std::string("Improper Number Of Vector Norms");	// (15)
             GAMMAerror(hdr,msg,noret);  break;
    case 16: msg = std::string("Improper Number Of Coordinates");	// (16)
             GAMMAerror(hdr,msg,noret);  break;
    case 20: msg = std::string("Negative Relaxation Rate Specified");// (20)
             GAMMAerror(hdr,msg,noret);  break;
    case 21: msg = std::string("Negative Exchange Rate Specified");	// (21)
             GAMMAerror(hdr,msg,noret);  break;
    case 22: msg = std::string("Negative Vector Norm Specified");	// (22)
             GAMMAerror(hdr,msg,noret);  break;
    case 29: msg = std::string("Cannot Set Relaxation Rates");	// (29)
             GAMMAerror(hdr,msg,noret);  break;
    case 30: msg = std::string("Cannot Set R1 Relaxation Rates");	// (30)
             GAMMAerror(hdr,msg,noret);  break;
    case 31: msg = std::string("Cannot Set R2 Relaxation Rates");	// (31)
             GAMMAerror(hdr,msg,noret);  break;
    case 32: msg = std::string("Cannot Set Spin Types");		// (32)
             GAMMAerror(hdr,msg,noret);  break;
    case 33: msg = std::string("Cannot Set Exchange Rates");		// (33)
             GAMMAerror(hdr,msg,noret);  break;
    case 34: msg = std::string("Cannot Set Vector Norms");		// (34)
             GAMMAerror(hdr,msg,noret);  break;
    case 35: msg = std::string("Cant Set Up Magnetiztion Vector(s)");// (35)
             GAMMAerror(hdr,msg,noret);  break;
    case 40: msg = std::string("Parameters Insufficient");		// (40)
             GAMMAerror(hdr,msg,noret);  break;
    case 50: msg = std::string("Cannot Determine Vector Offset");	// (50)
             GAMMAerror(hdr,msg,noret);  break;
    case 51: msg = std::string("Cannot Determine Isotope Type");	// (51)
             GAMMAerror(hdr,msg,noret);  break;
    case 52: msg = std::string("Cannot Determine R1 Relax. Rate");	// (52)
             GAMMAerror(hdr,msg,noret);  break;
    case 53: msg = std::string("Cannot Determine R2 Relax. Rate");	// (53)
             GAMMAerror(hdr,msg,noret);  break;
    case 54: msg = std::string("Cannot Determine Vector Components");// (54)
             GAMMAerror(hdr,msg,noret);  break;
    case 55: msg = std::string("Cannot Determine Spin Assignment");	// (55)
             GAMMAerror(hdr,msg,noret);  break;
    case 56: msg = std::string("Can't Determine # Of Sub-Vectors");	// (56)
             GAMMAerror(hdr,msg,noret);  break;
    default: GAMMAerror(hdr, eidx, noret); break;		//(-1)
    }
  }

void BlochSys::BSerror(int eidx, const std::string& pname, int noret) const
  {
  std::string hdr("Bloch Spin System");
  std::string msg;
  switch(eidx)
    {
    case 101:                                                         // (101)
      msg = std::string("Can't Find Parameters For ")
          + pname;
     GAMMAerror(hdr, msg, noret); break;
    case 102:                                                         // (102)
      msg = std::string("Can't Find ") + pname 
          + std::string(" In Parameter Set");
     GAMMAerror(hdr, msg, noret); break;
    default: GAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error   (-1)
    }
  }

volatile void BlochSys::BSfatal(int eidx) const
  {  
  BSerror(eidx, 1);				// Normal non-fatal error
  if(eidx) BSerror(0);				// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

volatile void BlochSys::BSfatal(int eidx, const std::string& pname) const
  {  
  BSerror(eidx, pname, 1);			// Normal non-fatal error
  if(eidx) BSerror(0);				// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii                   Bloch System Checking Functions
// ____________________________________________________________________________

bool BlochSys::CheckR1s(const std::vector<double>& R1s, bool warn) const
  {
  unsigned nm = Offsets.size();			// Number of vectors
  if(nm != R1s.size())
    {
    if(warn)
      {
      BSerror(10, 1);				// Dimensioning problem
      BSerror(11, 1);				// Improper # of R1 values
      }
    return false;
    }
  for(unsigned i=0; i<nm; i++)
    if(R1s[i] < 0)
      {
      if(warn) BSerror(20, 1);			// Negative rate set
      return false;
      }
  return true;
  }

bool BlochSys::CheckR2s(const std::vector<double>& R2s, bool warn) const
  {
  unsigned nm = Offsets.size();			// Number of vectors
  if(nm != R2s.size())
    {
    if(warn)
      {
      BSerror(10, 1);				// Dimensioning problem
      BSerror(12, 1);				// Improper # of R2 values
      }
    return false;
    }
  for(unsigned i=0; i<nm; i++)
    if(R2s[i] < 0)
      {
      if(warn) BSerror(20, 1);			// Negative rate set
      return false;
      }
  return true;
  }

bool BlochSys::CheckIsos(const std::vector<Isotope>& Is, bool warn) const
  {
  unsigned nm = Offsets.size();			// Number of vectors
  if(nm != Is.size())
    {
    if(warn)
      {
      BSerror(10, 1);				// Dimensioning problem
      BSerror(13, 1);				// Improper # of R2 values
      }
    return false;
    }
  return true;
  }

bool BlochSys::CheckKs(const std::vector<double>& Ks, bool warn) const
  {
  unsigned nm = Offsets.size();			// Number of vectors
  if(nm != Ks.size())
    {
    if(warn)
      {
      BSerror(10, 1);				// Dimensioning problem
      BSerror(14, 1);				// Improper # of K values
      }
    return false;
    }
  for(unsigned i=0; i<nm; i++)
    if(Ks[i] < 0)
      {
      if(warn) BSerror(21, 1);			// Negative rate set
      return false;
      }
  return true;
  }

bool BlochSys::CheckSpins(int ns1, int ns2, bool warn) const
  {
  if(ns1 != ns2)
    {
    if(warn) { BSerror(10,1); BSerror(9,1); }
    return false;
    }
  return true;
  }

bool BlochSys::CheckNorms(const std::vector<double>& norms, bool warn) const
  {
  int nm = Offsets.size();			// Number of vectors
  if(nm != size())
    {
    if(warn)
      {
      BSerror(10, 1);				// Dimensioning problem
      BSerror(15, 1);				// Improper # of norms
      }
    return false;
    }
  for(int i=0; i<nm; i++)
    if(norms[i] < 0)
      {
      if(warn) BSerror(22, 1);			// Negative norm set
      return false;
      }
  return true;
  }

bool BlochSys::CheckCoords(const coord_vec& MxMyMzs, bool warn) const
  {
  int nm = Offsets.size();			// Number of vectors
  if(nm != MxMyMzs.size())
    {
    if(warn)
      {
      BSerror(10, 1);				// Dimensioning problem
      BSerror(16, 1);				// Improper # of coordinates
      }
    return false;
    }
  return true;
  }

// ____________________________________________________________________________
// iii                     Bloch System Setup Functions
// ____________________________________________________________________________

/* These are protected functions because they allow specific aspects of the
   spin system to be set up without worrying about system consistency! There     
   are two ways to set up a Bloch system: A.) By specific magnetization vectors
   & B.) By a spin network.  While the former is more direct, the latter fits
   into a more generic GAMMA sheme for MR based simulations.

   Magnetization vectors: 1.) Number of vectors
                          2.) Frequency offset for each vector
                          3.) R1 & R2 rates for each vector
                          4.) Kex rates for each vector pair

   Spin system:           1.) A GAMMA spin system 
                              a. Number of spins
                              b. Chemical shifts
                              c. Scalar couplings
                          2.) R1 & R2 rates for each spin 
                          3.) Kex rates for each spin pair

   Since GAMMA spin systems are often used in programs, method B will be used
   preferentially. This will be indicated by parameter NSpins in typical spin
   system fashion. Failing to find NSpins, method A will be used & parameter
   NMagVecs will should indicate how many vectors are defined in the system. */

bool BlochSys::SetSystem(const ParameterSet& pset, int idx, bool warn)
  {
//               Filter Out Parameters With Prefix [idx]

  ParameterSet  subpset;                        // Working parameter set
  if(idx != -1) subpset = pset.strip(idx);      // Get only params with [#]
  else          subpset = pset;                 // Or use full pset

//                    Zero Current Bloch System

  Offsets.clear();				// Set chemical shift
  R1rates.clear();				// Set R1 rate
  R2rates.clear();				// Set R2 rate
  isotopes.clear();				// Set isotope
  Krates.clear();				// No exchange rates
  Spins.clear();				// No vector spins
  _M = MagVec();				// Zero mag vector

//               Read Number Of Sub-Vectors In System

  int nm;					// For # spins or vectors
  if(!GetNVects(subpset, nm, false))		// If # vectors not found
    {						// look for # of spins
    if(!GetNSpins(subpset, nm, false))		//   If # spins not found
      {						//   we are in trouble
      if(warn) BSerror(56, 1);			//     Warn no # mag vects
      if(warn) BSerror(102, "NMagVects", 1);	//     Warn no # mag vects
      if(warn) BSerror(102, "NSpins", 1);	//     Warn no # mag vects
      if(warn) BSerror(40, 1);			//     Warn not enough params
      return false;				//     Return we failed
      }
    }

//                    Read Magnetization Vector

  if(!_M.SetSubVects(pset,nm,true))		// Try to read |M> in
    {						// If that fails
    if(warn) BSerror(40, 1);			//   Warn insufficient params
    return false;				//   Return we failed
    }

//                    Read System Vector Quatities

  if(!SetVects(subpset, nm, warn))		// Set up nm vectors
    {						// If that fails
    if(warn) BSerror(40, 1);			//   Warn insufficient params
    return false;				//   Return we failed
    }
  SetExchange(subpset, nm, false);		// Set nex exchange rates
  return true;					// Didn't set things up
  }

bool BlochSys::GetNSpins(const ParameterSet& pset, int& ns, bool warn) const
  {
  std::string pstate, pname("NSpins");		// Parameter name
  ParameterSet::const_iterator item;		// A pix into parameter list
  item = pset.seek(pname);                      // Pix in param. list for NSpins
  if(item != pset.end())                        // If it was found in the list
    {
    (*item).parse(pname,ns,pstate);		//   Retrieve # of spins
    return true;				//   Return we found it
    }
  if(warn) BSerror(102, pname, 1);
  return false;
  }

bool BlochSys::GetNVects(const ParameterSet& pset, int& nm, bool warn) const
  {
  std::string pstate, pname("NMagVecs");	// Parameter name
  ParameterSet::const_iterator item;		// A pix into parameter list
  item = pset.seek(pname);                      // Pix in param. list for NSpins
  if(item != pset.end())                        // If it was found in the list
    {
    (*item).parse(pname,nm,pstate);		//   Retrieve # of spins
    return true;				//   Return we found it
    }
  if(warn) BSerror(102, pname, 1);
  return false;
  }

// ----------------------------------------------------------------------------
//        Set Up The Bloch System By Magnetization Vector Designations         
// ----------------------------------------------------------------------------

bool BlochSys::SetVects(const ParameterSet& pset, int N, bool warn)
  {
  double  W;					// For offset value
  Isotope I;					// For isotope type
  double  R1;					// For R1 rate
  double  R2;					// For R2 rate
  int     Sp;					// For spin assigment
  for(int i=0; i<N; i++)
    {
    if(!GetVect(pset, i, W, I, R1, R2, Sp))
      {
// sosi
      return false;
      }
    Offsets.push_back(W);			// Set offset value
    R1rates.push_back(R1);			// Set R1 rate
    R2rates.push_back(R2);			// Set R2 rate
    isotopes.push_back(I);			// Set isotope
    Spins.push_back(Sp);			// Set spin
    }
  return true;
  }

bool BlochSys::GetVect(const ParameterSet& pset, int i, double& W, Isotope& I,
                   double& R1, double& R2, int& Sp, bool warn) const
  {
  bool TF;
  TF =       GetW(pset,   i,  W, warn);		// Get offset frequency
             GetIso(pset, i,  I, false);	// Get isotope type
  TF = (TF & GetR1(pset,  i, R1, warn));	// Get R1 relaxation rate
  TF = (TF & GetR2(pset,  i, R2, warn));	// Get R2 relaxation rate
             GetSp(pset,  i, Sp, false);	// Get spin assignment
  return TF;
  }

bool BlochSys::GetW(const ParameterSet& pset, int i, double& v, bool warn) const
  {
  std::string pname = "v(" + Gdec(i) + ")";	// Offset parameter name
  ParameterSet::const_iterator item;		// A pix into parameter list
  item = pset.seek(pname);			// Pix in parameter list for v(i)
  std::string pstate;				// Temp string
  if(item != pset.end())			// Retrieve the offset value
    {
    (*item).parse(pname,v,pstate);
    v *= HZ2RAD;
    return true;
    }
  pname = std::string("M") + pname;
  item = pset.seek(pname);			// Pix in parameter list for v(i)
  if(item != pset.end())			// Retrieve the offset value
    {
    (*item).parse(pname,v,pstate);
    v *= HZ2RAD;
    return true;
    }
  if(warn)
    {
    BSerror(50, 1);
    BSerror(102, pname, 1);
    }
  return false;
  }

bool BlochSys::GetIso(const ParameterSet& pset, int i, Isotope& I, bool warn) const
  {
  std::string pname = "Iso(" + Gdec(i) + ")";	// Isotope parameter name
  ParameterSet::const_iterator item; 		// A pix into parameter list
  item = pset.seek(pname);			// Pix in parameter list for Iso(i)
  std::string Ival, pstate;			// Temp strings
  if(item != pset.end())			// Retrieve the offset value
    {
    (*item).parse(pname,Ival,pstate);
    I = Isotope(Ival);
    return true;
    }
  pname = std::string("M") + pname;
  item = pset.seek(pname);			// Pix in parameter list for MIso(i)
  if(item != pset.end())			// Retrieve the offset value
    {
    (*item).parse(pname,Ival,pstate);
    I = Isotope(Ival);
    return true;
    }
  if(warn)
    {
    BSerror(51, 1);
    BSerror(102, pname, 1);
    }
  I = DEFISO;
  return false;
  }

bool BlochSys::GetR1(const ParameterSet& pset, int i, double& R1, bool warn) const
  {
  std::string pname = "R1(" + Gdec(i) + ")";	// Offset parameter name
  ParameterSet::const_iterator item; 		// A pix into parameter list
  item = pset.seek(pname);			// Pix in parameter list for R1(i)
  std::string pstate;				// Temp string
  if(item != pset.end())			// Retrieve the offset value
    {
    (*item).parse(pname,R1,pstate);
    return true;
    }
  pname = "T1(" + Gdec(i) + ")";		// Offset parameter name
  item = pset.seek(pname);			// Pix in parameter list for T1(i)
  if(item != pset.end())			// Retrieve the offset value
    {
    (*item).parse(pname,R1,pstate);
    R1 = 1.0/R1;
    return true;
    }
  if(warn)
    {
    BSerror(52, 1);
    BSerror(102, pname, 1);
    }
  return false;
  }

bool BlochSys::GetR2(const ParameterSet& pset, int i, double& R2, bool warn) const
  {
  std::string pname = "R2(" + Gdec(i) + ")";	// Offset parameter name
  ParameterSet::const_iterator item; 		// A pix into parameter list
  item = pset.seek(pname);			// Pix in parameter list for R2(i)
  std::string pstate, pname1;			// Temp string
  if(item != pset.end())			// Retrieve the offset value
    {
    (*item).parse(pname,R2,pstate);
    return true;
    }
  pname1 = "T2(" + Gdec(i) + ")";		// Offset parameter name
  item = pset.seek(pname1);			// Pix in parameter list for T2(i)
  if(item != pset.end())			// Retrieve the offset value
    {
    (*item).parse(pname1,R2,pstate);
    R2 = 1.0/R2;
    return true;
    }
  pname1 = "LW(" + Gdec(i) + ")";		// Offset parameter name
  item = pset.seek(pname1);			// Pix in parameter list for T2(i)
  if(item != pset.end())			// Retrieve the offset value
    {
    (*item).parse(pname1,R2,pstate);
    R2 *= PI;
    return true;
    }
  if(warn)
    {
    BSerror(53, 1);
    BSerror(102, pname, 1);
    }
  return false;
  }

bool BlochSys::GetSp(const ParameterSet& pset, int i, int& Sp, bool warn) const
  {
  std::string pname = "Spin(" + Gdec(i) + ")";	// Parameter name
  ParameterSet::const_iterator item; 		// A pix into parameter list
  item = pset.seek(pname);			// Pix in list for MCoord(i)
  std::string pstate, pname1;			// Temp string
  if(item != pset.end())			// Retrieve the offset value
    {
    (*item).parse(pname1,Sp,pstate);
    return true;
    }
  if(warn)
    {
    BSerror(55, 1);
    BSerror(102, pname, 1);
    }
  Sp = i;
  return false;
  }

bool BlochSys::SetExchange(const ParameterSet& pset, int nm, bool warn)	
  {
  int nex = (nm*nm-nm)/2;			// Number of exchange rates 
  Krates  = std::vector<double> (nex, 0);	// Begin with no exchange
  ParameterSet::const_iterator item; 		// A pix into parameter list
  std::string pnst("Kex(");			// Parameter name start
  std::string pnfi(")");			// Parameter name end
  std::string pname, pstate;
  double Kval;
  int i, j, k;					// Mag. vector indices
  for(i=0, k=0; i<nm-1; i++)			// Loop over vector pairs
    for(j=i+1; j<nm; j++, k++)
      {
      pname = pnst+Gdec(i)+","+Gdec(j)+")";	//   Parameter Kex(i,j)
      item = pset.seek(pname);			//   Pix in list for Kex(i)
      if(item != pset.end())			//   Retrieve exchange value
        {
        (*item).parse(pname,Kval,pstate);
         Krates[k] = Kval;
        }
      }
  return true;
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A                  SPIN SYSTEM CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                              Simple Constructors
// ----------------------------------------------------------------------------

BlochSys::BlochSys(int nm)
  {
  Isotope H("1H");				// Default isotope
  int nex  = (nm*nm-nm)/2;			// Number of exchange rates 
  Offsets  = std::vector<double>  (nm, 0);	// Array of shifts
  isotopes = std::vector<Isotope> (nm, H);	// Array of isotopes
  R1rates  = std::vector<double>  (nm, 0);	// Array of R1 rates
  R2rates  = std::vector<double>  (nm, 0);	// Array of R2 rates
  Krates   = std::vector<double>  (nex,0);	// Array of exchange rates
  Spins    = std::vector<int>     (nm, 0);	// Array of spins
  _M       = MagVec(nm);			// Magnetization vector
  }

BlochSys::BlochSys(const BlochSys& sys)
  {
  Offsets  = sys.Offsets;			// Copy mag. vector offsets
  isotopes = sys.isotopes;			// Copy mag. vector isotopes
  R1rates  = sys.R1rates;			// Copy mag. vector R1 rates
  R2rates  = sys.R2rates;			// Copy mag. vector R2 rates
  Krates   = sys.Krates;			// Copy mag. vector K rates
  Spins    = sys.Spins;				// Copy mag. vector spins
  _M       = sys._M;				// Copy nagnetization vector
  }

// ----------------------------------------------------------------------------
//                   Constructors Using Single Vector Arguments
// ----------------------------------------------------------------------------

BlochSys::BlochSys(double w, double R1, double R2) 
  {
  Offsets  = std::vector<double>(1, w*HZ2RAD);	// Set chemical shift
  R1rates  = std::vector<double>(1, R1);	// Set R1 rate
  R2rates  = std::vector<double>(1, R2);	// Set R2 rate
  Isotope H("1H");				// Default isotope
  isotopes = std::vector<Isotope>(1,  H);	// Set isotope
  Krates   = std::vector<double>();		// No exchange rates
  Spins    = std::vector<int>     (1, 0);	// Associated spin
  _M       = MagVec(1);				// Magnetization vector
  }

// ----------------------------------------------------------------------------
//                  Constructors Using Multiple Vector Arguments
// ----------------------------------------------------------------------------

BlochSys::BlochSys(const std::vector<double>& SH,
                          const std::vector<double>& R1s, const std::vector<double>& R2s) 
  {
  Offsets = SH;					// Set chemical shifts
  if(!CheckR1s(R1s)) BSfatal(30);		// Cant set R1 values
  if(!CheckR2s(R2s)) BSfatal(31);		// Cant set R2 values
  R1rates = R1s;				// Set R1 rates
  R2rates = R2s;				// Set R2 rates
  Isotope H("1H");				// Default isotope
  int nm  = Offsets.size();			// Number of mag. vectors
  int nex  = (nm*nm-nm)/2;			// Number of exchange rates 
  isotopes = std::vector<Isotope> (nm,  H);	// Array of isotopes
  Krates   = std::vector<double>  (nex, 0);	// Array of exchange rates
  Spins    = std::vector<int>     (nm,  0);	// Array of spins
  _M       = MagVec(nm);			// Magnetization vector
  }

BlochSys::BlochSys(const std::vector<double>& SH, const std::vector<Isotope>& Is,
                          const std::vector<double>& R1s, const std::vector<double>& R2s) 
  {
  Offsets = SH;					// Set chemical shifts
  if(!CheckR1s(R1s)) BSfatal(30);		// Cant set R1 values
  if(!CheckR2s(R2s)) BSfatal(31);		// Cant set R2 values
  if(!CheckIsos(Is)) BSfatal(32);		// Cant set isotope values
  isotopes = Is;				// Set isotopes
  R1rates  = R1s;				// Set R1 rates
  R2rates  = R2s;				// Set R2 rates
  int nm  = Offsets.size();			// Number of mag. vectors
  int nex  = (nm*nm-nm)/2;			// Number of exchange rates 
  Krates   = std::vector<double>  (nex, 0);	// Array of exchange rates
  Spins    = std::vector<int>     (nm,  0);	// Array of spins
  _M       = MagVec(nm);			// Magnetization vector
  }

BlochSys::BlochSys(const std::vector<double>& SH, const std::vector<double>& R1s, 
                           const std::vector<double>& R2s, const std::vector<double>& Ks) 
  {
  Offsets = SH;					// Set chemical shifts
  if(!CheckR1s(R1s)) BSfatal(30);		// Cant set R1 values
  if(!CheckR2s(R2s)) BSfatal(31);		// Cant set R2 values
  if(!CheckR2s(R2s)) BSfatal(33);		// Cant set exchange values
  R1rates  = R1s;				// Set R1 rates
  R2rates  = R2s;				// Set R2 rates
  Isotope H("1H");				// Default isotope
  int nm  = Offsets.size();			// Number of mag. vectors
  isotopes = std::vector<Isotope> (nm, H);	// Array of isotopes
  Krates   = Ks;				// Set K  rates
  Spins    = std::vector<int>     (nm, 0);	// Array of spins
  _M       = MagVec(nm);			// Magnetization vector
  }

BlochSys::BlochSys(const std::vector<double>& SH, const std::vector<Isotope>& Is,
const std::vector<double>& R1s, const std::vector<double>& R2s, const std::vector<double>& Ks) 
  {
  Offsets = SH;					// Set chemical shifts
  if(!CheckR1s(R1s)) BSfatal(30);		// Cant set R1 values
  if(!CheckR2s(R2s)) BSfatal(31);		// Cant set R2 values
  if(!CheckIsos(Is)) BSfatal(32);		// Cant set isotope values
  if(!CheckR2s(R2s)) BSfatal(33);		// Cant set exchange values
  isotopes = Is;				// Set isotopes
  R1rates  = R1s;				// Set R1 rates
  R2rates  = R2s;				// Set R2 rates
  Krates   = Ks;				// Set K  rates
  int nm = Offsets.size();			// Number of vectors
  Spins    = std::vector<int>     (nm, 0);		// Array of spins
  _M       = MagVec(nm);			// Magnetization vector
  }

// ----------------------------------------------------------------------------
//                     Constructors Using Spin Systems
// ----------------------------------------------------------------------------

/* The spin system will be assumed to be weakly coupled! In this manner it is
   easier to determine which transitions are associated with which spins.    */

BlochSys::BlochSys(const spin_system& sys, const RBasic& Rs)
  {
  int ns       = sys.spins();			// Number of spins
  if(!CheckSpins(ns, Rs.spins())) BSfatal(29);	// Check spin count matches
  gen_op H     = How(sys);			// Isotropic Hamiltonian
  gen_op sigma = Fx(sys);			// System total Fx operator
  gen_op D;					// Detection operator
  double mag;					// |M> intensity
  TTable1D TT;					// Transitions table
  int i,j,k;					// Dummy indices
  std::vector<MagVec> Ms;				// Vector of |M>s
  MagVec M;					// Single vector |M>
  int nc, NC=0;					// Components 
  for(i=0; i<ns; i++)				// Loop over all spins
    {
    D = Fm(sys, i);				//   Detect F- spin i
    acquire1D Acq(D, H);			//   Set up acquisition
    TT = Acq.table(sigma);			//   Find all transitions  
    nc = TT.size();				//   Number of transitions
    NC += nc;					//   Total # of transitions
    M = MagVec(nc);				//   Initialize |M>
    for(j=0; j<TT.size(); j++)			//   Loop transitions
      {
      Offsets.push_back(TT.Fr(j));		//     Set offset
      isotopes.push_back(sys.isotope(i));	//     Set spin type
      R1rates.push_back(Rs.R1(i));		//     Set R1 rate
      R2rates.push_back(Rs.R2(i));		//     Set R2 rate
      mag = norm(TT.I(j));			//     Set norm
      M.put(mag, 3*j+2);			//     Sub-vector @ (0,0,mag)
      Spins.push_back(i);			//     Set spin
      }
    Ms.push_back(M);				// Store |M> this spin
    }
  _M = MagVec(NC);				// Concatonated |M>
  for(i=0, k=0; i<ns; i++)
    {
    nc = Ms[i].NComps();
    for(j=0; j<nc; j++, k++)
      _M.put(Ms[i].Mz(j), 3*k+2);
    }
  int nm  = Offsets.size();			// Number of mag. vectors
  int nex = (nm*nm-nm)/2;			// Number of exchange rates 
  Krates  = std::vector<double>  (nex, 0);		// Array of exchange rates
  }


BlochSys::BlochSys(const spin_system& sys, const matrix& Ks)
  {
  int ns       = sys.spins();			// Number of spins
  if(!CheckSpins(ns, Ks.rows())) BSfatal(33);	// Check spin count matches
  if(!CheckSpins(ns, Ks.cols())) BSfatal(33);	// Check spin count matches
  gen_op H     = How(sys);			// Isotropic Hamiltonian
  gen_op sigma = Fx(sys);			// System total Fx operator
  gen_op D;					// Detection operator
  TTable1D TT;					// Transitions table
  int i,j;					// Dummy indices
  double mag;
  for(i=0; i<ns; i++)				// Loop over all spisn
    {
    D = Fm(sys, i);				//   Detect F- spin i
    acquire1D Acq(D, H);			//   Set up acquisition
    TT = Acq.table(sigma);			//   Find all transitions  
    for(j=0; j<TT.size(); j++)
      {
      Offsets.push_back(TT.Fr(j));		//     Set offset
      isotopes.push_back(sys.isotope(i));	//     Set spin type
      R1rates.push_back(0.0);			//     Set R1 rate
      R2rates.push_back(0.0);			//     Set R2 rate
      mag = norm(TT.I(j));
// sosi
//      Vectors.push_back(coord(0.,0.,mag));	//     Set vector coords
      Spins.push_back(i);			//     Set spin
      }
    }
  int nm  = Offsets.size();			// Number of mag. vectors
  int I,J;					// Spin indices
  for(i=0; i<nm-1; i++)				// Loop over all vector pairs
    {						//   Get spin assoicated with
    I = Spins[i];				//   magnetization vector i
    for(j=i+1; j<nm; j++)
      {						//   Get spin associated with
      J = Spins[j];				//   magnetization vector j
      Krates.push_back(Ks.getRe(I,J));		//   Store exchange rate
      }
    }
  }

BlochSys::BlochSys(const spin_system& sys, const RBasic& Rs, const matrix& Ks)
  {
  int ns       = sys.spins();			// Number of spins
  if(!CheckSpins(ns, Rs.spins())) BSfatal(29);	// Check spin count matches
  if(!CheckSpins(ns, Ks.rows()))  BSfatal(33);	// Check spin count matches
  if(!CheckSpins(ns, Ks.cols()))  BSfatal(33);	// Check spin count matches
  gen_op H     = How(sys);			// Isotropic Hamiltonian
  gen_op sigma = Fx(sys);			// System total Fx operator
  gen_op D;					// Detection operator
  TTable1D TT;					// Transitions table
  int i,j;					// Dummy indices
  double mag;
  for(i=0; i<ns; i++)				// Loop over all spisn
    {
    D = Fm(sys, i);				//   Detect F- spin i
    acquire1D Acq(D, H);			//   Set up acquisition
    TT = Acq.table(sigma);			//   Find all transitions  
    for(j=0; j<TT.size(); j++)
      {
      Offsets.push_back(TT.Fr(j));		//     Set offset
      isotopes.push_back(sys.isotope(i));	//     Set spin type
      R1rates.push_back(Rs.R1(i));		//     Set R1 rate
      R2rates.push_back(Rs.R2(i));		//     Set R2 rate
      mag = norm(TT.I(j));
// sosi
//      Vectors.push_back(coord(0.,0.,mag));	//     Set vector coords
      Spins.push_back(i);			//     Set spin
      }
    }
  int nm  = Offsets.size();			// Number of mag. vectors
  int I,J;					// Spin indices
  for(i=0; i<nm-1; i++)				// Loop over all vector pairs
    {						//   Get spin assoicated with
    I = Spins[i];				//   magnetization vector i
    for(j=i+1; j<nm; j++)
      {						//   Get spin associated with
      J = Spins[j];				//   magnetization vector j
      Krates.push_back(Ks.getRe(I,J));		//   Store exchange rate
      }
    }
  }

// ----------------------------------------------------------------------------
//                     Constructors Using Other Objects
// ----------------------------------------------------------------------------

BlochSys::BlochSys(const TTable1D& TT, const std::string& Iso)
  {
  Isotope I(Iso);
  double mag;
  for(int j=0; j<TT.size(); j++)
    {
    Offsets.push_back(TT.Fr(j));		//     Set offset
    isotopes.push_back(I);			//     Set spin type
    R1rates.push_back(TT.R2(j));		//     Set R1 rate
    R2rates.push_back(TT.R2(j));		//     Set R2 rate
    mag = norm(TT.I(j));
// sosi
//    Vectors.push_back(coord(0.,0.,mag));	//     Set vector coords
    Spins.push_back(0);				//     Set spin
    }
  int nm  = Offsets.size();			// Number of mag. vectors
  int nex  = (nm*nm-nm)/2;			// Number of exchange rates 
  Krates   = std::vector<double>  (nex, 0);		// Array of exchange rates
  }

// ----------------------------------------------------------------------------
//                          Assignment and Destruction
// ----------------------------------------------------------------------------

BlochSys& BlochSys::operator= (const BlochSys& sys)
  {
  if(this == &sys) return *this;		// Return if self-assign
  Offsets  = sys.Offsets;			// Copy mag. vector offsets
  isotopes = sys.isotopes;			// Copy mag. vector isotopes
  R1rates  = sys.R1rates;			// Copy mag. vector R1 rates
  R2rates  = sys.R2rates;			// Copy mag. vector R2 rates
  Krates   = sys.Krates;			// Copy mag. vector K rates
  Spins    = sys.Spins;				// Copy mag. vector spins
  _M       = sys._M;				// Copy mag. vector
  return *this;
  }

BlochSys::~BlochSys () { }

// ____________________________________________________________________________
// B                 Magnetization Vector Offset Access
// ____________________________________________________________________________

/*
double BlochSys::Offset(int i) const
  {
  return Offsets[i];				// Return offset in 1/sec
  }

void BlochSys::Offset(double w, int i)
  {
  Offsets[i] = w;
  }
*/

// ____________________________________________________________________________
// C              Magnetization Vector Isotope & Spin Access
// ____________________________________________________________________________

int BlochSys::NIso() const
  {
  int nm = isotopes.size();			// Number of vectors
  int NI = 0;					// No isotopes so far
  if(!nm) return NI;				// If none, no isotopes
  std::vector<Isotope> UI;				// Unique isotopes
  UI.push_back(isotopes[0]);			// 1st one is unique!
  NI++; 					// We have at least one
  int i, j;					// Temp indices
  bool unique = true;				// Flag for unique isotope
  Isotope I;					// Temp isotope
  for(i=1; i<nm; i++)				// Loop remaining vectors
    {
    I = isotopes[i];				//   Get the isotope type
    unique = true;				//   Assume it is unique
    for(j=0; j<NI && unique; j++)		//   Compare with others
      { if(I == UI[j]) unique=false; }		//   to see if really unique
    if(unique)					//   If it is indeed unique
      { UI.push_back(isotopes[nm]); NI++; }	//   store in list & count
    }
  return NI;					// Return number of isotopes
  }						// (that are unique)

int BlochSys::IsoMaxLength() const
  {
  int l, len=0;					// Start with no length
  int nm = isotopes.size();			// Number of vectors
  for(int i=0; i<nm; i++)			// Loop magnetization vectors
    {
    l = ((isotopes[i]).symbol()).length();	//   Length of i vector symbol
    len = gmax(len, l);				//   Set for maximum
    }
  return len;
  }

int BlochSys::NSpins() const
  {
  int nm = size();				// Number of vectors
  int NS = 0;					// No spins so far
  if(!nm) return NS;				// If none, no spins
  std::vector<int> US;				// Unique spins
  US.push_back(Spins[0]);			// 1st one is unique!
  NS++; 					// We have at least one
  int i, j;					// Temp indices
  bool unique = true;				// Flag for unique spin
  int I;					// Temp spin
  for(i=1; i<nm; i++)				// Loop remaining vectors
    {
    I = Spins[i];				//   Get the spin type
    unique = true;				//   Assume it is unique
    for(j=0; j<NS && unique; j++)		//   Compare with others
      { if(I == US[j]) unique=false; }		//   to see if really unique
    if(unique)					//   If it is indeed unique
      { US.push_back(Spins[nm]); NS++; }	//   store in list & count
    }
  return NS;					// Return number of spins
  }						// (that are unique)


// ____________________________________________________________________________
// D              Magnetization Vector Relaxation Rate Access
// ____________________________________________________________________________

// ____________________________________________________________________________
// E                Magnetization Vector Exchange Rate Access
// ____________________________________________________________________________

double BlochSys::R1(int i) const { return R1rates[i];     }
double BlochSys::T1(int i) const { return 1.0/R1rates[i]; }
double BlochSys::R2(int i) const { return R2rates[i];     }
double BlochSys::T2(int i) const { return 1.0/R2rates[i]; }
double BlochSys::LW(int i) const { return R2rates[i]/PI;  }

double BlochSys::MaxExchange() const
  {
  double K = 0;				// Begin with 0 exchange rate
  int nm = Offsets.size();		// Number of vectors defined
  int nex = (nm*nm - nm)/2;		// Number of exchange rates
  for(int i=0; i<nex; i++)
    K = gmax(K, Krates[i]);
  return K;
  }

// ____________________________________________________________________________
// F                Magnetization Vector Component Access
// ____________________________________________________________________________

std::vector<double> BlochSys::Norms() const                        { return _M.Norms(); }
void           BlochSys::Norms(const std::vector<double>& Ns) { _M.Norms(Ns); }

double BlochSys::Norm(int i) const { return _M.Norm(i); }
void   BlochSys::Norm(double nv, int i) { _M.Norm(nv,i); } 

// ____________________________________________________________________________
// G                         Bloch Equation Arrays
// ____________________________________________________________________________

/* For each magnetization vector we have an associated 3x3 block in these 
   array.  In the case of the H array (B fields) and R array (relaxation) the
   blocks appear as follows.

                                                        
       [        0          -w + w       g*B  * sin(phi) ]       [ R   0   0  ]
       |                     0   rf        1            |       |  2         |
       |                                                |       |            |
   H = |     w - w             0       -g*B  * cos(phi) |   R = | 0   R   0  |
       |      0   rf                       1            |       |      2     |
       |                                                |       |            |
       [ -gB * sin(phi)   gB * cos(phi)       0         |       | 0   0   R  |
       [    1               1                           ]       [          1 ]
            
   We used the following definitions. The rf-field strength is given by 
   gB1 = gamma*B1. The field offset is wrf and the field phase is phi. The
   vector offset is w0 and the relaxation times are R1 = 1/T1, R2=1/T2.      */

// sosi - remove these if not used at all....

matrix BlochSys::H() const { return B(); }
matrix BlochSys::H(double gamB1, double wrf, double phi) const { return B(gamB1, wrf, phi); }

matrix BlochSys::B() const
  {
  int nm = Offsets.size();				// Number of vectors
  matrix Hmx(nm*3, nm*3, complex0);			// Initial matrix
  double w;						// Rot. frame frequency
  int ix3 = 0;						// For base block index
  for(int i=0; i<nm; i++, ix3+=3)			// Loop over each vector
    {							// & fill its 3x3 block
    w = Offsets[i];					//   Vector offset
    Hmx.put( w, ix3,   ix3+1);				//   <0|H|1>=  w
    Hmx.put(-w, ix3+1, ix3  );				//   <1|H|0>= -w
    }
  return Hmx;
  }

matrix BlochSys::B(double gamB1, double wrf, double phi) const
  {
  int nm = Offsets.size();				// Number of vectors
  matrix Hmx(nm*3, nm*3, complex0);			// Initial matrix
  double delw;						// Rot. frame frequency
  double sinphi = sin(phi);				// RF-field sin(phase)
  double cosphi = cos(phi);				// RF-field cos(phase)
  int ix3 = 0;						// For base block index
  double gB1 = gamB1*HZ2RAD;				// Set rf strengh 1/s
  wrf *= HZ2RAD;					// Set rf freqeuncy in 1/s
  for(int i=0; i<nm; i++, ix3+=3)			// Loop over each vector
    {							// & fill its 3x3 block
    delw = Offsets[i]-wrf;				//   Rot. frame offset
    Hmx.put( delw,       ix3,   ix3+1);			//   <0|H|1>= delw
    Hmx.put( gB1*sinphi, ix3,   ix3+2);			//   <0|H|2>= gB1*S(phi)
    Hmx.put(-delw,       ix3+1, ix3  );			//   <1|H|0>=-delw
    Hmx.put(-gB1*cosphi, ix3+1, ix3+2);			//   <1|H|2>=-gB1*C(phi)
    Hmx.put(-gB1*sinphi, ix3+2, ix3  );			//   <2|H|0>=-gB1*S(phi)
    Hmx.put( gB1*cosphi, ix3+2, ix3+1);			//   <2|H|1>= gB1*C(phi)
    }
  return Hmx;
  }

matrix BlochSys::R() const
  {
  int nm = Offsets.size();				// Number of vectors
  matrix Rmx(nm*3, nm*3, d_matrix_type);		// Initial matrix
  int ix3 = 0;						// For base block index
  for(int i=0; i<nm; i++, ix3+=3)			// Loop over each vector
    {							// & fill its 3x3 block
    Rmx.put(R2rates[i], ix3,   ix3  );			//   <0|R|0> = R2
    Rmx.put(R2rates[i], ix3+1, ix3+1);			//   <1|R|1> = R2
    Rmx.put(R1rates[i], ix3+2, ix3+2);			//   <2|R|2> = R1
    }
  return Rmx;
  }

matrix BlochSys::K() const
  {
  int nm = Offsets.size();				// Number of vectors
  matrix Kmx(nm*3, nm*3, complex0, h_matrix_type);	// Exchange matrix
  int i,j,k;						// Dummy indices
  for(i=0, k=0; i<nm-1; i++)				// Loop first vectors
    for(j=i+1; j<nm; j++,k++)				// Loop second vectors
      {
      Kmx.put(Krates[k],    3*i,  3*i);		// Diagonal <Mx1|Mx1>
      Kmx.put(Krates[k],    3*i+1,3*i+1);	// Diagonal <My1|My1>
      Kmx.put(Krates[k],    3*i+2,3*i+2);	// Diagonal <Mz1|Mz1>
      Kmx.put_h(-Krates[k], 3*i,  3*j);		// OffDiag  <Mx1|Mx2>
      Kmx.put_h(-Krates[k], 3*i+1,3*j+1);	// OffDiag  <My1|My2>
      Kmx.put_h(-Krates[k], 3*i+2,3*j+2);	// OffDiag  <Mz1|Mz2>
      Kmx.put(Krates[k],    3*j,  3*j);		// Diagonal <Mx2|Mx2>
      Kmx.put(Krates[k],    3*j+1,3*j+1);	// Diagonal <My2|My2>
      Kmx.put(Krates[k],    3*j+2,3*j+2);	// Diagonal <Mz2|Mz2>
      }
  return Kmx;
  }

// ____________________________________________________________________________
// H                     Bloch Equation Column Vectors
// ____________________________________________________________________________

/* These will magnetization vectors. Detection operators will be row vectors */

      MagVec  BlochSys::Meq() const { return MagVec(size()); }
const MagVec& BlochSys::Mo()  const { return _M; }

MagVec BlochSys::Mx() const
  {
  int nm = size();				// Number of vectors
  MagVec M(nm);					// Empty vector
  for(int i=0, j=0; i<nm; i++)			// Loop over vector elements
    {M.put(1,j++); M.put(0,j++); M.put(0,j++);}
  return M;
  }

MagVec BlochSys::My() const
  {
  int nm = size();				// Number of vectors
  MagVec M(nm);					// Empty vector
  for(int i=0, j=0; i<nm; i++)			// Loop over vector elements
    {M.put(0,j++); M.put(1,j++); M.put(0,j++);}
  return M;
  }

MagVec BlochSys::Mz() const
  {
  int nm = size();				// Number of vectors
  MagVec M(nm);					// Empty vector
  for(int i=0, j=0; i<nm; i++)			// Loop over vector elements
    {M.put(0,j++); M.put(0,j++); M.put(1,j++);}
  return M;
  }

MagVec BlochSys::Mss(const matrix& L, const matrix& R) const
  { return inv(L)*R*Meq(); }

MagVec BlochSys::Mss(const matrix& L, const matrix& R, const col_vector& Meq) const
  { return inv(L)*R*Meq; }

// ____________________________________________________________________________
// I                       Bloch Equation Row Vectors
// ____________________________________________________________________________

row_vector BlochSys::DetectMu(int u) const
  {
  int nm = size();				// Number of vectors
  row_vector D(3*nm, complex0);			// Vector for Mx detection
  int jst = 0;
  switch(u)
    {
    case 0:  jst = 0; break;
    case 1:  jst = 1; break;
    case 2:  jst = 2; break;
    default: jst = 0; break;
    }
  for(int i=0, j=jst; i<nm; i++,j+=3)		// Loop over Mxi positions
    D.put(complex1, j);				// & set vector to detect
  return D;
  }

row_vector BlochSys::DetectMu(int k, int u) const
  {
  int nm = Offsets.size();			// Number of vectors
  row_vector D(3*nm, complex0);			// Vector for Mx detection
  int jst = 0;
  switch(u)
    {
    case 0:  jst = 0; break;
    case 1:  jst = 1; break;
    case 2:  jst = 2; break;
    default: jst = 0; break;
    }
  for(int i=0, j=jst; i<nm; i++, j+= 3)			// Loop over Mxi positions
    if(Spins[i] == k) D.put(complex1, j);	// & set vector to detect
  return D;					// Only if Mxi assocated 
  }						// with spin k though

row_vector BlochSys::DetectMu(const std::string& I, int u) const
  {
  int nm = Offsets.size();			// Number of vectors
  int bd = 3*nm;				// Bloch dimension
  row_vector D(bd, complex0);			// Vector for Mx detection
  Isotope II(I);
  int ist = 0;
  switch(u)
    {
    case 0:  ist = 0; break;
    case 1:  ist = 1; break;
    case 2:  ist = 2; break;
    default: ist = 0; break;
    }
  for(int i=ist; i<bd; i+= 3)			// Loop over Mxi positions
    if(isotopes[i] == II) D.put(complex1, i);	// & set vector to detect
  return D;					// Only if Mxi assocated 
  }						// with spin k though

row_vector BlochSys::DetectMx()                const { return DetectMu(0); }
row_vector BlochSys::DetectMy()                const { return DetectMu(1); }
row_vector BlochSys::DetectMz()                const { return DetectMu(2); }

row_vector BlochSys::DetectMx(int i)           const { return DetectMu(i, 0); }
row_vector BlochSys::DetectMy(int i)           const { return DetectMu(i, 1); }
row_vector BlochSys::DetectMz(int i)           const { return DetectMu(i, 2); }

row_vector BlochSys::DetectMx(const std::string& I) const { return DetectMu(I, 0); }
row_vector BlochSys::DetectMy(const std::string& I) const { return DetectMu(I, 1); }
row_vector BlochSys::DetectMz(const std::string& I) const { return DetectMu(I, 2); }

// ____________________________________________________________________________
// J                    Bloch System Auxiliary Functions
// ____________________________________________________________________________

int BlochSys::size() const { return Offsets.size(); }

// ____________________________________________________________________________
// K                  Bloch System Parameter Set Functions
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//            Functions To Make A Spin System From A Parameter Set
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// L                     Bloch System Input Functions
// ____________________________________________________________________________

bool BlochSys::read(const std::string& filein, int idx, int warn)
  {
  ParameterSet pset;                    // Declare a parameter set
  if(!pset.read(filein, warn?1:0))      // Read in pset from file
    {
    if(warn)
      {
      BSerror(40, filein, 1);		// Problems with file filein
      if(warn>1) BSfatal(101,filein);	// Can't read from file filein
      else       BSerror(101,1);
      }
    return false;
    }
  return read(pset, idx, warn);
  }

bool BlochSys::read(const ParameterSet& pset,int idx, int warn)
  {
  if(SetSystem(pset, idx, warn?true:false))
    return true;
  if(warn)                                      // Looks like we can't read
    {                                           // the system at all
    std::string sl("");				//   String for index
    if(idx != -1) sl = std::string(Gdec(idx));	//   Set it if one exists
                  BSerror(77, sl, 1);
    if(warn > 1)  BSfatal(5);
    else          BSerror(78, sl);
    }
  return false;
  }

std::string BlochSys::ask_read(int argc, char* argv[], int argn)
  {
  std::string filename;                              // Name of spin system file
  query_parameter(argc, argv, argn,             // Get filename from command
       "\n\tBloch system filename? ", filename); // Or ask for it
  read(filename);                               // Read system from filename
  return filename;
  }

std::string BlochSys::ask_read(int argc, char* argv[], int argn, const std::string& def)
  {
  std::string msg = "\n\tBloch system filename ["	// Query we will ask if
             + def + "]? ";                     // it is needed
  std::string filename = def;                        // Name of spin system file
  ask_set(argc,argv,argn,msg,filename);         // Or ask for it
  read(filename);                               // Read system from filename
  return filename;                              // Return filename
  }

 
// ____________________________________________________________________________
// M                        SYSTEM STANDARD I/O FUNCTIONS
// ____________________________________________________________________________

        // Input                sys	: A Bloch system (this)
        //                      ostr    : Output stream
        // Output               none    : Spin system info is sent
        //                                into the output stream

std::ostream& BlochSys::print(std::ostream& ostr) const
  {
  int ns = Offsets.size();			// Get number mag. vectors
  std::string hdr("Bloch Spin System");		// Output header
  int hl;					// Header length
  if(!ns)					// If there are no vectors
    {						// just write system is empty
    hdr = "Empty " + hdr;
    hl = hdr.length();
    ostr << "\n" << std::string(40-hl/2, ' ') << hdr;
    return ostr;
    }

// ----------------------------------------------------------------------------
//                Output The Vector Offsets, Relaxation Rates
// ----------------------------------------------------------------------------

  hl = hdr.length();				// Main header length
  ostr << "\n" << std::string(40-hl/2,' ') << hdr;	// Output main header
  ostr << "\n" << std::string(40-hl/2,' ')		// Output main header underline
       << std::string(hdr.length(), '_');

  hdr = "Magnetization Vectors";		// Magnetization vector header
  hl = hdr.length();				// Mag. vector header length
  ostr << "\n\n";
  ostr << std::string(40-hl/2,' ') << hdr;		// Output mag vec. header
  int    cgap = 1;				// Gap between columns
  std::string cspc = std::string(cgap, ' ');		// String for gap
  int NI = NIso();
  std::vector<std::string> chdrs;				// Column headers
  chdrs.push_back(std::string("|M>"));		// 1st column header
  if(NI > 1)
      chdrs.push_back(std::string("Types"));		// 2st column header
  chdrs.push_back(std::string("  Offsets  "));	// 3rd column header
  chdrs.push_back(std::string(" R1 Rates "));	// 4th column header
  chdrs.push_back(std::string(" R2 Rates "));	// 5th column header
  chdrs.push_back(std::string(" T1 Times "));	// 6th column header
  chdrs.push_back(std::string(" T2 Times "));	// 7th column header
  chdrs.push_back(std::string(" Linewidths "));	// 8th column header
  int nc = chdrs.size();			// Number of coluns
  int i, llen=0;				// Column index, line length
  for(i=0; i<nc-1; i++)			// Determine line length
    llen += chdrs[i].size() + cgap;
  llen += chdrs[i].size();
  std::string spc = std::string(40-llen/2, ' ');	// Spacer to center line
  ostr << "\n\n" << spc;			// Start of header line
  for (i=0; i<nc; i++)
    {
    ostr << chdrs[i];
    if(i<nc-1) ostr << cspc;
    }
  ostr << "\n" << spc;			// Start of header underline
  for (i=0; i<nc; i++)
    {
    ostr << std::string(chdrs[i].length(), '-');
    if(i<nc-1) ostr << cspc;
    }
  int vlen = 1;				// Set spacing on 1st column
  if(ns>10)   vlen +=1;
  if(ns>100)  vlen +=1;
  if(ns>1000) vlen +=1;
  std::string vst((chdrs[0].length()-vlen)/2, ' ');
  std::string vfi(chdrs[0].length()-vlen-vst.length(), ' ');

    int tlen = IsoMaxLength();			// Set spacing on 2nd column
  std::string tst((chdrs[1].length()-tlen)/2, ' ');
  std::string tfi(chdrs[1].length()-tlen-tst.length(), ' ');

  double abss, ss;
  double R1, R2;
  double T1, T2, LW;
  for(i=0; i<ns; i++)
    {
    ostr << "\n" << spc;
    ostr << vst << Gdec(i, vlen)        << vfi << cspc;
    if(NI > 1)
      ostr << tst << isotopes[i].symbol() << tfi << cspc;

    ss = Offsets[i]*RAD2HZ;
    abss = fabs(ss);
    if(abss > 1.e9)       ostr << Gform("%7.2f", ss*1.e-9) << " GHz" << cspc;
    else if(abss > 1.e6)  ostr << Gform("%7.2f", ss*1.e-6) << " MHz" << cspc;
    else if(abss > 1.e3)  ostr << Gform("%7.2f", ss*1.e-3) << " KHz" << cspc;
    else                  ostr << Gform("%7.2f", ss)       << " Hz " << cspc;

    R1 = R1rates[i];
    if(fabs(R1) > 1.e7)      ostr << Gform("%7.3f/ns", R1*1.e-9);
    else if(fabs(R1) > 1.e6) ostr << Gform("%7.3f/us", R1*1.e-6);
    else if(fabs(R1) > 1.e3) ostr << Gform("%7.3f/ms", R1*1.e-3);
    else                     ostr << Gform("%7.3f/s ", R1);
    ostr << cspc;

    R2 = R2rates[i];
    if(fabs(R2) > 1.e9)      ostr << Gform("%7.3f/ns", R2*1.e-9);
    else if(fabs(R2) > 1.e6) ostr << Gform("%7.3f/us", R2*1.e-6);
    else if(fabs(R2) > 1.e3) ostr << Gform("%7.3f/ms", R2*1.e-3);
    else                     ostr << Gform("%7.3f/s ", R2);
    ostr << cspc;

    if(R1rates[i])
      {
      T1 = 1.0/R1;
      if(fabs(T1) < 1.e-9)      ostr << Gform("%7.3f ns", T1*1.e9);
      else if(fabs(T1) < 1.e-6) ostr << Gform("%7.3f us", T1*1.e6);
      else if(fabs(T1) < 1.e-3) ostr << Gform("%7.3f ms", T1*1.e3);
      else                      ostr << Gform("%7.3f s ", T1);
      ostr << cspc;
      }
    else
      ostr << std::string(chdrs[5].length(), ' ') << cspc;
    if(R2rates[i])
      {
      T2 = 1.0/R2;
      if(fabs(T2) < 1.e-9)      ostr << Gform("%7.3f ns", T2*1.e9);
      else if(fabs(T2) < 1.e-6) ostr << Gform("%7.3f us", T2*1.e6);
      else if(fabs(T2) < 1.e-3) ostr << Gform("%7.3f ms", T2*1.e3);
      else                      ostr << Gform("%7.3f s ", T2);
      ostr << cspc;
      }
    else
      ostr << std::string(chdrs[6].length(), ' ') << cspc;
    LW = fabs(R2rates[i]/PI);
    if(fabs(LW) > 1.e11)     ostr << Gform("%7.3f GHz", LW*1.e-9);
    else if(fabs(LW) > 1.e8) ostr << Gform("%7.3f MHz", LW*1.e-6);
    else if(fabs(LW) > 1.e5) ostr << Gform("%7.3f KHz", LW*1.e-3);
    else                     ostr << Gform("%7.3f Hz ", LW);
    }

// ----------------------------------------------------------------------------
//                       Output The Vector Coordinates
// ----------------------------------------------------------------------------

  hdr = "Vector Components";			// Coordinates header
  hl = hdr.length();				// Exchange header length
  ostr << "\n\n\n";				// Begin with new lines
  ostr << std::string(40-hl/2,' ') << hdr;		// Output exchange header
  int NS = NSpins();
  chdrs.clear();				// Column headers
  chdrs.push_back(std::string("|M>"));		// 1st column header
  if(NS > 1)
    chdrs.push_back(std::string("Spins"));		// 2st column header
  chdrs.push_back(std::string("    Mx    "));	// 3rd column header
  chdrs.push_back(std::string("    My    "));	// 4th column header
  chdrs.push_back(std::string("    Mz    "));	// 5th column header
  chdrs.push_back(std::string("   |M>|   "));	// 6th column header
  nc = chdrs.size();				// Number of columns
  llen=0;					// Line length
  for(i=0; i<nc-1; i++)				// Determine |M> col width
    llen += chdrs[i].size() + cgap;
  llen += chdrs[i].size();
  spc = std::string(40-llen/2, ' ');			// Spacer to center line
  ostr << "\n\n" << spc;			// Start of header line
  for (i=0; i<nc; i++)
    {
    ostr << chdrs[i];
    if(i<nc-1) ostr << cspc;
    }
  ostr << "\n" << spc;				// Start of header underline
  for (i=0; i<nc; i++)
    {
    ostr << std::string(chdrs[i].length(), '-');
    if(i<nc-1) ostr << cspc;
    }
  tlen = Gdec(NS).length();			// Set spacing on 2nd column
  tst = std::string((chdrs[1].length()-tlen)/2, ' ');
  tfi = std::string(chdrs[1].length()-tlen-tst.length(), ' ');


  for(i=0; i<ns; i++)
    {
    ostr << "\n" << spc;
    ostr << vst << Gdec(i, vlen)        << vfi << cspc;
    if(NS > 1)
      ostr << tst << Spins[i] << tfi << cspc;

    ostr << Gform("%10.4f", _M.Mx(i))   << cspc;
    ostr << Gform("%10.4f", _M.My(i))   << cspc;
    ostr << Gform("%10.4f", _M.Mz(i))   << cspc;
    ostr << Gform("%10.4f", _M.Norm(i)) << cspc;
    }

// ----------------------------------------------------------------------------
//                   Output The Vector Exchange Rates (If Any)
// ----------------------------------------------------------------------------

  if(MaxExchange())					// If any exchange processes
    {							// defined we output them too
    hdr = "Exchange Rates (1/sec)";			// Exchange rate header
    hl = hdr.length();					// Exchange header length
    ostr << "\n\n" << std::string(40-hl/2,' ')<<hdr << "\n";	// Output exchange header

    std::vector<std::string> chdrs;				// Column headers
    chdrs.push_back(std::string("|M>"));			// 1st column header
    chdrs.push_back(std::string(" | "));			// 2nd column header

    int    cw   = 8;					// Column width
    int    cgap = 1;					// Gap between columns
    int    llen = (ns-1)*cw + (ns-2)*cgap
                + chdrs[0].length() + chdrs[1].length();
    std::string spc(40-llen/2, ' ');
    std::string cspc = std::string(cgap, ' ');			// String for gap
    ostr << "\n" << spc << chdrs[0] << chdrs[1];	// 1st column header
    int i, j, k;
    for(i=1; i<ns; i++)
      {
      ostr << Gdec(i, cw);
      if(i < ns-1) ostr << cspc;
      }
    ostr << "\n" << spc << std::string("--- | ");		// 1st column header
    for(i=1; i<ns; i++)
      {
      ostr << std::string(cw, '-');
      if(i < ns-1) ostr << std::string(cgap, ' ');
      }
    for(i=0, k=0; i<ns-1; i++)
      {
      ostr << "\n" << spc << Gdec(i, chdrs[0].length());
      ostr << chdrs[1];
      for(j=1; j<=i; j++)
        ostr << std::string(cw, ' ') << cspc;
      for(j=i+1; j<ns; j++, k++)
        {
        ostr << Gform("%8.2f", Krates[k]);
        if(j < ns-1) ostr << cspc;
        }
      
      }
    }
  return ostr;
  }

std::ostream& operator<< (std::ostream& ostr, const BlochSys& sys)
  { return sys.print(ostr); }
    
#endif								// BlochSys.cc
