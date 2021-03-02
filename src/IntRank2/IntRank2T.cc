/* IntRank2T.cc *************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Rank 2 Interaction Spin Tensors            Implementation	**
**									**
**      Scott Smith                                                     **
**      Copyright © 1997                                                **
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
** A variable of type IntRank2T represents the spin tensor associated 	**
** with a rank 2 interaction defined for a particular nuclear spin or   **
** nuclear spin pair.  The spin tensor components are stored in         **
** in irreducible spherical format.                                     **
**                                                                      **
*************************************************************************/

#ifndef   IntRank2T_cc_			// Is file already included?
#  define IntRank2T_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <IntRank2/IntRank2T.h>		// Include interface definition
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Matrix/matrix.h>		// Include matrices
#include <Basics/Isotope.h>		// Include isotopes
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <HSLib/SpinOpSng.h>		// Include single spin operators
#include <Basics/StringCut.h>		// Include GAMMA string parsing
#include <ESRLib/CubicIon.h>		// Include cubic ion electrons
#if defined(_MSC_VER) || defined(__SUNPRO_CC)				// If we are using MSVC++
 #pragma warning (disable : 4786)       //   Kill STL namelength warnings
#endif

using std::string;			// Using libstdc++ strings
using std::list;			// Using libstdc++ STL lists
using std::cout;			// Using libstdc++ standard output
using std::vector;			// Using libstdc++ STL vectors
using std::ostream;			// Using libstdc++ output streams

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i             CLASS INTERACTION SPIN TENSOR ERROR HANDLING
// ____________________________________________________________________________

/*       Input                IST     : Rank 2 spin tensor (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pn      : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */
 
void IntRank2T::ISTerror(int eidx, int noret) const
  {
  string hdr("Rank 2 Spin Tensor");
  switch(eidx)
    {
    case 0: GAMMAerror(hdr,"Program Aborting.....",        noret); break;// (0)
    case 2: GAMMAerror(hdr,"Problems During Construction", noret); break;// (2)
    case 8: GAMMAerror(hdr,"Unknown Spin Interaction Type",noret); break;// (8)
    case 9: GAMMAerror(hdr,"Unknown Spin Pair Interaction",noret); break;// (9)
    case 10:GAMMAerror(hdr,"I Must Be A Multiple Of 1/2",  noret); break;// (10)
    case 11:GAMMAerror(hdr,"I & S Must Be Multiple Of 1/2",noret); break;// (11)
    case 12:GAMMAerror(hdr,"Interaction Requires I>1/2",   noret); break;// (12)
    case 14:GAMMAerror(hdr,"No Quantum # Exists, I=0,S=1", noret); break;// (14)
    case 15:GAMMAerror(hdr,"G Tensor Must Have Electron",  noret); break;// (15)
    case 16:GAMMAerror(hdr,"e-/e- or nuc./nuc. Hyperfine", noret); break;// (16)
    case 17:GAMMAerror(hdr,"e-/nucleus Dipolar Interact.", noret); break;// (17)
    case 18:GAMMAerror(hdr,"Cannot Alter Spin Tensor",     noret); break;// (18)
    case 19:GAMMAerror(hdr,"Quantum # < 1/2 Specified?",   noret); break;// (19)
    case 20:GAMMAerror(hdr,"Negative Quantum # Specified?",noret); break;// (20)
    case 21:GAMMAerror(hdr,"Quantum # Not Multiple Of 1/2",noret); break;// (21)
    case 22:GAMMAerror(hdr,"Quantum Number < 1 Specified?",noret); break;// (22)
    case 23:GAMMAerror(hdr,"Unreasonably Large Quantum #", noret); break;// (23)
    case 30:GAMMAerror(hdr,"Spin Hilbert Space MisMatch",  noret); break;// (30)
    case 31:GAMMAerror(hdr,"Composite Space Tensor Bad",   noret); break;// (31)
    case 49:GAMMAerror(hdr,"Bad Tensor Component Index",   noret); break;// (49)
    case 50:GAMMAerror(hdr,"Component Index Is m=[-2,2]",  noret); break;// (50)
    case 51:GAMMAerror(hdr,"Setting Quantum # To 1/2",     noret); break;// (51)
    case 52:GAMMAerror(hdr,"Setting Quantum # To 0",       noret); break;// (52)
    case 53:GAMMAerror(hdr,"Improper Spin Designation",    noret); break;// (53)
    case 60:GAMMAerror(hdr,"Problems Setting Quantum #",   noret); break;// (60)
    case 61:GAMMAerror(hdr,"Problems Setting I Quantum #", noret); break;// (61)
    case 62:GAMMAerror(hdr,"Problems Setting S Quantum #", noret); break;// (62)
    case 63:GAMMAerror(hdr,"Cannot Set From Parameters",   noret); break;// (63)
    case 64:GAMMAerror(hdr,"Spins Have The Same Index",    noret); break;// (64)
    case 65:GAMMAerror(hdr,"Negative Spin Index In Pair!", noret); break;// (65)
    case 66:GAMMAerror(hdr,"Trouble Setting Spin Tensor",  noret); break;// (66)
    case 70:GAMMAerror(hdr,"Spin Must Be An Electron",     noret); break;// (70)
    case 71:GAMMAerror(hdr,"Spin Cannot Be An Electron",   noret); break;// (71)
    case 79:GAMMAerror(hdr,"Cant Set Interaction Tensor",  noret); break;// (79)
    case 80:GAMMAerror(hdr,"Cannot Set Single Spin Tensor",noret); break;// (80)
    case 81:GAMMAerror(hdr,"Cannot Set Spin Pair Tensor",  noret); break;// (81)
    case 82:GAMMAerror(hdr,"Electron Paired With Nucleus", noret); break;// (82)
    case 83:GAMMAerror(hdr,"Disallowed Spin Pairing",      noret); break;// (83)
    case 84:GAMMAerror(hdr,"Electron Paired With Electron",noret); break;// (84)
    case 85:GAMMAerror(hdr,"Nucleus Paired With Nucleus",  noret); break;// (85)
    }
  }

volatile void IntRank2T::ISTfatal(int eidx) const
  {
  ISTerror(eidx, 1);				// Output error message
  if(eidx) ISTerror(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }


/* Err Ind.     Default Message            Err Ind.     Default Message
   -------- -----------------------------  -------- ---------------------------
      1     Problems With File pname          2     Cannot Read Parameter pname
      3     Invalid Use of Function pname
      4     Use Of Deprecated Function pname
      5     Please Use Class Member Funciton pname
   default  Unknown Error - pname                                           */

void IntRank2T::ISTerror(int eidx, const string& pname, int noret) const 
  {
  string hdr("Rank 2 Spin Tensor");
  string msg;
  switch(eidx)
    {
    case 8:							// (8)
      msg = string("Unknown Single Spin Interaction Type Specified");
      GAMMAerror(hdr, msg, noret); break;
    case 9:							// (9)
      msg = string("Problems During ")
          + pname + string(" Interaction Spin Tensor Construction");
      GAMMAerror(hdr, msg, noret); break;
    case 10:                                                    // (10)
      msg = string("Unknown Single Spin Interaction Type, ")
          + pname + string(" Specified");
      GAMMAerror(hdr, msg, noret); break;
    case 11:							// (11)
      msg = string("Unknown Spin Pair Interaction Type, ")
          + pname + string(" Specified");
      GAMMAerror(hdr, msg, noret); break;
    case 12:							// (12)
      msg = string("Odd Spin Quantum Value Of ")
          + pname + string(" Specified?");
      GAMMAerror(hdr, msg, noret); break;
    case 13:							// (13)
      msg = string("Problems Adjusting ")
          + pname + string(" Interaction Spin Tensor");
      GAMMAerror(hdr, msg, noret); break;
    case 14:							// (14)
      msg = string("Perhaps You Should Use A ")
          + pname + string(" Interaction Instead");
      GAMMAerror(hdr, msg, noret); break;
    case 15:							// (15)
      msg = string("Improper Spin Quantum Numbers ") + pname;
      GAMMAerror(hdr, msg, noret); break;
    case 65:							// (65)
      msg = string("Spin Pair Indices Are ") + pname;
      GAMMAerror(hdr, msg, noret); break;
    case 79:							// (79)
      msg = string("Spin Tensor Interaction Index Is ") + pname;
      GAMMAerror(hdr, msg, noret); break;
    case 80:							// (80)
      msg = string("Spin Tensor Spin Index Is ") + pname;
      GAMMAerror(hdr, msg, noret); break;
    case 81:							// (81)
      msg = string("Spin Tensor Spin Indices Are ") + pname;
      GAMMAerror(hdr, msg, noret); break;
    case 108:                                                   // (108)
     msg = string("Specified Spin Paring Of ")  + pname;
     GAMMAerror(hdr, msg, noret); break;
    case 109:                                                   // (109)
     msg = string("Specified Spin Of Type ")  + pname;
     GAMMAerror(hdr, msg, noret); break;
    case 110:                                                   // (110)
     msg = string("Use Of Unknown Isotope Type ")  + pname;
     GAMMAerror(hdr, msg, noret); break;
    default: GAMMAerror(hdr, eidx, pname, noret); break;
    }
  }

volatile void IntRank2T::ISTfatal(int eidx, const string& pname) const
  {
  ISTerror(eidx, pname, eidx);			// Output error message
  if(eidx) ISTerror(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii                INTERACTION SPIN TENSOR SETUP FUNCTIONS
// ____________________________________________________________________________

/* The functions in this section are PRIVATE because the value(s) of I (and S)
   MUST be set prior to their use. Thus we protect the user from trying to 
   generate the tensor components before a spin Hilbert space has been set.

   As a further note, since all non-null construction of rank 2 interaction
   spin tensors use these functions, any linked lists of interactions will
   depend on this section and this section only.                             */

// ----------------------------------------------------------------------------
//        Set Up Single Spin - Field Interaction Spin Tensor Components
// ----------------------------------------------------------------------------

/* We set the 5 spherical spin tensor components applicable to interactions
   involving a single spin interacting with a static magnetic field. Two such
   cases would be shift anisotropy (NMR) & g electron (EPR) interactions. Since
   this class maintains "scaled" spin tensors, all single spin - field  
   interactions will share the same spin tensors. The five arrays will only
   differ between single spin interactions if the spin quantum value differs
   (or if the field orientation is set differently). Typically, the static
   magnetic field is taken to be a normalized vector pointing along +z.
   
                                            +
                                         m  |
                              T    = (-1)  T
                               2,m          2,-m

                     1/2
                  [4]                      1
            T   = |-| * I         T    = - - I           T    = 0
             2,0  [6]    z         2,1     2  +           2,2 

		                  SINGLE SPIN - FIELD INTERACTIONS
           Input	IST	: Interaction Spin Tensor (this)
           Output       none    : Interaction Spin Tensor spherical
                                  spin components are generated
           Note                 : MAKE SURE Ival & Sval ARE SET FIRST!
	   Note			: These are for the spin with field
				  ala shift ansiotropy and electron G
				  interactions. Field along +z, norm 1
	   Note			: For G interactions Ival is usually 2
	  			  (Ival = 2*0.5 + 1, where I = 1/2 e-)       */

void IntRank2T::setSPF()
  {
  matrix IM = Im(Ival);                         // The operator I-
  matrix IP = Ip(Ival);                         // The operator I+
  matrix IZ = Iz(Ival);                         // The operator Iz
  T0  = (2./sqrt(6.))*IZ;			// T20  = 2*Iz/sqrt(6)
  T1  = -0.5*IP;				// T21  =-1/2(Bz * I+)
  Tm1 =  0.5*IM;				// T2m1 = 1/2(Bz * I-)
  Tm2 = 0.0*Ie(Ival);				// T2   = 0
  T2  = Tm2;					// T2m2 = 0
  }

// ----------------------------------------------------------------------------
//            Set Up Single Spin Interaction Spin Tensor Components
// ----------------------------------------------------------------------------

/* Here we set the 5 spherical spin tensor components applicable to single spin
   interactions. This currently includes only quadrupolar interactions. However
   since this class maintains "scaled" spin tensors, all single spin 
   interactions would share the same spin tensors. The five arrays will only
   differ between single spin interactions if the spin quantum value differs.                         
                                              +
                                           m  |
                                T    = (-1)  T
                                 2,m          2,-m
            1/2
         [1]      2                      1                          1  2
   T   = |-| * [3I - I(I+1)]     T   = - -(I I + I I )       T    = - I
    2,0  [6]      z               2,1    2  + z   z +         2,2   2  + 

	  				    SINGLE SPIN-SELF INTERACTION 
           Input                IST	: Interaction Spin Tensor (this)
           Output               none    : Interaction Spin Tensor spherical
                                          spin components are generated
           Note				: MAKE SURE Ival & Sval ARE SET FIRST!
	   Note				: These are for the spin with itself
					  ala quadrupolar interactions       */

void IntRank2T::setSPQ() 
  {
  double I = double(Ival-1)/2.0;		// Get spin quantum number
  matrix IE = Ie(Ival);                         // The operator 1
  matrix IM = Im(Ival);                         // The operator I-
  matrix IP = Ip(Ival);                         // The operator I+
  matrix IZ = Iz(Ival);                         //            2  
  T0  = (3.*IZ*IZ-(I*(I+1))*IE)/sqrt(6.);	// T20  = [3Iz -I(I+1)]/sqrt(6)
  T1  = -0.5*(IP*IZ + IZ*IP);			// T2m1 = -(I+Iz + IzI+)/2
  Tm1 = 0.5*(IM*IZ + IZ*IM);			// T2m1 =  (I-Iz + IzI-)/2
  T2  = 0.5*IP*IP;				// T22  =  (I+I+)/2
  Tm2 = 0.5*IM*IM;				// T2m2 =  (I-I-)/2
  }


// ----------------------------------------------------------------------------
//            Set Up Spin Pair Interaction Spin Tensor Components
// ----------------------------------------------------------------------------

/* Here we set the 5 spherical spin tensor components applicable to spin-spin
   interactions. These currently include hyperfine (elecron-nucleon) and
   dipolar (electron-electron and nucleon-nucleon) interactions.  Since this
   class maintains "scaled" spin tensors, all spin-spin interactions share the
   same spin tensors. The five arrays will only differ between interactions if
   the spin quantum values differ.                         
                                              +
                                           m  |
                                T    = (-1)  T
                                 2,m          2,-m
           1/2
        [1]                              1                            1   
  T   = |-| * [3I S - I.S]       T   = - -(I S + I S )         T    = - I S
   2,0  [6]      z z              2,1    2  + z   z +           2,2   2  + +  

	  				      SPIN PAIR INTERACTIONS
           Input                IST	: Interaction Spin Tensor (this)
           Output               none    : Interaction Spin Tensor spherical
                                          spin components are generated
           Note				: MAKE SURE Ival & Sval ARE SET FIRST!
	   Note				: Hyperfine interactions will typically
	  				  have either Ival and/or Sval of 2
	  				  (Ival = 2*0.5 + 1, where I = 1/2 e-)
	   Note				: These do NOT 't add anything to
	  				  the linked list Dip_HF_list        */

// sosi - In setSPSP the zero should use IP & IM since we have to do both now!
 
void IntRank2T::setSPSP()
  {
  matrix IE = Ie(Ival);                         // The operator Ie
  matrix SE = Ie(Sval);                         // The operator Se
  matrix IM = tensor_product(Im(Ival), SE);     // The operator I-
  matrix IP = tensor_product(Ip(Ival), SE);     // The operator I+
  matrix IZ = tensor_product(Iz(Ival), SE);     // The operator Iz
  matrix IY = tensor_product(Iy(Ival), SE);     // The operator Iy
  matrix IX = tensor_product(Ix(Ival), SE);     // The operator Ix
  matrix SM = tensor_product(IE, Im(Sval));     // The operator S-
  matrix SP = tensor_product(IE, Ip(Sval));     // The operator S+
  matrix SZ = tensor_product(IE, Iz(Sval));     // The operator Sz
  matrix SY = tensor_product(IE, Iy(Sval));     // The operator Sy
  matrix SX = tensor_product(IE, Ix(Sval));     // The operator Sx
  T0  =  (2.*IZ*SZ - IX*SX - IY*SY)/sqrt(6.);	// T20  = (3IzSz-I.S)/sqrt(6)
  T1  = -0.5*(IP*SZ + IZ*SP);               	// T21  =-1/2(I+Sz + IzS+)
  Tm1 =  0.5*(IM*SZ + IZ*SM);               	// T2m1 = 1/2(I-Sz + IzS-)
  T2  =  0.5*IP*SP;                         	// T22  = 1/2*I+*S+
  Tm2 =  0.5*IM*SM;                         	// T2m2 = 1/2*I-*S-
  }


// ____________________________________________________________________________
// iii           RANK 2 SPIN TENSOR FROM PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

/* These functions allow spin tensors to be set up from external ASCII files
   and/or GAMMA parameters sets. Since the parameters which define a rank 2
   spin tensor are often inter-related, we keep these private to insure that
   tensor integrity is maintained. Since objects of type IntRank2T are not
   intended to be used directly, these functions are typically used by classes
   derived from us.
  
           Input                IST	: Rank 2 spin tensor (this)
                                pset    : A parameter set
                                Pbase   : Parameter base name
                                idx     : Index value
                                warn    : Warning output label
                                           0 = no warnings
                                          !0 = warnings
           Output               IsoI	: Tensor spin isotope type
                                          obtained from parameters in pset
           			I	: Tensor spin quantum number
                                          obtained from parameters in pset
	   Note				: No isotope restrictions are herein
	  				  considered

   Look for I value. Allowed parameters are the following:

    Iso(idx)         - Spin I type (2H, 131Xe, ...)                                
    Ion(idx)	     - Magnetic Ion type (Ce3+, Ho3+,....)
    Pbase,  Pbase(#) - Spin I quantum number (0.5, 1.5, 2.0, 2.5,..)         */

// ----------------------------------------------------------------------------
//                          Single Spin Spin Tensors
// ----------------------------------------------------------------------------

bool IntRank2T::getIso(const ParameterSet& pset, std::string& IS,
                                                        int i, bool warn) const
  {
  ParameterSet::const_iterator item;         // A pix into parameter list
  string pname("Iso");				// Parameter name base Iso
  string Nidx = "";                             // Name addition per index
  if(i >= 0) Nidx				// If index exists then set up
    += string("(")+Gdec(i)+string(")");		// the param. name to append
  pname += Nidx;				// Parameter name Iso(i)
  string pstate;				// Dummy string for comment
  item = pset.seek(pname);			// Seek parameter in pset
  IS = "";					// Empty any isotope label
  if(item != pset.end())			// If Iso(i) found, parse
    { 						//   info, set I value from 
    (*item).parse(pname,IS,pstate);		//   type, return we found it
    return true;
    }
  pname = string("Ion") + Nidx;			// This is parameter name
  item = pset.seek(pname);			// Seek parameter in pset
  if(item != pset.end())			// If Ion(i) found, parse
    { 						//   info, set I value from
    (*item).parse(pname,IS,pstate);		//   type, return we found it
    CubicIon X;					//   An empty cubic ion
    X.initialize();				//   Insure cubic ions set up
    return true;
    }
  if(warn)					// If we just cannot find it
    {						// issue warnings as needed
    ISTerror(2, pname, 1);			//   Cant find Ion(i)
    pname = string("Iso") + Nidx;		//   Alternate parameter name
    ISTerror(2, pname, 1);			//   Cant find Iso(i)
    }						//   These are never fatal
  return false;					// Return we failed
  }

bool IntRank2T::scanIso(const ParameterSet& pset, std::string& IS,
                                                        int i, bool warn) const
  {
  ParameterSet::const_iterator item;         // A pix into parameter list
  string pname("Iso");				// Parameter name base Iso
  string Nidx = "";                             // Name addition per index
  if(i >= 0) Nidx				// If index exists then set up
    += string("(")+Gdec(i)+string(")");		// the param. name to append
  pname += Nidx;				// Parameter name Iso(i)
  string pstate;				// Dummy string for comment
  cout << "\n\t\t\tLooking For Parameter " 
       << pname;
  item = pset.seek(pname);			// Seek parameter in pset
  IS = "";					// Empty any isotope label
  if(item != pset.end())			// If Iso(i) found, parse
    { 						//   info, set I value from 
    (*item).parse(pname,IS,pstate);		//   type, return we found it
    cout << " .... OK";
    return true;
    }
  cout << " .... Not Found";
  pname = string("Ion") + Nidx;			// This is parameter name
  cout << "\n\t\t\tLooking For Parameter " 
       << pname;
  if(item != pset.end())			// If Ion(i) found, parse
    { 						//   info, set I value from
    cout << " .... OK";
    (*item).parse(pname,IS,pstate);		//   type, return we found it
    return true;
    }
  cout << " .... Not Found";
  if(warn)					// If we just cannot find it
    {						// issue warnings as needed
    ISTerror(2, pname, 1);			//   Cant find Ion(i)
    pname = string("Iso") + Nidx;		//   Alternate parameter name
    ISTerror(2, pname, 1);			//   Cant find Iso(i)
    }						//   These are never fatal
  return false;					// Return we failed
  }

bool IntRank2T::getIqn(const ParameterSet& pset, const string& Pbase,
                                          double& qn, int idx, bool warn) const
  {
  qn = 0;					// Start with no quantum #
  ParameterSet::const_iterator item;         // A pix into parameter list
  string Nidx = "";                             // Name addition per index
  if(idx >= 0) Nidx                             // If index exists then set up
    += string("(")+Gdec(idx)+string(")");	// the parameter name to append
  string pstate;				// string for parameter comment
  string pname = Pbase +  Nidx;			// Parameter name for I(#)
  item = pset.seek(pname);			// Seek parameter in pset
  if(item != pset.end())			// If Pname(idx) was found,
    {						// parse info & set qn value
    (*item).parse(pname,qn,pstate);		// then return we have found it
    return true;
    } 						//   info, set I value from 
  if(warn)					// If we just cannot find it
    {						// issue warnings as needed
    ISTerror(2, pname, 1);			//   Cant find Pbase(idx)
    ISTerror(52, 1); 				//   Problems set quantum #
    }						//   These are never fatal
  return false;					// Return we failed
  }  

bool IntRank2T::scanIqn(const ParameterSet& pset, const string& Pbase,
                                          double& qn, int idx, bool warn) const
  {
  qn = 0;					// Start with no quantum #
  ParameterSet::const_iterator item;         // A pix into parameter list
  string Nidx = "";                             // Name addition per index
  if(idx >= 0) Nidx                             // If index exists then set up
    += string("(")+Gdec(idx)+string(")");	// the parameter name to append
  string pstate;				// string for parameter comment
  string pname = Pbase +  Nidx;			// Parameter name for I(#)
  cout << "\n\t\t\tLooking For Parameter " 
       << pname;
  item = pset.seek(pname);			// Seek parameter in pset
  if(item != pset.end())			// If Pname(idx) was found,
    {						// parse info & set qn value
    cout << " .... OK";
    (*item).parse(pname,qn,pstate);		// then return we have found it
    return true;
    } 						//   info, set I value from 
  cout << " .... Not Found";
  if(warn)					// If we just cannot find it
    {						// issue warnings as needed
    ISTerror(2, pname, 1);			//   Cant find Pbase(idx)
    ISTerror(52, 1); 				//   Problems set quantum #
    }						//   These are never fatal
  return false;					// Return we failed
  }  

bool IntRank2T::SpinCheck(const string& II, bool warn) const
  {
  if(Isotope::known(II)) 			// Make sure that string II
    return true; 				// is a valid isotope in GAMMA
  if(warn)					// Not a valid isotope? Then
    {						// we issue warnings if set
    ISTerror(110, II, 1);			//   Unknown isotope type
    ISTerror(53,1); 				//   Bad spin type specified
    }						//   This is never fatal
  return false;					// Failure, don't know type
  }

bool IntRank2T::SpinCheck(const Isotope& II, bool epr, bool warn) const
  {
  if(!epr &&  !II.electron()) return true;	// For nmr, nuclei OK
  if( epr &&   II.electron()) return true;	// For epr, electrons OK
  if(warn)					// Isotope wrong type,warn
    epr?ISTerror(70,1):ISTerror(71,1);		// Must/Must Not Be Electron
  return false;					// Failure, wrong spin type
  }

bool IntRank2T::SpinCheck(double Iqn, bool quad, bool warn) const
  {
  int twoI = int(2.0*Iqn);              // 2*I: 1/2->1, 1->2, 3/2->3, ......
  if(twoI<1)				// Insure minimum I value is 1/2
    {					// If less than 1/2 we issue 
    if(warn)				// warnings if needed
      {
      string oddI = Gform("%5.1f",Iqn);	//   String for I value
      ISTerror(19, 1);			//   Quantum # < 1/2 specified
      ISTerror(12, oddI, 1);		//   Odd spin quantum value
      }
    return false;			// We failed because I<1/2
    }	
  if(twoI > 40)				// Insure maximum I values is 20
    {					// (which is way way big)
    if(warn)				// If I more than 20 we issue
      {					// warnings if needed
      string oddI = Gform("%5.1f",Iqn);	//   String for I value
      ISTerror(23, 1);			//   Unreasonably large #
      ISTerror(12, oddI, 1);		//   Odd spin quantum value
      }
    return false;			// We failed because I>20
    }	
  if(2.0*Iqn - double(twoI))		// Insure I is an integer 
    { 					// multiple of 1/2
    if(warn)				// If not multiple of 1/2 we
      {					// issue warnings if needed
      string oddI = Gform("%5.1f",Iqn);	//   String for I value
      ISTerror(21, 1);			//   Not multiple of 1/2
      ISTerror(12, oddI, 1);		//   Odd spin quantum value
      }
    return false;			// We failed because I != m*1/2
    }	
  if(quad && Iqn <= 0.5)		// If quadrupolar need I > 1/2
    {
    if(warn)				// If not > 1/2 we issue 
      {					// warnings as needed
      string oddI = Gform("%5.1f",Iqn);	//   String for I value
      ISTerror(12, 1);			//   Not > 1/2
      ISTerror(12, oddI, 1);		//   Odd spin quantum value
      }
    return false;
    }
  return true;
  }

// ----------------------------------------------------------------------------
//                          Spin Pair Spin Tensors
// ----------------------------------------------------------------------------

bool IntRank2T::getIsos(const ParameterSet& pset, int idxI, int idxS, 
                                       string& II, string& IS, bool warn) const
  {
  bool TFI=getIso(pset,II,idxI,warn);		// Try & read spin I type
  bool TFS=getIso(pset,IS,idxS,warn);		// Try & read spin S type
  return TFI&TFS;
  }

bool IntRank2T::scanIsos(const ParameterSet& pset, int idxI, int idxS, 
                                       string& II, string& IS, bool warn) const
  {
  bool TFI=scanIso(pset,II,idxI,warn);		// Try & read spin I type
  bool TFS=scanIso(pset,IS,idxS,warn);		// Try & read spin S type
  return TFI&TFS;
  }

bool IntRank2T::getIqns(const ParameterSet& pset, const string& Pbase,
                               double& Iqn, double& Sqn, int i, bool warn) const
  {
  Iqn = 0.0;					// No I spin quantum number
  Sqn = 0.0;					// No S spin quantum number
  string p = Pbase + string("Iqn");		// Spin I qn parameter 
  bool TFI=getIqn(pset,p,Iqn,i,warn);		// Try & read spin I qn
  p = Pbase + string("Sqn");			// Spin S qn parameter 
  bool TFS=getIqn(pset,p,Sqn,i,warn);		// Try & read spin S qn
  return TFI&TFS;
  }

bool IntRank2T::getIqns(const ParameterSet& pset, const string& Pbase,
                  double& Iqn, double& Sqn, int idxI, int idxS, bool warn) const
  {
  Iqn = 0.0;					// No I spin quantum number
  Sqn = 0.0;					// No S spin quantum number
  string p = Pbase + string("qn");		// Spin I qn parameter 
  bool TFI=getIqn(pset,p,Iqn,idxI,warn);	// Try & read spin I qn
  bool TFS=getIqn(pset,p,Sqn,idxS,warn);	// Try & read spin S qn
  return TFI&TFS;
  }

bool IntRank2T::scanIqns(const ParameterSet& pset, const string& Pbase,
                               double& Iqn, double& Sqn, int i, bool warn) const
  {
  string p = Pbase + string("Iqn");		// Spin I qn parameter 
  cout << "\n\t\t\tLooking For Parameter " << p;
  bool TFI=getIqn(pset,p,Iqn,i,warn);		// Try & read spin I type
  if(TFI) cout << " .... OK";
  else    cout << " .... Not Found";
  p = Pbase + string("Sqn");			// Spin S qn parameter 
  cout << "\n\t\t\tLooking For Parameter " << p;
  bool TFS=getIqn(pset,p,Sqn,i,warn);		// Try & read spin I type
  if(TFS) cout << " .... OK";
  else    cout << " .... Not Found";
  return TFI&TFS;
  }

bool IntRank2T::SpinCheck(const string& II, const string& IS, bool warn) const
  {
  if(!SpinCheck(II, warn)) return false; 	// Insure II is valid isotope
  if(!SpinCheck(IS, warn)) return false; 	// Insure IS is valid isotope
  return true;
  }

bool IntRank2T::SpinCheck(const Isotope& II, const Isotope& IS, 
                                                    bool epr, bool warn) const
  {
  if( epr &&  II.nepair(IS)) return true;	// For epr e-/nucleus pair OK
  if(!epr && !II.nepair(IS)) return true;	// For nmr e-/e- or n/n OK
  if(warn)
    {
    if(epr &&  II.electron()) ISTerror(84, 1);	//   EPR: electron w/ electron
    if(epr && !II.electron()) ISTerror(85, 1);	//   EPR: nucleus w/ nucleus
    if(!epr)                  ISTerror(82, 1);	//   NMR: electron w/ nucleus
    ISTerror(83, 1);				//   Disallowed spin pairing
    string pn = II.symbol() + " With "		//   String for spins paired
                            + IS.symbol();
    ISTerror(108,pn,1); 			//   Specified spin pairing
    }
  return false;
  }

bool IntRank2T::SpinCheck(double Iqn, double Sqn, bool warn) const
  {
  if(SpinCheck(Iqn, false, warn)		// Make sure Iqn is valid 
  && SpinCheck(Sqn, false, warn)) return true;	// Make sure Sqn is valid
  if(warn)
    {
    string oddIS = Gform("%5.1f",Iqn)
                 + string(" & ")
                 + Gform("%5.1f",Sqn);
    ISTerror(15, oddIS, 1);
    }
  return false;
  }


bool IntRank2T::SpinCheck(int idxI, int idxS, bool warn) const
  {
  if(idxI<0 || idxS<0) 				// Insure neither spin pair
    {						// index is negative
    if(warn)					// If either is negative we
      {						// issue warnings and fail
      ISTerror(65, 1);				//   Negative index in pair
      string pname=Gdec(idxI)+" & "+Gdec(idxS);	//   String for indices
      ISTerror(65, pname, 1);			//   These spin pair indices
      }
    return false;				// We have failed
    }
  if(idxI != idxS) return true;			

  if(warn)
    {
    ISTerror(64, 1);				// Both spins have same index
    string pname=Gdec(idxI)+" & "+Gdec(idxS);	// String for indices
    ISTerror(65, 1);				// These spin pair indices
    }
  return false;
  }








bool IntRank2T::setTIso(const ParameterSet& pset,
                                           int ttype, int idx, int warn)
  {
  string II;
  bool TF=getIso(pset,II,idx,warn?true:false);	// Try & read spin type
  if(!TF)					// If can't read the spin type
    {						// issue warnings, fail (quit?)
    if(warn) (warn>1)?ISTfatal(63):ISTerror(63);//   Cant set from parameters
    return false;				//   We were unsuccessful
    }
  if(!SpinCheck(II))				//  First insure we know this
    { if(warn>1) ISTfatal(63); return false; }	//    Cant set from parameters
  Isotope IS = Isotope(II);			//  If we do, make an isotope
  switch(ttype)					//  and try & construct ourself
    {
    default:
    case 0:					// 0: Spin-Field No e- (CSA)
      if(!SpinCheck(IS,0,warn?1:0))		//   Cannot Be An Electron
        { 
        if(warn>1) ISTfatal(63);		//    Cant set from parameters
        else       ISTerror(63,1);		//    Specified spin type II
        return false;
        } 
      *this = IntRank2T(IS.qn());
      return true;
      break;
    case 1:					// 1: Spin-Field    e- (G)
      if(!SpinCheck(IS,1,warn?1:0))		//   Must Be An Electron
        { 
        if(warn>1) ISTfatal(63);		//    Cant set from parameters
        else       ISTerror(63,1);		//    Specified spin type II
        return false;
        } 
      *this = IntRank2T(IS.qn());
      return true;
      break;
    case 2:					// 2: Spin-Self  No e- (Quad)
      if(!SpinCheck(IS,0,warn?1:0))		//   Cannot Be An Electron
        { 
        if(warn>1) ISTfatal(63);		//    Cant set from parameters
        else       ISTerror(63,1);		//    Specified spin type II
        return false;
        } 
      double qn = 0;
      *this = IntRank2T(IS.qn(),false);
      if(qn<1.0)				//    Cannot have I<1
        {
        if(warn)   ISTerror(22,1);		//    Cannot have I<1
        if(warn>1) ISTfatal(63);		//    Cant set from parameters
        else       ISTerror(109,II,1);		//    Specified spin type II
        return false;
          }
      return true;
      break;
    }
  return false;
  }

bool IntRank2T::setTIsos(const ParameterSet& pset,
                                int ttype, int idxI, int idxS, int warn)
  {
  string II, IS;				// Strings for isotope names
  if(!getIsos(pset,idxI,idxS,II,IS,warn?true:false))	// Try tp get isotope names
    return false;				// Return false if failure

  Isotope ISI = Isotope(II);			//  Make an isotope for spin I
  Isotope ISS = Isotope(IS);			//  Make an isotope for spin S
  switch(ttype)
    {
    default:
    case 0:					// 0: Spin-Spin n/n or e-/e-
      if(ISI.nepair(ISS))
        { 
        if(warn)   ISTerror(71,1);		//    Cant Be Electron-Nucleus
        if(warn>1) ISTfatal(63);		//    Cant set from parameters
        else
          {
          ISTerror(109,II,1);			//    Specified spin type II
          ISTerror(109,IS,1);			//    Specified spin type IS
          }
        return false;
        } 
     break;
    case 1:					// 1: Spin-Spin n/e- pair
      if(!ISI.nepair(ISS))
        { 
        if(warn)   ISTerror(70,1);		//    Must Be Electron-Nucleus
        if(warn>1) ISTfatal(63);		//    Cant set from parameters
        else
          {
          ISTerror(109,II,1);			//    Specified spin type II
          ISTerror(109,IS,1);			//    Specified spin type IS
          }
       return false;
       } 
     break;
     }
  *this = IntRank2T(ISI.qn(), ISS.qn());
  return true;
  }


bool IntRank2T::setTqn(const ParameterSet& pset, const string& PBase,
                                           int ttype, int idx, int warn)
  {
  double QN;
  bool TF = getIqn(pset,PBase,QN,idx,warn?true:false);	// Try & read spin type
  if(!TF) return false;					// Quit if no qn found
  if(!SpinCheck(QN, warn?true:false)) return false;	// Quit if bad qn found
  switch(ttype)					//  If quantum # valid, try
    { 						//  & construct ourself
    default:
    case 0:					// 0: Spin-Field No e- (CSA)
    case 1:					// 1: Spin-Field    e- (G)
      *this = IntRank2T(QN);
      return true;
      break;
    case 2:					// 1: Spin-Field    e- (G)
      if(QN<1.0)
        {
        string oddI = Gform("%5.1f", QN);
        if(warn)   ISTerror(22,1);		// Cannot have I<1
        if(warn>1) ISTfatal(63);		// Cant set from parameters
        else       ISTerror(12, oddI, 1);	// Odd mz value found
        return false;
        }
      *this = IntRank2T(QN, false);
      return true;
      break;
    }
  return false;
  }

bool IntRank2T::setTqns(const ParameterSet& pset, const string& PBase,
                                                         int idx, int warn)
  {
  string pname = PBase + string("I");		// 1st Spin Parameter 
  double QNI;
  bool TF = getIqn(pset,PBase,QNI,idx,warn?true:false);	// Try & read I spin qn
  if(!TF) return false;			// Quit if no qn found
  if(!SpinCheck(QNI, warn?true:false)) return false;	// Quit if bad qn found
  pname = PBase + string("S");			// 2nd Spin Parameter 
  double QNS;
  TF = getIqn(pset,PBase,QNS,idx,warn?true:false);	// Try & read S spin qn
  if(!TF) return false;			// Quit if no qn found
  if(!SpinCheck(QNS, warn?true:false)) return false;	// Quit if bad qn found
  *this = IntRank2T(QNI, QNS);
  return true;
  }

bool IntRank2T::setTqns(const ParameterSet& pset, const string& PBase,
                                               int idxI, int idxS, int warn)
  {
  double QNI;
  bool TF = getIqn(pset,PBase,QNI,idxI,warn?true:false);// Try & read I spin qn
  if(!TF) return false;			// Quit if no qn found
  if(!SpinCheck(QNI, warn?true:false)) return false;	// Quit if bad qn found
  double QNS;
  TF = getIqn(pset,PBase,QNS,idxS,warn?true:false);// Try & read S spin qn
  if(!TF) return false;			// Quit if no qn found
  if(!SpinCheck(QNS, warn?true:false)) return false;	// Quit if bad qn found
  *this = IntRank2T(QNI, QNS);
  return true;
  }

bool IntRank2T::setT(const ParameterSet& pset, const string& PBase,
                                           int ttype, int idx, int warn)
  {
  if(setTIso(pset,        ttype, idx, warn?1:0)) return true;
  if(setTqn( pset, PBase, ttype, idx, warn?1:0)) return true;
  if(warn)
    {
    ISTerror(80,1);				// Cant set 1 spin T
    string pname = Gdec(idx);
    if(warn > 1) ISTfatal(80, pname);
    else         ISTerror(80, pname,1);
    }
  return false;
  }

bool IntRank2T::setTs(const ParameterSet& pset, const string& PBase,
                                                int ttype, int idxI, int warn)
  {
  if(setTqns( pset, PBase, idxI, warn?1:0)) return true;
  if(warn)
    {
    ISTerror(79,1);
    string pname = Gdec(idxI);
    if(warn > 1) ISTfatal(79, pname);
    else         ISTerror(79, pname,1);
    }
  return false;
  }

bool IntRank2T::setTs(const ParameterSet& pset, const string& PBase,
                                      int ttype, int idxI, int idxS, int warn)
  {
  if(setTIsos(pset,        ttype, idxI, idxS, warn?1:0)) return true;
  if(setTqns( pset, PBase,        idxI, idxS, warn?1:0)) return true;
  if(warn)
    {
    ISTerror(81,1);				// Cant set spin pair T
    string pname = Gdec(idxI)+" & "+Gdec(idxS);
    if(warn > 1) ISTfatal(81, pname);
    else         ISTerror(81, pname,1);
    }
  return false;
  }

// ____________________________________________________________________________
// iv           RANK 2 SPIN TENSOR HILBERT SPACE TRANSFORMATIONS
// ____________________________________________________________________________

/* These functions will take a spin tensor component (matrix) defined in the
   single spin or spin pair Hilbert space and convert it into an array defined
   in any specified composite Hilbert space. Individual Hilbert spaces in the
   specified composite must match those of the spins involved in the tensor.

   For a single spin interaction the conversion is quite simple. We just take
   direct products with identity matrices in the single spin Hilbert space, 
   except when the spin is ourself in which case we use our own component.
   For an N spin composite Hilbert space of dimension chs, we have

  Mx(chs) = Ie(hs )*Ie(hs )* ... *Ie(hs   )*Mx(hs )*Ie(hs   )* .... *Ie(hs   )
                 0       1             i-1       i       i+1              N-1

   For a spin pair interaction their is no simple product that will do the job,
   so the conversion is a little more messy. We have to map the composite space
   matrix elements <I|Mx(chs)|J> back to the spin pair elements <i|Mx(hs)|j>,
   or vice-versa depending on ones point of view. We can do this by first 
   mapping the basis functions of the spin pair                            */

//matrix HSorder(const vector<int>& HSs, int i, int j)

vector<int> IntRank2T::HSorder(const vector<int>& HSs, int i)
  {
  int ns = HSs.size();				// Number of spins
  int j;					// Spin index
  int cHS = 1;					// Composite Hilbert space
  int HSo = 1;					// Hilbert space < spin i
  int HSi = HSs[i];				// Hilbert space   spin i
  int HSf = 1;					// Hilbert space > spin i
  for(j=0; j<i; j++)    HSo *= HSs[j];		// Set space before spin i
  for(j=i+1; j<ns; j++) HSf *= HSs[j];		// Set space after spin i 
  cHS = HSo * HSi * HSf;			// Set full composiste space
  if(HSo == 1) HSo = 0;				// Cannot have HS of 1
  if(HSf == 1) HSf = 0;				// Cannot have HS of 1
  vector<int> ssmap(cHS);
  int a, b, c;					// Subspace basis indices
  int B;					// Composite space index
  for(a=0, B=0; a<HSo || !a; a++)		// Loop space < spin i 
    for(b=0; b<HSi; b++)			// Loop space   spin i
      for(c=0; c<HSf || !c; c++, B++)		// Loop space > spin i
        ssmap[B] = b;				// Set basis mapping b->B
  return ssmap;
  }

matrix IntRank2T::blow_up(const matrix& mx, const vector<int>& HSs,
                                                            int i, int j) const
  {
  if(HSs[i] != Ival)				// First check that spin i's
    {						// Hilbert space matches us
    ISTerror(30, 1);				//   Hilbert space mismatch
    ISTfatal(31);				//   Bad composite space tensor
    }
  matrix bmx;					// The composite space array
  int ns = int(HSs.size());			// # spins in composite space
  int k;					// Temp. spin index
  if(Sval == 0)					// Single Spin Tensor
    {
    for(k=0; k<ns; k++)				//   Loop over spins for the
      {						//   composite space, taking
      if(i == k)				//   direct products with 
        bmx = tensor_product(bmx, mx);		//   ourself and I's
      else
        bmx = tensor_product(bmx, Ie(HSs[k])); 
      }
    }
  else						// Spin Pair Tensor
    {
    if(HSs[j] != Sval)				//   Check that spin j's
      {						//   Hilbert space matches us
      ISTerror(30, 1);				//     Hilbert space mismatch
      ISTfatal(31);				//     Bad comp. space tensor
      }
    if(ns == 2) return mx;			//   If only spins, we are it
//    matrix_type mt = mx.stored_type();		//   Original array storage
//    need to set matrix type too
    int chs=1;					//   Need composite space 
    for(k=0; k<ns; k++) chs *= HSs[k]; 		//   dimension
    bmx=matrix(chs,chs,complex0);		//   Set 0 matrix in comp. sp.

//        Map Basis Functions For Spins i & j As Well As Both Together

    vector<int> imap = HSorder(HSs, i);		//   Get mapping of i to chs
    vector<int> jmap = HSorder(HSs, j);		//   Get mapping of j to chs
    int A, B;					//   Comp. space basis indices
    int a, b;					//   Spin space basis indices
    vector<int> cmap(chs);			//   For composite space map
    for(A=0; A<chs; A++)			//   Loop comp. space basis
      cmap[A] = imap[A]*HSs[j] + jmap[A];

//          Map Basis Functions As Not Involved With Spins i & j

    vector<int> xmap(chs, 0);			//   For composite space map
    int sf = 1;
    for(k=0; k<ns; k++)
      {
      if(k!=i && k!=j)
        {
        imap = HSorder(HSs,k);
        for(A=0; A<chs; A++)
          xmap[A] += sf*imap[A];
        sf *= 100;
        }
      }

//       Now Map Matrix Elements From Two Spin Space To Composite Space
//      (Basis Functions Without i & j Must Be Identical Or No Element)

    for(A=0; A<chs; A++)			// Loop composite space basis
      {						// functions (for array rows)
      a = cmap[A];
      for(B=0; B<chs; B++)
        {
        b = cmap[B];
        if(xmap[A] == xmap[B]) 
          bmx.put(mx.get(a,b), A, B);
        }
      }
    }
  return bmx;					// Return composite space mx
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A          INTERACTION SPIN TENSOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

IntRank2T::IntRank2T()
  { Ival=0; Sval=0; }			// No I or S spin Hilbert space 

IntRank2T::IntRank2T(const IntRank2T &IST1)
  {
  Ival = IST1.Ival;			// Copy the I Hilbert space
  Sval = IST1.Sval;			// Copy the S Hilbert space
  T0   = IST1.T0;			// Spherical spin component m=0
  T1   = IST1.T1;			// Spherical spin component m=1
  Tm1  = IST1.Tm1;			// Spherical spin component m=-1
  T2   = IST1.T2;			// Spherical spin component m=2
  Tm2  = IST1.Tm2;			// Spherical spin component m=-2
  }

// ----------------------------------------------------------------------------
//     Direct Single Spin Constructors That Need Multiple Input Parameters
// ----------------------------------------------------------------------------

/* These constructors handle spin tensors applicable in single spin and
   spin-field interactions. An example of a single spin interaction is the
   nuclear quadrupolar interaction.  Examples of spin-field interactions are
   nuclear shift anistropy and electron g interactions. For these constructors
   we need to know about the spin Hilbert space of only on spin particle. To
   obtain the Hilbert space one either uses a spin isotope type declaration
   which implicitly knows the spin or one directly specifies the value.

           Input                IST	: Interaction Spin Tensor (this)
	  			IsoS	: Spin S isotope type
	                   or
	  			Iqn	: Spin I quantum number
                                ttype   : Spin tensor type
				 	   true = Spin Field (SA, G)
                                           false = Spin Self (Quad)
        			warn    : Warning output label
                                           0 = no warnings
                                           1 = non-fatal warnings
                                          >1 = warnings & exit
	   Output		none    : Spin Tensor constructed            */

IntRank2T::IntRank2T(const string& IsoI, bool ttype, int warn)
  {
  Isotope II(IsoI);				// Form isotope of type IsoI
  *this = IntRank2T(II, ttype , warn);		// Use overloaded constructor
  }						// Here qn = 0.5, 1.0, 1.5, ..

IntRank2T::IntRank2T(const Isotope& II, bool ttype, int warn)
  {
  double IIqn = II.qn();			// Get spin quantum number
  *this = IntRank2T(IIqn, ttype, warn);		// Use overloaded constructor
  }						// Here qn = 0.5, 1.0, 1.5, ..

IntRank2T::IntRank2T(double Iqn, bool ttype, int warn)
  {
  Ival = int(2.0*Iqn + 1);			// Set the I Hilbert space
  Sval = 0;					// Set the S Hilbert space
  if(!SpinCheck(Iqn,warn?true:false))		// Insure spin quantum number
    {						// is even multiple of 1/2
    if(ttype) ISTerror(9, "Single Spin",  1);	// Construction problems
    else      ISTerror(9, "Single Field", 1);	// Construction problems
    ISTfatal(0);				// Aborting....
    }
  if(ttype) setSPF();  				// Spin With Field Interaction
  else      setSPQ(); 				// Spin With Itself Interaction
  }

// ----------------------------------------------------------------------------
//     Direct Spin Pair Constructors That Need Multiple Input Parameters
// ----------------------------------------------------------------------------

/* These constructors handle spin tensors applicable in spin-spin interactions.
   One such interaction is dipole-dipole (nulcear-nuclear or electron-electron)
   and another would be hyperfine (nuclear-electron). These thus need to know
   about the spin Hilbert spaces of two spin particles. To obtain the Hilbert
   space one the two spins we either use spin isotope type declarations which
   implicitly knows the spins or one directly specifies their values.

           Input                IST	: Interaction Spin Tensor (this)
	  			IsoI	: Spin I isotope type
	  			IsoS	: Spin S isotope type
                            or
	  			Iqn	: Spin I quantum number
	  			Sqn	: Spin S quantum number
	   Output		none    : Spin Tensor constructed            */


IntRank2T::IntRank2T(const string& IsoI, const string& IsoS)
  {
  Isotope II(IsoI);				// Form isotope of type IsoI
  Isotope IS(IsoS);				// Form isotope of type IsoS
  *this = IntRank2T(II,IS);			// Use overloaded constructor
  }

IntRank2T::IntRank2T(const Isotope& II, const Isotope& IS)
  { *this=IntRank2T(II.qn(),IS.qn()); } 	// Use overloaded constructor

IntRank2T::IntRank2T(double Iqn, double Sqn)
  {
  Ival = int(2.0*Iqn + 1);			// Set the I Hilbert space
  Sval = int(2.0*Sqn + 1);			// Set the S Hilbert space
  if(!SpinCheck(Iqn, Sqn, 1))			// Insure reasonable Iqn & Sqn
    { 						// values, both >0 & n*1/2
    ISTerror(3, 1);
    ISTfatal(11);
    }
  setSPSP();
  }

// ----------------------------------------------------------------------------
//            Rank 2 Interaction Spin Tensor Assignment & Desctruction
// ----------------------------------------------------------------------------

IntRank2T::~IntRank2T () { }

void IntRank2T::operator= (const IntRank2T &IST1) 
  {
  Ival = IST1.Ival;			// Copy spin I Hilbert space
  Sval = IST1.Sval;			// Copy spin S Hilbert space
  T0   = IST1.T0;			// Spherical spin component m=0
  T1   = IST1.T1;			// Spherical spin component m=1
  Tm1  = IST1.Tm1;			// Spherical spin component m=-1
  T2   = IST1.T2;			// Spherical spin component m=2
  Tm2  = IST1.Tm2;			// Spherical spin component m=-2
  }


// ____________________________________________________________________________
// B       INTERACTION SPIN TENSOR SPIN TENSOR I & S VALUE ACCESS
// ____________________________________________________________________________

        // Input                IST     : Interaction Spin Tensor
        // Output               I,S     : Quantum value of I/S (e.g. 0.5, 1.5)
        //                      IV,SV   : Hilbert space of I/S (e.g. 2, 6) 
        //                      HS      : Spin Hilbert space of IST
 
double IntRank2T::Izval() const     { return double(Ival-1)/2.0; }
double IntRank2T::Szval() const     { return Sval?double(Sval-1)/2.0:0.0; }
int    IntRank2T::IV()    const     { return Ival; }
int    IntRank2T::SV()    const     { return Sval; }
int    IntRank2T::HS()    const     { return Sval?Ival*Sval:Ival; }
double IntRank2T::qn(bool TF) const { return TF?Izval():Szval(); }

// ____________________________________________________________________________
// C       INTERACTION SPIN TENSOR SPIN TENSOR COMPONENT ACCESS
// ____________________________________________________________________________

/* These functions provide direct access to individual spin tensor spherical
   components.  Note that if the spin tensor has not been initialized each
   component will be a NULL matrix.  The index spans m = [-2, 2]. Note that the
   functions will return arrays in the single or spin pair Hilbert space 
   unless a composite Hilbert space is provided.                             */

// ----------------------------------------------------------------------------
//        Single Spin Or Spin Pair Hilbert Space Spin Tensor Components
// ----------------------------------------------------------------------------
 
matrix IntRank2T::T20()           const { return T0;  }
matrix IntRank2T::T21()           const { return T1;  }
matrix IntRank2T::T2m1()          const { return Tm1; }
matrix IntRank2T::T22()           const { return T2;  }
matrix IntRank2T::T2m2()          const { return Tm2; }
matrix IntRank2T::T2m(int m)      const { return Tcomp(m); }  
matrix IntRank2T::Tcomp(int comp) const
  {
  switch(comp)
    {
    case 0:  return T0;  break;
    case 1:  return T1;  break;
    case -1: return Tm1; break;
    case 2:  return T2;  break;
    case -2: return Tm2; break;
    default:
      {
      ISTerror(49,1);
      ISTerror(50,1);
      ISTfatal(0);
      }
    }
  return T0;
  }

// ----------------------------------------------------------------------------
//               Composite Hilbert Space Spin Tensor Components
// ----------------------------------------------------------------------------

matrix IntRank2T::T20(const  vector<int>& HSs,int i,int j) const
  { return blow_up(T0, HSs, i, j); }
matrix IntRank2T::T21(const  vector<int>& HSs,int i,int j) const
  { return blow_up(T1, HSs, i, j); }
matrix IntRank2T::T2m1(const vector<int>& HSs,int i,int j) const
  { return blow_up(Tm1, HSs, i, j); }
matrix IntRank2T::T22(const  vector<int>& HSs,int i,int j) const
  { return blow_up(T2, HSs, i, j); }
matrix IntRank2T::T2m2(const vector<int>& HSs,int i,int j) const
  { return blow_up(Tm2, HSs, i, j); }
matrix IntRank2T::T2m(int m, const vector<int>&HSs, int i, int j) const
  { return blow_up(Tcomp(m), HSs, i, j); }  
matrix IntRank2T::CompSpace(const matrix& mx, 
                             const vector<int>&HSs, int i, int j) const
  { return blow_up(mx, HSs, i, j); }  

// ____________________________________________________________________________
// D           INTERACTION SPIN TENSOR SPIN CHECKING FUNCTIONS
// ____________________________________________________________________________

        // Input                IST	: Interaction Spin Tensor (this)
        //                      Iqn     : Spin quantum value of I
        // Output               TF      : True if I is an even multiple of 1/2

//------------------------------------------------------------------------------
//                        Isotope Designation Checks             
//------------------------------------------------------------------------------

// ____________________________________________________________________________
// P                          STANDARD I/O FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//     Functions That Generate Strings to Simplify and Modularize Printing
//-----------------------------------------------------------------------------

string IntRank2T::StringI() const
  {
  string SI;
  if(!Sval) SI = string("Spin Quantum Number:   ");
  else      SI = string("Spin I Quantum Number: ");
  int intI = int(Izval());                      // Get I as an integer
  double x = Izval()/double(intI);              // Get remainder of I/int(I)
  if(fabs(x-1.0)>1.e-6) intI = 1;               // Flag if I is whole number
  else                  intI = 0;               // or if I isn't whole
  if(intI)
    SI += string(Gform("%3i",int(2*Izval())))
       +  string("/2      ");
  else
    SI += string(Gform("%3i", int(Izval())))
       +  string("        ");
  return SI;
  }

string IntRank2T::StringS() const
  {
  string SS("Spin S Quantum Number: ");
  int intS = int(Szval());                      // Get S as an integer
  double x = Szval()/double(intS);              // Get remainder of S/int(S)
  if(fabs(x-1.0)>1.e-6) intS = 1;               // Flag if S is whole number
  else                  intS = 0;               // or if S isn't whole
  if(intS)
    SS += string(Gform("%3i",int(2*Szval())))
       +  string("/2      ");
  else
    SS += string(Gform("%3i", int(Szval())))
       +  string("        ");
  return SS;
  }

string IntRank2T::StringIS() const
  {
  int intI = int(Izval());                      // Get I as an integer
  double x = Izval()/double(intI);              // Get remainder of I/int(I)
  if(fabs(x-1.0)>1.e-6) intI = 1;               // Flag if I is whole number
  else                  intI = 0;               // or if I isn't whole

  int intS = int(Szval());                      // Get S as an integer
  double y = Szval()/double(intS);              // Get remainder of S/int(S)
  if(fabs(y-1.0)>1.e-6) intS = 1;               // Flag if S is whole number
  else                  intS = 0;               // or if S isn't whole

  string SIS("Spin Quantum Numbers:  ");
  if(intI)
    SIS += string(Gform("%3i",int(2*Izval())))
        +  string("/2, ");
  else
    SIS += string(Gform("%3i", int(Izval())))
        +  string(",   ");
  if(intS)
    SIS += string(Gform("%1i",int(2*Szval())))
        +  string("/2 ");
  else
    SIS += string(Gform("%1i", int(Szval())))
        +  string("   ");
  return SIS;
  }

//-----------------------------------------------------------------------------
//  Functions That Generate string Arrays to Simplify and Modularize Printing
//-----------------------------------------------------------------------------

/* These functions return spin tensor components in string format.  This is
   done to facilitate printing, in particular printing of spin tensors within
   rank 2 interactions.
 
	   Input		IST	: Interaction Spin Tensor
	  			M       : Ang. momentum component [0,4]
                                m       : Ang. momentum component [-2,2]
           Output               TSS     : Pointer to array of hs strings
	  				  where hs is the spin pair
	  				  Hilbert space dimension

     				[ x.x, x.x, x.x]
  		        T     = [ x.x, x.x, x.x]
  		    	 2,m	[ x.x, x.x, x.x]
  				[ x.x, x.x, x.x]

              M = { 0, 1, ..., 4 } <===> m = { 0, 1, -1, 2, -2 }             */

vector<string> IntRank2T::T2Strings(int m) const
  {
  switch(m)
    {
    case  0: return TStrings(0); break;
    case  1: return TStrings(1); break;
    case -1: return TStrings(2); break;
    case  2: return TStrings(3); break;
    case -2: return TStrings(4); break;
    default:                     break;
    }
  return TStrings(0);				// Should never get here
  }

vector<string> IntRank2T::TStrings(int M) const
  {
  int mvals[5]    = { 0, 1, -1, 2, -2 };	// Actual m values
  matrix T  = Tcomp(mvals[M]);			// Spin tensor component
  int hs    = T.cols();				// Spin array dimension
  vector<string> TTS(hs);			// An array of hs strings
  if(!hs) return  TTS; 				// Quit if null array

  string mlabs[5] = { "0 ", "1 ", "-1",		// Labels for m
                            "2 ", "-2" };
  int Treal = T.is_real();			// Real vs. complex flag
  int j;
  for(int i=0; i<hs; i++)			// Loop all array columns
    {
    if(i == hs/2-1)    TTS[i] = string("T    = ");
    else if(i == hs/2) TTS[i] = string(" 2,")
                              + mlabs[M]
                              + string("  ");
    else               TTS[i] = string("       ");
    TTS[i] += string("[");			// Begin array row i
    if(Treal)					// Loop over the columns
      {						// and add in elements
      for(j=0; j<hs; j++)
        TTS[i] += Gform("%7.3f", T.getRe(i,j));
      }
    else
      {
      for(j=0; j<hs; j++)
        TTS[i] += string("(")
                     +  string(Gform("%7.3f", T.getRe(i,j)))
                     +  string(", ")
                     +  string(Gform("%7.3f", T.getIm(i,j)))
                     +  string(")");
      }
    TTS[i] += string(" ]");			// End array row i
    }
  return TTS;
  }


//-----------------------------------------------------------------------------
//     Functions That Generate Ouput Of The Rank 2 Interaction Spin Tensor
//-----------------------------------------------------------------------------

ostream& IntRank2T::print(ostream& ostr, int fflag) const

        // Input                IST	: Interaction Spin Tensor (this)
        //                      ostr	: Output stream
	//			fflag   : Format flag
	//				    0 - Basic Parameters
	//				   !0 - Full output
        // Output               none    : Interaction Spin Tensor parameters 
        //                                placed into the output stream

// sosi - this should probably use the above strings once they have been tested

  {
  if(!Ival)
    {
    string hdr("Empty Interaction Spin Tensor");
    int hl = hdr.length();
    string spacer = string(40 - hl/2, ' ');
    ostr << "\n\n" << spacer << hdr << "\n";
    return ostr;
    }

//	       Now Output The Spin Tensor Spherical Components
//     These Will Be Printed Lines Which Will Look Like The Following 
//		      (Repeated For All 5 m Values)
//
// 				[ x.x, x.x, x.x]
//		        T     = [ x.x, x.x, x.x]
//		    	 2,m	[ x.x, x.x, x.x]
//				[ x.x, x.x, x.x]

  int hs = T0.rows();				// Get Hilbert space
  int wid = hs*6 + 11;				// Width of printed array
  int aperl = 70/wid;				// Arrays per line
  int mind[5] = { 0, 1, -1, 2, -2 };		// Access indices 
  string mlabs[5] = {"0 ","1 ","-1","2 ","-2"}; // Labels for m
  int reim[5];					// Flag if real arrays 
  int i, j, m, n;				// Array indices
  for(m=0; m<5; m++)				// Set real vs. complex
    reim[m] = (Tcomp(mind[m])).is_real();	// array flags
  string spacer = string((80-aperl*wid)/2, ' ');
  int left = 5;
  for(m=0; m<5; m+=aperl)
    {

//   Print All Array Rows Above Where the Label Will Be Written 

    for(i=0; i<hs/2-1; i++)			// Print array upper half
      {						// Index i is the row
      ostr << spacer;				//   Begin with space
      for(n=m; n<m+aperl && n<5; n++)		//   Loop m values if multiple
        {					//   arrays per line
        ostr << "         [";			//	Start array n, line i
        if(reim[n]) for(j=0; j<hs; j++) ostr << Gform("%7.3f", Tcomp(mind[n]).getRe(i,j));
        else for(j=0; j<hs; j++)        ostr << Tcomp(mind[n]).get(0,j);
        ostr << " ]";				//	End array n line i
        }					//   Next m if multi-array line 
      ostr << "\n";				// End row i
      }

//       Print Two Array Rows Where the Label Will Be Written 

    ostr << spacer;
    for(n=m; n<m+aperl && n<5; n++)		// Do row with T label
      {
      ostr << "  T    = [";
      if(reim[n]) for(j=0; j<hs; j++) ostr << Gform("%7.3f", Tcomp(mind[n]).getRe(i,j));
      else      for(j=0; j<hs; j++) ostr << Tcomp(mind[n]).get(0,i);
      ostr << " ]";
      left--;
      }
    ostr << "\n";
    ostr << spacer;
    i++;
    for(n=m; n<m+aperl && n<5; n++)		// Do row with m value
      {
      ostr << "   2," << mlabs[n] << "  [";
      if(reim[n]) for(j=0; j<hs; j++) ostr << Gform("%7.3f", Tcomp(mind[n]).getRe(i,j));
      else      for(j=0; j<hs; j++) ostr << Tcomp(mind[n]).get(0,i);
      ostr << " ]";
      }
    ostr << "\n";

//      Print All Array Rows Below Where the Label Will Be Written 

    for(i++; i<hs; i++)				// Loop lower rows
      {
      ostr << spacer;
      for(n=m; n<m+aperl && n<5; n++)
        {
        ostr << "         [";
        if(reim[n]) for(j=0; j<hs; j++) ostr << Gform("%7.3f", Tcomp(mind[n]).getRe(i,j));
        else      for(j=0; j<hs; j++) ostr << Tcomp(mind[n]).get(0,j);
        ostr << " ]";
        }
      ostr << "\n";
      }
    if(left < aperl) spacer = string((80-left*wid)/2, ' ');
    ostr << "\n\n";
    }
  return ostr;
  }

ostream& operator<< (ostream& out, const IntRank2T& IST)
  { return IST.print(out); }


// ____________________________________________________________________________
// F               RANK 2 TENSOR LIST/VECTOR SUPPORT FUNCTIONS
// ____________________________________________________________________________

/* Spin systems will often maintain linked lists of rank 2 spin tensors which
   coincide with the spin, spin-field, and spin-spin interactions present.
   Use of C++ STL list and vector classes can be used for such purposes if
   the ususal comparison functions are defined. Below are the definitions.

        // Input                IST     : Interaction Spin Tensor (this)
        //                      IST1    : Another interaction spin tensor
        // Output               T/F     : TRUE if IST is 1st or == IST1
        ///F_list ==                    - Equality
        ///F_list !=                    - Inequality                         */

bool IntRank2T::operator==(const IntRank2T &IST1) const
  { return (IST1.Ival==Ival && IST1.Sval==Sval); } 

bool IntRank2T::operator!=(const IntRank2T &IST1) const
  { return (IST1.Ival!=Ival || IST1.Sval!=Sval); } 

bool IntRank2T::operator<(const IntRank2T &IST1) const
  { 
  if(Ival < IST1.Ival) return true;
  else if(Ival == IST1.Ival && Sval < IST1.Sval) return true;
  return false;
  }

bool IntRank2T::operator>(const IntRank2T &IST1) const
  {
  if(Ival > IST1.Ival) return true;
  else if(Ival == IST1.Ival && Sval > IST1.Sval) return true;
  return false;
  }

int IntRank2T::lessthan(double Iqn, double Sqn) const

        // Input                IST     : Interaction Spin Tensor (this)
        //                      Iqn,Sqn : Spin quantum values of I&S
        // Output               TF      : True if Ist(Iqn,Sqn) spin
        //                                quantum values that are below
        //                                those input (Iqn,Sqn)

  {
  int newI = int(2.0*Iqn + 1);			// To compare Ival's
  if     (Ival < newI) return 1;		// True if IST(Iqn) < Iqn
  else if(Ival > newI) return 0;		// False if IST(Iqn) > Iqn
  else 						// Check S quantum values
    {						// if the I values are equal
    if(!Sqn) return 0;				// False if I's match & no S's
    int newS = int(2.0*Sqn + 1);		// S value to compare
    if(Sval < newS) return 1;			// True: I's match, IST(Sqn)<Sqn
    }
  return 0;
  }

 
int IntRank2T::morethan(double Iqn, double Sqn) const
 
        // Input                IST     : Interaction Spin Tensor (this)
        //                      Iqn,Sqn : Spin quantum values of I&S
        // Output               TF      : True if Ist(Iqn,Sqn) spin
        //                                quantum values that are above
        //                                those input (Iqn,Sqn)

  {
  int newI = int(2.0*Iqn + 1);			// To compare Ival's
  if     (Ival > newI) return 1;		// True if IST(Iqn) > Iqn
  else if(Ival < newI) return 0;		// False if IST(Iqn) < Iqn
  else 						// Check S quantum values
    {						// if the I values are equal
    if(!Sqn) return 0;				// False if I's match & no S's
    int newS = int(2.0*Sqn + 1);		// S value to compare
    if(Sval > newS) return 1;			// True: I's match, IST(Sqn)>Sqn
    }
  return 0;
  }
 
int IntRank2T::equalto(double Iv, double Sv) const
  {
  if(Iv != Ival) return 0;			// False if I's don't match
  if(!Sv)        return 1;			// True is I's match, no S's
  if(Sv != Sval) return 0;			// False if S's don't match
  return 1;					// True, both I's and S's match
  }
 
#endif						// IntRank2T.cc


/* Final Comments:

    1. This class knows about isotopes.  That is (almost exclusively) to allow
       for simple constructors from isotope type(s).  In turn, that forces
       us to perform isotope type checking to insure that electons are used
       in electron interactions and nuclei in nuclear spin interactions. 
       That's a pain, so I may dumb this class down before it is released.
    2. Object of type IntRank2T are not much use in and of themselves.  They
       are meant to be combined with spatial tensors (class IntRank2A) and
       interaction constants to treat spin and spin-spin interactions.  This
       is done in the derived class IntRank2.  Unless there is a specific
       need to deal directly with this class, users should use objects of type
       IntRank2 (or specific interaction classes derived from IntRank2) in 
       their GAMMA programs.
    3. This class defines the enumeration used for specification of particular
       interaction types.  However the type is NOT part of the class structure.
       That is, no spin tensor knows what type it is.  Rather the types are
       used during construction so that we generate the proper spin tensor
       components.                                                           */

