/* sys_dynamic.cc ***********************************************-*-c++-*-
**									**
**	                        G A M M A				**
**								 	**
**	Dynamic Spin System                     Implementation		**
**								 	**
**	Copyright (c) 1991, 1992, 1993				 	**
**	Scott Smith						 	**
**	Eidgenoessische Technische Hochschule	 			**
**	Labor fuer physikalische Chemie				 	**
**	8092 Zuerich / Switzerland				 	**
**								 	**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**								 	**
** The class sys_dynamic defines a collection nuclear spins with	**
** specific environments and motional properties.			**
**									**
** NOTE: This class is sheduled for a dramatic change before GAMMA 4.0	**
**       is released.  Rather than using class coord_vec as a base and	**
**       space_T classes for spin and spin pair tensor tracking, this	**
**       will switch over to use the newer Rank2Int class(es).  That	**
**       should be more efficient and coincide exactly with the new	**
**       solid state modules due to be in GAMMA soon.			**
**									**
*************************************************************************/

#ifndef _sys_dynamic_cc_		// Is file already included?
#define _sys_dynamic_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

const double mu     = 1.e-7;		// [mu/(4*pi)] J-sec -C  -m
const double rt6po5 = 1.941625913;	// sqrt[6*pi/5]

#include <GamGen.h>			// Include OS specific stuff
#include <LSLib/sys_dynamic.h>		// Include dynamic spin system def
#include <HSLib/SpinSystem.h>		// Includes spin system knowledge
#include <Level1/SpaceT.h>		// Includes spatial tensor knowledge
#include <Level1/coord.h>		// Includes coordinate knowledge
#include <Level1/coord_vec.h>		// Includes coordinate vector knowledge
#include <Matrix/matrix.h>		// Includes matrix knowledge
#include <Matrix/row_vector.h>		// Includes row vector knowledge
#include <Basics/Gconstants.h>		// Includes GAMMA constants
#include <Basics/Gutils.h>  		// Includes query knowledge
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Basics/StringCut.h>		// Include string parsing

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i              CLASS DYNAMIC SPIN SYSTEM ERROR HANDLING
// ____________________________________________________________________________


        // Input                dsys	: Dynamic spin system (this)
        //                      eidx    : Error index
	//			pname   : string for error message
        //                      noret	: Flag for return (0=no)
        // Output               none  	: Error message
        //                                Program execution stopped

void sys_dynamic::DSerror(int eidx, int noret) const
  {
  std::string hdr("Anisotropic System");
  std::string msg;
  switch(eidx)
    {
    case 0: GAMMAerror(hdr, 0, noret); break;   // Program Aborting        (0)
    case 3: GAMMAerror(hdr, 3, noret); break;   // !Construct from pset    (3)
    case 5: GAMMAerror(hdr, 5, noret); break;   // !Write to param. file   (5)
    case 9: msg = std::string("No Correlation Times Have Been Specified");
            GAMMAerror(hdr, msg, noret); break; //                         (9)
    case 10: msg = std::string("Inappropriate Correlation Time Specified");
            GAMMAerror(hdr, msg, noret); break; //                         (10)
    case 11: msg = std::string("Warning - No Quad. Coupling Allowed for I=1/2");
            GAMMAerror(hdr, msg, noret); break; //                         (11)
    case 12: msg = std::string("Warning - Setting Asymmetry of a Zero Tensor");
            GAMMAerror(hdr, msg, noret); break; //                         (12)
    case 13: msg = std::string("Cannot Determine # Spins In The System");
            GAMMAerror(hdr, msg, noret); break; //                         (13)
    case 14: msg = std::string("Sorry, Dipolar Tensor Operation Not Allowed Yet");
            GAMMAerror(hdr, msg, noret); break; //                         (14)
    case 15: msg = std::string("Attempted Dipole Access of Spin with Itself");
            GAMMAerror(hdr, msg, noret); break; //                         (15)
    case 21: msg = std::string("Can't Read System from Parameter File");
            GAMMAerror(hdr, msg, noret); break; //                         (21)
    case 22: msg = std::string("Unreasonable Internuclear Distance");
            GAMMAerror(hdr, msg, noret); break; //                         (22)
    case 30: msg = std::string("Problem(s) With Defined Exchange Process");
            GAMMAerror(hdr, msg, noret); break; //                         (30)
    case 31: msg = std::string("Negative Mutual Exchange Rate Specified");
            GAMMAerror(hdr, msg, noret); break; //                         (31)
    case 32: msg = std::string("Specified Heternuclear Spin Exchange");
            GAMMAerror(hdr, msg, noret); break; //                         (32)
    case 33: msg = std::string("Exchange Between Improperly Indexed Spins");
            GAMMAerror(hdr, msg, noret); break; //                         (33)
    case 34: msg = std::string("Accessed Exchange Process NonExistant");
            GAMMAerror(hdr, msg, noret); break; //                         (34)
    case 35: msg = std::string("Cannot Get/Set Exchange Rate");
            GAMMAerror(hdr, msg, noret); break; //                         (35)
    case 36: msg = std::string("Mutual Exchange Of Heteronuclei Disallowed!");
            GAMMAerror(hdr, msg, noret); break; //                         (36)
    case 37: msg = std::string("Invalid Mutual Exchange Process Declared!");
            GAMMAerror(hdr, msg, noret); break; //                         (37)
    case 44: msg = std::string("Warning - Unspecified Shifts Set to 0 Hz");
            GAMMAerror(hdr, msg, noret); break; //                         (44)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }

void sys_dynamic::DSerror(int eidx, const std::string& pname, int noret) const
  {
  std::string hdr("Anisotropic System");
  std::string msg;
  switch (eidx)
    {
    case 1:  GAMMAerror(hdr, 1, pname, noret);  break;	// File Problems  (1)
    case 2:  GAMMAerror(hdr, 2, pname, noret);  break;	// !Read Parameter(2)
    case 5:  msg = "Warning - No Chemical Shift Information: Spin " + pname;
             GAMMAerror(hdr, msg, noret);  break; 	// Bad Parameter  (5)
    case 6:  msg = "Parameter " + pname + " Is The Culprit!";
             GAMMAerror(hdr, msg, noret);  break;	// Bad Parameter  (6)
    case 7:  msg = "There Is No Exchange Process " + pname;
             GAMMAerror(hdr, msg, noret);  break;	// Bad Parameter  (7)
    default: GAMMAerror(hdr, -1, pname, noret); break;	// Unknown Error  (-1)
    }
  }

volatile void sys_dynamic::DSfatal(int eidx) const
  {
  DSerror(eidx, 1);				// Output error message
  if(eidx) DSerror(0,1);			// Write we're aborting
  GAMMAfatal();					// Clean exit from program
  }

volatile void sys_dynamic::DSfatal(int eidx, const std::string& pname) const
  {
  DSerror(eidx,pname,1);			// Output error message
  if(eidx) DSerror(0,1);			// Write we're aborting
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii             CLASS DYNAMIC SPIN SYSTEM SETUP FUNCTIONS
// ____________________________________________________________________________

/* These are protected functions because they allow specific aspects of the
   spin system to be set up without worrying about system consistency!       */



	// Input		sys      : Spin system (this)
	// 			pset     : A parameter set
	//			idx	 : Parameter index value used for
	//				   prefix [#] in input names
	//			warn	 : Warning output level
	//					0 = no warnings
	//					1 = warnings
	//				       >1 = fatal warnings
	// Output		TF	 : Spin system is filled with
	//				   parameters in pset
	// Note				 : Three things are gleaned from
	//				   the parameter set from spin_sys
	//				   1.) The number of spins
	//				   2.) An isotope type for each spin
	//				   3.) An optional spin system name
	// 				   Three more parameters are taken
	//				   from pset to complete spin_system
	//				   4.) The chemical shifts
	//				   5.) The coupling constants
	//				   6.) A spectrometer frequency
	// Note				 : Functions which place a spin_system
	//				   into a parameter set must contain
	//				   add the information read here

int sys_dynamic::setSsys(const ParameterSet& pset, int idx, int warn)
  {
//	Use Only Parameters With Prefix [#], So Clip [#] From Names First

  ParameterSet  subpset;			// Working parameter set
  if(idx != -1) subpset=pset.strip(idx);	// Get only params with [#]
  else          subpset=pset;			// Or use full pset

//	          Now Set Up The Spin System Directly

  int TF=1;					// Track if we read OK
  int ns = getSpins(subpset, warn?1:0);		// Get the number of spins
  if(ns<=0)					// If there are no spins
    {
    if(warn) DSerror(13,1);
    else     DSerror(13);       
    return 0;
    }
  *this = sys_dynamic(ns); 			// Set system for ns spins
  setIs(subpset);				// Set the isotope types
  setName(subpset);				// Read in the system name
  setBasis(matrix(HS(), HS(), i_matrix_type));	// Set up default basis (matrix)
  setOm(subpset);				// Set spectrometer frequency
  if(coord_vec::read(pset,-1,0)) setDip();	// Set the dipolar tensors 
  SetCSA(subpset);				// Set shift & shift anisotropy
  setJs(subpset); 				// Set J couplings (spin_system)
  setQuad(subpset);				// Set quadrupolar parameters
  setRand(subpset);				// Set random field parameters
  setTaus(subpset);				// Set the correlation time(s) 
  setKs(subpset, false);			// Set the exchange processes
  if(electrons())				// Read G's & J's if electrons
    {						// present in the spin system
    setGs(subpset);				// 	Set g-factors
    setAs(subpset);				// 	Set hyperfine couplings
    } 
  return TF;
  }


// ____________________________________________________________________________
// iii             CLASS DYNAMIC SPIN SYSTEM CHECKING FUNCTIONS
// ____________________________________________________________________________

bool sys_dynamic::CheckExch(int p, bool warn) const
  {
  if(p>=0 && p<int(MExs.size())) return true;
  if(warn) DSerror(34, 1);
  return false;
  }

bool sys_dynamic::CheckExch(ExchProcM& XP, bool warn) const
  {
  int np = XP.NComps();				// # of exchange processes
  if(!np) return true;				// If none, XP is OK
  Isotope I = isotope(XP.Comp(0));		// First exchanging spin type
  for(int i=1; i<np; i++)			// Loop through all others
    if(isotope(XP.Comp(i)) != I)		// exchanging and return false
      { if(warn) DSerror(36,1); return false; }	// if any are not the same type
  return true;
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A              DYNAMIC SPIN SYSTEM CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

sys_dynamic::sys_dynamic() : spin_system(0), coord_vec(0) { Taus.xyz(0,0,0); }

sys_dynamic::sys_dynamic(int ns) : spin_system(ns), coord_vec(ns)
  {
  Taus.xyz(0,0,0);				// Set correlation times to zero
  if(ns)					// If spins set, allocate space
    {
    shift_As = std::vector<space_T> (ns);	// Shift tensors for each spin
    quad_As  = std::vector<space_T> (ns);	// Quadrupole tensors, each spin
    if(ns>1) 					// Dipole tensors, per spin pair
      {
      int nd = (ns*ns-ns)/2;			//   Number of dipoles
      dip_As = std::vector<space_T> (nd);	//   Set up spatial tensors
      }
    rand_As = std::vector<double>(ns+1, 0.0);	// Rand. Fld."tensors", per spin
    } 						// Insure RF "tensors" zeroed
  }

sys_dynamic::sys_dynamic(const sys_dynamic &dsys1)
            : spin_system(dsys1), coord_vec(dsys1)
  {
  shift_As = dsys1.shift_As;			// Copy the shift tensors
  dip_As   = dsys1.dip_As;			// Copy the dipolar tensors
  quad_As  = dsys1.quad_As;			// Copy the quadrupolar tensors
  rand_As  = dsys1.rand_As;			// Copy random field "tensors"
  Taus     = dsys1.Taus;			// Copy the correlation times
  MExs     = dsys1.MExs;			// Copy mutual exchange processes
  }

sys_dynamic& sys_dynamic::operator= (const sys_dynamic& dsys)
  {
  if(this == &dsys) return *this;	// Avoid self copying
  spin_system::operator=(dsys);		// Equate the spin system part
  coord_vec::operator=(dsys);		// Equate the coord vector part
  shift_As = dsys.shift_As;		// Copy chemical shift tensors
  dip_As   = dsys.dip_As;		// Copy the dipolar tensors
  quad_As  = dsys.quad_As;		// Copy the quadrupolar tensors
  rand_As  = dsys.rand_As;		// Copy random field "tensors"
  Taus = dsys.Taus;			// Copy correlation times 
  MExs = dsys.MExs;			// Copy mutual exchange processes
  return *this;
  }

sys_dynamic::~sys_dynamic () { }

// ____________________________________________________________________________
// B                      CHEMICAL SHIFT MANIPULATIONS
// ____________________________________________________________________________

// ************* Chemical Shift Manipulations in Hertz & PPM ******************
 
/* These functions allow users to either get or set the isotropic chemical
   shift of a spin or set of spins. These functions could be handled entirely
   by the base class where it not for the shift tensors in this class.  Since
   the shift tensors MUST coincide with the isotropic values, sys_dynamic must
   mirror the shift functions of class spin_system.  There are a few important
   items worth taking note of.
    
   1.) Electrons CANNOT have their "shift" set with these functions.  Users
       must set the electron g factor to get the ESR semi-equivalent for e-.
   2.) GAMMA assumes system shifts are stored relative to a base Larmor
       frequency particular for the associated spin's isotope type.  For a
       positive gyromagnetic ratio, a spin with a positive shift is deshielded
       relative to what the base Larmor frequency is.
   3.) Shifts in PPM may not be set/obtained unless a spectrometer field
       strength has been specified (Omega).

       Function    Return   Arguments              Results
      ----------   ------   ---------   --------------------------------------
       shift       void      spin,w     Set spin shift to w spin = [0,ns)
       shift       double     spin      Retrieve shift of spin in Hz
       lab_shift   double     spin      Retrieve shift of spin in Hz, lab frame
       PPM         spin      spin,PPM   Set spin shift to PPM
       PPM         void       spin      Retrieve shift to spin in PPM
       shifts      double       w       Set all nuclear spin shifts to w in Hz
       maxShift    double               Retrieve largest shift in Hz
       maxShift    double      iso      Get largest shift in Hz, iso spins only
       minShift    double               Retrieve smallest shift in Hz
       minShift    double      iso      Get smallest shift Hz, iso spins only
       medianShift double               Get middle of shifts in Hz           */

void sys_dynamic::shifts(double nu)
  {
  spin_system::shifts(nu);			// Set isotopic shifts
  int s=spins();
  for (int i=0; i<s; i++)			// Set tensor values too
    {						//  (these are in PPM)
    nu = spin_system::PPM(i);				
    if((shift_As[i]).exists())
    (shift_As[i]).iso(nu);
    }
  }

void sys_dynamic::shift(double nu, int spin) { shift(spin,nu); }
void sys_dynamic::shift(int spin, double nu)
  {
  spin_system::shift(spin, nu);			// Set isotropic shift value
  if((shift_As[spin]).exists())			// Also set the tensor value
    {						//  (this is in PPM)
    nu = spin_system::PPM(spin);
    (shift_As[spin]).iso(nu);
    }
  }

        // Input                nu     : chemical shift offset value (Hertz)
	//			spin   : Spin label
	// Output 		none   : Chemical shifts modified by offset
	//				 for all spins of same isotope type
	//				 as the input spin

double sys_dynamic::shift(int spin) const { return spin_system::shift(spin); }
void   sys_dynamic::offsetShifts(double nu, int spin)
  {
  spin_system::offsetShifts(nu, spin);		// Offset all isotropic shifts
  int s=spins();				// Number of spins
  Isotope I = isotope(spin);			// First spin type
  for(int j=0; j<s; j++)			// Offset all tensor values
    if((shift_As[j]).exists() && I==isotope(j)) //  (these are done in PPM)
      {
      nu = spin_system::PPM(j);
      (shift_As[j]).iso(nu);
      }
  }


void sys_dynamic::offsetShifts(double nu, const std::string& Iso)

        // Input                nu     : chemical shift offset value (Hertz)
	//			spin   : spin label
	// Output 		none   : chemical shifts modified by offset
	//				 for all spins of same isotope type
	//				 as the input spin

    {
    spin_system::offsetShifts(nu, Iso);		// Offset all isotropic shifts
    int s=spins();
    for(int i=0; i<s; i++)			// Offset all tensor values
      if((shift_As[i]).exists()			//  (these are done in PPM)
	              && Iso == symbol(i))
        {
        nu = spin_system::PPM(i);
        (shift_As[i]).iso(nu);
        }
    }

// ******************* Chemical Shift Manipulations in PPM ********************


void sys_dynamic::PPM(int spin, double freq)

	// Input		spin   : spin of spin system [0, nspins-1]
	// 			freq   : chemical shift of spin (PPM)
	// Output		none   : set chemical shift of spin (PPM)

  {
  spin_system::PPM(spin, freq);		// Set all isotropic value
  if((shift_As[spin]).exists())		// Also set the tensor value
  (shift_As[spin]).iso(freq);
  }


// ********************* Chemical Shift Tensor Functions **********************


double sys_dynamic::delz(int spin) const

	// Input		spin   : spin of spin system [0, nspins-1]
	// Output		none   : delz component of chemical shift
	//				 spatial tensor of spin (PPM)

  { if((shift_As[spin]).exists()) return shift_As[spin].delz(); return 0; }



	// Input		spin   : spin of spin system [0, nspins-1]
	//			delzz  : value of delzz in PPM
	// Output		none   : delz component of chemical shift
	//				 spatial tensor of spin set to
	//				 input value of delzz

void sys_dynamic::delz(int spin, double delzz)
    {
    double Aiso = 0;
    if((shift_As[spin]).exists())	// If CSA tensor exists set delzz
      shift_As[spin].delz(delzz);
    else				// If CSA tensor does not exist
      {				        // construct it appropriately
      if(Omega() != 0.0)
        Aiso = spin_system::PPM(spin);	// Keep all in PPM				
      shift_As[spin] =  A2(Aiso, delzz);
      }
    return;
    }


double  sys_dynamic::Ceta(int i)  const       { return shift_As[i].eta(); }
void    sys_dynamic::Ceta(int i, double ceta) { shift_As[i].eta(ceta); }
space_T sys_dynamic::TC(int spin) const       { return shift_As[spin]; }
void    sys_dynamic::TC(const space_T& A, int spin)
// sosi - should check the tensor rank too.
  {
  check_spin(spin);
  shift_As[spin] = A;			// Set shift tensor
  PPM(spin, A.iso());			// Set isotropic shift
  }


row_vector sys_dynamic::xiC_vector( ) const

	// Input		dsys  : A dynamic system (this)
	// Return		dximx : A matrix of CSA interaction
	//				constants (rad/sec)
	// Note			      : The PPM units of delzz are
	//				counteracted by the MHz units
	//				on Omega

/*		      1/2
	  CSA   [6*pi]
	xi    = |----| * gamma * B * del  (i) = K * Omega  * del  (i)
	  i     [ 5  ]        i   0     zz               i      zz           */

  {
  double delzz, xii;
  double Omi;				// Spectrometer frequency (MHz)
  int ns = spins();
  row_vector xivec(ns);
  for(int i=0; i<spins(); i++)
    {
    Omi = Omega(i);			// Frequency of spin i (MHz)
    delzz = delz(i);			// delzz of spin i (PPM)
    xii = rt6po5*PIx2*Omi*delzz;	// 2*pi converts Hz to rad/sec
    xivec.put(xii, i);			// Store the Xi value for spin i
    }
  return xivec;
  }

double sys_dynamic::xiC(int i) const

	// Input		dsys  : A dynamic system (this)
	//			i     : Spin index
	// Return		xii   : The gamma CSA interaction
	//				constant (rad/sec)
	// Note			      : The PPM units of delzz are
	//				counteracted by the MHz units
	//				on Omega
	//
	//		 1/2
	//   CSA   [6*pi]
	// xi    = |----| * gamma * B * del  (i) = K * Omega  * del  (i)
	//   i     [ 5  ]        i   0     zz               i      zz

  {
  double Omi = Omega(i);		// Frequency of spin i (MHz)
  double delzz = delz(i);		// delzz of spin i (PPM)
  double xii = rt6po5*PIx2*Omi*delzz;	// 2*pi converts Hz to rad/sec
  return xii;
  }

	// Input		dsys  : A dynamic system (this)
	// Return		TF    : Return true if any shift tensors
	//				are present in the system

bool sys_dynamic::CSA( ) const
  {
  bool TF = false;			// Set for no shift tensors
  if(!shift_As.size()) return TF;	// Should be NULL if non
  int ns = spins();			// Loop over the system spins
  for(int i=0; i<ns && !TF; i++)	// and set TF true if any exist
    if((shift_As[i]).exists()) TF=true;
  return TF;
  }

// ____________________________________________________________________________
// C                        DIPOLAR MANIPULATIONS
// ____________________________________________________________________________


// -------------------- Spin Coordinate Access Functions ----------------------
 
void sys_dynamic::coords(const coord_vec& cvec, double cutoff)

        // Input                cvec   : A vector of coordinates
        //                      cutoff : Dipolar cutoff distance
        // Output               none   : Sets the spin coordinates
        //                               and generates dipolar tensors
        //                               for all spin pairs
        // Note                        : If the distance between two
        //                               spins is > cutoff then the dipolar
        //                               tensor is left zero.

  {
// sosi this is new. probably should keep distances >0
cutoff=5;
  int ns = spins();			// Loop over the system spins
  for(int i=0; i<ns; i++)		// Set all the spin coordinates
    put(cvec.get(i), i);		// assuming cvec units are meters
  setDip();				// Fill in all the dipolar tensors
  }
        // Output (Coord)	TF	: Return true if any spin coordinates
        //				  are present in the system

bool sys_dynamic::Coord() const { return coord_vec::size()?true:false; } 

void   sys_dynamic::DCC(int i, int j, double nu) { Ddelz(i, j, nu); }
double sys_dynamic::DCC(int i, int j)   const    { return Ddelz(i, j); }
double sys_dynamic::Ddelz(int i, int j) const
                                  { return dip_As[dipole(i,j)].delz(); }
void   sys_dynamic::Ddelz(int i, int j, double delzz)
  {
  if(spins() !=2) DSfatal(14);
  dip_As[dipole(i,j)] = A2(0,delzz);
  }


double sys_dynamic::Deta(int i, int j) const

	// Input		i      : Spin of spin system [0, nspins-1]
	// 			j      : Spin of spin system [0, nspins-1]
	// Output		none   : Return asymmetry of dipolar
	//				 coupling spatial tensor of spin pair

  { return dip_As[dipole(i,j)].eta(); }


void sys_dynamic::Deta(int i, int j, double Deta)

	// Input                i      : Spin of spin system [0, nspins-1]
	// 			j      : Spin of spin system [0, nspins-1]
	// 			eta    : Dipolar asymmetry of spin pair
	// Output		none   : Set dipolar asymmetry of spin pair

  { dip_As[dipole(i,j)].eta(Deta); }


space_T sys_dynamic::AD(int i, int j) const { return dip_As[dipole(i,j)]; }

	// Input		dsys   : Dynamic spin system (this)
	//      		i      : Spin of spin system [0, nspins-1]
	// 			j      : Spin of spin system [0, nspins-1]
	// 			spin   : Spin index
	// Output		sphT   : Spin dipolar spatial tensor
	//				 between spin1 and spin2


space_T sys_dynamic::AD(int dip) const { return dip_As[dip]; }

	// Input		dsys   : Dynamic spin system (this)
	//      		dip    : Dipole index
	// Output		sphT   : Dipolar spatial tensor for the dipole


int sys_dynamic::dipoles() const { return spinpairs(); }

	// Input		dsys     : Dynamic spin system (this)
	// Output		int 	 : Number of dipoles (spin pairs)
	//				   in the dynamic spin system

int sys_dynamic::dipole(int spin1, int spin2) const

	// Input		dsys     : Dynamic spin system (this)
	//      		spin1    : Spin of spin system [0, nspins-1]
	// 			spin2    : Spin of spin system [0, nspins-1]
	// Output		int 	 : Dipole index for the spin pair

   {
   if(spin1 == spin2) DSfatal(15);// Can't do spin with itself
   int s1 = spin1;			// Copy spin indices
   int s2 = spin2;
   if(spin1 > spin2)			// Insure s1 < s2
     {
     s1 = spin2;
     s2 = spin1;
     }
   int dip = 0;				// Set the dipole count to zero
   for(int i=0; i<s1; i++)		// Add dipoles from previous rows
     dip += spins()-1-i;
   dip += s2 - s1 -1;			// Add previous dipoles this row
   return dip;
   }


matrix sys_dynamic::xiD_matrix() const
//return dximx(spins(), spins(), complex0)

	// Input		dsys  : A dynamic system (this)
	// Return		dximx : A matrix of dipolar interaction
	//				constants (unit rad/sec)
	//					  -1      -2
	// Note			      : 1T = 1 J-C  -sec-m	

/*	  	          1/2
	  	    [6*pi]     mu 
	       -2 * |----|   * --- * hbar * gamma  * gamma
	  D         [ 5  ]     4pi               i        j
	xi   = ____________________________________________
	  ij		           3
	  		          r
	  		           ij                                        */
  {
  matrix dximx(spins(), spins(), complex0);
  const double K = -2*rt6po5*mu*HBAR;
  double xiij, rij, r3;
  double gi, gj;
  int npts = spins();
  for(int i=0; i<npts-1; i++)
    {					//		      -1  -1
    gi = gamma(i);			// gamma of spin i sec  -T
    for(int j=i+1; j<npts; j++)
      {					//		      -1  -1
      gj = gamma(j);			// gamma of spin j sec  -T
      rij = Rad(get(i), get(j));	// dist. between i & j (meters)
      r3 = rij*rij*rij;
      xiij = K*gi*gj/r3;
      dximx.put(xiij, i, j);
      dximx.put(xiij, j, i);
      }
    }
  return dximx;
  }




	// Input		dsys  : A dynamic system (this)
	// Return		TF    : Return true if any dipolar tensors
	//				are present in the system

bool sys_dynamic::Dip( ) const
    {
    bool TF = false;			// Set for no dipolar tensors
    if(!dip_As.size()) return TF;	// Should be NULL if none exist
    int ns = spins();			// Loop over the system spins
    for(int i=0; i<ns-1 && !TF; i++)	// and set TF true if any exist
      for(int j=i+1; j<ns && !TF; j++)	//
        if((dip_As[dipole(i,j)]).exists()) TF = true;
   return TF;
   }

// ____________________________________________________________________________
// D                       QUADRUPOLAR MANIPULATIONS
// ____________________________________________________________________________

/* These functions allow access to the system quadrupolar interactions. 
   Functions exist to obtain/set the quadrupolar coupling constant(s),

	  		 1/2   QCC   	        1/2  del  (i)
	     Q     [6*pi]         i       [6*pi]        zz
	   xi    = |----| * ----------  = |----|  * ----------
	     i     [ 5  ]   2I (2I -1)    [ 5  ]    2I (2I -1)
	                      i   i                   i   i      

             Function                  Action
             --------    ----------------------------------------------
               QCC       Get/Set Specified Spin QCC Value In Hz
               Qdelz     Get/Set Specified Spin QCC Value In Hz (= QCC)      */

void   sys_dynamic::QCC(int spin, double nu) { Qdelz(spin, nu); }
double sys_dynamic::QCC(int spin)   const    { return Qdelz(spin); }
double sys_dynamic::Qdelz(int spin) const
  { if((quad_As[spin]).exists()) return quad_As[spin].delz(); return 0; }


void sys_dynamic::Qdelz(int spin, double delzz)
  {
  if((quad_As[spin]).exists())		// If Q tensor exists set delzz
      quad_As[spin].delz(delzz);	// which is stored in kHz
  else if(qn(spin) > 0.5)		// If Q tensor does not exist and I>1/2
      quad_As[spin] = A2(0,delzz);	// construct one with only delzz !=0
  else DSerror(11);			// Warn of quad. setting on I=1/2 spin
  }


double sys_dynamic::Qeta(int spin) const
  { if((quad_As[spin]).exists()) return quad_As[spin].eta(); return 0; }

void sys_dynamic::Qeta(int spin, double Qeta)

	// Input		spin   : Spin of spin system [0, nspins-1]
	// 			eta    : Quadrupolar asymmetry of spin
	// Output		none   : Set quadrupolar asymmetry of spin

  {
  if((quad_As[spin]).exists())	// If Q tensor exists set delzz
    quad_As[spin].eta(Qeta);
  else if(qn(spin) > 0.5)	// If Q tensor does not exist and I>1/2
    DSerror(12);		// warn that setting eta of zero tensor
  else DSerror(11);		// Warn of quadrupolar setting I=1/2 spin
  return;
  }


space_T sys_dynamic::TQ(int spin) const { return quad_As[spin]; }
void    sys_dynamic::TQ(const space_T& A, int i) {check_spin(i); quad_As[i]=A;}

row_vector sys_dynamic::xiQ_vector( )
//return xivec(spins(), complex0)

	// Input		dsys  : A dynamic system (this)
	// Return		dximx : A matrix of quadrupolar interaction
	//				constants (rad/sec)

  {
  row_vector xivec(spins(), complex0);
  double I=0.5, delzz, xii, Ifact;
  for(int i=0; i<spins(); i++)
    {
    xii = 0;					// Set the Xi value to zero
    I = qn(i);					// Get the spin quantum number
    if(I > 0.5)					// Only mess with I>1/2 spins
      {
      delzz = QCC(i);				// delzz of spin i (Hz)
      Ifact = 2.0*I*(2.0*I - 1.0);		// I based denomenator
      xii = rt6po5*PIx2*delzz/Ifact;		// 2*pi converts Hz to rad/sec
      }
    xivec.put(xii, i);				// Store the Xi value for spin i
    }
  return xivec;
  }


double sys_dynamic::xiQ(int i)

	// Input		dsys  : A dynamic system (this)
	//			i     : Spin index
	// Return		xii   : Quadrupolar int. constant (rad/sec)

  {
  double xii = 0;
  double I = qn(i);
  double delzz=0;
  double Ifact=0;
  if(I > 0.5)
    {
    delzz = QCC(i);				// delzz of spin i (Hz)
    Ifact = 2.0*I*(2.0*I - 1.0);		// I based denomenator
    xii = rt6po5*PIx2*delzz/Ifact;		// 2*pi converts Hz to rad/sec
    }
  return xii;
  }




        // Input                dsys  : A dynamic system (this)
        // Return               TF    : Return true if ANY quadrupolar tensors
        //                              are present in the system

bool sys_dynamic::Quad( ) const
  {
  if(!quad_As.size()) return false;	// Return 0 if no quad. tensors at all
  bool TF = false;			// Set for no quadrupolar tensors
  int ns = spins();			// Loop over the system spins
  for(int i=0; i<ns && !TF; i++)	// and set TF true if any exist
      if((quad_As[i]).exists()) TF=true;
  return TF;
  }


// ____________________________________________________________________________
// E                       RANDOM FIELD MANIPULATIONS
// ____________________________________________________________________________

/* These functions allow access to the system random field interactions. Note
   that such interactions are more of a "catch all" for generating relaxation
   effects.  They are treated as rank 1 interactions and basically scaled to
   produce a specified single quantum transition linewidth constribution.
   In GAMMA, these "tensors" are just double precision numbers which represent
   SQT linewidths, i.e. relaxation under a random field interaction should
   produce a the specified linewidth(s).  Similarly, the random field inter-
   action constant is given by
	  	     	                            1/2
	               R              [    R       ]
	             xi    = 2 * pi * |LWhh / J(w )| 
	               i              [    i     i ]                         */

double sys_dynamic::TR(int i) const {if(rand_As.size()) return rand_As[i]; return 0;}
double sys_dynamic::tauR()    const { return rand_As[spins()]; }

row_vector sys_dynamic::xiR_vector( ) const

	// Input		dsys  : A dynamic system (this)
	// Return		dximx : A matrix of random field interaction
	//				constants (rad/sec)
	// Note			      : In the random field treatment, the
	//				xi value is defined in terms of single
	//				quantum transition linewidths

  {
  row_vector xivec(spins(), complex0);
  for(int i=0; i<spins(); i++) xivec.put(xiR(i), i);
  return xivec;
  }


double sys_dynamic::xiR(int i) const
  {
  double w0  = Omega(i)*PIx2*1.e6; 	// SQT frequency in rad/sec
  double tau = tauR();			// Get system random field tau value
  if(!tau) tau = 1.e-15; 		// If 0, set extreme narrowing (1 fsec)
  double Jw0 = tau/(1+tau*tau*w0*w0);	// Get general J(w0) = 4pi*Jred(w0)
  return sqrt(PI*rand_As[i]/(2.0*Jw0)); // Return xii value
  }

// ____________________________________________________________________________
//                           PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//        Functions To Make A Parameter Set From A Dynamic Spin System
// ----------------------------------------------------------------------------

sys_dynamic::operator ParameterSet( )

	// Input		dsys  : Dynamic spin system
	//  			pset  : Parameter set
	// Output		pset  : Parameter set with only
	//			        dynamic spin system parameters

   { ParameterSet pset; pset += *this; return pset; }


void operator+= (ParameterSet& pset, sys_dynamic &dsys)

	// Input		dsys  : Dynamic spin system
	//  			pset  : Parameter set
	// Output		pset  : Parameter set with only
	//			        dynamic spin system parameters

  {
  pset += (spin_system)dsys;		// First add spin_system parameters
  pset += (coord_vec)dsys;		// Also add coord_vec parameters

   std::string pname;			// Now add unique sys_dynamic parameters
   std::string pstate;
   SinglePar par;
   int ns = dsys.spins();
   for (int i=0; i<ns; i++)		// Add the shift tensors
     {
     if((dsys.shift_As[i]).exists())
       {
       pstate = std::string("Chemical Shift Tensor");
       pname = std::string("CS_T(");
       pname += std::string(Gdec(i));
       pname += std::string(")");
       par = dsys.shift_As[i].param(pname,pstate);
       pset.push_back(par);
       }
     }
   pstate = "Correlation Times";	// Correlation times (nsecs)
   pname = "Taus";
   coord tpersec = 1.e9*dsys.taus();
   par = tpersec.param(pname,pstate);
   pset.push_back(par);
   return;
   } 

// ----------------------------------------------------------------------------
//        Functions To Make A Dynamic Spin System From A Parameter Set
// ----------------------------------------------------------------------------


int sys_dynamic::setCoords(const ParameterSet& pset, int mand)

	// Input		dsys	: Dynamic spin system (this)
        //                      pset	: A parameter set
	//			mand	: Flag whether coords mandatory
        // Output               T/F     : True if spin coordiantes are
        //                                found in the pset.  The system
        //                                spin coordinates are set from
        //                                the parameters in pset
	// Note				: These be mandatory if mand!=0

  {
  std::string pnames[6];			// Begin possible param. names
  pnames[3] = "Coord(";				// - for input in Angstroms
  pnames[4] = "Coordn(";			// - for input in nanometers
  pnames[5] = "Coordm(";			// - for input in meters
  SinglePar par;				// A single parameter
  ParameterSet::const_iterator item;	// A pix into parameter list
  int ccount = 0;				// Parameter count
  coord pt;					// A coordinate point
  int i, j;
  for(i=0; i<spins(); i++)			// Retrieve the spin coordinates
    {
    for(j=0; j<3; j++)				// Construct possible parameter
      pnames[j] = pnames[3+j] + Gdec(i) + ")";	// names for this spin
    j = 0;
    par = SinglePar(pnames[j]);
    while((j<3) && (pset.seek(par)==pset.end()))
       j++;
    par = SinglePar(pnames[j]);
    item = pset.seek(par);		// Pix parameter list for Omega 
    switch(j)
      {
      case 0:				// Coord input in Angstroms,
        pt = coord(*item);		// keep in meters
        put(pt*1.e-10, i);
        ccount++;
        break;
      case 1:				// Coord input in nanometers,
        pt = coord(*item);		// keep in meters
        put(pt*1.e-9, i);
        ccount++;
        break;
      case 2:				// Coord input in meters
        pt = coord(*item);
        put(pt,i);
        ccount++;
        break;
      case 3:				// No coordinate specified
       default: 			// for spin i
        if(mand)
          {
          DSerror(2, pnames[0]); 	// Can't read in Coord(i)
          DSfatal(3);		// Can't construct system 
          }
       put(0,0,0,i);
       }
     }
  return ccount;
  }
 

void sys_dynamic::setDip( )

	// Input		dsys	: Dynamic spin system (this)
        // Output               none    : Spin system dipolar
	//				  relaxation values are set from
	//				  coordinates in dsys

// Note: Unlike shielding and quadrupolar tensors which are read in directly
//       from their PAS components in the spin system file (using the class
//       space_T), dipolar tensors may be set from spin coordinates instead.

  {
  const double K = mu*HBAR;			// muo*hbar/(4*pi)
  matrix the_mx = coord_vec::thetas();		// Get the dipole theta angles
  matrix phi_mx = coord_vec::phis();		// Get the dipole phi angles
  matrix d_mx = coord_vec::distances();		// Get the dipole distances
  int dip=0;					// Dipole counter
  double rad = 0.0;				// For dipolar distance(m)
  double dcc = 0.0;				// For dipolar coupling constant
  double alpha=0.0, beta=0.0;			// For dipole Euler angles
  int ns = spins();
  for(int i=0, j=0; i<ns-1; i++)		// Loop over possible dipoles
    for(j=i+1; j<ns; j++)
      {
      rad = d_mx.getRe(i,j);			// Dipole internuclear distance
      if((rad) < 1.e-11)			// Check radius is reasonable
           DSfatal(22);			//   (cutoff at 0.1 Angstrom)
      alpha = phi_mx.getRe(i,j);		// Dipole polar ang. over from x
      beta = the_mx.getRe(i,j);			// Dipole polar ang. down from z
      dcc = gamma(i)*gamma(j)*K/(rad*rad*rad);	// Dipole coupling constant
      dip_As[dip] = A2(0.0, dcc, 0., 		// Set dip. tensors, each pair
                             alpha, beta, 0.);	// with Diso, eta, & chi all 0
// sosik - this A2 constructor is not setting the angles properly!
      dip++;					// Increment the dipole count
      }
// Read in the Dipolar Tensors -------------------------
// sosi - not sure if I'll ever implement this.
//  Better done in a solids computation?
  }
 

void sys_dynamic::SetCSA(const ParameterSet& pset)

	// Input		dsys	: Dynamic spin system (this)
        //                      pset	: A parameter set
        // Output               none    : Spin system shift & shift anisotropy
	//				  relaxation values are set from
	//				  parameters in pset

  {
  SinglePar par;
  ParameterSet::const_iterator item;
  std::string pname;
  std::string pstate;
  std::string sval;
  int i, j;

			// Read in the Chemical Shift Information --------------

   SetFlags(0);					// Set all spin flags false
   std::string pnames[10];
   pnames[3] = std::string("CS_T(");			// Start of shift tensor parameter
   pnames[4] = std::string("PPM(");			// Start of PPM shift parameter
   pnames[5] = std::string("v(");			// Start of Hz shift parameter
   SinglePar pars[5];
   double pdatad;
   int count = 0;
   int ns = spins();				// Number of spins in system
   for(i=0; i<ns; i++)				// Loop through all system spins
     {
     shift(i, 0.0);				// Initialize the shift to 0 Hz
     for(j=0; j<3; j++)				// Construct possible parameter
       pnames[j] = pnames[3+j] + Gdec(i) + ")";	// names for the shift of this spin
     pars[0] = SinglePar(pnames[0]);
     pars[1] = SinglePar(pnames[1]);
     pars[2] = SinglePar(pnames[2]);
     if(pset.seek(pars[0]) != pset.end())
       {					// Allowed only if Omega
       item = pset.seek(pars[0]);
       shift_As[i] = space_T(*item);		// Set shift tensors, each spin
       spin_system::PPM(i,(shift_As[i]).iso());	// Assumed the tensor is in PPM
       SetFlag(i,1);
       if(!Omega())
         {
         DSerror(4, pnames[0]);		// Can't do without Omega
         DSfatal(3);			// Can't construct from pset
         }
       }
     else if(pset.seek(pars[1]) != pset.end())	// Then attempt to read
       { 					// isotropic value in PPM
       par = SinglePar(pnames[1]);
       item = pset.seek(par);
       (*item).parse(pname,pdatad,pstate);
       spin_system::PPM(i, pdatad);
       SetFlag(i,1);
       if(!Omega())
         {
         DSerror(4, pnames[1]);		// Disallowed without Omega
         DSfatal(3);			// Can't construct from pset
         }
       }
     else if(pset.seek(pars[2]) != pset.end())	// Then attempt to read
       { 					// isotropic value in Hz
       item = pset.seek(pars[2]);
       (*item).parse(pname,pdatad,pstate);
       spin_system::shift(i, pdatad);
       SetFlag(i,1);
       }
     else					// Count the number of
       count++; 				// spins with no shift set
     }
   if(count)					// Loop all the spins in system,
     {						// output warnings for all with
     count = 0;					// no shifts specified
     for(i=0; i<ns; i++)
       if(!GetFlag(i))
         {
         count++;
         if(count==1) DSerror(5, Gdec(i));	// No shift info for spin: i
         else
           {
           std::cout << "; Spin " << i;		// Neither for this spin i
           if(count == 5) count = 0;
           }
         }
     DSerror(44);			// Warning, Unspecified Shifts at 0 Hz
     }
  SetFlags(0);				// Leave system with all spin flags false
  }


void sys_dynamic::setQuad(const ParameterSet& pset)

	// Input		dsys	: Dynamic spin system (this)
        //                      pset	: A parameter set
        // Output               none    : Spin system quadrupolar
	//				  relaxation values are set from
	//				  parameters in pset

  {
  double I;
  int Qtensor = 0;
  std::string qnames[8];
  qnames[4] = "Q_T(";				// Start name for quad tensor
  qnames[5] = "QCC(";				// Start name quad. couple (Hz) 
  qnames[6] = "QCCk(";				// Start name quad. couple (kHz)
  qnames[7] = "QCCM(";				// Start name quad. couple (MHz)
  int ns = spins();				// Spins in the system
  SinglePar pars[8], par;
  ParameterSet::const_iterator item;
  std::string pname, pstate;
  double pdatad;
  for(int i=0, j=0; i<ns; i++)			// Loop all the spins in system
    {
    I = qn(i);					// Get the spin quantum number
    if(I > 0.5)				// Only do spins with I>1/2
      {
      QCC(i, 0.0);				// Initialize QCC to zero Hz
      for(j=0; j<4; j++)			// Construct possible
        { 					// parameter names and
        qnames[j] = qnames[j+4] + Gdec(i) + ")";// parameters
        pars[j] = SinglePar(qnames[j]);
        }
//      if(pset.seek(SinglePar(qnames[0])))	// First try to read QCC tensor
      if(pset.seek(pars[0]) != pset.end())	// First try to read QCC tensor
        {
        par = SinglePar(qnames[0]);
        item = pset.seek(par);
        quad_As[i] = space_T(*item);		// Set quad. tensors, each spin
        QCC(i,((quad_As[i]).delz())*1.e3);	// Assumed that the tensor
        } 					// QCC is read in kHz!
//      else if(pset.seek(SinglePar(qnames[1]))) // Next try to read QCC in Hz
      else if(pset.seek(pars[1]) != pset.end())	// Next try to read QCC in Hz
        {
        par = SinglePar(qnames[1]);
        item = pset.seek(par);
        (*item).parse(pname,pdatad,pstate);
        QCC(i, pdatad);			// Set the QCC constant
        Qtensor = 1;				// Flag to look for E.A.'s
        } 					// and eta value
//      else if(pset.seek(SinglePar(qnames[2]))) // Next try to read the QCC
      else if(pset.seek(pars[2]) != pset.end())	// Next try to read the QCC
        { 					// in kHz
        par = SinglePar(qnames[2]);
        item = pset.seek(par);
        (*item).parse(pname,pdatad,pstate);
        QCC(i, pdatad*1.e3);			// Set the QCC constant
        Qtensor = 1;				// Flag to look for E.A.'s
        } 					// and eta value
//      else if(pset.seek(SinglePar(qnames[3]))) // Next try to read QCC in MHz
      else if(pset.seek(pars[3]) != pset.end())	// Next try to read QCC in MHz
        {
        par = SinglePar(qnames[3]);
        item = pset.seek(par);
        (*item).parse(pname,pdatad,pstate);
        QCC(i, pdatad*1.e6);			// Set the QCC constant
        Qtensor = 1;				// Flag to look for E.A.'s
        } 					// and eta value
      if(Qtensor)
        {}					// Look for Euler angles
      } 					// and eta value here
    Qtensor = 0;				// Unset the flag to look
    } 						// for E.A.'s and eta value
  } 


void sys_dynamic::setRand(const ParameterSet& pset)

	// Input		dsys	: Dynamic spin system (this)
        //                      pset	: A parameter set
        // Output               none    : Spin system random field 
	//				  relaxation values are set from
	//				  parameters in pset

  {
  std::string rnames[4];
  double rflw=0, rftau=0;
  rnames[0] = "RFlw";				// Rdm Field lwhh (Hz) all spins
  rnames[1] = "RFlw(";				// Rdm Field lwhh (Hz) one spin
  rnames[2] = "RFtau";				// Rdm Field corr. time all spins
  SinglePar par(rnames[0]);			// Working single parameter
  ParameterSet::const_iterator item;		// Working pix into parameter set
  int i, ns=spins();
  std::string pname, pstate;
  if(pset.seek(par) != pset.end())
    {
    par = SinglePar(rnames[0]);
    item = pset.seek(par);
    (*item).parse(pname,rflw,pstate);
    for(i=0; i<ns; i++)				// Loop all the spins in system
      rand_As[i] = rflw;			// Assumed linewidth is in Hz
    }
  else
    for(i=0; i<ns; i++)				// Loop all the spins in system
      {
      rnames[3] = rnames[1] + Gdec(i) + ")";
      par = SinglePar(rnames[3]);
      if(pset.seek(par) != pset.end()) 		// Try to read random field lwhh
        {
        par = SinglePar(rnames[3]);
        item = pset.seek(par);
        (*item).parse(pname,rflw,pstate);
        rand_As[i] = rflw;			// Assumed that linewidth in Hz
        }
      else rand_As[i] = 0.0;
      }
  par = SinglePar(rnames[2]);
  if(pset.seek(par) != pset.end()) 		// Try read rand.f. tau (psec)
     {
     par = SinglePar(rnames[2]);
     item = pset.seek(par);
     (*item).parse(pname,rftau,pstate);
     rand_As[ns] = rftau*1.e-12;		// Store correlation time in sec
     }
  else rand_As[ns] = 0.0;			// Store correlation time in sec
  }
 

bool sys_dynamic::setKs(const ParameterSet& pset, bool warn)
  {
  MExs.clear(); 				// We don't know any processes
  ExchProcM pro;				// Single exchange processes
  int i=0;                                      // Exchange process index
  while(pro.read(pset, i, false))               // Try & read exchange process i
    {
    if(!CheckExch(pro))	DSfatal(37); 		// Insure process is valid
    MExs.push_back(pro);			// If we can, store it in procs
    i++;                                        // Increment to next process
    }
  if(!i)                                        // If we didn't read any
    { if(warn) DSerror(63,1); return false; }	// Warn and return our failure
  return true;                                  // Return we read em OK
  }

 
void sys_dynamic::setTaus(const ParameterSet& pset, int mand)

	// Input		dsys	: Dynamic spin system (this)
        //                      pset	: A parameter set
	//			mand	: Flag whether taus mandatory
        // Output               none    : Spin system correlation times
	//				  are from the parameters in pset
	// Note				: These are mandatory if mand != 0

  {
  std::string names[4] = { "Taus", "Tausp",	// Possible parameter names
                      "Ds", "DsHz" };	// (ns, ps, rad/sec, Hz)
  ParameterSet::const_iterator item;	// A pix into parameter list
  int i = 0;				// Parameter name index
  int found = 0;			// Flag if any tau parameter found
  for(i=0; i<4 && !found; i++)		// Look in pset for any tau parameter
    {
    item = pset.seek(names[i]); 	// Parameter to look for (Taus)
    if(item != pset.end()) found++;	// If item exists, we've found a tau
    }
  if(!found)				// If we couldn't find any tau
    {					// related parameters, then
    if(mand)				// 1) If mandatory, we must abort!
      {
      DSerror(9);			// 	No taus specified
      DSfatal(3);			// 	Can't construct from pset
      }
    Taus.xyz(0,0,0);			// 2) If not mandatory, set to 0
    return;				// and be done with it
    }
  i--;					// We found parameters, this index
  switch(i)				// Now we'll set the taus.  They
    {					// are all stored in sec!
    case 0:				// Tau input in nsec
      Taus = coord(*item)*1.e-9; break;
    case 1:				// Tau input in psec, keep in sec
      Taus = coord(*item)*1.e-12; break;
    case 2:				// Diffusion constants input (10**8/sec)
      Taus = coord(*item)*1.e+12;
      Taus.invert(); break;
    case 3:				// Diffusion constants input (10**8 Hz)
      Taus = coord(*item)*3.76988e+8; 	// 37.69911 = 2*pi*6.
      Taus.invert(); break;
    default: DSfatal(0);		// We should never reach here!
    }
  double tmin = taux();
  if(tauy() < tmin) tmin = tauy();
  if(tauz() < tmin) tmin = tauz();
  if(tmin <= 0)				// (Not hit if no taus set)
    {
    DSerror(10);			// Bad correlation time specified
    DSfatal(3);			// Can't construct from pset
    }
  }


void sys_dynamic::operator= (const ParameterSet& pset)

	// Input		dsys     : Dynamic spin system (this)
	// 			pset     : A parameter set
	// Output		none	 : Dynamic spin system filled with
	//				   parameters in pset
	// Point 1			 : Four things are gleaned from
	//				   the parameter set from spin_system
	//				     1.) The number of spins
	//				     2.) An isotope type for each spin
	//				     3.) An optional spin system name
	//				     4.) A spectrometer frequency
	// Point 2			 : Two things normally in spin_system
	//				   here read only if tensor info missing
	//				     5.) Chemical shifts (perhaps)
	// 				     6.) Coupling constants (perhaps)
	// Point 3			 : Data read from coord_vec pset is
	//				     7.) Spin coordinates
	// Point 4			 : Parameters read specifically for the
	//				   dynamic system are
	//				     8.) Tensor information
	//				     9.) Correlation times
	// Note				 : Functions placing  a dynamic system
	//				   into a parameter set must contain
	//				   all the information read here
	// Note				 : This function can't use corresponding
	//				   function of base class spin_system
	//				   due to the tensor treatment of shifts
	//				   and J-couplings.  It does mimic the
	//				   spin_system function to some extent.

  {
  int ns = getSpins(pset);              // Get the number of spins (spin_sys)
  *this = sys_dynamic(ns); 		// Set system for ns spins
  setIs(pset);                          // Get the isotope types (spin_sys)
  setName(pset);                        // Get the system name (spin_sys)
  setOm(pset);                          // Get spect. frequency (spin_system)
  int ccount = setCoords(pset);		// Get the spin coordinates
  if(ccount) setDip();			// Set the dipolar tensors 
  SetCSA(pset);				// Set shift & shift anisotropy
  setJs(pset); 				// Set scalar couplings (spin_system)
  setQuad(pset);			// Set quadrupolar parameters
  setRand(pset);			// Set random field parameters
  setTaus(pset);			// Set the correlation time(s) 
  setKs(pset, false);			// Set the exchange processes
  return;
  } 

// ____________________________________________________________________________
//                     DYNAMIC SYSTEM OUTPUT FUNCTIONS
// ____________________________________________________________________________


void sys_dynamic::write(const std::string &filename)

	// Input		dsys   	 : Dynamic spin system (this)
	//			filename : Output file name
	// Output		none 	 : Dynamic spin system is written as a 
	//				   parameter set to file filename

  {
  std::ofstream ofstr(filename.c_str());	// Open filename for input
  if(!ofstr.good())			// If file bad then exit
    {
    DSerror(1, filename);		// Problems with file filename
    DSfatal(5);			// Can't write to file filename
    }
  ParameterSet pset;			// Declare a parameter set
  pset += (*this);			// Add in spin system parameters
  pset.write(ofstr);			// Write pset to output stream
  }


// ____________________________________________________________________________
//                     DYNAMIC SYSTEM INPUT FUNCTIONS
// ____________________________________________________________________________



	// Input		dsys     : Dynamic spin system (this)
	// 			filename : Input filename
	//			idx	 : Parameter index vlaues used for
	//				   prefix [#] in input names
	//			warn	 : Warning output level
	//				      0 = no warnings
	//				      1 = warnings
	//				     >1 = fatal warnings 
	// Output		TF	 : Dynamic spin system filled with
	//				   parameters read from file
	//				   TRUE if successful
	// Note				 : File filename should be an ASCII
	//				   file with known sys parameters
	// Input		dsys     : Dynamic spin system (this)
	// 			pset	 : A parameter set
	//			idx	 : Parameter index vlaues used for
	//				   prefix [#] in input names
	//			warn	 : Warning output level
	//				      0 = no warnings
	//				      1 = warnings
	//				     >1 = fatal warnings 
	// Output		TF	 : Dynamic spin system filled with
	//				   parameters read from file
	//				   TRUE if successful


int sys_dynamic::read(const std::string &filename, int idx, int warn)
  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(filename, warn?1:0))	// Try and read in pset
    {					// If we cannot read the file then
    if(warn)				// we'll issue warnings as desired
      {
      DSerror(1,filename,1);		//	Problems with file
      if(warn>1) DSfatal(21);	//	This is a fatal problem
      else       DSerror(21,1);	//	Or perhaps it isn't so bad
      }
    return 0;
    }
  return read(pset, idx, warn);	// Use overload to fill up system
  }

int sys_dynamic::read(const ParameterSet& pset, int idx, int warn)
  {
  int TF = setSsys(pset, idx, warn?1:0);// Use overload to read system
  if(!TF)				// If we couldn't set system from pset
    {					// then we'll issue warnings as desired
    if(warn)
      {
      DSerror(8, 1);			//	Problems with pset
      if(warn>1) DSfatal(21);	//	This is a fatal problem
      else       DSerror(21,1);	//	Or perhaps it isn't so bad
      }
    return 0;
    }
  return TF;
  }



	// Input		dsys    : Dynamic spin system (this)
        //                      argc    : Number of arguments
        //                      argv    : Vector of argc arguments
        //                      argn    : Argument index
        // Output               void    : The parameter argn of array argc
        //                                is used to supply a filename
        //                                from which the spin system is read
        //                                If the argument argn is not in argv,
        //                                the user is asked to supply a filename
        // Note                         : The file should be an ASCII file
        //                                containing recognized sys parameters
        // Note                         : The spin system is modifed (filled)

std::string sys_dynamic::ask_read(int argc, char* argv[], int argn)
  {
  std::string filename;				// Name of spin system file  
  query_parameter(argc, argv, argn,		// Get filename from command
       "\n\tSpin system filename? ", filename); // line or ask for them
  read(filename);		           	// Read system from filename
  return filename;
  }

std::string sys_dynamic::ask_read(int argc, char* argv[], int argn, const std::string& def)
  {
  std::string msg = "\n\tSpin system filename ["     // Query we will ask if
             + def + "]? ";                     // it is needed
  std::string filename = def;                        // Name of spin system file
  ask_set(argc,argv,argn,msg,filename);         // Or ask for it
  read(filename);                               // Read system from filename
  return filename;                              // Return filename
  }

// ____________________________________________________________________________
//                 DYNAMIC SYSTEM CORRELATION TIME FUNCTIONS
// ____________________________________________________________________________

/* These functions allow access to the system correlation times.  Functions
   exist to obtain all three in a coordinate { taux, tauy, tauz } or obtain
   them individually where each ordinate is returned as a double.  The values
   must be set individually. The input/output units are always in seconds.   */

coord  sys_dynamic::taus() const     { return Taus;     }
double sys_dynamic::taux() const     { return Taus.x(); }
double sys_dynamic::tauy() const     { return Taus.y(); }
double sys_dynamic::tauz() const     { return Taus.z(); }
void   sys_dynamic::taux(double tau) { Taus.x(tau);     }
void   sys_dynamic::tauy(double tau) { Taus.y(tau);     }
void   sys_dynamic::tauz(double tau) { Taus.z(tau);     }

// ____________________________________________________________________________
//                         EXCHANGE RATE MANIPULATIONS
// ____________________________________________________________________________
 
///F_list Kex                    - Exchange rate array (1/sec)

double sys_dynamic::Kex(int p) const
  {
  if(!CheckExch(p)) { DSerror(7, Gdec(p), 1); DSfatal(35); }
  return  MExs[p].Kex();
  }

void sys_dynamic::Kex(double K, int p)
  {
  if(!CheckExch(p)) { DSerror(7, Gdec(p), 1); DSfatal(35); }
  MExs[p].Kex(K);
  }

matrix sys_dynamic::Kex() const { return Kmx; }
void   sys_dynamic::Kex_zero()  { matrix X; Kmx = X; }
 
        // Input                dsys     : Dynamic spin system (this)
        // Output               void     : All mutual exchange process
        //                                 removed.
void sys_dynamic::Kex(int i, int j, double K)
 
        // Input                dsys     : Dynamic spin system (this)
        //                      i,j      : Spin pair indices
        //                      K        : Exchange rate (1/sec)
        // Output               void     : Exchange rate between spins
        //                                 i & j is set to K
        // Note                          : This is mutual exchange, so
        //                                 i & j must be same isotopes
        // Note                          : K must be non-negative
 
  {
  if(K<0)				// Insure exchange rate non-negative
    {
    DSerror(30);			//	Exchange process problems
    DSfatal(31);			//	Negative exchange rate!
    }
  if(symbol(i) != symbol(j))		// Insure exchange is homonuclear
    {
    DSerror(30);			//	Exchange process problems
    DSfatal(32);			//	Exchange between heteronuclei
    }
  return;
  }

 
void sys_dynamic::Kex(int N, int* Is, double K)
 
        // Input                dsys     : Dynamic spin system (this)
        //                      N        : Number of spins
        //                      Is       : Array of spin indices
        //                      K        : Exchange rate (1/sec)
        // Output               void     : Exchange process between N
        //                               : spins of indices Is set to K
        // Note                          : This is mutual exchange of a
        //                                 cyclical nature -
        //                                 Is[0]<->Is[1]<->...<->I[N-1]<->Is[0]
        // Note                          : All spins must be same isotopes
        // Note                          : K must be non-negative
 
  {
// sosi - lost this function, must now rebuild.
if(Is == NULL) N =7;
K =3;
  return;
  }

const std::vector<ExchProcM>& sys_dynamic::MExProcs() const
  { return MExs; }

// ______________________________________________________________________
//                       STANDARD I/O FUNCTIONS
// ______________________________________________________________________

std::vector<std::string> sys_dynamic::PtStrings(int cw1, int cw2, int digs) const
  {
  std::vector<std::string> StrArray;			// Make a string vector
  if(spins() < 2) return StrArray;			// No strings if no spins
  if(!coord_vec::size()) return StrArray;		// No strings if no coords
  if(!max_x() && !max_y() && max_z()) return StrArray;	// No strings if no coords
  int ns = spins();					// Spins in the system
  std::string NA(cw2, '-');				// Use this if no coordinate
  std::string fmt = std::string("%") + Gdec(cw2-4)	// Set # format string
                  + std::string(".")  + Gdec(digs)
                  + std::string("f");
  std::string line0, lineX, lineY, lineZ;		// Strings to add in
  double absx, absy, absz;				// For coordinates 
  std::string un;					// For units
  double sf;						// For scaling factors
  int i;
  line0 = CenterString("Spin Coords", cw1);		// Header row 1st column
  lineX = CenterString("X", cw1) + std::string(" :");	// X row First column
  lineY = CenterString("Y", cw1) + std::string(" :");	// X row First column
  lineZ = CenterString("Z", cw1) + std::string(" :");	// X row First column

  for(i=0; i<ns; i++)					// Loop spins & add coords
    {
    absx = fabs(coord_vec::x(i));
    absy = fabs(coord_vec::y(i));
    absz = fabs(coord_vec::z(i));
    if        (absx <= 1.e-13 && absy <= 1.e-13 && absz <= 1.e-13) { sf=1.0;   un = " A  "; }
    else if   (absx >= 1.e-6  || absy >= 1.e-6  || absz >= 1.e-6)  { sf=1.e6;  un = " um "; }
    else if   (absx >= 1.e-9  || absy >= 1.e-9  || absz >= 1.e-9)  { sf=1.e9;  un = " nm "; }
    else if   (absx >= 1.e-10 || absy >= 1.e-10 || absz >= 1.e-10) { sf=1.e10; un = " A  "; }
    else                                                           { sf=1.e12; un = " pm "; }

    lineX += std::string(" ") + Gform(fmt.c_str() + un, coord_vec::x(i)*sf);
    lineY += std::string(" ") + Gform(fmt.c_str() + un, coord_vec::y(i)*sf);
    lineZ += std::string(" ") + Gform(fmt.c_str() + un, coord_vec::z(i)*sf);
    }
 
  StrArray.push_back(line0);
  StrArray.push_back(lineX);
  StrArray.push_back(lineY);
  StrArray.push_back(lineZ);

  return StrArray;					// Return the string array
  }

std::vector<std::string> sys_dynamic::AQStrings(int cw1, int cw2, int digs) const
  {
  std::vector<std::string> StrArray;			// Make a string vector
  int ns = spins();					// Spins in the system
  if(!ns) return StrArray;				// No strings if no spins
 
  if(!Quad()) return StrArray;				// No string if no quads
  std::string blanks(cw2, ' ');				// Filler for non-existing #s
  std::string szero = CenterString("0", cw2);		// Print for zero numbers
  std::string NA(cw2, '-');				// Use this if no quadrupole
  std::string fmt = std::string("%") + Gdec(cw2-4)	// Set # format string
                  + std::string(".")  + Gdec(digs)
                  + std::string("f");

//  std::string line0 = "Quad Tensors";
  std::string line1 =  std::string("Quad QCC  ") + std::string(" :");
  std::string line2 =  std::string("Quad eta  ") + std::string(" :");
  std::string line3 =  std::string("Quad alpha") + std::string(" :");
  std::string line4 =  std::string("Quad beta ") + std::string(" :");
  std::string line5 =  std::string("Quad gamma") + std::string(" :");
  std::string deg(" deg ");
  std::string ubk("     ");
  std::string Qun;
  double qcc, absqcc, scf;

  std::string b(" ");
  for(int i=0; i<ns; i++)				// Loop spins & add QA's
    {
    if((quad_As[i]).exists())
      {
      qcc = quad_As[i].delz();
      absqcc = fabs(qcc);
           if(absqcc >= 1.e6) { scf = 1.e-6; Qun = std::string(" MHz "); }
      else if(absqcc >= 1.e3) { scf = 1.e-3; Qun = std::string(" kHz "); }
      else                    { scf = 1.0;   Qun = std::string(" Hz  "); }
      line1 += b + Gform(fmt.c_str() + Qun, qcc * scf);
      line2 += b + Gform(fmt.c_str() + ubk, quad_As[i].eta());
      line3 += b + Gform(fmt.c_str() + deg, quad_As[i].alpha());
      line4 += b + Gform(fmt.c_str() + deg, quad_As[i].beta());
      line5 += b + Gform(fmt.c_str() + deg, quad_As[i].gamma());
      }
    else if(electron(i) || qn(i) < 1.0)
      {
      line1 += b + NA;
      line2 += b + NA;
      line3 += b + NA;
      line4 += b + NA;
      line5 += b + NA;
      }
    else
      {
      line1 += b + szero;
      line2 += b + szero;
      line3 += b + szero;
      line4 += b + szero;
      line5 += b + szero;
      }
    } 

//  StrArray.push_back(line0);
  StrArray.push_back(line1);
  StrArray.push_back(line2);
  StrArray.push_back(line3);
  StrArray.push_back(line4);
  StrArray.push_back(line5);
  return StrArray;					// Return the string array
  }

std::ostream& sys_dynamic::printAC(std::ostream& ostr) const

	// Input		dsys     : Dynamic spin system (this)
	// 			ostr	 : Output stream
	// Output		none	 : Dynamic spin system 
	//				   shift anisotropy spatial tensors
	//				   sent to the output stream

  {
  int ns = spins();				// Number of system spins
  int shield = 0;				// Flag if shielding tensors
  int i; 					// Temp indices, dipole index
  for(i=0; i<ns && !shield; i++) 		// Check if SA tensors exist
    if((shift_As[i]).exists())
      shield = 1;
  if(!shield) return ostr;			// Exit if no SA tensors
  ostr << "\nShift Tensors"; 			// Write out the shift tensors
  ostr << " (PPM)";  
  std::string blanks10("          ");		// Filler for non-existing #s
  std::string szero("      0   ");			// Print for zero numbers
  std::string sthdrs[6] = { "Aiso     ", "CSA      ",// Dipolar tensor component
                       "eta      ", "alpha    ",// headers
                       "beta     ", "gamma    " }; 
  double cmpval=0;				// Tensor component index
  for(int cmp=0; cmp<6; cmp++)
    {
    ostr << "\n" << sthdrs[cmp] << ":";		// Write component header
    for(i=0; i<ns; i++)				// Loop through spins
      {
      if((shift_As[i]).exists())
        {
        if(cmp==0)      cmpval=shift_As[i].iso();
        else if(cmp==1) cmpval=1.5*shift_As[i].delz();
        else if(cmp==2) cmpval=shift_As[i].eta();
        else if(cmp==3) cmpval=shift_As[i].alpha();
        else if(cmp==4) cmpval=shift_As[i].beta();
        else if(cmp==5) cmpval=shift_As[i].gamma();
        ostr << Gform("%10.2f", cmpval);
        }
      else ostr << szero;
      }
    }
  ostr << "\n";
  return ostr;
  }


std::ostream& sys_dynamic::printAD(std::ostream& ostr) const

	// Input		dsys     : Dynamic spin system (this)
	// 			ostr	 : Output stream
	// Output		none	 : Dynamic spin system 
	//				   dipolar spatial tensors
	//				   sent to the output stream

  {
  int ns = spins();				// Number of system spins
  int dipolars = 0;				// Flag if Dipolar Tensors
  int k, l, dip=0;				// Temp indices, dipole index
  for(k=0; k<ns-1 && !dipolars; k++)		// Check if dipoles tensors exist
    for(l=k+1; l<ns && !dipolars; l++, dip++)
      if((dip_As[dip]).exists())
        dipolars = 1;
  if(!dipolars) return ostr;			// Exit if no dipolar tensors
  ostr << "\nDipolar Tensors";			// Write out the dipolar
  ostr << " (delz in rad/sec)";			// tensors title 
  int diptot = 0;				// Need to count per line
  int doiso=0, doeta=0, dogam=0;		// Count each non-zero cmpt
  dip = 0;					// For dipole indexing
  for(k=0; k<ns-1; k++)				// Loop over all the dipoles
    {						// and look at the Aiso,
    for(l=k+1; l<ns; l++, dip++)		// eta, and gamma values.
      {						// Normally these are zero
      if(dip_As[dip].iso())   doiso++;		// for a dipolar tensor so
      if(dip_As[dip].eta())   doeta++;		// we don't need to print
      if(dip_As[dip].gamma()) dogam++;		// them.  But we check just
      dip ++;					// in case.....
      }
    }
  std::string blanks10("          ");		// Filler for non-existing #s
  std::string sthdrs[6] = { "Aiso     ", "delz(DCC)",// Dipolar tensor component
                       "eta      ", "alpha    ",// headers
                       "beta     ", "gamma    " }; 
  int stcmpi = 0;				// Flag to print Aiso
  if(doiso) stcmpi=1;				// If no Aiso, don't print
  int stcmps = 5;				// Number of components
  if(dogam) stcmps++;				// Adjust to print gamma
  double cmpval=0;				// Tensor component index

  for(k=0; k<ns-1; k++)				// Begin with 1st spin index
    {
    ostr <<"\nSpin "<< Gdec("%d2", k)<<"  :";
    for(int cmp=stcmpi; cmp<stcmps; cmp++)
      {
      if(cmp==2 && !doeta) break;
      dip = 0;
      ostr << "\n" << sthdrs[cmp] << ":";
      for(l=0; l<ns; l++)
        {
        if(l<=k) ostr << blanks10;
        else
          {
          if(cmp==0)      cmpval=dip_As[diptot+dip].iso();
          else if(cmp==1) cmpval=dip_As[diptot+dip].delz();
          else if(cmp==2) cmpval=dip_As[diptot+dip].eta();
          else if(cmp==3) cmpval=dip_As[diptot+dip].alpha();
          else if(cmp==4) cmpval=dip_As[diptot+dip].beta();
          else if(cmp==5) cmpval=dip_As[diptot+dip].gamma();
          ostr << Gform("%10.2f", cmpval);
          dip ++;
          }
        }
      }
    diptot += dip;			// Update dipole index
    ostr << "\n";
    }
  return ostr;
  }


std::ostream& sys_dynamic::printAQ(std::ostream& ostr) const

	// Input		dsys     : Dynamic spin system (this)
	// 			ostr	 : Output stream
	// Output		none	 : Dynamic spin system 
	//				   quadrupolar spatial tensors
	//				   sent to the output stream

  {
  int ns = spins();				// Number of system spins
  int quad = 0;					// Flag if quad. tensors
  int i; 					// Temp index
  for(i=0; i<ns && !quad; i++) 			// Check if quad tensors exist
    if((quad_As[i]).exists())
      quad = 1;
  if(!quad) return ostr;			// Exit if no quad tensors
  ostr << "\nQuadrupolar Tensors (delz in KHz)";// Write out the Quadrupolar tensors
  std::string blanks10("          ");		// Filler for non-existing #s
  std::string szero("      0   ");			// Print for zero numbers
  std::string sthdrs[5] = { "delz(QCC)",		// Quadrupolar tensor component
                       "eta      ", "alpha    ",// headers
                       "beta     ", "gamma    " }; 
  double cmpval=0;				// Tensor component index
  for(int cmp=0; cmp<5; cmp++)
    {
    ostr << "\n" << sthdrs[cmp] << ":";		// Write component header
    for(i=0; i<ns; i++)				// Loop through spins
      {
      if((quad_As[i]).exists())
        {
        if(cmp==0)      cmpval = (quad_As[i].delz()) * 1.e-3;
        else if(cmp==1) cmpval=quad_As[i].eta();
        else if(cmp==2) cmpval=quad_As[i].alpha();
        else if(cmp==3) cmpval=quad_As[i].beta();
        else if(cmp==4) cmpval=quad_As[i].gamma();
        ostr << Gform("%10.2f", cmpval);
        }
      else ostr << szero;
      }
    }
  ostr << "\n";
  return ostr;
  }


std::ostream& sys_dynamic::printARDM(std::ostream& ostr) const

	// Input		dsys     : Dynamic spin system (this)
	// 			ostr	 : Output stream
	// Output		none	 : Dynamic spin system 
	//				   random field spatial tensors
	//				   sent to the output stream

  {
  int ns = spins();				// Number of system spins
  int rand = 0;					// Flag if random field tensors
  int i; 					// Temp index
  for(i=0; i<ns && !rand; i++) 			// Check rand.fld. tensors exist
    if(rand_As[i]) rand = 1;
  if(!rand) return ostr;			// Exit if no rand.fld. tensors 
  ostr << "\nRandom Field Parameters";		// Write out Random Field stuff
  ostr << "\nRFlwhh   :";			// Write rand.fld. SQT linewidths
  for(i=0; i<ns; i++)
    {
    if(rand_As[i]) ostr << Gform("%10.2f", rand_As[i]);
    else           ostr << "      0   ";
    }
  ostr << "\nRFtau    :"; 			// Effective correlation time
  if(rand_As[ns])
    ostr << rand_As[ns]*1.e12 << " picoseconds";
  else ostr << " Extreme Narrowing";
  ostr << "\n";
  return ostr;
  }



	// Input		dsys     : Dynamic spin system (this)
	// 			ostr	 : Output stream
	// Output		none	 : Dynamic spin system 
	//				   rotational correlation times
	// Note				: Taus are internally stored
	//				  in seconds!

std::ostream& sys_dynamic::printTaus(std::ostream& ostr) const
  {
  if(!Taus.x() && !Taus.y() && !Taus.z())	// Exit if no taus defined
    return ostr;
  double taumin = Taus.x();			// Find smallest tau value
  if(taumin > Taus.y()) taumin = Taus.y();	// in order to set units
  if(taumin > Taus.z()) taumin = Taus.z();	// for reasonable printing

  double taufact=1.0;				// Print scaling factor
  std::string U("sec");				// Print units
  if(taumin < 1.e-12)				// See if print in femto
    {						// If so, must scale by
    taufact = 1.e15;				// 10^15 and use fsec units
    U = " fs";
    }
  else if(taumin < 1.e-9)			// See if print in pico
    {						// If so, must scale by
    taufact = 1.e12;				// 10^12 and use psec units
    U = " ps";
    }
  else if(taumin < 1.e-6)			// See if print in nano
    {						// If so, must scale by
    taufact = 1.e9;				// 10^9 and use nsec units
    U = " ns";
    }
  else if(taumin < 1.e-3)			// See if print in micro
    {						// If so, must scale by
    taufact = 1.e6;				// 10^6 and use usec units
    U = " us";
    }
  else if(taumin < 1)				// See if print in milli
    {						// If so, must scale by
    taufact = 1.e3;				// 10^3 and use msec units
    U = " ms";
    }

  std::string T[3] = { "x", "y", "z" };		// Tau axis labels
  ostr << "\nCorrelation Times";		// Write out header
  for(int i=0; i<3; i++)			// Print all three correlation
    {						// times
    ostr << "\ntau(" << T[i] << ")   :"
         << Gform("%10.4f", Taus.get(i)*taufact)
         << U;
//        << setw(10) << setprecision(4)
//        << Taus.get(i) * taufact << U;

    }
  ostr << "\n";
  return ostr;
  }

// ----------------------------------------------------------------------
//          Print Non-Mutual Exchange Processes (If Any Defined)
// ----------------------------------------------------------------------

std::ostream& sys_dynamic::printEX(std::ostream& ostr) const
  {
  int np = MExs.size();				// # of exchange processes
  if(!np) return ostr;				// Exist if no exchange
  std::string hdr = "Mutual Exchange Processes";	// Header
  int hl = hdr.length();			// Header length
  std::string spacer = std::string(40-hl/2, ' ');		// Print header centered
  ostr << "\n" << spacer << hdr << "\n";
  int i, len=0;					// Maximum string length
  std::vector<std::string> xps;				// Exchange process strings
  for(i=0; i<np; i++)				// Loop exchange processes
    {
    xps.push_back(MExs[i].ExchStr());		//   Store process string
    len = gmax(len, int(xps[i].length()));	//   Track longest string
    }
  std::string RS("     Rate    ");		// Title for rate column
  std::string RH("-------------");		// U.L. for rate column
  std::string RU("/sec");			// Rate output units
  std::string EXSS("Exchanging Spins");		// Title for process column
  std::string EXSH("----------------");		// U.L. for process column
  std::string SB("    ");			// Spacer between columns
  if(len>16)					// If process strings are
    {						// longer than colunn title
    EXSS = std::string((len-16)/2, ' ') + EXSS;	// adjust column width to
    EXSH = std::string((len-16)/2, '-') + EXSH;	// match
    EXSS += std::string(len-EXSS.length(), ' ');
    EXSH += std::string(len-EXSH.length(), '-');
    }
  std::string EXSP("");				// If process strings short
  if(len<16) EXSP = std::string((16-len)/2, ' ');// buffer to center in column

  hdr = RS + SB + EXSS;				// This is header, all columns
  hl = hdr.length();				// Length of the header
  spacer = std::string(40-hl/2, ' ');		// To center the header
  ostr << "\n" << spacer << hdr;		// Print header centered
  ostr << "\n" << spacer << RH << SB << EXSH;	// Print U.L.s centered
    
  double K;					// An exchange process
  for(i=0; i<np; i++)				// Loop all pricesses
    {
    K = MExs[i].Kex();				//   Get exchange rate
    ostr << "\n" << spacer << Gform("%9.3f", K)	//   Output line for this
         << RU << SB << EXSP << xps[i];		//   exchange process
    }
  ostr << "\n";					// Added line for space
  return ostr;					// We are done
  }


// ----------------------------------------------------------------------
//                      Print Entire Spin System
// ----------------------------------------------------------------------

std::ostream& sys_dynamic::print(std::ostream& ostr) const
  {
  spin_system::print(ostr);		// Write spin_system parameters
  if(!spins()) return ostr;		// Exit now if no spins
  if(coord_vec::size() && spins()>1)	// If coordinates present print
    {					// them if multiple spins
    coord_vec::printf(ostr, 1);
    ostr << "\n";
    }
  printAD(ostr);			// Print dipolar tensors
  printAC(ostr); 			// Print shift anisotropy tensors
  printAQ(ostr);			// Print quadurpolar tensors
  printARDM(ostr);			// Print random field tensors
  printTaus(ostr); 			// Print correlation times
  printEX(ostr); 			// Print mutual exchange processes
  return ostr;
  }

std::ostream& operator<< (std::ostream& out, const sys_dynamic& dsys)
  { return dsys.print(out); }


	// Input		dsys     : Dynamic spin system (this)
	// 			out      : Output stream
	//			full	 : A print flag
	// Output		none	 : Dynamic spin system dipoles 
	//				   sent to the output stream

std::ostream& sys_dynamic::print_D(std::ostream& out, int full)
  {
  coord_vec dipv;				// Coordinate vector for info
  int ns = spins();				// Number of spins in the system
  if(full)					// If full flag is set, then
    {
    full = 1;					//	Make sure full is 1
    dipv = vectors_f();				//	Fill dipv with values for
    }						//	all spin pairs (i-j & j-i)
  else						// If full flag not set, then
    dipv = vectors();				//	Values in dipv for (i->j, not j->i)

  out << "\nDipoles";				// Begin outputting the info
  int kr=0, kt=0, kp=0;
  int i, j;
  for(i=0; i<ns-1+full; i++)
    {
    out << "\nRadius " << Gdec("%d2", i) << ":";//	Output distance between spins
    for(j=0; j<ns; j++)				//	This loop does spin i with all
      if (full || i<j)				//	j spins greater than index i
        {					//      if full flag is not set or with all
	out << Gform("%10.2f", dipv.x(kr));	//	j spins if full flag has been set
        kr++;
        }
      else
	out << "          ";

    out << "\nTheta  " << Gdec("%d2", i) << ":";//	Output spherical angle theta (degrees)
    for(j=0; j<ns; j++)				//	Again, loop does spin i will spins j>i
    if(full || i<j)				//	unless full flag is set then i with all
      {						//	spins j is output
      out << Gform("%10.2f", dipv.y(kt)*180.0);
      kt++;
      }
    else
      out << "          ";

    out << "\nPhi    " << Gdec("%d2", i) << ":";//	Output the spherical angle phi (degrees)
    for(j=0; j<ns; j++)				//	Loop does spin i with spins j>i unless
      if(full || i<j)				//	full flag is set, the does spin i with
        {					//	all spins j
	out << Gform("%10.2f", dipv.z(kp)*180.0);
        kp++;
        }
      else
	out << "          ";
    }
  out << "\n";
  return out;
  }

#endif						// sys_dynamic.cc
