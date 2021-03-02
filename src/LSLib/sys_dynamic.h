/* sys_dynamic.h ************************************************-*-c++-*-
**									**
** 	                          G A M M A				**
**								 	**
**	Dynamic Spin System                         Interface		** 
**								 	**
**	Copyright (c) 1992, 1993				 	**
**	Scott Smith 							**
**	Eidgenossische Technische Hochschule			 	**
**	Labor fuer physikalische Chemie				 	**
**	8092 Zuerich / Switzerland				 	**
**								 	**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**                                                                      **
** The class sys_dynamic defines a collection nuclear spins with        **
** specific environments and motional properties.                       **
**                                                                      **
*************************************************************************/

///Chapter Class Dynamic Spin System (sys_dynamic)
///Section Overview
///Body    The class 
///Section Available Dynamic Spin System Functions

#ifndef   Sys_dynamic_h_		// Is this file already included?
#  define Sys_dynamic_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Include system specifics
#include <HSLib/SpinSystem.h>		// Base isotropic spin system 
#include <Level1/coord.h>		// Know about coordinates
#include <Level1/coord_vec.h>		// Know coordinate vectors
#include <Basics/ParamSet.h>		// Know parameter sets
#include <Matrix/matrix.h>		// Know matrices 
#include <Level1/SpaceT.h>		// Know spatial tensors
#include <Matrix/row_vector.h>		// Know row vectors
#include <Level1/ExProcessM.h>		// Know mutual exchange processes
#include <vector>			// Know libstdc++ STL vectors

class sys_dynamic: public spin_system, public coord_vec 
  {
  std::vector<space_T> shift_As;	// Array of shift tensors
  std::vector<space_T>   dip_As;	// Array of dipolar coupling tensors
  std::vector<space_T>  quad_As;	// Array of quad. coupling tensors
  std::vector<double>   rand_As;	// Array of random relaxation values
  std::vector<ExchProcM>   MExs;	// Mutual exchange processes
  coord Taus;				// Correlation times
  matrix Kmx;				// Array of mutual exchange rates
  
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
//                  CLASS DYNAMIC SPIN SYSTEM ERROR HANDLING
// ____________________________________________________________________________

        // Input                dsys    : Dynamic spin system (this)
        //                      eidx    : Error index
        //                      pname   : String for error message
	//			noret	: Flag for return (0=return)
        // Output               none    : Error message output 
        //                                Program execution stopped if fatal

         void DSerror(int eidx,                           int noret=0) const;
         void DSerror(int eidx, const std::string& pname, int noret=0) const;
volatile void DSfatal(int eidx)                                        const;
volatile void DSfatal(int eidx, const std::string& pname)              const;

// ____________________________________________________________________________
// ii              CLASS DYNAMIC SPIN SYSTEM SETUP FUNCTIONS
// ____________________________________________________________________________


virtual int setSsys(const ParameterSet& pset,int idx=-1,int wrn=2);

        // Input                dsys    : Dynamic spin system (this)
	// 			pset     : A parameter set
	//			idx	 : Parameter index value used for
	//				   prefix [#] in input names
	//			wrn	 : Warning output level
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


// ____________________________________________________________________________
// iii             CLASS DYNAMIC SPIN SYSTEM CHECKING FUNCTIONS
// ____________________________________________________________________________

bool CheckExch(int p,         bool warn=true) const;
bool CheckExch(ExchProcM& XP, bool warn=true) const;
     
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A             DYNAMIC SPIN SYSTEM CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

///Center Basic Functions
///F_list =		      - Assignment

MSVCDLC              sys_dynamic();
MSVCDLC              sys_dynamic(int spins);
MSVCDLC              sys_dynamic(const sys_dynamic &dsys);
MSVCDLL sys_dynamic& operator= (const  sys_dynamic &dsys);
MSVCDLC virtual      ~sys_dynamic ();

// ______________________________________________________________________
// B                   CHEMICAL SHIFT MANIPULATIONS
// ______________________________________________________________________

///Center Spin System Chemical Shift Functions 
///F_list shifts		- Set all chemical shifts
///F_list shift			- Set/Retrieve a specific chemical shift
///F_list offsetShifts		- Offset all chemical shifts
///F_list PPM			- Set/Retrieve chemical shift in PPM
///F_list delz			- Shift tensor delzz value

// ************** Chemical Shift Manipulations in Hertz ***************

	// Input		shift	: Chemical shift (Hz)
	//			spin	: Spin of spin system [0, nspins-1]
        //	                offset	: Chemical shift offset value (Hz)
	//			iso	: Isotope type
	// Output (shifts)	none	: Sets all chemical shifts to shift
	// Output (shift)	none	: Set/Get shift of spin (Hz)
	// Output (offsetShifts)none	: All spins of same isotope type as spin
	//				  have their shifts modified by offset
	// Output (offsetShifts)none	: All spins of same isotope type as iso
	//				  have their shifts modified by offset

MSVCDLL void   shifts(double shift=0);
MSVCDLL void   shift(double shift, int spin);
MSVCDLL void   shift(int spin, double shift);
MSVCDLL double shift (int spin) const;
MSVCDLL void   offsetShifts(double offset, int spin); 
MSVCDLL void   offsetShifts(double offset, const std::string& iso); 

// *************** Chemical Shift Manipulations in PPM ****************

	// Input		int    : spin of spin system [0, nspins-1]
	// 			double : chemical shift of spin (PPM)
	// Output		none   : default chemical shift of spin (PPM)

MSVCDLL void PPM (int, double);

// ***************** Chemical Shift Tensor Functions ******************

	// Input		dsys     : Dynamic spin system (this)
	// 			spin     : Spin index [0, nspins-1]
	//			delzz    : value of delzz in PPM
        //                      ceta     : CSA asymmetry of spin [0,1]
	//			A        : A spatial tensor
	// Output (delz)	none     : delz component of chemical shift
	//				   spatial tensor of spin (PPM)
	// Output (delz)	none     : delz component of chemical shift
	//				   spatial tensor of spin set to
	//				   input value of delzz
        // Output (Ceta)        none     : Return asymmetry of CSA tensor 
        // Output (Ceta)        none     : Set shift anisotropy asymmetry of spin 
	// Output (Ceta)	sphT	 : Spatial tensor for spin i shift
	// Output (TC)		void	 : The spatial chemical shift
	//				   tensor of spin i is set to A
	// Output (TC)		A	 : Returns chemcial shift A of spin
	// Output (xiC_vector)	dximx    : A matrix of CSA interaction
	//				   constants (rad/sec)
	// Output (xiC)		xii      : The gamma CSA interaction
	//				   constant (rad/sec)
        // Output (CSA)         TF       : Return true if any shift tensors
        //                                 are present in the system
	// Note			         : The PPM units of delzz are
	//			   	   counteracted by the MHz units
	//				   on Omega
	//		 1/2
	//   CSA   [6*pi]
	// xi    = |----| * gamma * B * del  (i) = K * Omega  * del  (i)
	//   i     [ 5  ]        i   0     zz               i      zz

MSVCDLL double     delz(int spin) const;
MSVCDLL void       delz(int spin, double delzz);
MSVCDLL double     Ceta(int spin) const; 
MSVCDLL void       Ceta(int spin, double ceta); 
MSVCDLL space_T    TC(int   spin) const;
MSVCDLL void       TC(const space_T& A, int spin);
MSVCDLL row_vector xiC_vector()   const; 
MSVCDLL double     xiC(int spin)  const;
MSVCDLL bool       CSA()          const;

// ______________________________________________________________________
//                         DIPOLAR MANIPULATIONS
// ______________________________________________________________________

///Center Dipole Related Functions

// ----------------- Spin Coordinate Access Functions -------------------
 
        // Input                cvec   : A vector of coordinates
	//			cutoff : Dipolar cutoff distance
        // Output               none   : Sets the spin coordinates
        //                               and generates dipolar tensors
        //                               for all spin pairs
        // Note                        : If the distance between two
        //                               spins is > cutoff then the dipolar
        //                               tensor is left zero.

MSVCDLL void coords(const coord_vec& cvec, double cutoff=5.e-10);

        // Output (Coord)	TF	: Return true if any spin coordinates
        //				  are present in the system

MSVCDLL bool Coord( ) const;
 
// ------------------ Dipolar Tensor Access Functions -------------------

	///F_list AD                     - Dipolar spatial tensor
	///F_list dipoles	         - Number of dipoles in system
	///F_list xiD_matrix	         - Matrix of dipolar interaction constants

	// Input		dsys	: A dynamic system (this)
	//			spin1	: Spin of spin system [0, nspins-1]
	// 			spin2	: Spin of spin system [0, nspins-1]
	//			delzz	: Value of dipolar delzz in Hz
	// 			nu	: Dipolar coupling of spin pair (Hertz)
	//      		dip	: Dipole index in spin system
	// Output (DCC)		none	: Set dipolar coupling of spin pair
	// Output (DCC)		DCC	: Return dipolar coupling constant
	//				  of the specified spin pair
	// Output (Ddelz)	none	: Return delz component of dipolar
	//				  coupling spatial tensor of spin pair (Hz)
	// Output (Ddelz)	none	: Set delz component of dipolar
	//			   	  coupling spatial tensor of spin pair
	// Output (Deta)	none	: Return asymmetry of dipolar
	//				  coupling spatial tensor of spin pair
	// 			eta	: Dipolar asymmetry of spin pair
	// Output (Deta)	none	: Set dipolar asymmetry of spin pair
	// Output (AD)		sphT	: Get spatial dipolar spatial tensor
	//				  between spin1 and spin2
	// Output (AD)		sphT	: Set spatial dipolar spatial tensor
	//				  for dipole dip
	// Output (dipoles)	int	: Number of dipoles (spin pairs)
	//				  in the dynamic spin system
	// Output (dipole)	int	: Dipole index for the spin pair
	// Output (xiD_matrix)	dximx	: A matrix of dipolar interaction
	//				  constants (unit rad/sec)
        // Output (DIP)		TF	: Return true if any dipolar tensors
        //                                are present in the system
	// Note				: DCC is defined to be the same as the
	//			          dipolar tensor delzz value
	// Note				: Normally there is NO asymmetry in a
	//				  dipolar interaction
	//					    -1      -2
	// Note				: 1T = 1 J-C  -sec-m	
	//
	//		     1/2
	//	            [6*pi]   mu 
	//	       -2 * |----| * --- * hbar * gamma  * gamma
	//	  D         [ 5  ]   4pi               i        j
	//	xi   = __________________________________________
	//	  ij		           3
	//			          r
	//			           ij


MSVCDLL double  DCC(int   spin1, int spin2) const;
MSVCDLL void    DCC(int   spin1, int spin2, double nu);
MSVCDLL double  Ddelz(int spin1, int spin2) const;
MSVCDLL void    Ddelz(int spin1, int spin2, double delzz);
MSVCDLL double  Deta(int  spin1, int spin2) const;
MSVCDLL void    Deta(int  spin1, int spin2, double Deta);
MSVCDLL space_T AD(int    spin1, int spin2) const;
MSVCDLL space_T AD(int dip) const;
MSVCDLL int     dipoles() const;
MSVCDLL int     dipole(int spin1, int spin2) const;
MSVCDLL matrix  xiD_matrix() const;
MSVCDLL bool    Dip( ) const;

// ____________________________________________________________________________
//                         QUADRUPOLAR MANIPULATIONS
// ____________________________________________________________________________

	// Input		dsys	: A dynamic system (this)
	//			spin	: Spin of spin system [0, nspins-1]
	// 			nu	: Quadrupolar coupling of spin (Hertz)
	//			delzz	: Value of quadrupolar delzz in Hz
	// 			Qeta    : Quadrupolar asymmetry of spin
        //                      A	: A quadrupolar spatial tensor
	// Output (QCC)		none	: Set quadrupolar coupling of spin
	// Output (QCC)		QCC	: Return quadrupolar coupling constant
	//				  of the specified spin
	// Output (Qdelz)	none	: Return delz component of quadrupolar
	//				  coupling spatial tensor of spin (Hz)
	// Output (Qdelz)	none	: Set delz component of quadrupolar
	//			   	  coupling spatial tensor of spin
	// Output (Qeta)	none	: Return asymmetry of quadrupolar
	//				  coupling spatial tensor of spin
	// Output (Qeta)	none	: Set quadrupolar asymmetry of spin [0,1]
	// Output (TQ)		sphT	: Spin quadrupolar spatial tensor
        // Output (TQ)          void	: The spatial quadrupolar
        //                                tensor of spin i is set to A
	// Output (xiQ_vector)	xii	: Quadrupolar interaction constants (rad/sec)
	// Output (xiQ)		xii	: Quadrupolar interaction constant  (rad/sec)
        // Output (Quad)	TF	: Return true if any quadrupolar tensors
        //                                are present in the system
	// Note				: QCC is defined to be the same as the
	//			          quadrupolar tensor delzz value
	// Note				: delz is defined to be the same as the
	//				  quadrupolar coupling constant
	//
	//		 1/2   QCC   	        1/2  del  (i)
	//   Q     [6*pi]         i       [6*pi]        zz
	// xi    = |----| * ----------  = |----|  * ----------
	//   i     [ 5  ]   2I (2I -1)    [ 5  ]    2I (2I -1)
	//                    i   i                   i   i

MSVCDLL double     QCC(int   spin) const;
MSVCDLL void       QCC(int   spin, double nu);
MSVCDLL double     Qdelz(int spin) const;
MSVCDLL void       Qdelz(int spin, double delzz);
MSVCDLL double     Qeta(int  spin) const;
MSVCDLL void       Qeta(int  spin, double Qeta);
MSVCDLL space_T    TQ(int    spin) const;
MSVCDLL void       TQ(const space_T& A, int spin);
MSVCDLL row_vector xiQ_vector( );
MSVCDLL double     xiQ(int spin);
MSVCDLL bool       Quad() const;

// ____________________________________________________________________________
//                       RANDOM FIELD MANIPULATIONS
// ____________________________________________________________________________

	// Input		dsys      : Dynamic spin system (this)
	// 			spin      : Spin index
	// Output (TR)		LWhh	  : Spin random field "spatial tensor"
	//				    is the SQT linewidth at half-height
	// Output (tauR)	tauR	  : Effective field correlation time
	// Output (xiR_vector)  xiR_vector: Rand. Fld. interact. const. (rad/sec)
	// Return		dximx     : A matrix of random field interaction
	//				    constants (rad/sec)
	// Note			          : In the random field treatment, the
	//				    xi value is defined in terms of single
	//				    quantum transition linewidths
	//
	//		                  1/2
	//   R              [    R       ]
	// xi    = 2 * pi * |LWhh / J(w )| 
	//   i              [    i     i ]

MSVCDLL double     TR(int spin)  const;
MSVCDLL double     tauR()        const;
MSVCDLL row_vector xiR_vector( ) const;
MSVCDLL double     xiR(int spin) const;

// ____________________________________________________________________________
//                 SCALAR COUPLING CONSTANT MANIPULATIONS
// ____________________________________________________________________________


// All Handled by class spin_system

// ____________________________________________________________________________
//                           PARAMETER SET FUNCTIONS
// ____________________________________________________________________________
///Center Parameter Set Functions

// ----------------------------------------------------------------------------
//        Functions To Make A Parameter Set From A Dynamic Spin System
// ----------------------------------------------------------------------------

	///F_list =		       - Conversion
	///F_list +=		       - Unary Addition

	// Input		dsys   : Dynamic spin system
	//  			pset   : Parameter set
	// Output ParameterSet	pset   : Parameter set with
	//			         only sys_dynamic parameters
	// Output +=		pset   : Parameter set with
	//			         only system parameters

MSVCDLL             operator ParameterSet( );
MSVCDLL friend void operator+= (ParameterSet& pset, sys_dynamic &dsys);

// ----------------------------------------------------------------------------
//        Functions To Make A Dynamic Spin System From A Parameter Set
// ----------------------------------------------------------------------------


MSVCDLL int setCoords(const ParameterSet& pset, int mand=0);

        // Input                dsys    : Dynamic spin system (this)
        //                      pset    : A parameter set
        //                      mand    : Flag whether coords mandatory
        // Output               T/F	: True if spin coordiantes are
	//				  found in the pset.  The system
	//				  spin coordinates are set from
	//				  the parameters in pset
        // Note                         : These be mandatory if mand!=0


MSVCDLL void setDip( );

        // Input                dsys    : Dynamic spin system (this)
        // Output               none    : Spin system dipolar
        //                                relaxation values are set from
        //                                coordinates in dsys

 
MSVCDLL void SetCSA(const ParameterSet& pset);

        // Input                dsys    : Dynamic spin system (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system shift & shift anisotropy
        //                                relaxation values are set from
        //                                parameters in pset

 
MSVCDLL void setQuad(const ParameterSet& pset);

        // Input                dsys    : Dynamic spin system (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system quadrupolar
        //                                relaxation values are set from
        //                                parameters in pset


 
MSVCDLL void setRand(const ParameterSet& pset);

        // Input                dsys    : Dynamic spin system (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system random field 
        //                                relaxation values are set from
        //                                parameters in pset


MSVCDLL void setTaus(const ParameterSet& pset, int mand=0);

        // Input                dsys    : Dynamic spin system (this)
        //                      pset    : A parameter set
        //                      mand    : Flag whether taus mandatory
        // Output               none    : Spin system correlation times
        //                                are from the parameters in pset
        // Note                         : These are mandatory if mand != 0


MSVCDLL bool setKs(const ParameterSet& pset, bool warn=true);

        // Input                dsys    : Dynamic spin system (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system exchange processes
        //                                are from the parameters in pset


MSVCDLL void operator= (const ParameterSet& pset);

	// Input		dsys   : Dynamic spin system (this)
	// 			pset   : A parameter set
	// Output		none   : Dynamic spin system filled with
	//				 parameters from pset
	// Note			       : Three things are gleaned from
	//				 the parameter set for a sys_dynamic
	//				 1.) The number of atoms
	//				 2.) The atomic weight for each atom
	//				 3.) An optional sys_dynamic name
	// Note			       : Functions which place a sys_dynamic
	//				 into a parameter set must contain
	//				 all the information read here


MSVCDLL virtual void write(const std::string &filename);

	// Input		dsys   	 : Dynamic spin system (base)
	//			filename : Output file name
	// Output		none 	 : Dynamic spin system is written as a 
	//				   parameter set to file filename
	///F_list write		         - Write system to disk file



// ____________________________________________________________________________
//                     DYNAMIC SYSTEM INPUT FUNCTIONS
// ____________________________________________________________________________

	// Input		dsys     : Dynamic spin system (this)
	// 			filename : Input filename
	// 			pset	 : A parameter set
	//			idx	 : Parameter index values used for
	//				   prefix [#] in input names
	//			warn	 : Warning output level
	//				      0 = no warnings
	//				      1 = warnings
	//				     >1 = fatal warnings 
	// Output		TF	 : Dynamic spin system filled with
	//				   parameters read from file
	//				   or int parameter set pset
	//				   TRUE if successful
	// Note				 : File filename should be an ASCII
	//				   file with known sys parameters

MSVCDLL virtual int read(const std::string &filename, int idx=-1, int warn=2);
MSVCDLL virtual int read(const ParameterSet& pset,    int idx=-1, int warn=2);

MSVCDLL virtual std::string ask_read(int argc, char* argv[], int argn);
MSVCDLL virtual std::string ask_read(int argc, char* argv[], int argn, const std::string& def);

        // Input                dsys    : Dynamic spin system (this)
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

// ____________________________________________________________________________
//                 DYNAMIC SYSTEM CORRELATION TIME FUNCTIONS
// ____________________________________________________________________________

///Center Correlation Time Functions
///F_list taus		         - Correlation times in seconds
///F_list tauy		         - Y-axis correlation time in seconds
///F_list tauz		         - Z-axis correlation time in seconds
///F_list taux		         - X-axis correlation time in seconds

	// Output		TAUS 	 : Coordinate containing the
	//				   three correlation times (sec)
	// Output		tau 	 : X-axis correlation time (sec)
	//			tau      : X-axis correlation time (sec)
	// Output		tau 	 : Z-axis correlation time (sec)

	// Output			: X-axis correlation time set
	// Output			: Y-axis correlation time set
	// Output			: Z-axis correlation time set

MSVCDLL coord  taus() const;
MSVCDLL double taux() const;
MSVCDLL double tauy() const;
MSVCDLL double tauz() const;

MSVCDLL void taux(double tau);
MSVCDLL void tauy(double tau);
MSVCDLL void tauz(double tau);

// ____________________________________________________________________________
//                          EXCHANGE RATE MANIPULATIONS
// ____________________________________________________________________________

///Center Spin System Exchange Rate Functions 
///F_list Kex		         - Exchange rate array (1/sec)


	// Input		dsys     : Dynamic spin system (this)
        //                      i,j      : Spin pair indices
        //                      K        : Exchange rate (1/sec)
        //                      N        : Number of spins
        //                      Is       : Array of spin indices
        //                      K        : Exchange rate (1/sec)
	// Output		Kex      : Array containing exchange
	//				   rate information
        // Output               void     : All mutual exchange process
        //                                 removed.
        // Output               void     : Exchange rate between spins
        //                                 i & j is set to K
        // Note                          : This is mutual exchange, so
        //                                 i & j must be same isotopes
        // Note                          : K must be non-negative
        // Output               void     : Exchange process between N
        //                               : spins of indices Is set to K
        // Note                          : This is mutual exchange of a
        //                                 cyclical nature - 
        //                                 Is[0]<->Is[1]<->...<->I[N-1]<->Is[0]
        // Note                          : All spins must be same isotopes
        // Note                          : K must be non-negative

MSVCDLL double Kex(int p) const;
MSVCDLL void   Kex(double K, int p);
MSVCDLL matrix Kex() const;
MSVCDLL void   Kex_zero();
MSVCDLL void   Kex(int i, int j, double K);
MSVCDLL void   Kex(int N, int* Is, double K);

MSVCDLL const std::vector<ExchProcM>& MExProcs() const;


// ____________________________________________________________________________
//                            STANDARD I/O FUNCTIONS
// ____________________________________________________________________________

///Center Dynamic Spin System I/O Functions

MSVCDLL std::vector<std::string> PtStrings(int w1=10, int w2=12, int digs=2) const;
MSVCDLL std::vector<std::string> AQStrings(int w1=10, int w2=12, int digs=2) const;

MSVCDLL std::ostream& printAC(std::ostream& ostr) const;

        // Input                dsys     : Dynamic spin system (this)
        //                      ostr     : Output stream
        // Output               none     : Dynamic spin system 
        //                                 shift anisotropy spatial tensors
        //                                 sent to the output stream


MSVCDLL std::ostream& printAD(std::ostream& ostr) const;
 
        // Input                dsys     : Dynamic spin system (this)
        //                      ostr     : Output stream
        // Output               none     : Dynamic spin system
        //                                 dipolar spatial tensors
        //                                 sent to the output stream


MSVCDLL std::ostream& printAQ(std::ostream& ostr) const;

        // Input                dsys     : Dynamic spin system (this)
        //                      ostr     : Output stream
        // Output               none     : Dynamic spin system
        //                                 quadrupolar spatial tensors
        //                                 sent to the output stream


MSVCDLL std::ostream& printARDM(std::ostream& ostr) const;

        // Input                dsys     : Dynamic spin system (this)
        //                      ostr     : Output stream
        // Output               none     : Dynamic spin system 
        //                                 random field spatial tensors
        //                                 sent to the output stream


MSVCDLL std::ostream& printTaus(std::ostream& ostr) const;

        // Input                dsys     : Dynamic spin system (this)
        //                      ostr     : Output stream
        // Output               none     : Dynamic spin system 
        //                                 rotational correlation times


MSVCDLL std::ostream& printEX(std::ostream& ostr) const;

        // Input                dsys     : Dynamic spin system (this)
        //                      ostr     : Output stream
        // Output               none     : Dynamic spin system
        //                                 exchange processes



	// Input		dsys     : Dynamic spin system (this)
	// 			out      : Output stream
	// Output		none	 : Dynamic spin system parameters set to
	//				   output stream
	///F_list print		         - Write system to output stream

// ----------------------------------------------------------------------
//                      Print Entire Spin System
// ----------------------------------------------------------------------

        // Input                out      : Output stream;
        //                      dsys     : Dynamic spin system to write
        // Output               	 : Modifies output stream
	///F_list << 		         - Standard Output

MSVCDLL virtual std::ostream& print(std::ostream& out) const;
MSVCDLL friend  std::ostream& operator<< (std::ostream& out, const sys_dynamic& dsys);

MSVCDLL std::ostream& print_D(std::ostream& out, int full=0);

	// Input		dsys     : Dynamic spin system (this)
	// 			out      : Output stream
	//			full	 : print flag
	// Output		none	 : Dynamic spin system dipoles 
	//				   sent to the output stream
	///F_list print		         - Write dipoles to output stream
};

#endif						// sys_dynamic.h
