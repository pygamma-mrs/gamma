/* SpinSystem.h *************************************************-*-c++-*-
**									**
**	                        G A M M A				**
**								 	**
**	Spin System                                Interface		**
**						 			**
**	Copyright (c) 1990, 1991, 1992				 	**
**	Scott Smith						 	**
**	Eidgenoessische Technische Hochschule			 	**
**	Labor fuer physikalische Chemie				 	**
**	8092 Zuerich / Switzerland				 	**
**						 			**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**								 	**
** The class spin system defines the number of spins, their chemical	**
** shifts, scalar coupling constants, isotope types, and input/output	**
** routines.  In addition, there are spin flags, a field strength, and	**
** spin pair flags.  A default isotope type and a default spectrometer	**
** frequency is also maintained herein.					**
**									**
*************************************************************************/

///Chapter Class Spin System
///Section Overview
///Body    None
///Section Available Spin System Functions

#ifndef   Spin_system_h_			// Is file already included?
#  define Spin_system_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/SpinSys.h>		// Includes base class
#include <Basics/ParamSet.h>		// Includes parameter sets
#include <Basics/Gconstants.h>		// Include GAMMA1H definition
#include <Basics/Gutils.h>		// Include GAMMA errors/queries
#include <string>			// Include libstdc++ strings
#include <vector>			// Include libstdc++ STL vectors

class spin_system: public spin_sys
{
private:

  std::vector<double> cshifts;		// Isotropic chemical shifts
  std::vector<double> gfacts;		// Isotropic g-factors
  std::vector<double> Jcouplings;	// Isotropic J coupling constants
  std::vector<double> Acouplings;	// Isotropic A coupling constants
  std::vector<int>    _spflags;		// Flags to tag spin pairs
  double Omega1H;			// Spectrometer (1H) frequency (MHz)
  static double DefOm;			// Default spectrometer frequency

  
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                  CLASS SPIN SYSTEM ERROR HANDLING
// ____________________________________________________________________________

        // Input                sys     : Spin system (this)
        //                      eidx    : Error index
        //                      noret   : Flag for linefeed (0=linefeed)
        //                      pname   : string in message
     
void SYSTerror(int eidx, int noret=0) const;
void SYSTerror(int eidx, const std::string& pname, int noret=0) const;
volatile void SYSTfatality(int eidx) const;
volatile void SYSTfatality(int eidx, const std::string& pname) const;

// ____________________________________________________________________________
// ii                   CLASS SPIN SYSTEM SETUP FUNCTIONS
// ____________________________________________________________________________


virtual int setSsys(const ParameterSet& pset,int idx=-1,int wrn=2);

	// Input		sys      : Spin system (this)
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

protected:

bool setOm(const ParameterSet& pset);

        // Input                sys     : A spin system (this)
        //                      pset    : A parameter set
        // Output               TF      : Spin system field strength
        //                                is set from parameter in pset
        //                                True if set from pset, False
        //                                if set to a default value

// ____________________________________________________________________________
// iii              CLASS SPIN SYSTEM SPIN CHECKING FUNCTIONS
// ____________________________________________________________________________

bool checkNotE(int spin, int warn=1) const;

        // Input                sys     : Spin system (this)
        // 			spin	: Spin index
	//			warn    : Flag if fatal error can occur
	//				      0 = no warnings
	//				      1 = non-fatal warnings
	//				     >1 = fatal error
        // Output               T/F	: TRUE if spin isn't electon
	//				  FALSE otherwise


bool checkNotN(int spin, int warn=1) const;

        // Input                sys     : Spin system (this)
        // 			spin	: Spin index
	//			warn    : Flag if fatal error can occur
	//				      0 = no warnings
	//				      1 = non-fatal warnings
	//				     >1 = fatal error
        // Output               T/F	: TRUE if spin isn't nucleus
	//				  FALSE otherwise


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ____________________________________________________________________________
// A                   SPIN SYSTEM CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________
///Center Spin System Algebraic
///F_list spin_system		- Constructor
///F_list ~			- Destructor
///F_list =			- Assignment	

MSVCDLC               spin_system(int spins=0);
MSVCDLC               spin_system(const spin_system &sys);
MSVCDLL virtual       ~spin_system ();
MSVCDLL spin_system&  operator= (const spin_system &sys);

// ____________________________________________________________________________
// B                CHEMICAL SHIFT AND G-FACTOR MANIPULATIONS
// ____________________________________________________________________________

///Center Spin System Chemical Shift Functions 
///F_list shifts	- Set all chemical shifts
///F_list shift		- Set/Retrieve a specific chemical shift
///F_list maxShift	- Retrieve largest chemical shift
///F_list minShift  	- Retrieve smallest chemical shift
///F_list medianShift 	- Retrieve average chemical shift
///F_list offsetShifts	- Offset all chemical shifts
///F_list PPM		- Set/Retrieve chemical shift in PPM
///Center Spin System Electron G-Factor Functions 
///F_list gfactor	- Set/Retrieve a spin's g-factor
///F_list eshift	- Retrieve e- shift relative to free e- at set Bo
///F_list eshift_lab	- Retrieve e- shift relative to free e- at set Bo L.Fr.

// ************* Chemical Shift Manipulations in Hertz & PPM ******************

MSVCDLL virtual void   shifts(double shift=0);
MSVCDLL virtual void   shift(int, double);
MSVCDLL virtual double shift(int)                  const;
MSVCDLL         double maxShift()                  const;
MSVCDLL         double maxShift(const std::string& Iso) const;
MSVCDLL         double minShift()                  const;
MSVCDLL         double minShift(const std::string& Iso) const;
MSVCDLL         double medianShift()               const;
MSVCDLL         double lab_shift(int)              const;	// Typically ~10^8 !
MSVCDLL virtual void   offsetShifts(double OF,int i=0); 
MSVCDLL virtual void   offsetShifts(double OF,const std::string& Iso); 
MSVCDLL virtual void   PPM (int, double);
MSVCDLL         double PPM (int) const;

//**************************** g-factors For EMR ******************************

/* These functions allow users to obtain or set electon g-factors.  These are
   NOT applicable to nuclear spin! The g-factor is a unitless quantity that
   will typically have a value near that of g for a free electron: 2.00232.
   Frequency <-> Field conversions can be made based on this value using the
   conversion factor Gauss2Hz which is planks constant (6.6262e-27 erg-s/cycle)
   over the Bohr magneton (9.2471e-21 erg/G). Gauss2Hz = 0.714474e-6 G/Hz.
   GAMMA also can track an "electron shift" which is the electon resonance at
   a particular filed strength relative to a base Larmor frequency of a free
   electron.  Note that e- has a negative gyromagnetic ratio, thus electrons
   with a positive shift are shielded!                                       */

MSVCDLL double gfactor(int    spin) const;
MSVCDLL void   gfactor(int    spin, double g);
MSVCDLL double eshift(int     spin) const;
MSVCDLL double lab_eshift(int spin) const;
MSVCDLL double efield(int     spin) const;
MSVCDLL double efield_lab(int spin) const;

// ____________________________________________________________________________
// C           SCALAR & HYPERFINE COUPLING CONSTANT MANIPULATIONS
// ____________________________________________________________________________

///Center Spin System Scalar Coupling Constant Functions 
///F_list Js		        - Set all coupling constants in Hz
///F_list J		        - Set/Retrieve coupling constant in Hz
///Center Spin System Hyperfine Coupling Constant Functions 
///F_list As		        - Set all coupling constants in Gauss
///F_list A		        - Set/Retrieve hyperfine constant in Hz

// ----------------------- Scalar Coupling Functions --------------------------

        // Input                sys     : A spin system (this)
        //			Jval	: Scalar coupling value (Hz)
	// 			int	: Spin index [0, nspins)
	//       		int	: Spin index [0, nspins)
	// Output Js		none	: Sets all J couplings to Jval
	//        J(i,j,d)      none	: Sets (i,j) J coupling to Jval
	//        J(d,i,j)      none	: Sets (i,j) J coupling to Jval
	//        J(i,j)	double	: Returns (i,j) J coupling in Hz

MSVCDLL virtual void   Js(double Jval=0);
MSVCDLL virtual void   J(int, int, double);
MSVCDLL virtual void   J(double, int, int);
MSVCDLL         double J(int, int) const;


// --------------------- Hyperfine Coupling Functions -------------------------

        // Input                sys     : A spin system (this)
	// 			Aval	: Hyperfine coupling in Gauss
	// 			int	: Spin index [0, nspins)
	//       		int	: Spin index [0, nspins)
	// Output As		none	: Sets all HF couplings to Aval
	//        A(i,j,d)      none	: Sets (i,j) HF coupling to Aval
	//        A(d,i,j)      none	: Sets (i,j) HF coupling to Aval
	//        A(i,j)	double	: Returns (i,j) HF coupling in G
	//        AHz(i,j)	double	: Returns (i,j) HF coupling in Hz

MSVCDLL virtual void   As(double Aval=0);
MSVCDLL virtual void   A(int, int, double);
MSVCDLL virtual void   A(double, int, int);
MSVCDLL         double A(int, int)   const;
MSVCDLL         double AHz(int, int) const;

// ____________________________________________________________________________
// D                SPECTROMETER FREQUENCY & FIELD MANIPULATIONS
// ____________________________________________________________________________

///Center Spin System Spectrometer Frequency Functions 
///F_list Omega		     - Set/Retrieve spectrometer frequency
///F_list OmegaAdjust	     - Set spectrometer frequency, PPM shifts static 
///F_list Bo		     - Set spectrometer field strength in Tesla
 
/* These functions allow the user to set or obtain the spectrometer field
   strength.  Typically this is done based on the proton Larmor frequency,
   using the function Omega.  To simulate working on a 500 MHz NMR spectrometer
   you simply set the system Omega to 500.  Additional functions allow users
   to set the field strength based on a resonance frequency of ANY isotope
   (if 1H doesn't suit your fancy).  Last but not least you can access, but not
   set, the field strength in Tesla using the function Bo. Function Bo will of
   course always return a non-negative number.
 
   --> IMPORTANT NOTE: (Re-)Setting the field strength using an Omega function
       DOES NOT CHANGE EXISTING SYSTEM SHIFTS IN HZ!!!  That means that it does
       change the system shifts in PPM.  In contrast, the function OmegaAdjust
       KEEPS THE SYSTEM SYSTEM SHIFTS IN PPM CONSTANT. For example, say you
       which to do a field study by looping through a simulation at different
       field strengths ---> use OmegaAdjust to change the field strength so
       so that you system shifts will be constant on a PPM scale but properly
       change on a Hz scale to reflect altered field.
 
       Function      Return   Arguments            Results
      ----------     ------   ---------   -------------------------------------
       Omega         void        Om       Set 1H Larmor frequency to Om (MHz)
       Omega         void      Om, iso    Set iso Larmor frequency to Om (MHz)
       Omega         double     spin      Get Larmor frequency to spin
       Omega         double     iso       Get Larmor frequency for isotope iso
       Bo            double               Get Field Strength (Tesla)
       OmegaAdjust   double      Om       Set 1H Larmor freq. to Om (MHz)    */

MSVCDLL void   Omega(double Om);				// Om : 1H spect. freq.  (MHz)
MSVCDLL void   Omega(double Om, const std::string& iso);	// Om : iso spect. freq. (MHz)
MSVCDLL double Omega(int spin=-1) const;			// Return Larmor of spin (MHz)
MSVCDLL double Omega(const std::string& iso) const;	// Return Larmor of iso  (MHz)
MSVCDLL double Bo() const;					// Return Field (Gauss)
MSVCDLL void   OmegaAdjust(double Om);			// Reset 1H Larmor, PPM static
MSVCDLL void   FieldAdjust(double B);			// Reset Field Strength

/* These functions work but are now deprecated, they're identical to Omega. */
        
MSVCDLL void   spectrometer_frequency(double freq);
MSVCDLL double spectrometer_frequency() const;

// ____________________________________________________________________________
// E                      SPIN PAIR FLAG FUNCTIONS
// ____________________________________________________________________________

// --------- Functions Which Get/Set Spin Pair Flags Within The System --------

MSVCDLL void spflags(int TF);

        // Input                   TF : TRUE/FALSE status
        // Output                none : All spin pairs  have their flags set
        //                              to TRUE/FALSE
        ///F_list flags               - Set all spin pair flags TRUE or FALSE


MSVCDLL void spflag(int spin1, int spin2, int TF);

        // Input               spin1 : Spin in the system
        //                     spin2 : Another spin in the system
        //                      int  : TRUE/FALSE status 
        // Output               none : Flag for spin pair spin1-spin2 has
	//			       its spin pair flag set to TF
        ///F_list flag               - Set or retrieve specific spin pair flag T/F status.

 
MSVCDLL int spflag(int spin1, int spin2) const;
 
        // Input               spin1 : Spin in the system
        //                     spin2 : Another spin in the system
        // Output               int  : Current T/F status of spin
	//			       pair spin1-spin2

// ____________________________________________________________________________
// F                   SPIN SYSTEM ASSOCIATED FUNCTIONS
// ____________________________________________________________________________
///Center Spin System Associated Functions 

MSVCDLL double center(int spin=0);

	// Input	sys	: Spin system(this)
	// 		spin	: Spin number
	// Return	centerf : Spin system frequency center based
	//			  on the shifts of all spins of the
	//			  same iosotope type as input spin
	///F_list center        - Retrieve spin system spectrum center frequency

// ____________________________________________________________________________
// G			  Nyquist Frequency Functions
// ____________________________________________________________________________

	// Input	sys	: Spin system (this)
	//		<***>   : Selectivity
	//		fact    : Extension Factor (%)
	//		lwhh    : Anticipated half-height linewidth
	// Return	Nyqf    : Approximate Nyquist frequency needed to
	//			  insure proper quadrature acquisitions
	//			  without foldover.  This is based on the
	//			  spin system shifts and coupling constants
	//			  for the spins of the same isotope type as
	//			  the input spin
	// Note		     	: Returned value is in Hertz
	// Note		     	: For broad peaks due to relaxation or
	//			  severe apodization this can return value 
	//			  too small if no lwhh is set.  Compensation
	//			  can be made externally by adjusting the
	//			  returned value.

	// Overloads	<***>	: Selectivity 
	// 1. int	 spin   : Based on this spin's isotope type
	// 2. string	 iso    : Based on this isotope designation
	// 3. Isotope	 iso    : Based on this isotope


MSVCDLL double Nyquist(int spin,               double fact, double lwhh) const;
MSVCDLL double Nyquist(const std::string& iso, double fact, double lwhh) const;
MSVCDLL double Nyquist(const Isotope& iso,     double fact, double lwhh) const;

// ____________________________________________________________________________
// H                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

///Center Spin System Parameter Set Functions 
///F_list =		      - Conversion of system to parameter set
///F_list +=		      - Add system parameters to parameter set

// ----------------------------------------------------------------------------
//            Functions To Make A Parameter Set From A Spin System
// ----------------------------------------------------------------------------
 
MSVCDLL operator ParameterSet( ) const;

	// Input		ss    : A spin system (this)
	// Output		pset  : A parameter set with
	//			        only spin system parameters

MSVCDLL friend void operator+= (ParameterSet& pset, const spin_system &ss);

	// Input		ss    : spin system
	//  			pset  : parameter set
	// Output		pset  : parameter set with
	//			        only spin system parameters

MSVCDLL virtual void PSetAdd(ParameterSet& pset, int idx=-1) const;

        // Input                ss      : A spin system
        //                      pset    : Parameter set
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        // Output               void    : Spin system parameters are
        //                                are added to the parameter set
        //                                with interaction prefix idx
        // Note                         : The parameter names & types
        //                                output here MUST match those used
        //                                in setting spin systems
        //                                from parameters sets

// ----------------------------------------------------------------------------
//            Functions To Make A Spin System From A Parameter Set
// ----------------------------------------------------------------------------

MSVCDLL void setJs(const ParameterSet& pset);

        // Input                ssys    : A spin system (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system scalar couplings
        //                                are set from parameters in pset


MSVCDLL void setAs(const ParameterSet& pset);

        // Input                ssys    : A spin system (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system hyperfine couplings
        //                                are set from parameters in pset
        // Note                         : This restricts the setting of A(i,j)
        //                                to have j>i for all spin pairs
        // Note                         : Scalar couplings are NOT read
        //                                in with this function
 

MSVCDLL void setShifts(const ParameterSet& pset);
        // Input                ssys    : A spin system (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system isotropic shifts
        //                                are set from parameters in pset
        // Note                         : If the system spectrometer
        //                                frequency is not present, PPM
        //                                shifts are disabled


MSVCDLL void setGs(const ParameterSet& pset);

        // Input                ssys    : A spin system (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system isotropic g-factors
        //                                are set from parameters in pset
        //                              : A field strength must have been
        //                                specified before this function


MSVCDLL virtual void operator= (const ParameterSet& pset);

	// Input		sys   : Spin system (this)
	// 			pset  : A parameter set
	// Output		none  : Spin system filled with
	//				parameters n pset
	// Note			      : Functions which place a spin_sys
	//				into a parameter set must contain
	//				add the information read here
 
// ----------------------------------------------------------------------------
//     Functions To Output Isotropic System To ASCII From A Parameter Set
// ----------------------------------------------------------------------------

        // Input                ss      : Spin system (this)
        //                      filename: Output file name
        //                      ofstr   : Output file stream 
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : Spin system is written as a
        //                                parameter set to file filename
	//				  or into output filestream ostr
        ///F_list write                 - Write system to file (as  pset).
        // Note                         : This depends on function PSetAdd!

MSVCDLL virtual int write(const std::string &filename, int idx=-1, int warn=2) const;
MSVCDLL virtual int write(std::ofstream& ofstr,        int idx=-1, int warn=2) const; 

// ____________________________________________________________________________
// I                        SPIN SYSTEM INPUT FUNCTIONS
// ____________________________________________________________________________
 
///F_list read		- Read spin system from disk file
///F_list ask_read	- Ask for file, read spin system from file


	// Input		sys      : Spin system (this)
	// 			filename : Input filename
	//			idx	 : Parameter index value used for
	//				   prefix [#] in input names
	//			warn	 : Warning output level
	//					0 = no warnings
	//					1 = warnings
	//				       >1 = fatal warnings
	// Output		none	 : Spin system is filled with
	//				   parameters read from file
	// Note			 	 : The file should be an ASCII file
	//				   containing recognized sys parameters


MSVCDLL virtual int read(const std::string& fn,    int idx=-1, int warn=2);
MSVCDLL virtual int read(const ParameterSet& pset, int idx=-1, int warn=2);

	// Input		sys      : Spin system (this)
	// 			pset	 : A parameter set
	//			idx	 : Parameter index value used for
	//				   prefix [#] in input names
	//			warn	 : Warning output level
	//					0 = no warnings
	//					1 = warnings
	//				       >1 = fatal warnings
	// Output		TF	 : Spin system is filled with
	//				   parameters in pset
	//				   TRUE if system filled properly


MSVCDLL virtual std::string ask_read(int argc, char* argv[], int argn);
MSVCDLL virtual std::string ask_read(int argc, char* argv[], int argn,
                                                    const std::string& def);

	// Input		sys     : A spin system (this)
	//			argc	: Number of arguments
	//			argv    : Vecotr of argc arguments
	//			argn    : Argument index
	// Output		void    : The parameter argn of array argc
	//				  is used to supply a filename
	//				  from which the spin system is read
	//				  If the argument argn is not in argv,
	//				  the user is asked to supply a filename
        //                                The set filename is returned 
	// Note			 	: The file should be an ASCII file
	//				  containing recognized sys parameters
	// Note			 	: The spin system is modifed (filled)


// ____________________________________________________________________________
// I                       STANDARD I/O FUNCTIONS
// ____________________________________________________________________________

	///F_list print		         - Print to an output stream
        ///F_list <<                     - Send spin system to an output stream

        // Input                sys     : A spin system (this)
        //                      ostr    : Output stream
        // Output printvs	none    : Spin system chemical shifts
        //                                sent into the output stream
        // Note   printvs		: Printed in both Hz and PPM if
        //                                the field strength is set
        // Output printGs	none    : Spin system g-factors are
        //                                sent into the output stream
        // Output printJs       none    : Spin system scalar coupling info
        //                                is sent into the output stream
	// Output printAs	none    : Spin system hyperfine coupling info
        //                                is sent into the output stream
        // Output printO	none    : Spin system field strength
        //                                info sent to the output stream
 

MSVCDLL std::ostream& printvs(std::ostream& ostr) const;
MSVCDLL std::ostream& printGs(std::ostream& ostr) const;
MSVCDLL std::ostream& printJs(std::ostream& ostr) const;
MSVCDLL std::ostream& printAs(std::ostream& ostr) const;
MSVCDLL std::ostream& printO(std::ostream&  ostr) const;
 
        // Input                out      : output stream
        // Output               none	 : modifies output stream

MSVCDLL virtual std::ostream& print(std::ostream& out, bool hdr=true) const;
MSVCDLL friend  std::ostream& operator<<(std::ostream& out, const spin_system& sys);

//-----------------------------------------------------------------------------
//                    Strings Used For Generic Output Functions
//-----------------------------------------------------------------------------

MSVCDLL virtual std::vector<std::string> SYSStrings(int w1=10,int w2=12,int w3=1) const;
MSVCDLL std::vector<std::string> VStrings(int   colwd=12, int digs=2) const;
MSVCDLL std::vector<std::string> PPMStrings(int colwd=12, int digs=2) const;
MSVCDLL std::vector<std::string> GFStrings(int  colwd=12, int digs=2) const;
MSVCDLL std::vector<std::string> BeStrings(int  colwd=12, int digs=2) const;
MSVCDLL std::vector<std::string> JStrings(int   colwd=12, int digs=2) const;
MSVCDLL std::vector<std::string> AStrings(int   colwd=12, int digs=2) const;
MSVCDLL std::vector<std::string> OmStrings(int  colwd=12, int digs=2) const;

};
 
#endif							// SpinSystem.h
