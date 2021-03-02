/* SpinSys.h ****************************************************-*-c++-*-
**									**
** 	                           G A M M A				**
**									**
**	Basic Spin System                            Interface		**
**								 	**
**	Copyright (c) 1991, 1992, 1993, 1994				**
**	Scott Smith and Tilo Levante			 		**
**	Eidgenoessische Technische Hochschule	 			**
**	Labor fuer physikalische Chemie		 			**
**	8092 Zuerich / Switzerland		 			**
**							 		**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**									**
** Description								**
**									**
** Class spin_sys defines a quantity which contains a count of the	**
** number of spins and their isotope types.  As it is assumed that	**
** complex spin system classes will be derived from spin_sys, the	**
** following functions have been declared virtual.			**
**									**
**      ~       The spin_sys destructor					**
**      print   The spin_sys printing function				**
**									**
** This allows any any derived classes to overwrite the defaults of	**
** these functions.							**
**									**
*************************************************************************/

///Chapter Class Spin Sys
///Section Overview
///Body    The class Spin Sys defines the basic physical
///        attributes of a system of spins. These essential
///        quantities are the number of spins and their
///        associated spin angular momentum. Functions are
///        provided for simplified access to all spin sys
///        quantities and for disk I/O.
///Section Available Spin Sys Functions

#ifndef  SpinSys_h_			    // Is this file already included?
#define  SpinSys_h_ 			// If no, then remember it

#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Basics/Isotope.h>		// Include Isotope knowledge
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Matrix/row_vector.h>		// Include row vectors
#include <Matrix/matrix.h>		// Include matrices
#include <HSLib/Basis.h>		// Include bases
#include <string>			// Include stdlibc++ strings
#include <vector>			// Include STL vectors
#include <list>				// Inlcude STL list class

#if defined(_MSC_VER) || defined(__SUNPRO_CC)	// Microsoft VC++ does not handle vectors
  typedef std::vector<int> flagvec;		// of bool due to lack of defined < & >
#else						// operators.  So we used int vectors.
  typedef std::vector<bool> flagvec;
#endif
  
class spin_sys
  {
  int nspins;                      // Number of spins in the system
  flagvec spinflags;               // Flags to tag spins (used by SpinOp)
  std::vector<Isotope> Isotopes;   // The spin isotopes in the system
  std::string sysname;             // The spin system name - max. 50 characters
  static int _warn;		   // Flag whether warning messages
  static std::string DefIso;       // Default isotope type
  basis bsmx;			   // A default basis matrix defined by spin HS

 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                     CLASS SPIN SYSTEM ERROR HANDLING
// ____________________________________________________________________________

        // Input                sys     : Spin system (this)
        //                      eidx    : Error index
        //                      noret   : Flag for linefeed (0=linefeed)
        //                      pname   : string in message

void error(int eidx, int noret=0) const;
void error(int eidx, const std::string& pname, int noret=0) const;
volatile void fatality(int eidx) const;

// ____________________________________________________________________________
// ii                 CLASS SPIN SYSTEM SETUP FUNCTIONS
// ____________________________________________________________________________

protected:

virtual int setSsys(const ParameterSet& pset, int idx=-1, int warn=2);

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
	//				   the parameter set for a spin_sys
	//				   1.) The number of spins
	//				   2.) An isotope type for each spin
	//				   3.) An optional spin system name
	// Note				 : Functions which place a spin_sys
	//				   into a parameter set must contain
	//				   add the information read here


void setBasis(const matrix& mx);

	// Input		sys	: Spin system (this)
	// 			mx	: A default basis matrix
	// Output		void	: The system basis matrix is set
	// Note				: There is NO dimension checking


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
public:

        // Input                sys     : Spin system (this)
        //                      spin    : Spin index
        //                      spin1   : First spin index
        //                      spin2   : Second spin index
        //                      warn    : Flag if fatal error can occur
        //                                    0 = no warnings
        //                                    1 = non-fatal warnings
        //                                   >1 = fatal error
        // Output check_spin   T/F     : TRUE if 0 <= spin < nspins
        // Output check_spins  T/F     : TRUE if 0 <= spin(1,2) < nspins
        // Note                         : If die is set non-zero, a program
        //                                abort will be sent

int check_spin(int  spin,             int die=1) const;
int check_spins(int spin1, int spin2, int die=1) const;

// ____________________________________________________________________________
// A                SPIN SYSTEM CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

///Center Spin Sys Algebraic
///F_list spin_sys	      - Constructor
///F_list =	              - Assignment

MSVCDLC                   spin_sys();
MSVCDLC                   spin_sys(int    spins);
MSVCDLC                   spin_sys(const  spin_sys& sys);
MSVCDLC virtual           ~spin_sys();
MSVCDLL virtual spin_sys& operator=(const spin_sys& sys);

// ____________________________________________________________________________
// B                              COMPARISONS
// ____________________________________________________________________________

///F_list ==	              - Equality
///F_list !=	              - Inequality

MSVCDLL int operator==(const spin_sys& sys) const;
MSVCDLL int operator!=(const spin_sys &sys) const;

// ____________________________________________________________________________
// C                   BASIC SPIN SYSTEM MANIPULATIONS
// ____________________________________________________________________________

///Center Basic Functions
///F_list spins	              - Number of spins
///F_list spinpairs	      - Number of spin pairs
///F_list HS	              - Retrieve spin or spin system Hilbert space

        // Input                sys   : Spin system (this)
	// Input		spin  : Spin in the spin system
	// Output spins		int   : Number of spins in the system
	// Output spinpairs	int   : Number of spin pairs in the system
	// Output HS		int   : Size of spin system Hilbert Space
	// Output HS(i)     	int   : Size of specific spin Hilbert Space

MSVCDLL int spins()     const;
MSVCDLL int spinpairs() const;
MSVCDLL int HS()        const;
MSVCDLL int HS(int spin) const;

// ____________________________________________________________________________
// D          SPIN ANGULAR MOMENTUM AND ISOTOPE MANIUPLATIONS
// ____________________________________________________________________________

// sosi - the gamma(Iso) function is completely unnecessary! it is a member
//        function of class isotope I think

// ----------------------------------------------------------------------------
//                Functions To Get & Set Spin Isotope Values
// ----------------------------------------------------------------------------

/* Note that in GAMMA one may specify an isotope using the spin symbol as a
   string.  Some examples: 1H, 2H, 13C, 14N, 131Xe, .....  There is only one
   exception currently and that is an electron whose symbol is taken to be
   e- (rather than 0e).

///F_list isotope	- Set or retrieve spin isotope type
///F_list symbol	- Retrieve spin isotope type as a string, e.g. 19F.
///F_list qn		- Retrieve spin angular momentum of a spin or the system
///F_list element	- Retrieve spin element type as a string, e.g. Carbon.
///F_list momentum	- Retrieve spin angular momentum as a string.
///F_list gamma	  	- Retrieve gyromagnetic ration of a spin.

Function  Arguments                                Result
========  ==========  =========================================================
isotope    i, name    Sets spin isotope i to that specified by name (e.g. 23Na)
isotope    i, Iso     Sets spin isotope i to that specified by Iso
isotope      i        Returns the isotope type of spin i (as Isotope)
weight       i        Returns the atomic weight of spin i (in amu)
symbol       i        Returns the symbol for spin i (e.g. 23Na), default DefIso
qn           i        Returns I of spin i in units of hbar (e.g. 0.5, 1.5, ...)
qn                    Returns total Iz of the system in units of hbar
element      i        Returns name of spin i (e.g. Uranium)     
momentum     i        Returns spin angular momentum of spin i as string (1/2)
momentum              Returns total s.a.m. of the system as a string
gamma        i        Returns gyromagnetic ratio of spin i in rad/(sec*T)
gamma       Iso       Return gyromagnetic ration of isotope I 

   Note that spin indices span [0, nspins-1]                                 */

MSVCDLL virtual        void        isotope(int, const std::string&);
MSVCDLL virtual        void        isotope(int, const Isotope&);
MSVCDLL virtual const  Isotope&    isotope(int)                  const;
MSVCDLL                double      weight(int)                   const;
MSVCDLL                std::string symbol(int spin)              const;
MSVCDLL                double      qn(int spin)                  const;
MSVCDLL                double      qn()                          const;
MSVCDLL                std::string element(int spin)             const;
MSVCDLL                std::string momentum(int spin)            const;
MSVCDLL                std::string momentum()                    const;
MSVCDLL                double      gamma(int spin)               const;
MSVCDLL                double      gamma(const std::string& iso) const;
MSVCDLL         const  std::vector<Isotope>& IsoVec()            const;

// ----------------------------------------------------------------------------
//        Functions To Access Spin States & Default Basis Functions
// ----------------------------------------------------------------------------

/* In GAMMA there is a default basis associated with each spin system.  The
   spin system is taken to form a composite Hilbert spin space of dimension  

                                 -----
                            HS =  | |  2I + 1
                                  | |    i
                                i spins
   
   and there are HS spin basis functions.  The default basis for a single
   spin is simply ordered { alpha, beta, gamma, delta, ...... } and the
   default basis of a multi-spin system is obtained from direct products of
   individual spin bases.  The most familar example is for spin 1/2 systems,
   where (using a=alpha, b=beta, ....) the basis functions will be

      2 spins: aa  ab ba bb        3 spins: aaa aab aba abb aab bab bba bbb

   Note that these are not total Fz ordered, they are ordered according to
   the direct product of indiviual spin basis values!  The functions below
   allow users to obtain information regarding these basis functions with
   respect to indiviual spins and the basis as a whole.                      */

///F_list HSvect	- Get vector of spin Hilbert spaces
///F_list qState	- Get state vector of a particular quantum state
///F_list qStates	- Get array of all quantum states
///F_list qnState	- Get total Fz of a state
///F_list qnState	- Get total Fz of all states
///F_list qnDist	- Retrieve a statistic over the quantum numbers
///F_list CoherDist	- Retrieve a statistic over the different coherences

MSVCDLL std::vector<int> HSvect() const;
MSVCDLL row_vector       qState(int state) const;

	// Input                sys	   : Spin system (this)
	// 			state	   : State index [0, HS)
	// Output               row_vector : Vector with elements for each spin
	//				     having its quantum value in state
	// Note				   : If sys were two spin 1/2 returns
	//				        state  return vector 
	//					  0       0.5,  0.5
	//					  1       0.5, -0.5
	//					  2      -0.5,  0.5
	//					  3      -0.5, -0.5
	
MSVCDLL matrix qStates() const;

	// Input                sys	: Spin system (this)
	// Output 		matrix	: Matrix (HS x ns) where each 
	//				  row i constains qState(i)
	// Note				: If sys were two spin 1/2 returns
	//					       [  0.5,  0.5 ]
	//					  mx = |  0.5, -0.5 ]
	//					       | -0.5,  0.5 ]
	//					       [ -0.5, -0.5 ]


MSVCDLL double qnState(int state) const;

	// Input                sys	: Spin system (this)
	// 			state	: State index [0, HS)
	// Output		double	: Quantum # of a state (total Fz)
	// Note				: If sys were two spin 1/2 returns
	//				        state    return
	//					  0         1
	//					  1         0
	//					  2         0
	//					  3        -1


MSVCDLL col_vector qnStates() const;

	// Input                sys	: Spin system (this)
	// Output		cvec	: A column vector of HS elements
	//				  containing all states total Fz
	// Note				: If sys were two spin 1/2 returns
	//				                       t
	//					[ 1, 0, 0, -1 ]


MSVCDLL row_vector qnDist() const;

	// Input                sys	: Spin system (this)
	// Output		rvec    : Row vector with (2*qn()+1) elements
	//				  containing the number of states
	//				  with the quantum number q = qn()-i
	//				  where i is the index in the vector


MSVCDLL row_vector CoherDist() const;

	// Output		row_vector : Vector of size (4*qn()+1). Each
	//				     element contains the number of
	//				     transitions per coherence order
	//				     The ZQC is in the middle.


/*        Functions To Get Quick Spin System Information      

   These are helper functions often used in GAMMA programs for TF tests of
   the system.  Allows users to quickly find is the system is heteronuclear,
   whether all spin 1/2, whether it contains electrons, etc.                 */

///F_list homonuclear         - Tests if sys homonuclear
///F_list heteronuclear       - Tests if sys heteronuclear
///F_list spinhalf            - Tests if sys has all spins I=1/2
///F_list electrons           - Number of electrons in sys
///F_list enpair              - Test electron/nucleus pair
///F_list isotopes            - No. unique isotopes

/*    Function   Arguments Output                     Comments
   ------------- --------- ------  --------------------------------------------
    homonuclear    ---      bool   True if all spins in system homonuclear
   heteronuclear   ---      bool   True if any spin types s in system differ
     electron       i       bool   True if spin is an electron
     nucleon       i       bool   True if spin is an nucleon
     spinhalf      ---      bool   True if all spin in system have I=1/2
     electrons     ---      int    Returns the number of electron in system
      nepair       i,j      bool   True if spin i & j are a nucleus/e- pair
      enpair       i,j      bool   True if spin i & j are a nucleus/e- pair
     pairidx       i,j      int    Index (dipolar) for spin pair i&j
     isotopes      ---      int    Number of unique isotopes in system
     isotopes       i      string  Isotope type of (isotope) index i
     isotopes     string    bool   True if isotope of type input in system   */

MSVCDLL bool        homonuclear()                  const;
MSVCDLL bool        heteronuclear()                const;
MSVCDLL bool        electron(int i)                const;
MSVCDLL bool        nucleon(int i)                 const;
MSVCDLL bool        spinhalf()                     const;
MSVCDLL int         electrons()                    const;
MSVCDLL int         nucleons()                     const;
MSVCDLL bool        nepair(int i,  int j)          const;
MSVCDLL bool        enpair(int i,  int j)          const;
MSVCDLL bool        eepair(int i,  int j)          const;
MSVCDLL bool        nnpair(int i,  int j)          const;
MSVCDLL int         pairidx(int i, int j)          const;
MSVCDLL int         isotopes()                     const;
MSVCDLL std::string isotopes(int idx)              const;
MSVCDLL bool        isotopes(const std::string& I) const;

// ____________________________________________________________________________
// E                           SPIN FLAG FUNCTIONS
// ____________________________________________________________________________

/* --------------- Functions Which Get/Set Flags Within The System ------------

   For each spin in the system there exists an On-Off "flag" which can be
   used in functions that act on a set of chosen spins.  They do not have
   any other quality that affects the internal workings of the class. 

   Function     Argument(s)   Output               Comments

  SetFlag          i,TF        void   Sets spin i's flag to TF
  SetFlags          TF         void   Sets all spin flags to TF
  SetFlags      isoin,TF       void   Sets flags for spins of type isoin to TF
  SetFlags        iso,TF       void   Sets flags for isotopes of type iso to TF
  GetFlag           i          bool   Gets flag for spin i
  GetFlags                    vector  Gets spin system spin flags
  GetFlags         TF         vector  Get vector of flags for sys all set to TF
  GetFlags      i,TF,DTF      vector  Get vector of sys flags, all DTF but i=TF
  GetFlags    isoin,TF,DTF    vector  Get vector of sys flags, isoin?TF:DTF  */

MSVCDLL void    SetFlag(int spin, bool TF);
MSVCDLL void    SetFlags(bool TF);
MSVCDLL void    SetFlags(const std::string& isoin, bool TF);
MSVCDLL void    SetFlags(const Isotope& Iso, bool TF);
MSVCDLL bool    GetFlag(int i) const;
MSVCDLL flagvec GetFlags() const;
MSVCDLL flagvec GetFlags(bool TF) const;
MSVCDLL flagvec GetFlags(int spin, bool TF, bool DefTF=0) const;
MSVCDLL flagvec GetFlags(const std::string& isoin,bool TF,bool DTF=0) const;
MSVCDLL flagvec GetFlags(const Isotope& isoin,bool TF,bool DTF=0) const;

// ____________________________________________________________________________
// F                          SPIN SYSTEM NAME
// ____________________________________________________________________________

	// Input		spin sys : spin system
	//			i        : spin index
	//			sysname  : spin system name
	// Output		none     : internal spin system name set
	//			           or name of system returned (i<0)
	//                                 or name of spin i returned (i>=0)
	///F_list name	                 - Set or retrieve spin system name.

MSVCDLL       void         name(const std::string& name);
MSVCDLL const std::string& name(int i=-1) const;

        // Input                sys     : A spin system (this)
        //                      warnf   : Warning flag
        // Output               none    : The function sets/gets the intermal
        //                                value of _warn to warnf
        // Note                         : _warnf determines wheteher some
        //                                non-fatal problems issue a warning
        //                                message or not

MSVCDLL void warnings(int warnf);
MSVCDLL int  warnings() const;


MSVCDLL std::string IsoDefault();
 
        // Input                sys     : A spin system (this)
        // Output               none    : The function returns a string
        //                                of the default isotope type

 
MSVCDLL void IsoDefault(const std::string& DI);
 
        // Input                sys     : A spin system (this)
	//			DI	: string for an isotope type
        // Output               none    : The function sets the intermal
        //                                default isotope type to DI
 
// ____________________________________________________________________________
// G                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//        Functions To Make A Parameter Set From A Solid Spin System
// ----------------------------------------------------------------------------

///Center Spin Sys I/O
///F_list =		- Assignment of a spin sys to/from a pset
///F_list +=		- Addition of a spin sys to a pset
///F_list write		- Write system to file (as  pset).

	// Input		sys	: A spin system (this)
	//			pset    : A parameter set
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
	// Output ()		pset	: A parameter set with
	//				  only spin system parameters
	// Output +=		pset	: Parameter set with spin system
	//				  parameters added in
        // Output PSetAdd       void    : Spin system parameters are
        //                                are added to the parameter set
        //                                with interaction prefix idx
        // Note                         : The parameter names & types
        //                                output here MUST match those used
        //                                in setting spin systems
        //                                from parameters sets

MSVCDLL              operator ParameterSet( ) const;
MSVCDLL friend  void operator+=(ParameterSet& pset, const spin_sys &ss);
MSVCDLL virtual void PSetAdd(ParameterSet& pset, int idx=-1) const;

// ----------------------------------------------------------------------------
//        Functions To Make A Solid Spin System From A Parameter Set
// ----------------------------------------------------------------------------


MSVCDLL int getSpins(const ParameterSet& pset, int warn=0) const;
 
        // Input                sys	: A spin system (this)
        //                      pset    : A parameter set
	//			warn	: Warning level
	//					0 = no warnings
	//					1 = warnings
	//				       >1 = fatal warnings
        // Output               ns      : Number of spins specified by 
        //                                parameter NSpins in pset
        // Note                         : This does NOT set the number of spins
	// Note				: Return -1 if # of spins not found
 

MSVCDLL void setName(const ParameterSet& pset);

        // Input                sys     : A spin system (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system name is
        //                                set from parameter in pset
        // Note                         : This is not a required parameter

     
MSVCDLL void setIs(const ParameterSet& pset);
 
        // Input                sys     : A spin system (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system isotope types
        //                                are set from parameters in pset
 

MSVCDLL void operator=(const ParameterSet& pset);

	// Input		sys      : Spin system (this)
	// 			pset     : A parameter set
	// Output		none	 : Spin system filled with
	//				   parameters in pset
	// Note				 : Three things are gleaned from
	//				   the parameter set for a spin_sys
	//				   1.) The number of spins
	//				   2.) An isotope type for each spin
	//				   3.) An optional spin system name
	// Note				 : Functions which place a spin_sys
	//				   into a parameter set must add in
	//				   in the same information read here


// ----------------------------------------------------------------------------
//      Functions To Output Spin System To ASCII From A Parameter Set
// ----------------------------------------------------------------------------

MSVCDLL virtual int write(const std::string &filename, int ix=-1, int wn=2) const;
 
        // Input                ss      : Spin system (base)
        //                      filename: Output file name 
        //                      ix	: Parameter index value used for 
        //                                prefix [#] in output names 
        //                      wn 	: Warning level 
        // Output               none    : Spin system is written as a
        //                                parameter set to file filename

 
MSVCDLL virtual int write(std::ofstream& ofstr, int idx=-1, int warn=2) const;
 
        // Input                ss      : Spin system (base)
        //                      ofstr   : Output file stream
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : Spin system is written as a
        //                                parameter set to output filestream
        // Note                         : This depends on function PSetAdd!


// ____________________________________________________________________________
//                          SPIN SYSTEM INPUT FUNCTIONS
// ____________________________________________________________________________

///F_list read		- Read spin sys from disk file
///F_list ask_read	- Ask for file, read spin sys from file

	// Input		sys      : Spin system (this)
	// 			filename : Input filename
	// 			pset	 : A parameter set
	//			idx	 : Parameter index value used for
	//				   prefix [#] in input names
	//			warn	 : Warning output level
	//					0 = no warnings
	//					1 = warnings
	//				       >1 = fatal warnings
	// Output		none	 : Spin system is filled with
	//				   parameters read from file
	//				   or from pset
	//				   TRUE if system filled properly
	// Note			 	 : The file should be an ASCII file
	//				   containing recognized sys parameters

MSVCDLL virtual int read(const std::string& filename, int idx=-1, int warn=2);
MSVCDLL virtual int read(const ParameterSet& pset,    int idx=-1, int warn=2);

MSVCDLL virtual std::string ask_read(int argc, char* argv[], int argn);
MSVCDLL virtual std::string ask_read(int argc, char* argv[], int argn,
                                                    const std::string& def);

	// Input		sys     : A basic spin system (this)
	//			argc	: Number of arguments
	//			argv    : Vector of argc arguments
	//			argn    : Argument index
	// Output		filename: The parameter argn of array argc
	//				  is used to supply a filename
	//				  from which the spin system is read
	//				  If the argument argn is not in argv,
	//				  the user is asked to supply a filename
	//				  The file name is returned
	// Note			 	: The file should be an ASCII file
	//				  containing recognized sys parameters
	// Note			 	: The spin system is modifed (filled)


// ____________________________________________________________________________
// H                     DEFAULT BASIS FUNCTIONS
// ____________________________________________________________________________

///F_list get_basis	- Spin sys default basis (as a matrix)

MSVCDLL basis get_basis() const;

	// Input		sys      : Spin system
	// Output		bs	 : The default basis defined by
	//				   the spin system Hilbert space
	// Note				 : The basis is output in matrix
	//				   format since actual bases are
	//				   not known to GAMMA on this level 

// ____________________________________________________________________________
// I                     SPIN SPACE MAPPING FUNCTIONS
// ____________________________________________________________________________

/* These functions map the system's single spin and spin pair basis functions
   as well as their coherences to those of the full system.  This is quite 
   useful in switching between reduced and full Hilbert spaces.              */

//-----------------------------------------------------------------------------
//                           Single Spin Mappings
//-----------------------------------------------------------------------------


MSVCDLL matrix BasisMap1() const;

        // Input                sys      : Spin system (this)
        // Output               bmap     : An ns x hs array of basis mappings.
        //                                 For each spin is a row the length
        //                                 of the full Hilbert space.  Each
        //                                 column (system basis function) will
        //                                 contain the spin's corresponding
        //                                 basis function # in the REAL part.
        //                                 The imaginary part will contain the
        //                                 total basis Fz - spin's mz

 
MSVCDLL matrix TransitionMap1() const;

        // Input                sys      : Spin system (this)
        // Output               tmap     : Array (ns x hs) of coherence maps
        //                                 For each spin is a row the length
        //                                 of the full Liouville space.  Each
        //                                 column (system coherences) will
        //                                 contain the spin's corresponding
        //                                 coherence # in the REAL part.
        //                                 The imaginary part will contain the
        //                                 total basis Fz - spin's mz
 
//-----------------------------------------------------------------------------
//                              Spin Pair Mapping
//-----------------------------------------------------------------------------
 
MSVCDLL matrix BasisMap2() const;
 
        // Input                sys      : Spin system (this)
        // Output               bmap     : A nd x hs array of basis mappings.
        //                                 Each spin pair is a row the length
        //                                 of the full Hilbert space.  Each
        //                                 column (system basis function) will
        //                                 contain the spin pair's corresponding
        //                                 basis function # in the REAL part.
        //                                 The imaginary part will contain the
        //                                 total basis Fz - spin pair's mz


MSVCDLL matrix TransitionMap2() const;

        // Input                sys      : Spin system (this)
        // Output               tmap     : An nd x hs array of coherence maps.
        //                                 Each spin pair has a row the length
        //                                 of the full Liouville space.  Each
        //                                 column (system coherences) will
        //                                 contain the spin pair's corresponding
        //                                 coherence # in the REAL part.

// ______________________________________________________________________
// J                       STANDARD I/O FUNCTIONS
// ______________________________________________________________________

///F_list print			- Send spin system to output stream.
///F_list <<			- Send spin system to an output stream
///F_list printstrings		- Vector of strings for printing system

 
//-----------------------------------------------------------------------------
//                        Commonly Used Output Functions
//-----------------------------------------------------------------------------
 
        // Input                out      : output stream;
        // Output               none     : modifies output stream

MSVCDLL virtual std::ostream& print(std::ostream& out, bool hdr=true) const;
MSVCDLL friend  std::ostream& operator<<(std::ostream& out, const spin_sys& sys);

MSVCDLL virtual std::vector<std::string> printstrings() const;
     
        // Input                sys      : Spin system (this)
        // Output               SV       : string vector containing strings
        //                                 which can be used for printing
        //                                 the spin system

//-----------------------------------------------------------------------------
//                    Strings Used For Generic Output Functions
//-----------------------------------------------------------------------------

/* These functons return vectors of strings that can be used in functions that
   print out spin system information. The strings are of a specified width so
   that they can easily form nice columns when printed. The value of colwd sets
   the width of the srings returned. Function SYSStrings will return strings
   for printing the entire spin system. Each string in that case will appear
   as
                |<--w1-->|_:w3|<--w2-->|w3|<--w2-->|w3|<--w2-->|......
   e.g.         Isotope    :      1H          13C         2H

   where the column widths have default and minimal values built in.         */

MSVCDLL virtual std::vector<std::string> SYSStrings(int w1=10,int w2=5,int w3=1) const;
MSVCDLL std::vector<std::string> SIStrings(int  colwd=10) const;
MSVCDLL std::vector<std::string> SYMStrings(int colwd=10) const;
MSVCDLL std::vector<std::string> SAMStrings(int colwd=10) const;

};
 
#endif							// SpinSys.h
