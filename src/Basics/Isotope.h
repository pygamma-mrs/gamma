/* Isotope.h ****************************************************-*-c++-*-
**									**
**                              G A M M A				**
**									**
**	Isotope 		           Interface Definition		**
**							`	 	**
**	Copyright (c) 1999					 	**
**	Tilo Levante, Scott A. Smith				 	**
**	Eidgenoessische Technische Hochschule			 	**
**	Labor fur physikalische Chemie			 		**
**	8092 Zurich / Switzerland			 		**
**						 			**
**      $Header: $
**						 			**
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**							 		**
** Class Isotope is the working spin isotope object in GAMMA programs.	**
** Each isotope object allows for access to many intrinsic spin isotope	**
** properties (see class IsotopeData).  The class also manages a list	**
** of spin isotopes available to GAMMA programs, each Isotope being a	**
** pointer to a particular entry in the isotope list.			**
**									**
*************************************************************************/

#ifndef   GIsotope_h_				// Is file already included?
#  define GIsotope_h_ 1				// If no, then remember it 
#  if defined(GAMPRAGMA)			// If it is the GNU compiler
#    pragma interface				// then this is the interface
#  endif

#include <Basics/IsotopeData.h>			// Include Isotope Data (list entry)
#include <Basics/ParamSet.h>			// Include GAMMA parameter sets
#include <fstream>				// Include file stream handling
#include <vector>				// Include libstdc++ STL vectors
#include <GamGen.h>				// Know MSVCDLL (__declspec)

class Isotope
  {

protected:					// Give derived classes access

  static std::vector<IsotopeData> Isotopes;	// Pointer to isotopes list
  static double                   RELFRQ1H;	// Relative 1H Larmor frequency
         int                      Iso;		// An isotope index
   
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// i                      CLASS ISOTOPE ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   I      : Isotope (this)
                               eidx    : Error index
                               noret   : Flag for linefeed (0=linefeed)
                               pname   : Added error message for output
       Output                  none    : Error message output
                                         Program execution stopped if fatal

  The following error messages use the defaults set in the Gutils package

               Case                          Error Message

               (0)                     Program Aborting.....
               (1)                     Problems With Input File Stream
               (2)                     Problems With Output File Stream
               default                 Unknown Error                        */

         void Isoerror(int eidx,                           int noret=0) const;
         void Isoerror(int eidx, const std::string& pname, int noret=0) const;
volatile void Isofatal(int eidx)                                        const;
volatile void Isofatal(int eidx, const std::string& pname)              const;
 
// ____________________________________________________________________________
// ii               CLASS ISOTOPE DEALINGS WITH ISOTOPE LISTS
// ____________________________________________________________________________
 
/* The code for this next function was generated from a program which parsed
   the original GAMMA file Isotopes.list.  Users should keep in mind that
   since our Isotopes.list may be "incomplete" in that we don't have
   all values for all isotopes, we'll use a default value for anything
   that we don't know.  Since most values in the list are positive, the
   default will be negative!  In those rare cases where an isotope does
   have some negetive values (15N for example) we'll just hope that they
   won't have the ISODEFVAL I've picked.  If they do, just change the
   number for ISODEFVAL to any negative number not in the list and we'll
   be just fine and still know what an "undefined" value is......           */
 
/* const double ISODEFVAL = -1.1e6;                     // Default value    */
 
 
void set_Isotope_list();
 
        // Input                I : A dummy isotope (this)
        // Output            void : Fills isotopes list with all spins
 
void SetRel1HF();

        // Input                I : A dummy isotope (this)
        // Output            void : Sets the relative 1H frequency
        // Note                   : Must have filled isotopes list first!

/* GAMMA's isotopes, rather than having a known gyromagnetic ratio, have
   a relative Larmor frequency from which gyromagnetic ratio may be
   calculated.  The base relative Larmor frequency used in GAMMA is that
   of a proton, RELFRQ1H, and to date will be 400.13 (MHz) in GAMMA.
   The gyromagnetic ratio of isotope type i relates to its relative
   frequency RelFreq according to
 
        gamma  = RelFreq*GAMMA1H/RELFRQ1H = Omega * gamma  / Omega
             i                                   i       1H       1H
 
   GAMMA has the DEFAULT value GAMMA1H = 2.67519e8 (1/T-sec) set in its
   constants module (Gutils).  It sets the default proton relative frequency
   herein because that value (usually 400.13) is imported into GAMMA from a
   stardard isotopes list.  Were that list to change, perhaps the listed base
   frequencys would change. But that is O.K. because it will be set here
   accordingly.
 
   Note however that this does NOT affect any GAMMA programs that wish to set
   alternate field strengths through spin systems.  The relative frequencies
   in GAMMA's isotope list are Larmor frequencies but used ONLY to determine
   gamma values. Larmor frequencies in MR simulations are set independently
   by the user in typical GAMMA programs.                                    */


// ____________________________________________________________________________
// iii              Class Isotope Private Parameter Set Functions
// ____________________________________________________________________________

bool SetIsotope(const ParameterSet& pset, int idx=-1, bool warn=true);

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ____________________________________________________________________________
// A                    ISOTOPE CONSTRUCTION, ASSIGNMENT
// ____________________________________________________________________________

MSVCDLC          Isotope();
MSVCDLC          Isotope(const         Isotope& I);
MSVCDLC          Isotope(const    std::string&  I);
MSVCDLL Isotope& operator= (const      Isotope& I);
MSVCDLC virtual  ~Isotope();

        // Input               I : An isotope (this)
        // Output           none : Isotope is destructed (do nothing)
        // Note                  : The isotopes list, Isotopes, does use
        //                         additional memory but will NOT be
        //                         destroyed until a program is terminated!

// ____________________________________________________________________________
// C                        ISOTOPE ACCESS FUNCTIONS
// ____________________________________________________________________________

        // Input                I : An Isotope (this)
	// Note			  : All functions use IsotopesData 
	//			    except "gamma" which is calculated based
	//			    on stored relative isotope frequencies

/*    Function        Return      Units                Example(s)
      --------     -------------  -----   -------------------------------
         qn        mz   (double)   hbar   0.5, 1.0, 1.5, ....
         HS        2I+1 (int)      none   2, 3, 4, .....
      momentum     mz   (string)   none   1/2, 1, 3/2, .....
       symbol           (string)   none   1H, 2H, 13C, 19F, ....
       name             (string)   none   Hydrogen, Lithium, Carbon, ...
      element           (string)   none   H, Li, F, C, ......
       number           (int)      none   1<-H, 3<-Li, 6<-C, .....
       mass             (int)      amu    1<-1H, 2<-2H, 13<-13C, ....
      weight            (double)   g/m    1.00866<-1H, 7.016<-7Li, ....
       gamma            (double)  1/T-s   2.67519*10^8 
 relative_frequency    	(double)   MHz    400.13, 155.503, ...                */

      MSVCDLL double       qn()                 const;
      MSVCDLL int          HS()                 const;
      MSVCDLL std::string  momentum()           const;
const MSVCDLL std::string& symbol()             const;
const MSVCDLL std::string& name()               const;
const MSVCDLL std::string& element()            const;
      MSVCDLL int          number()             const;
      MSVCDLL int          mass()               const;
      MSVCDLL double       weight()             const;
      MSVCDLL double       gamma()              const;
      MSVCDLL double       receptivity()        const;
      MSVCDLL double       relative_frequency() const;
      MSVCDLL bool         electron()           const;

// ____________________________________________________________________________
// D                        ISOTOPE I/O FUNCTIONS
// ____________________________________________________________________________

/* These are function that will read/write an Isotope designation from/to
   an external ASCII file.                                                   */

virtual MSVCDLL bool read(const std::string& filename, int idx=-1, int warn=2);
virtual MSVCDLL bool read(const ParameterSet& pset,    int idx=-1, int warn=2);

/* These are functions that output formatted information concerning the
   Isotope to a specified output stream of file. 

           Input                I       : An isotope (this) 
                                ostr    : Output stream 
           Output               none    : Isotope info placed into the
                                          output stream ostr                 */

MSVCDLL std::vector<std::string> printStrings(bool hdr=true) const;
virtual MSVCDLL std::ostream& print(std::ostream& ostr) const;
friend  MSVCDLL std::ostream& operator<< (std::ostream&    ostr, const Isotope& I);

// ____________________________________________________________________________
// E                         ISOTOPE LIST FUNCTIONS
// ____________________________________________________________________________

/* Allow users to search for a specific isotope in the list or find whether
   a particular isotope is contained in the list.  Remember that each Isotope
   only contains a pointer to the Isotopes list with an index for which 
   Isotope in the list it actually is. 
 
           Input                I : An isotope (this) 
                               ID : A single isotope
                           symbol : A string designating an isotope
           Output               i : Index of isotope in list if it exists.
                                    If id does NOT exist, returns is -1
                               TF : TRUE if isotope of type "symbol" known
	   Note			  : Will return FALSE if symbol not found
	  			    due to no isotope list                    */

virtual MSVCDLL int                      seek(const   IsotopeData& ID);
virtual MSVCDLL bool                     exists(const std::string& symbol);
static  MSVCDLL bool                     known(const  std::string& symbol);
static  MSVCDLL int                      size();
static  MSVCDLL std::vector<std::string> PrintListStrings();
static  MSVCDLL void                     PrintList(std::ostream& ostr, bool hdr=true);

// ____________________________________________________________________________
// F                   Isotope List Expansion Functions
// ____________________________________________________________________________

/* These functions allow users to add additional spins into the working GAMMA
   isotopes list.  This will primarily be used in EPR/ESR/EMR simulations when
   electrons with I>1/2 are required. For example, this will occur when dealing
   with an electron about a metal center when the surroundings effecitvely
   cause it to be split. To add in an additional spin one simply constructs a
   new isotope type then adds/subtracts it with these functions.  Note that one
   may NOT remove any "standard" isotopes in GAMMA nor may one add an isotope
   that already exists.                                                      */

static MSVCDLL bool AddIsotope(const IsotopeData& ID, int warn=2);

// ____________________________________________________________________________
// G                Isotope Container Support Functions
// ____________________________________________________________________________

/* These functions allow for comparisons between spin isotopes. Such functions
   are required under some compilers if one desired to build containers of
   spin isotopes (STL lists, vectors, ....).  Equal spin isotopes will simply
   point to the same entry in the isotopes list.  For sorting purposes we
   go by the spin Hilbert space associated with the isotope.                 */

virtual MSVCDLL bool operator== (const Isotope& I) const;
virtual MSVCDLL bool operator!= (const Isotope& I) const;
virtual MSVCDLL bool operator<  (const Isotope& I) const;
virtual MSVCDLL bool operator>  (const Isotope& I) const;

// ____________________________________________________________________________
// H                    Isotope Auxiliary Functions
// ____________________________________________________________________________

/* These are just handy to have available for other parts of GAMMA           */

MSVCDLL bool nepair(const Isotope& S) const;
MSVCDLL bool enpair(const Isotope& S) const;
MSVCDLL bool eepair(const Isotope& S) const;
MSVCDLL bool nnpair(const Isotope& S) const;


};						// End class Isotope


#endif								// Isotope.h

