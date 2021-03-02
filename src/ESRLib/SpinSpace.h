/* SpinSpace.h **************************************************-*-c++-*-
**									                                    **
** 	                           G A M M A				                **
**									                                    **
**	Spin Hilbert Space                                 Interface		**
**								 	                                    **
**	Copyright (c) 2000				                                    **
**	Scott A. Smith                                                      **
**  National High Magnetic Field Laboratory                           	**
**  1800 E. Paul Dirac Drive                                          	**
**  Box 4005                                                          	**
**  Tallahassee, FL 32306                                             	**
**							 		                                    **
**      $Header: $
**								 	                                    **
*************************************************************************/

/*************************************************************************
**									                                    **
** Description								                            **
**									                                    **
** Class SpinSpace defines a composite spin Hilbert space. This will be **
** a vector of individual spin spaces in addition with quantity which contains a count of the	**
** number of spins and their isotope types.  As it is assumed that	**
** complex spin system classes will be derived from SpinSpace, the	**
** following functions have been declared virtual.			**
**									**
**      ~       The SpinSpace destructor					**
**      print   The SpinSpace printing function				**
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

#ifndef _SpinSpace_h_		// Is this file already included?
#define _SpinSpace_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface		// Then this is the interface
#endif

#include <Basics/Isotope.h>        // Include Isotope knowledge
#include <Basics/ParamSet.h>       // Include parameter sets
#include <Matrix/row_vector.h>     // Include row vectors
#include <Matrix/matrix.h>         // Include matrices
#include <string>                  // Include stdlibc++ strings
#include <vector>                  // Include STL vectors
#include <list>                    // Inlcude STL list class
#ifdef _MSC_VER                    // Microsoft VC++ does not handle vectors
  typedef std::vector<int> flagvec;// of bool due to lack of defined < & >
#else                              // operators.  So we used int vectors.
  typedef std::vector<bool> flagvec;
#endif

class SpinSpace
  {
  flagvec spinflags;               // Flags to tag spaces active/inactive
  std::vector<int> *HSs;           // The spin Hilbert spaces
  static int _warn;	               // Flag whether warning messages
  matrix bsmx;                     // A default basis matrix defined by spin HS

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                     CLASS SPIN SPACE ERROR HANDLING
// ____________________________________________________________________________

/*      Input               HSSp    : Hilbert spin space (this)
                            eidx    : Error index
                            noret   : Flag for linefeed (0=linefeed)
                            pname   : string in message
        Output              none    : Error message output
                                      Program execution stopped if fatal     */

void SpinSpace::error(int eidx, int noret=0) const;
void SpinSpace::error(int eidx, const std::string& pname, int noret=0) const;
volatile void SpinSpace::fatality(int eidx) const;

// ____________________________________________________________________________
// ii                   CLASS SPIN SPACE VIRTUAL COPYING
// ____________________________________________________________________________

void virtualCopy()   { references()++; }
void virtualDelete() { references()--; if(references() <= 0) delete mx; }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ____________________________________________________________________________
// A                  SPIN SPACE CONTAINER SUPPORT FUNCTIONS
// ____________________________________________________________________________

///Center Spin Sys Algebraic
///F_list SpinSpace	      - Constructor
///F_list =	              - Assignment

/* These constructors set up a new spin Hilbert space. There are only a couple
   of ways to make a spin Hilbert space:

         Input Arguments                   Resulting System
      ---------------------    -------------------------------------------
               -               Empty spin Hilbert space (HS=0)
              HSSp             Spin Hilbert space duplicate

    The destructor here does nothing except delete the ion index because the
    ions list, although using some small amount of memory, must remain until
    all Cubic Systems in a  program are deleted.  This will not be done
    until program completion.                                                */

        SpinSpace::SpinSpace();
        SpinSpace::SpinSpace(const SpinSpace& HSSp);
virtual SpinSpace::~SpinSpace();
virtual SpinSpace& operator=(const SpinSpace& HSSp);

// ____________________________________________________________________________
// B                         SPIN SPACE COMPARISONS
// ____________________________________________________________________________

// ____________________________________________________________________________
// C                       SPIN SPACE ACCESS FUNCTIONS
// ____________________________________________________________________________

///Center Basic Functions
///F_list size             - Number of spins
///F_list spins	           - Number of spins
///F_list HS	           - Spin or spin system Hilbert space
///F_list qn               - Spin or spin system quantum number

int         SpinSpace::size()             const;
int         SpinSpace::spins()            const;
int         SpinSpace::HS()               const;
int         SpinSpace::HS(int spin)       const;
double      SpinSpace::qn(int spin)       const;
double      SpinSpace::qn()               const;
std::string SpinSpace::momentum(int spin) const;
std::string SpinSpace::momentum()         const;

    // Input		int     : spin label [0, nspins-1]
	// Output		double  : I(spin) in units of hbar
    //                              : default value is 0.5
	// Input		int     : spin label [0, nspins-1]
	// Output		char *  : momentum (e.g. 1/2, 3/2)
    //                                default is 1/2


/*        Functions To Access Spin States & Default Basis Functions

   In GAMMA there is a default basis associated with each spin system.  The
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

      2 spins: aa  aß ba bb        3 spins: aaa aab aba abb aab bab bba bbb

   Note that these are not total Fz ordered, they are ordered according to
   the direct product of indiviual spin basis values!  The functions below
   allow users to obtain information regarding these basis functions with
   respect to indiviual spins and the basis as a whole.                      */


// ____________________________________________________________________________
//                          SPIN SPACE INPUT FUNCTIONS
// ____________________________________________________________________________

///F_list read		- Read spin sys from disk file
///F_list ask_read	- Ask for file, read spin sys from file

// ____________________________________________________________________________
// H                     DEFAULT BASIS FUNCTIONS
// ____________________________________________________________________________

///F_list get_basis	- Spin sys default basis (as a matrix)

matrix get_basis() const;

	// Input		sys      : Spin system
	// Output		bs	 : The default basis defined by
	//				   the spin system Hilbert space
	// Note				 : The basis is output in matrix
	//				   format since actual bases are
	//				   not known to GAMMA on this level

// ______________________________________________________________________
// J                       STANDARD I/O FUNCTIONS
// ______________________________________________________________________

///F_list print			- Send spin system to output stream.
///F_list <<			- Send spin system to an output stream

/* These functions allow users to write spin Hilbert spaces in formatted
   ASCII to an output stream.

 Function  Arguments                      Result
---------- --------- -----------------------------------------------------
  print     ostream  Writes spin Hilbert space in ASCII to output stream
   <<       ostream  Writes spin Hilbert space in ASCII to output stream */


virtual std::ostream& print(std::ostream& out) const;
friend  std::ostream& operator<<(std::ostream& out, const SpinSpace& sys);

// ____________________________________________________________________________
// E                Hilbert Spin Space Container Support Functions
// ____________________________________________________________________________

/* Aside for providing basic tests as to whether two spaces are equivalent or
   not, these operators are necessary if any STL container classes are to
   be used based on Hilbert spin spaces (e.g. list<SpinSpace> or
   vector<SpinSpace>)                                                        */

///F_list ==	              - Equality
///F_list !=	              - Inequality

bool SpinSpace::operator==(const SpinSpace& HSSp) const;
bool SpinSpace::operator!=(const SpinSpace& HSSp) const;
bool SpinSpace::operator< (const SpinSpace& HSSp) const;
bool SpinSpace::operator> (const SpinSpace& HSSp) const;
};

#endif							                                 // SpinSpace.h
