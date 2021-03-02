/* BlochSys.h ***************************************************-*-c++-*-
**									**
**	                        G A M M A				**
**								 	**
**	Bloch Spin System                            Interface		**
**						 			**
**      Copyright (c) 2001                                              **
**      Dr. Scott A. Smith						**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**						 			**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**								 	**
** The class BlochSys defines a spin system used as a basis for 	**
** simulations based on the phenomenological Bloch equations. The	**
** Bloch equations account for spin magnetization evolution under	**
** a static Bo field, an applied rf-field (B1), simplistic relaxation	**
** (T1 & T2), and possibly exchange between magnetizaiton vectors.	**
**									**
** Bloch spin system tracks any number of spins and their associated	**
** magnetization vectors. To each spin, and all of its associated	**
** magnetizationv vectors, it assigns single longitudinal & transverse	**
** relaxation rates. Between any two spins (and their magnetization	**
** vectors) it assigns an exchange rate.				**
**									**
** Based on this Bloch spin system, functions can easily built that	**
** generate magnetization vectors, relaxation matrices, exchange	**
** matrices, and evolution matrices with/without an rf-field present.	**
**									**
*************************************************************************/

///Chapter Class Spin System
///Section Overview
///Body    None
///Section Available Spin System Functions

#ifndef   BlochSys_h_			// Is file already included?
#  define BlochSys_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Basics/ParamSet.h>		// Includes parameter sets
#include <Basics/Gconstants.h>		// Include GAMMA1H definition
#include <Basics/Gutils.h>		// Include GAMMA errors/queries
#include <Basics/Isotope.h>		// Include GAMMA isotopes
#include <HSLib/SpinSystem.h>		// Include isotropic systems
#include <Level1/coord.h>		// Include GAMMA coordinates
#include <Level1/coord_vec.h>		// Include GAMMA cordinate vectors
#include <Level2/RelaxBas.h>		// Include basic relaxation
#include <Level2/TrnsTable1D.h>		// Include transition tables
#include <Bloch/MagVec.h>		// Include magnetization vectors
#include <string>			// Include libstdc++ strings
#include <vector>			// Include libstdc++ STL vectors

class BlochSys
  {
  std::vector <double>  Offsets;	// Mag. vector offsets
  std::vector <Isotope> isotopes;	// Mag. vector associated isotope
  std::vector <double>  R1rates;	// Mag. vector R1 rates
  std::vector <double>  R2rates;	// Mag. vector R2 rates
  std::vector <double>  Krates;		// Mag. vector exchange rates
  std::vector <int>     Spins;		// Mag. vector associated spin
  MagVec                _M;		// Magnetization vector
  
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                    Bloch System Error Handling
// ____________________________________________________________________________


        // Input                sys     : Spin system (this)
        //                      ei	: Error index
        //                      nr	: Flag for linefeed (0=linefeed)
        //                      pn	: string in message
     
         void BSerror(int ei,                        int nr=0) const;
         void BSerror(int ei, const std::string& pn, int nr=0) const;
volatile void BSfatal(int ei) const;
volatile void BSfatal(int ei, const std::string& pn) const;

// ____________________________________________________________________________
// ii                   Bloch System Checking Functions
// ____________________________________________________________________________

bool CheckR1s(const    std::vector<double>& R1s, bool warn=true) const;
bool CheckR2s(const    std::vector<double>& R2s, bool warn=true) const;
bool CheckIsos(const   std::vector<Isotope>& Is, bool warn=true) const;
bool CheckKs(const     std::vector<double>&  Ks, bool warn=true) const;
bool CheckSpins(int    ns1, int ns2,        bool warn=true) const;
bool CheckNorms(const  std::vector<double>&  Ns, bool warn=true) const;
bool CheckCoords(const coord_vec&       Ms, bool warn=true) const;

// ____________________________________________________________________________
// iii                     Bloch System Setup Functions
// ____________________________________________________________________________

/* These are protected functions because they allow specific aspects of the
   spin system to be set up without worrying about system consistency!       */

bool SetSystem(const ParameterSet& pset, int idx=-1, bool warn=true);

bool GetNSpins(const ParameterSet& pset,int& ns,bool wn=true) const;
bool GetNVects(const ParameterSet& pset,int& nm,bool wn=true) const;
bool SetVects(const  ParameterSet& pset,int   N,bool wn=true);

// ----------------------------------------------------------------------------
//         These Read All Single Magnetization Vector Parameters
// ----------------------------------------------------------------------------


bool GetVect(const ParameterSet& pset, int i, double& v, Isotope& I,
                        double& R1, double& R2, int& Sp, bool warn=true) const;

bool GetW(const   ParameterSet& pset, int i, 
                                              double& v, bool warn=true) const;
bool GetIso(const ParameterSet& pset, int i,
                                             Isotope& I, bool warn=true) const;
bool GetR1(const  ParameterSet& pset, int i, 
                                             double& R1, bool warn=true) const;
bool GetR2(const  ParameterSet& pset, int i, 
                                             double& R2, bool warn=true) const;
bool GetSp(const  ParameterSet& pset, int i,
                                                int& Sp, bool warn=true) const;

// ----------------------------------------------------------------------------
//         These Read All Magnetization Vector Pair Parameters
// ----------------------------------------------------------------------------

bool SetExchange(const ParameterSet& pset, int nm, bool warn=true);

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ____________________________________________________________________________
// A                   SPIN SYSTEM CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

///Center Spin System Algebraic
///F_list BlochSys		- Constructor
///F_list ~			- Destructor
///F_list =			- Assignment	

// ----------------------------------------------------------------------------
//                          Simple Constructors
// ----------------------------------------------------------------------------

MSVCDLC BlochSys(int spins=0);
MSVCDLC BlochSys(const BlochSys& sys);

// ----------------------------------------------------------------------------
//                   Constructors Using Single Vector Arguments
// ----------------------------------------------------------------------------

MSVCDLC BlochSys(double w, double R1, double R2);

// ----------------------------------------------------------------------------
//                  Constructors Using Multiple Vector Arguments
// ----------------------------------------------------------------------------

MSVCDLC BlochSys(const std::vector<double>& SH, 
               const std::vector<double>& R1s, const std::vector<double>& R2s);

MSVCDLC BlochSys(const std::vector<double>& SH, const std::vector<Isotope>& Is,
               const std::vector<double>& R1s, const std::vector<double>& R2s);

MSVCDLC BlochSys(const std::vector<double>& SH, 
               const std::vector<double>& R1s, const std::vector<double>& R2s,
                                                const std::vector<double>& Ks);

MSVCDLC BlochSys(const std::vector<double>& SH, const std::vector<Isotope>& Is,
               const std::vector<double>& R1s, const std::vector<double>& R2s,
                                                const std::vector<double>& Ks);

// ----------------------------------------------------------------------------
//                     Constructors Using Spin Systems
// ----------------------------------------------------------------------------

MSVCDLC BlochSys(const spin_system& sys, const RBasic& Rs);
MSVCDLC BlochSys(const spin_system& sys, const matrix& Ks);
MSVCDLC BlochSys(const spin_system& sys, const RBasic& Rs, const matrix& Ks);

// ----------------------------------------------------------------------------
//                     Constructors Using Other Objects
// ----------------------------------------------------------------------------

MSVCDLC BlochSys(const TTable1D& TT, const std::string& Iso=DEFISO);

// ----------------------------------------------------------------------------
//                       Assignment and Destruction
// ----------------------------------------------------------------------------

MSVCDLC ~BlochSys();
MSVCDLL BlochSys& operator= (const BlochSys& sys);

// ____________________________________________________________________________
// B                 Magnetization Vector Offset Access
// ____________________________________________________________________________

// ____________________________________________________________________________
// C                Magnetization Vector Isotope Access
// ____________________________________________________________________________

MSVCDLL int NIso()         const;
MSVCDLL int IsoMaxLength() const;

MSVCDLL int NSpins()       const;

// ____________________________________________________________________________
// D                Magnetization Vector Exchange Rate Access
// ____________________________________________________________________________

MSVCDLL double R1(int i) const;
MSVCDLL double T1(int i) const;
MSVCDLL double R2(int i) const;
MSVCDLL double T2(int i) const;
MSVCDLL double LW(int i) const;

MSVCDLL double MaxExchange() const;

// ____________________________________________________________________________
// F                Magnetization Vector Component Access
// ____________________________________________________________________________

MSVCDLL std::vector<double> Norms() const;
MSVCDLL void                Norms(const std::vector<double>& Ns);

MSVCDLL double Norm(int i) const;
MSVCDLL void   Norm(double nv, int i);

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

MSVCDLL matrix H()                                       const;
MSVCDLL matrix H(double gamB1, double w=0, double phi=0) const;

MSVCDLL matrix B()                                       const;
MSVCDLL matrix B(double gamB1, double w=0, double phi=0) const;
MSVCDLL matrix R()                                       const;
MSVCDLL matrix K()                                       const;

// ____________________________________________________________________________
// H                     Bloch Equation Column Vectors
// ____________________________________________________________________________

MSVCDLL       MagVec  Meq() const;
MSVCDLL const MagVec& Mo() const;
MSVCDLL       MagVec  Mx() const;
MSVCDLL       MagVec  My() const;
MSVCDLL       MagVec  Mz() const;
MSVCDLL       MagVec  Mss(const matrix& L, const matrix& R) const;
MSVCDLL       MagVec  Mss(const matrix& L, const matrix& R,
                                                  const col_vector& Meq) const;

// ____________________________________________________________________________
// I                       Bloch Equation Row Vectors
// ____________________________________________________________________________

MSVCDLL row_vector DetectMu(                 int u) const;
MSVCDLL row_vector DetectMu(int k,           int u) const;
MSVCDLL row_vector DetectMu(const std::string& I, int u) const;

MSVCDLL row_vector DetectMx() const;
MSVCDLL row_vector DetectMy() const;
MSVCDLL row_vector DetectMz() const;

MSVCDLL row_vector DetectMx(int i) const;
MSVCDLL row_vector DetectMy(int i) const;
MSVCDLL row_vector DetectMz(int i) const;

MSVCDLL row_vector DetectMx(const std::string& I) const;
MSVCDLL row_vector DetectMy(const std::string& I) const;
MSVCDLL row_vector DetectMz(const std::string& I) const;

// ____________________________________________________________________________
// J                    Bloch System Auxiliary Functions
// ____________________________________________________________________________

MSVCDLL int size() const;

// ____________________________________________________________________________
// K                  Bloch System Parameter Set Functions
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//            Functions To Make A Spin System From A Parameter Set
// ----------------------------------------------------------------------------


/*
virtual int write(const std::string &filename, int idx=-1, int warn=2) const;
virtual int write(std::ofstream& ofstr, int idx=-1, int warn=2) const; 
*/

// ____________________________________________________________________________
// L                     Bloch System Input Functions
// ____________________________________________________________________________

///F_list read		- Read Bloch system from disk file
///F_list ask_read	- Ask for file, read Bloch system from file

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


MSVCDLL bool read(const std::string& fn,    int idx=-1, int warn=2);
MSVCDLL bool read(const ParameterSet& pset, int idx=-1, int warn=2);

MSVCDLL std::string ask_read(int argc, char* argv[], int argn);
MSVCDLL std::string ask_read(int argc, char* argv[], int argn,
                                                       const std::string& def); 

// ____________________________________________________________________________
// M                   Bloch System Standard I/O Functions
// ____________________________________________________________________________

MSVCDLL std::ostream& print(std::ostream& out) const;
MSVCDLL friend std::ostream& operator<<(std::ostream& out, const BlochSys& sys);

};

#endif						// BlochSys.h
