/* relax_WBRexch.h ***********************************************
**								**
** 	                    G A M M A				**
**								**
**    WBR Relaxation & Exchange 		   Interface 	**
**						 		**
**	Copyright (c) 1996			 		**
**	Scott A. Smith					 	**
**      National High Magnetic Field Laboratory			**
**      1800 E. Paul Dirac Drive				**
**      Box 4005						**
**      Tallahassee, FL 32306					**
**								**
**      $Header: $
**								**
*****************************************************************/

/*****************************************************************
**							 	**
** 	Description					 	**
**							 	**
**  This class contains the parameters and functionality that	**
**  may be used control the computation of relaxation and 	**
**  exchange effects in GAMMA.  It is meant to be used along	**
**  with dynamic spin systems when accessing the various	**
**  functions provides in the "relax_*" GAMMA modules.		**
**  						 		**
**  There is nothing provided in the class which cannot be done	**
**  in GAMMA directly, however its use may afford computational	**
**  efficiency as well as code simplicity and flexibility.	**
**						 		**
*****************************************************************/

///Chapter Wangness, Bloch & Redfield Relaxation + Exchange
///Section Overview
///Body The ...
///Section Available WBR & Exchange Conrol Functions

#ifndef   Relax_WBREx_h_		// Is this file already included?
#  define Relax_WBREx_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <string>			// Know stdlibc++ string class
#include <Basics/ParamSet.h>		// Know GAMMA parameter sets
#include <LSLib/sys_dynamic.h>		// Know anisotropic systems 
#include <HSLib/GenOp.h>		// Know general operators
#include <LSLib/SuperOp.h>		// Know superoperators
#include <Level1/SpaceT.h>		// Know spatial tensors
#include <Level1/SpinT.h>		// Know spin tensors
#include <Matrix/matrix.h>		// Know GAMMA matrices
#include <Matrix/complex.h>		// Know GAMMA complex numbers

//
// Added function declaration here before the class definition.
//
// The friend declarations that are included in the class
// are no longer sufficient (in gcc 4.3.3) as function declaration.
// 

// KY -- come back to this later - this function is kiling the build 
// MSVCDLL void ask_relax(int argc, char* argv[], int& argn,
//                               super_op& R, const sys_dynamic& sys, gen_op& H);

MSVCDLL void REXijkl(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS=0);

MSVCDLL void REXijk(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS=0);

MSVCDLL void REXkij(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS=0);

MSVCDLL void REXij(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS=0);
 
MSVCDLL void REXmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w, int hs,
              double* taus, double* c1s, double* c2s, double xi1xi2,
               double w0, double w1, double w2, int l, int level, int autoc, int DFS=0, int het=0);

MSVCDLL void REXmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w, int hs,
              double* taus, double* c1s, double* c2s, double xi1xi2,
         double w0, double w1, double w2, int l, int level, int autoc, int DFS, int het,
                               gen_op& Fz11, double W11, gen_op& Fz12, double W12,
                               gen_op& Fz21, double W21, gen_op& Fz22, double W22);

MSVCDLL void REXrfijkl(super_op &LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS=0);

MSVCDLL void REXrfijk(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS=0);
MSVCDLL void REXrfkij(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS=0);

MSVCDLL void REXrfij(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS=0);

MSVCDLL void REXrfmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, matrix* J12,
                  double* J, double* w, int rank, int level, int autoc, int DFS=0, int het=0);

MSVCDLL void REX_4(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12);

MSVCDLL complex REX_4(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
		 	         int rank, int a, int b, int aa, int bb);

MSVCDLL void REX_4(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12,
                               gen_op& Fz11, double W11, gen_op& Fz12, double W12,
                               gen_op& Fz21, double W21, gen_op& Fz22, double W22);

MSVCDLL complex REX_4(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
                        int rank, int a, int b, int aa, int bb,
                               gen_op& Fz1, double W1, gen_op& Fz2, double W2);

MSVCDLL void REX_3(super_op& LOp, double* w, int rank, gen_op* T1s, gen_op* T2s,
					 matrix& J12, double cutoff=1.e-2);

MSVCDLL void REXrf_4(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix* J12a);

MSVCDLL complex REXrf_4(int hs, gen_op* T1s, gen_op* T2s, matrix* J12,
		 	         int rank, int a, int b, int aa, int bb);

MSVCDLL void REXrf_3(super_op& LOp, double* w, int rank, gen_op* T1s,
                                         gen_op* T2s, matrix* J12, double cutoff=1.e-6);




class WBRExch 
  {
  int DD; 	                      // Flag for Dipolar Relaxation
  int CC; 	                      // Flag for CSA Relaxation
  int QQ; 	                      // Flag for Quadrupolar Relaxation
 
  int DDdfs; 	                      // Flag for Dipolar Relaxation DFS
  int CCdfs; 	                      // Flag for CSA Relaxation DFS
  int QQdfs; 	                      // Flag for Quadrupolar Relaxation DFS
 
  int DC; 	                      // Flag for Dipolar-CSA Cross Correlation
  int DCdfs; 	                      // Flag for Dipolar-CSA Cross Correlation DFS

  int DQ; 	                      // Flag for Dipolar-Quad Cross Correlation
  int DQdfs; 	                      // Flag for Dipolar-Quad Cross Corrlation DFS

  int QC; 	                      // Flag for Quadrupolar-CSA Cross Correlation
  int QCdfs; 	                      // Flag for Quadrupolar-CSA Cross Correlation DFS

  int level; 	                      // Level of Relaxation Computation
  int type; 	                      // Type of Relaxation Computation



// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
 
// ____________________________________________________________________________
// i                RELAXATION & EXCHANGE CONTROL ERROR HANDLING
// ____________________________________________________________________________
 
 
void WBRerror(int eidx, int noret=0) const;
 
        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message
 
 
void WBRerror(int eidx, const std::string& pname, int noret=0) const;
 
        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      eidx    : Flag for error type
        //                      pname   : string included in message
        // Output               none    : Error message
        //                      noret   : Flag for return (0=return)
 
 
volatile void WBRfatality(int eidx) const;

        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      eidx    : Flag for error type
        // Output               none    : Error message
        //                                Program execution stopped
 
 
volatile void WBRfatality(int eidx, const std::string& pname) const;
 
        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      eidx    : Flag for error type
        //                      pname   : string included in message
        // Output               none    : Error message
        //                                Program execution stopped

 
// ____________________________________________________________________________
// ii            RELAXATION & EXCHANGE CONSTRUCTION AUXILIARY FUNCTIONS
// ____________________________________________________________________________


 void assign(const ParameterSet& pset, int DF=1, int CF=1, int QF=1);
 
        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      pset     : A parameter set
        //                      DF       : Flag for possible dipolar relaxation
        //                      CF       : Flag for possible CSA relaxation
        //                      QF       : Flag for possible quadrupolar relaxation
        // Output               none     : WBRE filled with relaxation and 
        //                                 exchange control parameters
        // Note                          : Several parameters are gleaned from
        //                                 the parameter set for a WBRE
        //                                  1.) Rlevel: The level of the computation
        //                                  2.) Rtype:  The type of the computation
        //                                  3.) RDD:    Flag for dipolar relaxation
        //                                  4.) RDDdfs: Flag for dipolar DFS terms
        //                                  5.) RCC:    Flag for CSA relaxation
        //                                  6.) RCCdfs: Flag for CSA DFS terms
        //                                  7.) RQQ:    Flag for Quadrupolar relaxation
        //                                  8.) RQQdfs: Flag for Quad DFS terms
        //                                  9.) RDC:    Flag for Dip-CSA cross correlation
        //                                 10.) RDCdfs: Flag for Dip-CSA DFS terms
        //                                 11.) RDQ:    Flag for Dip-Quad cross correlation
        //                                 12.) RDQdfs: Flag for Dip-Quad DFS terms
        //                                 13.) RQC:    Flag for Quad-CSA cross correlation 
        //                                 14.) RQCdfs: Flag for Quad-CSA DFS terms
        // Note                          : Functions which place a WBRE
        //                                 into a parameter set must contain
        //                                 add the information read here
 


  public:

     
// --------------------------------------------------------------------------------
// ------------------------------ PUBLIC FUNCTIONS --------------------------------
// --------------------------------------------------------------------------------

// ________________________________________________________________________________ 
// A    RELAXATION & EXCHANGE CONTROL CONSTRUCTION, ASSIGNMENT, DESTRUCTION
// ________________________________________________________________________________


MSVCDLC WBRExch();

        // Input                none  :
        // Output               none  : WBR Relaxation & Exchange constructor
        ///F_list spin_system         - Constructor


MSVCDLC WBRExch(const WBRExch& WBRE);
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        // Output               none  : New WBRExch is constructed
        //                              which is identical to WBRE


MSVCDLC ~WBRExch();
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        // Output               none  : WBRExch is deleted


MSVCDLL WBRExch& operator= (const WBRExch &WBRE);

        // Input                WBRE1 : Relaxation/exchange controls (this)
        //                      WBRE  : Relaxation/exchange controls
        // Output               none  : WBRE set identical to WBRE1
        ///F_list =                   - Relaxation/exchange controls assignment
 

// ________________________________________________________________________________
// B             RELAXATION & EXCHANGE CONTROL ACCESS FUNCTIONS
// ________________________________________________________________________________

// --------------------------------------------------------------------------------
//                          Computation Level Access
// --------------------------------------------------------------------------------

 
MSVCDLL void Level(int i);
 
        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for relaxation computation level
        // Output               none  : WBRE level set
 

MSVCDLL int Level() const;
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        // Output               level : Relaxation computation level
 
// --------------------------------------------------------------------------------
//                          Computation Type Access
// --------------------------------------------------------------------------------
 
 
MSVCDLL void Type(int i);
 
        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for relaxation computation type
        // Output               none  : WBRE type set

 
MSVCDLL int Type() const;
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        // Output               type  : Relaxation computation type

   
// --------------------------------------------------------------------------------
//                     Dipolar Relaxation Control Access
// --------------------------------------------------------------------------------
 

MSVCDLL void Dip(int i=1);
 
        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for dipolar relaxation
        // Output               none  : WBRE set to include/neglect
	//				dipole-dipole relaxation effects

MSVCDLL void DipDFS(int i=1);
 
        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for dipolar DFS
        // Output               none  : WBRE set to include/neglect
	//				dipole-dipole dynamic frequency
	//				shift relaxation effects


MSVCDLL void DipCSA(int i=1);
 
        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for dipolar-CSA cross correlation 
        // Output               none  : WBRE set to include/neglect
	//				dipole-CSA cross correlation effects


MSVCDLL void DipCSADFS(int i=1);
 
        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for dipolar-CSA cross correlation 
	//				dynamic frequency shift effects
        // Output               none  : WBRE set to include/neglect
	//				dipole-CSA cross correlation DFS effects


MSVCDLL void DipQuad(int i=1);
 
        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for dipolar-Quad cross correlation 
        // Output               none  : WBRE set to include/neglect
	//				dipole-Quad cross correlation effects


MSVCDLL void DipQuadDFS(int i=1);
 
        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for dipolar-Quad cross correlation 
	//				dynamic frequency shift effects
        // Output               none  : WBRE set to include/neglect
	//				dipole-Quad cross correlation DFS effects

 
// --------------------------------------------------------------------------------
//                     Shift Anisotropy Relaxation Control Access
// --------------------------------------------------------------------------------
 

MSVCDLL void CSA(int i=1);
 
        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for CSA relaxation
        // Output               none  : WBRE set to include/neglect
	//				CSA-CSA relaxation effects


MSVCDLL void CSADFS(int i=1);
 
        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for CSA DFS
        // Output               none  : WBRE set to include/neglect
	//				CSA-CSA dynamic frequency
	//				shift relaxation effects


MSVCDLL void CSADip(int i=1);
 
        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for dipolar-CSA cross correlation 
        // Output               none  : WBRE set to include/neglect
	//				dipole-CSA cross correlation effects


MSVCDLL void CSADipDFS(int i=1);

        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for dipolar-CSA cross correlation 
	//				dynamic frequency shift effects
        // Output               none  : WBRE set to include/neglect
	//				dipole-CSA cross correlation DFS effects



MSVCDLL void CSAQuad(int i=1);

        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for CSA-Quad cross correlation 
        // Output               none  : WBRE set to include/neglect
	//				CSA-Quad cross correlation effects


MSVCDLL void CSAQuadDFS(int i=1);

        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for CSA-Quad cross correlation 
	//				dynamic frequency shift effects
        // Output               none  : WBRE set to include/neglect
	//				CSA-Quad cross correlation DFS effects

 
// --------------------------------------------------------------------------------
//                     Quadrupolar Relaxation Control Access
// --------------------------------------------------------------------------------
 

MSVCDLL void Quad(int i=1);
 
        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for quadrupolar relaxation
        // Output               none  : WBRE set to include/neglect
	//				quadrupole-quadrupole relaxation effects


MSVCDLL void QuadDFS(int i=1);
 
        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for quadrupolar DFS
        // Output               none  : WBRE set to include/neglect
	//				quadrupole-quadrupole dynamic frequency
	//				shift relaxation effects


MSVCDLL void QuadDip(int i=1);

        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for dipolar-Quad cross correlation 
        // Output               none  : WBRE set to include/neglect
	//				dipole-Quad cross correlation effects


MSVCDLL void QuadDipDFSQuad(int i=1);
 
        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for dipolar-Quad cross correlation 
	//				dynamic frequency shift effects
        // Output               none  : WBRE set to include/neglect
	//				dipole-Quad cross correlation DFS effects


MSVCDLL void QuadCSA(int i=1);

        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for CSA-Quad cross correlation 
        // Output               none  : WBRE set to include/neglect
	//				CSA-Quad cross correlation effects


MSVCDLL void QuadCSADFS(int i=1);

        // Input                WBRE  : Relaxation/exchange controls (this)
	//			i     : Flag for CSA-Quad cross correlation 
	//				dynamic frequency shift effects
        // Output               none  : WBRE set to include/neglect
	//				CSA-Quad cross correlation DFS effects

 
// --------------------------------------------------------------------------------
//                         Mutual Exchange Control Access
// --------------------------------------------------------------------------------

// --------------------------------------------------------------------------------
//                       Non-Mutual Exchange Control Access
// --------------------------------------------------------------------------------


// ________________________________________________________________________________
// C             RELAXATION & EXCHANGE PARAMETER AUXILIARY FUNCTIONS
// ________________________________________________________________________________


MSVCDLL double LWhh(const sys_dynamic& sys, const std::string& Iso);
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      sys   : Dynamic spin system
        //                      Iso   : A spin isotope label
        // Output               LW    : Maximum linewidth at half-height (Hz)
        //                              estimated under simple assumptions
        // Note			      : Here, the motion is assumed to be diffusive and
        //                           	the system moving as an isotropic top.


// ________________________________________________________________________________
// D             RELAXATION & EXCHANGE PARAMETER SET FUNCTIONS
// ________________________________________________________________________________

/*  Function/Operator                         Result
    -----------------   ------------------------------------------------------
    ParameterSet        Convert parameters into a parameter set
       +=               Add parameters to a parameter set
    SetZero             Set all flags to zero
    SetLevel            Set level from parameter Rlevel in pset (default 4)
    SetType             Set type from parameter Rtype in pset   (default 0)
    SetDip              Set dipolar flags from RDD & RDDdfs in pset
    SetSA               Set shift anisotropy flags from RCC & RCCdfs in pset
    SetQuad             Set quadrupolar flags from RQQ & RQQdfs in pset
    SetQuad             Set quadrupolar flags from RQQ & RQQdfs in pset
    SetQuad             Set quadrupolar flags from RQQ & RQQdfs in pset
    SetQuad             Set quadrupolar flags from RQQ & RQQdfs in pset
    SetQuad             Set quadrupolar flags from RQQ & RQQdfs in pset
 
    Again, note that operations which return TF will trigger class matrix to
    take appropriate action when they return FALSE.  This usually means that
    the operation cannot be handled by a Hermitian array so the array type
    is changed to be more generic.  The the operation is reattempted.       */
 
MSVCDLL             operator ParameterSet( ) const;
MSVCDLL friend void operator+= (ParameterSet& pset, const WBRExch& WBRE);

MSVCDLL void SetZero();
MSVCDLL void SetLevel(const ParameterSet& pset);
MSVCDLL void SetType(const  ParameterSet& pset);
MSVCDLL void SetDip(const   ParameterSet& pset);
MSVCDLL void SetSA(const    ParameterSet& pset);
MSVCDLL void SetQuad(const  ParameterSet& pset);
MSVCDLL void SetDCX(const   ParameterSet& pset);

        // Input                WBRE    : Relaxation/exchange controls (this)
        //                      pset    : A parameter set
        // Output               none    : WBRE dipolar-SA cross correlation
	//				 flags are set from parameters in pset
        // Note                         : Parameters are as follows:
        //                                RDC    - dipolar-SA cross correlation flag
        //                                RDCdfs - dipolar-SA cross correlation  DFS flag


MSVCDLL void SetDQX(const ParameterSet& pset);

        // Input                WBRE    : Relaxation/exchange controls (this)
        //                      pset    : A parameter set
        // Output               none    : WBRE dipolar-quad cross correlation flags
        //                                are set from parameters in pset
        // Note                         : Parameters are as follows:
        //                                RDQ    - dipolar-quad cross correlation flag
        //                                RDQdfs - dipolar-quad cross correlation  DFS flag
        // Note                         : The dipolar and quad flags should be
        //                                set before this is called


MSVCDLL void SetQCX(const ParameterSet& pset);

        // Input                WBRE    : Relaxation/exchange controls (this)
        //                      pset    : A parameter set
        // Output               none    : WBRE quad-SA cross correlation flags
        //                                are set from parameters in pset
        // Note                         : Parameters are as follows:
        //                                RQC    - quad-SA cross correlation flag
        //                                RQCdfs - quad -SA cross correlation  DFS flag
        // Note                         : The quad and SA flags should be
        //                                set before this is called
 
 
MSVCDLL WBRExch& operator= (const ParameterSet& pset);
 
        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      pset     : A parameter set
        // Output               none     : WBRE filled with relaxation and
        //                                 exchange control parameters
        // Note                          : Three things are gleaned from
        //                                 the parameter set for a WBRE
        //                                 1.) The level of the computation
        //                                 2.) The type of the computation
        //                                 3.) Flag for dipolar relaxation
        // Note                          : Functions which place a WBRE
        //                                 into a parameter set must contain
        //                                 add the information read here



// --------------------------------------------------------------------------------
//                    Prepare Quadrupolar Relaxation Components
//                   (Only If Quadrupolar Effects Are Included)
// --------------------------------------------------------------------------------

//void prepCSA(const sys_dynamic& sys, matrix& Xis, spin_T* Ts, space_T* As)

        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      sys   : Dynamic spin system
        //                      Xis   : Matrix for CSA interacation constants
        //                      Ts    : NULL array for CSA spin tensor compts.
        //                      As    : NULL array for CSA space tensor compts.
        // Output               void  : WBR Relaxation and Exchange
        //                              superoperator

MSVCDLL void prepQuad(const sys_dynamic& sys, matrix& Xis, spin_T* Ts, space_T* As) const;

        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      sys   : Dynamic spin system
        //                      Xis   : Matrix for Quad interacation constants
        //                      Ts    : NULL array for Quad spin tensor compts.
        //                      As    : NULL array for Quad space tensor compts.
        // Output               void  : WBR Relaxation and Exchange
        //                              superoperator

// ________________________________________________________________________________
// E                     RELAXATION & EXCHANGE OUTPUT FUNCTIONS
// ________________________________________________________________________________

 
MSVCDLL std::ostream& printDip(std::ostream& ostr) const;
 
        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      ostr     : An output stream
        // Output               ostr     : The output stream modified by
        //                                 the WBRE Dipolar parameters


MSVCDLL std::ostream& printSA(std::ostream& ostr) const;
 
        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      ostr     : An output stream
        // Output               ostr     : The output stream modified by
        //                                 the WBRE Shift Anisotropy parameters
 

MSVCDLL std::ostream& print(std::ostream& ostr) const;

        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      ostr     : An output stream
        // Output               ostr     : The output stream modified by
        //                                 the information in WBRE
 

MSVCDLL std::ostream& printQ(std::ostream& ostr) const;

        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      ostr     : An output stream
        // Output               ostr     : The output stream modified by
        //                                 the WBRE Quadrupolar parameters

 
MSVCDLL friend std::ostream &operator << (std::ostream &ostr, const WBRExch& WBRE);
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      ostr  : Output stream
        // Output               none  : WBRE info is sent to the output stream


// ________________________________________________________________________________
// E                   RELAXATION & EXCHANGE INPUT FUNCTIONS
// ________________________________________________________________________________


MSVCDLL void read(const std::string& filename);

        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      filename : Input filename
        // Output               none     : WBRE parameters filled with
        //                                 parameters read from file
        // Note                          : The file should be an ASCII file
        //                                 containing recognized WBR parameters
 
 
MSVCDLL void read(const std::string& filename, const sys_dynamic& sys);

        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      filename : Input filename
        //                      sys      : A dynamic spin system
        // Output               none     : WBRE parameters filled with
        //                                 parameters read from file
        //                                 under the constraints set by the     
        //                                 input spin system
        // Note                          : The file should be an ASCII file
        //                                 containing recognized WBR parameters


MSVCDLL void ask_read(int argc, char* argv[], int argn);

        // Input                WBRE    : Relaxation/exchange controls (this)
        //                      argc    : Number of arguments
        //                      argv    : Vector of argc arguments
        //                      argn    : Argument index
        // Output               void    : The parameter argn of array argc
        //                                is used to supply a filename
        //                                from which WBRE parameters are read
        //                                If the argument argn is not in argv,
        //                                the user is asked to supply a filename
        // Note                         : The file should be an ASCII file
        //                                containing recognized WBRE parameters
        // Note                         : The WBRE is modifed (filled)


MSVCDLL void ask_read(int argc, char* argv[], int argn, const sys_dynamic& sys);
 
        // Input                WBRE    : Relaxation/exchange controls (this)
        //                      argc    : Number of arguments
        //                      argv    : Vector of argc arguments
        //                      argn    : Argument index
        //                      sys     : A dynamic spin system
        // Output               void    : The parameter argn of array argc
        //                                is used to supply a filename
        //                                from which WBRE parameters are read
        //                                under the constraints set by the      
        //                                input spin system
        //                                If the argument argn is not in argv,
        //                                the user is asked to supply a filename
        // Note                         : The file should be an ASCII file
        //                                containing recognized WBRE parameters
        // Note                         : The WBRE is modifed (filled)
 

// ________________________________________________________________________________
// F              INTERACTIVE RELAXATION & EXCHANGE SETUP FUNCTIONS
// ________________________________________________________________________________


MSVCDLL void ask(int argc, char* argv[], int& argn);

        // Input                sys      : Spin system (this)
        //                      ostr     : An output stream
        // Output               ostr     : The output stream modified by


// ________________________________________________________________________________
// G              RELAXATION & EXCHANGE SUPEROPERATOR PRODUCTION
// ________________________________________________________________________________


MSVCDLL super_op REX(const sys_dynamic& sys, gen_op& Ho, int fext=0);
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      sys   : Dynamic spin system
        //                      Ho    : Static Hamiltonian
 	//                      fext  : Signals use of external functions
        // Output               LOp   : WBR Relaxation and Exchange
        //                              superoperator
        // Note                         Computed in the eigenbasis of Ho
 
 
MSVCDLL super_op REXrf(const sys_dynamic& sys, gen_op& Heff, double Wrf, int fext=0);
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      sys   : Dynamic spin system
        //                      Heff  : Effective static Hamiltonian
	//			Wrf   : RF field offset
 	//                      fext  : Signals use of external functions
        // Output               LOp   : WBR Relaxation and Exchange
        //                              superoperator
        // Note                         Computed in the eigenbasis of Heff
 

// ______________________________________________________________________
// ************************ Class WBRExch LeftOvers *********************
// ______________________________________________________________________


// KY - this function is causing problems in the build...
// MSVCDLL friend void ask_relax(int argc, char* argv[], int& argn,
//                               super_op& R, const sys_dynamic& sys, gen_op& H);

        // Input        argc    : Number of command line arguments
        //              argv    : Command line arguments
        //              argn    : Initial command line argument
        //                        for relaxation parameters
        //              R       : Relaxation superoperator
        //              sys     : Dynamic spin system
        //              H       : Hamiltonian
        // Output       none    : Function is void.  The relaxation
        //                        matrix R has different effects added
        //                        in depending upon user requests


//*************************************************************
//*************************************************************
//*************************************************************
//*************************************************************
 
MSVCDLL friend void REXijkl(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS);

MSVCDLL friend void REXijk(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS);

MSVCDLL friend void REXkij(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS);

MSVCDLL friend void REXij(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS);
 
MSVCDLL friend void REXmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w, int hs,
              double* taus, double* c1s, double* c2s, double xi1xi2,
               double w0, double w1, double w2, int l, int level, int autoc, int DFS, int het);

// ---------------- Routing to Specific Level Computation ---------------
// ------------------ This Is For Heteronuclear Systems -----------------

MSVCDLL friend void REXmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w, int hs,
              double* taus, double* c1s, double* c2s, double xi1xi2,
         double w0, double w1, double w2, int l, int level, int autoc, int DFS, int het,
                               gen_op& Fz11, double W11, gen_op& Fz12, double W12,
                               gen_op& Fz21, double W21, gen_op& Fz22, double W22);

 


//*************************************************************
//*************************************************************
//*************************************************************
//*************************************************************

MSVCDLL friend void REXrfijkl(super_op &LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS);


 
// --------------- Spin-Pair with Spin Functions -------------------

MSVCDLL friend void REXrfijk(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS);


MSVCDLL friend void REXrfkij(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS);



// -------------------- Spin with Spin Functions -------------------

MSVCDLL friend void REXrfij(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS);


MSVCDLL friend void REXrfmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, matrix* J12,
                  double* J, double* w, int rank, int level, int autoc, int DFS, int het);



//*************************************************************
//*************************************************************
//*************************************************************
//*************************************************************



// ---------------- Level 4 via Element Calculation ---------------------

MSVCDLL friend void REX_4(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12);

MSVCDLL friend complex REX_4(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
		 	         int rank, int a, int b, int aa, int bb);


// ---------------- Level 4 via Element Calculation ---------------------
// ------------ Multiple Rotating Frame (Heteronuclear) -----------------

MSVCDLL friend void REX_4(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12,
                               gen_op& Fz11, double W11, gen_op& Fz12, double W12,
                               gen_op& Fz21, double W21, gen_op& Fz22, double W22);
//friend void REX_4(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12,
//                               gen_op& Fz1, double W1, gen_op& Fz2, double W2);

 
MSVCDLL friend complex REX_4(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
                        int rank, int a, int b, int aa, int bb,
                               gen_op& Fz1, double W1, gen_op& Fz2, double W2);
 

// ---------------- Level 3 via Element Calculation ---------------------
//                  (Applies Secular Approximation)

MSVCDLL friend void REX_3(super_op& LOp, double* w, int rank, gen_op* T1s, gen_op* T2s,
					 matrix& J12, double cutoff);


//*************************************************************
//*************************************************************
//*************************************************************
//*************************************************************

// ----------------- Level 4 via Element Calculation --------------------

MSVCDLL friend void REXrf_4(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix* J12a);

MSVCDLL friend complex REXrf_4(int hs, gen_op* T1s, gen_op* T2s, matrix* J12,
		 	         int rank, int a, int b, int aa, int bb);

// ---------------- Level 3 via Element Calculation ---------------------
//                  (Applies Secular Approximation)

MSVCDLL friend void REXrf_3(super_op& LOp, double* w, int rank, gen_op* T1s,
                                         gen_op* T2s, matrix* J12, double cutoff);


  };


#endif 						/* __RELAX_WBREx_H__ */
