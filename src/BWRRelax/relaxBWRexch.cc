/* relaxWBRexch.cc ******************************************************
**									**
**                               G A M M A				**
**									**
**    WBR Relaxation & Exchange           	   Implementation	**
**									**
**      Copyright (c) 1996						**
**      Scott A. Smith							**
**      National High Magnetic Field Laboratory				**
**      1800 E. Paul Dirac Drive					**
**      Box 4005							**
**      Tallahassee, FL 32306						**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
** Description							 	**
**									**
**  This class contains the parameters and functionality that may be	**
**  used to control the computation of relaxation and exchange effects	**
**  in GAMMA.  It is meant to be used along with dynamic spin systems	**
**  when accessing the various functions provided in the "relax_*"	**
**  GAMMA modules.							**
**									**
**  There is nothing provided in the class which cannot be done	in	**
**  GAMMA directly, however its use may afford computational efficiency	**
**  as well as code simplicity and flexibility.				**
**									**
*************************************************************************/

#ifndef _relax_WBREx_cc_		// Is this file already included?
#define _relax_WBREx_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation          	// This is the implementation
#endif

#include <BWRRelax/relaxBWRexch.h>	// Include the header
#include <Basics/ParamSet.h>		// Include parameter sets
#include <stdlib.h>
#include <string>			// Include strings
#include <LSLib/sys_dynamic.h>		// Include dynamic spin systems 
#include <HSLib/GenOp.h>		// Include operators 
#include <LSLib/SuperOp.h>		// Include superoperators
#include <BWRRelax/relaxNMR.h>		// Knowledge of relaxtion routines
#include <HSLib/SpinOpCmp.h>		// Knowledge of spin operators
#include <BWRRelax/relaxanalyze.h>	// Knowledge of analysis routines
#include <Basics/Gconstants.h>		// Know GAMMA1H value
#include <Level1/nmr_tensor.h>		// Know common MR spin tensors


                                                                                
// ________________________________________________________________________________
// i                RELAXATION & EXCHANGE CONTROL ERROR HANDLING
// ________________________________________________________________________________

                                                                                
void WBRExch::WBRerror(int eidx, int noret) const

        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message

  {
  std::cout<< "\nClass WBRExch: ";
  switch(eidx)
    {
    case 0:								// (0)
      std::cout << "Program Aborting.....";
      break;
    case 9:                                                             // (9)
       std::cout << "Setting Relaxation Computation Level To 4";
      break;
    case 10:                                                             // (10)
       std::cout << "Setting Relaxation Computation Type To 0";
      break;
    case 11:                                                             // (11)
      std::cout << "Setting Dipolar Dynamic Frequency Shift Flag To 1";
      break;
    case 12:                                                             // (12)
      std::cout << "Setting Dipolar Relaxation Flag To 1";
      break;
    case 13:                                                             // (13)
      std::cout << "Setting Shift Anisotropy Dynamic Frequency Shift Flag To 1";
      break;
    case 14:                                                             // (14)
      std::cout << "Setting Shift Anisotropy Relaxation Flag To 1";
      break;
    case 15:                                                             // (15)
      std::cout << "Setting Quadrupolar Dynamic Frequency Shift Flag To 1";
      break;
    case 16:                                                             // (16)
      std::cout << "Setting Quadrupolar Relaxation Flag To 1";
      break;
    case 17:                                                             // (17)
      std::cout << "Setting Dipole-SA Cross-Correlation Dynamic Frequency Shift Flag To 1";
      break;
    case 18:                                                             // (18)
      std::cout << "Setting Dipole-SA Cross-Correlation Flag To 1";
      break;
    case 19:                                                             // (19)
      std::cout << "Setting Dip-Quad Cross-Correlation Dynamic Frequency Shift Flag To 1";
      break;
    case 20:                                                             // (20)
      std::cout << "Setting Dipole-Quadrupole Cross-Correlation Flag To 1";
      break;
    case 21:                                                             // (21)
      std::cout << "Setting SA-Quad Cross-Correlation Dynamic Frequency Shift Flag To 1";
      break;
    case 22:                                                             // (22)
      std::cout << "Setting Shift Anisotropy - Quadrupole Cross-Correlation Flag To 1";
      break;
    case 31:                                                             // (31)
      std::cout << "Dipolar DFS Effects Disallowed Without Dip-Dip Relaxation";
      break;
    case 32:                                                             // (32)
      std::cout << "Dipole-CSA Cross-Corr. Disallowed Without Dip-Dip Relaxation";
      break;
    case 33:                                                             // (33)
      std::cout << "Dipole-CSA Cross-Corr. Disallowed Without CSA-CSA Relaxation";
      break;
    case 34:                                                             // (34) 
      std::cout << "Dipole-CSA Cross-Corr. DFS Disallowed If No Dip-Dip Relaxation"; 
      break; 
    case 35:                                                             // (35)
      std::cout << "Dipole-CSA Cross-Corr. DFS Disallowed If No CSA-CSA Relaxation"; 
      break;
    case 36:                                                             // (36)
      std::cout << "Dipole-CSA Cross-Corr. DFS Disallowed If No Dipole-CSA Terms"; 
      break;
    case 37:                                                             // (37)
      std::cout << "Dipole-Quad Cross-Corr. Disallowed If No Dip-Dip Relaxation";
      break;
    case 38:                                                             // (38)
      std::cout << "Dipole-Quad Cross-Corr. Disallowed If No Quad-Quad Relaxation";
      break;
    case 39:                                                             // (39)
      std::cout << "Dipole-Quad Cross-Corr. DFS Disallowed If No Dip-Dip Relaxation";
      break;
    case 40:                                                             // (40)
      std::cout << "Dipole-Quad Cross-Corr. DFS Disallowed If No Quad-Quad Relaxation";
      break;
    case 41:                                                             // (41)
      std::cout << "Dipole-Quad Cross-Corr. DFS Disallowed If No Dipole-Quad Terms";
      break;
    case 42:                                                             // (42)
      std::cout << "Shift Anisotropy DFS Effects Disallowed If No CSA-CSA Relaxation";
      break;
    case 43:                                                             // (43)
      std::cout << "Quadrupole-CSA Cross-Corr. Disallowed If No Quad-Quad Relaxation";
      break;
    case 44:                                                             // (44)
      std::cout << "Quadrupole-CSA Cross-Corr. Disallowed If No CSA-CSA Relaxation";
      break;
    case 45:                                                             // (45)
      std::cout << "CSA-Quad Cross-Corr. DFS Disallowed If No Quad-Quad Relaxation";
      break;
    case 46:                                                             // (46)
      std::cout << "CSA-Quad Cross-Corr. DFS Disallowed If No CSA-CSA Relaxation";
      break;
    case 47:                                                             // (47)
      std::cout << "CSA-Quad Cross-Corr. DFS Disallowed If No CSA-Quad Cross Terms";
      break;
    case 48:                                                             // (48)
      std::cout << "Quadrupolar DFS Effects Disallowed If No Quad-Quad Relaxation";
      break;
    case 50:                                                             // (50)
      std::cout << "Problems Constructing WBR Relaxation & Exchange Controls\n";
      break;
    case 60:                                                             // (60)
       std::cout << "\n\tWarning: Unable to Obtain Proper Ho(lab) Eigenbasis!";
      break;
    default:
      std::cout<<"Unknown error";
    }
  if(!noret) std::cout << ".\n";
  }  


void WBRExch::WBRerror(int eidx, const std::string& pname, int noret) const

        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      eidx    : Flag for error type
        //                      pname   : string included in message
        // Output               none    : Error message

  {
  std::cout<< "\nClass WBRExch: ";
  switch(eidx)
    {
    case 1:                                                             // (1)
      std::cout << "Problems with File " << pname;
      break;
    case 2:                                                             // (2)
      std::cout << "Cannot Read Parameter " << pname << " From Parameter Set or File";
      break;
     }
  if(!noret) std::cout << ".\n";
  }  


volatile void WBRExch::WBRfatality(int eidx) const

        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      eidx    : Flag for error type
        // Output               none    : Error message
        //                                Program execution stopped

  {                                                                 
  WBRerror(eidx);			// Output the error message
  if(eidx) WBRerror(0);			// State that its fatal
  GAMMAfatal();				// Clean exit from program
  }


volatile void WBRExch::WBRfatality(int eidx, const std::string& pname) const

        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      eidx    : Flag for error type
        //                      pname   : string included in message
        // Output               none    : Error message
        //                                Program execution stopped

  {                                                                 
  WBRerror(eidx, pname, 1);		// Output the error message
  if(eidx) WBRerror(0);			// State that its fatal
  GAMMAfatal();				// Clean exit from program
  }

 
// ____________________________________________________________________________
// ii            RELAXATION & EXCHANGE CONSTRUCTION AUXILIARY FUNCTIONS
// ____________________________________________________________________________
 
      
void WBRExch::assign(const ParameterSet& pset, int DF, int CF, int QF)
 
        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      pset	 : A parameter set
	//			DF       : Flag for possible dipolar relaxation
	//			CF       : Flag for possible CSA relaxation
	//			QF       : Flag for possible quad. relaxation
        // Output               none     : WBRE filled with relaxation and 
	//				   exchange control parameters
        // Note                          : Several parameters are gleaned from
        //                                 the parameter set for a WBRE
        //                                  1.) Rlevel: Level of computation
        //                                  2.) Rtype:  Type of computation
        //                                  3.) RDD:    Flag dip. relaxation
        //                                  4.) RDDdfs: Flag dipolar DFS
        //                                  5.) RCC:    Flag CSA relaxation
        //                                  6.) RCCdfs: Flag CSA DFS terms
        //                                  7.) RQQ:    Flag Quad. relaxation
        //                                  8.) RQQdfs: Flag Quad DFS terms
        //                                  9.) RDC:    Flag Dip-CSA x-corr.
        //                                 10.) RDCdfs: Flag Dip-CSA DFS
        //                                 11.) RDQ:    Flag Dip-Quad x-corr.
        //                                 12.) RDQdfs: Flag Dip-Quad DFS
        //                                 13.) RQC:    Flag Quad-CSA x-corr.
        //                                 14.) RQCdfs: Flag Quad-CSA DFS
        // Note                          : Functions which place a WBRE
        //                                 into a parameter set must contain
        //                                 add the information read here
 
   {
   SetLevel(pset);                              // Set R computation level
   SetType(pset);                               // Set R computation type
   if(DF) SetDip(pset);				// Set Dip R controls if desired
   else   Dip(0);				// else insure no Dip relaxation
   if(CF) SetSA(pset);				// Set SA R controls if desired
   else	  CSA(0);				// else insure no SA relaxation
   if(QF) SetQuad(pset);			// Set quad R controls if desired
   else	  Quad(0);				// else insure no quad relaxation
   if(DD && CC) SetDCX(pset);			// Set Dip-SA cross-corr controls
   else						// if desired else insure no dip-SA
     {						// cross correlation
     DC = 0;
     DCdfs = 0; 
     }
   if(DD && QQ) SetDQX(pset);			// Set Dip-Quad cross-corr controls
   else						// if desired, else insure no
     {						// dip-quad cross correlation
     DQ = 0;
     DQdfs = 0; 
     }
   if(QQ && CC) SetDQX(pset);			// Set Quad-SA cross corr controls
   else						// if desired else insure no
     {						// quad-SA cross correlation
     QC = 0;
     QCdfs = 0; 
     }
   }
 
                                                                                
// --------------------------------------------------------------------------------
// ------------------------------ PUBLIC FUNCTIONS --------------------------------
// --------------------------------------------------------------------------------
 
// ________________________________________________________________________________
// A    RELAXATION & EXCHANGE CONTROL CONSTRUCTION, ASSIGNMENT, DESTRUCTION
// ________________________________________________________________________________
 

WBRExch::WBRExch()

        // Input                none  :
        // Output               none  : WBR Relaxation & Exchange constructor
 
  {
  DD = 1;			// Set for dipolar relaxation
  CC = 1;			// Set for CSA relaxation
  QQ = 1;			// Set for quadrupolar relaxation

  DDdfs = 1;			// Set for dipolar DFS
  CCdfs = 1;			// Set for CSA DFS 
  QQdfs = 1;			// Set for quadrupolar DFS 

  DC = 1;			// Set for dip-CSA cross-correlation
  DCdfs = 1;			// Set for dip-CSA CC DFS

  DQ = 1;			// Set for dip-quad cross-correlation
  DQdfs = 1;			// Set for dip-quad CC DFS

  QC = 1;			// Set for quad-CSA cross-correlation
  QCdfs = 1;			// Set for quad-CSA CC DFS

  level = 4;			// Set for level 4 computation
  type = 0;			// Set for full cross-correlation
  }
 

WBRExch::WBRExch(const WBRExch& WBRE)
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        // Output               none  : New WBRExch is constructed
        //                              which is identical to WBRE
 
  {
  DD = WBRE.DD;			// Copy dipolar relaxation flag
  CC = WBRE.CC;			// Copy SA relaxation flag
  QQ = WBRE.QQ;			// Copy quadrupolar relaxation flag

  DDdfs = WBRE.DDdfs;		// Copy dipolar DFS flag
  CCdfs = WBRE.CCdfs;		// Copy SA DFS flag
  QQdfs = WBRE.QQdfs;		// Copy quadrupolar DFS flag

  DC = WBRE.DC;
  DCdfs = WBRE.DCdfs;

  DQ = WBRE.DQ;
  DQdfs = WBRE.DQdfs;

  QC = WBRE.QC;
  QCdfs = WBRE.QCdfs;

  level = WBRE.level;
  type = WBRE.type;
  }


WBRExch::~WBRExch() {}
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        // Output               none  : WBRExch is deleted
 


WBRExch& WBRExch::operator= (const WBRExch &WBRE)

        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      WBRE  : Coordinate
        // Output               none  : pt set identical to pt1

  {
  DD = WBRE.DD;
  CC = WBRE.CC;
  QQ = WBRE.QQ;

  DDdfs = WBRE.DDdfs;
  CCdfs = WBRE.CCdfs;
  QQdfs = WBRE.QQdfs;

  DC = WBRE.DC;
  DCdfs = WBRE.DCdfs;

  DQ = WBRE.DQ;
  DQdfs = WBRE.DQdfs;

  QC = WBRE.QC;
  QCdfs = WBRE.QCdfs;

  level = WBRE.level;
  type = WBRE.type;
  return *this;
  }


// ________________________________________________________________________________
// B             RELAXATION & EXCHANGE CONTROL ACCESS FUNCTIONS
// ________________________________________________________________________________

// --------------------------------------------------------------------------------
//			    Computation Level Access
// --------------------------------------------------------------------------------

void WBRExch::Level(int i)
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      i     : Flag for relaxation computation level
        // Output               none  : WBRE level set

   {
   if(i>4 || i<-4) i=4;
   level = i;
   }


int WBRExch::Level() const { return level; }
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        // Output               level : Relaxation computation level

// --------------------------------------------------------------------------------
//			    Computation Type Access
// --------------------------------------------------------------------------------


void WBRExch::Type(int i) { type = i; }
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      i     : Flag for relaxation computation type
        // Output               none  : WBRE type set
	// Note                       : Values of type are as follows

        //                                 0 = Auto- & Cross-Correlation
        //                                 + = Auto-Correlation only
        //                                 - = Cross-Correlation only


int WBRExch::Type() const { return type; }
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        // Output               type : Relaxation computation type


// --------------------------------------------------------------------------------
//		       Dipolar Relaxation Control Access
// --------------------------------------------------------------------------------


void WBRExch::Dip(int i)
   {
   if(!i)
     {
     DD = 0;				// No dipolar relaxation
     DDdfs=0;				// No dipolar dynamic frequency shifts
     DC = 0;				// No dipolar-CSA cross correlation
     DCdfs=0;				// No dipolar-CSA dyn. freq. shifts
     DQ = 0;				// No dipolar-Quad cross correlation
     DQdfs=0;				// No dipolar-Quad dyn. freq. shifts
     }
   else DD = 1;
   }


void WBRExch::DipDFS(int i)
   {
   if(!i) DDdfs = 0; 			// No dipolar dynamic frequency shifts
   else if(!DD) WBRerror(31);		// No dipolar DFS if no DD relaxation
   else DDdfs=1;			// Set for dipolar dyn. freq. shifts
   }


void WBRExch::DipCSA(int i)
   {
   if(!i)
     {
     DC = 0;				// No dipolar-CSA cross correlation
     DCdfs=0;				// No dipolar-CSA dyn. freq. shifts
     }
   else if(!DD) WBRerror(32);		// No dip-CSA DFS if no DD relaxation
   else if(!CC) WBRerror(33);		// No dip-CSA DFS if no CSA-CSA relax.
   else DC=1;				// Set for dipole-CSA cross correlation
   }


void WBRExch::DipCSADFS(int i)
   {
   if(!i) DCdfs=0;			// No dipolar-CSA dyn. freq. shifts
   else if(!DD) WBRerror(34);		// No dip-CSA DFS if no DD relaxation
   else if(!CC) WBRerror(35);		// No dip-CSA DFS if no CSA-CSA relax.
   else if(!DC) WBRerror(36);		// No dip-CSA DFS if no Dip-CSA x-terms
   else DCdfs=1;			// Set for dipole-CSA cross correlation
   }


void WBRExch::DipQuad(int i)
   {
   if(!i)
     {
     DQ = 0;				// No dipolar-Quad cross correlation
     DQdfs=0;				// No dipolar-Quad dyn. freq. shifts
     }
   else if(!DD) WBRerror(37);		// No dip-Quad DFS if no DD relaxation
   else if(!QQ) WBRerror(38);		// No dip-Quad DFS if no Q-Q relax.
   else DQ=1;				// Set for dipole-Quad cross correlat.
   }


void WBRExch::DipQuadDFS(int i)
   {
   if(!i) DQdfs=0;			// No dipolar-Quad dyn. freq. shifts
   else if(!DD) WBRerror(39);		// No dip-Quad DFS if no DD relaxation
   else if(!QQ) WBRerror(40);		// No dip-Quad DFS if no Q-Q relax.
   else if(!DQ) WBRerror(41);		// No dip-Quad DFS if no D-Q x-terms
   else DQdfs=1;			// Set for dipole-Quad cross correlat.
   }

// --------------------------------------------------------------------------------
//		       Shift Anisotropy Relaxation Control Access
// --------------------------------------------------------------------------------


void WBRExch::CSA(int i)
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      i     : Flag for relaxation computation type
        // Output               none  : Sets the parameters affecting the
	//				flags related to CSA relaxation
	// Note                       : Doesn't affect DFS or the
	//				cross-correlation terms

  {
  if(!i)				// Here if flagging NO CSA
    {
    CC = 0;				// No CSA relaxation
    CCdfs=0;				// No CSA dynamic freq. shifts
    QC = 0;				// No quad-CSA cross correlation
    QCdfs=0;				// No quad-CSA dyn. freq. shifts
    }
  else CC = 1;				// Set for CSA, NO DFS or CSA-X
  }


void WBRExch::CSADFS(int i)
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      i     : Flag for relaxation computation type
        // Output               none  : Sets the parameters affecting the
	//				flags related to CSA DFS relaxation
	// Note                       : Doesn't affect cross-correlation terms

  {
  if(!i) CCdfs = 0; 			// No CSA dynamic frequency shifts
  else if(!CC) WBRerror(42);		// No CSA DFS if no CSA relaxation
  else CCdfs=1;				// Set for CSA dyn. freq. shifts
  }

void WBRExch::CSADip(int i)    { DipCSA(i); }
void WBRExch::CSADipDFS(int i) { DipCSADFS(i); }


void WBRExch::CSAQuad(int i)
  {
  if(!i)
    {
    QC = 0;				// No CSA-Quad cross correlation
    QCdfs=0;				// No CSA-Quad dyn. freq. shifts
    }
  else if(!QQ) WBRerror(43);		// No CSA-Quad DFS if no Quad-Quad relaxation
  else if(!CC) WBRerror(44);		// No CSA-Quad DFS if no CSA-CSA relaxation
  else QC=1;				// Set for dipole-Quad cross correlation
  }


void WBRExch::CSAQuadDFS(int i)
  {
  if(!i) QCdfs=0;			// No CSA-Quad dyn. freq. shifts
  else if(!QQ) WBRerror(45);		// No CSA-Quad DFS if no Quad-Quad relaxation
  else if(!CC) WBRerror(46);		// No CSA-Quad DFS if no CSA-CSA relaxation
  else if(!QC) WBRerror(47);		// No CSA-Quad DFS if no CSA-Quad cross terms
  else QCdfs=1;				// Set for CSA-Quad cross correlation
  }

// --------------------------------------------------------------------------------
//		       Quadrupolar Relaxation Control Access
// --------------------------------------------------------------------------------


void WBRExch::Quad(int i)
  {
  if(!i)
    {
    QQ = 0;				// No quadrupolar relaxation
    QQdfs=0;				// No quadrupolar dynamic frequency shifts
    DQ = 0;				// No dipolar-Quad cross correlation
    DQdfs=0;				// No dipolar-Quad dyn. freq. shifts
    QC = 0;				// No quadrupolar-CSA cross correlation
    QCdfs=0;				// No quadrupolar-CSA dyn. freq. shifts
    }
  else QQ = 1;
  }


void WBRExch::QuadDFS(int i)
  {
  if(!i) QQdfs = 0; 			// No Quad dynamic frequency shifts
  else if(!QQ) WBRerror(48);		// No Quad DFS if no Quad relaxation
  else QQdfs=1;				// Set for Quad dyn. freq. shifts
  }


void WBRExch::QuadDip(int i) { DipQuad(i); }

void WBRExch::QuadDipDFSQuad(int i) { DipQuadDFS(i); }

void WBRExch::QuadCSA(int i) { CSAQuad(i); }

void WBRExch::QuadCSADFS(int i) { CSAQuadDFS(i); }
 

// --------------------------------------------------------------------------------
//		           Mutual Exchange Control Access
// --------------------------------------------------------------------------------

// --------------------------------------------------------------------------------
//		         Non-Mutual Exchange Control Access
// --------------------------------------------------------------------------------

 
// ________________________________________________________________________________
// C             RELAXATION & EXCHANGE PARAMETER AUXILIARY FUNCTIONS
// ________________________________________________________________________________


double WBRExch::LWhh(const sys_dynamic& sys, const std::string& Iso)
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      sys   : Dynamic spin system
        //                      Iso   : A spin isotope label
        // Output               LW    : Maximum linewidth at half-height (Hz)
        //                              estimated under simple assumptions
        // Note                       : Here, the motion is assumed to be diffusive and
        //                              the system moving as an isotropic top.

  {
  double lwhh = 0;
  if(DD) lwhh += LWhh_DD_max(sys, Iso); 
  if(CC) lwhh += LWhh_CC_max(sys, Iso); 
  if(QQ) lwhh += LWhh_QQ_max(sys, Iso); 
  return lwhh;
  }
 

// ________________________________________________________________________________
// D             RELAXATION & EXCHANGE PARAMETER SET FUNCTIONS
// ________________________________________________________________________________

//---------------------------------------------------------------------------------
//                Parameter Set Generation from WBR Parameters
//---------------------------------------------------------------------------------

WBRExch::operator ParameterSet( ) const

        // Input                WBRE  : Relaxation/exchange controls (this)
        // Output               pset  : Parameter set with
        //                              only WBRE parameters

  { ParameterSet pset; pset += *this; return pset; }     // Add in WBRE parameters


void operator+= (ParameterSet& pset, const WBRExch& WBRE)

        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      pset  : A parameter set
        // Output               pset  : Parameter set with
        //                              WBRE parameters added to it

   {
// sosi - this is still not complete
   std::string pname;
   std::string pdata;
   std::string pstate;
 
   SinglePar par;
   pname = std::string("Rlevel");			// Add Computation Level
   pstate = std::string("Relaxation Computation Level");
   pdata = WBRE.Level();
   par = SinglePar(pname, 0, pdata, pstate);
   pset.push_back(par);
 
   pname = std::string("Rtype");			// Add Computation Type
   pstate = std::string("Relaxation Computation Type");
   pdata = WBRE.Type();
   par = SinglePar(pname, 0, pdata, pstate);
   pset.push_back(par);
 
   return;
   }

//---------------------------------------------------------------------------------
//                  WBR Specification From Parameter Set
//---------------------------------------------------------------------------------

void WBRExch::SetZero()
 
        // Input                WBRE	: Relaxation/exchange controls (this)
        // Output               none	: All relaxation and exchange flags
	//				  are set to zero

  {
  DD = 0;				// No dipolar relaxation
  DDdfs = 0;				// No dipolar DFS 
  CC = 0;				// No shift anisotropy relaxation
  CCdfs = 0;				// No shift anisotropy DFS 
  QQ = 0;				// No quadrupolar relaxation
  QQdfs = 0;				// No quadrupolar DFS 
  DC = 0;				// No dip-shift cross correlation
  DCdfs = 0;				// No dip-shift DFS 
  DQ = 0;				// No dip-quad cross correlation
  DQdfs = 0;				// No dip-quad DFS 
  QC = 0;				// No quad-shift cross correlation
  QCdfs = 0;				// No quad-shift DFS 
  }


void WBRExch::SetLevel(const ParameterSet& pset)
 
        // Input                WBRE	: Relaxation/exchange controls (this)
        //                      pset	: A parameter set
        // Output               none	: WBRE level is set from parameter
	//				  Rlevel in pset
	// Note				: Default level is 4

   {
   std::string pstate, pname("Rlevel");		// Relax comput. level
   SinglePar par(pname);			// Parameter for this name
   ParameterSet::const_iterator item;	// Pix in parameter list
   item = pset.seek(par); 			// Pix in p. list for Rlevel
   if(item != pset.end())			// Retrieve the comput. level
     (*item).parse(pname,level,pstate);
   else
     {
     level = 4;
     WBRerror(2, pname, 1);			// Cannot read parameter
     WBRerror(9);				// Set level to 4
     }
   }


void WBRExch::SetType(const ParameterSet& pset)
 
        // Input                WBRE	: Relaxation/exchange controls (this)
        //                      pset	: A parameter set
        // Output               none	: WBRE type is set from parameter
	//				  Rtype in pset
	// Note				: Default type is 0

  {
  std::string pstate, pname("Rtype");		// Relaxation computation type 
  SinglePar par(pname);				// Parameter for this name
  ParameterSet::const_iterator item;		// Pix in parameter list
  item = pset.seek(par); 			// Pix in p. list for Rlevel
  if(item != pset.end()) 			// Retrieve computation type
    (*item).parse(pname,type,pstate);
  else
    {
    type= 0;
    WBRerror(2, pname, 1);			// Cannot read parameter
    WBRerror(10);				// Set type to 0
    }
  }


void WBRExch::SetDip(const ParameterSet& pset)
 
        // Input                WBRE	: Relaxation/exchange controls (this)
        //                      pset	: A parameter set
        // Output               none	: WBRE dipolar relaxation flags are
	//				  is set from parameters in pset
	// Note				: Parameters are as follows:
	//				  RDD    - dipolar relaxation flag
	//				  RDDdfs - dipolar DFS flag

  {
  std::string pstate, pname("RDD");			// Dipolar relaxation flag
  SinglePar par(pname);				// Parameter for this name
  ParameterSet::const_iterator item;		// Pix in parameter list
  item = pset.seek(par); 			// Pix in p. list for RDD
  int ival;					// Dummy iteger for value
  if(item != pset.end())			// Retrieve dipolar relax. flag
    {
    (*item).parse(pname,ival,pstate);
    if(ival) DD = 1;
    else     DD = 0;
    pname = std::string("RDDdfs");			// Dipolar DFS relaxation flag
    par = SinglePar(pname);			// Pix in p. list for Rlevel
    item = pset.seek(par); 			// Pix in p. list for RDDdfs
    if(item != pset.end())			// Retrieve dipolar relax. flag
      {
      (*item).parse(pname,ival,pstate);
      if(ival) DDdfs = 1;
      else     DDdfs = 0;
      }
    else
      {
      DDdfs = 1; 
      WBRerror(2, pname, 1);			// Cannot read parameter
      WBRerror(11, pname);			// Setting Dip DFS flag to 1
      }
    }
  else
    {
    DD = 1;
    DDdfs = 1; 
    WBRerror(2, pname, 1);			// Cannot read parameter
    WBRerror(12, pname);			// Setting Dip relaxation flag to 1
    WBRerror(11, pname);			// Setting Dip DFS flag to 1
    }
  }


void WBRExch::SetSA(const ParameterSet& pset)
 
        // Input                WBRE	: Relaxation/exchange controls (this)
        //                      pset	: A parameter set
        // Output               none	: WBRE SA relaxation flags are
	//				  is set from parameters in pset
	// Note				: Parameters are as follows:
	//				  RCC    - SA relaxation flag
	//				  RCCdfs - SA DFS flag

  {
  std::string pstate, pname("RCC");			// SA relaxation flag
  SinglePar par(pname);				// Parameter for this name
  ParameterSet::const_iterator item;		// Pix in parameter list
  item = pset.seek(par); 			// Pix in p. list for RCC
  int ival;					// Dummy integer for value
  if(item != pset.end())			// Retrieve the CSA relaxation flag
    {
    (*item).parse(pname,ival,pstate);
    if(ival) CC = 1;
    else     CC = 0;
    pname = std::string("RCCdfs");			// CSA DFS relaxation flag
    par = SinglePar(pname);			// Pix in p. list for Rlevel
    item = pset.seek(par);
    if(item != pset.end())			// Retrieve the dipolar relaxation flag
      {
      (*item).parse(pname,ival,pstate);
      if(ival) CCdfs = 1;
      else     CCdfs = 0;
      }
    else
      {
      CCdfs = 1; 
      WBRerror(2, pname, 1);			// Cannot read parameter
      WBRerror(13, pname);			// Setting SA DFS flag to 1
      }
    }
  else
    {
    CC = 1;
    CCdfs = 1; 
    WBRerror(2, pname, 1);			// Cannot read parameter
    WBRerror(14, pname, 1);			// Setting SA relax flag to 1
    WBRerror(13, pname);			// Setting SA DFS flag to 1
    }
  }



void WBRExch::SetQuad(const ParameterSet& pset)
 
        // Input                WBRE	: Relaxation/exchange controls (this)
        //                      pset	: A parameter set
        // Output               none	: WBRE quad relaxation flags are
	//				  is set from parameters in pset
	// Note				: Parameters are as follows:
	//				  RQQ    - quadrupolar relaxation flag
	//				  RQQdfs - quadrupolar DFS flag

  {
  std::string pstate, pname("RQQ");			// Quad relaxation flag
  SinglePar par(pname);				// Parameter for this name
  ParameterSet::const_iterator item;		// Pix in parameter list
  item = pset.seek(par); 			// Pix in p. list for RQQ
  int ival;					// Dummy integer for value
  if(item != pset.end())			// Retrieve the Quad relaxation flag
    {
    (*item).parse(pname,ival,pstate);
    if(ival) QQ = 1;
    else     QQ = 0;
    pname = std::string("RQQdfs");			// Quadrupolar DFS relaxation flag
    par = SinglePar(pname);			// Parameter for this name
    item = pset.seek(par);			// Pix in pset for RQQdfs
    if(item != pset.end())			// Retrieve the dipolar relaxation flag
      {
      (*item).parse(pname,ival,pstate);
      if(ival) QQdfs = 1;
      else     QQdfs = 0;
      }
    else
      {
      QQdfs = 1; 
      WBRerror(2, pname, 1);			// Cannot read parameter
      WBRerror(15, pname);			// Setting Quad DFS flag to 1
      }
    }
  else
    {
    QQ = 1;
    QQdfs = 1; 
    WBRerror(2, pname, 1);			// Cannot read parameter
    WBRerror(16, pname, 1);			// Setting Quad relax flag to 1
    WBRerror(15, pname);			// Setting Quad DFS flag to 1
    }
  }


void WBRExch::SetDCX(const ParameterSet& pset)
 
        // Input                WBRE	: Relaxation/exchange controls (this)
        //                      pset	: A parameter set
        // Output               none	: WBRE dipolar-SA cross correlation flags
	//				  are set from parameters in pset
	// Note				: Parameters are as follows:
	//				  RDC    - dipolar-SA cross correlation flag
	//				  RDCdfs - dipolar-SA cross correlation  DFS flag
	// Note				: The dipolar and SA flags should be
	//				  set before this is called

  {
  if(!DD || !CC) return;			// Quit if no Dip or no SA
  std::string pstate, pname("RDC");			// Dip-CSA relaxation flag
  SinglePar par(pname);				// Parameter for this name
  ParameterSet::const_iterator item;		// Pix in parameter list
  item = pset.seek(par); 			// Pix in p. list for RDC
  int ival;
  if(item != pset.end())			// Retrieve the Dip-CSA relaxation flag
    {
    (*item).parse(pname,ival,pstate);
    if(ival) DC = 1;
    pname = std::string("RDCdfs");			// Dip-CSA DFS relaxation flag
    par = SinglePar(pname);			// Parameter for this name
    item = pset.seek(par);			// Pix in pset for RDCdfs
    if(item != pset.end())			// Retrieve the dipolar relaxation flag
      {
      (*item).parse(pname,ival,pstate);
      if(ival) DCdfs = 1;
      }
    else
      {
      DCdfs = 1; 
      WBRerror(2, pname, 1);			// Cannot read parameter
      WBRerror(17, pname);			// Set DCdfs to 1
      }
    }
  else
    {
    DC = 1;
    DCdfs = 1; 
    WBRerror(2, pname, 1);			// Cannot read parameter
    WBRerror(18, pname, 1);			// Setting DC Xcorr flag to 1
    WBRerror(17, pname);			// Set DCdfs to 1
    }
  }



void WBRExch::SetDQX(const ParameterSet& pset)
 
        // Input                WBRE	: Relaxation/exchange controls (this)
        //                      pset	: A parameter set
        // Output               none	: WBRE dipolar-SA cross correlation flags
	//				  are set from parameters in pset
	// Note				: Parameters are as follows:
	//				  RDC    - dipolar-SA cross correlation flag
	//				  RDCdfs - dipolar-SA cross correlation  DFS flag
	// Note				: The dipolar and SA flags should be
	//				  set before this is called

  {
  if(!DD || !QQ) return;			// Quit if no Dip or no quad
  std::string pstate, pname("RDQ");			// Dip-Quad relaxation flag
  SinglePar par(pname);				// Parameter for this name
  ParameterSet::const_iterator item;		// Pix in parameter list
  item = pset.seek(par); 			// Pix in p. list for RDQ
  int ival;
  if(item != pset.end())			// Get Dip-Quad relaxation flag
    {
    (*item).parse(pname,ival,pstate);
    if(ival) DQ = 1;
    pname = std::string("RDQdfs");			// Dip-Quad DFS relaxation flag
    par = SinglePar(pname);			// Parameter for this name
    item = pset.seek(par);			// Pix in pset for RDQdfs
    if(item != pset.end())			// Get DQ DFS relaxation flag
      {
      (*item).parse(pname,ival,pstate);
      if(ival) DQdfs = 1;
      }
    else
      {
      DQdfs = 1; 
      WBRerror(2, pname, 1);			// Cannot read parameter
      WBRerror(19, pname);			// Set DQdfs to 1
      }
    }
  else
    {
    DQ = 1;
    DQdfs = 1; 
    WBRerror(2, pname, 1);			// Cannot read parameter
    WBRerror(20, pname, 1);			// Setting DQ Xcorr flag to 1
    WBRerror(19, pname);			// Set DQdfs to 1
    }
  }


void WBRExch::SetQCX(const ParameterSet& pset)
 
        // Input                WBRE	: Relaxation/exchange controls (this)
        //                      pset	: A parameter set
        // Output               none	: WBRE quad-SA cross correlation flags
	//				  are set from parameters in pset
	// Note				: Parameters are as follows:
	//				  RQC    - quad-SA cross correlation flag
	//				  RQCdfs - quad -SA cross correlation  DFS flag
	// Note				: The quad and SA flags should be
	//				  set before this is called

  {
  if(!QQ || !CC) return;			// Quit if no Quad or no SA
  std::string pstate, pname("RQC");			// Quad-CSA relaxation flag
  SinglePar par(pname);				// Parameter for this name
  ParameterSet::const_iterator item;		// Pix in parameter list
  item = pset.seek(par); 			// Pix in p. list for RQC
  int ival;
  if(item != pset.end())			// Retrieve the Quad-CSA relaxation flag
    {
    (*item).parse(pname,ival,pstate);
    if(ival) QC = 1;
    pname = std::string("RQCdfs");			// Quad-CSA DFS relaxation flag
    par = SinglePar(pname);			// Parameter for this name
    item = pset.seek(par);			// Pix in pset for RQCdfs
    if(item != pset.end())			// Retrieve the QC DFS flag
      {
      (*item).parse(pname,ival,pstate);
      if(ival) QCdfs = 1;
      }
    else
      {
      QCdfs = 1; 
      WBRerror(2, pname, 1);			// Cannot read parameter
      WBRerror(21, pname);			// Set QCdfs to 1
      }
    }
  else
    {
    QC = 1;
    QCdfs = 1; 
    WBRerror(2, pname, 1);			// Cannot read parameter
    WBRerror(22, pname, 1);			// Setting QC Xcorr flag to 1
    WBRerror(21, pname);			// Set QCdfs to 1
    }
  }
 

WBRExch& WBRExch::operator= (const ParameterSet& pset)
 
        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      pset	 : A parameter set
        // Output               none     : WBRE filled with relaxation and 
	//				   exchange control parameters
        // Note                          : Three things are gleaned from
        //                                 the parameter set for a WBRE
        //                                 1.) The level of the computation
        //                                 2.) The type of the computation
        //                                 3.) Flag for dipolar relaxation
        // Note                          : Functions which place a WBRE
        //                                 into a parameter set must contain
        //                                 add the information read here
 
   {
   SetZero();					// Zero all flags
   SetLevel(pset);				// Set R computation level
   SetType(pset);				// Set R computation type
   SetDip(pset);				// Set dipolar relaxation flags
   SetSA(pset);					// Set SA relaxation flags
   SetQuad(pset);				// Set quadrupolar relaxation flags 
   SetDCX(pset);				// Set Dip-SA cross-corr flags
   SetDQX(pset);				// Set Dip-SA cross-corr flags
   SetQCX(pset);				// Set Quad-SA cross-corr flags
   return *this;
   }


// ________________________________________________________________________________
// E                RELAXATION & EXCHANGE OUTPUT FUNCTIONS
// ________________________________________________________________________________


std::ostream& WBRExch::printDip(std::ostream& ostr) const

        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      ostr     : An output stream
        // Output               ostr     : The output stream modified by
        //                                 the WBRE Dipolar parameters

  {
  ostr << "\n\t** Dipolar Relaxation: ";		// Main Header
  if(!DD)						// Dipole relaxation?
    {							// If no, then exit
    ostr << "No";
    ostr.flush();
    return ostr;
    }
  std::string spacer(40, ' ');  
  ostr << "         Yes";				// If yes see which types
  if(DDdfs) ostr << "\n" << spacer << "- with DFS";	// Here if D-D DFS
  else      ostr << "\n" << spacer << "- no DFS";
  if(!DC) ostr << "\n" << spacer << "- no Dipole-CSA"; 	// D-SA cross-correlation?
  else							// If there is this type
    {							// then output that and see
    ostr << "\n" << spacer << "- with Dipole-CSA";	// if D-SA DFS included
    if(DCdfs) ostr << ", with DFS ";
    else      ostr << ", no DFS ";
    }
  if(!DQ) ostr << "\n" << spacer << "- without Dipole-Quad";
  else
    {
    ostr << "\n" << spacer << "- with Dipole-Quad";
    if(DQdfs) ostr << ", with DFS";
    else      ostr << ", no DFS ";
    }
  return ostr;
  }


std::ostream& WBRExch::printSA(std::ostream& ostr) const

        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      ostr     : An output stream
        // Output               ostr     : The output stream modified by
        //                                 the WBRE Shift Anisotropy parameters

  {
  ostr << "\n\t** Shift Anisotropy Relaxation: ";
  if(!CC)						// Shift anisotropy relaxation?
    {							// If no, then exit
    ostr << "No";
    ostr.flush();
    return ostr;
    }
  std::string spacer(40, ' ');  
  ostr << "Yes";
  if(CCdfs) ostr << "\n" << spacer << "- with DFS";	// Here if C-C DFS
  else      ostr << "\n" << spacer << "- no DFS";
  if(!DC) ostr << "\n" << spacer << "- no Dipole-CSA"; 	// D-SA cross-correlation?
  else
    {
    ostr << "\n" << spacer << "- with Dipole-CSA";	// if D-SA DFS included
    if(DCdfs) ostr << ", with DFS ";
    else      ostr << ", no DFS ";
    }
  if(!QC) ostr << "\n" << spacer << "- no Quad-CSA"; 	// Q-SA cross-correlation?
  else
    {
    ostr << "\n" << spacer << "- with Quad-CSA";	// if Q-SA DFS included
    if(QCdfs) ostr << ", with DFS ";
    else      ostr << ", no DFS ";
    }
  return ostr;
  }


std::ostream& WBRExch::printQ(std::ostream& ostr) const

        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      ostr     : An output stream
        // Output               ostr     : The output stream modified by
        //                                 the WBRE Quadrupolar parameters

  {
  ostr << "\n\t** Quadrupolar Relaxation:      ";
  if(!QQ)						// Quadrupolar relaxation?
    {							// If no, then exit
    ostr << "No";
    ostr.flush();
    return ostr;
    }
  std::string spacer(40, ' ');  
  ostr << "Yes";
  if(QQdfs) ostr << "\n" << spacer << "- with DFS";	// Here if Q-Q DFS
  else      ostr << "\n" << spacer << "- no DFS";
  if(!DQ) ostr << "\n" << spacer << "- no Dipole-Quad";	// D-Q cross-correlation?
  else
    {
    ostr << "\n" << spacer << "- with Dipole-Quad";	// if D-Q DFS included
    if(DQdfs) ostr << ", with DFS ";
    else      ostr << ", no DFS ";
    }
  if(!QC) ostr << "\n" << spacer << "- no Quad-CSA"; 	// Q-SA cross-correlation?
  else
    {
    ostr << "\n" << spacer << "- with Quad-CSA";	// if Q-SA DFS included
    if(QCdfs) ostr << ", with DFS ";
    else      ostr << ", no DFS ";
    }
  return ostr;
  }


std::ostream& WBRExch::print(std::ostream& ostr) const

        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      ostr     : An output stream
        // Output               ostr     : The output stream modified by
        //                                 the WBRE parameters

  {
  ostr << "\n\tSpecified Relaxation and Exchange Effects";
  ostr << "\n\t-----------------------------------------";
  ostr << "\n\n\t* The Computation Level is " << level;  
  if(type > 0)      ostr << "\n\t* Only Auto-Correlation Terms Are Included\n";
  else if(type < 0) ostr << "\n\t* Only Cross-Correlation Terms Are Included\n";
  else              ostr << "\n\t* Both Auto- and Cross-Correlation Are Terms Included\n";
  printDip(ostr);				// Output dipolar contributions
  printSA(ostr);				// Output SA contributions
  printQ(ostr);					// Output quad contributions
// print G(ostr)				// Output g-tensor contributions
// print HF(ostr)				// Output hyperfine contributions
// print RF(ostr)				// Output random field contributions
// print K(ostr)				// Output exchange contributions
  return ostr;
  }
 
 
std::ostream &operator << (std::ostream &ostr, const WBRExch& WBRE)
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      ostr  : Output stream
        // Output               none  : WBRE info is sent
        //                              to the output stream
 
  { return WBRE.print(ostr); }


// ________________________________________________________________________________
// F                RELAXATION & EXCHANGE INPUT FUNCTIONS
// ________________________________________________________________________________


void WBRExch::read(const std::string& filename)

        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      filename : Input filename
        // Output               none     : WBRE parameters filled with
        //                                 parameters read from file
        // Note                          : The file should be an ASCII file
        //                                 containing recognized WBR parameters

   {
   ParameterSet pset;                // Declare a parameter set
   std::ifstream inp(filename.c_str());	// Open filename for input
   if(!inp.good())                      // If file bad then exit
     {   
     WBRerror(1, filename, 1);		//	Problems with file
     WBRfatality(50);			//	Can't construct WBR
     }   
   SinglePar par;
   while(par.read(inp))                 // Read all file parameters
     {   
     if(!pset.contains(par))            // Add them to the parameter list if
       pset.push_back(par);                   // not already included
     }
   (*this) = pset;                      // Fill up spin_sys with parameters
   return;
   }
 

void WBRExch::read(const std::string& filename, const sys_dynamic& sys)

        // Input                WBRE     : Relaxation/exchange controls (this)
        //                      filename : Input filename
	//			sys      : A dynamic spin system
        // Output               none     : WBRE parameters filled with
        //                                 parameters read from file
	//				   under the constraints set by the	
	//				   input spin system
        // Note                          : The file should be an ASCII file
        //                                 containing recognized WBR parameters

   {
//		Read Input File For WBRExch Parameters

   ParameterSet pset;                // Declare a parameter set
   std::ifstream inp(filename.c_str()); 	// Open filename for input
   if(!inp.good())                      // If file bad then exit
     {   
     WBRerror(1, filename, 1);		//	Problems with file
     WBRfatality(50);			//	Can't construct WBR
     }   
   SinglePar par;			// Working parameter
   while(par.read(inp))                 // Read all file parameters
     {   
     if(!pset.contains(par))            // Add them to the parameter list if
       pset.push_back(par);                   // not already included
     }
//		Set WBRExch Parameters From Those in File

   (*this) = pset;                      // Fill up spin_sys with parameters

//	Filter Specified WBRExch Parameters Based on Info in Spin System

   int DipT = 0;			// Flag for no dipolar relaxation
   int CSAT = 0;			// Flag for no CSA relaxation
   int QuadT = 0;			// Flag for no quad. relaxation
   int i=0, j=0;
   int ns = sys.spins();
   for(i=0; i<ns; i++)
     {
     if(sys.delz(i)) CSAT++;		// Check if any CSA tensors exist
     if(sys.QCC(i)) QuadT++;		// Check if any quad tensors exist
     }
   for(i=0; i<ns-1; i++)
     for(j=i+1; j<ns; j++)
       if(sys.DCC(i,j)) DipT++;		// Check if any dipolar tensors exist
// sosi - should issue warnings here if requested terms are not possible
//        because they are not in the spin system parameters

   if(!DipT)				// No dipolar couplings in spin system
     {					// Insure none are set in WBRE
     DD = 0;				// No dipolar relaxation
     DDdfs = 0;				// No dipolar DFS 
     DC = 0;				// No dip-shift cross correlation
     DCdfs = 0;				// No dip-shift DFS 
     DQ = 0;				// No dip-quad cross correlation
     DQdfs = 0;				// No dip-quad DFS 
     }
   if(!CC)				// No shift tensors in spin system
     {					// Insure none are set in WBRE
     CC = 0;				// No shift anisotropy relaxation
     CCdfs = 0;				// No shift anisotropy DFS 
     DC = 0;				// No dip-shift cross correlation
     DCdfs = 0;				// No dip-shift DFS 
     QC = 0;				// No quad-shift cross correlation
     QCdfs = 0;				// No quad-shift DFS 
     }
   if(!QQ)				// No quadrupolar tensors in spin sys
     {					// Insure none are set in WBRE
     QQ = 0;				// No quadrupolar relaxation
     QQdfs = 0;				// No quadrupolar DFS 
     DQ = 0;				// No dip-quad cross correlation
     DQdfs = 0;				// No dip-quad DFS 
     QC = 0;				// No quad-shift cross correlation
     QCdfs = 0;				// No quad-shift DFS 
     }

   return;
   }
 
 
void WBRExch::ask_read(int argc, char* argv[], int argn)
 
        // Input                WBRE    : Relaxation/exchange controls (this)
        //                      argc    : Number of arguments
        //                      argv    : Vector of argc arguments
        //                      argn    : Argument index
        // Output               void    : The parameter argn of array argc
        //                                is used to supply a filename
        //                                from which WBRE parameters are read
        //                                If the argument argn is not in argv,
        //                                user is asked to supply a filename
        // Note                         : The file should be an ASCII file
        //                                containing recognized WBRE parameters
        // Note                         : The WBRE is modifed (filled)
 
  {
  std::string filename;                            // Name of spin system file
  query_parameter(argc, argv, argn,           // Get filename from command
   "\n\tRelaxation/Exchange Controls Filename? "
                                  , filename); 	// Or ask for it
  read(filename);                             // Read parameters from filename
  }
 
 
void WBRExch::ask_read(int argc, char* argv[], int argn, const sys_dynamic& sys)
 
        // Input                WBRE    : Relaxation/exchange controls (this)
        //                      argc    : Number of arguments
        //                      argv    : Vector of argc arguments
        //                      argn    : Argument index
	//			sys     : A dynamic spin system
        // Output               void    : The parameter argn of array argc
        //                                is used to supply a filename
        //                                from which WBRE parameters are read
	//				  under the constraints set by the	
	//				  input spin system
        //                                If the argument argn is not in argv,
        //                                the user is asked to supply a filename
        // Note                         : The file should be an ASCII file
        //                                containing recognized WBRE parameters
        // Note                         : The WBRE is modifed (filled)
 
    {
    std::string filename;                            // Name of spin system file
    query_parameter(argc, argv, argn,           // Get filename from command
     "\n\tRelaxation/Exchange Controls Filename? "
                                  , filename); // Or ask for it
    int DD = 0;					// Flag for no dipolar relax.
    int CC = 0;					// Flag for no CSA relaxation
    int QQ = 0;					// Flag for no quad. relax.
    int i=0, j=0;
    int ns = sys.spins();
    for(i=0; i<ns; i++)
      {
      if(sys.delz(i)) CC++;			// Check if any CSA tensors 
      if(sys.QCC(i)) QQ++;			// Check if any quad tensors
      }
    for(i=0; i<ns-1; i++)
      for(j=i+1; j<ns; j++)
        if(sys.DCC(i,j)) DD++;			// Check if any dipolar tensors
    read(filename, sys);			// Read parameters from filename
    }

                                                                       
// ________________________________________________________________________________
// F              INTERACTIVE RELAXATION & EXCHANGE SETUP FUNCTIONS
// ________________________________________________________________________________
 

void WBRExch::ask(int argc, char* argv[], int& argn)

        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      ostr     : An output stream
        // Output               ostr     : The output stream modified by

    {
    std::string yn;

//			  Pure Dipolar Relaxation

    query_parameter(argc, argv, argn,                   // Dipolar relaxation
         "\n\tInclude Dipolar Relaxation (y/n)? ", yn);
    argn++;
    if(yn=="y")
      {
      DD=1;						// Set for D-D relax.
      query_parameter(argc, argv, argn,			// Dipolar DFS relax.
         "\n\tInclude Dipolar Dynamic Frequency Shifts (y/n)? ", yn);
      argn++;
      if(yn=="y") DDdfs=1;
      else        DDdfs=0;
      }
    else
      {
      DD=0;
      DDdfs=0;
      }

//			     Pure CSA Relaxation

    query_parameter(argc, argv, argn,                   // CSA relaxation
     "\n\tInclude CSA Relaxation (y/n)? ", yn);
    argn++;
    if(yn=="y")
      {
      CC=1;
      query_parameter(argc, argv, argn,			// CSA DFS relaxation effects
         "\n\tInclude CSA Dynamic Frequency Shifts (y/n)? ", yn);
      argn++;
      if(yn=="y") CCdfs=1;
      else        CCdfs=0;
      }
    else
      {
      CC=0;
      CCdfs=0;
      }

//			     Dipole-CSA Cross Correlation

    if(DD && CC)
      {
      query_parameter(argc, argv, argn,                   // Dipole-CSA relaxation
           "\n\tInclude Dip-CSA Cross Correlation (y/n)? ", yn);
      argn++;
      if(yn=="y")
        {
        DC=1;
        query_parameter(argc, argv, argn,		// Dip-CSA DFS effects
           "\n\tInclude Dip-CSA Dynamic Frequency Shifts (y/n)? ", yn);
        if(yn=="y") DCdfs=1;
        else        DCdfs=0;
        }
      else
        {
        DC=0;
        DCdfs=0;
        }
      }

//			     Pure Quad Relaxation

    query_parameter(argc, argv, argn,                   // Quad relaxation
     "\n\tInclude Quadrupolar Relaxation (y/n)? ", yn);
    argn++;
    if(yn=="y")
      {
      QQ=1;
      query_parameter(argc, argv, argn,			// Quad DFS relaxation effects
         "\n\tInclude Quad Dynamic Frequency Shifts (y/n)? ", yn);
      argn++;
      if(yn=="y") QQdfs=1;
      else        QQdfs=0;
      }
    else
      {
      QQ=0;
      QQdfs=0;
      }

//			     Dipole-Quad Cross Correlation

    if(DD && QQ)
      {
      query_parameter(argc, argv, argn,                   // Dipole-Quad relaxation
           "\n\tInclude Dip-Quad Cross Correlation (y/n)? ", yn);
      argn++;
      if(yn=="y")
        {
        DQ=1;
        query_parameter(argc, argv, argn,		// Dip-Quad DFS effects
           "\n\tInclude Dip-Quad Dynamic Frequency Shifts (y/n)? ", yn);
        if(yn=="y") DQdfs=1;
        else        DQdfs=0;
        }
      else
        {
        DQ=0;
        DQdfs=0;
        }
      }

//			     Quad-CSA Cross Correlation

    if(DD && CC)
      {
      query_parameter(argc, argv, argn,                   // Quad-CSA relaxation
           "\n\tInclude Quad-CSA Cross Correlation (y/n)? ", yn);
      argn++;
      if(yn=="y")
        {
        QC=1;
        query_parameter(argc, argv, argn,		// Quad-CSA DFS effects
           "\n\tInclude Quad-CSA Dynamic Frequency Shifts (y/n)? ", yn);
        if(yn=="y") QCdfs=1;
        else        QCdfs=0;
        }
      else
        {
        QC=0;
        QCdfs=0;
        }
      }
    }

                                                                       
// ________________________________________________________________________________
// G           INTERACTIVE RELAXATION & EXCHANGE SUPEROPERATOR PRODUCTION
// ________________________________________________________________________________
 
 
super_op WBRExch::REX(const sys_dynamic& sys, gen_op& Ho, int fext)
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      sys   : Dynamic spin system
        //                      Ho    : Static Hamiltonian
//                      fext    : Signals external functions
        // Output               LOp   : WBR Relaxation and Exchange
	//				superoperator
        // Note                         Computed in the eigenbasis of Ho
 
   {
   int hs = sys.HS();                           // Total system Hilbert space
   int ls = hs*hs;                              // Total system Liouville space
   Ho.set_EBR();                                // Insure Ho in its eigenbasis (Hz)
   matrix mx(ls, ls, complex0);			// Construct zero superoperator
   super_op LOp(mx, Ho.get_basis());		// The returned superoperator
   super_op LOp0(mx, Ho.get_basis());
   super_op LOpDFS(mx, Ho.get_basis());		// For dynamic shift superoperators
   super_op LOpX(mx, Ho.get_basis());
   super_op LOpXdfs(mx, Ho.get_basis());

// --------------------------------------------------------------------------------
//                      Prepare the Spectral Density Components
//        (5 Tau Values, A Chi Value, & Lab. Frame Transition Frequencies)
// --------------------------------------------------------------------------------

   double taus[5];                              // Get 5 taus for spec. densities
   taust(taus, sys.taus());			// for a rigid asymmetric top
   double chi = chit(sys.taus());               // Get the system chi value
   double *w;                                // For system energy levels (LAB)
   w = new double[hs];
   if(abs(level) > 1)                           // For higher level computations
     {                                          // where spect.dens. are evaluated
     gen_op HZ = Hz(sys);			// Zeeman Hamiltonian (lab frame)
     HZ.Op_base(Ho, 1.e-7);			// Put Hz into basis of Ho
     if(!HZ.test_EBR()) WBRerror(60);		// Insure the result is proper
     for(int k=0; k<hs; k++)			// Calculate trans. frequencies
       w[k] = Re(HZ.get(k,k)) + Re(Ho.get(k,k));// in the laboratory frame
     }   

// --------------------------------------------------------------------------------
//                     Prepare Dipolar Relaxation Components
//		      (Only If Dipolar Effects Are Included)
// --------------------------------------------------------------------------------

   int ns = sys.spins();			// Total number of spins
   int i, j;
   matrix xiDs; 				// Array for Dipolar int. constants
   spin_T *TDip = NULL;				// Spin tensors for each dipole
   space_T *ADip = NULL;			// Spatial tensors for each dipole
   if(DD)
     {
     xiDs = xiD(sys);				// Dipolar interaction constants
						// Prepare the Spin Tensors
     int ns = sys.spins();			// Total number of spins
     int ndip = sys.dipoles();			// Total number of dipoles
     TDip = new spin_T[ndip];			// Space for the dipolar spin tensors
     int ij=0;
     for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
       for(j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
         if(Re(xiDs.get(i,j)))
           {
           TDip[ij] = T_D(sys,i,j);
           ij++;
           }
     ADip = new space_T[ndip];			// Space for dipolar spatial tensors
     ij=0;
     for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
       for(j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
         if(Re(xiDs.get(i,j)))
           ij++;
     }

// --------------------------------------------------------------------------------
//                       Prepare CSA Relaxation Components
//	             (Only If Shielding Effects Are Included)
// --------------------------------------------------------------------------------

   matrix xiCs; 				// Matrix for CSA interact. consts.
   spin_T *TCSA = NULL;				// Spin tensors for each spin
   space_T *ACSA = NULL;			// Spatial tensors for each spin
   if(CC)					// Generate values only if CSA
     {						// is to be considered 
     xiCs = xiCSA(sys);				//   Matrix of Xi's for CSA
     TCSA = new spin_T[ns];			//   Space for CSA spin tensors
     ACSA = new space_T[ns];			//   Space for CSA spatial tensors
     for(i=0; i<ns; i++)			//   Set CSA applicable tensors
       {					//   Note that we skip any spins
       if(Re(xiCs.get(i,i)))			//   having no shift anisotropy
         {
         ACSA[i] = sys.TC(i);
         TCSA[i] = T_CS2(sys,i);
         }
       }
     }


// --------------------------------------------------------------------------------
//                    Prepare Quadrupolar Relaxation Components
//	             (Only If Quadrupolar Effects Are Included)
// --------------------------------------------------------------------------------

   matrix xiQs;					// Matrix Quad. interact. constants
   spin_T *TQuad = NULL;			// Spin tensors for each spin
   space_T *AQuad = NULL;			// Spatial tensors for each spin
   if(QQ)					// Generate values only if Quad
     {						// interactions are considered 
     xiQs = xiQ(sys);				//   Matrix of Xi's for quadrupoles
     TQuad = new spin_T[ns];			//   Space for Quad spin tensors
     AQuad = new space_T[ns];			//   Space for Quad spatial tensors
     for(i=0; i<ns; i++)			//   Set Quad applicable tensors
       {					//   Note that we skip any spins
       if(xiQs.getRe(i,i))			//   having no quadrupole moment
         {
         AQuad[i] = sys.TQ(i);
         TQuad[i] = T_Q(sys,i);
         }
       }
     }

// --------------------------------------------------------------------------------
//                  Determine Dipolar Relaxation Superoperator
// --------------------------------------------------------------------------------
 
   if(DD)
     {
if(fext)
  {
std::cout << "\n\t\tadding external pure dipolar"; 
     Rijkl(LOp, sys, Ho, w, xiDs, xiDs, ADip, ADip,	// Dipole-dipole relaxation superop
                  TDip, TDip, taus, chi, type, level);
     if(DDdfs)
       {
std::cout << "\n\t\tadding external dipolar DFS"; 
       LOpDFS = LOp0;
       Rijklds(LOpDFS, sys, Ho, w, xiDs,xiDs,ADip,ADip,	// Dip-Dip DFS relaxation superop
                  TDip, TDip, taus, chi, type, level);
       LOp += complexi*LOpDFS;
       }
  }
else
  {
//std::cout << "\n\t\tadding internal pure dipolar"; 
//if(DDdfs) std::cout << " with DFS"; 
     REXijkl(LOp, sys, Ho, w, xiDs, xiDs, ADip, ADip,	// Dipole-dipole relaxation superop
             TDip, TDip, taus, chi, type, level, DDdfs);
  }
     }

// --------------------------------------------------------------------------------
//          Determine Chemical Shift Anisotropy Relaxation Superoperator
// --------------------------------------------------------------------------------

   if(CC)
     {
if(fext)
  {
std::cout << "\n\t\tadding pure CSA"; 
     Rij(LOp, sys, Ho, w, xiCs, xiCs, ACSA, ACSA,	// CSA-CSA relaxation superop
                  TCSA, TCSA, taus, chi, type, level);
     if(CCdfs)
       {
//std::cout << "\n\t\tadding CSA DFS"; 
       LOpDFS = LOp0;
       Rijds(LOpDFS, sys, Ho, w, xiCs,xiCs, ACSA, ACSA,	// CSA-CSA DFS relaxation superop
                  TCSA, TCSA, taus, chi, type, level);
       LOp += complexi*LOpDFS;
       }
  }
else
  {
//std::cout << "\n\t\tadding pure CSA"; 
//if(CCdfs) std::cout << " with DFS"; 
     REXij(LOp, sys, Ho, w, xiCs, xiCs, ACSA, ACSA,	// CSA-CSA relaxation superop
                  TCSA, TCSA, taus, chi, type, level, CCdfs);
  }
     }

// --------------------------------------------------------------------------------
//                Determine Quadrupolar Relaxation Superoperator
// --------------------------------------------------------------------------------

   if(QQ)
     {
if(fext)
  {
//std::cout << "\n\t\tadding pure Quad"; 
     Rij(LOp, sys, Ho, w, xiQs, xiQs, AQuad, AQuad,	// Quad-Quad relaxation superop
                  TQuad, TQuad, taus, chi, type, level);
     if(QQdfs)
       {
//std::cout << "\n\t\tadding Quad DFS"; 
       LOpDFS = LOp0;
       Rijds(LOpDFS, sys, Ho, w, xiQs, xiQs, AQuad, AQuad,	// Quad-Quad DFS relaxation superop
                  TQuad, TQuad, taus, chi, type, level);
       LOp += complexi*LOpDFS;
       }
  }
else
  {
//for(int a=0, aaa=0; a<hs; a++)
//for(int aa=0; aa<hs; aa++, aaa++)
//for(int b=0, bbb=0; b<hs; b++)
//for(int bb=0; bb<hs; bb++, bbb++)
//if(norm(LOp.get(aaa,bbb)) > 1.e-4)
//RQQel(sys, Ho, TQuad, TQuad, a, aa,b,bb);
//std::cout << "\n\t\tadding pure Quad"; 
//if(QQdfs) std::cout << " with DFS"; 
     REXij(LOp, sys, Ho, w, xiQs, xiQs, AQuad, AQuad,	// Quad-Quad relaxation superop
           TQuad, TQuad, taus, chi, type, level, QQdfs);
  }
     }

// --------------------------------------------------------------------------------
//            Determine Dipole & CSA Cross Correlation Superoperator
// --------------------------------------------------------------------------------

   if(DC)
     {
if(fext)
  {
std::cout << "\n\t\tadding Dipole-CSA"; 
     Rijk(LOpX, sys, Ho, w, xiDs, xiCs, ADip, ACSA,
                  TDip, TCSA, taus, chi, 0, level);
     if(level == 4)
       Rkij(LOpX, sys, Ho, w, xiCs, xiDs, ACSA, ADip,
                    TCSA, TDip, taus, chi, 0, level);
     if(DCdfs)
       {
std::cout << "\n\t\tadding Dipole-CSA DFS"; 
       Rijkds(LOpXdfs, sys, Ho, w, xiDs, xiCs, ADip, ACSA,
                       TDip, TCSA, taus, chi, 0, level);
       if(level == 4)
         Rkijds(LOpXdfs, sys, Ho, w, xiCs, xiDs, ACSA, ADip,
                         TCSA, TDip, taus, chi, 0, level);
       }
  }
else
  {
//std::cout << "\n\t\tadding Dipole-CSA"; 
//if(DCdfs) std::cout << " with DFS"; 
     REXijk(LOpX, sys, Ho, w, xiDs, xiCs, ADip, ACSA,
                    TDip, TCSA, taus, chi, 0, level, DCdfs);
     if(level == 4)
       REXkij(LOpX, sys, Ho, w, xiCs, xiDs, ACSA, ADip,
                    TCSA, TDip, taus, chi, 0, level, DCdfs);
  }
     }

// --------------------------------------------------------------------------------
//          Determine Dipole & Quadrupolar Cross Correlation Superoperator
// --------------------------------------------------------------------------------

   if(DQ)
     {
if(fext)
  {
//std::cout << "\n\t\tadding Dipole-Quad"; 
     Rijk(LOpX, sys, Ho, w, xiDs, xiQs, ADip, AQuad,
                  TDip, TQuad, taus, chi, 0, level);
     if(level == 4)
       Rkij(LOpX, sys, Ho, w, xiQs, xiDs, AQuad, ADip,
                    TQuad, TDip, taus, chi, 0, level);
     if(DQdfs)
       {
//std::cout << "\n\t\tadding Dipole-Quad DFS"; 
       Rijkds(LOpXdfs, sys, Ho, w, xiDs, xiQs, ADip, AQuad,
                    TDip, TQuad, taus, chi, 0, level);
       if(level == 4)
         Rkijds(LOpXdfs, sys, Ho, w, xiQs, xiDs, AQuad, ADip,
                      TQuad, TDip, taus, chi, 0, level);
       }
  }
else
  {
//std::cout << "\n\t\tadding Dipole-Quad"; 
//if(DQdfs) std::cout << " with DFS"; 
     REXijk(LOpX, sys, Ho, w, xiDs, xiQs, ADip, AQuad,
                  TDip, TQuad, taus, chi, 0, level, DQdfs);
     if(level == 4)
       REXkij(LOpX, sys, Ho, w, xiQs, xiDs, AQuad, ADip,
                    TQuad, TDip, taus, chi, 0, level, DQdfs);
  }
     }

// --------------------------------------------------------------------------------
//          Determine CSA & Quadrupolar Cross Correlation Superoperator
// --------------------------------------------------------------------------------

if(!fext)
{  if(QC)
     {
     REXij(LOpX, sys, Ho, w, xiCs, xiQs, ACSA, AQuad,
                  TCSA, TQuad, taus, chi, 0, level, QCdfs);
     if(level == 4)
       REXij(LOpX, sys, Ho, w, xiQs, xiCs, AQuad, ACSA,
                    TQuad, TCSA, taus, chi, 0, level, QCdfs);
     }
}
else
{  if(QC)
     {
     Rij(LOpX, sys, Ho, w, xiCs, xiQs, ACSA, AQuad,
                  TCSA, TQuad, taus, chi, 0, level);
     if(level == 4)
       Rij(LOpX, sys, Ho, w, xiQs, xiCs, AQuad, ACSA,
                    TQuad, TCSA, taus, chi, 0, level);
     if(QCdfs)
       {
       Rijds(LOpXdfs, sys, Ho, w, xiCs, xiQs, ACSA, AQuad,
                    TCSA, TQuad, taus, chi, 0, level);
       if(level == 4)
         Rijds(LOpXdfs, sys, Ho, w, xiQs, xiCs, AQuad, ACSA,
                      TQuad, TCSA, taus, chi, 0, level);
       }
    }
}


// --------------------------------------------------------------------------------
//               Scale The (Summed) Cross Correlation Superoperator
// --------------------------------------------------------------------------------

// Unless Level 4 Calculation Has Been Applied We Must Double The Cross-Correlation
// Superoperator.  For Level 4 LOpX = LOp(mu1,mu2) + LOp(mu2,mu1), Whereas In Other
// Computations We Need Just Use LOpX = 2 * LOp(mu1,mu2) = 2 * LOp(mu2,mu1).

   double fact = 1.0;					// Set Level 4 factor
   if(DC || DQ || QC)					// Only do this of cross
     {							// correlation is included
     if(level != 4) fact = 2.0;				// Set Level != 4 factor
     if((DCdfs || DQdfs || QCdfs) && fext)		// 1st add in dynamic shift
       LOpX += complexi*LOpXdfs;			// cross terms to superop
     LOp += fact*LOpX; 					// Now scale the superop
     }

//for(int a=0, aaa=0; a<hs; a++)
//for(int aa=0; aa<hs; aa++, aaa++)
//for(int b=0, bbb=0; b<hs; b++)
//for(int bb=0; bb<hs; bb++, bbb++)
//if(norm(LOp.get(aaa,bbb)) > 1.e-4)
//{
//if(DD)RDDel(sys, Ho, TDip, TDip, a, aa,b,bb);
//if(CC)RDDel(sys, Ho, TCSA, TCSA, a, aa,b,bb);
//if(QQ)RDDel(sys, Ho, TDip, TDip, a, aa,b,bb);
//}

//                 Now Compare with Imposed Secular Approximation

int debug=0;
if(debug)
{
std::cout << "\n\n\t\tTesting Oscillating Terms"; 
   if(sys.heteronuclear())
     {
     int ni = sys.isotopes();			// Isotope types in system
     gen_op X, *HZs, *FZs;			// For selective Zeeman, Fz's
     HZs = new gen_op[ni];
     FZs = new gen_op[ni];
     std::string *types;
     types = new std::string[ni];				// For isotope labels
     int I;
     for(I=0; I<ni; I++)			// Loop the isotope types
       {
       types[I] = sys.isotopes(I);		// Store isotope label
       X = Hz(sys,types[I]);			// Isotope selective Zeeman
       X.Op_base(Ho, 1.e-7);			// Put this into basis of Ho
       HZs[I] = X;				// Store isotope selective Zeeman
       X = Fz(sys,types[I]);			// Isotope selective Zeeman
       X.Op_base(Ho, 1.e-7);			// Put this into basis of Ho
       FZs[I] = X;				// Store isotope selective Zeeman
       }
     int hs = sys.HS();				// Get Hilbert space size
//int NOUT = 5*hs*hs;
int np = 0;
int full=0;
double sum1, sum2;
double x1, x2, x3;
     int aaa=0, bbb=0;
     double W1, W2, delW; 
     complex Rel;
     std::cout << "\n\tR ";
if(full) std::cout << "W1 W2 ";
     std::cout << "delW ";
     for(int a=0; a<hs; a++)			// Sum over transition a-aa
       for(int aa=0; aa<hs; aa++)
         {
         W1 = w[a] - w[aa];			// Transition a-aa frequency
         bbb = 0;
         for(int b=0; b<hs; b++)		// Sum over transition b-bb
           for(int bb=0; bb<hs; bb++)
             {
             W2 = w[b] - w[bb];		// Transition b-bb frequency
             delW = W1-W2;
             Rel = LOp.get(aaa,bbb);		// Get the current LOp element
//if(np<NOUT && norm(Rel)>1.e-4)
if(norm(Rel)>1.e-4)
{
             std::cout << "\n" << Rel;
if(full)
  {
             std::cout << "  " << W2;
             std::cout << "  " << W1;
  }
             std::cout << "  " << delW << "   ";
             sum1=0;
             for(I=0; I<ni; I++)
               {
               x1 = Re((FZs[I]).get(a,a)); 
               x2 = Re((FZs[I]).get(aa,aa));
               x3 = sys.Omega(types[I]);
//               x3 = Re((HZs[I]).get(a,aa));
               sum1 += (x1-x2)*x3; 
if(full)
  {
               std::cout << "  " << x1;
               std::cout << ", " << x2;
               std::cout << "  " << x3;
               std::cout << "  " << (x1-x2)*x3;;
  }
               }
             std::cout << "  " << sum1;
             sum2=0;
             for(I=0; I<ni; I++)
               {
               x1 = Re((FZs[I]).get(b,b)); 
               x2 = Re((FZs[I]).get(bb,bb));
               x3 = sys.Omega(types[I]);
//               x3 = Re((HZs[I]).get(b,bb));
               sum2 += (x1-x2)*x3; 
if(full)
  {
               std::cout << "  " << x1;
               std::cout << ", " << x2;
               std::cout << "  " << x3;
               std::cout << "  " << (x1-x2)*x3;;
  }
               }
             std::cout << "  " << sum2;
             std::cout << "  " << sum1-sum2;
delete [] HZs;
delete [] FZs;
delete [] types;
}
             bbb++;
np++;
             }
         aaa++;
         }
     }

}
   delete [] w;
   return LOp;
   }


/// sosixxxx


// --------------------------------------------------------------------------------
//                    Prepare Quadrupolar Relaxation Components
//	             (Only If Quadrupolar Effects Are Included)
// --------------------------------------------------------------------------------

//void WBRExch::prepCSA(const sys_dynamic& sys, matrix& Xis, spin_T* Ts, space_T* As)
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      sys   : Dynamic spin system
        //                      Xis   : Matrix for CSA interacation constants
	//			Ts    : NULL array for CSA spin tensor compts.
	//			As    : NULL array for CSA space tensor compts.
        // Output               void  : WBR Relaxation and Exchange
	//				superoperator

void WBRExch::prepQuad(const sys_dynamic& sys, matrix& Xis, spin_T* Ts, space_T* As) const
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      sys   : Dynamic spin system
        //                      Xis   : Matrix for Quad interacation constants
	//			Ts    : NULL array for Quad spin tensor compts.
	//			As    : NULL array for Quad space tensor compts.
        // Output               void  : WBR Relaxation and Exchange
	//				superoperator
 
   {
   int ns = sys.spins();			// Get number of spins
   Xis = xiQ(sys);				// Matrix of Xi's for quadrupoles
   Ts  = new spin_T[ns];			// Space for Quad spin tensors
   As  = new space_T[ns];			// Space for Quad spatial tensors
   for(int i=0; i<ns; i++)				// Set Quad applicable tensors
     {						//   Note that we skip any spins
     if(Xis.getRe(i,i))			//   having no quadrupole moment
       {
       As[i] = sys.TQ(i);
       Ts[i] = T_Q(sys,i);
       }
     }
   }


/// sosixxxx


 
super_op WBRExch::REXrf(const sys_dynamic& sys, gen_op& Heff, double Wrflab, int fext)
 
        // Input                WBRE  : Relaxation/exchange controls (this)
        //                      sys   : Dynamic spin system
        //                      Heff  : Effective static Hamiltonian
//                      fext    : Signals external functions
        // Output               LOp   : WBR Relaxation and Exchange
	//				superoperator
        // Note                         Computed in the eigenbasis of Heff
 
   {
   int hs = sys.HS();                           // Total system Hilbert space
   int ls = hs*hs;                              // Total system Liouville space
   Heff.set_EBR();				// Insure Heff is in its eigenbasis (Hz)
   matrix mx(ls, ls, complex0);			// Construct zero superoperator
   super_op LOp(mx, Heff.get_basis());
   super_op LOp0(mx, Heff.get_basis());
   super_op LOpDFS(mx, Heff.get_basis());
   super_op LOpX(mx, Heff.get_basis());
   super_op LOpXdfs(mx, Heff.get_basis());

//                      Prepare the Spectral Density Components
//        (5 Tau Values, A Chi Value, & Lab. Frame Transition Frequencies)

   double taus[5];                              // Get the 5 taus for spectral densities
   taust(taus, sys.taus());			// applicable to a rigid asymmetric top
   double chi = chit(sys.taus());               // Get the system chi value
   double *w;                                // Set up for system energy levels (LAB)
   w = new double[hs];
   if(abs(level) > 1)                           // Needed for higher level computations
     {                                          // where spectral densities are evaluated
// have to check this crap
     Heff.eigvals(w);                           // Frequencies w in field rotating frame
//     gen_op HZ = Hz(sys);			// Zeeman Hamiltonian (lab frame)
//     HZ.Op_base(Heff, 1.e-7);			// Put Hz into basis of Heff
//     if(!HZ.test_EBR()) WBRerror(60);		// Insure the result is proper
//     for(int k=0; k<hs; k++)			// Calculate the transition frequencies
//       w[k] = Re(HZ.get(k,k))+Re(Heff.get(k,k));// in the laboratory frame
     }   

//                  Prepare Dipolar Relaxation Components
//                 (Only If Dipolar Effects Are Included)
   
   int ns = sys.spins();			// Total number of spins
   int i, j;
   matrix xiDs; 				// Array for Dipolar int. constants
   spin_T *TDip = NULL;				// Spin tensors for each dipole
   space_T *ADip = NULL;			// Spatial tensors for each dipole
   if(DD)
     {
     xiDs = xiD(sys);				// Dipolar interaction constants
						// Prepare the Spin Tensors
     int ns = sys.spins();			// Total number of spins
     int ndip = sys.dipoles();			// Total number of dipoles
     TDip = new spin_T[ndip];			// Space for the dipolar spin tensors
     int ij=0;
     for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
       for(j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
         if(Re(xiDs.get(i,j)))
           {
           TDip[ij] = T_D(sys,i,j);
           ij++;
           }
     ADip = new space_T[ndip];			// Space for dipolar spatial tensors
     ij=0;
     for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
       for(j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
         if(Re(xiDs.get(i,j)))
           {
           ij++;
           }
     }

//                  Prepare CSA Relaxation Components
//              (Only If Shielding Effects Are Included)
/// sosixxxx

   matrix xiCs; 				// Matrix for CSA interaction constants
   spin_T *TCSA = NULL;				// Spin tensors for each spin
   space_T *ACSA = NULL;			// Spatial tensors for each spin
   if(CC)
     {
     xiCs = xiCSA(sys);				// Matrix of Xi's for CSA
						// Prepare CSA Spin Tensors
     TCSA = new spin_T[ns];			// Space for CSA spin tensors
     for(i=0; i<ns; i++)			// Set them all to CSA spin tensors
       if(Re(xiCs.get(i,i)))
         TCSA[i] = T_CS2(sys,i);
						// Prepare the CSA Spatial Tensors
     ACSA = new space_T[ns];			// Space for CSA spatial tensors
     for(i=0; i<ns; i++)			// Set them all to CSA space tensors
       if(Re(xiCs.get(i,i)))
         ACSA[i] = sys.TC(i);
     }


//                  Prepare Quadrupolar Relaxation Components
//                 (Only If Quadrupolar Effects Are Included)

   matrix xiQs;					// Matrix for quad interaction constants
   spin_T *TQuad = NULL;			// Spin tensors for each spin
   space_T *AQuad = NULL;			// Spatial tensors for each spin
   if(QQ)
     {
     xiQs = xiQ(sys);				// Matrix of Xi's for quadrupoles
     TQuad = new spin_T[ns];			// Space for Quad spin tensors
     AQuad = new space_T[ns];			// Space for Quad spatial tensors
     for(i=0; i<ns; i++)			// Set them all to Quad spin tensors
       if(xiQs.getRe(i,i))
         {
         AQuad[i] = sys.TQ(i);
         TQuad[i] = T_Q(sys,i);
         }
     }

//             Determine Dipolar Relaxation Superoperator
 
   if(DD)
     {
if(fext)
  {
std::cout << "\n\t\tadding pure dipolar"; 
     Rrfijkl(LOp, sys, Heff, w, Wrflab, xiDs, xiDs, ADip, ADip,	// Dipole-dipole relaxation superop
                  TDip, TDip, taus, chi, type, level);
     if(DDdfs)
       {
std::cout << "\n\t\tadding dipolar DFS"; 
       LOpDFS = LOp0;
       Rrfijklds(LOpDFS,sys,Heff,w,Wrflab,xiDs,xiDs,ADip,ADip,	// Dip-Dip DFS relaxation superop
                  TDip, TDip, taus, chi, type, level);
       LOp += complexi*LOpDFS;
       }
  }
else
  {
std::cout << "\n\t\tadding pure dipolar"; 
//if(DDdfs) std::cout << " with DFS";
     REXrfijkl(LOp, sys, Heff, w, Wrflab, xiDs, xiDs, ADip, ADip,	// Dipole-dipole relaxation superop
                      TDip, TDip, taus, chi, type, level, DDdfs);
  }
     }

//          Determine Chemical Shift Anisotropy Relaxation Superoperator

   if(CC)
     {
if(fext)
  {
std::cout << "\n\t\tadding pure CSA"; 
     Rrfij(LOp, sys, Heff, w, Wrflab, xiCs, xiCs, ACSA, ACSA,	// CSA-CSA relaxation superop
                  TCSA, TCSA, taus, chi, type, level);
     if(CCdfs)
       {
std::cout << "\n\t\tadding CSA DFS"; 
       LOpDFS = LOp0;
       Rrfijds(LOpDFS,sys,Heff,w,Wrflab,xiCs,xiCs,ACSA,ACSA,	// CSA-CSA DFS relaxation superop
                  TCSA, TCSA, taus, chi, type, level);
       LOp += complexi*LOpDFS;
       }
  }
else
  {
std::cout << "\n\t\tadding pure CSA"; 
//if(CCdfs) std::cout << " with DFS";
     REXrfij(LOp, sys, Heff, w, Wrflab, xiCs, xiCs, ACSA, ACSA,	// CSA-CSA relaxation superop
                  TCSA, TCSA, taus, chi, type, level, CCdfs);
  }
     }

//          Determine Quadrupolar Relaxation Superoperator

   if(QQ)
     {
if(fext)
  {
std::cout << "\n\t\tadding pure quadrupolar"; 
     Rrfij(LOp, sys, Heff, w, Wrflab, xiQs, xiQs, AQuad, AQuad,	// Quad-Quad relaxation superop
                  TQuad, TQuad, taus, chi, type, level);
     if(QQdfs)
       {
std::cout << "\n\t\tadding quad DFS"; 
       LOpDFS = LOp0;
       Rrfijds(LOpDFS,sys,Heff,w,Wrflab,xiQs,xiQs,AQuad,AQuad,	// Quad-Quad DFS relaxation superop
                  TQuad, TQuad, taus, chi, type, level);
       LOp += complexi*LOpDFS;
       }
  }
else
  {
std::cout << "\n\t\tadding pure quadrupolar"; 
//if(QQdfs) std::cout << " with DFS";
     REXrfij(LOp, sys, Heff, w, Wrflab, xiQs, xiQs, AQuad, AQuad,	// Quad-Quad relaxation superop
                  TQuad, TQuad, taus, chi, type, level, QQdfs);
  }
     }

//          Determine Dipole & CSA Cross Correlation Superoperator

   if(DC)
     {
if(fext)
  {
std::cout << "\n\t\tadding Dipole-CSA"; 
     Rrfijk(LOpX, sys, Heff, w, Wrflab, xiDs, xiCs, ADip, ACSA,
                  TDip, TCSA, taus, chi, 0, level);
     if(level == 4)
       Rrfkij(LOpX, sys, Heff, w, Wrflab, xiCs, xiDs, ACSA, ADip,
                    TCSA, TDip, taus, chi, 0, level);
     if(DCdfs)
       {
std::cout << "\n\t\tadding Dipole-CSA DFS"; 
       Rrfijkds(LOpXdfs, sys, Heff, w, Wrflab, xiDs, xiCs, ADip, ACSA,
                       TDip, TCSA, taus, chi, 0, level);
       if(level == 4)
         Rrfkijds(LOpXdfs, sys, Heff, w, Wrflab, xiCs, xiDs, ACSA, ADip,
                         TCSA, TDip, taus, chi, 0, level);
       }
  }
else
  {
std::cout << "\n\t\tadding Dipole-CSA"; 
//if(DCdfs) std::cout << " with DFS"; 
     REXrfijk(LOpX, sys, Heff, w, Wrflab, xiDs, xiCs, ADip, ACSA,
                  TDip, TCSA, taus, chi, 0, level, DCdfs);
     if(level == 4)
       REXrfkij(LOpX, sys, Heff, w, Wrflab, xiCs, xiDs, ACSA, ADip,
                    TCSA, TDip, taus, chi, 0, level, DCdfs);
  }
     }

//          Determine Dipole & Quadrupolar Cross Correlation Superoperator

   if(DQ)
     {
if(fext)
  {
std::cout << "\n\t\tadding Dipole-Quad"; 
     Rrfijk(LOpX, sys, Heff, w, Wrflab, xiDs, xiQs, ADip, AQuad,
                  TDip, TQuad, taus, chi, 0, level);
     if(level == 4)
       Rrfkij(LOpX, sys, Heff, w, Wrflab, xiQs, xiDs, AQuad, ADip,
                    TQuad, TDip, taus, chi, 0, level);
     if(DQdfs)
       {
std::cout << "\n\t\tadding Dipole-Quad DFS"; 
       Rrfijkds(LOpXdfs, sys, Heff, w, Wrflab, xiDs, xiQs, ADip, AQuad,
                    TDip, TQuad, taus, chi, 0, level);
       if(level == 4)
         Rrfkijds(LOpXdfs, sys, Heff, w, Wrflab, xiQs, xiDs, AQuad, ADip,
                      TQuad, TDip, taus, chi, 0, level);
       }
  }
else
  {
std::cout << "\n\t\tadding Dipole-Quad"; 
if(DQdfs) std::cout << " with DFS"; 
     REXrfijk(LOpX, sys, Heff, w, Wrflab, xiDs, xiQs, ADip, AQuad,
                  TDip, TQuad, taus, chi, 0, level, DQdfs);
     if(level == 4)
       REXrfkij(LOpX, sys, Heff, w, Wrflab, xiQs, xiDs, AQuad, ADip,
                    TQuad, TDip, taus, chi, 0, level, DQdfs);
  }
     }

//          Determine CSA & Quadrupolar Cross Correlation Superoperator

   if(QC)
     {
if(fext)
  {
std::cout << "\n\t\tadding CSA-Quad"; 
     Rrfij(LOpX, sys, Heff, w, Wrflab, xiCs, xiQs, ACSA, AQuad,
                  TCSA, TQuad, taus, chi, 0, level);
     if(level == 4)
       Rrfij(LOpX, sys, Heff, w, Wrflab, xiQs, xiCs, AQuad, ACSA,
                    TQuad, TCSA, taus, chi, 0, level);
     if(QCdfs)
       {
std::cout << "\n\t\tadding CSA-Quad DFS"; 
       Rrfijds(LOpXdfs, sys, Heff, w, Wrflab, xiCs, xiQs, ACSA, AQuad,
                    TCSA, TQuad, taus, chi, 0, level);
       if(level == 4)
         Rrfijds(LOpXdfs, sys, Heff, w, Wrflab, xiQs, xiCs, AQuad, ACSA,
                      TQuad, TCSA, taus, chi, 0, level);
       }
  }
else
  {
std::cout << "\n\t\tadding CSA-Quad"; 
if(QCdfs) std::cout << " with DFS";
     REXrfij(LOpX, sys, Heff, w, Wrflab, xiCs, xiQs, ACSA, AQuad,
                  TCSA, TQuad, taus, chi, 0, level, QCdfs);
     if(level == 4)
       REXrfij(LOpX, sys, Heff, w, Wrflab, xiQs, xiCs, AQuad, ACSA,
                    TQuad, TCSA, taus, chi, 0, level, QCdfs);
  }
     }

   double fact = 1.0;
   if(DC || DQ || QC)
     {
     if(level != 4) fact = 2.0;
     if(DCdfs || DQdfs || QCdfs)
       LOpX += complexi*LOpXdfs;
     LOp += fact*LOpX; 
     }
   delete [] w;
   return LOp;
   }


// ________________________________________________________________________________
//                              Class WBRExch LeftOvers 
// ________________________________________________________________________________

// KY - come back to this later - this function is killing the build 
//  void ask_relax(int argc, char* argv[], int& argn,
//                            super_op& R, const sys_dynamic& sys, gen_op& H)

//         // Input        argc    : Number of command line arguments
//         //              argv    : Command line arguments
//         //              argn    : Initial command line argument
//         //                        for relaxation parameters
//         //              R       : Relaxation superoperator
//         //              sys     : Dynamic spin system
//         //              H       : Hamiltonian
//         // Output       none    : Function is void.  The relaxation
//         //                        matrix R has different effects added
//         //                        in depending upon user requests

//     {
//     std::string dd="n", cc="n", dc="n", dddfs="n", ccdfs;
//     query_parameter(argc, argv, argn,                   // Dipolar relaxation
//      "\n\tInclude Dipolar Relaxation (y/n)? ", dd);
//     argn++;
//     if(dd=="y")
//       {
//       query_parameter(argc, argv, argn,			// Dipolar DFS relaxation effects
//          "\n\tInclude Dipolar Dynamic Frequency Shifts (y/n)? ", dddfs);
//       argn++;
//       }
//     query_parameter(argc, argv, argn,                   // Dipolar relaxation
//      "\n\tInclude CSA Relaxation (y/n)? ", cc);
//     argn++;
//     if(cc=="y")
//       {
//       query_parameter(argc, argv, argn,			// Dipolar DFS relaxation effects
//          "\n\tInclude CSA Dynamic Frequency Shifts (y/n)? ", ccdfs);
//       argn++;
//       }
//     if((dd=="y") && (cc=="y"))
//       {
//       query_parameter(argc, argv, argn,                   // Dipole-CSA relaxation
//            "\n\tInclude DCX Relaxation (y/n)? ", dc);
//       argn++;
//       }
//     if(dd == "y")
//       {
//       std::cout << "\n\tComputing The Dipolar Relaxation Matrix";
//       R += RDD(sys, H);
//       }
//     if(cc == "y")
//       {
//       std::cout << "\n\tComputing The CSA Relaxation Matrix";
//       R += RCC(sys, H);
//       }
//     if(dc == "y")
//       {
//       std::cout << "\n\tComputing The Dipole-CSA Relaxation Matrix";
//       R += RDCX(sys, H);
//       }
//     complex icmplx(0,1);                                // z = 0 + 1i
//     if(dddfs == "y")
//       {
//       std::cout << "\n\tComputing The Dipolar DFS Relaxation Matrix";
//       R += icmplx*RDDds(sys,H);
//       }
//     if(ccdfs == "y")
//       {
//       std::cout << "\n\tComputing The CSA DFS Relaxation Matrix";
//       R += icmplx*RCCds(sys,H);
//       }
//     return;
//     }

//*************************************************************************
//*************************************************************************
//*************************************************************************
//*************************************************************************
//*************************************************************************

// --------- Spin-Pair with Spin-Pair Interaction Functions -------------

void REXijkl(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi,
                                         int type, int level, int DFS)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			xi1s  : Interaction constant, mu1
	//			xi2s  : Interaction constant, mu2
	//			A1    : Spatial tensor, mu1
	//			A2    : Spatial tensor, mu2
	//			T1    : Spin tensor, mu1
	//			T2    : Spin tensor, mu2
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = Auto & Cross Correlation
	//				   + = Auto Correlation Only
	//				   - = Cross Correlation Only
	//			level : Relaxation treatment level
	//			DFS   : Dynamic frequency shift flag
	//				   + = Normal & DFS
	//				   0 = Normal, No DFS
	//				   - = Only DFS

/*	             Two Spin-Pair Rank 2 Mechanisms

                                              >=i     >=k
                      --- ---             --- --- --- ---
                      \   \               \   \   \   \
               LOp += /   /    R        = /   /   /   /   R
  	              --- ---   mu1,mu2   --- --- --- ---  ij,kl
                      mu1 mu2              i   j   k   l                 */

   {
   int het = sys.heteronuclear();
int rank=2;
matrix theta = sys.thetas();			// Get dipole theta values (radians) 
matrix phi = sys.phis(); 			// Get dipole phi values (radians)
double alphaij, betaij;			// Two Euler angles for ij dipole
double alphakl, betakl;			// Two Euler angles for kl dipole
   double xi1, xi2, xi1xi2=0;			// For specific interaction constants
   coord EA1, EA2;				// For specific Euler angles
   double c1s[5];				// Set up 5 coefficients interaction 1
   double c2s[5];				// Set up 5 coefficients interaction 2
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[5];
   T2s = new gen_op[5];
   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
   int m;					// z-component of ang. momentum
   int ij = 0;					// Sum over all dipolar pairs
   int kl = 0;					// and build relaxation matrix
//sosi - added for heteronuclear cutoff (not yet implemented)
   int ijhet=0;					// Flag if ij is heteronuclear
   int klhet=0;					// Flag if kl is heteronuclear
   for(int i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)
       {
//sosi - added for heteronuclear cutoff (not yet implemented) 
       ijhet = 0;
       if(sys.element(i) != sys.element(j))
         ijhet = 1;
       xi1 = Re(xi1s.get(i,j));			// Dipolar interaction constant i&j
alphaij = Re(phi.get(i,j));		// Alpha Euler angle: dipole ij
betaij = Re(theta.get(i,j));		// Beta Euler angle: dipole kl
//       Jcoeffs(c1s, alphaij, betaij, 0.0, chi);	// Set the 5 coefficients for ij
EA1.xyz(alphaij, betaij, 0.0);
         Jcoeffs(c1s, EA1, chi);		// Set the 5 coefficients for ij
       for(m=-2; m<3; m++)			// Put spin tensor for i,j into a
         {					// vector of operators in basis of Ho
         T1s[m+2] = gen_op(T1[ij].component(2,m));
         T1s[m+2].Op_base(Ho);
         }
       kl = 0;					// dipole kl count to zero
       for(int k=0; k<ns-1; k++)		// Sum spin pairs k&l: dipole kl
         for(int l=k+1; l<ns; l++)
           {
//sosi - added for heteronuclear cutoff (not yet implemented)
           klhet = 0;
           if(sys.element(k) != sys.element(l))
             klhet = 1;
//if(klhet != ijhet)
//std::cout << "\n\tSpin Pairs " << kl << " and " << ij 
//     << " Shouldn't Contribute!, one homo, one hetero";
//else if(ijhet)
//std::cout << "\n\tSpin Pairs " << kl << " and " << ij 
//     << " Maybe Shouldn't Contribute, both hetero!";
           if((ij == kl) && (type >= 0))	// Auto-correlation term
             {
             xi1xi2 = xi1*xi1;			// xi(ij)*xi(ij)
             if(abs(level) <2)
               {
               w0 = 0.0;			// Need only zero, single, & double
               w1 = sys.gamma(i)/GAMMA1H; 	// quantum transition frequencies
               w1 *= sys.Omega()*1.0e6;
               w2 = 2.0*w1;
                }
//sosi - new stuff
//             if(sys.element(i) == sys.element(j))
               REXmumu(LOp,T1s,T1s,w,hs,taus,
                     c1s,c1s,xi1xi2,w0,w1,w2,rank,level,1,DFS,het);
//             else
//               {
//               double RotFW = sys.Omega(j)*1.e6;
//               gen_op Fz1 = Fz(sys,j);
//               Fz1.Op_base(Ho, 1.e-7);			// Put Hz into basis of Ho
//               REXmumu(LOp,T1s,T1s,w,hs,taus,c1s,c1s,xi1xi2,
//                        w0,w1,w2,rank,level,1,DFS,het,Fz1,RotFW,Fz1,RotFW);
//               }
             }
           else if((ij != kl) && (type <= 0))	// Cross-correlation term
             {
             xi2 = Re(xi2s.get(k,l));		// Dipolar interaction constant k&l
             xi1xi2 = xi1*xi2;			// xi(ij)*xi(kl)
alphakl = Re(phi.get(k,l));	// Get the alpha Euler angle kl 
betakl = Re(theta.get(k,l));	// Get the beta Euler angle kl
//             Jcoeffs(c2s, alphakl, betakl, 	// Set the 5 coefficients for kl
//				   0.0,  chi);  // where gammakl = 0 dipolar
EA2.xyz(alphakl, betakl, 0.0);
             Jcoeffs(c2s, EA2, chi);		// Set the 5 coefficients for kl
             for(m=-2; m<3; m++)		// Put spin tensor for k,l into a
               {				// vector of operators in basis of Ho
               T2s[m+2] = gen_op(T2[kl].component(2,m));
               T2s[m+2].Op_base(Ho);
               }
             if(abs(level) <2)
               {
               wi = sys.gamma(i)/GAMMA1H;	// Need only zero, single and double
               wi *= sys.Omega()*1.0e6; 	// quantum transition frequencies
               wj = sys.gamma(j)/GAMMA1H;
               wj *= sys.Omega()*1.0e6; 	// !! Don't know which to use here!!
               w0 = wi-wj;
               w1 *= wi;
               w2 = wi+wj;
               }
// sosi - new stuff
//             if(!het)
               REXmumu(LOp,T1s,T2s,w,hs,taus,
                 c1s,c2s,xi1xi2,w0,w1,w2,rank,level,0, DFS,het);
//             else
//               {
//               gen_op Fz21 = Fz(sys,j);
//               gen_op Fz2 = Fz(sys,l);
//              double RotFW1 = 0;
//             double RotFW2 = 0;
//            if(sys.element(i) != sys.element(j))
//             {
//            RotFW1 = sys.Omega(j)*1.e6;
//           Fz1.Op_base(Ho, 1.e-7);			// Put Hz into basis of Ho
//          }
//       else RotFW1 = 0;
//               if(sys.element(k) != sys.element(l))
//                 {
//                 RotFW2 = sys.Omega(l)*1.e6;
//                 Fz2.Op_base(Ho, 1.e-7);			// Put Hz into basis of Ho
//                 }
//               REXmumu(LOp,T1s,T1s,w,hs,taus,c1s,c1s,xi1xi2,
//                 w0,w1,w2,rank,level,1,DFS,het,Fz1,RotFW1,Fz2,RotFW2);
//               }
             }
           kl++;				// Increment second dipole
           }
       ij++;					// Increment first dipole
       }
   gen_op Op;					// Psuedo destruction of operators
   if(T1s)
     {
     for(int ii=0; ii<5; ii++)
       T1s[ii] = Op;
     T1s = NULL;
     }
   if(T2s)
     {
     for(int jj=0; jj<5; jj++)
       T2s[jj] = Op;
     T2s = NULL;
     }
   return;
   if(A1 != NULL) type = 0;			// Compiler likes this used
   if(A2 != NULL) type = 0;			// Compiler likes this used
   }


// --------------- Spin-Pair with Spin Functions -------------------

void REXijk(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi,
                                         int type, int level, int DFS)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			xi1s  : Interaction constant, mu1
	//			xi2s  : Interaction constant, mu2
	//			A1    : Spatial tensor, mu1
	//			A2    : Spatial tensor, mu2
	//			T1    : Spin tensor, mu1
	//			T2    : Spin tensor, mu2
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = Auto & Cross Correlation
	//				   + = Auto Correlation Only
	//				   - = Cross Correlation Only
	//			level : Relaxation treatment level

/*	       Mechanism 1 = Spin-Pair; Mechanism 2 = Spin

                      --- ---             --- --- ---
                      \   \               \   \   \
               LOp += /   /    R        = /   /   /   R
  	              --- ---   mu1,mu2   --- --- ---  ij,k
                      mu1 mu2              i   j   k                     */

   {
   int het = sys.heteronuclear();
   int rank=2;
   matrix theta = sys.thetas();			// Dipole theta values (rads) 
   matrix phi = sys.phis(); 			// Dipole phi values (radians)
   double alphaij, betaij;			// Two Euler angles ij dipole

   double xi1, xi2, xi1xi2;			// Specific intera. constants
   coord EA2;					// For specific Euler angles
   double c1s[5];				// For 5 coeff. interaction 1
   double c2s[5];				// For 5 coeff. interaction 2
   gen_op *T1s;					// Compiler hated gen_op T1s[5]
   gen_op *T2s;					// For spin tensor components
   T1s = new gen_op[5];
   T2s = new gen_op[5];
   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
   int m;					// z-component of ang. momentum
   double cutoff = 0;
   int ij = 0;					// Sum over all dipolar pairs
   for(int i=0; i<ns-1; i++)			// Sum spin pairs i&j:dipole ij
     for(int j=i+1; j<ns; j++, ij++)
       {
       xi1 = Re(xi1s.get(i,j));			// Dipolar int. constant i&j
       alphaij = Re(phi.get(i,j));		// Alpha Euler angle: dipole ij
       betaij = Re(theta.get(i,j));		// Beta Euler angle: dipole kl
       Jcoeffs(c1s, alphaij, betaij, 0.0, chi);	// Set 5 coefficients for ij
       for(m=-2; m<3; m++)			// Put spin tensor for i,j into
         {					// vector of ops in basis of Ho
         T1s[m+2] = gen_op(T1[ij].component(2,m));
         T1s[m+2].Op_base(Ho);
         }
       for(int k=0; k<ns; k++)			// Sum over spins k
         {
         xi2 = xi2s.getRe(k,k);		// Spin k (mu2) int. constant
         xi1xi2 = xi1*xi2;			// xi(mu1,i)*xi(mu2,j)
         if(fabs(xi1xi2) > cutoff)		// Only add !=0 contributions
           {
           EA2 = (A2[k]).PASys_EA();		// Spin k (mu2) Euler angles
           Jcoeffs(c2s, EA2, chi);		// Set the 5 coefficients for k
           for(m=-2; m<3; m++)			// Put spin tensor for k into a
             {					// vector of ops in basis of Ho
             T2s[m+2] = gen_op(T2[k].component(2,m));
             T2s[m+2].Op_base(Ho);
             }
           if(abs(level) <2)
             {
             wi = sys.gamma(i)/GAMMA1H; 	// Need zero, single, & double
             wi *= sys.Omega()*1.0e6;	 	// quantum trans. frequencies
             wj = sys.gamma(j)/GAMMA1H;
             wj *= sys.Omega()*1.0e6;
             w0 = wi-wj;
             w1 *= wi;				// !! Which to use here??!!
             w2 = wi+wj;
             }
// sosi = the new stuff here
//           if(!het)
             REXmumu(LOp, T1s, T2s, w, hs, taus,
               c1s,c2s,xi1xi2,w0,w1,w2,rank,level,0,DFS,het);
//           else
//             {
//             double RotFW11 = sys.Omega(i)*1.e6;
//             gen_op Fz11 = Fz(sys,sys.symbol(i));
//             double RotFW12 = 0;
//             gen_op Fz12 = Fz11;
//             if(sys.symbol(i) != sys.symbol(j))
//               {
//               RotFW12 = sys.Omega(j)*1.e6;
//               Fz12 = Fz(sys,sys.symbol(j));
//               }
//             double RotFW21 = sys.Omega(k)*1.e6;
//             gen_op Fz21 = Fz(sys,sys.symbol(k));
//             double RotFW22 = 0;
//             gen_op Fz22 = Fz21;
//             Fz11.Op_base(Ho, 1.e-7);
//             Fz12.Op_base(Ho, 1.e-7);
//             Fz21.Op_base(Ho, 1.e-7);
//             Fz22.Op_base(Ho, 1.e-7);
//             REXmumu(LOp,T1s,T2s,w,hs,taus,c1s,c2s,xi1xi2,
//             w0,w1,w2,rank,level,1,DFS,het,
//             Fz11,RotFW11,Fz12,RotFW12,
//             Fz21,RotFW21,Fz22,RotFW22);
//             }
           }
         } 					// Increment spin k (mu2)
       }					// Increment pair ij (mu1)
   return;
   if(A1 != NULL) type = 0;			// Compiler likes these used
   }


void REXkij(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi,
                                         int type, int level, int DFS)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			xi1s  : Interaction constant, mu1
	//			xi2s  : Interaction constant, mu2
	//			A1    : Spatial tensor, mu1
	//			A2    : Spatial tensor, mu2
	//			T1    : Spin tensor, mu1
	//			T2    : Spin tensor, mu2
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = Auto & Cross Correlation
	//				   + = Auto Correlation Only
	//				   - = Cross Correlation Only
	//			level : Relaxation treatment level

/*	       Mechanism 1 = Spin; Mechanism 2 = Spin-Pair

                      --- ---             --- --- ---
                      \   \               \   \   \
               LOp += /   /    R        = /   /   /   R
  	              --- ---   mu1,mu2   --- --- ---  k,ij
                      mu1 mu2              k   i   j                     */

   {
   int het = sys.heteronuclear();
   int rank=2;
   matrix theta = sys.thetas();			// Get dipole theta values (radians) 
   matrix phi = sys.phis(); 			// Get dipole phi values (radians)
   double alphaij, betaij;			// Two Euler angles for ij dipole

   double xi1, xi2, xi1xi2;			// For specific interaction constants
   coord EA1 ;					// For specific Euler angles
   double c1s[5];				// Set up 5 coefficients interaction 1
   double c2s[5];				// Set up 5 coefficients interaction 2
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[5];
   T2s = new gen_op[5];
   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
   int m;					// z-component of ang. momentum
   int ij=0;					// z-component of ang. momentum
   double cutoff = 0;
   for(int k=0; k<ns; k++)			// Sum over spins k (mu1)
     {
     xi1 = Re(xi1s.get(k,k));			// Get spin k (mu1) interaction constant
     if(fabs(xi1) > cutoff)			// Only add non-trivial contributions
       {
       EA1 = (A1[k]).PASys_EA(); 		// Get spin k (mu1) space tensor Euler angles
       Jcoeffs(c1s, EA1, chi);			// Set the 5 coefficients for k
       for(m=-2; m<3; m++)			// Put spin tensor for i (mu1) into a
         {					// vector of operators in basis of Ho
         T1s[m+2] = gen_op(T1[k].component(2,m));
         T1s[m+2].Op_base(Ho);
         }
       ij = 0;					// Set dipole count to zero
       for(int i=0; i<ns-1; i++)		// Sum spin pairs i&j: dipole ij (mu2)
         for(int j=i+1; j<ns; j++)
           {
           xi2 = Re(xi2s.get(i,j));		// Dipolar interaction constant k&l
           xi1xi2 = xi1*xi2;			// xi(k)*xi(ij)
           alphaij = Re(phi.get(i,j));		// Alpha Euler angle: dipole ij
           betaij = Re(theta.get(i,j));		// Beta Euler angle: dipole kl
           Jcoeffs(c2s,alphaij,betaij,0.0,chi);	// Set the 5 coefficients for ij
           for(m=-2; m<3; m++)			// Put spin tensor for i,j into a
             {					// vector of operators in basis of Ho
             T2s[m+2] = gen_op(T2[ij].
                            component(2,m));
             T2s[m+2].Op_base(Ho);
             }
           if(abs(level) <2)
             {
             wi = sys.gamma(i)/GAMMA1H;	// Need only zero, single and double
             wi *= sys.Omega()*1.0e6;	 	// quantum transition frequencies
             wj = sys.gamma(j)/GAMMA1H;
             wj *= sys.Omega()*1.0e6;	 	// !! Don't know which to use here!!
             w0 = wi-wj;
             w1 *= wi;
             w2 = wi+wj;
             }
// sosi = the new stuff here
//           if(!het)
             REXmumu(LOp, T1s, T2s, w, hs, taus,
               c1s,c2s,xi1xi2,w0,w1,w2,rank,level,0,DFS,het);
//           else
//             {
//             double RotFW11 = sys.Omega(k)*1.e6;
//             gen_op Fz11 = Fz(sys,sys.symbol(k));
//             double RotFW12 = 0;
//             gen_op Fz12 = Fz11;
//             double RotFW21 = sys.Omega(i)*1.e6;
//             gen_op Fz21 = Fz(sys,sys.symbol(i));
//             double RotFW22 = 0;
//             gen_op Fz22 = Fz21;
//             if(sys.symbol(i) != sys.symbol(j))
//               {
//               RotFW22 = sys.Omega(j)*1.e6;
//               Fz22 = Fz(sys, sys.symbol(j));
//               }
//             Fz11.Op_base(Ho, 1.e-7);			// Put Hz into basis of Ho
//             Fz12.Op_base(Ho, 1.e-7);			// Put Hz into basis of Ho
//             Fz21.Op_base(Ho, 1.e-7);			// Put Hz into basis of Ho
//             Fz22.Op_base(Ho, 1.e-7);			// Put Hz into basis of Ho
//             REXmumu(LOp,T1s,T2s,w,hs,taus,c1s,c2s,xi1xi2,
//             w0,w1,w2,rank,level,1,DFS,het,
//             Fz11,RotFW11,Fz12,RotFW12,Fz21,RotFW21,Fz22,RotFW22);
//             }
           ij++;
           }					// Increment second dipole (mu2)
       }
     }						// Increment spin (mu1)
   gen_op Op;					// Psuedo destruction of operators
   if(T1s)
     {
     for(int ii=0; ii<5; ii++)
       T1s[ii] = Op;
     T1s = NULL;
     }
   if(T2s)
     {
     for(int jj=0; jj<5; jj++)
       T2s[jj] = Op;
     T2s = NULL;
     }
   return;
   if(A2 != NULL) type = 0;			// Compiler likes these used
   }



// -------------------- Spin with Spin Functions -------------------

void REXij(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi,
                                           int type, int level, int DFS)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			xi1s  : Interaction constant, mu1
	//			xi2s  : Interaction constant, mu2
	//			A1    : Spatial tensor, mu1
	//			A2    : Spatial tensor, mu2
	//			T1    : Spin tensor, mu1
	//			T2    : Spin tensor, mu2
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = Auto & Cross Correlation
	//				   + = Auto Correlation Only
	//				   - = Cross Correlation Only
	//			level : Relaxation treatment level
	//			DFS   : Dynamic frequency shift flag
	//				   + = Normal & DFS
	//				   0 = Normal, No DFS

/*	          Two Single Spin Rank 2 Mechanisms

                      --- ---             --- ---
                      \   \               \   \
               LOp += /   /    R        = /   /   R
  	              --- ---   mu1,mu2   --- ---  i,j
                      mu1 mu2              i   j                         */

   {
   int het = sys.heteronuclear();
   int rank=2;
   double cutoff = 1e-12;
   double xi1, xi2, xi1xi2;			// For specific interaction constants
   coord EA1, EA2;				// For specific Euler angles
   double c1s[5];				// Set up 5 coefficients interaction 1
   double c2s[5];				// Set up 5 coefficients interaction 2
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[5];
   T2s = new gen_op[5];
   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
   int m;					// z-component of ang. momentum
   for(int i=0; i<ns; i++)			// Sum over spins i (mu1)
     {
     xi1 = Re(xi1s.get(i,i));			// Get spin i (mu1) interaction constant
     if(fabs(xi1) > cutoff)			// Only add non-trivial contributions
       {
       EA1 = (A1[i]).PASys_EA();	 		// Get spin i (mu1) space tensor Euler angles
       Jcoeffs(c1s, EA1, chi);			// Set the 5 coefficients for i
       for(m=-2; m<3; m++)			// Put spin tensor for i (mu1) into a
         {					// vector of operators in basis of Ho
         T1s[m+2] = gen_op(T1[i].component(2,m));
         T1s[m+2].Op_base(Ho);
         }
       for(int j=0; j<ns; j++)			// Sum over spins j
         {
         if((i==j) && (type >= 0))		// Auto-correlation term
           {
           xi1xi2 = xi1*xi1;			// xi(mu1,i)*xi(mu1,i)
           if(abs(level) <2)
             {
             w0 = 0.0;				// Need only zero, single, & double
             w1 = sys.gamma(i)/GAMMA1H;	// quantum transition frequencies
             w1 *= sys.Omega()*1.0e6;
             w2 = 2.0*w1;
             }
           if(fabs(xi1xi2) > cutoff)		// Only add non-trivial contributions
             REXmumu(LOp, T1s, T1s, w, hs, taus,
                c1s,c1s,xi1xi2,w0,w1,w2,rank,level,1, DFS,het);
           }
         else if((i!=j) && (type<=0))		// Cross-correlation term
           {
           xi2 = Re(xi2s.get(j,j));		// Get spin j (mu2) interaction constant
           xi1xi2 = xi1*xi2;			// xi(mu1,i)*xi(mu2,j)
           if(fabs(xi1xi2) > cutoff)		// Only add non-trivial contributions
             {
             EA2 = (A2[j]).PASys_EA();		// Get spin j (mu2) space tensor Euler angles
             Jcoeffs(c2s, EA2, chi);		// Set the 5 coefficients for j
             for(m=-2; m<3; m++)		// Put spin tensor for j into a
               {				// vector of operators in basis of Ho
               T2s[m+2] = gen_op(T2[j].component(2,m));
               T2s[m+2].Op_base(Ho);
               }
             if(abs(level) <2)
               {
               wi = sys.gamma(i)/GAMMA1H; 	// Need only zero, single, & double
               wi *= sys.Omega()*1.0e6; 	// quantum transition frequencies
               wj = sys.gamma(j)/GAMMA1H;
               wj *= sys.Omega()*1.0e6;
               w0 = wi-wj;
               w1 *= wi;			// !! Don't know which to use here!!
               w2 = wi+wj;
               }
// sosi - current reasoning is that these terms never need a secular approx
het=0;
             REXmumu(LOp,T1s,T2s,w,hs,taus,
               c1s,c2s,xi1xi2,w0,w1,w2,rank,level,0,DFS,het);
             }
           }
         }					// Increment second spin
       }
     }						// Increment first spin
   return;
   }


// ---------------- Routing to Specific Level Computation ---------------
// ----------------- This Must Be For Homonuclear Systems ---------

   void REXmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w, int hs,
              double* taus, double* c1s, double* c2s, double xi1xi2,
         double w0, double w1, double w2, int l, int level, int autoc, int DFS, int het)

	// Input                LOp   : Superoperator
	//                      T1s   : Spin tensor components of mu1
	// 			T2s   : Spin tensor components of mu2
	//			w     : Transisiton frequencies
	//			hs    : Hilbert space size
	//			taus  : 5 effective correlation times
	//			c1s   : J coefficients of mu1
	//			c2s   : J coefficients of mu2
	//			xi1xi2: Interaction constant product 
	//			w0    : Zero quantum transition frequency
	//			w1    : Single quantum transition frequency
	//			w2    : Double quantum transition frequency
	//			l     : Rank of the relaxation mechanisms
	//			level : Relaxation treatment level
	//			autoc : Flag for auto correlation vs cross correlation
	//			DFS   : Dynamic frequency shift flag
	//				   + = Normal & DFS
	//				   0 = Normal, No DFS
	//				   - = Only DFS
	// Output		LOp   : Relaxation superoperator
	// Note			      :	Computed in the basis T1s and T2s

//                            LOp += R   (Level)
//                                    1,2  

   {
//sosi
w0=0;
w1=0;
w2=0;
autoc=0;
double secapp = 1.e6;
   matrix J12;					// Used for J's levels > 1
//   complex J0=complex0,J1=complex0,J2=complex0;	// Used for J's levels 0 & 1
   switch(level)
     {
     case 4:					// Level 4 mu1-mu2: element by element
       if(DFS >= 0)				// Get all reduced spectral densities
         {
         J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// These are the normal (Lorentzian) ones
         if(DFS)
           J12 += complexi*Q_reduced(w,hs,taus,c1s,c2s,1);	// These are the normal (Lorentzian) ones
         }
       else					// These are the imaginary DFS ones
         J12 = complexi*Q_reduced(w,hs,taus,c1s,c2s,1);
       J12 *= complex(xi1xi2);			// Scale by the interaction constants
// sosi - added this as a temporary patch to avoid the heteronuclear non-secular terms
       if(!het) REX_4(LOp, l, T1s, T2s, J12);
       else REX_3(LOp, w, l, T1s, T2s, J12, secapp); 
       break;
//     case -4:					// Level 4 mu1-mu2: double commutator
//       if(DFS >= 0)				// Get all reduced spectral densities
//         J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// These are the normal (Lorentzian) ones
//       if(DFS)					// These are the imaginary DFS ones
//         J12 += complexi*Q_reduced(w,hs,taus,c1s,c2s,1);
//       J12 *= complex(xi1xi2);			// Scale by the interaction constants
//       R_4s(LOp, l, T1s, T2s, J12);
//       break;
     case 3:					// Level 3 mu1-mu2: element by element
       if(DFS >= 0)				// Get all reduced spectral densities
         {
//std::cout << "\n\tComputing Level 3";
         J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// These are the normal (Lorentzian) ones
         if(DFS)
{
           J12 += complexi*Q_reduced(w,hs,taus,c1s,c2s,1);	// These are the normal (Lorentzian) ones
//std::cout << ", with DFS"; 
}
         }
       else					// These are the imaginary DFS ones
         J12 = complexi*Q_reduced(w,hs,taus,c1s,c2s,1);
       J12 *= complex(xi1xi2);			// Scale by the interaction constants
       REX_3(LOp, w, l, T1s, T2s, J12);
//       REX_3(LOp, w, l, T1s, T2s, J12, secapp);
       break;
//     case -3:					// Level 3 mu1-mu2: double commutator
//       if(DFS >= 0)				// Get all reduced spectral densities
//         J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// These are the normal (Lorentzian) ones
//       if(DFS)					// These are the imaginary DFS ones
//         J12 += complexi*Q_reduced(w,hs,taus,c1s,c2s,1);
//       J12 *= complex(xi1xi2);			// Scale by the interaction constants
//       R_3s(LOp, w, l, T1s, T2s, J12);
//       break;
//     case 2:					// Level 2 mu1-mu2: element by element
//       if(DFS >= 0)				// Get all reduced spectral densities
//         J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// These are the normal (Lorentzian) ones
//       if(DFS)					// These are the imaginary DFS ones
//         J12 += complexi*Q_reduced(w,hs,taus,c1s,c2s,1);
//       J12 *= complex(xi1xi2);			// Scale by the interaction constants
//       R_2(LOp, l, T1s, T2s, J12);
//       break;
//     case -2:					// Level 2 mu1-mu2: double commutator
//       if(DFS >= 0)				// Get all reduced spectral densities
//         J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// These are the normal (Lorentzian) ones
//       if(DFS)					// These are the imaginary DFS ones
//         J12 += complexi*Q_reduced(w,hs,taus,c1s,c2s,1);
//       J12 *= complex(xi1xi2);			// Scale by the interaction constants
//       R_2s(LOp, l, T1s, T2s, J12);
//       break;
//     case 1: 					// Level 1 mu1-mu2: double commutator
//       if(DFS >= 0)				// Get all reduced spectral densities
//         {
//         J0 = J_reduced(taus, c1s, c1s, w0, 1);	// J(0)
//         J1 = J_reduced(taus, c1s, c1s, w1, 1);	// J(w0)
//         J2 = J_reduced(taus, c1s, c1s, w2, 1);	// J(2w0)
//         }
//       if(DFS)					// These are the imaginary DFS ones
//         {
//         J0 += complexi*Q_reduced(taus, c1s, c1s, w0, 1);	// J(0)
//         J1 += complexi*Q_reduced(taus, c1s, c1s, w1, 1);	// J(w0)
//         J2 += complexi*Q_reduced(taus, c1s, c1s, w2, 1);	// J(2w0)
//         }
//       if(autoc)
//         R_AC_1(T1s, LOp, l,
//                     J0*xi1xi2,J1*xi1xi2,J2*xi1xi2);
//       else
//         R_CC_1(T1s, T2s, LOp, l,
//                   J0*xi1xi2,J1*xi1xi2,J2*xi1xi2);
//       break;
//     case 0:					// Level 0 mu1-mu2: element by element
//       if(DFS >=0)				// These are the imaginary DFS ones
//         J0 = J_reduced(taus,c1s,c2s,0.0,1);	// Need J(0) only, extreme narrowing
//       if(DFS)					// These are the imaginary DFS ones
//         J0 = complexi*Q_reduced(taus,c1s,c2s,0.0,1);	// Need J(0) only, extreme narrowing
//       if(fabs(xi1xi2*J0) > 1.e-6)
//         R_0(LOp,l,T1s,T2s,complex(xi1xi2*J0));
//       break;
//     default:					// Level 0 mu1-mu2: double commutator
//       if(DFS >=0)				// These are the imaginary DFS ones
//         J0 = J_reduced(taus,c1s,c2s,0.0,1);	// Need J(0) only, extreme narrowing
//       if(DFS)					// These are the imaginary DFS ones
//         J0 = complexi*Q_reduced(taus,c1s,c2s,0.0,1);	// Need J(0) only, extreme narrowing
//       if(fabs(xi1xi2*J0) > 1.e-6)
//         if(autoc)
//           R_AC_0(T1s, LOp, l, xi1xi2*J0);
//         else
//           R_CC_0(T1s,T2s,LOp,l,xi1xi2*J0);
//       break;
     }
   return;
   }


// ---------------- Routing to Specific Level Computation ---------------
// ------------------ This Is For Heteronuclear Systems -----------------

   void REXmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w, int hs,
              double* taus, double* c1s, double* c2s, double xi1xi2,
         double w0, double w1, double w2, int l, int level, int autoc, int DFS, int het,
                                        gen_op& Fz11, double W11, gen_op& Fz12, double W12,
                                        gen_op& Fz21, double W21, gen_op& Fz22, double W22)

	// Input                LOp   : Superoperator
	//                      T1s   : Spin tensor components of mu1
	// 			T2s   : Spin tensor components of mu2
	//			w     : Transisiton frequencies
	//			hs    : Hilbert space size
	//			taus  : 5 effective correlation times
	//			c1s   : J coefficients of mu1
	//			c2s   : J coefficients of mu2
	//			xi1xi2: Interaction constant product 
	//			w0    : Zero quantum transition frequency
	//			w1    : Single quantum transition frequency
	//			w2    : Double quantum transition frequency
	//			l     : Rank of the relaxation mechanisms
	//			level : Relaxation treatment level
	//			autoc : Flag for auto correlation vs cross correlation
	//			DFS   : Dynamic frequency shift flag
	//				   + = Normal & DFS
	//				   0 = Normal, No DFS
	//				   - = Only DFS
	//			m1    : Selective Fz values for T1s
	//			W1    : Rotating frame difference for T1's
	//			m2    : Selective Fz values for T2s
	//			W2    : Rotating frame difference for T2's
	// Output		LOp   : Relaxation superoperator
	// Note			      :	Computed in the basis T1s and T2s

//                            LOp += R   (Level)
//                                    1,2  

   {
w0=0;
w1=0;
w2=0;
autoc=0;
het=0;
   double secapp = 1.e6;
   matrix J12;					// Used for J's levels > 1
   switch(level)
     {
     case 4:					// Level 4 mu1-mu2: element by element
       if(DFS >= 0)				// Get all reduced spectral densities
         {
         J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// These are the normal (Lorentzian) ones
         if(DFS)
           J12 += complexi*Q_reduced(w,hs,taus,c1s,c2s,1);	// These are the normal (Lorentzian) ones
         }
       else					// These are the imaginary DFS ones
         J12 = complexi*Q_reduced(w,hs,taus,c1s,c2s,1);
       J12 *= complex(xi1xi2);			// Scale by the interaction constants
       REX_4(LOp,l,T1s,T2s,J12,Fz11,W11,Fz12,W12,
                               Fz21,W21,Fz22,W22);
       break;
// sosi the rest doesn't yet work anyway.
     case 3:					// Level 3 mu1-mu2: element by element
     default:
       if(DFS >= 0)				// Get all reduced spectral densities
         {
         J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// These are the normal (Lorentzian) ones
         if(DFS)
           J12 += complexi*Q_reduced(w,hs,taus,c1s,c2s,1);	// These are the normal (Lorentzian) ones
         }
       else					// These are the imaginary DFS ones
         J12 = complexi*Q_reduced(w,hs,taus,c1s,c2s,1);
       J12 *= complex(xi1xi2);			// Scale by the interaction constants
       REX_3(LOp, w, l, T1s, T2s, J12, secapp);
       break;
     }
   return;
   }


//*************************************************************************
//*************************************************************************
//*************************************************************************
//*************************************************************************
//*************************************************************************


// --------- Spin-Pair with Spin-Pair Interaction Functions -------------

   void REXrfijkl(super_op &LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			w     : Transition frequencies
	//			xi1s  : Interaction constant, mu1
	//			xi2s  : Interaction constant, mu2
	//			A1    : Spatial tensor, mu1
	//			A2    : Spatial tensor, mu2
	//			T1    : Spin tensor, mu1
	//			T2    : Spin tensor, mu2
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = Auto & Cross Correlation
	//				   + = Auto Correlation Only
	//				   - = Cross Correlation Only
	//			level : Relaxation treatment level
	//			DFS   : Dynamic frequency shift flag
	//				   + = Normal & DFS
	//				   0 = Normal, No DFS
	//				   - = Only DFS

/*	             Two Spin-Pair Rank 2 Mechanisms

                      --- ---             --- --- --- ---
                      \   \               \   \   \   \
               LOp += /   /    R        = /   /   /   /   R
  	              --- ---   mu1,mu2   --- --- --- ---  ij,kl
                      mu1 mu2              i   j   k   l                 */

   {
int het = sys.heteronuclear();
int rank=2;
matrix theta = sys.thetas();			// Get dipole theta values (radians) 
matrix phi = sys.phis(); 			// Get dipole phi values (radians)
double alphaij, betaij;			// Two Euler angles for ij dipole
double alphakl, betakl;			// Two Euler angles for kl dipole
   double xi1, xi2, xi1xi2=0;			// Specific interaction constants
   coord EA1, EA2;				// Specific Euler angles
   double c1s[5];				// Set up 5 coefficients interaction 1
   double c2s[5];				// Set up 5 coefficients interaction 2
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[5];
   T2s = new gen_op[5];
//   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
   int m;					// z-component of spin ang. momentum
   int ij = 0;                                  // Sum over all dipolar pairs
   int kl = 0;                                  // and build relaxation matrix

//
//		    Prepare the Interaction Constants
//		     and Spectral Density Components

  double mWrf=0;
  double J[3];
  matrix* J12 = NULL;                          // Matrices of spectral densities
  if(abs(level) > 1)                           // Needed for higher level computations
    {
    J12 = new matrix[5];
    Heff.eigvals(w);                           // Frequencies w in field rotating frame
    }

//	      Start the Relaxation Superoperator Computation

   for(int i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)
       {
       xi1 = Re(xi1s.get(i,j));			// 1st interaction constant, i&j
alphaij = Re(phi.get(i,j));		// Alpha Euler angle: dipole ij
betaij = Re(theta.get(i,j));		// Beta Euler angle: dipole kl
//      Jcoeffs(c1s, alphaij, betaij, 0.0, chi);	// Set the 5 coefficients for ij
EA1.xyz(alphaij, betaij, 0.0);
       Jcoeffs(c1s, EA1, chi);			// Set the 5 coefficients for ij
       for(m=-2; m<3; m++)			// Put spin tensor for i,j into a
         {					// vector of operators in basis of Fz
         T1s[m+2] = gen_op((T1[ij].component(2,m)));
         T1s[m+2].Op_base(Heff);
         }
       kl = 0;
       for(int k=0; k<ns-1; k++)		// Sum spin pairs k&l: dipole kl
         for(int l=k+1; l<ns; l++)
           {
           if((ij == kl) && (type >= 0))	// Auto-correlation term
             {
             xi1xi2 = xi1*xi1;			// xi(ij)*xi(ij)
             if(abs(level) > 1)
               for(m=-rank; m<=rank; m++)	// Set up all the reduced spectral
                 {				// density functions J(w-m*Wrf) 
                 mWrf = double(m)*Wrflab;
                 if(DFS>=0)
                   {
                   J12[m+2]=J_red_shft(w,mWrf,
 		             hs,taus,c1s,c1s,1);
                   if(DFS>0)
                     J12[m+2]+= complexi*Q_red_shft(w,mWrf,
 		                   hs,taus,c1s,c1s,1);
                   }
                 else
                   J12[m+2] = complexi*Q_red_shft(w,mWrf,
 		                   hs,taus,c1s,c1s,1);
                 J12[m+2] *= complex(xi1xi2);	// Scale by the interaction constants
                 }
//              else
//                 {				// density functions J(w-m*Wrf) 
//                 w0 = 0.0;			// For Level < 2, need only zero,
//                 w1 = sys.gamma(i)/GAMMA1H; 	// single, and double quantum
//                 w1 *= sys.Omega()*1.0e6; 	// transition frequencies
//                 J[0] = xi1xi2*J_reduced(taus,	// J(0)
//                                 c1s, c1s,w0,1);
//                 J[1] = xi1xi2*J_reduced(taus,	// J(w0)
// 				 c1s, c1s,w1,1);
//                 J[2] = xi1xi2*J_reduced(taus,	// J(2w0)
//         		      c1s,c1s,2.0*w2,1);
//                 }
             REXrfmumu(LOp,T1s,T1s,J12,J,w,rank,level,1, DFS, het);
             }
           else if((ij != kl) && (type <= 0))	// Cross-correlation term
             {
             xi2 = Re(xi2s.get(k,l));		// Dipolar interaction constant k&l
             xi1xi2 = xi1*xi2;			// xi(ij)*xi(kl)
alphakl = Re(phi.get(k,l));	// Get the alpha Euler angle kl 
betakl = Re(theta.get(k,l));	// Get the beta Euler angle kl
//            Jcoeffs(c2s, alphakl, betakl, 	// Set the 5 coefficients for kl
//			   0.0,  chi);  // where gammakl = 0 dipolar
             EA2.xyz(alphakl, betakl, 0.0);
             Jcoeffs(c2s, EA2, chi); 		// Set the 5 coefficients for kl
             for(m=-2; m<3; m++)		// Put spin tensor for k,l into a
               {				// vector of operators in basis of Fz
               T2s[m+2] = gen_op((T2[kl].component(2,m)));
               T2s[m+2].Op_base(Heff);
               if(abs(level) > 1)
                 {
                 mWrf = double(m)*Wrflab;
                 if(DFS>=0)
                   {
                   J12[m+2]=J_red_shft(w,mWrf,
 		             hs,taus,c1s,c2s,1);
                   if(DFS>0)
                     J12[m+2]+= complexi*Q_red_shft(w,mWrf,
 		                   hs,taus,c1s,c2s,1);
                   }
                 else
                   J12[m+2] = complexi*Q_red_shft(w,mWrf,
 		                   hs,taus,c1s,c2s,1);
                 J12[m+2] *= complex(xi1xi2);	// Scale by the interaction constants
                 }
               }
//             if(abs(level) < 2)
//               {
//               wi = sys.gamma(i)/GAMMA1H;	// Need only zero, single and double
//               wi *= sys.Omega()*1.0e6; 	// quantum transition frequencies
//               wj = sys.gamma(j)/GAMMA1H;
//               wj *= sys.Omega()*1.0e6; 	// !! Don't know which to use here!!
//               w0 = wi-wj;
//               w1 *= wi;
//               w2 = wi+wj;
//               J[0] = xi1xi2*J_reduced(taus,	// J(0)
//                                 c1s, c2s,w0,1);
//               J[1] = xi1xi2*J_reduced(taus,	// J(w0)
// 				 c1s, c2s,w1,1);
//               J[2] = xi1xi2*J_reduced(taus,	// J(2w0)
//         		      c1s,c2s,2.0*w2,1);
//               }
             REXrfmumu(LOp,T1s,T2s,J12,J,w,rank,level,0, DFS, het);
             }
           kl++;				// Increment second dipole
           }
       ij++;					// increment first dipole
       }

     gen_op Op;					// Psuedo destruction of operators
     if(T1s)
       {
       for(int ii=0; ii<5; ii++)
         T1s[ii] = Op;
       T1s = NULL;
       }
     if(T2s)
       {
       for(int jj=0; jj<5; jj++)
         T2s[jj] = Op;
        T2s = NULL;
       }
   return;
   if(A1==NULL && A2==NULL) type=0;		// Compiler likes these used
   }



// --------------- Spin-Pair with Spin Functions -------------------

void REXrfijk(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			w     : Transition frequencies
	//			Wrflab: Lab frame rf-field frequency
	//			xi1s  : Interaction constant, mu1
	//			xi2s  : Interaction constant, mu2
	//			A1    : Spatial tensor, mu1
	//			A2    : Spatial tensor, mu2
	//			T1    : Spin tensor, mu1
	//			T2    : Spin tensor, mu2
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = Auto & Cross Correlation
	//				   + = Auto Correlation Only
	//				   - = Cross Correlation Only
	//			level : Relaxation treatment level
	//			DFS   : Dynamic frequency shift flag
	//				   + = Normal & DFS
	//				   0 = Normal, No DFS
	//				   - = Only DFS

/*	       Mechanism 1 = Spin-Pair; Mechanism 2 = Spin

                      --- ---             --- --- ---
                      \   \               \   \   \
               LOp += /   /    R        = /   /   /   R
  	              --- ---   mu1,mu2   --- --- ---  ij,k
                      mu1 mu2              i   j   k                     */

   {
int het = sys.heteronuclear();
   int rank=2;
// sosi: these three variables are still dipole specific!!
matrix theta = sys.thetas();			// Get dipole theta values (radians) 
matrix phi = sys.phis(); 			// Get dipole phi values (radians)
double alphaij, betaij;				// Two Euler angles for ij dipole
   double xi1, xi2, xi1xi2=0;			// For specific interaction constants
   coord EA1, EA2;				// For specific Euler angles
   double c1s[5];				// Set up 5 coefficients interaction 1
   double c2s[5];				// Set up 5 coefficients interaction 2
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[5];
   T2s = new gen_op[5];
//   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
//   int ls = hs*hs;				// Total system Liouville space
   int m;					// z-component of ang. momentum
   double cutoff = 1e-12;
   double mWrf=0;
   double J[3];
   matrix* J12 = NULL;
   if(abs(level) > 1)                           // Needed for higher level computations
     {
     J12 = new matrix[5];
     Heff.eigvals(w);                           // Frequencies w in field rotating frame
     }
   int ij = 0;					// Sum over all dipolar pairs
   for(int i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)
       {
       xi1 = Re(xi1s.get(i,j));			// Dipolar interaction constant i&j
alphaij = Re(phi.get(i,j));		// Alpha Euler angle: dipole ij
betaij = Re(theta.get(i,j));		// Beta Euler angle: dipole kl
//       Jcoeffs(c1s, alphaij, betaij, 0.0, chi);	// Set the 5 coefficients for ij
EA1.xyz(alphaij, betaij, 0.0);
       Jcoeffs(c1s, EA1, chi);
       for(m=-2; m<3; m++)			// Put spin tensor for i,j into a
         {					// vector of operators in basis of Heff
         T1s[m+2] = gen_op(T1[ij].component(2,m));
         T1s[m+2].Op_base(Heff);
         }
       for(int k=0; k<ns; k++)			// Sum over spins k
         {
         xi2 = Re(xi2s.get(k,k));		// Get spin k (mu2) interaction constant
         xi1xi2 = xi1*xi2;			// xi(mu1,i)*xi(mu2,j)
         if(fabs(xi1xi2) > cutoff)		// Only add non-trivial contributions
           {
           EA2 = (A2[k]).PASys_EA();		// Get spin k (mu2) space tensor Euler angles
           Jcoeffs(c2s, EA2, chi);		// Set the 5 coefficients for k
           for(m=-2; m<3; m++)			// Put spin tensor for k into a
             {					// vector of operators in basis of Heff
             T2s[m+2] = gen_op(T2[k].component(2,m));
             T2s[m+2].Op_base(Heff);
             }
           if(abs(level) > 1)
             for(m=-rank; m<=rank; m++)		// Set up all the reduced spectral
               {				// density functions J(w-m*Wrf) 
               mWrf = double(m)*Wrflab;
               if(DFS >= 0)
                 {
                 J12[m+2]=J_red_shft(w,mWrf,
 		             hs,taus,c1s,c2s,1);
                 if(DFS > 0)
                   J12[m+2] += complexi*Q_red_shft(w,mWrf,
 		                 hs,taus,c1s,c2s,1);
                 }
               else
                 J12[m+2] = complexi*Q_red_shft(w,mWrf,
 		                 hs,taus,c1s,c2s,1);
               J12[m+2] *= complex(xi1xi2);	// Scale by the interaction constants
               }
//           else
//             {
//             wi = sys.gamma(i)/GAMMA1H; 	// Need only zero, single, & double
//             wi *= sys.Omega()*1.0e6;	 	// quantum transition frequencies
//             wj = sys.gamma(j)/GAMMA1H;
//             wj *= sys.Omega()*1.0e6;
//             w0 = wi-wj;
//             w1 *= wi;				// !! Don't know which to use here!!
//             w2 = wi+wj;
//             }
// sosi - added this to insure that heteronuclear terms are not included
//if(sys.element(i) == sys.element(k) && sys.element(j) == sys.element(k))
           REXrfmumu(LOp,T1s,T2s,J12,J,w,rank,level,0,DFS, het);
           }
         } 					// Increment spin k (mu2)
       }					// Increment pair ij (mu1)
   return;
   if(A1==NULL) type=0;				// Compiler likes these used
   }


void REXrfkij(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			w     : Transition frequencies
	//			Wrflab: Lab frame rf-field frequency
	//			xi1s  : Interaction constant, mu1
	//			xi2s  : Interaction constant, mu2
	//			A1    : Spatial tensor, mu1
	//			A2    : Spatial tensor, mu2
	//			T1    : Spin tensor, mu1
	//			T2    : Spin tensor, mu2
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = Auto & Cross Correlation
	//				   + = Auto Correlation Only
	//				   - = Cross Correlation Only
	//			level : Relaxation treatment level
	//			DFS   : Dynamic frequency shift flag
	//				   + = Normal & DFS
	//				   0 = Normal, No DFS
	//				   - = Only DFS

/*	       Mechanism 1 = Spin; Mechanism 2 = Spin-Pair

                      --- ---             --- --- ---
                      \   \               \   \   \
               LOp += /   /    R        = /   /   /   R
  	              --- ---   mu1,mu2   --- --- ---  k,ij
                      mu1 mu2              k   i   j                     */

   {
int het = sys.heteronuclear();
// sosi: these three variables are still dipole specific!!
   int rank=2;
matrix theta = sys.thetas();			// Get dipole theta values (radians) 
matrix phi = sys.phis(); 			// Get dipole phi values (radians)
double alphaij, betaij;			// Two Euler angles for ij dipole

   double xi1, xi2, xi1xi2;			// For specific interaction constants
   coord EA1, EA2;				// For specific Euler angles
   double c1s[5];				// Set up 5 coefficients interaction 1
   double c2s[5];				// Set up 5 coefficients interaction 2
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[5];
   T2s = new gen_op[5];
//   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
//   int ls = hs*hs;				// Total system Liouville space
   int m;					// z-component of ang. momentum
   int ij=0;					// z-component of ang. momentum
   double cutoff = 0;
//   double cutoff = 1e-12;
   double mWrf=0;
   double J[3];
   matrix* J12 = NULL;
   if(abs(level) > 1)                           // Needed for higher level computations
     {
     J12 = new matrix[5];
     Heff.eigvals(w);                           // Frequencies w in field rotating frame
     }
   for(int k=0; k<ns; k++)			// Sum over spins k (mu1)
     {
     xi1 = Re(xi1s.get(k,k));			// Get spin k (mu1) interaction constant
     if(fabs(xi1) > cutoff)			// Only add non-trivial contributions
       {
       EA1 = (A1[k]).PASys_EA(); 		// Get spin k (mu1) space tensor Euler angles
       Jcoeffs(c1s, EA1, chi);			// Set the 5 coefficients for k
       for(m=-2; m<3; m++)			// Put spin tensor for i (mu1) into a
         {					// vector of operators in basis of effH
         T1s[m+2] = gen_op(T1[k].component(2,m));
//         T1s[m+2] = gen_op((T1[k].component(2,m)).matrix());
//         T1s[m+2] = gen_op(T1[k].component(2,m),
//                               sys.get_basis());
         T1s[m+2].Op_base(Heff);
         }
       ij = 0;					// Set dipole count to zero
       for(int i=0; i<ns-1; i++)		// Sum spin pairs i&j: dipole ij (mu2)
         for(int j=i+1; j<ns; j++)
           {
           xi2 = Re(xi2s.get(i,j));		// Dipolar interaction constant k&l
           xi1xi2 = xi1*xi2;			// xi(k)*xi(ij)
alphaij = Re(phi.get(i,j));		// Alpha Euler angle: dipole ij
betaij = Re(theta.get(i,j));		// Beta Euler angle: dipole kl
//Jcoeffs(c2s,alphaij,betaij,0.0,chi);	// Set the 5 coefficients for ij
EA2.xyz(alphaij, betaij, 0.0);
           Jcoeffs(c2s,EA2,chi);	// Set the 5 coefficients for ij
           for(m=-2; m<3; m++)			// Put spin tensor for i,j into a
             {					// vector of operators in basis of Heff
             T2s[m+2] = gen_op(T2[ij].component(2,m));
             T2s[m+2].Op_base(Heff);
             }
           if(abs(level) > 1)
             for(m=-rank; m<=rank; m++)		// Set up all the reduced spectral
               {				// density functions J(w-m*Wrf) 
               mWrf = double(m)*Wrflab;
               if(DFS >= 0)
                 {
                 J12[m+2]=J_red_shft(w,mWrf,
 		             hs,taus,c1s,c2s,1);
                 if(DFS > 0)
                   J12[m+2]+=complexi*Q_red_shft(w,mWrf,
 		                hs,taus,c1s,c2s,1);
                 }
               else
                 J12[m+2]=complexi*Q_red_shft(w,mWrf,
 		                hs,taus,c1s,c2s,1);
               J12[m+2] *= complex(xi1xi2);	// Scale by the interaction constants
               }
//           else
//             {
//             wi = sys.gamma(i)/GAMMA1H;	// Need only zero, single and double
//             wi *= sys.Omega()*1.0e6;	 	// quantum transition frequencies
//             wj = sys.gamma(j)/GAMMA1H;
//             wj *= sys.Omega()*1.0e6;	 	// !! Don't know which to use here!!
//             w0 = wi-wj;
//             w1 *= wi;
//             w2 = wi+wj;
//             }
// sosi - added this to insure that heteronuclear terms are not included
//if(sys.element(i) == sys.element(k) && sys.element(j) == sys.element(k))
           REXrfmumu(LOp,T1s,T2s,J12,J,w,rank,level,0,DFS,het);
           ij++;
           }					// Increment second dipole (mu2)
       }
     }						// Increment spin (mu1)
   gen_op Op;					// Psuedo destruction of operators
   if(T1s)
     {
     for(int ii=0; ii<5; ii++)
       T1s[ii] = Op;
     T1s = NULL;
     }
   if(T2s)
     {
     for(int jj=0; jj<5; jj++)
       T2s[jj] = Op;
     T2s = NULL;
     }
   return;
   if(A2==NULL) type=0;				// Compiler likes these used
   }



// -------------------- Spin with Spin Functions -------------------

void REXrfij(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			w     : Transition frequencies
	//			xi1s  : Interaction constant, mu1
	//			xi2s  : Interaction constant, mu2
	//			A1    : Spatial tensor, mu1
	//			A2    : Spatial tensor, mu2
	//			T1    : Spin tensor, mu1
	//			T2    : Spin tensor, mu2
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = Auto & Cross Correlation
	//				   + = Auto Correlation Only
	//				   - = Cross Correlation Only
	//			level : Relaxation treatment level
	//			DFS   : Dynamic frequency shift flag
	//				   + = Normal & DFS
	//				   0 = Normal, No DFS
	//				   - = Only DFS

/*	          Two Single Spin Rank 2 Mechanisms

                      --- ---             --- ---
                      \   \               \   \
               LOp += /   /    R        = /   /   R
  	              --- ---   mu1,mu2   --- ---  i,j
                      mu1 mu2              i   j                         */

   {
int het = sys.heteronuclear();
   double xi1, xi2, xi1xi2;			// For specific interaction constants
   coord EA1, EA2;				// For specific Euler angles
   double c1s[5];				// Set up 5 coefficients interaction 1
   double c2s[5];				// Set up 5 coefficients interaction 2
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[5];
   T2s = new gen_op[5];
//   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
   int m;					// z-component of ang. momentum
   double cutoff = 1e-12;
   int rank = 2;
   double mWrf=0;
   double J[3];
   matrix* J12 = NULL;
   if(abs(level) > 1)                           // Needed for higher level computations
     {
     J12 = new matrix[5];
     Heff.eigvals(w);                           // Frequencies w in field rotating frame
     }

   for(int i=0; i<ns; i++)			// Sum over spins i (mu1)
     {
     xi1 = Re(xi1s.get(i,i));			// Get spin i (mu1) interaction constant
     if(fabs(xi1) > cutoff)			// Only add non-trivial contributions
       {
       EA1 = (A1[i]).PASys_EA();		// Get spin i (mu1) space tensor Euler angles
       Jcoeffs(c1s, EA1, chi);			// Set the 5 coefficients for i
       for(m=-2; m<3; m++)			// Put spin tensor for i (mu1) into a
         {					// vector of operators in basis of Heff
         T1s[m+2] = gen_op(T1[i].component(2,m));
         T1s[m+2].Op_base(Heff);
         }
       for(int j=0; j<ns; j++)			// Sum over spins j
         {
         if((i==j) && (type >= 0))		// Auto-correlation term
           {
           xi1xi2 = xi1*xi1;			// xi(mu1,i)*xi(mu1,i)
           if(abs(level) > 1)
             for(m=-rank; m<=rank; m++)		// Set up all the reduced spectral
               {				// density functions J(w-m*Wrf)
               mWrf = double(m)*Wrflab;
               if(DFS >= 0)
                 {
                 J12[m+2]=J_red_shft(w,mWrf,
                 hs,taus,c1s,c1s,1);
                 if(DFS)
                   J12[m+2] += complexi*Q_red_shft(w,mWrf,
                                 hs,taus,c1s,c1s,1);
                 }
               else 
               J12[m+2] = complexi*Q_red_shft(w,mWrf,
                                 hs,taus,c1s,c1s,1);
               J12[m+2] *= complex(xi1xi2);	// Scale by the interaction constants
               }
//           if(abs(level) <2)
//             {
//             w0 = 0.0;				// Need only zero, single, & double
//             w1 = sys.gamma(i)/GAMMA1H;	// quantum transition frequencies
//             w1 *= sys.Omega()*1.0e6;
//             J[0] = xi1xi2*J_reduced(taus,	// J(0)
//                                 c1s, c1s,w0,1);
//             J[1] = xi1xi2*J_reduced(taus,	// J(w0)
//                                 c1s, c1s,w1,1);
//             J[2] = xi1xi2*J_reduced(taus,	// J(2w0)
//                                 c1s,c1s,2.0*w2,1);
//             }
           if(fabs(xi1xi2) > cutoff)		// Only add non-trivial contributions
             REXrfmumu(LOp,T1s,T1s,J12,J,w,rank,level,1, DFS,het);
           }
         else if((i!=j) && (type<=0))		// Cross-correlation term
           {
           xi2 = Re(xi2s.get(j,j));		// Get spin j (mu2) interaction constant
           xi1xi2 = xi1*xi2;			// xi(mu1,i)*xi(mu2,j)
           if(fabs(xi1xi2) > cutoff)		// Only add non-trivial contributions
             {
             EA2 = (A2[j]).PASys_EA();		// Get spin j (mu2) space tensor Euler angles
             Jcoeffs(c2s, EA2, chi);		// Set the 5 coefficients for j
             for(m=-2; m<3; m++)		// Put spin tensor for j into a
               {				// vector of operators in basis of Heff
               T2s[m+2] = gen_op(T2[j].component(2,m));
               T2s[m+2].Op_base(Heff);
               if(abs(level) > 1)
                 {
                 mWrf = double(m)*Wrflab;
                 if(DFS >= 0)
                   {
                   J12[m+2]=J_red_shft(w,mWrf,
                   hs,taus,c1s,c2s,1);
                   if(DFS)
                     J12[m+2] += complexi*Q_red_shft(w,mWrf,
                                   hs,taus,c1s,c2s,1);
                   }
                 else 
                   J12[m+2] = complexi*Q_red_shft(w,mWrf,
                                 hs,taus,c1s,c2s,1);
                 J12[m+2] *= complex(xi1xi2);   // Scale by the interaction constants
                 }
               }
//             if(abs(level) <2)
//               {
//               wi = sys.gamma(i)/GAMMA1H; 	// Need only zero, single, & double
//               wi *= sys.Omega()*1.0e6; 	// quantum transition frequencies
//               wj = sys.gamma(j)/GAMMA1H;
//               wj *= sys.Omega()*1.0e6;
//               w0 = wi-wj;
//               w1 *= wi;			// !! Don't know which to use here!!
//               w2 = wi+wj;
//               J[0] = xi1xi2*J_reduced(taus,	// J(0)
//                                 c1s, c2s,w0,1);
//               J[1] = xi1xi2*J_reduced(taus,	// J(w0)
//                                 c1s, c2s,w1,1);
//               J[2] = xi1xi2*J_reduced(taus,	// J(2w0)
//                                 c1s,c2s,2.0*w2,1);
//               }
// sosi - added this to insure that heteronuclear terms are not included
//if(sys.element(i) == sys.element(j))
             REXrfmumu(LOp,T1s,T2s,J12,J,w,rank,level,0,DFS,het);
             }
           }
         }					// Increment second spin
       }
     }						// Increment first spin
   return;
   }



// ---------------- Routing to Specific Level Computation ---------------

   void REXrfmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, matrix* J12,
               double* J, double* w, int rank, int level, int autoc, int DFS, int het)

	// Input                LOp   : Superoperator
	//                      T1s   : Spin tensor components of mu1
	// 			T2s   : Spin tensor components of mu2
	//			J12   : Array of reduced spectral density
	//				matricies scaled by xi1*xi2
	//			J     : Vector containing J0, J1, & J2
	//				scaled by xi1*xi2
	//			rank  : Rank of the two interactions
	//			level : Relaxation treatment level
	//			autoc : Flag auto- vs cross- correlation
	// Output		void  : Relaxation superoperator LOp altered
        //                      DFS   : Dynamic frequency shift flag
        //                                 + = Normal & DFS
        //                                 0 = Normal, No DFS
        //                                 - = Only DFS

	// Note			      :	Computed in the basis T1s and T2s

//                            LOp += R   (Level)
//                                    1,2  
   {
double* X=J;
X=NULL;
autoc=0;
DFS=0;
//   double cutoff = 1.e-6;			// Cutoff for 0 spectral density
   double secapr = 1.e6; 			// Cutoff for secular approximation
   switch(level)
     {
     case 4:					// Level 4 mu1-mu2: element by element
     if(!het) REXrf_4(LOp, rank, T1s, T2s, J12);
     else REXrf_3(LOp, w, rank, T1s, T2s, J12, secapr);
       break;
//     case -4:					// Level 4 mu1-mu2: double commutator
//       Rrf_4s(LOp, rank, T1s, T2s, J12);
//       break;
     case 3:					// Level 3 mu1-mu2: element by element
       REXrf_3(LOp, w, rank, T1s, T2s, J12, secapr);
       break;
//     case -3:					// Level 3 mu1-mu2: double commutator
//       Rrf_3s(LOp, w, rank, T1s, T2s, J12);
//       break;
//     case 2:					// Level 2 mu1-mu2: element by element
//       Rrf_2(LOp, rank, T1s, T2s, J12);
//       break;
//     case -2:					// Level 2 mu1-mu2: double commutator
//       Rrf_2s(LOp, rank, T1s, T2s, J12);
//       break;
//     case 1: 					// Level 1 mu1-mu2: double commutator
//       if(autoc)
//         R_AC_1(T1s,LOp,rank,J[0],J[1],J[2]);
//       else
//         R_CC_1(T1s,T2s,LOp,rank,J[0],J[1],J[2]);
//       break;
//     case 0:					// Level 0 mu1-mu2: element by element
//         Rrf_0(LOp,rank,T1s,T2s,complex(J[0]));
//       break;
//     default:					// Level 0 mu1-mu2: double commutator
//       if(fabs(J[0]) > cutoff)
//         if(autoc)
//           R_AC_0(T1s, LOp, rank, J[0]);
//         else
//           R_CC_0(T1s,T2s,LOp,rank,J[0]);
//       break;
     }
   return;
   }



//*************************************************************************
//*************************************************************************
//*************************************************************************
//*************************************************************************
//*************************************************************************


// ________________________________________________________________________________
//                RELAXATION SUPEROPERATORS VIA ELEMENT CALCULATIONS
// ________________________________________________________________________________

// ---------------- Level 4 via Element Calculation ---------------------
// -------------- Single Rotating Frame (Homonuclear) -------------------

void REX_4(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12)

	// Input		LOp   : Current relaxation superoperator
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	// Return		void  : Level 4 Relaxation superop for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R4 
	//				          out	   in     12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases


//                  <a,a'| LOp |b,b'> += <a,a'| R4   |b,b'>
//		                                  1,2

// where a, a', b, b' are basis function indices.

   {
   int hs = T1s[0].dim();			// Get Hilbert space size
   int aaa=0, bbb=0;
   complex Rel;
   for(int a=0; a<hs; a++)			// Sum over transition a-aa
     for(int aa=0; aa<hs; aa++)
       {
       bbb = 0;
       for(int b=0; b<hs; b++)			// Sum over transition b-bb
         for(int bb=0; bb<hs; bb++)
           {
           Rel = LOp.get(aaa,bbb);		// Get the current LOp element
           Rel += REX_4(hs, T1s, T2s, J12,	// Add the R component
	                   rank, a, b, aa, bb);
           LOp.put(aaa, bbb, Rel);		// Put in the new element
           bbb++;
           }
       aaa++;
       }
   return;
   }


complex REX_4(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
		 	         int rank, int a, int b, int aa, int bb)

	// Input		hs    : Spin system Hilbert space
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//			rank  : Rank of the two interactions
	//			a, aa : 1st coherence (transition) indices
	//			b, bb : 2nd coherence (transition) indices
	// Output		Rel   : Level 4 Relaxation superoperator element
	//
	//				         <a,a'|R4  |b,b'>
	//					         12
	//				which will connect density operator elements
	//				<a|sigma|a'> and <b|sigma|b'> to each other.
	//
	// Note			      :	T1s, T2s, and LOp assumed in proper bases
	// Note			      : This routine uses NO secular approximation.
	//				It should be valid when both T1 & T2 involve
	//				spins in the same rotating frame. The T1
	//				spin(s) can be in a different rotating frame
	//				than the T2 spin(s) though.

//			The Following Formula is Implemented Herein
//			-------------------------------------------

//                     --- [ ---
//                     \   | \                 m       m 
// <a,a'|R4   |b,b'> = /   | /   delta     <a|T |g><b|T |g> J  (w  )
//	   1,2	       --- | ---      a',b'    1       2     12  gb
//		        m  [  g
//
//                               m        m                       m        m
//                         - <a|T |b><a'|T |b'> J  (w    ) - <b'|T |a'><b|T |a> J  (w  )
//	     	                 1        2      12  b'a'         1        2     12  ab
//
//                                                   ---                                     ]
//                                                   \               m        m              |
//                                                +  /   delta   <g|T |a'><g|T |b'> J  (b'g) |
//	     	                                     ---      a,b    1        2      12      |
//                                                    g                                      ]
//
// where index m sums over angular momentum components of the spin tensors, index g over basis
// functions, & spectral density function values, J, include the interaction constants xi1 & xi2.

// NOTE:  Indices a' & b' in the above equations are mapped int aa and bb in the code below.
//        It is valid only for interactions of spins existing in the same rotating frame. 

   {
   complex Rel = 0;
   complex J2 = J12.get(bb,aa);
   complex J3 = J12.get(a,b);
   int k=0, g=0;
   for(int m = -rank; m<=rank; m++)			// Sum over all tensor components
     {
     Rel-=J2*T1s[k].get(a,b)*conj(T2s[k].get(aa,bb));	// Add in terms RII
     Rel-=J3*T1s[k].get(bb,aa)*conj(T2s[k].get(b,a));	// Add in terms RIII
//     Rel-=J2*T1s[k].get(a,b)*T2s[k].get(aa,bb);		// Add in terms RII
//     Rel-=J3*T1s[k].get(bb,aa)*T2s[k].get(b,a);		// Add in terms RIII
     for(g=0; g<hs; g++)
       {
       if(aa == bb)					// Add in terms RI over gamma sum
         Rel += J12.get(g,b)*
		T1s[k].get(a,g)*conj(T2s[k].get(b,g));
//		T1s[k].get(a,g)*T2s[k].get(b,g);
       if(a == b) 					// Add in terms RIV over gamma sum
         Rel += J12.get(bb,g)*
              T1s[k].get(g,aa)*conj(T2s[k].get(g,bb));
//              T1s[k].get(g,aa)*T2s[k].get(g,bb);
       }
     k++;
     }
   return Rel;
   }


// ---------------- Level 4 via Element Calculation ---------------------
// ------------ Multiple Rotating Frame (Heteronuclear) -----------------


void REX_4(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12,
                      gen_op& Fz11, double W11, gen_op& Fz12, double W12,
                      gen_op& Fz21, double W21, gen_op& Fz22, double W22)

	// Input		LOp   : Current relaxation superoperator
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//			m1    : Selective Fz values for T1s
	//			W1    : Rotating frame difference for T1's
	//			m2    : Selective Fz values for T2s
	//			W2    : Rotating frame difference for T2's
	// Return		void  : Level 4 Relaxation superop for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R4 
	//				          out	   in     12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases


//                  <a,a'| LOp |b,b'> += <a,a'| R4   |b,b'>
//		                                  1,2

// where a, a', b, b' are basis function indices.

   {
double cutoff = 200.0;
   int hs = T1s[0].dim();			// Get Hilbert space size
   int aaa=0, bbb=0;
   complex Rel;
   double delwaa, delwbb; 
   for(int a=0; a<hs; a++)			// Sum over transition a-aa
     for(int aa=0; aa<hs; aa++)			//   (Rows of R: aaa)
       {
       delwaa = Re(Fz11.get(a,a)-Fz11.get(aa,aa))*W11
              + Re(Fz12.get(a,a)-Fz12.get(aa,aa))*W12;
       bbb = 0;
       for(int b=0; b<hs; b++)			// Sum over transition b-bb
         for(int bb=0; bb<hs; bb++)		//   (Cols of R: bbb)
           {
           delwbb = Re(Fz11.get(b,b)-Fz21.get(bb,bb))*W21
                  + Re(Fz12.get(b,b)-Fz22.get(bb,bb))*W22;
	   if(fabs(delwaa-delwbb) < cutoff)	// Apply secular approximation
             {
             Rel = LOp.get(aaa,bbb);		// Get the current LOp element
             Rel += REX_4(hs,T1s,T2s,J12,
	                        rank,a,b,aa,bb);
             LOp.put(aaa, bbb, Rel);		// Put in the new element
             }
//           Rel += REX_4(hs, T1s, T2s, J12,	// Add the R component
//	             rank, a, b, aa, bb,
//                             Fz1, W1, Fz2, W2);
//           LOp.put(aaa, bbb, Rel);		// Put in the new element
           bbb++;
           }
       aaa++;
       }
   return;
   }

complex REX_4(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
	                int rank, int a, int b, int aa, int bb,
                               gen_op& Fz1, double W1, gen_op& Fz2, double W2)

	// Input		hs    : Spin system Hilbert space
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//			rank  : Rank of the two interactions
	//			a, b  : 1st transition indices
	//			aa, bb: 2nd transition indices
	//			m1    : Selective Fz values for T1s
	//			W1    : Rotating frame difference for T1's
	//			m2    : Selective Fz values for T2s
	//			W2    : Rotating frame difference for T2's
	// Output		Rel   : Level 4 Relaxation superoperator element
	//
	//				         <a,a'|R4  |b,b'>
	//					         12
	//				which will connect density operator elements
	//				<a|sigma|a'> and <b|sigma|b'> to each other.
	//
	// Note			      :	T1s, T2s, and LOp assumed in proper bases
	// Note			      : This routine uses a MINIMAL secular approximation.
	//				It should be valid when both T1 & T2 involve
	//				spins in the same rotating frame. The T1
	//				spin(s) can be in a different rotating frame
	//				than the T2 spin(s) though.

//			The Following Formula is Implemented Herein
//			-------------------------------------------

//                     --- [ --- 
//                     \   | \                                 m       m 
// <a,a'|R4   |b,b'> = /   | /   del             delta     <a|T |g><b|T |g> J  (w  )
//	   1,2	       --- | ---    (ag,1),(gb,2)     a',b'    1       2     12  gb
//		        m  [  g
//
//                            m        m                                       m        m
//    - del               <a|T |b><a'|T |b'> J  (w    ) - del      del    <b'|T |a'><b|T |a> J  (w  )
//	   (ab,1),(a'b',2)    1        2      12  b'a'       b'a',1   ba,2     1        2     12  ab
//
//                                                   ---                                     ]
//                                                   \               m        m              |
//                                                +  /   delta   <g|T |a'><g|T |b'> J  (b'g) |
//	     	                                     ---      a,b    1        2      12      |
//                                                    g                                      ]
//
// where index m sums over angular momentum components of the spin tensors, index g over basis
// functions, & spectral density function values, J, include the interaction constants xi1 & xi2.
//
// Secular approximations are made when the Kronecker deltas are applied.  The Kronecker
// delta are expressed above in shorthand notation and shown expanded below.
//
//              del              = del
//	           (ag,1),(gb,2)      (m    a-m    g)*/\W , -(m    g-m    )*/\W
//                                      z,1S   z,1S   -- 1     z,2S   z,2S  -- 2
//
// where 1S indicates all spins in the second rotating frame of tensor component 1 and W
// indicates the difference in the two rotating frames associated with that tensor.     1
// If a particular tensor, say T1, were associated with only a single rotating frame
// then its "rotating frame difference frequency", W  should be zero.
//                                                  1

// NOTE:  Indices a' & b' in the above equations are mapped int aa and bb in the code below.
//        It is valid only for interactions of spins existing in multiple rotating frames
//	  assuming there is one rotating frame per isotope type.

   {
   complex Rel = 0;
   complex J2 = J12.get(bb,aa);
   complex J3 = J12.get(a,b);
   int k=0, g=0;
   int K2=0, K3=0;
   if(!(Re(Fz1.get(a,a)-Fz1.get(b,b))*W1-Re(Fz2.get(aa,aa)-Fz2.get(bb,bb))*W2))
     K2=1;
   if(!(Re(Fz1.get(bb,bb)-Fz1.get(aa,aa))*W1-Re(Fz2.get(b,b)-Fz2.get(a,a))*W2))
     K3=1;
   for(int m = -rank; m<=rank; m++)			// Sum over all tensor components
     {
     if(K2)
       Rel-=J2*T1s[k].get(a,b)*T2s[k].get(aa,bb);	// Add in terms RII
     if(K3)
       Rel-=J3*T1s[k].get(bb,aa)*T2s[k].get(b,a);	// Add in terms RIII
     for(g=0; g<hs; g++)
       {
       if(aa == bb)					// Add in terms RI over gamma sum
         {
         if(!(Re(Fz1.get(a,a)-Fz1.get(g,g))*W1-Re(Fz2.get(b,b)-Fz2.get(g,g))*W2))
           Rel += J12.get(g,b)*
		  T1s[k].get(a,g)*T2s[k].get(b,g);
         }
       if(a == b) 					// Add in terms RIV over gamma sum
         {
         if(!(Re(Fz1.get(g,g)-Fz1.get(aa,aa))*W1-Re(Fz2.get(g,g)-Fz2.get(bb,bb))*W2))
           Rel += J12.get(bb,g)*
              T1s[k].get(g,aa)*T2s[k].get(g,bb);
         }
       }
     k++;
     }
   return Rel;
   }

// ---------------- Level 3 via Element Calculation ---------------------
//                  (Applies Secular Approximation)

void REX_3(super_op& LOp, double* w, int rank, gen_op* T1s, gen_op* T2s,
					 matrix& J12, double cutoff)

	// Input		LOp   : Current relaxation superoperator
	//			w     : Vector of system energies
	//				in rad/sec in the LAB frame
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//			cutoff: Secular approximation cutoff value (Hz)
	// Return		void  : Level 3 Relaxation superop for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R3
	//				          out	   in     12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases


//   <a,a'| LOp |b,b'> += <a,a'| R3   |b,b'> == del         * <a,a'| R4   |b,b'>
//		                   1,2             w   ,w              1,2
//                                                  aa'  bb'

// where a, a', b, b' are basis function indices.

   {
   int hs = T1s[0].dim();			// Get Hilbert space size
   double delwaa, delwbb;
   int aaa=0, bbb=0;
   complex Rel;
   for(int a=0; a<hs; a++)			// Sum over transition a-aa
     for(int aa=0; aa<hs; aa++)
       {
       delwaa = w[a] - w[aa];			// Transition a-aa frequency
       bbb = 0;
       for(int b=0; b<hs; b++)			// Sum over transition b-bb
         for(int bb=0; bb<hs; bb++)
           {
           delwbb = w[b] - w[bb];		// Transition b-bb frequency
	   if(fabs(delwaa-delwbb) < cutoff)	// Apply secular approximation
             {
             Rel = LOp.get(aaa,bbb);
             Rel += REX_4(hs,T1s,T2s,J12,
	                        rank,a,b,aa,bb);
             LOp.put(aaa, bbb, Rel);		// Add to the relaxation matrix element
             }
           bbb++;
           }
       aaa++;
       }
   return;
   }

//*************************************************************************
//*************************************************************************
//*************************************************************************
//*************************************************************************
//*************************************************************************


// ----------------- Level 4 via Element Calculation --------------------

void REXrf_4(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix* J12)

	// Input		LOp   : Current relaxation superoperator
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st mu
	// 			T2s   : Spin tensor components, 2nd mu
	//			J12   : Matrices of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//				evaluated at m based frequencies
	//			rank  : Rank of the 2 interactions
	// Return		void  : Level 4 relax. superop for interactions
	//				1 & 2 added to input superoperator LOp
	//
	//				LOp    = LOp   + R4
	//				   out	    in     12
	//
	// Note				T1s, T2s, & LOp assumed in proper bases


	//                  <a,a'| LOp |b,b'> += <a,a'|R4   |b,b'>
	//                                               1,2
   
	// where a, a', b, b' are basis function indices.

   {
   int hs = T1s[0].dim();			// Get Hilbert space size
   int aaa=0, bbb;
   complex Rel;
   for(int a=0; a<hs; a++)			// Sum over transition a-aa
     for(int aa=0; aa<hs; aa++)
       {
       bbb = 0;
       for(int b=0; b<hs; b++)			// Sum over transition b-bb
         for(int bb=0; bb<hs; bb++)
           {
           Rel = LOp.get(aaa,bbb);		// Get the current LOp element
           Rel += REXrf_4(hs,T1s,T2s,J12,		// Add the R component
			       rank,a,b,aa,bb);
           LOp.put(aaa, bbb, Rel);		// Put in the new element
           bbb++;
           }
       aaa++;
       }
   return;
   }


complex REXrf_4(int hs, gen_op* T1s, gen_op* T2s, matrix* J12,
		 	         int rank, int a, int b, int aa, int bb)


	// Input		hs    : Spin system Hilbert space
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrices of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//				evaluated at m based frequencies
	//				times the two interaction constants
	//			rank  : Rank of the 2 interactions
	//			a, b  : 1st transition indices
	//			aa, bb: 2nd transition indices
	// Output		Rel   : Level 4 relaxation superop element
	//				for interactions 1 & 2 for component m
	//
	//				    <a,aa|R4  |b,bb>
	//					    12
	//
	// Note				T1s, T2s assumed in proper bases

//                    --- [ ---
//                    \   | \                 m       m
// <a,a'|R   |b,b'> = /   | /   delta     <a|T |g><b|T |g> J  (w  - mW  )
//        1,2         --- | ---      a',b'    1       2     12  gb    rf
//                     m  [  g
//
//                               m        m
//                         - <a|T |b><a'|T |b'> J  (w    - mW  )
//                               1        2      12  b'a'    rf
//
//                                         m        m
//                                  - <b'|T |a'><b|T |a> J  (w  - mW  )
//                                         1        2     12  ab    rf
//
//                            ---                                            ]
//                            \               m        m                     |
//                         +  /   delta   <g|T |a'><g|T |b'> J  (w   - mW  ) |
//                            ---      a,b    1        2      12  b'g    rf  |
//                             g                                             ]
   
// where a, a', b, b' and g are basis function indices.  The index m sums over
// components of angular momentum and the spectral densities are scaled by the
// two interaction constants.


   {
   complex Rel(0,0);
   int k=0, g=0;
   for(int m=-rank; m<=rank; m++)		// Sum over tensor components
     {
     Rel -= J12[k].get(bb,aa)*			// Add in terms RII
             T1s[k].get(a,b)*
              conj(T2s[k].get(aa,bb));
     Rel -= J12[k].get(a,b)* 			// Add in terms RIII
             T1s[k].get(bb,aa)*
              conj(T2s[k].get(b,a));
     for(g=0; g<hs; g++)
       {
       if(aa == bb)				// Add terms RI over gamma sum
         Rel += J12[k].get(g,b)*
		T1s[k].get(a,g)*
                conj(T2s[k].get(b,g));
       if(a == b) 				// Add terms RIV over gamma sum
         Rel += J12[k].get(bb,g)*
                T1s[k].get(g,aa)*
                conj(T2s[k].get(g,bb));
       }
     k++;
     }
   return Rel;
   }


// ---------------- Level 3 via Element Calculation ---------------------
//                  (Applies Secular Approximation)

void REXrf_3(super_op& LOp, double* w, int rank, gen_op* T1s,
                                         gen_op* T2s, matrix* J12, double cutoff)

	// Input		LOp   : Current relaxation superoperator
	//			w     : Vector system enery levels of Heff
	//				(rad/sec) in the LAB frame
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st mu
	// 			T2s   : Spin tensor components, 2nd mu
	//			J12   : Matrices of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//				evaluated at m based frequencies
        //                      cutoff: Secular approximation cutoff value (Hz)
	// Return		void  : Level 3 relax. superop for interactions
	//				1 & 2 added to input superoperator LOp
	//
	//				LOp    = LOp   + R3
	//				   out	    in     12
	//
	// Note				T1s, T2s, & LOp assumed in
	//				proper bases (eigenbasis of Heff)


// <a,a'| LOp |b,b'> += <a,a'| R3   |b,b'> == del           * <a,a'|R4   |b,b'>
//                               1,2             w    , w             1,2
//                                                a,a'   b,b'
   
// where a, a', b, b' are basis function indices.


   {
   int hs = T1s[0].dim();			// Get Hilbert space size
   double delwaa, delwbb;
   int aaa=0, bbb=0;
   complex Rel;
   for(int a=0; a<hs; a++)			// Sum over transition a-aa
     for(int aa=0; aa<hs; aa++)
       {
       delwaa = w[a] - w[aa];			// Transition a-aa frequency
       bbb = 0;
       for(int b=0; b<hs; b++)			// Sum over transition b-bb
         for(int bb=0; bb<hs; bb++)
           {
           delwbb = w[b] - w[bb];		// Transition b-bb frequency
	   if(fabs(delwaa-delwbb) < cutoff)	// Apply secular approximation
             {
             Rel = LOp.get(aaa,bbb);
             Rel += REXrf_4(hs,T1s,T2s,J12,
			       rank,a,b,aa,bb);
             LOp.put(aaa, bbb, Rel);		// Add to the relax. mx. element
             }
           bbb++;
           }
       aaa++;
       }
   return;
   }

#endif 							/* __RELAX_WBREx_CC__ */
