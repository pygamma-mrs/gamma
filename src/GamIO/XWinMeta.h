/* XWinMeta.h ***************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWinMeta                                   Interface		**
**                                                                      **
**      Copyright (c) 1999                                              **
**      Scott Smith                                                     **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**                                                                      **
**      $Header: 
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**      Description                                                     **
**                                                                      **
** The XWin* files provide an interface to Bruker XWinNMR (uxnmr) data  **
** sets. This class embodies a Bruker parameter file, meta, which is    **
** uset to control plotting and output of NMR spectra within XWinNMR.	**
** Yes, there just are not enough parameter files associated with a 	**
** spectrum so this joins the list.  Its an ASCII file, and has nothing	**
** to do with GAMMA. In fact, it will only be used output in a crude	**
** sense. Such meta files are loaded with parameters which are not	**
** particularly useful to GAMMA. We are not needed when reading XWinNMR	**
** data and there will be nothing changable in meta files from within	**
** GAMMA when they are written. In fact, were it not for XWinNMR 	**
** complaining about their absense, I wouldn't deal with them at all.	**
** But, being forced to, GAMMA will just write the same damn meta file	**
** out whenever and wherever XWinNMR wants one. Below is the XWinNMR	**
** 1D data set directory hierarchy:					**
**									**
**                          __ acqu  (changable parameter file)		**
**			   / 						**
**                        /___ acqus (static parameter file)		**
**			 /						**
**  expname -- expnum --< ---- fid (binary data)			**
**			 \						**
**			  \___ pdata -- 1 -- proc, procs, meta, outd	**
**									**
*************************************************************************/


#ifndef   XWinMeta_H_ 				// Is file already included?
#  define XWinMeta_H_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the implementation
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <string>				// Include libstdc++ strings
#include <GamIO/XWinPSet.h>			// Include Bruier parameter sets

class XWinMeta : public XWinPSet

// ----------------------------------------------------------------------------
// ------------------------------- STRUCTURE ----------------------------------
// ----------------------------------------------------------------------------
 
  {
  std::string  mname;			// Parameter file name
  int     oldflag;			// Flag for old format 
  std::string _ANGLE1;
  std::string _ANGLE2;
  std::string _ANGLE3;
  std::string _ASPA;
  std::string _AXCELL;
  std::string _AXDLBL;
  std::string _AXDTIC;
  std::string _AXGDIS;
  std::string _AXGLEN;
  std::string _AXLABL;
  std::string _AXLDIR;
  std::string _AXLIM1;
  std::string _AXLIM2;
  std::string _AXLOOP;
  std::string _AXFLAG;
  std::string _AXFORM;
  std::string _AXSTRG;
  std::string _ARWIDTH;
  std::string _ARHEIGHT;
  std::string _CONTRI;
  std::string _D_TYPE;
  std::string _DISABLE;
  std::string _FILENAMES;
  std::string _INTENSITY;
  std::string _LINEBF;
  std::string _HUE;
  std::string _MAGIC;
  std::string _MODIFY;
  std::string _OBJECTS;
  std::string _OWNER;
  std::string _PENCOL;
  std::string _PIC_ANZ;
  std::string _PICNUMB;
  std::string _PICS;
  std::string _PLANENR;
  std::string _PLDHEI;
  std::string _PLDIWI;
  std::string _PLLABF;
  std::string _PLMSHI;
  std::string _PLROTA;
  std::string _PLSTRG;
  std::string _SATURATION;
  std::string _SCALFLG;
  std::string _SHADOW;
  std::string _SOFTCLIP;
  std::string _TXMODE;
  std::string _TXALLI;
  std::string _TXFNST;
  std::string _TXSOSP;
  std::string _TXCELL;
  std::string _TYPNAM;
  std::string _TYPK;
  std::string _TYPNR;
  std::string _XHIGHEST; 
  std::string _XLEAST;
  std::string _XLENGTH; 
  std::string _XORIGIN; 
  std::string _XSCALE; 
  std::string _YHIGHEST;
  std::string _YLEAST;
  std::string _YLENGTH;
  std::string _YORIGIN; 
  std::string _YSCALE; 
  std::string _ZHIGHEST;
  std::string _ZLEAST;
  std::string _ZORIGIN;
  std::string _ZSCALE;


 
// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________       
// i                      XWinNMR Meta File Error Handling
// ____________________________________________________________________________

/* These functions take care of any errors encountered when reading, writing,   
   and setting parameters in Bruker acquisition parameter files.
 
	Input		AcqPar  : XWinNMR acqusition parameters (this)
			eidx    : Error index
                        pname   : Additional error message
                        noret   : Flag for linefeed (0=linefeed)
        Output          void    : An error message is output                 */

         void XWinMetaerror(int    eidx, int noret=0) const;
         void XWinMetaerror(int    eidx, const std::string& pname, int noret=0) const;
volatile void XWinMetafatality(int eidx) const;
volatile void XWinMetafatality(int eidx, const std::string& pname) const;

 
// ____________________________________________________________________________
// ii                    XWinNMR Meta File Setup Functions
// ____________________________________________________________________________

void SetBase(int old);
void SetAxis(int axis, int old);
void SetLine(int axis, int old);
void SetPic(int axis,  int old);
void SetSpec(int axis, int old);
void SetText(int axis, int old);
void SetPeak(int pk,   int old);
void SetImag(int pk,   int old);

public:
// ____________________________________________________________________________ 
// A             XWinMeta Parameter File Constructors, Destructor
// ____________________________________________________________________________
 
/* These are the constructors of the class handling Bruker XWinNMR output
   parameter files.  This doesn't do anything in particular, it is the write
   functions that perform the work.                                          */

XWinMeta();
XWinMeta(const std::string& name);
virtual ~XWinMeta();                                                             

// ____________________________________________________________________________
// B                  XWinMeta Parameter Access Functions
// ____________________________________________________________________________
 
/* These functions allow users to set any more important parameters that 
   these files contain.                                                      */

void OldFlag(int of);
 
//             WE DONT CARE ABOUT ANY OF THESE (at least not yet)

// ____________________________________________________________________________
// C                        XWinMeta Input Functions
// ____________________________________________________________________________

/* These functions will read in the 1D print/plot parameters from an XWinNMR
   parameter file, typically named meta.  By design, the Bruker parameter
   file is initially read into a GAMMA parameter set so that ALL parameters in
   the file are stored (class XWinPSet).  Subsequently, the parameters in the
   parameter set are parsed (herein) to obtain values of consequence to GAMMA
   and these are explicitly maintained variables in this class.

    Function                               Purpose
  ____________________________________________________________________________
     read              Read in parameter set for class object.  This is done
                       special since Bruker format is NOT in GAMMA parameter
                       format so it must be parsed appropriately.
    parsePSet          Converts parameters in internal pset to specific values
     getPar            Returns parameter found in the parameter set herein
                       This function is inherited from base class XWinPSet. */

//virtual int read(const string& filein, int warn=1);
//virtual int read(int warn=1);
//        int parsePSet(int warn=1);
//int getPar(const string& n,string& v,int i=-1, int w=0) const;

//        WE DONT WANT TO INPUT ANY OF THIS STUFF (at least not yet)

// ____________________________________________________________________________
// D                       XWinMeta Output Functions
// ____________________________________________________________________________

/* These function allow for output of NMR parameters directly into a Bruker
   XWinNMR ASCII parameter file (meta).                                      */
 
int write(const std::string& name, int warn=2, int of=-1);
int write(int warn=2);
void writeBase(std::ofstream& ofstr) const;
void writeInit(std::ofstream& ofstr,  int ef=1) const;
void writeLutab(std::ofstream& ofstr, int ef=1) const;
void writeLine(std::ofstream& ofstr,  int ef=1) const;
void writeAxis(std::ofstream& ofstr,  int ef=1) const;
void writePic(std::ofstream& ofstr,   int ef=1) const;
void writeSpec(std::ofstream& ofstr,  int ef=1) const;
void writeText(std::ofstream& ofstr,  int ef=1) const;
void writePeak(std::ofstream& ofstr,  int ef=1) const;
void writeImag(std::ofstream& ofstr,  int ef=1) const;
void writeFirst(std::ofstream& ofstr)           const;
void writeDraw(std::ofstream& ofstr)            const;
void writeXYZ(std::ofstream& ofstr)             const;



// ____________________________________________________________________________
// E                    XWinMeta Standard Output Functions
// ____________________________________________________________________________

/* These function allow for a quick output of the parameter file contents.
   They don't have anything to do with output while running XWinNMR, rather
   users can just glance at meta parameters or store then in a small file.   */

//ostream& printPset(ostream& O) const;              INHERITED
//virtual ostream& print(ostream& ostr, int full=0, int hdr=1) const;
//friend ostream& operator<< (ostream& O, const XWinMeta& P);
};

#endif 								// XWinMeta.h

                                                                                



