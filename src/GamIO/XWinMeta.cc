/* XWinMeta.cc **************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWinMeta                                  Implementation	**
**                                                                      **
**      Copyright (c) 1999                                              **
**      Scott Smith                                                     **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/
/*************************************************************************
**                                                                      **
**      Description                                                     **
**                                                                      **
** The XWin* files provide an interface to Bruker XWinNMR (uxnmr) data  **
** sets. This class embodies a Bruker parameter file, meta, which seems **
** to control access to NMR parameters within XWinNMR. Yet another      **
** ASCII file, this has nothing to do with GAMMA and will only be used  **
** for output. Typically this file will be dir/1/pdata/meta where dir   **
** is the directory base that contains XWinNMR data.  Note that without **
** the existence of an meta file one cannot display status parameters   **
** under the "Output" menu in XWinNMR.                                  **
**                                                                      **
*************************************************************************/

#ifndef   XWinMeta_CC_			// Is file already included?
#  define XWinMeta_CC_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <GamIO/XWinMeta.h>		// Include the interface
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <string>			// Include libstdc++ strings
#include <iostream>			// Include libstdc++ io streams
#ifndef _MSC_VER                        // If not using MSVC++ then we
 #include <sys/time.h>			// Include time and date access
 #include <unistd.h>			// Include POSIX getcwd function
#else                                   // and if using MSVC++ then I guess
 #include <time.h>			// Include time and date access
#endif

using std::string;			// Using libstdc++ strings
using std::ofstream;			// Using libstdc++ output file streams

// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      XWinNMR Meta File Error Handling
// ____________________________________________________________________________
 
/* These functions take care of any errors encountered when reading, writing,
   and setting parameters in Bruker acquisition parameter files.

        Input		eidx    : Error index
        		pname   : Additional error message
        		noret	: Flag for linefeed (0=linefeed)
        Output		void    : An error message is output                 */
 
void XWinMeta::XWinMetaerror(int eidx, int noret) const
  {
  string hdr("XWinNMR Output Parameter File");
  string msg;
  switch (eidx)
    {
    case 25: GAMMAerror(hdr,"Cannot Write Parameters To File",noret);break; // (25)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }  

    
void XWinMeta::XWinMetaerror(int eidx, const string& pname, int noret) const
  {                                                                             
  string hdr("XWinNMR Output Parameter File");
  string msg;
  switch(eidx)
    {
    case 21:msg = string("Cannot Write To ") + pname;
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (21)
    default: GAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error  (-1)
    }
  }

volatile void XWinMeta::XWinMetafatality(int eidx) const
  {                                                                 
  XWinMetaerror(eidx, 1);			// Normal non-fatal error
  if(eidx) XWinMetaerror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program                          // Quit program
  }

volatile void XWinMeta::XWinMetafatality(int eidx, const string& pname) const
  {                                                                 
  XWinMetaerror(eidx, pname, 1);		// Normal non-fatal error
  if(eidx) XWinMetaerror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }
 
// ____________________________________________________________________________
// ii                    XWinNMR Meta File Setup Functions
// ____________________________________________________________________________

 
void XWinMeta::SetBase(int old)
  {
  switch(old)
    {
    default:
    case 0:					// This is for newer meta files
      _OWNER    = string("sosi");
      _MAGIC    = string("0x8001");
      _OBJECTS  = "18";
      _PICS     = "7";
      _ASPA     = "51";
      _ARWIDTH  = "278.57144";
      _ARHEIGHT = "180";
      _XLENGTH  = "278.57144";
      _YLENGTH  = "200";
      break;
    case 1:					// This is for oldermeta files
      _OWNER    = string("bg");
      _MAGIC    = string("0x8002");
      _OBJECTS  = "28";
      _PICS     = "12";
      _ASPA     = "411";
      _ARWIDTH  = "270";
      _ARHEIGHT = "205";
      _XLENGTH  = "270";
      _YLENGTH  = "205";
      break;
    }
  }


void XWinMeta::SetAxis(int axis, int old)
  {
  _AXLOOP = "0";
  switch(old)
    {
    default:
    case 0:					// This is for newer meta files
      {
      _PENCOL   = string("0x1");
      _HUE      = string("0xdc");
      switch(axis)
        {
        default:
        case 1: 				//    Setting for 1st Axis
          _TYPNAM   = string("<yaxis>");
          _PICNUMB  = string("2");
          _XORIGIN  = "0";
          _YORIGIN  = "0.5";
          _XLEAST   = "0";
          _YLEAST   = "-3.4028235e+38";
          _XHIGHEST = "0.2";
          _YHIGHEST = "3.4028235e+38";
          _ANGLE1   = "0";
          _AXDLBL   = "1";
          _AXLABL   = "0";
          _AXFLAG   = "4";
          _AXFORM   = "0";
          _AXLDIR   = "1";
          _AXDTIC   = "2";
          _AXLIM1   = "9.5";
          _AXLIM2   = "-3";
          _AXCELL   = "0";
          _AXGDIS   = "0";
          _AXGLEN   = "0";
          _AXSTRG   = string("");
          break;
        case 2:					// Setting for 2nd Axis
          _TYPNAM   = string("<xaxis>");
          _PICNUMB  = string("3");
          _HUE      = string("0");
          _PENCOL   = string("0x1");
          _XORIGIN  = "0";
          _YORIGIN  = "0";
          _XLEAST   = "-3.4028235e+38";
          _YLEAST   = "0";
          _XHIGHEST = "3.4028235e+38";
          _YHIGHEST = "1";
          _ANGLE1   = "0";
          _AXDLBL   = "0";
          _AXLABL   = "1";
          _AXFLAG   = "1";
          _AXFORM   = "1";
          _AXLDIR   = "1";
          _AXSTRG   = string("<ppm>");
          _AXDTIC   = "0.5";
          _AXLIM1   = "0";
          _AXLIM2   = "0";
          _AXCELL   = "-4";
          _AXGDIS   = "0";
          _AXGLEN   = "0";
          break;
        }					// End new meta axis switch
      }						// End new meta file settings
      break;
    case 1:					// This is for older meta files
      {
      _PENCOL   = string("0x6");
      switch(axis)
        {
        default:
        case 1: 				//    Setting for 1st Axis
          _TYPNAM   = string("<f1axis>");
          _PICNUMB  = string("3");
          _XORIGIN  = "0";
          _YORIGIN  = "0";
          _XLEAST   = "-14";
          _YLEAST   = "-3.4028235e+38";
          _XHIGHEST = "2.5";
          _YHIGHEST = "3.4028235e+38";
          _ANGLE1   = "90";
          _AXDLBL   = "0";
          _AXLABL   = "1";
          _AXFLAG   = "4";
          _AXFORM   = "1";
          _AXLDIR   = "0";
          _AXDTIC   = "0.1";
          _AXLIM1   = "0";
          _AXLIM2   = "1";
          _AXCELL   = "-4";
          _AXGDIS   = "0.5";
          _AXGLEN   = "14";
          _AXSTRG   = string("<ppm>");
          _HUE      = string("0x82");
          break;
        case 2:					// Setting for 2nd Axis
          _TYPNAM   = string("<f2axis>");	// in older meta files
          _PICNUMB  = string("7");
          _HUE      = string("0x82");
          _XORIGIN  = "0";
          _YORIGIN  = "-1.5";
          _XLEAST   = "-3.4028235e+38";
          _YLEAST   = "0";
          _XHIGHEST = "3.4028235e+38";
          _YHIGHEST = "15.5";
          _ANGLE1   = "0";
          _AXDLBL   = "0";
          _AXLABL   = "1";
          _AXFLAG   = "1";
          _AXFORM   = "1";
          _AXLDIR   = "1";
          _AXSTRG   = string("<ppm>");
          _AXDTIC   = "0.1";
          _AXLIM1   = "0";
          _AXLIM2   = "1";
          _AXCELL   = "-4";
          _AXGDIS   = "0.5";
          _AXGLEN   = "14";
          break;
        }					// End old meta axis switch
      }						// End old meta axis settings
    }						// End switch new/old meta
  }

void XWinMeta::SetLine(int line, int old)
  {
  if(old) line *= -1;
  switch(line)
    {
    default:
    case -1:					// First line parameters old
    case 1:					// First line parameters new
      _TYPNAM   = string("<alabel>");
      _PICNUMB  = string("0");
      _DISABLE  = string("0x1");
      _SHADOW   = "0";
      _PLANENR  = string("0");
      _MODIFY   = string("0");
      _YORIGIN  = "0";
      _XSCALE   = "10";
      _YSCALE   = "10";
      _XLEAST   = "-3.4028235e+38";
      _YLEAST   = "-3.4028235e+38";
      _XHIGHEST = "3.4028235e+38";
      _YHIGHEST = "3.4028235e+38";
      _FILENAMES = string("");
      _LINEBF   = string("");
      break;
    case 2:
      _TYPNAM   = string("<slabel>");
      _PICNUMB  = string("5");
      _PLANENR  = string("0x7");
      _YORIGIN  = "0.5";
      _YSCALE   = "4.7044828e-08";
      _XLEAST   = "0";
      _YLEAST   = "0";
      _XHIGHEST = "0";
      _YHIGHEST = "0";
      _FILENAMES = string("<1r>");
      _LINEBF   = string("");
      _DISABLE  = string("0");
      _SHADOW   = "0";
      _MODIFY   = string("0x1");
      _XSCALE   = "1";
      break;
    case -2:
      _TYPNAM   = string("<frame>");
      _PICNUMB  = string("10");
      _PLANENR  = string("0x5");
      _YORIGIN  = "0";
      _YSCALE   = "1";
      _XLEAST   = "-3.4028235e+38";
      _YLEAST   = "-3.4028235e+38";
      _XHIGHEST = "3.4028235e+38";
      _YHIGHEST = "3.4028235e+38";
      _FILENAMES = string("<int2d>");
      _LINEBF   = string("(0..9)\n10\n10\n0\n10\n0\n0\n10\n0\n10\n10");
      _DISABLE  = string("0");
      _SHADOW   = "0";
      _MODIFY   = string("0x1");
      _XSCALE   = "1";
      break;
    }
  }

void XWinMeta::SetPic(int pic, int old)
  {
  if(old) pic *= -1;
  switch(pic)
    {
    default:
    case 1:					// First pic parameters new
      _TYPNAM   = string("<1spctrl>");
      _PICNUMB  = string("1");    
      _SOFTCLIP = string("0x1");
      _DISABLE  = string("0");    
      _PLANENR  = string("0x1");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _INTENSITY  = string("0xff");
      _XORIGIN  = "0";
      _YORIGIN  = "0";
      _XSCALE   = "10";
      _YSCALE   = "10";
      _XLEAST   = "0";
      _YLEAST   = "0";
      _XHIGHEST = "16";
      _YHIGHEST = "17";
      _ANGLE1   = "0";
      _PIC_ANZ  = "2";
      break;
    case -1:					// First pic parameters old
      _TYPNAM   = string("<2spctrl>");
      _PICNUMB  = string("1");    
      _SOFTCLIP = string("0x1");
      _DISABLE  = string("0");    
      _PLANENR  = string("0x1");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _INTENSITY  = string("0xff");
      _XORIGIN  = "0";
      _YORIGIN  = "0";
      _XSCALE   = "10";
      _YSCALE   = "10";
      _XLEAST   = "0";
      _YLEAST   = "0";
      _XHIGHEST = "21.5";
      _YHIGHEST = "20.5";
      _ANGLE1   = "0";
      _PIC_ANZ  = "2";
      break;
    case 2:                                     // Second pic parameters new
      _TYPNAM   = string("<ypic>");
      _PICNUMB  = string("2");
      _SOFTCLIP = string("0x1");
      _DISABLE  = string("0x1");
      _PLANENR  = string("0x1");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _INTENSITY  = string("0xff");
      _XORIGIN  = "0";
      _YORIGIN  = "3.5";
      _XSCALE   = "1";
      _YSCALE   = "1";
      _XLEAST   = "0";
      _YLEAST   = "-2.5";
      _XHIGHEST = "0.2";
      _YHIGHEST = "10";
      _ANGLE1   = "0";
      _PIC_ANZ  = "1";
      break;
    case -2:                                    // Second pic parameters old
      _TYPNAM   = string("<f1clip>");
      _PICNUMB  = string("2"); 
      _SOFTCLIP = string("0x1"); 
      _DISABLE  = string("0");
      _PLANENR  = string("0x1");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _INTENSITY  = string("0xff");
      _XORIGIN  = "19";
      _YORIGIN  = "1.5";
      _XSCALE   = "1"; 
      _YSCALE   = "1";
      _XLEAST   = "0";
      _YLEAST   = "-2.5";
      _XHIGHEST = "14";
      _YHIGHEST = "19";
      _ANGLE1   = "-90"; 
      _PIC_ANZ  = "1";
      break;
    case 3:					// Third pic parameters new
      _TYPNAM   = string("<region>");
      _PICNUMB  = string("3");    
      _SOFTCLIP = string("0");
      _DISABLE  = string("0");    
      _PLANENR  = string("0x1");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _INTENSITY  = string("0xff");
      _XORIGIN  = "348.3031";
      _YORIGIN  = "0";
      _XSCALE   = "-7.8445568";
      _YSCALE   = "1";
      _XLEAST   = "-3.4028235e+38";
      _YLEAST   = "-3.4028235e+38";
      _XHIGHEST = "3.4028235e+38";
      _YHIGHEST = "3.4028235e+38";
      _ANGLE1   = "0";
      _PIC_ANZ  = "2";
      break;
    case -3:					// Third pic parameters old
      _TYPNAM   = string("<f1reg>");
      _PICNUMB  = string("3");    
      _SOFTCLIP = string("0");
      _DISABLE  = string("0");    
      _PLANENR  = string("0x1");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _INTENSITY  = string("0xff");
      _XORIGIN  = "13.195404";
      _YORIGIN  = "0";
      _XSCALE   = "-1.4836284";
      _YSCALE   = "1";
      _XLEAST   = "-3.4028235e+38";
      _YLEAST   = "-3.4028235e+38";
      _XHIGHEST = "3.4028235e+38";
      _YHIGHEST = "3.4028235e+38";
      _ANGLE1   = "0";
      _PIC_ANZ  = "2";
      break;
    case 4:                                     // Fourth pic parameters new
      _TYPNAM   = string("<filoff>");
      _PICNUMB  = string("4");
      _SOFTCLIP = string("0");
      _DISABLE  = string("0");
      _PLANENR  = string("0x1");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _INTENSITY  = string("0xff");
      _XORIGIN  = "0";
      _YORIGIN  = "3.5";
      _XSCALE   = "-0.0029305699";
      _YSCALE   = "1";
      _XLEAST   = "0";
      _YLEAST   = "-3.4028235e+38";
      _XHIGHEST = "2.1474836e+09";
      _YHIGHEST = "3.4028235e+38";
      _ANGLE1   = "0";
      _PIC_ANZ  = "3";
      break;
    case -4:					// Third pic parameters old
      _TYPNAM   = string("<f1filo>");
      _PICNUMB  = string("4");    
      _SOFTCLIP = string("0");
      _DISABLE  = string("0");    
      _PLANENR  = string("0x1");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _INTENSITY  = string("0xff");
      _XORIGIN  = "0";
      _YORIGIN  = "14.5";
      _XSCALE   = "-0.0029305699";
      _YSCALE   = "1e-07";
      _XLEAST   = "0";
      _YLEAST   = "-3.4028235e+38";
      _XHIGHEST = "2.1474836e+09";
      _YHIGHEST = "3.4028235e+38";
      _ANGLE1   = "0";
      _PIC_ANZ  = "2";
      break;
    case 5:                                    // Fifth pic parameters old
      _TYPNAM   = string("<data>");
      _PICNUMB  = string("5");
      _SOFTCLIP = string("0x1");
      _DISABLE  = string("0");
      _PLANENR  = string("0x1");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _INTENSITY  = string("0xff");
      _XORIGIN  = "0";
      _YORIGIN  = "0";
      _XSCALE   = "1";
      _YSCALE   = "1";
      _XLEAST   = "0";
      _YLEAST   = "0";
      _XHIGHEST = "2.1474836e+09";
      _YHIGHEST = "10";
      _ANGLE1   = "0";
      _PIC_ANZ  = "3";
      break;
    case -5:					// Third pic parameters old
      _TYPNAM   = string("<f1pclip>");
      _PICNUMB  = string("5");    
      _SOFTCLIP = string("0x1");
      _DISABLE  = string("0");    
      _PLANENR  = string("0x1");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _INTENSITY  = string("0xff");
      _XORIGIN  = "0";
      _YORIGIN  = "0";
      _XSCALE   = "1";
      _YSCALE   = "1";
      _XLEAST   = "0";
      _YLEAST   = "-5000000";
      _XHIGHEST = "2.1474836e+09";
      _YHIGHEST = "45000000";
      _ANGLE1   = "0";
      _PIC_ANZ  = "2";
      break;
    case 6:                                     // Sixth Pic parameters new
      _TYPNAM   = string("<tipic>");
      _PICNUMB  = string("6");
      _SOFTCLIP = string("0x1");
      _DISABLE  = string("0");
      _PLANENR  = string("0x1");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _INTENSITY  = string("0xff");
      _XORIGIN  = "0";
      _YORIGIN  = "180";
      _XSCALE   = "10";
      _YSCALE   = "10";
      _XLEAST   = "0";
      _YLEAST   = "0";
      _XHIGHEST = "25";
      _YHIGHEST = "2";
      _ANGLE1   = "0";
      _PIC_ANZ  = "1";
      break;
    case -6:					// Third pic parameters old
      _TYPNAM   = string("<f2clip>");
      _PICNUMB  = string("6");    
      _SOFTCLIP = string("0x1");
      _DISABLE  = string("0");    
      _PLANENR  = string("0x1");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _INTENSITY  = string("0xff");
      _XORIGIN  = "5";
      _YORIGIN  = "1.5";
      _XSCALE   = "1";
      _YSCALE   = "1";
      _XLEAST   = "0";
      _YLEAST   = "-1.5";
      _XHIGHEST = "14";
      _YHIGHEST = "19";
      _ANGLE1   = "0";
      _PIC_ANZ  = "1";
      break;
    case 7:                                     // Seventh Pic parameters new
      _TYPNAM   = string("<parpic>");
      _PICNUMB  = string("7");
      _SOFTCLIP = string("0x1");
      _DISABLE  = string("0");
      _PLANENR  = string("0x4");
      _SATURATION = string("0");
      _HUE        = string("0");
      _INTENSITY  = string("0");
      _XORIGIN  = "210";
      _YORIGIN  = "180";
      _XSCALE   = "10";
      _YSCALE   = "10";
      _XLEAST   = "0";
      _YLEAST   = "-18";
      _XHIGHEST = "6.8571429";
      _YHIGHEST = "0";
      _ANGLE1   = "0";
      _PIC_ANZ  = "1";
      break;
    case -7:                                     // Seventh Pic parameters new
      _TYPNAM   = string("<f2reg>");
      _PICNUMB  = string("7");
      _SOFTCLIP = string("0");
      _DISABLE  = string("0");
      _PLANENR  = string("0x1");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _INTENSITY  = string("0xff");
      _XORIGIN  = "13.195434";
      _YORIGIN  = "0";
      _XSCALE   = "-1.4836284";
      _YSCALE   = "1";
      _XLEAST   = "-3.4028235e+38";
      _YLEAST   = "-3.4028235e+38";
      _XHIGHEST = "3.4028235e+38";
      _YHIGHEST = "3.4028235e+38";
      _ANGLE1   = "0";
      _PIC_ANZ  = "3";
      break;
    case -8:					// Third pic parameters old
      _TYPNAM   = string("<f2filo>");
      _PICNUMB  = string("8");    
      _SOFTCLIP = string("0");
      _DISABLE  = string("0");    
      _PLANENR  = string("0x1");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _INTENSITY  = string("0xff");
      _XORIGIN  = "0";
      _YORIGIN  = "14.5";
      _XSCALE   = "-0.0029305699";
      _YSCALE   = "1e-07";
      _XLEAST   = "0";
      _YLEAST   = "-3.4028235e+38";
      _XHIGHEST = "2.1474836e+09";
      _YHIGHEST = "3.4028235e+38";
      _ANGLE1   = "0";
      _PIC_ANZ  = "2";
      break;
    case -9:					// Third pic parameters old
      _TYPNAM   = string("<f2pclip>");
      _PICNUMB  = string("9");    
      _SOFTCLIP = string("0x1");
      _DISABLE  = string("0");    
      _PLANENR  = string("0x1");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _INTENSITY  = string("0xff");
      _XORIGIN  = "0";
      _YORIGIN  = "0";
      _XSCALE   = "1";
      _YSCALE   = "1";
      _XLEAST   = "0";
      _YLEAST   = "-5000000";
      _XHIGHEST = "2.1474836e+09";
      _YHIGHEST = "45000000";
      _ANGLE1   = "0";
      _PIC_ANZ  = "2";
      break;
    case -10:					// Third pic parameters old
      _TYPNAM   = string("<cntclip>");
      _PICNUMB  = string("10");    
      _SOFTCLIP = string("0x1");
      _DISABLE  = string("0");    
      _PLANENR  = string("0x1");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _INTENSITY  = string("0xff");
      _XORIGIN  = "0";
      _YORIGIN  = "13.195404";
      _XSCALE   = "1";
      _YSCALE   = "-1.4836284";
      _XLEAST   = "-0.54229641";
      _YLEAST   = "-0.54231644";
      _XHIGHEST = "8.8940287";
      _YHIGHEST = "8.8940086";
      _ANGLE1   = "0";
      _PIC_ANZ  = "3";
      break;
    case -11:					// Third pic parameters old
      _TYPNAM   = string("<tipic>");
      _PICNUMB  = string("11");    
      _SOFTCLIP = string("0x1");
      _DISABLE  = string("0");    
      _PLANENR  = string("0x1");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _INTENSITY  = string("0xff");
      _XORIGIN  = "20";
      _YORIGIN  = "180";
      _XSCALE   = "10";
      _YSCALE   = "10";
      _XLEAST   = "0";
      _YLEAST   = "0";
      _XHIGHEST = "25";
      _YHIGHEST = "2";
      _ANGLE1   = "0";
      _PIC_ANZ  = "1";
      break;
    case -12:					// Third pic parameters old
      _TYPNAM   = string("<parpic>");
      _PICNUMB  = string("12");    
      _SOFTCLIP = string("0x1");
      _DISABLE  = string("0");    
      _PLANENR  = string("0x4");
      _SATURATION = string("0");
      _HUE        = string("0");
      _INTENSITY  = string("0");
      _XORIGIN  = "200";
      _YORIGIN  = "180";
      _XSCALE   = "10";
      _YSCALE   = "10";
      _XLEAST   = "0";
      _YLEAST   = "-18";
      _XHIGHEST = "6.8571429";
      _YHIGHEST = "0";
      _ANGLE1   = "0";
      _PIC_ANZ  = "1";
      break;
    }
  }

void XWinMeta::SetSpec(int spec, int old)
  {
  if(old) spec *= -1;
  switch(spec)
    {
    default:
    case 1:					// First spec parameters new
      _TYPNAM     = string("<spctrm>");
      _PICNUMB    = string("5");
      _DISABLE    = string("0");
      _PLANENR    = string("0x7");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _PENCOL     = string("0");
      _XORIGIN    = "0";
      _YORIGIN    = "0.5";
      _YSCALE     = "4.7044828e-08";
      _YLEAST     = "0";
      _YHIGHEST   = "0";
      _FILENAMES   = string("<1r>");
      _D_TYPE     = string("0"); 
      _SCALFLG    = "10";
      break;
    case -1:					// First spec parameters new
      _TYPNAM     = string("<f1projp>");
      _PICNUMB    = string("5");
      _DISABLE    = string("0");
      _PLANENR    = string("0x4");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _PENCOL     = string("0");
      _XORIGIN    = "0";
      _YORIGIN    = "0";
      _YSCALE     = "1";
      _YLEAST     = "0";
      _YHIGHEST   = "0";
      _FILENAMES   = string("<p2r2>");
      _D_TYPE     = string("0"); 
      _SCALFLG    = "4.5";
      break;
    case 2:					// Second spec parameters new
      _TYPNAM     = string("<intgrl>");
      _PICNUMB    = string("5");
      _DISABLE    = string("0");
      _PLANENR    = string("0x6");
      _SATURATION = string("0"); 
      _HUE        = string("0x55");
      _PENCOL     = string("0x4");
      _XORIGIN    = "0";
      _YORIGIN    = "1";
      _YSCALE     = "5.4840832e-10";
      _YLEAST     = "-3.4028235e+38";
      _YHIGHEST   = "3.4028235e+38";
      _FILENAMES   = string("<int>");
      _D_TYPE     = string("0x1");
      _SCALFLG    = "12";
      break;
    case -2:					// Second spec parameters old
      _TYPNAM     = string("<f1projn>");
      _PICNUMB    = string("5");
      _DISABLE    = string("0x1");
      _PLANENR    = string("0x4");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _PENCOL     = string("0");
      _XORIGIN    = "0";
      _YORIGIN    = "0";
      _YSCALE     = "1";
      _YLEAST     = "0";
      _YHIGHEST   = "0";
      _FILENAMES   = string("<n2r2>");
      _D_TYPE     = string("0"); 
      _SCALFLG    = "4.5";
      break;
    case -3:					// Third spec parameters old 
      _TYPNAM     = string("<f2projp>");
      _PICNUMB    = string("9");
      _DISABLE    = string("0");
      _PLANENR    = string("0x4");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _PENCOL     = string("0");
      _XORIGIN    = "0";
      _YORIGIN    = "0";
      _YSCALE     = "1";
      _YLEAST     = "0";
      _YHIGHEST   = "0";
      _FILENAMES   = string("<p2r1>");
      _D_TYPE     = string("0"); 
      _SCALFLG    = "4.5";
      break;
    case -4:					// Fourth spec parameters old 
      _TYPNAM     = string("<f2projn>");
      _PICNUMB    = string("9");
      _DISABLE    = string("0x1");
      _PLANENR    = string("0x4");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _PENCOL     = string("0");
      _XORIGIN    = "0";
      _YORIGIN    = "0";
      _YSCALE     = "1";
      _YLEAST     = "0";
      _YHIGHEST   = "0";
      _FILENAMES   = string("<n2r1>");
      _D_TYPE     = string("0"); 
      _SCALFLG    = "4.5";
      break;
    }
  }


void XWinMeta::SetText(int txt, int old)
  {
  if(old) txt *= -1;
  switch(txt)
    {
    default:
    case 1:					// First spec parameters new
      _TYPNAM   = string("<title>");
      _PICNUMB  = string("6");
      _PLANENR  = string("0x4");
      _MODIFY   = string("0x1");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _INTENSITY  = string("0xff");
      _PENCOL     = string("0"); 
      _YLEAST   = "0";
      _XHIGHEST = "25";
      _YHIGHEST = "2";
      _FILENAMES = string("<title>");
      _TXMODE   = "97";
      _TXALLI   = "5";
      _TXFNST   = "0";
      _TXSOSP   = "33";
      _TXCELL   = "-6";
      break;
    case -1:					// First spec parameters old
      _TYPNAM   = string("<title>");
      _PICNUMB  = string("11");
      _PLANENR  = string("0x4");
      _MODIFY   = string("0x1");
      _SATURATION = string("0");
      _HUE        = string("0x73");
      _INTENSITY  = string("0xc0");
      _PENCOL     = string("0x9"); 
      _YLEAST   = "0";
      _XHIGHEST = "25";
      _YHIGHEST = "2";
      _FILENAMES = string("<title>");
      _TXMODE   = "97";
      _TXALLI   = "1";
      _TXFNST   = "5";
      _TXSOSP   = "33";
      _TXCELL   = "-4";
      break;
    case 2:					// Second spec parameters new
      _TYPNAM   = string("<parmtr>");
      _PICNUMB  = string("7");
      _PLANENR  = string("0");
      _MODIFY   = string("0");
      _SATURATION = string("0xff");
      _HUE        = string("0x80");
      _INTENSITY  = string("0xff");
      _PENCOL     = string("0"); 
      _YLEAST   = "-18";
      _XHIGHEST = "6.8571429";
      _YHIGHEST = "0";
      _FILENAMES = string("<parm.txt>");
      _TXMODE   = "0";
      _TXALLI   = "3";
      _TXFNST   = "0";
      _TXSOSP   = "0";
      _TXCELL   = "-4";
      break;
    case -2:					// Second spec parameters old
      _TYPNAM   = string("<parmtr>");
      _PICNUMB  = string("12");
      _PLANENR  = string("0");
      _MODIFY   = string("0");
      _SATURATION = string("0");
      _HUE        = string("0x73");
      _INTENSITY  = string("0xc0");
      _PENCOL     = string("0x9"); 
      _YLEAST   = "-18";
      _XHIGHEST = "6.8571429";
      _YHIGHEST = "0";
      _FILENAMES = string("<parm.txt>");
      _TXMODE   = "0";
      _TXALLI   = "5";
      _TXFNST   = "0";
      _TXSOSP   = "0";
      _TXCELL   = "-4";
      break;
    }
  }

void XWinMeta::SetPeak(int pk, int old)
  {
  _TYPK      = string("<PEAKLBL>");
  _TYPNR     = "12";
  _SOFTCLIP  = string("0");
  _SHADOW    = "0";
  _MODIFY    = string("0x1");
  _INTENSITY = string("0xff");
  _XORIGIN  = "0";
  _ZORIGIN  = "0";
  _YSCALE   = "1";
  _ZSCALE   = "1";
  _ZLEAST   = "0";
  _ZHIGHEST = "0";
  _ANGLE1   = "0";
  _ANGLE2   = "0";
  _ANGLE3   = "0";
  if(old) pk *= -1;
  switch(pk)
    {
    default:
    case 1:                                     // First peak parameters new
      _TYPNAM     = string("<peaklb>");
      _PICNUMB    = string("4");
      _PLANENR    = string("0x7");
      _DISABLE    = string("0");
      _SATURATION = string("0");
      _HUE        = string("0xdc");
      _PENCOL     = string("0x1"); 
      _YORIGIN    = "0.5";
      _YSCALE     = "1e-07";
      _XLEAST     = "-3.4028235e+38";
      _YLEAST     = "-5000000";
      _XHIGHEST   = "3.4028235e+38";
      _YHIGHEST   = "1.16e+08";
      _FILENAMES   = string("<peaks>");
      _PLROTA     = "0";
      _PLDIWI     = "7";
      _PLLABF     = string("");
      _PLSTRG     = string("<Hz>");
      _PLDHEI     = "1e+08";
      _PLMSHI     = "130.49713";
      break;
    case -1:                                     // First peak parameters old 
      _TYPNAM     = string("<f1peak>");
      _PICNUMB    = string("4");
      _DISABLE    = string("0x1");
      _PLANENR    = string("0x7");
      _SATURATION = string("0");
      _HUE        = string("0x73");
      _PENCOL     = string("0x5"); 
      _YORIGIN    = "0";
      _XLEAST     = "-3.4028235e+38";
      _YLEAST     = "-5000000";
      _YSCALE     = "1";
      _XHIGHEST   = "3.4028235e+38";
      _YHIGHEST   = "70000000";
      _FILENAMES   = string("<peak2p>");
      _PLROTA     = "1";
      _PLDIWI     = "4";
      _PLLABF     = string("<label2p>");
      _PLSTRG     = string("<ppm>");
      _PLDHEI     = "50000000";
      _PLMSHI     = "689.992";
      break;
    case -2:                                     // Second peak parameters old 
      _TYPNAM     = string("<f2peak>");
      _PICNUMB    = string("8");
      _PLANENR    = string("0x7");
      _DISABLE    = string("0x1");
      _HUE        = string("0x73");
      _SATURATION = string("0");
      _PENCOL     = string("0x5"); 
      _YORIGIN    = "0";
      _XLEAST     = "-3.4028235e+38";
      _YLEAST     = "-5000000";
      _YSCALE     = "1";
      _XHIGHEST   = "3.4028235e+38";
      _YHIGHEST   = "70000000";
      _FILENAMES   = string("<peak1p>");
      _PLROTA     = "0";
      _PLDIWI     = "4";
      _PLLABF     = string("<label1p>");
      _PLSTRG     = string("<ppm>");
      _PLDHEI     = "50000000";
      _PLMSHI     = "689.99188";
      break;
    }
  }


void XWinMeta::SetImag(int im, int old)
  {
  _TYPK       = string("<IMAG>");
  _TYPNR      = "8";
  _PICNUMB    = string("10");
  _PLANENR    = string("0xff");
  _PENCOL     = string("0");
  _SOFTCLIP   = string("0");
  _SHADOW     = "0";
  _MODIFY     = string("0x1");
  _INTENSITY  = string("0xff");
  _XORIGIN    = "0";
  _ZORIGIN    = "0";
  _XSCALE     = "1";
  _ZSCALE     = "1";
  _ZLEAST     = "0";
  _ZHIGHEST   = "0";
  _XLEAST     = "0";
  _YLEAST     = "0";
  _XHIGHEST   = "0";
  _YHIGHEST   = "0";
  _ANGLE1     = "0";
  _ANGLE2     = "0";
  _ANGLE3     = "0";
  if(old) im *= -1;
  switch(im)
    {
    default:
    case 1:                                      // First imag parameters new
    case -1:                                     // First imag parameters old
      _TYPNAM     = string("<contour>");
      _DISABLE    = string("0");
      _PENCOL     = string("0");
      _FILENAMES   = string("<2rr>");
      _CONTRI     = "0";
      break;
    case -2:                                     // Second peak parameters ol
      _TYPNAM     = string("<contadd>");
      _DISABLE    = string("0x1");
      _FILENAMES   = string("</u/data/demo/nmr/H5/1/pdata/1/2rr>");
      _CONTRI     = "-1";
      break;
    }
  }


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A             XWinMeta Parameter File Constructors, Destructor
// ____________________________________________________________________________
 
/* These are the constructors of the class handling Bruker XWinNMR output
   parameter files.  This doesn't do anything in particular, it is the write
   functions that perform the work.                                          */

XWinMeta::XWinMeta()
  {
  mname = string("meta");		// Set a default file name
  oldflag = 0;				// Set for newer format
  }

XWinMeta::XWinMeta(const string& name)
  {
  mname = name;				// Set a default file name
  oldflag = 0;				// Set for newer format
  }

XWinMeta::~XWinMeta()  { }                                                         


// ____________________________________________________________________________
// B                  XWinMeta Parameter Access Functions
// ____________________________________________________________________________

/* These functions allow users to set any more important parameters that
   these files contain.                                                      */

void XWinMeta::OldFlag(int of)
  { of?oldflag=1:oldflag=0; }

// ____________________________________________________________________________
// C                         XWinMeta Output Functions
// ____________________________________________________________________________
 
/* These function allow for output of NMR parameters directly into a Bruker
   XWinNMR ASCII parameter file (meta).                                      */

int XWinMeta::write(const string& dname, int warn, int of)
  { mname = dname; if(of>=0) OldFlag(of); return write(warn); }

int XWinMeta::write(int warn)
  {
  ofstream ofstr(mname.c_str());		// Open ASCII file for output
  if(!ofstr.good())				// If file bad then exit
     {
     if(warn)
       {
       XWinMetaerror(2, 1);			// Output filestream problems
       XWinMetaerror(25,1);			// Can't write parameters out
       if(warn==1) XWinMetaerror(21,mname,1);	// Can't write parameters out
       else        XWinMetafatality(21,mname);	// with out without die
       } 
     return 0;
     }  
  string nn("##");			// Bruker ## parameter line start
  string nns("##$");			// Bruker ##$ parameter line start
  string ss("$$");			// Bruker $$ parameter line start
  SetBase(oldflag);			// Set base meta parameters
  writeBase(ofstr);			// Output base parameters
  writeInit(ofstr);			// Output initial block
  if(oldflag) writeLutab(ofstr);	// Output LUTAB block (old)
  SetLine(1,oldflag);			// Set 1st line parameters
  writeLine(ofstr);			// Write 1st line
  SetPic(1, oldflag);			// Set 1st pic parameters
  writePic(ofstr);			// Write 1st pic
  SetPic(2, oldflag);			// Set 2nd pic parameters
  writePic(ofstr);			// Write 2nd pic
  if(oldflag)				// Older meta different
    { 					// order now than newer meta
    SetPic(3, oldflag);			// Set 3rd pic parameters
    writePic(ofstr);			// Write 3rd pic (old only)
    SetAxis(1, oldflag);		// Set up axis parameters
    writeAxis(ofstr);			// Output axis parameters
    SetPic(4, oldflag);			// Set 4th pic parameters
    writePic(ofstr);			// Write 4th pic
    SetPic(5, oldflag);			// Set 5th pic parameters
    writePic(ofstr);			// Write 5th pic
    SetSpec(1, oldflag);		// Set 1st spectrum parameters
    writeSpec(ofstr);			// Write 1st spectrum parameters
    SetSpec(2, oldflag);		// Set 2nd spectrum parameters
    writeSpec(ofstr);			// Write 2nd spectrum parameters
    SetPeak(1, oldflag);		// Set 1st peak parameters
    writePeak(ofstr);			// Write 1st peak parameters
    SetPic(6, oldflag);			// Set 6th pic parameters
    writePic(ofstr);			// Write 6th pic
    SetPic(7, oldflag);			// Set 7th pic parameters
    writePic(ofstr);			// Write 7th pic
    SetAxis(2, oldflag);		// Set up axis parameters
    writeAxis(ofstr);			// Output axis parameters
    SetPic(8, oldflag);			// Set 8th pic parameters
    writePic(ofstr);			// Write 8th pic
    SetPic(9, oldflag);			// Set 9th pic parameters
    writePic(ofstr);			// Write 9th pic
    SetSpec(3, oldflag);		// Set 3rd spectrum parameters
    writeSpec(ofstr);			// Write 1st spectrum parameters
    SetSpec(4, oldflag);		// Set 4th spectrum parameters
    writeSpec(ofstr);			// Write 2nd spectrum parameters
    SetPeak(2, oldflag);		// Set 2nd peak parameters
    writePeak(ofstr);			// Write 2nd peak parameters
    SetPic(10, oldflag);		// Set 10th pic parameters
    writePic(ofstr);			// Write 10th pic
    SetLine(2,oldflag);			// Set 2nd line parameters
    writeLine(ofstr);			// Write 2nd line
    SetImag(1,oldflag);			// Set 1st imag block
    writeImag(ofstr);			// Write 1st image 
    SetImag(2,oldflag);			// Set 2nd imag block
    writeImag(ofstr);			// Write 2nd image
    SetPic(11, oldflag);		// Set 11th pic parameters
    writePic(ofstr);			// Write 11th pic
    SetText(1, oldflag);		// Set 1st text parameters
    writeText(ofstr);			// Write 1st text block
    SetPic(12, oldflag);		// Set 12th pic parameters
    writePic(ofstr);			// Write 12th pic
    SetText(2, oldflag);		// Set 2nd text parameters
    writeText(ofstr, 0);		// Write 2nd text block
    }
  else
    {
    SetAxis(1, oldflag);		// Set up axis parameters
    writeAxis(ofstr);			// Output axis parameters
    SetPic(3, oldflag);			// Set 3rd pic parameters
    writePic(ofstr);			// Write 3rd pic
    SetAxis(2, oldflag);		// Set up axis parameters
    writeAxis(ofstr);			// Output axis parameters
    SetPic(4, oldflag);			// Set 1st pic parameters
    writePic(ofstr);			// Write 1st pic
    SetPic(5, oldflag);			// Set 1st pic parameters
    writePic(ofstr);			// Write 1st pic
    SetSpec(1, oldflag);		// Set 1st spectrum parameters
    writeSpec(ofstr);			// Write 1st spectrum parameters
    SetLine(2, oldflag);		// Set 2nd line parameters
    writeLine(ofstr);			// Write 2nd line
    SetSpec(2, oldflag);		// Set 2nd spectrum parameters
    writeSpec(ofstr);			// Write 2nd spectrum parameters
    ofstr << nns << "TYPK= <INTLBL>\n";
    ofstr << nns << "TYPNR= 13\n";
    ofstr << nns << "TYPNAM= <intlbl>\n";
    ofstr << nns << "PICNUMB= 4\n";
    ofstr << nns << "SOFTCLIP= 0\n";
    ofstr << nns << "DISABLE= 0\n";
    ofstr << nns << "SHADOW= 0\n";
    ofstr << nns << "PLANENR= 0x6\n";
    ofstr << nns << "MODIFY= 0x1\n";
    ofstr << nns << "SATURATION= 0xff\n";
    ofstr << nns << "HUE= 0x80\n";
    ofstr << nns << "INTENSITY= 0xff\n";
    ofstr << nns << "PENCOL= 0\n";
    ofstr << nns << "XORIGIN= 0\n";
    ofstr << nns << "YORIGIN= 0\n";
    ofstr << nns << "ZORIGIN= 0\n";
    ofstr << nns << "XSCALE= 1\n";
    ofstr << nns << "YSCALE= 1\n";
    ofstr << nns << "ZSCALE= 1\n";
    ofstr << nns << "XLEAST= -3.4028235e+38\n";
    ofstr << nns << "YLEAST= -2.5\n";
    ofstr << nns << "ZLEAST= 0\n";
    ofstr << nns << "XHIGHEST= 3.4028235e+38\n";
    ofstr << nns << "YHIGHEST= 0\n";
    ofstr << nns << "ZHIGHEST= 0\n";
    ofstr << nns << "ANGLE1= 0\n";
    ofstr << nns << "ANGLE2= 0\n";
    ofstr << nns << "ANGLE3= 0\n";
    ofstr << nns << "FILENAME= <int>\n";
    ofstr << nns << "ILDIWI= 6\n";
    ofstr << nns << "ILROTA= 0\n";
    ofstr << nns << "ILCELL= -4\n";
    ofstr << nns << "ILCONT= 0\n";
    ofstr << nns << "ILSCAL= 0\n";
    ofstr << nns << "ILINVL= 0.2\n";
    ofstr << nns << "ILINSH= 0.5\n";
    ofstr << "\n";
    SetPeak(1, oldflag);		// Set 1st peak parameters
    writePeak(ofstr);		// Write 1st peak parameters
    SetPic(6, oldflag);			// Set 6th pic parameters
    writePic(ofstr);			// Write 6th pic
    SetText(1, oldflag);			// Set 1st text parameters
    writeText(ofstr);			// Write 6th pic
    SetPic(7, oldflag);			// Set 1st pic parameters
    writePic(ofstr);			// Write 1st pic
    SetText(2, oldflag);			// Set 2nd text parameters
    writeText(ofstr, 0);			// Write 6th pic
    }
  ofstr << nn << "END= ";
  return 1;
  }


void XWinMeta::writeBase(ofstream& ofstr) const
  {
  string nn("##");
  string nns("##$");
  ofstr << nn << "TITLE= UXNMR Metafile\n";
  ofstr << nn << "DATATYPE= UXNMR Metafile Objects\n";
  ofstr << nn << "ORIGIN= UXNMR, Bruker Analytische Messtechnik GmbH\n";
  ofstr << nn  << "OWNER= "             << _OWNER       << "\n";
  ofstr << nns << "MAGIC NUMBER= "      << _MAGIC       << "\n";
  ofstr << nns << "VERSION= "           << "<930315.0>" << "\n";
  ofstr << nns << "NUMBER OF OBJECTS= " << _OBJECTS     << "\n";
  ofstr << nns << "NUMBER OF PICS= "    << _PICS        << "\n";
  ofstr << nns << "ASPA SIZE= "         << _ASPA        << "\n";
  }

void XWinMeta::writeInit(ofstream& ofstr, int ef) const
  {
  string nns("##$");
  ofstr << nns << "TYPK= "         << "<INIT>"   << "\n";
  ofstr << nns << "TYPNR= "        << "0"        << "\n";
  ofstr << nns << "TYPNAM= "       << "<survey>" << "\n";
  ofstr << nns << "ORIENT= "       << "0"        << "\n";
  ofstr << nns << "PAPER_FLG= "    << "0"        << "\n";
  ofstr << nns << "INIFIT= "       << "0x1"      << "\n";
  ofstr << nns << "EDGE_COL= "     << "0"        << "\n";
  ofstr << nns << "PAPER_FORMAT= " << "0"        << "\n";
  ofstr << nns << "ARWIDTH= "   << _ARWIDTH  << "\n";
  ofstr << nns << "ARHEIGHT= "  << _ARHEIGHT << "\n";
  ofstr << nns << "XHLOLEFT= "  << "0"       << "\n";
  ofstr << nns << "YHLOLEFT= "  << "0"       << "\n";
  ofstr << nns << "XHUPRIGHT= " << "1"       << "\n";
  ofstr << nns << "YHUPRIGHT= " << "1"       << "\n";
  ofstr << nns << "XLENGTH= "   << _XLENGTH  << "\n";
  ofstr << nns << "YLENGTH= "   << _YLENGTH  << "\n";
  ofstr << nns << "REDUCTION= " << "1"       << "\n";
  ofstr << nns << "NOARPRES= "  << "0"       << "\n";
  if(ef) ofstr << "\n";
  }


void XWinMeta::writeLutab(ofstream& ofstr, int ef) const
  {
  string nns("##$");
  ofstr << nns << "TYPK= <LUTAB>\n";
  ofstr << nns << "TYPNR= 14\n";
  ofstr << nns << "TYPNAM= <lutab>\n";
  ofstr << nns << "PICNUMB= 0\n";
  ofstr << nns << "SOFTCLIP= 0\n";
  ofstr << nns << "DISABLE= 0\n";
  ofstr << nns << "SHADOW= 0\n";
  ofstr << nns << "PLANENR= 0x5\n";
  ofstr << nns << "MODIFY= 0x1\n";
  ofstr << nns << "SATURATION= 0\n";
  ofstr << nns << "HUE= 0x82\n";
  ofstr << nns << "INTENSITY= 0xff\n";
  ofstr << nns << "PENCOL= 0x6\n";
  ofstr << nns << "XORIGIN= 627\n";
  ofstr << nns << "YORIGIN= 4\n";
  ofstr << nns << "ZORIGIN= 0\n";
  ofstr << nns << "XSCALE= 1\n";
  ofstr << nns << "YSCALE= 1\n";
  ofstr << nns << "ZSCALE= 1\n";
  ofstr << nns << "XLEAST= 0\n";
  ofstr << nns << "YLEAST= 0\n";
  ofstr << nns << "ZLEAST= 0\n";
  ofstr << nns << "XHIGHEST= 60\n";
  ofstr << nns << "YHIGHEST= 250\n";
  ofstr << nns << "ZHIGHEST= 0\n";
  ofstr << nns << "ANGLE1= 0\n";
  ofstr << nns << "ANGLE2= 0\n";
  ofstr << nns << "ANGLE3= 0\n";
  ofstr << nns << "LUT_THRES= 190.60001\n";
  ofstr << nns << "LUT_WINDOW= 6.999999\n";
  ofstr << nns << "LUT_GREY= 2\n";
  ofstr << nns << "LUT_MODE= 1\n";
  if(ef) ofstr << "\n";
  }


void XWinMeta::writeLine(ofstream& ofstr, int ef) const
  {
  string nns("##$");
  ofstr << nns << "TYPK= "       << "<LINE>" << "\n";
  ofstr << nns << "TYPNR= "      << "6"      << "\n";
  ofstr << nns << "TYPNAM= "     << _TYPNAM  << "\n";
  ofstr << nns << "PICNUMB= "    << _PICNUMB << "\n";
  ofstr << nns << "SOFTCLIP= "   << "0"      << "\n";
  ofstr << nns << "DISABLE= "    << _DISABLE << "\n";
  ofstr << nns << "SHADOW= "     << _SHADOW  << "\n";
  ofstr << nns << "PLANENR= "    << _PLANENR << "\n";
  ofstr << nns << "MODIFY= "     << _MODIFY  << "\n";
  ofstr << nns << "SATURATION= " << "0xff"   << "\n";
  ofstr << nns << "HUE= "        << "0x80"   << "\n";
  ofstr << nns << "INTENSITY= "  << "0xff"   << "\n";
  ofstr << nns << "PENCOL= "     << "0"      << "\n";
  ofstr << nns << "XORIGIN= "    << "0"      << "\n";
  ofstr << nns << "YORIGIN= "    << _YORIGIN << "\n";
  ofstr << nns << "ZORIGIN= "    << "0"      << "\n";
  ofstr << nns << "XSCALE= "     << _XSCALE  << "\n";
  ofstr << nns << "YSCALE= "     << _YSCALE  << "\n";
  ofstr << nns << "ZSCALE= "     << "1"      << "\n";
  ofstr << nns << "XLEAST= "     << _XLEAST  << "\n";
  ofstr << nns << "YLEAST= "     << _YLEAST  << "\n";
  ofstr << nns << "ZLEAST= "     << "0"      << "\n";
  ofstr << nns << "XHIGHEST= "   << _XHIGHEST<< "\n";
  ofstr << nns << "YHIGHEST= "   << _YHIGHEST<< "\n";
  ofstr << nns << "ZHIGHEST= "   << "0"      << "\n";
  ofstr << nns << "ANGLE1= "     << "0"      << "\n";
  ofstr << nns << "ANGLE2= "     << "0"      << "\n";
  ofstr << nns << "ANGLE3= "     << "0"      << "\n";
  if(_FILENAMES.length())
    ofstr << nns << "FILENAME= " << _FILENAMES<< "\n";
  ofstr << nns << "LIDIMS= "     << "2"      << "\n";
  ofstr << nns << "LIMARK= "     << "0"      << "\n";
  ofstr << nns << "LIFILE= "     << "0"      << "\n";
  ofstr << nns << "LIFONT= "     << "0"      << "\n";
  if(_LINEBF.length())
    ofstr << nns << "LINEBF= "   << _LINEBF  << "\n";
  ofstr << nns << "LIMKSI= "     << "0"      << "\n";
  ofstr << nns << "LIZCOO= "     << "0"      << "\n";
  ofstr << nns << "LICELL= "     << "0"      << "\n";
  ofstr << nns << "LISLNT= "     << "0"      << "\n";
  if(ef) ofstr << "\n";
  }


void XWinMeta::writeAxis(ofstream& ofstr, int ef) const
  {
  string nns("##$");
  ofstr << nns << "TYPK= "       << "<AXIS>"  << "\n";
  ofstr << nns << "TYPNR= "      << "4"       << "\n";
  ofstr << nns << "TYPNAM= "     << _TYPNAM   << "\n";
  ofstr << nns << "PICNUMB= "    << _PICNUMB  << "\n";
  ofstr << nns << "SOFTCLIP= "   << "0"       << "\n";
  ofstr << nns << "DISABLE= "    << "0"       << "\n";
  ofstr << nns << "SHADOW= "     << "0"       << "\n";
  ofstr << nns << "PLANENR= "    << "0x5"     << "\n";
  ofstr << nns << "MODIFY= "     << "0x1"     << "\n";
  ofstr << nns << "SATURATION= " << "0"       << "\n";
  ofstr << nns << "HUE= "        << _HUE      << "\n";
  ofstr << nns << "INTENSITY= "  << "0xff"    << "\n";
  ofstr << nns << "PENCOL= "     << _PENCOL   << "\n";
  ofstr << nns << "XORIGIN= "    << _XORIGIN  << "\n";
  ofstr << nns << "YORIGIN= "    << _YORIGIN  << "\n";
  ofstr << nns << "ZORIGIN= "    << "0"       << "\n";
  ofstr << nns << "XSCALE= "     << "1"       << "\n";
  ofstr << nns << "YSCALE= "     << "1"       << "\n";
  ofstr << nns << "ZSCALE= "     << "1"       << "\n";
  ofstr << nns << "XLEAST= "     << _XLEAST   << "\n";
  ofstr << nns << "YLEAST= "     << _YLEAST   << "\n";
  ofstr << nns << "ZLEAST= "     << "0"       << "\n";
  ofstr << nns << "XHIGHEST= "   << _XHIGHEST << "\n";
  ofstr << nns << "YHIGHEST= "   << _YHIGHEST << "\n";
  ofstr << nns << "ZHIGHEST= "   << "0"       << "\n";
  ofstr << nns << "ANGLE1= "      << _ANGLE1 << "\n";
  ofstr << nns << "ANGLE2= "      << "0"      << "\n";
  ofstr << nns << "ANGLE3= "      << "0"      << "\n";
  ofstr << nns << "AXDLBL= "      << _AXDLBL  << "\n";
  ofstr << nns << "AXLABL= "      << _AXLABL  << "\n";
  ofstr << nns << "AXFLAG= "      << _AXFLAG  << "\n";
  ofstr << nns << "AXFORM= "      << _AXFORM  << "\n";
//  if(!oldflag)
    ofstr << nns << "AXLDIR= "      << _AXLDIR  << "\n";
  ofstr << nns << "AXGRFLAG= "    << "0"      << "\n";
  if(_AXSTRG.length())
    ofstr << nns << "AXSTRG= "    << _AXSTRG  << "\n";
  ofstr << nns << "AXLTIC= "      << "0.2"    << "\n";
  ofstr << nns << "AXDTIC= "      << _AXDTIC  << "\n";
  ofstr << nns << "AXLIM1= "      << _AXLIM1  << "\n";
  ofstr << nns << "AXLIM2= "      << _AXLIM2  << "\n";
  ofstr << nns << "AXCELL= "      << _AXCELL  << "\n";
  ofstr << nns << "AXGDIS= "      << _AXGDIS  << "\n";
  ofstr << nns << "AXGLEN= "      << _AXGLEN  << "\n";
  ofstr << nns << "AXLOOP= "      << _AXLOOP  << "\n";
  if(ef) ofstr << "\n";
  }

void XWinMeta::writePic(ofstream& ofstr, int ef) const
  {
  string nns("##$");
  ofstr << nns << "TYPK= "       << "<PIC>"     << "\n";
  ofstr << nns << "TYPNR= "      << "1"         << "\n";
  ofstr << nns << "TYPNAM= "     << _TYPNAM     << "\n";
  ofstr << nns << "PICNUMB= "    << _PICNUMB    << "\n";
  ofstr << nns << "SOFTCLIP= "   << _SOFTCLIP   << "\n";
  ofstr << nns << "DISABLE= "    << _DISABLE    << "\n";
  ofstr << nns << "SHADOW= "     << "0"         << "\n";
  ofstr << nns << "PLANENR= "    << _PLANENR    << "\n";
  ofstr << nns << "MODIFY= "     << "0x1"       << "\n";
  ofstr << nns << "SATURATION= " << _SATURATION << "\n";
  ofstr << nns << "HUE= "        << _HUE        << "\n";
  ofstr << nns << "INTENSITY= "  << _INTENSITY  << "\n";
  ofstr << nns << "PENCOL= "     << "0"         << "\n";
  ofstr << nns << "XORIGIN= "    << _XORIGIN    << "\n";
  ofstr << nns << "YORIGIN= "    << _YORIGIN    << "\n";
  ofstr << nns << "ZORIGIN= "    << "0"         << "\n";
  ofstr << nns << "XSCALE= "     << _XSCALE     << "\n";
  ofstr << nns << "YSCALE= "     << _YSCALE     << "\n";
  ofstr << nns << "ZSCALE= "     << "1"         << "\n";
  ofstr << nns << "XLEAST= "     << _XLEAST     << "\n";
  ofstr << nns << "YLEAST= "     << _YLEAST     << "\n";
  ofstr << nns << "ZLEAST= "     << "0"         << "\n";
  ofstr << nns << "XHIGHEST= "   << _XHIGHEST   << "\n";
  ofstr << nns << "YHIGHEST= "   << _YHIGHEST   << "\n";
  ofstr << nns << "ZHIGHEST= "   << "0"         << "\n";
  ofstr << nns << "ANGLE1= "     << _ANGLE1     << "\n";
  ofstr << nns << "ANGLE2= "     << "0"         << "\n";
  ofstr << nns << "ANGLE3= "     << "0"         << "\n";
  ofstr << nns << "PIC_ANZ= "    << _PIC_ANZ    << "\n";
  ofstr << nns << "PIC_PRJ= "    << "0"         << "\n";
  ofstr << nns << "PIC_CSC= "    << "0"         << "\n";
  ofstr << nns << "PIC_ZSL= "    << "0"         << "\n";
  if(ef) ofstr << "\n";
  }

void XWinMeta::writeSpec(ofstream& ofstr, int ef) const
  {
  string nns("##$");
  ofstr << nns << "TYPK= "       << "<SPEC>"    << "\n";
  ofstr << nns << "TYPNR= "      << "2"         << "\n";
  ofstr << nns << "TYPNAM= "     << _TYPNAM     << "\n";
  ofstr << nns << "PICNUMB= "    << _PICNUMB    << "\n";
  ofstr << nns << "SOFTCLIP= "   << "0"         << "\n";
  ofstr << nns << "DISABLE= "    << _DISABLE    << "\n";
  ofstr << nns << "SHADOW= "     << "0"         << "\n";
  ofstr << nns << "PLANENR= "    << _PLANENR    << "\n";
  ofstr << nns << "MODIFY= "     << "0x1"       << "\n";
  ofstr << nns << "SATURATION= " << _SATURATION << "\n";
  ofstr << nns << "HUE= "        << _HUE        << "\n";
  ofstr << nns << "INTENSITY= "  << "0xff"      << "\n";
  ofstr << nns << "PENCOL= "     << _PENCOL     << "\n";
  ofstr << nns << "XORIGIN= "    << _XORIGIN    << "\n";
  ofstr << nns << "YORIGIN= "    << _YORIGIN    << "\n";
  ofstr << nns << "ZORIGIN= "    << "0"         << "\n";
  ofstr << nns << "XSCALE= "     << "1"         << "\n";
  ofstr << nns << "YSCALE= "     << _YSCALE     << "\n";
  ofstr << nns << "ZSCALE= "     << "1"         << "\n";
  ofstr << nns << "XLEAST= "     << "0"         << "\n";
  ofstr << nns << "YLEAST= "     << _YLEAST     << "\n";
  ofstr << nns << "ZLEAST= "     << "0"         << "\n";
  ofstr << nns << "XHIGHEST= "   << "0"         << "\n";
  ofstr << nns << "YHIGHEST= "   << _YHIGHEST   << "\n";
  ofstr << nns << "ZHIGHEST= "   << "0"         << "\n";
  ofstr << nns << "ANGLE1= "     << "0"         << "\n";
  ofstr << nns << "ANGLE2= "     << "0"         << "\n";
  ofstr << nns << "ANGLE3= "     << "0"         << "\n";
  ofstr << nns << "FILENAME= "   << _FILENAMES   << "\n";
  ofstr << nns << "D_TYPE= "     << _D_TYPE     << "\n";
  ofstr << nns << "INCREMENT= "  << "1"         << "\n";
  ofstr << nns << "S_OFFSET= "   << "0"         << "\n";
  ofstr << nns << "MARKER= "     << "0"         << "\n";
  ofstr << nns << "PROCESS= "    << "0"         << "\n";
  ofstr << nns << "INTEG= "      << "0"         << "\n";
  ofstr << nns << "STEPS= "      << "0"         << "\n";
  ofstr << nns << "SPFLAG= "     << "0"         << "\n";
  ofstr << nns << "SPCELL= "     << "0"         << "\n";
  ofstr << nns << "SP_SHUFFLE= " << "0"         << "\n";
  ofstr << nns << "SCALFLG= "    << _SCALFLG    << "\n";
  ofstr << nns << "SPASSF= "     << "0"         << "\n";
  ofstr << nns << "SPASSW= "     << "0"         << "\n";
  if(ef) ofstr << "\n";
  }

void XWinMeta::writeText(ofstream& ofstr, int ef) const
  {
  string nns("##$");
  ofstr << nns << "TYPK= "       << "<TEXT>"    << "\n";
  ofstr << nns << "TYPNR= "      << "3"         << "\n";
  ofstr << nns << "TYPNAM= "     << _TYPNAM     << "\n";
  ofstr << nns << "PICNUMB= "    << _PICNUMB    << "\n";
  ofstr << nns << "SOFTCLIP= "   << "0"         << "\n";
  ofstr << nns << "DISABLE= "    << "0"         << "\n";
  ofstr << nns << "SHADOW= "     << "0"         << "\n";
  ofstr << nns << "PLANENR= "    << _PLANENR    << "\n";
  ofstr << nns << "MODIFY= "     << _MODIFY     << "\n";
  ofstr << nns << "SATURATION= " << _SATURATION << "\n";
  ofstr << nns << "HUE= "        << _HUE        << "\n";
  ofstr << nns << "INTENSITY= "  << _INTENSITY  << "\n";
  ofstr << nns << "PENCOL= "     << _PENCOL     << "\n";
  ofstr << nns << "XORIGIN= "    << "0"         << "\n";
  ofstr << nns << "YORIGIN= "    << "0"         << "\n";
  ofstr << nns << "ZORIGIN= "    << "0"         << "\n";
  ofstr << nns << "XSCALE= "     << "1"         << "\n";
  ofstr << nns << "YSCALE= "     << "1"         << "\n";
  ofstr << nns << "ZSCALE= "     << "1"         << "\n";
  ofstr << nns << "XLEAST= "     << "0"         << "\n";
  ofstr << nns << "YLEAST= "     << _YLEAST     << "\n";
  ofstr << nns << "ZLEAST= "     << "0"         << "\n";
  ofstr << nns << "XHIGHEST= "   << _XHIGHEST   << "\n";
  ofstr << nns << "YHIGHEST= "   << _YHIGHEST   << "\n";
  ofstr << nns << "ZHIGHEST= "   << "0"         << "\n";
  ofstr << nns << "ANGLE1= "     << "0"         << "\n";
  ofstr << nns << "ANGLE2= "     << "0"         << "\n";
  ofstr << nns << "ANGLE3= "     << "0"         << "\n";
  ofstr << nns << "FILENAME= "   << _FILENAMES   << "\n";
  ofstr << nns << "TXMODE= "     << _TXMODE     << "\n";
  ofstr << nns << "TXALLI= "     << _TXALLI     << "\n";
  ofstr << nns << "TXFNST= "     << _TXFNST     << "\n";
  ofstr << nns << "TXFNAL= "     << "0"         << "\n";
  ofstr << nns << "TXALFG= "     << "0"         << "\n";
  ofstr << nns << "TXSOSP= "     << _TXSOSP     << "\n";
  ofstr << nns << "TXFILE= "     << "1"         << "\n";
  ofstr << nns << "TXCELL= "     << _TXCELL     << "\n";
  ofstr << nns << "TXSLNT= "     << "0"         << "\n";
  if(ef) ofstr << "\n";
  }

void XWinMeta::writePeak(ofstream& ofstr, int ef) const
  {
  string nns("##$");
  writeFirst(ofstr);		// Write 1st block
  writeDraw(ofstr);		// Write drawing block
  writeXYZ(ofstr);		// Write {X,Y,Z} block
  ofstr << nns << "FILENAME= "   << _FILENAMES   << "\n";
  ofstr << nns << "PLROTA= "     << _PLROTA     << "\n";
  ofstr << nns << "PLDIWI= "     << _PLDIWI     << "\n";
  ofstr << nns << "PLMARK= "     << "0"         << "\n";
  ofstr << nns << "PLSIGN= "     << "0"         << "\n";
  ofstr << nns << "PLOCHA= "     << "0"         << "\n";
  if(_PLLABF.length())
    ofstr << nns << "PLLABF= "   << _PLLABF     << "\n";
  ofstr << nns << "PLSTRG= "     << _PLSTRG     << "\n";
  ofstr << nns << "PLCELL= "     << "-4"        << "\n";
  ofstr << nns << "PLOFFS= "     << "0"         << "\n";
  ofstr << nns << "PLSCAL= "     << "0"         << "\n";
  ofstr << nns << "PLDHEI= "     << _PLDHEI     << "\n";
  ofstr << nns << "PLMSHI= "     << _PLMSHI     << "\n";
  ofstr << nns << "PLINVL= "     << "2000000"   << "\n";
  ofstr << nns << "PLINSH= "     << "5000000"   << "\n";
  ofstr << nns << "PLMINT= "     << "0"         << "\n";
  if(ef) ofstr << "\n";
  }

void XWinMeta::writeImag(ofstream& ofstr, int ef) const
  {
  string nns("##$");
  writeFirst(ofstr);		// Write 1st block
  writeDraw(ofstr);		// Write drawing block
  writeXYZ(ofstr);		// Write {X,Y,Z} block
  ofstr << nns << "FILENAME= "   << _FILENAMES   << "\n";
  ofstr << nns << "B_MAX= "      << "2147483647" << "\n";
  ofstr << nns << "B_MIN= "      << "0"         << "\n";
  ofstr << nns << "UT= "         << "0"         << "\n";
  ofstr << nns << "LT= "         << "0"         << "\n";
  ofstr << nns << "CONULE= "     << "3"         << "\n";
  ofstr << nns << "CONUCO= "     << "4"         << "\n";
  ofstr << nns << "CONTRI= "     << _CONTRI     << "\n";
  ofstr << nns << "CONCOL= "     << "(0..9)"    << "\n";
  ofstr << "0 115 192 0 9\n"
        << "0 213 255 0 8\n"
        << "0 170 255 0 7\n"
        << "255 128 255 0 0\n"
        << "0 220 255 0 1\n"
        << "0 32 255 0 2\n"
        << "0 54 255 0 3\n"
        << "0 85 255 0 4\n"
        << "0 120 255 0 5\n"
        << "0 130 255 0 6\n";
  ofstr << nns << "CONTOUR= "    << "0"         << "\n";
  ofstr << nns << "SMX_DTYPE= "  << "0"         << "\n";
  ofstr << nns << "SCALE_TYP= "  << "1"         << "\n";
  ofstr << nns << "SMX_X1DIM= "  << "0"         << "\n";
  ofstr << nns << "SMX_X2DIM= "  << "0"         << "\n";
  ofstr << nns << "SMX_X3DIM= "  << "0"         << "\n";
  if(ef) ofstr << "\n";
  }

void XWinMeta::writeFirst(ofstream& ofstr) const
  {
  string nns("##$");
  ofstr << nns << "TYPK= "       << _TYPK       << "\n";
  ofstr << nns << "TYPNR= "      << _TYPNR      << "\n";
  ofstr << nns << "TYPNAM= "     << _TYPNAM     << "\n";
  }

void XWinMeta::writeDraw(ofstream& ofstr) const
  {
  string nns("##$");
  ofstr << nns << "PICNUMB= "    << _PICNUMB    << "\n";
  ofstr << nns << "SOFTCLIP= "   << _SOFTCLIP   << "\n";
  ofstr << nns << "DISABLE= "    << _DISABLE    << "\n";
  ofstr << nns << "SHADOW= "     << _SHADOW     << "\n";
  ofstr << nns << "PLANENR= "    << _PLANENR    << "\n";
  ofstr << nns << "MODIFY= "     << _MODIFY     << "\n";
  ofstr << nns << "SATURATION= " << _SATURATION << "\n";
  ofstr << nns << "HUE= "        << _HUE        << "\n";
  ofstr << nns << "INTENSITY= "  << _INTENSITY  << "\n";
  ofstr << nns << "PENCOL= "     << _PENCOL     << "\n";
  }

void XWinMeta::writeXYZ(ofstream& ofstr) const
  {
  string nns("##$");
  ofstr << nns << "XORIGIN= "    << _XORIGIN    << "\n";
  ofstr << nns << "YORIGIN= "    << _YORIGIN    << "\n";
  ofstr << nns << "ZORIGIN= "    << _ZORIGIN    << "\n";
  ofstr << nns << "XSCALE= "     << _XSCALE     << "\n";
  ofstr << nns << "YSCALE= "     << _YSCALE     << "\n";
  ofstr << nns << "ZSCALE= "     << _ZSCALE     << "\n";
  ofstr << nns << "XLEAST= "     << _XLEAST     << "\n";
  ofstr << nns << "YLEAST= "     << _YLEAST     << "\n";
  ofstr << nns << "ZLEAST= "     << _ZLEAST     << "\n";
  ofstr << nns << "XHIGHEST= "   << _XHIGHEST   << "\n";
  ofstr << nns << "YHIGHEST= "   << _YHIGHEST   << "\n";
  ofstr << nns << "ZHIGHEST= "   << _ZHIGHEST   << "\n";
  ofstr << nns << "ANGLE1= "     << _ANGLE1     << "\n";
  ofstr << nns << "ANGLE2= "     << _ANGLE2     << "\n";
  ofstr << nns << "ANGLE3= "     << _ANGLE3     << "\n";
  }
#endif							// XWinMeta.cc\n";
