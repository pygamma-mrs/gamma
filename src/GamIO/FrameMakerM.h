/* FrameMakerM.h *************************************************-*-c++-*-
**									**
**	                          G A M M A 				**
**								 	**
**	FrameMaker MIF Functions		Interface		**
**						 			**
**	Copyright (c) 1991, 1992				 	**
**	Tilo Levante, Scott Smith					**
**	Eidgenoessische Technische Hochschule	 			**
**	Labor fuer physikalische Chemie				 	**
**	8092 Zurich / Switzerland				 	**
**								 	**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**								 	**
** The FrameMaker Library provides functions to create different 	**
** structures readable in the desktop publishing program FrameMaker.	**
** This module is part of the FrameMaker set in GAMMA and contains the	**
** functions that output common MIF constructs.				**
**								 	**
*************************************************************************/

#ifndef   GFrameMakerM_h_		// Is this file already included?
#  define GFrameMakerM_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <GamIO/FrameMakerP.h>		// Know about FM parameters

// ____________________________________________________________________________
//                   FRAMEMAKER MIF BASIC CONSTRUCTS 
// ____________________________________________________________________________


MSVCDLL int FM_ID();

	// Input	none	: none
	// Output	int	: The internal FM object ID
	//			  is incremented and returned
  

MSVCDLL void FM_Begin(std::ostream &out);

	// Input	ostr    : An output stream
	// Output	void	: The output stream is modified
	//			  to contain a comment about its origin


MSVCDLL void FM_End(std::ostream &out);

	// Input	ostr    : An output stream
	// Output	void	: The output stream is modified
	//			  to contain a comment about its end


MSVCDLL void FM_Group (std::ostream &out, int ID, int GroupID=0);

	// Input	out       : output stream
	// 		ID        : ID number
	// 		GroupID   : Group ID number
	// Output	none      : Function is void.  Puts the commands
	//			    to group all graphic objects with
	//			    i.d. ID together with i.d. GroupID
 
// ____________________________________________________________________________
// B                 FRAMEMAKER MIF ANCHORED FRAME CONSTUCTS
// ____________________________________________________________________________
 

MSVCDLL void FM_AFrames_Begin(std::ostream &out);
 
        // Input        ostr    : An output stream
        // Output       void    : The output stream is modified
        //                        to contain an anchored frame start


MSVCDLL void FM_AFrames_End(std::ostream &out);

        // Input        ostr    : An output stream
        // Output       void    : The output stream is modified
        //                        to contain an anchored frame end
 

MSVCDLL void FM_AFrame_Set(std::ostream &out, double xsize, double ysize, int ID=11);

        // Input        ostr    : An output stream                       
        //              xsize   : Frame horizontal dimensions in cm
        //              ysize   : Frame vertical dimensions in cm
        //              ID      : Anchored frame ID 
        // Output       void    : The output stream is modified
        //                        to contain anchored frame settings
 

MSVCDLL void FM_AFrame_Set(std::ostream &out, const FMPar& FMP, int ID=11);

        // Input	out	: Open output stream
	//		FMP	: FrameMaker plot parameters
	//		ID	: FrameMaker ID value
	// Output	void	: The output stream is modified to
	//			  contain the start of an anchored
	//			  in FrameMaker MIF format
	// Note			: The anchored frame height & width
	//			  in cm is taken from FMP
	// Note			: The anchored frame will have ID


MSVCDLL void FM_AFrame_End(std::ostream &out);
 
        // Input        ostr    : An output stream 
        // Output       void    : The output stream is modified 
        //                        to contain end of anchored frame
 
// ____________________________________________________________________________
// C                    FRAMEMAKER MIF TEXT FLOW CONSTUCTS
// ____________________________________________________________________________


MSVCDLL void FM_TextFlow_Set(std::ostream &out);


MSVCDLL void FM_TextFlow_End(std::ostream &out);


MSVCDLL void FM_Paragraph_Set(std::ostream &out);


MSVCDLL void FM_Paragraph_End (std::ostream &out, int ID=11);


MSVCDLL void FM_ParaText_End(std::ostream &out);


MSVCDLL void FM_TextLine(std::ostream &out, int ID, double x, double y, char z, 
                                  double size=0.0, int align=0, int angle=0);

	// Input	out       : output stream
	// 		ID        : Text group ID number (0 = no group)
	//		x	  : Text x-alignment coordinate
	//		y	  : Text y-alignment coordinate
	//		z	  : Character to be put in TextLine
	// 		size 	  : Font size in points
	// 		align     : Arrow on line (TRUE = yes, FALSE = no)
	// 		angle     : Angle flag (0=center, <0=Left, >0=Right)
	// Output	none      : Function is void.  Puts the commands
	//			    to produce a FrameMaker TextLine in MIF format
	// Note			  : Although MIF character format(s) can be specified
	//			    here it is left for user to do so within FM itself,
	//			    except for the font size.
	// Note			  : Additional generic object data can be included here
	//			    as specified in the Frame manual "MIF Reference"
	//			    such as color and pen pattern.  Not in currently


MSVCDLL void FM_TextLine(std::ostream &out, int ID, int al, double x, double y, double z);


MSVCDLL void FM_TextLine(std::ostream &out, int ID, double x, double y, 
                      std::string text, double size=0.0, int align=0, int angle=0);

	// Input	out       : output stream
	// 		ID        : Text group ID number (0 = no group)
	//		x	  : Text x-alignment coordinate
	//		y	  : Text y-alignment coordinate
	//		text	  : string to be put in TextLine
	// 		size 	  : Font size in points
	// 		arrow     : Arrow on line (TRUE = yes, FALSE = no)
	// 		angle     : Angle flag (0=center, <0=Left, >0=Right)
	// Output	none      : Function is void.  Puts the commands
	//			    to produce a FrameMaker TextLine in MIF format
	// Note			  : Although MIF character format(s) can be specified
	//			    here it is left for user to do so within FM itself,
	//			    except for the font size.
	// Note			  : Additional generic object data can be included here
	//			    as specified in the Frame manual "MIF Reference"
	//			    such as color and pen pattern.  Not in currently

// ____________________________________________________________________________
// D                    FRAMEMAKER MIF TABLE CONSTUCTS
// ____________________________________________________________________________


MSVCDLL void FM_Tbl_Begin(std::ostream &out);

	// Input	out       : output stream
	// Output	none      : Function is void.  The output stream
	//			    is modified to contain a MIF declaration
	//			    of a table for FrameMaker


MSVCDLL void FM_Tbl_End(std::ostream &out);

	// Input	out       : output stream
	// Output	none      : Function is void.  The output stream
	//			    is modified to contain a MIF declaration
	//			    of a table end for FrameMaker


MSVCDLL void FM_TblBody_Begin(std::ostream &out);


MSVCDLL void FM_TblBody_End(std::ostream &out);


MSVCDLL void FM_Tbl_Title(std::ostream &out);


// ____________________________________________________________________________
// E                       FRAMEMAKER MIF GRAPIC OBJECTS 
// ____________________________________________________________________________


MSVCDLL void FM_Ellipse(std::ostream &out, double x, double y, double rx, double ry,
	      int GID=0, int fill=7, int color=0, int pen=0, double width=0.0);

	// Input	out       : output stream
	// 		group     : Line group ID number (0 = no group)
	// 		arrow     : Arrow on line (TRUE = yes, FALSE = no)
	// 		width     : Line width in points
	//		fill      : Ellipse fill pattern (7=solid white)
	//		x	  : Ellipse x coordinate
	//		y	  : Ellipse y coordinate
	//		rx	  : Ellipse width (cm)
	//		ry	  : Ellipse height (cm)
	// Output	none      : Function is void.  Puts the commands
	//			    to produce a FrameMaker Ellipse in MIF
	//			    format.


MSVCDLL void FM_Line(std::ostream &out, int group, int arrow, double width,
	                double x1, double y1,double x2, double y2, int pen=0);

	// Input	out       : output stream
	// 		group     : Line group ID number (0 = no group)
	// 		arrow     : Arrow on line (TRUE = yes, FALSE = no)
	// 		width     : Line width in points
	// 		x1,y1     : Initial coordinates
	// 		x2,y2     : Final coordinates
	// 		pen       : Pen type
	// Output	none      : Function is void.  Puts the commands
	//			    to produce a FrameMaker Line in MIF
	//			    format.  Line is an arrow from (x1,y1) to
	//			    (x2,y2), the coordinates relative to the
	//			    current Frame.


MSVCDLL void FM_PolyLine(std::ostream &out, const row_vector &vx, int ID, int fill,
				     int points, int pen = 0, int width = 1);

	// Input	out       : output stream
	// 		vx        : vector of x,y values (Re, Im)
	// 		ID        : PolyLine ID number
	// 		fill      : PolyLine fill type
	// 		points    : Number of points in PolyLine
	// 		pen       : Pen thickness
	// 		width     : Pen width
	// Output	none      : Function is void.  Puts the commands
	//			    to produce a FrameMaker PolyLine.



MSVCDLL void FM_Rectangle(std::ostream &out, int ID,
	 		   double x, double y, double width, double height);

	// Input	out       : output stream
	// 		ID        : Rectangle ID number
	//		x	  : Rectangle x coordinate
	//		y	  : Rectangle y coordinate
	//		width	  : Rectangle width
	//		height	  : Rectangle height
	// Output	none      : Function is void.  Puts the commands
	//			    to produce a FrameMaker Rectangle.



MSVCDLL void FM_Polygon(std::ostream &out, int ID, double x, double y,
					 double rx, int sides=0, int fill=7);

	// Input	out       : output stream
	// 		ID        : Polygon ID number
	//		x	  : Polygon x coordinate
	//		y	  : Polygon y coordinate
	//		rx	  : Polygon radius
	//		sides     : Polygon sides
	//		fill      : Polygon fill (default to solid white)
	// Output	none      : Function is void.  Puts the commands
	//			    to produce a FrameMaker Rectangle.
	// Note			  : The maximum allowed is an octagon
	//			    as set in this function by "maxsides"


#endif 						// FrameMakerM.h
