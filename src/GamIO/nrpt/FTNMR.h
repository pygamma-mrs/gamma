/* FTNMR.h **********************************************-*-c++-*-
**                                                              **
**                            G A M M A                         **
**                                                              **
**      FTNMR                            Interface		**
**                                                              **
**      Copyright (c) 1990                                      **
**      Scott Smith                                             **
**      Eidgenoessische Technische Hochschule                   **
**      Labor fuer physikalische Chemie                         **
**      8092 Zurich / Switzerland                               **
**                                                              **
**      $Header: $
**                                                              **
*****************************************************************/

/*****************************************************************
**                                                              **
**  Description                                                 **
**                                                              **
**  The FTNMR module provides functions for the input and       **
**  output of GAMMA data to/from the program FTNMR.  FTNMR is   **
**  an old program from Hare Research, Inc. for working up      **
**  NMR spectral data.  It has now been replaced by Felix,      **
**  and owned by Biosym. If you use FTNMR have fun with this    **
**  stuff, but don't expect any upgrades to these routines      **
**  unless you do them yourself.  These won't be supported      **
**                                                              **
*****************************************************************/

#ifndef _FTNMR_h_		// Is file already included?
#define _FTNMR_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface		// This is the interface
#endif

#include <Deprecated/block_1D.h>
#include <Deprecated/block_2D.h>

void FTNMR_smx(const char *filename, block_2D &BLK, int rc=0);

	// Input		filename : output file name
	//			BLK	 : 2D-data block
	//			rc       : flag for real versus complex
	// Output		none     : disk file filename is produced
	//				 : in FTNMR ".smx" format


block_2D FTNMR_smx(const char *filename);

	// Input		filename : output file name
	// Output		BLK	 : 2D-data block read from 
	// 				   disk file filename which
	//				   exists in FTNMR ".smx" format


void FTNMR_dat(const char *filename, block_1D &BLK, int rc=0);

	// Input		filename : output file name
	//			BLK	 : 1D-data block
	//			rc       : flag for real versus complex
	// Output		none     : disk file filename is produced
	//				 : in FTNMR ".dat" format


//block_1D FTNMR_dat(const char *filename);

	// Input		filename : output file name
	// Output		BLK	 : 1D-data block read from 
	// 				   disk file filename which
	//				   exists in FTNMR ".dat" format


#endif /* __FTNMR_H__ */
