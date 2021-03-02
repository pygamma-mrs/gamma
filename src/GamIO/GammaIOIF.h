/* GammaIOIF.h **************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**	GAMMA Level 1 Computation Library 		Interface 	**
**								 	**
**      Scott Smith                                                     **
**	Copyright (c) 1999						**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**						 			**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**								 	**
** This file contains all include statements associated with the GAMMA	**
** modules for interfacing with know external programs.  The programs	**
** are generally either plotting programs, magnetic resonance processing**
** programs, etc. This includes the modules with support output to 	**
** FrameMaker (document processing) and MATLAB (math manipulations)	**
** In addition there is a module for reading in lists of values from	**
** simple ASCII files.							**
**								 	**
*************************************************************************/

#ifndef __GammaIOIF_H__				// Is this already included?
#define __GammaIOIF_H__				// If no, then include it


//                          Gnuplot Support

#include <GamIO/Ggnuplot.h>			// Gnuplot interface
#include <GamIO/GgnuplotC.h>			// Gnuplot controls
#include <GamIO/Ggnuplot1D.h>			// Gnuplot X vs y plots
#include <GamIO/GgnuplotSph.h>			// Gnuplot spherical 3D plot

//                  MATLAB MAT Version 4/5 Support

#include <GamIO/BinIOBase.h>			// Binary I/O Byte Swapping
#include <GamIO/ML4Tag.h>			// MATLAB MAT Headers/Tags
#include <GamIO/ML4DElem.h>			// MATLAB MAT Data Element
#include <GamIO/ML5Hdr.h>			// MATLAB MAT Main Headers
#include <GamIO/ML5Tag.h>			// MATLAB MAT Sub-Element Tags
#include <GamIO/ML5SubE.h>			// MATLAB MAT Sub-Elements
#include <GamIO/ML5AF.h>			// MATLAB MAT Array Flags (SE)
#include <GamIO/ML5DA.h>			// MATLAB MAT Dimen. Array (SE)
#include <GamIO/ML5AN.h>			// MATLAB MAT Array Name (SE)
#include <GamIO/ML5Reals.h>			// MATLAB MAT Reals Array (SE)
#include <GamIO/ML5Imags.h>			// MATLAB MAT Imags Array (SE)
#include <GamIO/ML5DElem.h>			// MATLAB MAT Data Element
#include <GamIO/MatLabFile.h>			// MATLAB MAT File

//                    FrameMaker MIF Support

#include <GamIO/FrameMaker.h>			// FrameMaker interface
#include <GamIO/FrameMakerC.h>			// FrameMaker MIF constructs
#include <GamIO/FrameMakerM.h>
#include <GamIO/FrameMakerP.h>
#include <GamIO/FrameMakerS.h>			// FrameMaker Stack Plots
#include <GamIO/FrameMakerSph.h>		// FrameMaker Spherical Plots

//                    FrameMaker MIF Support

#include <GamIO/Gascii.h>			// Include ASCII file interface

//            Felix (Molecular Simulations) Support

//#include <GamIO/Felix.h>			// Include Felix files

//                    Bruker XWinNMR Support

#include <GamIO/XWin1D.h>			// 1D spectra I/O 
#include <GamIO/XWin2D.h>			// 2D spectra I/O
#include <GamIO/XWinAcqPar.h>			// Acquisition parameters 1D
#include <GamIO/XWinAcqu2s.h>			// Acquisition parameters 2D
#include <GamIO/XWinAcqus.h>			//
#include <GamIO/XWinFid.h>			//
#include <GamIO/XWinMeta.h>			// Bruker XWinNMR Meta Files
#include <GamIO/XWinOutd.h>			//
#include <GamIO/XWinPSet.h>			//
#include <GamIO/XWinProc2s.h>			// Bruker 2D procpar
#include <GamIO/XWinProcPar.h>			// Bruker procpar files
#include <GamIO/XWinProcs.h>			// Bruker 1D procpar
#include <GamIO/XWinSer.h>			// 2D Serial files
#include <GamIO/XWinSpec.h>

#endif 						// GammaIOIF.h
