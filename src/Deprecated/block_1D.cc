/*************************************************************************
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      Block_1D                                      Implementation    **
**                                                                      **
**      Copyright (c) 1990, 1991, 1992                                  **
**      Scott Smith                                                     **
**      Eidgenoessische Technische Hochschule                           **
**      Labor fuer physikalische Chemie                                 **
**      8092 Zurich / Switzerland                                       **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** This class specifies the structure and properties of a 1-dimensional **
** block of simulated (or any other) data.  The data block is a row     **
** vector in RAM memory and has	an associated parameter set.  Functions **
** are provided for basic algebraic manipulations, specification of	**
** parameters, and general I/O.  Many external programs deal with data	**
** sets coupled to a block of parameters and this class allows users    **
** to more readily support I/O to these programs.  One example would be **
** the Felix module which provides functions for input/output of GAMMA  **
** data to/from the program Felix from Molecular Simulations (MSI),     **
** nee Biosym, nee Hare Research, Inc.                                  **
**                                                                      **
*************************************************************************/

#ifndef   block_1D_cc_			// Is file already included?
#  define block_1D_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Deprecated/block_1D.h>	// Include class interface
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <assert.h>
#include <iostream>
using namespace std;
 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      CLASS BLOCK_1D ERROR HANDLING
// ____________________________________________________________________________


void block_1D::Blk1DError(int eidx, int noret) const
      
	// Input		BLK 	: A 1D data block (this)
        //                      eidx    : Error index   
        //                      noret   : Flag for linefeed (0=linefeed)
        // Output               none    : Error message output 
 
/* The following error messages use the defaults set in the Gutils package

                Case                          Error Message                

                (0)                     Program Aborting.....
                (9)                     Problems During Construction
                default                 Unknown Error                        */
 
  {
  std::string hdr("Block 1D");
  switch (eidx)
    {
    case 1:  GAMMAerror(hdr,"Mixing Blocks of Differing Dimensions",noret);	break;//(1)
    case 2:  GAMMAerror(hdr,"Problems With Block-Block Operation",noret);	break;//(2)
    case 5:  GAMMAerror(hdr,"Block Data Point NOT Accessible",noret);		break;//(5)
    case 6:  GAMMAerror(hdr,"Accessed Data Point EXCEEDS Block Size",noret);	break;//(6)
    case 50: GAMMAerror(hdr,"Error in Block-Block Addition",noret);		break;//(50)
    case 51: GAMMAerror(hdr,"Error in Block-Block Subtraction",noret);		break;//(51)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }


void block_1D::Blk1DError(int eidx, int i, int j, int noret) const
     
	// Input		BLK 	: A 1D data block (this)
	//                      eidx    : Error index
        //                      i,j     : Block indices
        //                      noret   : Flag for linefeed (0=linefeed)
        // Output               none    : Error message output

  {
  std::cout << "\nClass Block_1D: ";
  switch (eidx)
    {
    case 100:							// (100)
      std::cout << "Size = " << i << " and Current Point Accessed = " << j;
      break;
    case 101:							// (101)
      std::cout << "Block Data Size = " << i << " mixed with Block Data Size " << j;
      break;
    default:
      std::cout << "Unkown Error, number " << eidx;
    }
  if(!noret) std::cout << "\n";
  }


volatile void block_1D::Blk1DFatal(int eidx)

	// Input		BLK 	: A 1D data block (this)
        //                      eidx    : Error index
        // Output               none    : Error message output
        //                                Program execution stopped

  {
  Blk1DError(eidx,1);
  if(eidx) Blk1DError(0);
  GAMMAfatal();					// Clean exit from program
  }
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                 CLASS BLOCK_1D CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

	
block_1D::block_1D(int dim) : row_vector(dim,0), ParameterSet()

	// Input			dim  : integer
	// Output			none : creates an empty block_1D
	//				       of dimension specified by dim

  { set_defaults(1); }			// Set Block Parameter Defaults

// ____________________________________________________________________________
// B                  CLASS BLOCK_1D FUNCTIONS, BLOCK WITH BLOCK
// ____________________________________________________________________________

	
block_1D operator + (const block_1D &BLK1, const block_1D &BLK2)

	// Input		BLK1, BLK2 : two block_1Ds
	// Return		BLK        : block_1D which is the addition
        //	                             of the two input block_1Ds,
	//			             BLK =  BLK1 + BLK2

  { 
  block_1D BLK(BLK1);				// Construct form BLK1
  BLK.block_1D::operator+= (BLK2);		// Add in BLK2
  return BLK;
  }


void block_1D::operator += (const block_1D &BLK2)

	// Input		BLK1 : a block_1D
	// Return		BLK  : block_1D which is the input
        //	                       input block_1D BLK1 added to itself
	//			       BLK = BLK + BLK1

  {
//                 First Attempt To Add the Data Values

  if(BLK2.elements())				// Add BLK2 data if exists
    {
    if(!elements()) row_vector::operator=(BLK2);// Data is BLK2 if BLK1 empty
    else					// If BLK2 & BLK1 have data
      {						// then
      if(Blk1DCheck(BLK2))		//   Insure dimensions match
        row_vector::operator+= (BLK2);		//   Add data vectors
      else
        {
        Blk1DError(2);			// error in BLK with BLK operation
        Blk1DFatal(50);			// fatality during BLK + BLK
        }
      }
    }
// sosi - kludge for gamma 4.0
//  ParameterSet::operator += (BLK2);			// add BLK2 parameters
  return;
  }


block_1D operator - (const block_1D &BLK1, const block_1D &BLK2)

	// Input		BLK1 : a data block_1D
	//			BLK2 : a data block_1D
	// Return		BLK  : block_1D which is the subtraction
        //	                       of the two input block_1Ds,
	//			       BLK = BLK1 - BLK2

  {
  block_1D BLK(BLK1);
  BLK.block_1D::operator-= (BLK2);
  return BLK;
  }

block_1D operator - (const block_1D &BLK1)

	// Input		BLK1  : a block_1D
	// Return		BLK   : block_1D which is the negative
  //	                of the input block_1D, BLK = - BLK1

  {
  block_1D BLK;
	
#ifdef _MSC_VER
	// *DCT* This call leads to infinite recursion and will not compile on windows.
	cout << "\nFATAL: Call to this deprecated funtion leads to infinite recursion.\n";
	assert(0);
#else
	BLK.row_vector::operator= (-BLK1);	// vector negation
  BLK.ParameterSet::operator= (BLK1);		// copy parameter set BLK1 into BLK
#endif

  return BLK;
  }



void block_1D::operator -= (const block_1D &BLK2)

	// Input		BLK1 : a block_1D
	// 			BLK2 : a block_1D
	// Return		none : block_1D BLK1 is modified to have
        //	                       input block_1D BLK2 added to itself
	//			       BLK1 = BLK1 - BLK2

  {
  if(BLK2.elements())				// if BLK1 has data
    {
    if(!elements())				// if BLK1 has no data
      row_vector::operator= (-BLK2);		// vector negation
    else
      {
      if(Blk1DCheck(BLK2))		// insure dimensions match
        row_vector::operator-= (BLK2);		// subtract vectors
      else
        {
        Blk1DError(2);			// error in BLK with BLK operation
        Blk1DFatal(51);			// fatality during BLK - BLK
        }
      }
    }
  }


void block_1D::operator = (const block_1D &BLK)

	// Input		this,BLK  : two block_1D
	// Return		copies BLK to this, and destroys the old
        //                      contents of this block_1D

  {
  row_vector::operator= (BLK);	// Copy BLK vector
  ParameterSet::operator= (BLK);	// Copy BLK parameter set
  }


void block_1D::operator = (const matrix &mx)
// sosi - dimension check needed?

	// Input		this,BLK  : two block_1D
	// Return		copies BLK to this, and destroys the old
        //                      contents of this block_1D

  {
  row_vector::operator= (mx);	// Copy BLK vector
  ParameterSet::operator= (ParameterSet());	// Copy BLK parameter set
  }


// ______________________________________________________________________
//                CLASS BLOCK_1D FUNCTIONS, BLOCK WITH SCALAR
// ______________________________________________________________________

	
block_1D operator * (const block_1D &BLK1, const complex &z)

	// Input		BLK1 : a data block
	//			z    : a complex number
	// Return		BLK  : a new data block which is the input
	//			       block multiplied by scalar z

  { 
  block_1D BLK(BLK1);
  BLK *= z;
  return BLK;
  }				// Uses class row_vector

	
block_1D operator * (const block_1D &BLK1, double d)

	// Input		BLK1 : a data block
	//			z    : a complex number
	// Return		BLK  : a new data block which is the input
	//			       block multiplied by scalar z

  {
  block_1D BLK(BLK1);
  BLK *= d;
  return BLK;
  }				// Uses class row_vector


block_1D operator * (const complex &z, const block_1D &BLK1)

	// Input		BLK1 : a data block
	//			z    : a complex number
	// Return		BLK  : a new data block which is the input
	//			       block multiplied by scalar z

  {
  block_1D BLK(BLK1);
  BLK *= z;
  return BLK;
  }				// Uses class row_vector


block_1D operator * (double d, const block_1D &BLK1)

	// Input		BLK1 : a data block
	//			z    : a complex number
	// Return		BLK  : a new data block which is the input
	//			       block multiplied by scalar z

  {
  block_1D BLK(BLK1);
  BLK *= d;
  return BLK;
  } 				// Uses class row_vector


block_1D operator / (const block_1D &BLK1, const complex &z)

	// Input		BLK  : a data block
	//			z    : a complex number
	// Return		BLK1 : a new data block which is the input
	//			       block divided by scalar z

  {
  block_1D BLK(BLK1);
  BLK /= z;
  return BLK;
  }				// Uses class row_vector


// ______________________________________________________________________
//             CLASS BLOCK_1D FUNCTIONS, BLOCK WITH PARAMETER SET
// ______________________________________________________________________


//void operator += (block_1D &BLK, const ParameterSet &pset)

	// Input		BLK  : a data block
	//			pset : a parameter set
	// Return		none : the input data block is returned
	//			       with parameters from pset added
     
//  { (ParameterSet)BLK += pset;	}		// Add parameters in pset
// sosi: this does not work!, I don't know why yet.


// ______________________________________________________________________
//                  CLASS BLOCK_1D DATA SET FUNCTIONS
// ______________________________________________________________________


void product(block_1D &BLK1, const block_1D &BLK2)

	// Input		BLK1 : data block
	// 			BLK2 : data block
	// Return		none : BLK1 multiplied point by point by BLK2

  {
  BLK1 = block_1D(product((row_vector&)BLK1, (const row_vector&)BLK2));
//  block_1D BLKprod = block_1D(product((row_vector&)BLK1, (const row_vector&)BLK2));
//  BLK1.row_vector::operator= row_vector(BLKprod);
// sosi kludge
//  BLK1.ParameterSet::operator += ParameterSet(BLK2);	// add BLK2 parameters
  return;
  }


// ______________________________________________________________________
//                  CLASS BLOCK_1D PARAMETER SET FUNCTIONS
// ______________________________________________________________________

// *****************  Entire Parameter Set Manipulations ****************

void block_1D::set_defaults(int level)

	// Input		BLK  : data block(this)
	//			level: level to which defaults are set
	// Return		none : BLK parameter set defaults set

  {
  if(level == 0) clear();		// Zero entire parameter set
  std::string pname("Dim");
  std::string pstate("Data Set Dimension of Block");
  SinglePar par(pname, elements(), pstate);
  push_back(par);			// Add # dimension to pset
  }


// ______________________________________________________________________
//                        CLASS BLOCK_1D CHECKS
// ______________________________________________________________________


int block_1D::Blk1DCheck(int pt) const

  {
  if(pt >= elements())			// Insure point not beyond vector
    {
    Blk1DError(6);			// Point exceeds block data size
    Blk1DError(100, elements(), pt);
    return 0;
    }
  return 1;
  }


int block_1D::Blk1DCheck(const block_1D &BLK2) const

  {
  if(elements() != BLK2.elements())	// Dimension equivalence check
    {
    Blk1DError(1,1);			// Dimension mismatch
    Blk1DError(101, elements(), BLK2.elements());
    return 0;
    }
  return 1;
  }



// ____________________________________________________________________________ 
// i                      CLASS BLOCK_1D ERROR HANDLING
// ____________________________________________________________________________
 
 
void block_1D::write_pset(const std::string& filename) const

	// Input		BLK  : A 1D data block (this)
	//			file : Name of disk file to produce
	// Return		none : A file is produced on disk containing
	//			       the parameter set in BLK

  { ParameterSet::write(filename); }


std::ostream& operator << (std::ostream& ostr, const block_1D &BLK)

	// Input		BLK  : data block
	// Return		none : data contained in BLK sent
	//			       to output stream

  { return ostr << (const row_vector&)BLK; }


void print_pset(const block_1D &BLK)

	// Input		BLK  : data block
	// Return		none : parameter set contained in BLK sent
	//			       to standard I/O	

  {
  std::cout << "\n\n\t\t\t\tCurrent Block_1D Parameter Set\n\n";
  std::cout << (const ParameterSet&)BLK;
  }


void print_dset(block_1D &BLK)

	// Input		BLK  : data block
	// Return		none : BLK data block is set to standard output

  {
  std::cout << "\n\n\t\t\t\tCurrent Block_1D Data Vector\n\n";
  int npts = BLK.elements();
  int pt = 0;
  int flag = 0;
  while (pt < npts)
    {
      std::cout << pt+1 << ".\t" << BLK(pt);
      flag = !flag;
      if (flag)
	std::cout << "\t";
      else
	std::cout << "\n";
      pt++;
    }
  }


#endif 							// Block_1D.cc
