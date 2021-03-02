/* FTNMR.cc *********************************************-*-c++-*-
**								**
**	                      G A M M A				**
**							 	**
**	FTNMR		                 Implementation  	**
**					      		 	**
**	Copyright (c) 1990				 	**
**	Scott Smith 	          		 		**
**	Eidgenoessische Technische Hochschule	 		**
**	Labor fuer physikalische Chemie		 		**
**	8092 Zurich / Switzerland		 		**
**						 		**
**      $Header: $
**						 		**
*****************************************************************/

/*****************************************************************
**							 	**
**  Description						 	**
**							 	**
**  The FTNMR module provides functions for the input and	**
**  output of GAMMA data to/from the program FTNMR.  FTNMR is	**
**  an old program from Hare Research, Inc. for working	up 	**
**  NMR spectral data.  It has now been	replaced by Felix,	**
**  and owned by Biosym. If you use FTNMR have fun with this    **
**  stuff, but don't expect any upgrades to these routines	**
**  unless you do them yourself.  These won't be supported	**
**							 	**
*****************************************************************/

#ifndef _FTNMR_cc_		// Is file already included?
#define _FTNMR_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#endif

#include <GamIO/FTNMR.h>

void FTNMR_error (int error)
     
	// Input		error : Error Flag
	// Output		none  : Error Message Output

{
  cout << "\nFTNMR: ";
  switch (error)
    {
    case 0:							// (0)
      cout << "Program Aborting.";
      break;
    case 1:							// (1)
      cout << "|yinc| too Large!";
      break;
    default:
      cout << "Unkown Error, number " << error;
    }
  cout << "\n";
}


volatile void FTNMR_fatality (int error)

	// Input		error : error flag
	// Output		none  : Stops Execution & Error Message Output

{
  FTNMR_error (error);
  if(error)
    FTNMR_error (0);
  GAMMAfatal();					// Clean exit from program
}


void FTNMR_dat(const char *filename, block_1D &BLK, int rc)

	// Input		filename : output file name
	//			BLK	 : 1D-data block
	//			rc       : flag for real versus complex
	// Output		none     : disk file filename is produced
	//				 : in FTNMR ".dat" format

{
  ofstream f_dat;				// File handle
  f_dat.open(filename, ios::out);		// Open file filename
  int n = BLK.elements();			// Number of elements
  int m = 2*sizeof(float)*n+sizeof(int);	// 2*n floating + 1 int
  f_dat.write((char*)&m, sizeof(m));			// Total storage needed
  f_dat.write((char*)&n, sizeof(n));			// Output the integer
  ((row_vector)BLK).write(f_dat);    	// Output BLK (vector typed)
  f_dat.write((char*)&m, sizeof(m));			// Again total needed
  f_dat.close();
  return;
  rc = 0;					// Compiler likes this used
}


void FTNMR_smx(const char *filename, block_2D &BLK, int rc)

	// Input		filename : output file name
	//			BLK	 : 2D-data block
	//			rc       : flag for real versus complex
	// Output		none     : disk file filename is produced
	//				 : in FTNMR ".smx" format

{
  int nrows = BLK.rows();
  int ncols = BLK.cols();
  ofstream fp(filename);		// Open the file for writing
  int crows;				// Write FTNMR smx header
  crows = nrows / 2;
  fp.write((char*)&crows,sizeof(int));
  fp.write((char*)&ncols, sizeof(int));
  crows = 0;
  int i;
  for(i=2; i<4096; i++)
    fp.write((char*)&crows,sizeof(int));

  int j, n, m;
  i = 0;				// Write FTNMR smx data
  while (i < nrows)
    {
    j = 0;
    while (j < ncols)
      {
      m = 0;
      while (m < 64)
        {
        n = 0;
        while (n < 64)
          {
          float f;
          f = Re(BLK(i+m,j+n));
          fp.write((char*)&f, sizeof(float));
          f = Re(BLK(i+m+1,j+n));
          fp.write((char*)&f, sizeof(float));
	  n++;
          }
        m+=2;
        }
      j += 64;
      }
    i += 64;
    }
  fp.close();
  return;
  rc = 0;					// Compiler likes this used
}


block_2D FTNMR_smx(const char *filename)

	// Input		filename : output file name
	// Output		BLK	 : 2D-data block read from 
	// 				   disk file filename which
	//				   exists in FTNMR ".smx" format

{
  int nrows=0, ncols=0;
  ifstream fp(filename);		// Open the file for reading

  int crows=0;				// Read FTNMR smx header
  fp.read((char*)&crows, sizeof(int));
  fp.read((char*)&ncols, sizeof(int));
  fp.seekg(sizeof(int)*4096);

  block_2D BLK(nrows, ncols);		// Read FTNMR smx data
  int i, j, n, m;
  i = 0;
  while (i < nrows)
    {
    j = 0;
    while (j < ncols)
      {
      m = 0;
      while (m < 64)
        {
        n = 0;
        while (n < 64)
          {
          float f;
          fp.read((char*)&f, sizeof(float));
          BLK.put(f,i+m,j+n); 
          fp.read((char*)&f, sizeof(float));
          BLK.put(f,i+m+1,j+n); 
	  n++;
          }
        m+=2;
        }
      j += 64;
      }
    i += 64;
    }
  fp.close();
  return BLK;
}




//block_1D FTNMR_dat(const char *filename)

	// Input		filename : output file name
	// Output		BLK	 : 1D-data block read from 
	// 				   disk file filename which
	//				   exists in FTNMR ".dat" format

//{
//  fstream f_dat;					// File handle
//  f_dat.open(filename, io_readonly);		// Open file filename
//}

#endif /* __FTNMR_CC__ */
