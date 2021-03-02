/* Isotope.cc ***************************************************-*-c++-*-
**									**
** 	                   G A M M A 					**
**									**
**	Isotope 			          Implementation	**
**						 			**
**	Copyright (c) 1990			 			**
**	Tilo Levante, Scott A. Smith		 			**
**	Eidgenoessische Technische Hochschule	 			**
**	Labor fur physikalische Chemie			 		**
**	8092 Zurich / Switzerland			 		**
**							 		**
**      $Header: $
**							 		**
*************************************************************************/

/*************************************************************************
**						 			**
**  Description							 	**
**						 			**
**  Isotope maintains a list of isotopes. It is used as a reference	**
**  for an Isotope. 							**
**									**
*************************************************************************/
 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

#ifndef   Isotope_cc_			// Is the file already included?
#  define Isotope_cc_ 1			// If no, then remember it 
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <GamGen.h>			// Include OS specific stuff
#include <Basics/Isotope.h>		// Include the interface
#include <Basics/IsotopeData.h>		// Include list elements
#include <Basics/StringCut.h>		// Include string cutting
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <Basics/Gconstants.h>		// Include GAMMA constants
#include <Basics/ParamSet.h>		// Include GAMMA parameter sets

#include <string>			// Include libstdc++ string class
#include <vector>			// Include libstdc++ STL vectors
#include <fstream>			// Include libstdc++ I/O streams
#include <cmath>			// Include max and min functions

using std::string;			// Using libstdc++ strings
using std::vector;			// Using libstdc++ vectors

/* Note that although Isotopes is a construct of class Isotope, it must still
   be initialized because it is static. This helps out in library linking    */

vector<IsotopeData> Isotope::Isotopes;		// Set empty list
double              Isotope::RELFRQ1H = 400.13;	// Set rel. Om(1H)

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      CLASS ISOTOPE ERROR HANDLING
// ____________________________________________________________________________
 
/*      Input			I      : Isotope (this)
                               eidx    : Error index
                               noret   : Flag for linefeed (0=linefeed)
                               pname   : Added error message for output
       Output                  none    : Error message output
                                         Program execution stopped if fatal

  The following error messages use the defaults set in the Gutils package

               Case                          Error Message

               (0)                     Program Aborting.....
               (1)                     Problems With Input File Stream
               (2)                     Problems With Output File Stream
               default                 Unknown Error                        */

void Isotope::Isoerror(int eidx, int noret) const
  {
  string hdr("Isotope");
  string msg;
  switch(eidx)
    {
    case 10: msg = string("Spin To Be Added Already Exists!");         // (10)
             GAMMAerror(hdr, msg, noret); break;
    case 33: msg = string("Not Specified From Given Parameters");      // (33)
             GAMMAerror(hdr,msg,noret); break;
    case 34: msg = string("Cannot Properly Glean From Param. Set");    // (34)
             GAMMAerror(hdr,msg,noret); break;
    case 41: msg = string("Can't Read Parameters From File");          // (41)
             GAMMAerror(hdr,msg,noret); break;
    case 42: msg = string("Can't Read Isotope From Parameter File");   // (42)
             GAMMAerror(hdr,msg,noret); break;
    default: GAMMAerror(hdr, eidx, noret); break;	// Unknown Error       (-1)
    }
  }  

void Isotope::Isoerror(int eidx, const string& pname, int noret) const
  {
  string hdr("Isotope");
  string msg;
  switch(eidx)
    {
    case  2: msg=string("Attempted Access of Unknown Isotope ")+pname;	// (2)
      break;
    case  3: msg=string("Cannot Find Parameter ") + pname;		// (3)
      break;
    case 10: msg=string("Cannot Add Another Spin Of Type ") + pname
                + string(" Into GAMMA"); break;				//(10)
    default: msg=string("Unknown Error - ") + pname; break;
    }
  GAMMAerror(hdr, msg, noret);
  }

volatile void Isotope::Isofatal(int eidx) const
  {                                                                 
  Isoerror(eidx, 1);				// Output error message
  if(eidx) Isoerror(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

volatile void Isotope::Isofatal(int eidx, const string& pname) const
  {                                                                 
  Isoerror(eidx, pname, eidx);			// Output error message
  if(eidx) Isoerror(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }
 
// ____________________________________________________________________________
// ii               CLASS ISOTOPE DEALINGS WITH ISOTOPE LISTS
// ____________________________________________________________________________

/* The code for this next function was generated from a program which parsed
   the original GAMMA file Isotopes.list.  Users should keep in mind that
   since our Isotopes.list may be "incomplete" in that we don't have all values
   for all isotopes, we'll use a default value for anything that we don't know.
   Since most values in the list are positive, the default will be negative!
   In those rare cases where an isotope does have some negative values 
   (15N for example) we'll just hope that they won't have the ISODEFVAL I've
   picked.  If they do, just change the number for ISODEFVAL to any negative
   number not in the list and we'll be just fine and still know what an
   "undefined" value is...... 

           Input                I : A dummy isotope (this)
           Output            void : Fills isotopes list with all spins
	   Note			  : The value of RELFRQ1H is set in this
				    routine. It assumes the first entry in
				    the arrays (below) is the proton.
	   Note			  : The value NIso is set herein too!       */

// const double ISODEFVAL = -1.1e6;			// Default value 

void Isotope::set_Isotope_list()
  {
  int NIso = 131;					// Set # isotopes we know

/*          An array of known isotope spin Hilbert spaces (2*I+1)           */

  int Spins[131] = { 2, 3, 2, 2, 1, 1, 3, 4, 5, 4, 7, 4, 1, 2, 3, 2, 1, 
                     6, 1, 2, 1, 4, 1, 4, 6, 6, 2, 2, 4, 4, 4, 4, 4, 8, 
                     8, 6, 8, 13, 8, 4, 6, 2, 8, 4, 4, 4, 6, 4, 4, 10, 4, 
                     2, 4, 4, 6, 4, 10, 2, 6, 10, 6, 6, 10, 6, 6, 2, 6, 2, 
                     2, 2, 2, 10, 10, 2, 2, 6, 8, 2, 2, 6, 2, 4, 8, 4, 4, 
                     11, 8, 6, 8, 8, 8, 8, 6, 6, 4, 4, 4, 6, 6, 8, 8, 2, 
                     2, 6, 8, 15, 8, 10, 8, 2, 6, 6, 2, 4, 4, 4, 2, 4, 2, 
                     4, 2, 2, 2, 10, 10, 4, 6, 4, 4, 8, 2 };

/*          An array of known isotope atomic numbers (1=H, 6=C, ...)        */

  int Numbers[131] = { 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 5, 5, 6, 6, 
                       7, 7, 8, 8, 8, 9, 10, 10, 10, 11, 12, 13, 14, 15, 
                       16, 17, 17, 19, 19, 20, 21, 22, 22, 23, 23, 24, 25, 26, 
                       27, 28, 29, 29, 30, 31, 31, 31, 32, 33, 34, 35, 37, 37, 
                       38, 39, 40, 41, 42, 42, 43, 44, 44, 45, 46, 47, 47, 48, 
                       48, 49, 49, 50, 50, 51, 51, 52, 52, 53, 54, 54, 55, 56, 
                       56, 57, 57, 59, 60, 60, 62, 62, 63, 63, 64, 64, 65, 66, 
                       66, 67, 68, 69, 70, 70, 71, 71, 72, 72, 73, 74, 75, 75, 
                       76, 76, 77, 77, 78, 79, 80, 80, 81, 81, 82, 83, 85, 89, 
                       90, 91, 91, 92, 0 };

/*     An array of known isotope atomic masses  (1=1H, 2=2H, 13=13C, ...)   */

  int Masses[131] = { 1, 2, 3, 3, 4, 6, 6, 7, 8, 9, 10, 
                      11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 
                      21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 
                      43, 45, 47, 49, 50, 51, 53, 55, 57, 59, 61, 
                      63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 85, 
                      87, 87, 89, 91, 93, 95, 97, 99, 99, 101, 103, 
                      105, 107, 109, 111, 113, 113, 115, 117, 119, 121, 123, 
                      123, 125, 127, 129, 131, 133, 135, 137, 138, 139, 141, 
                      143, 145, 147, 149, 151, 153, 155, 157, 159, 161, 163, 
                      165, 167, 169, 171, 173, 175, 176, 177, 179, 181, 183, 
                      185, 187, 187, 189, 191, 193, 195, 197, 199, 201, 203, 
                      205, 207, 209, 211, 227, 229, 231, 233, 235, 0 };

/*   An array of known isotope atomic weights (1H = 1.008665 g/mole, ...)   */

  double Weights[131] = { 1.008665, 2.0140, 3.01605, 3.01603, 4.00260, 
                          6.01888, 6.01512, 7.01600, -1.1e6, 9.01218, 
                          10.0129, 11.00931, 12.00000, 13.00335, 14.00307, 
                          15.00011, 15.99491, -1.1e6, -1.1e6, 18.99840, 
                          19.99244, 20.99395, 21.99138, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, -1.1e6, 
                          0.000001 };

/*        An array of known isotope receptivities (based on 13C = 1)        */

  double Recepts[131] = { 5.68E+3, 8.21E-3, -1.1e6, 3.26E-3, -1.1e6, 
                          -1.1e6, 3.58, 1.54E+3, -1.1e6, 7.88E+1, 
                          2.21E+1, 7.54E+2, -1.1e6, 1.00, 5.69, 
                          2.19E-2, -1.1e6, 6.11E-2, -1.1e6, 4.73E+3, 
                          -1.1e6, 3.59E-2, -1.1e6, 5.25E+2, 1.54, 
                          1.17E+3, 2.09, 3.77E+2, 9.73E-2, 2.02E+1, 
                          3.77, 2.69, 3.28E-2, 5.27E-2, 1.71E+3, 
                          8.64E-1, 1.18, 7.55E-1, 2.16E+3, 4.90E-1, 
                          9.94E+2, 4.19E-3, 1.57E+3, 2.41E-1, 3.65E+2, 
                          2.01E+2, 6.65E-1, 2.37E+2, 3.19E+2, 6.17E-1, 
                          1.43E+2, 2.98, 2.26E+2, 2.77E+2, 4.34E+1, 
                          2.77E+2, 1.07, 6.68E-1, 6.04, 2.74E+3, 
                          2.88, 1.84, 2.13E+3, 8.30E-1, 1.56, 
                          1.77E-1, 1.41, 1.95E-1, 2.76E-1, 6.93, 
                          7.6, 8.38E+1, 1.89E+3, 1.95E+1, 2.52E+1, 
                          5.20E+2, 1.11E+2, 8.90E-1, 1.25E+1, 5.30E+2, 
                          3.18E+1, 3.31, 2.69E+2, 1.83, 4.41, 
                          4.30E-1, 3.36E+2, 1.65E+3, 2.34, 3.70E-1, 
                          1.26, 0.59, 4.83E+2, 4.53E+1, 2.33E-1, 
                          4.84E-1, 3.31E+2, 4.47E-1, 1.58, 1.02E+3, 
                          6.59E-1, 3.21, 4.44, 1.22, 1.72E+2, 
                          5.47, 6.70E-1, 1.69E-1, 2.04E+2, 5.89E-2, 
                          2.80E+2, 4.90E+2, 1.14E-3, 2.13, 5.36E-2, 
                          1.16E-1, 1.91E+1, 1.43E-1, 5.42, 1.08, 
                          2.89E+2, 7.69E+2, 1.18E+1, 7.77E+2, -1.1e6, 
                          -1.1e6, -1.1e6, -1.1e6, -1.1e6, 4.95E-3, 
                          -1.1e6 };

/*  An array of known isotope relative frequencies (from 1H @ 400.130 MHz)  */

  double RelFreqs[131] = { 400.130, 61.424, 426.791, 304.811, -1.1e6, -1.1e6, 
                           58.883, 155.503, -1.1e6, 56.226, 42.990, 128.330, 
                           -1.1e6, 100.613, 28.913, -40.561, -1.1e6, 54.242, 
                           -1.1e6, 376.498, -1.1e6, 31.586, -1.1e6, 105.842, 
                           24.496, 104.262, 79.494, 161.977, 30.714, 39.205, 
                           32.635, 18.670, 10.247, 26.925, 97.200, 22.563, 
                           22.557, 39.893, 105.246, 22.615, 99.012, 12.956, 
                           94.939, 35.756, 106.146, 113.705, 25.036, 96.035, 
                           122.028, 13.957, 68.514, 76.313, 100.249, 108.063, 
                           38.633, 130.927, 17.342, 19.606, 37.196, 97.936, 
                           26.076, 26.625, 90.061, 18.462, 20.691, 12.744, 
                           18.310, 16.197, 18.622, 84.832, 88.741, 87.492, 
                           87.680, 142.574, 149.212, 95.755, 51.853, 104.714, 
                           126.241, 80.062, 110.668, 32.811, 52.481, 39.749, 
                           44.466, 52.793, 56.522, 117.202, 21.755, 13.384, 
                           16.517, 13.160, 99.236, 43.818, 15.281, 19.102, 
                           90.741, 13.180, 18.338, 82.079, 11.564, 33.095, 
                           70.475, 19.414, 45.643, 31.722, 12.484, 7.478, 
                           47.976, 16.669, 90.129, 91.038, 9.131, 31.070, 
                           6.874, 7.486, 85.876, 6.918, 71.667, 26.457, 
                           228.970, 231.219, 83.459, 64.297, -1.1e6, -1.1e6, 
                           -1.1e6, -1.1e6, -1.1e6, 7.162, -263376.60322 };

/*                       An array of known isotope names                     */

  string Names[87] = { "Hydrogen",   "Deuterium",  "Tritium",  "Helium", 
                            "Lithium",    "Beryllium",  "Boron",    "Carbon", 
                            "Nitrogen",   "Oxygen",     "Fluorine", "Neon", 
                            "Sodium",     "Magnesium",  "Aluminum", "Silicon", 
                            "Phosphorus", "Sulfur",     "Chlorine", "Potassium", 
                            "Calcium",    "Scandium",   "Titanium", "Vanadium", 
                            "Chromium",   "Manganese",  "Iron",     "Cobalt", 
                            "Nickel",     "Copper",     "Zinc",     "Gallium", 
                            "Germanium",  "Arsenic",    "Selenium", "Bromine", 
                            "Rubidium",   "Strontium",  "Yttrium",  "Zirconium", 
                            "Niobium",    "Molybdenum", "Technetium", "Ruthenium", 
                            "Rhodium",    "Palladium",  "Silver",   "Cadmium", 
                            "Indium",     "Tin",        "Antimony", "Tellurium", 
                            "Iodine",     "Xenon",      "Cesium",   "Barium", 
                            "Lanthanum",  "Praseodymium", "Neodymium", "Samarium", 
                            "Europium",   "Gadolinium", "Terbium",  "Dysprosium", 
                            "Holmium",    "Erbium",     "Thulium",  "Ytterbium", 
                            "Lutetium",   "Hafnium",    "Tantalum", "Tungsten", 
                            "Rhenium",    "Osmium",     "Iridium",  "Platinum", 
                            "Gold",       "Mercury",    "Thallium", "Lead", 
                            "Bismuth",    "Astatine",   "Actinium", "Thorium", 
                            "Protactinium", "Uranium",  "Electron" };

/*                       Known Isotope Indices For Previous Names Array   */

  int NamesIdx[131] = { 0, 1, 2, 3, 3, 3, 4, 4, 4, 5, 6, 
                        6, 7, 7, 8, 8, 9, 9, 9, 10, 11, 11, 
                        11, 12, 13, 14, 15, 16, 17, 18, 18, 19, 19, 
                        20, 21, 22, 22, 23, 23, 24, 25, 26, 27, 28, 
                        29, 29, 30, 31, 31, 32, 33, 34, 35, 35, 36, 
                        36, 37, 38, 39, 40, 41, 41, 42, 43, 43, 44, 
                        45, 46, 46, 47, 47, 48, 48, 49, 49, 50, 50, 
                        51, 51, 52, 53, 53, 54, 55, 55, 56, 56, 57, 
                        58, 58, 59, 59, 60, 60, 61, 61, 62, 63, 63, 
                        64, 65, 66, 67, 67, 68, 68, 69, 69, 70, 71, 
                        72, 72, 73, 73, 74, 74, 75, 76, 77, 77, 78, 
                        78, 79, 80, 81, 82, 83, 84, 84, 85, 86 };

  string Elements[85] = { "H",  "He", "Li", "Be",  "B",  "C",  "N",  "O", 
                               "F",  "Ne", "Na", "Mg", "Al", "Si",  "P",  "S", 
                               "Cl",  "K", "Ca", "Sc", "Ti",  "V", "Cr", "Mn", 
                               "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", 
                               "Se", "Br", "Rb", "Sr",  "Y", "Zr", "Nb", "Mo", 
                               "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", 
                               "Sb", "Te",  "I", "Xe", "Cs", "Ba", "La", "Pr", 
                               "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", 
                               "Tm", "Yb", "Lu", "Hf", "Ta",  "W", "Re", "Os", 
                               "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "At", 
                               "Ac", "Th", "Pa",  "U",  "e" };

  int IElements[131] = {  0,  0,  0,  1,  1,  1,  2,  2,  2,  3,  4, 
                          4,  5,  5,  6,  6,  7,  7,  7,  8,  9,  9, 
                         9, 10, 11, 12, 13, 14, 15, 16, 16, 17, 17, 
                         18, 19, 20, 20, 21, 21, 22, 23, 24, 25, 26, 
                         27, 27, 28, 29, 29, 30, 31, 32, 33, 33, 34, 
                         34, 35, 36, 37, 38, 39, 39, 40, 41, 41, 42, 
                         43, 44, 44, 45, 45, 46, 46, 47, 47, 48, 48, 
                         49, 49, 50, 51, 51, 52, 53, 53, 54, 54, 55, 
                         56, 56, 57, 57, 58, 58, 59, 59, 60, 61, 61, 
                         62, 63, 64, 65, 65, 66, 66, 67, 67, 68, 69, 
                         70, 70, 71, 71, 72, 72, 73, 74, 75, 75, 76, 
                         76, 77, 78, 79, 80, 81, 82, 82, 83, 84 };

  IsotopeData ID;				// Working isotope data
  for(int i=0; i<NIso; i++)			// Loop over all isotopes
    {
    ID._HS          = Spins[i]; 		//   Set spin Hilbert space
    ID._number      = Numbers[i]; 		//   Set the atomic number
    ID._mass        = Masses[i]; 		//   Set the atomic mass
    ID._weight      = Weights[i]; 		//   Set the atomic weights
    ID._receptivity = Recepts[i]; 		//   Set the receptivities
    ID._relfreq     = RelFreqs[i]; 		//   Set the relative frequency
    ID._name        = Names[NamesIdx[i]];	//   Set the isotope name 
    ID._element     = Elements[IElements[i]];	//   Set the isotope element 
    ID._symbol      = Gdec(ID._mass)+ID._element;
    ID._iselectron  = false;                    //   Set the type to nucleus
    if(ID._element == string("e"))		//   Exceptions occur for all
      {                                         //   electrons in this scheme
      ID._symbol = string("e-");		//     Use symbol e- (NOT 0e)
      ID._iselectron = true;                    //     Flag it is electron
      }
    Isotopes.push_back(ID);			//   Set isotope list entry
    }

  RELFRQ1H = RelFreqs[0];			// Setting Base 1H Frequency
  }						// [not using SetRel1HF()]




        // Input                I : A dummy isotope (this)
        // Output            void : Sets the relative 1H frequency
	// Note			  : Must have filled isotopes list first!

/* GAMMA's isotopes, rather than having a known gyromagnetic ratio, have
   a relative Larmor frequency from which gyromagnetic ratio may be
   calculated.  The base relative Larmor frequency used in GAMMA is that
   of a proton, RELFRQ1H, and to date will be 400.13 (MHz) in GAMMA. 
   The gyromagnetic ratio of isotope type i relates to its relative
   frequency RelFreq according to

        gamma  = RelFreq*GAMMA1H/RELFRQ1H = Omega * gamma  / Omega
             i                                   i       1H       1H

   GAMMA has the DEFAULT value GAMMA1H = 2.67519e8 (1/T-sec) set in its
   constants module (Gutils).  It sets the default proton relative frequency
   herein because that value (usually 400.13) is imported into GAMMA from a
   stardard isotopes list.  Were that list to change, perhaps the listed base
   frequencys would change. But that is O.K. because it will be set here
   accordingly.

   Note however that this does NOT affect any GAMMA programs that wish to set
   alternate field strengths through spin systems.  The relative frequencies
   in GAMMA's isotope list are Larmor frequencies but used ONLY to determine
   gamma values. Larmor frequencies in MR simulations are set independently
   by the user in typical GAMMA programs.                                    */

void Isotope::SetRel1HF()
  {
  Isotope I("1H");				// Use proton for base
  RELFRQ1H = (Isotopes[I.Iso]).rel_freq();	// Here is rel.frq. (400.13)
  }


// ____________________________________________________________________________
// iii              Class Isotope Private Parameter Set Functions
// ____________________________________________________________________________

/* This is used for reading an Isotope in from a GAMMA Parameter Set         */

bool Isotope::SetIsotope(const ParameterSet& pset, int idx, bool warn)
  {
  string pname("Iso");				// Isotope parameter name
  if(idx != -1)					// Add a suffix of (#) where
    pname += string("(")+Gdec(idx)+string(")");	// # = idx if idx is not -1
  ParameterSet::const_iterator item;		// A pix into parameter list
  item = pset.seek(pname);			// Look in parameter set
  if(item != pset.end())			// If the parameter is found
    {
    string sval, pstate;			// Parameter value & statement
    (*item).parse(pname,sval,pstate);		// Parse out value as string
    *this = Isotope(sval);			// Set out isotope type
    return true;				// Return we were successful
    }
  if(warn) Isoerror(3, pname, 1);
  return false;					// Return we failed our quest
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A                     ISOTOPE CONSTRUCTION, ASSIGNMENT
// ____________________________________________________________________________

/* These construct single spin isotopes.  Note that on the very first isotope
   made in a program a persistent list of spin isotopes is automatically
   generated.  This remains in memory until the program ends. Also note that
   "construction" doesn't do a whole lot except set the isotope to point
   to a member in the isotope list.                                          */

Isotope::Isotope()
  {
  IsotopeData Proton("1H");			// Set mock default isotope,
  if(!Isotopes.size()) set_Isotope_list();	// Fill Isotopes list if empty
  Iso = seek(Proton);				// Get 1H index in Isotopes list
  if(Iso < 0) Isofatal(2, "1H");	 	// If no 1H FATAL ERROR!
  }

Isotope::Isotope(const Isotope& I)
        :Iso(I.Iso) {}				// Only copy index. Cannot have
						// any isotope without list set!
Isotope::Isotope(const string& symbol)
  {
  IsotopeData i(symbol);			// Form data for this isotope 
  if(!Isotopes.size()) set_Isotope_list();	// Fill Isotopes list if empty
  Iso = seek(i);				// Get index into Isotopes list
  if(Iso < 0) Isofatal(2, symbol); 		// Fatal error, unknown isotope
  }

Isotope& Isotope::operator= (const Isotope& I)
{ if(this == &I) return *this; Iso = I.Iso; return *this; }

Isotope::~Isotope() { }

// ____________________________________________________________________________
// C                        ISOTOPE ACCESS FUNCTIONS
// ____________________________________________________________________________
 
        // Input                I : An Isotope (this)
        // Note                   : All functions use IsotopesData
	//			    except gamma which must be calculated
	//			    based on relative frequencies

/*
      Function        Return      Units                Example(s)
      --------     -------------  -----   -------------------------------
         qn        mz   (double)   hbar   0.5, 1.0, 1.5, ....
         HS        2I+1 (int)      none   2, 3, 4, .....
      momentum     mz   (string)   none   1/2, 1, 3/2, .....
       symbol           (string)   none   1H, 2H, 13C, 19F, ....
       name             (string)   none   Hydrogen, Lithium, Carbon, ...
      element           (string)   none   H, Li, F, C, ......
       number           (int)      none   1<-H, 3<-Li, 6<-C, .....
       mass             (int)      amu    1<-1H, 2<-2H, 13<-13C, ....
      weight            (double)   g/m    1.00866<-1H, 7.016<-7Li, ....
       gamma            (double)   1/T-s  2.67519*10^8
    receptivity         (double)   none   400.13, 155.503, ...
 relative_frequency     (double)   MHz    400.13, 155.503, ...                */

      double       Isotope::qn()          const { return (Isotopes[Iso]).qn(); }
      int          Isotope::HS()          const { return (Isotopes[Iso]).HS(); }
      string  Isotope::momentum()         const { return (Isotopes[Iso]).momentum();}
const string& Isotope::symbol()           const { return (Isotopes[Iso]).symbol(); }
const string& Isotope::name()             const { return (Isotopes[Iso]).name(); }
const string& Isotope::element()          const { return (Isotopes[Iso]).element(); }
      int          Isotope::number()      const { return (Isotopes[Iso]).number(); }
      int          Isotope::mass()        const { return (Isotopes[Iso]).mass(); }
      double       Isotope::weight()      const { return (Isotopes[Iso]).weight(); }
      double       Isotope::receptivity() const { return (Isotopes[Iso]).recept(); }
      bool         Isotope::electron()    const { return (Isotopes[Iso]).electron(); }

double Isotope::relative_frequency() const {return (Isotopes[Iso]).rel_freq();}
double Isotope::gamma() const
                      { return ((Isotopes[Iso]).rel_freq()/RELFRQ1H)*GAMMA1H; }

// ____________________________________________________________________________
// D                        ISOTOPE I/O FUNCTIONS
// ____________________________________________________________________________

/*  These read/write the isotope to from an ASCII file via Single Parameters */

bool Isotope::read(const string& filename, int idx, int warn)
   {
   ParameterSet pset;                   // Declare a parameter set
   if(!pset.read(filename, warn?1:0))	// Read in pset from file
     {
     Isoerror(1, filename,1);		// Filename problems
     if(warn > 1) Isofatal(41);		// Cant read from file, fatal
     else         Isoerror(41);		// or non-fatal warning
     return false;
     }
   if(!read(pset, idx, warn?1:0))	// User overloaded function
     {
     Isoerror(1, filename,1);		// Filename problems
     if(warn > 1) Isofatal(42);		// Cannot read isotope from file
     else         Isoerror(42);		// or non-fatal one
     return false;
     }
   return true;
   }

bool Isotope::read(const ParameterSet& pset, int idx, int warn)
  {
  bool TF=SetIsotope(pset,idx,warn?true:false);	// Try & set ourself up
  if(!TF)                                       // If we didn't handle
    {                                           // setting ourself properly
    if(warn)                                    // then we'll issue some
      {                                         // warnings if desired
                   Isoerror(33, 1);		//   Cant read parameters
      if(warn > 1) Isofatal(34);		//   Can't set from parameters
      else         Isoerror(34,1);		//   or a warning issued
      }
    return false;                               // Return that we failed
    }
  return TF;
  }

vector<string> Isotope::printStrings(bool hdr) const
  {
  vector<string> PStrings;
if(!hdr) hdr = true;				// Why is hdr an argument?
if(!hdr) return PStrings;			// Use just to stop Borland warnings
  IsotopeData ID = Isotopes[Iso];
  PStrings = ID.printStrings(true);
  string gs = string(" Gamma    ")
            + Gform("%10.4f", gamma() * 1.e-8) 
             + string(" 10^-8/T-s\n");
  PStrings.push_back(gs);
  return PStrings;
  }

std::ostream& Isotope::print(std::ostream& ostr) const
  {
  vector<string> PStrings = printStrings(true);
  unsigned ns = PStrings.size();		// No. of strings to print
  unsigned i;					// Looping index
  unsigned ml = 0;				// Maximum string length
  for(i=0; i<ns; i++)				// Loop over strings and find
    ml = gmax(PStrings[i].length(), ml);	// the longest one
  int cl = 40 - ml/2;
  string cstr;
  if(cl > 0) cstr = string(cl, ' ');
  for(i=0; i<ns; i++)
    ostr << cstr << PStrings[i] << std::endl;
  return ostr; 
  }

std::ostream& operator<< (std::ostream& ostr, const Isotope& I)
  { return I.print(ostr); }
 
// ____________________________________________________________________________
// E                          ISOTOPE LIST FUNCTIONS
// ____________________________________________________________________________

/* Allow users to search for a specific isotope in the list or find whether
   a particular isotope is contained in the list.  Remember that each Isotope
   only contains a pointer to the Isotopes list with an index for which
   Isotope in the list it actually is.

           Input                I : An isotope (this)
                               ID : A single isotope
                           symbol : A string designating an isotope
           Output               i : Index of isotope in list if it exists.
                                    If id does NOT exist, returns is -1
                               TF : TRUE if isotope of type "symbol" known
           Note                   : Will return FALSE if symbol not found
                                    due to no isotope list                    */
int Isotope::seek(const IsotopeData& I)
  {
  if(!Isotopes.size()) return -1;		// If Isotopes empty, fail
  int NISO = Isotopes.size();			// Get size of list, some
  for(int i=0; i<NISO; i++)			// may have been added!
    if(I.symbol() == (Isotopes[i]).symbol()) return i;
  return -1;
  }

bool Isotope::exists(const string& symbol)
  {
  IsotopeData i(symbol);			// Form data for this isotope 
  if(!Isotopes.size()) set_Isotope_list(); 	// Fill Isotopes list if empty
  return (seek(i)<0)?false:true;
  }

bool Isotope::known(const string& symbol)
  { Isotope X; return X.exists(symbol); }

int Isotope::size() { return int(Isotopes.size()); }

vector<string> Isotope::PrintListStrings()
  {
  vector<string> PStrings;

  Isotope X("1H");			// A temp dummy isotope
  int NISO = X.size();			// Get size of list, some

//                     Set Up Output Column Alignment

  int nl = (Gdec(NISO)).length();	// Length of printed index
  string nspc(nl+1, ' ');		// Spacer for index column

  string csp("  ");			// Space between each columns
  int cspl = csp.length();		// Length of space between cols

  string sy;				// String for symbol
  int syl, syi, syf;			// For symbol alignment

  string icol("     ");			// Space between isotopes

  int nml = 0;				// For Name column width
  int nmi, nmf;				// For Name column spacing
  int i;				// Dummy index
  for(i=0; i<NISO; i++)			// Find longest name
    nml = gmax(nml, int(((Isotopes[i]).name()).length()));
  string nms("Name");			// For Name column header
  nmi = (nml-4)/2;
  nmf = nml - nmi - 4;
  if(nmi > 0) nms = string(nmi, ' ') + nms;
  if(nmf > 0) nms += string(nmf, ' ');

  int ncols = 2;			// No. isotopes per row we print
  int next  = 0;			// For line end
  int cw  = nl + 1       + cspl;
  cw += 6            + cspl;
  cw += nms.length() + cspl;
  cw *= ncols;
  cw += (ncols-1)*icol.length();
  string cen(40-cw/2, ' ');

//		       Set Up Column Headers

  string pline = cen;			// Initial spacing (center)
  for(i=0; i<ncols; i++)
    {
    pline += "Indx"   + csp;		// Index
    pline += "Symbol" + csp;		// Symbol
    pline += nms      + csp;		// Name
    pline += icol;			// Next Isotope
    }
  PStrings.push_back(pline);

  pline = cen;				// Initial spacing (center)
  for(i=0; i<ncols; i++)
    {
    pline += string(nl+1, '=')         + string(cspl, ' ');
    pline += string(6, '=')            + string(cspl, ' ');
    pline += string(nms.length(), '=') + string(cspl, ' ');
    pline += icol;
    }
  PStrings.push_back(pline);

  pline = string("");
  for(i=0; i<NISO; i++)
    {
    sy = (Isotopes[i]).symbol();		// Set symbol string
    syl = sy.length();
    syi = (6-syl)/2;
    syf = 6-syi-syl;
    if(syi>0) sy  = string(syi, ' ') + sy;
    if(syf>0) sy += string(syf, ' ');
    
    nms = (Isotopes[i]).name();
    nms += string(nml - nms.length(), ' ');
    if(!next) pline = cen;
    pline += Gdec(i+1, nl) + "." + csp;		// Print index
    pline += sy                  + csp;		// Print symbol
    pline += nms                 + csp;
    next++;
    if(next >= ncols) { next=0; PStrings.push_back(pline); }
    else              pline += icol;
    }
  if(next) PStrings.push_back(pline);
  return PStrings;
  }


void Isotope::PrintList(std::ostream& ostr, bool hdr)
  {
  if(hdr)
    {
    string H("Currently Available Isotopes");
    ostr << string(40-H.length()/2, ' ') << H << "\n";
    }
  Isotope X("1H");
  int NISO = X.size();			// Get size of list, some

//                     Set Up Output Column Alignment

  int nl = 1;				// Length of printed index
  if(NISO >= 9)    nl++;
  if(NISO >= 99)   nl++;
  if(NISO >= 999)  nl++;
  if(NISO >= 9999) nl++;
  string nspc(nl+1, ' ');		// Spacer for index

  string sy;				// String for symbol
  int syl, syi, syf;			// For symbol alignment

  string csp("  ");			// Space between each columns
  int cspl = csp.length();		// Length of space

  string icol("     ");			// Space between isotopes

  int nml = 0;				// For Name column width
  int nmi, nmf;				// For Name column spacing
  int i;				// Dummy index
  for(i=0; i<NISO; i++)			// Find longest name
    nml = gmax(nml, int(((Isotopes[i]).name()).length()));
  string nms("Name");			// For Name column header
  nmi = (nml-4)/2;
  nmf = nml - nmi - 4;
  if(nmi > 0) nms = string(nmi, ' ') + nms;
  if(nmf > 0) nms += string(nmf, ' ');

  int ncols = 2;			// Number of columns
  int next  = 0;			// For line end
  int cw  = nl + 1       + cspl;
  cw += 6            + cspl;
  cw += nms.length() + cspl;
  cw *= ncols;
  cw += (ncols-1)*icol.length();
  string cen(40-cw/2, ' ');

//			Write Column Headers

  ostr << "\n" << cen;
  for(i=0; i<ncols; i++)
    {
    ostr << "Indx"   << csp;		// Index
    ostr << "Symbol" << csp;		// Symbol
    ostr << nms      << csp;		// Name
    ostr << icol;			// Next Isotope
    }

  ostr << "\n" << cen;
  for(i=0; i<ncols; i++)
    {
    ostr << string(nl+1, '=')         << string(cspl, ' ');
    ostr << string(6, '=')            << string(cspl, ' ');
    ostr << string(nms.length(), '=') << string(cspl, ' ');
    ostr << icol;
    }

  for(i=0; i<NISO; i++)
    {
    sy = (Isotopes[i]).symbol();			// Set symbol string
    syl = sy.length();
    syi = (6-syl)/2;
    syf = 6-syi-syl;
    if(syi>0) sy = string(syi, ' ') + sy;
    if(syf>0) sy += string(syf, ' ');
    
    nms = (Isotopes[i]).name();
    nms += string(nml - nms.length(), ' ');
    if(!next) ostr << "\n" << cen;
    ostr << Gdec(i+1, nl) << "." << csp;		// Print index
    ostr << sy                   << csp;		// Print symbol
    ostr << nms                  << csp;
    next++;
    if(next >= ncols) next = 0;
    else              ostr << icol;
    }
  }


// ____________________________________________________________________________
// F                   Isotope List Expansion Functions
// ____________________________________________________________________________

/* These functions allow users to add additional spins into the working GAMMA
   isotopes list.  This will primarily be used in EPR/ESR/EMR simulations when
   electrons with I>1/2 are required. For example, this will occur when dealing
   with an electron about a metal center when the surroundings effecitvely 
   cause it to be split. To add in an additional spin one simply constructs a
   new isotope type then adds/subtracts it with these functions.  Note that one
   may NOT remove any "standard" isotopes in GAMMA nor may one add an isotope
   that already exists..						     */

bool Isotope::AddIsotope(const IsotopeData& ID, int warn)
  {
  Isotope X;					// A dummy isotope
  if(!Isotopes.size()) X.set_Isotope_list(); 	// Fill Isotopes list if empty
  if(X.seek(ID) >= 0)				// If isotope already in list
    {						// we cannot add it again!
    if(warn)					//   If warnings set we output
      {						//   error messages, maybe quit
      X.Isoerror(10, 1);			//     Spin already exists
      string pname = ID.symbol();		//     Get isotope name
      if(warn>1) X.Isofatal(10,ID.symbol());	//     We cannot add another
      else       X.Isoerror(10,   ID.symbol());	//     Quit program if needed
      }
    return false;				//   Return we failed to add it
    }
  Isotopes.push_back(ID);			// Add ID to Isotopes List
  return true;					// Return added successfully
  }

// ____________________________________________________________________________
// G        Isotope Container Support Functions & (In)Equality
// ____________________________________________________________________________

/* These functions allow for comparisons between spin isotopes. Such functions
   are required under some compilers if one desired to build containers of
   spin isotopes (STL lists, vectors, ....).  Equal spin isotopes will simply
   point to the same entry in the isotopes list.  For sorting purposes we
   go by the spin Hilbert space associated with the isotope.		     */

bool Isotope::operator== (const Isotope& i) const { return (Iso==i.Iso); }
bool Isotope::operator!= (const Isotope& i) const { return (Iso!=i.Iso); }
bool Isotope::operator<  (const Isotope& I) const { return (HS() < I.HS()); } 
bool Isotope::operator>  (const Isotope& I) const { return (HS() > I.HS()); }

// ____________________________________________________________________________
// H                    Isotope Auxiliary Functions
// ____________________________________________________________________________

/* These are just handy to have available for other parts of GAMMA           */

bool Isotope::nepair(const Isotope& S) const { return enpair(S); }

bool Isotope::enpair(const Isotope& S) const
  {
  if(electron()) return S.electron()?false:true;
  else           return S.electron()?true:false;
  }

bool Isotope::eepair(const Isotope& S) const
  { if(electron())  return S.electron()?true:false; return false; }

bool Isotope::nnpair(const Isotope& S) const
{ if(!electron()) return S.electron()?false:true; return false; }


#endif						// Isotope.cc
