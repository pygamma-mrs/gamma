/* Gconstants.i */
   
#ifndef PI				// Define PI once
#define PI 3.14159265358979323846
#endif
 
#ifndef PI2                             // Define 2PI once so we'll
#define PI2 6.28318530717958647692      // always have it around
#endif
 
extern const double PIx2;	// 2*pi

extern const double DEG2RAD;            // Degrees -> radians
extern const double RAD2DEG;            // Radians -> degrees

extern const double HZ2RAD;             // Cycles/sec -> rad/sec
extern const double RAD2HZ;             // Rad/sec -> cycles/sec

extern const double HZ2GAUSS;		// cycles/sec -> Gauss
extern const double GAUSS2HZ;		// Gauss -> cycles/sec
extern const double GHZ2GAUSS;		// h/beta*1.e9
extern const double GAUSS2GHZ;		// Gauss -> cycles/sec * 1.e-9

extern const double      MU_E; 		// Electron magnetic moment  (J/T)
extern const double      BOHRMAG;	// Bohr magneton             (J/T)
extern const double      GFREE;		// Free electron g factor
extern const double      GAMMAe;	// Free e gyromagnetic ratio (1/T-sec)
extern const double      GAMMA1H;	// Proton gyromagnetic ratio (1/T-sec)	
extern const std::string DEFISO;	// Default spin isotope symbol

extern const double PLANCK;		// Plancks constant (h)         (J-sec)
extern const double HBAR;   		// Plancks constant (h/2PI)     (J-sec)
