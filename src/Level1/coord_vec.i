/* coord_vec.h */
// Swig interface file.

%{
#include "Level1/coord_vec.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

%rename(__add__)  coord_vec::operator+ const;
%rename(__iadd__) coord_vec::operator+=;
%rename(__sub__)  coord_vec::operator- const;
%rename(__isub__) coord_vec::operator-=;

%rename(__mul__)  coord_vec::operator* const;
%rename(__imul__) coord_vec::operator*=;
%rename(__div__)  coord_vec::operator/ const;
%rename(__idiv__) coord_vec::operator/=;

%rename(__assign__) coord_vec::operator=;



class coord_vec
{

coord *Pts;			// Pointer to Array of Coordinates
int    Npts;		// Number of Coordinates


private:

void CVerror(int eidx, int noret=0) const;
void CVerror(int eidx, const std::string& pname, int noret=0) const;
volatile void CVfatality(int eidx) const;

int SetNPoints(const ParameterSet& pset, int warn=2);
int SetCoords(const  ParameterSet& pset, int warn=-1);
 
void check(int index) const;


public:

coord_vec( );
coord_vec(int pts);
coord_vec (const coord_vec& cvec1);
coord_vec(const ParameterSet& pset, int idx=-1, int warn=2);
coord_vec(const row_vector& X, const row_vector& Y, const row_vector& Z);

virtual ~coord_vec ();

coord_vec& operator= (const coord_vec& cvec1);

coord_vec xrotate(double theta, int rad=0) const;
coord_vec yrotate(double theta, int rad=0) const;
coord_vec zrotate(double phi,   int rad=0) const;

coord_vec rotate(const matrix& Rmx) const;
coord_vec rotate(double alpha, double beta, double gamma) const;
coord_vec rotate(const coord& EA) const;

void rotate_ip(double alpha, double beta, double gamma);
void rotate_ip(const coord &EA);

coord_vec translate(double    delx, double dely=0, double delz=0) const;
void      translate_ip(double delx, double dely=0, double delz=0);
coord_vec translate(const     coord& delpt) const;
void      translate_ip(const  coord& delpt);

coord_vec trans_x(double delx) const;
coord_vec trans_y(double dely) const;
coord_vec trans_z(double delz) const;

void trans_x_ip(double delx);
void trans_y_ip(double dely);
void trans_z_ip(double delz);

row_vector project(int projx, int projy) const;
      
//friend    coord_vec  operator *  (double r, const coord_vec& cvec1);

coord_vec  operator *  (double r) const;
coord_vec& operator *= (double r);
coord_vec  operator /  (double r) const;
coord_vec& operator /= (double r);

//coord_vec operator +  (const coord& pt) const;
//void      operator += (const coord& pt);
//coord_vec operator -  (const coord& pt) const;
//void      operator -= (const coord& pt);

coord_vec  operator +  (const coord_vec& cv) const;
coord_vec& operator += (const coord_vec& cv);
coord_vec  operator -  (const coord_vec& cv) const;
coord_vec& operator -= (const coord_vec& cv);

int size() const;

double max_x() const;
double max_y() const;
double max_z() const;

coord maxima() const;

void maxima(double &x, double &y, double &z) const;

double max_R() const;

void max_R(int &maxi, double &maxR) const;

coord_vec vectors() const;

coord_vec vectors_f() const;

double distance(int pt1, int pt2, int Angs=0) const;

matrix distances(int Angs=0) const;
matrix thetas(int    deg=0)  const;
matrix phis(int      deg=0)  const;


// Tried "renaming this like this (at the top of this file):
//     %rename(__funcall__)  coord_vec::operator();
// But it did not work, but worked without this "renamed".

coord     operator() (int index)                    const;

void      put(const coord &pt1, int index);
void      put(double x, double y, double z, int index);
coord     get(int index)                            const;
double    x(int index)                              const;
double    y(int index)                              const;
double    z(int index)                              const;
coord_vec get_block(int index, int npts)            const;
void      put_block(int index, const coord_vec& cv) const;

//operator ParameterSet( ) const;
//friend void operator+= (ParameterSet& pset, const coord_vec& cvec);

void PSetAdd(ParameterSet& pset, int idx=-1) const;

//int operator= (const ParameterSet& pset);
//int write(const std::string &filename, int idx=-1, int warn=2) const;
//int write(std::ofstream& ofstr,        int idx=-1, int warn=2) const;

bool is_zero( ) const;
 
virtual std::string ask_read(int argc, char* argv[], int argn,
                                         int idx=-1, int warn=2);

virtual int read(const std::string& filename, int idx=-1, int warn=2);
virtual int read(const ParameterSet& pset,    int idx=-1, int warn=2);

//std::ostream& printf(std::ostream& out, int units = 0) const;
//std::ostream& printCylindrical(std::ostream& ostr, double sf=1) const;
//std::ostream& printSpherical(std::ostream& ostr, double sf=1) const;
//std::ostream& printCartesian(std::ostream& ostr, double sf=1) const;
//       std::ostream& print(std::ostream& ostr, int ncols=2, int N=-1) const;
//friend std::ostream& operator << (std::ostream& ostr, const coord_vec &cvec);
 
coord_vec Cart2Sph(int rad=1) const;
coord_vec Sph2Cart(int rad=1) const;
 
coord_vec Cart2Cyl(int rad=1) const;
coord_vec Cyl2Cart(int rad=1) const;

coord_vec Sph2Cyl(int rad=1) const;
coord_vec Cyl2Sph(int rad=1) const;

};

