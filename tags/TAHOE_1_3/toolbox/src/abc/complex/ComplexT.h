/* $Id: ComplexT.h,v 1.8 2002-09-12 16:35:36 paklein Exp $ */
/* created: PAK/AFLP (05/19/1997) */
#ifndef _COMPLEX_T_H_
#define _COMPLEX_T_H_

/* environment */
#include "Environment.h"

/* forward declarations */
#include "ios_fwd_decl.h"

namespace Tahoe {

template <class nTYPE> class nArrayT;

/** complex number class */
class ComplexT
{
public:

	/** \name constructors */
	/*@{*/
	ComplexT(void);  				/**< default needed to making arrays     */
	ComplexT(int i);                /**< type conversion: int -> ComplexT    */
	ComplexT(double re);			/**< type conversion: double -> ComplexT */
	ComplexT(double re, double im);
	/*@}*/

	/** \name Real and Imaginary parts */
	/*@{*/
	double Re(void) const;
	double Im(void) const;
	ComplexT& toZ(double re, double im); /**< returns reference to this */
	/*@}*/
	
	/** \name Conjugate */
	/*@{*/
	ComplexT& Conjugate(const ComplexT& z);
	ComplexT& Conjugate(void) ;
	/*@}*/
	
	/** \name Real and Imaginary parts of arrays
	 * must be dimensioned BEFORE call */
	/*@{*/
	static void z_to_Re(const nArrayT<ComplexT>& z, nArrayT<double>& d);
	static void z_to_Im(const nArrayT<ComplexT>& z, nArrayT<double>& d);
	static void ReIm_to_z(const nArrayT<double>& re, const nArrayT<double>& im, nArrayT<ComplexT>& z);
	/*@}*/
	
	/** \name Polar components */
	/*@{*/
	double Magnitude() const;
	double Angle() const;
	/*@}*/

	/** \name Assignment operators */
	/*@{*/
	ComplexT& operator=(const ComplexT& zRHS);
	ComplexT& operator=(double re);
	ComplexT& operator=(int i);
	/*@}*/

	/** \name Addition */
	/*@{*/
	friend ComplexT operator+(const ComplexT& z1 , const ComplexT& z2) ;
	friend ComplexT operator+(double re1, const ComplexT& z2);
	ComplexT& operator+=(const ComplexT& zRHS);
	/*@}*/

	/** \name Subtraction */
	/*@{*/
	friend ComplexT operator-(const ComplexT& z1, const ComplexT& z2) ;
	friend ComplexT operator-(double re1, const ComplexT& z2);
	ComplexT& operator-=(const ComplexT& zRHS);
	/*@}*/

	/** \name Multiplication */
	/*@{*/
	friend ComplexT operator*(const ComplexT& z1, const ComplexT& z2);
	friend ComplexT operator*(double re1, const ComplexT& z2);
	friend ComplexT operator*(const ComplexT& z2,double re1);
	ComplexT& operator*=(const ComplexT& zRHS);
	/*@}*/

	/** \name Division */
	/*@{*/
	friend ComplexT operator/(const ComplexT& z1, const ComplexT& z2);
	friend ComplexT operator/(double re1, const ComplexT& z2);
	friend ComplexT operator/(const ComplexT& z2,double re1);
	ComplexT& operator/=(const ComplexT& zRHS);
	/*@}*/
	
	/** \name I/O */
	/*@{*/
	friend ostream& operator<<(ostream& out, const ComplexT& z);
	friend istream& operator>>(istream&  in, ComplexT& z);
	/*@}*/

	/** \name other Math functions */
	/*@{*/
	friend ComplexT log(const ComplexT& z);
	ComplexT& log_of(const ComplexT& z);
	/*@}*/
	
private:

	double fRe;
	double fIm;

};

/* inline functions */

/* Constructors */
inline ComplexT::ComplexT(void): fRe(0.0), fIm(0.0) { }
inline ComplexT::ComplexT(int i): fRe(i), fIm(0.0) { }
inline ComplexT::ComplexT(double re): fRe(re), fIm(0.0) { }
inline ComplexT::ComplexT(double re, double im): fRe(re), fIm(im) { }

/* Accessors */
inline double ComplexT::Re(void) const { return(fRe); }
inline double ComplexT::Im(void) const { return(fIm); }
inline ComplexT& ComplexT::toZ(double re, double im)
{
	fRe = re;
	fIm = im;

	return (*this);
}

/* Assignment operators */
inline ComplexT& ComplexT::operator=(const ComplexT& zRHS)
{
	fRe = zRHS.fRe;
	fIm = zRHS.fIm;

	return (*this);
}

inline ComplexT& ComplexT::operator=(double re)
{
	fRe = re;
	fIm = 0.0;

	return (*this);
}

inline ComplexT& ComplexT::operator=(int i)
{
	fRe = i;
	fIm = 0.0;

	return (*this);
}

/* Conjugate */
inline ComplexT& ComplexT::Conjugate(const ComplexT& z)
{
  fRe = z.fRe;
  fIm =-z.fIm;
  return *this;
}

inline ComplexT& ComplexT::Conjugate(void)
{
  fIm = -fIm; 
  return *this;
}

/* Addition */
inline ComplexT operator+(const ComplexT& z1, const ComplexT& z2)
{
	return ( ComplexT(z1.fRe + z2.fRe, z1.fIm + z2.fIm) );
}

inline ComplexT operator+(double re1, const ComplexT& z2)
{
	return( z2 + re1 );
}

inline ComplexT& ComplexT::operator+=(const ComplexT& zRHS)
{
	fRe += zRHS.fRe;
	fIm += zRHS.fIm;

	return (*this);
}

/* Subtraction */
inline ComplexT operator-(const ComplexT& z1,const ComplexT& z2)
{
	return ( ComplexT(z1.fRe - z2.fRe, z1.fIm - z2.fIm) );
}

inline ComplexT operator-(double re1, const ComplexT& z2)
{
	return( z2 - re1 );
}

inline ComplexT& ComplexT::operator-=(const ComplexT& zRHS)
{
	fRe -= zRHS.fRe;
	fIm -= zRHS.fIm;

	return (*this);
}

/* Multiplication */
inline ComplexT operator*(const ComplexT& z1, const ComplexT& z2)
{
	return ( ComplexT(z1.fRe*z2.fRe - z1.fIm*z2.fIm,
	                  z1.fIm*z2.fRe + z1.fRe*z2.fIm) );
}

inline ComplexT operator*(double re1, const ComplexT& z2)
{
	return( ComplexT(re1*z2.fRe, re1*z2.fIm)  );
}

inline ComplexT operator*(const ComplexT& z1, double re2)
{
	return( ComplexT(z1.fRe*re2, z1.fIm*re2)  );
}

inline ComplexT& ComplexT::operator*=(const ComplexT& zRHS)
{
	double tr = fRe*zRHS.fRe -  fIm*zRHS.fIm;
	double ti = fIm*zRHS.fRe +  fRe*zRHS.fIm;
	
	fRe = tr;
	fIm = ti;

	return (*this);
}

/* Division */
inline ComplexT operator/(const ComplexT& z1, const ComplexT& z2)
{
	double m = z2.fRe*z2.fRe + z2.fIm*z2.fIm;
	return ( ComplexT( (z1.fRe*z2.fRe + z1.fIm*z2.fIm)/m,
	                   (z1.fIm*z2.fRe - z1.fRe*z2.fIm)/m) );
}

inline ComplexT operator/(double re1, const ComplexT& z2)
{
	double m = z2.fRe*z2.fRe + z2.fIm*z2.fIm;
	return( ComplexT( z2.fRe*re1/m, - z2.fIm*re1/m)  );
}

inline ComplexT operator/(const ComplexT& z1, double re2)
{
	return( ComplexT(z1.fRe/re2, z1.fIm/re2)  );
}

inline ComplexT& ComplexT::operator/=(const ComplexT& zRHS)
{
	
	double m  = zRHS.fRe*zRHS.fRe + zRHS.fIm*zRHS.fIm;
	double tr =(fRe*zRHS.fRe + fIm*zRHS.fIm)/m;
	double ti =(fIm*zRHS.fRe - fRe*zRHS.fIm)/m;
	
	fRe = tr;
	fIm = ti;

	return (*this);
}
}

#endif /* _COMPLEX_T_H_ */
