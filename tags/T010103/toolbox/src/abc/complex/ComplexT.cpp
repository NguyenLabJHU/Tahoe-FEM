/* $Id: ComplexT.cpp,v 1.13 2002-10-20 22:38:52 paklein Exp $ */
/* created: PAK/AFLP (05/19/1997) */

#include "ComplexT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include "toolboxConstants.h"
#include "nArrayT.h"

/* array behavior */

using namespace Tahoe;

namespace Tahoe {
const bool ArrayT<ComplexT>::fByteCopy = true;
} /* namespace Tahoe */

/*
* Real and Imaginary parts of arrays - must be dimensioned BEFORE call
*/
void ComplexT::z_to_Re(const nArrayT<ComplexT>& z, nArrayT<double>& d)
{
	/* dimension check */
	if ( z.Length() != d.Length() ) throw ExceptionT::kOutOfRange;
	
	ComplexT* pz = z.Pointer();
	double*   pd = d.Pointer();
	
	for (int i = 0; i < z.Length(); i++)
		*pd++ = (*pz++).Re();
}

void ComplexT::z_to_Im(const nArrayT<ComplexT>& z, nArrayT<double>& d)
{
	/* dimension check */
	if ( z.Length() != d.Length() ) throw ExceptionT::kOutOfRange;
	
	ComplexT* pz = z.Pointer();
	double*   pd = d.Pointer();
	
	for (int i = 0; i < z.Length(); i++)
		*pd++ = (*pz++).Im();
}

void ComplexT::ReIm_to_z(const nArrayT<double>& re, const nArrayT<double>& im,
	nArrayT<ComplexT>& z)	
{
	/* dimension check */
	if ( re.Length() != im.Length() || im.Length() != z.Length() ) throw ExceptionT::kOutOfRange;

	ComplexT* pz  = z.Pointer();
	double*   pre = re.Pointer();
	double*   pim = im.Pointer();
	
	for (int i = 0; i < z.Length(); i++)
		(pz++)->toZ(*pre++,*pim++);
}

/* Polar components */
double ComplexT::Magnitude() const
{
	return ( sqrt(fRe*fRe + fIm*fIm) );
}

double ComplexT::Angle() const
{
	return ( atan2(fIm,fRe) );
}

namespace Tahoe {
/* I/O */
ostream& operator<<(ostream& out, const ComplexT& z)
{
	out << setw(kDoubleWidth) << z.fRe << " + i ";
	out << setw(kDoubleWidth) << z.fIm;
	
	return (out);
}

istream& operator>>(istream& in, ComplexT& z)
{
	in >> z.fRe >> z.fIm;	
	return(in);
}

/* other Math functions */
ComplexT log(const ComplexT& z)
{
	return ( ComplexT ( ::log( z.Magnitude() ) ,z.Angle() ) );
}

} // namespace Tahoe

ComplexT& ComplexT::log_of(const ComplexT& z)
{

	extern double log(double);

	fRe = /*std::*/log( z.Magnitude() );
	fIm = z.Angle();

	return (*this);
}

