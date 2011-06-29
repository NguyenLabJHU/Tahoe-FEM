/* $Id: ComplexT.cpp,v 1.18 2005-07-29 03:09:32 paklein Exp $ */
/* created: PAK/AFLP (05/19/1997) */
#include "ComplexT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include "toolboxConstants.h"
#include "nArrayT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<ComplexT>::fByteCopy = true;
} /* namespace Tahoe */

const char caller[] = "ComplexT";

/*
* Real and Imaginary parts of arrays - must be dimensioned BEFORE call
*/
void ComplexT::z_to_Re(const nArrayT<ComplexT>& z, nArrayT<double>& d)
{
	/* dimension check */
	if ( z.Length() != d.Length() ) ExceptionT::OutOfRange(caller);
	
	const ComplexT* pz = z.Pointer();
	double*   pd = d.Pointer();
	
	for (int i = 0; i < z.Length(); i++)
		*pd++ = (*pz++).Re();
}

void ComplexT::z_to_Im(const nArrayT<ComplexT>& z, nArrayT<double>& d)
{
	/* dimension check */
	if ( z.Length() != d.Length() ) ExceptionT::OutOfRange(caller);
	
	const ComplexT* pz = z.Pointer();
	double*   pd = d.Pointer();
	
	for (int i = 0; i < z.Length(); i++)
		*pd++ = (*pz++).Im();
}

void ComplexT::ReIm_to_z(const nArrayT<double>& re, const nArrayT<double>& im,
	nArrayT<ComplexT>& z)	
{
	/* dimension check */
	if ( re.Length() != im.Length() || im.Length() != z.Length() ) ExceptionT::OutOfRange(caller);

	ComplexT* pz = z.Pointer();
	const double* pre = re.Pointer();
	const double* pim = im.Pointer();
	
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
	fRe = ::log( z.Magnitude() );
	fIm = z.Angle();

	return (*this);
}

