// $Id: NOX_Tahoe_Vector.cpp,v 1.1 2002-03-28 16:40:35 paklein Exp $
#include "NOX_Tahoe_Vector.H"
#include "dArrayT.h"

using namespace NOX;
using namespace NOX::Tahoe;

Vector::Vector(void) { fArray = new dArrayT; }

Vector::Vector(const dArrayT& source, CopyType type)
{
	switch (type) {

		case DeepCopy:
			fArray = new dArrayT(source); 
			break;

		case CopyShape: // dArrayT's have no "shape"
			DimensionTo(source);
			break;  
	}
}

Vector::~Vector() { delete fArray; }

Abstract::Vector& Vector::operator=(const dArrayT& source)
{
	*fArray = source;
	return *this;
}

Abstract::Vector& Vector::operator=(const Abstract::Vector& source)
{
	return operator=(dynamic_cast<const Vector&>(source));
}

Abstract::Vector& Vector::operator=(const Vector& source)
{
	*fArray = source;
	return *this;
}

Abstract::Vector& Vector::init(double value)
{
	*fArray = value;
	return *this;
}

Abstract::Vector& Vector::abs(const Abstract::Vector& base)
{
	return abs(dynamic_cast<const Vector&>(base));
}

Abstract::Vector& Vector::abs(const Vector& base)
{
	const dArrayT& src = base;
#if __option(extended_errorcheck)
	if (fArray->Length() != src.Length()) throw eSizeMismatch;
#endif
	for (int i = 0; i < src.Length(); i++)
		(*fArray)[i] = fabs(src[i]);
	return *this;
}

Abstract::Vector& Vector::reciprocal(const Abstract::Vector& base)
{
	return reciprocal(dynamic_cast<const Vector&>(base));
}

Abstract::Vector& Vector::reciprocal(const Vector& base)
{
	const dArrayT& src = base;
#if __option(extended_errorcheck)
	if (fArray->Length() != src.Length()) throw eSizeMismatch;
#endif
	for (int i = 0; i < src.Length(); i++)
		(*fArray)[i] = 1.0/src[i];
	return *this;
}

Abstract::Vector& Vector::scale(double alpha)
{
	*fArray *= alpha;
	return *this;
}

Abstract::Vector& Vector::update(double alpha, const Abstract::Vector& a, double gamma)
{
	return update(alpha, dynamic_cast<const Vector&>(a), gamma);
}

Abstract::Vector& Vector::update(double alpha, const Vector& a, double gamma)
{
	fArray->SetToCombination(alpha, a, gamma, *this);
	return *this;
}

Abstract::Vector& Vector::update(double alpha, const Abstract::Vector& a, 
				 double beta, const Abstract::Vector& b,
				 double gamma)
{
  	return update(alpha, dynamic_cast<const Vector&>(a), 
                   beta, dynamic_cast<const Vector&>(b), gamma);
}

Abstract::Vector& Vector::update(double alpha, const Vector& a, 
				 double beta, const Vector& b,
				 double gamma)
{
	fArray->SetToCombination(alpha, a, beta, b, gamma, *this);
	return *this;
}


Abstract::Vector* Vector::clone(CopyType type) const
{
	Vector* newVec = new Vector(*fArray, type);
	return newVec;
}

double Vector::norm(Abstract::Vector::NormType type) const
{
	double n = 0.0;
	switch (type) {
		case INF:
			n = fArray->AbsMax();
			break;
		case ONE:
			n = fArray->AbsSum();
			break;
		case TWO:
		default:
			n = fArray->Magnitude();
			break;
	}
	return n;
}

double Vector::norm(const Abstract::Vector& weights, Abstract::Vector::NormType type) const
{
	return norm(dynamic_cast<const Vector&>(weights), type);
}

double Vector::norm(const Vector& weights, Abstract::Vector::NormType type) const
{
	double n = 0.0;
	switch (type) {
		case INF:
		case ONE:
			cerr << "\n Vector::norm: type not supported: " << type << endl;
			throw eGeneralFail;
			break;
		case TWO:
		default:
		{
			const dArrayT& w = weights;
			const dArrayT& a = *fArray;
			for (int i = 0; i < w.Length(); i++)
				n += w[i]*a[i]*a[i];
			n = sqrt(n);
			break;
		}
	}
	return n;
}

double Vector::dot(const Abstract::Vector& y) const
{
	return dot(dynamic_cast<const Vector&>(y));
}

double Vector::dot(const Vector& y) const
{
	return dArrayT::Dot(*this, y);
}

int Vector::length() const
{
	return fArray->Length();
}

/* dimension the vector based on the source. Copies the shape of the source. */
void Vector::DimensionTo(const dArrayT& source)
{
	fArray->Dimension(source.Length());
	*fArray = 0.0;
}
