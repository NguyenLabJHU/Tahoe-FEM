// $Id: NOX_Tahoe_Vector.cpp,v 1.6 2002-10-20 22:49:48 paklein Exp $
#include "NOX_Tahoe_Vector.h"

/* optional */
#ifdef __NOX__

#include "dArrayT.h"

using namespace Tahoe;
using namespace NOX;

Vector::Vector(void) { fArray = new dArrayT; }

Vector::Vector(const dArrayT& source, CopyType type):
	fArray(NULL)
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

NOX::Abstract::Vector& Vector::operator=(const dArrayT& source)
{
	*fArray = source;
	return *this;
}

NOX::Abstract::Vector& Vector::operator=(const NOX::Abstract::Vector& source)
{
	return operator=(dynamic_cast<const Vector&>(source));
}

NOX::Abstract::Vector& Vector::operator=(const Vector& source)
{
	*fArray = source;
	return *this;
}

NOX::Abstract::Vector& Vector::init(double value)
{
	*fArray = value;
	return *this;
}

NOX::Abstract::Vector& Vector::abs(const NOX::Abstract::Vector& base)
{
	return abs(dynamic_cast<const Vector&>(base));
}

NOX::Abstract::Vector& Vector::abs(const Vector& base)
{
	const dArrayT& src = base;
#if __option(extended_errorcheck)
	if (fArray->Length() != src.Length()) throw ExceptionT::kSizeMismatch;
#endif
	for (int i = 0; i < src.Length(); i++)
		(*fArray)[i] = fabs(src[i]);
	return *this;
}

NOX::Abstract::Vector& Vector::reciprocal(const NOX::Abstract::Vector& base)
{
	return reciprocal(dynamic_cast<const Vector&>(base));
}

NOX::Abstract::Vector& Vector::reciprocal(const Vector& base)
{
	const dArrayT& src = base;
#if __option(extended_errorcheck)
	if (fArray->Length() != src.Length()) throw ExceptionT::kSizeMismatch;
#endif
	for (int i = 0; i < src.Length(); i++)
		(*fArray)[i] = 1.0/src[i];
	return *this;
}

NOX::Abstract::Vector& Vector::scale(double alpha)
{
	*fArray *= alpha;
	return *this;
}

NOX::Abstract::Vector& Vector::update(double alpha, const NOX::Abstract::Vector& a, double gamma)
{
	return update(alpha, dynamic_cast<const Vector&>(a), gamma);
}

NOX::Abstract::Vector& Vector::update(double alpha, const Vector& a, double gamma)
{
	fArray->SetToCombination(alpha, a, gamma, *this);
	return *this;
}

NOX::Abstract::Vector& Vector::update(double alpha, const NOX::Abstract::Vector& a, 
				 double beta, const NOX::Abstract::Vector& b,
				 double gamma)
{
  	return update(alpha, dynamic_cast<const Vector&>(a), 
                   beta, dynamic_cast<const Vector&>(b), gamma);
}

NOX::Abstract::Vector& Vector::update(double alpha, const Vector& a, 
				 double beta, const Vector& b,
				 double gamma)
{
	fArray->SetToCombination(alpha, a, beta, b, gamma, *this);
	return *this;
}


NOX::Abstract::Vector* Vector::clone(NOX::CopyType type) const
{
	Vector* newVec = new Vector(*fArray, type);
	return newVec;
}

double Vector::norm(NOX::Abstract::Vector::NormType type) const
{
	double n = 0.0;
	switch (type) {
		case NOX::Abstract::Vector::INF:
			n = fArray->AbsMax();
			break;
		case NOX::Abstract::Vector::ONE:
			n = fArray->AbsSum();
			break;
		case NOX::Abstract::Vector::TWO:
		default:
			n = fArray->Magnitude();
			break;
	}
	return n;
}

double Vector::norm(const NOX::Abstract::Vector& weights, NOX::Abstract::Vector::NormType type) const
{
	return norm(dynamic_cast<const Vector&>(weights), type);
}

double Vector::norm(const Vector& weights, NOX::Abstract::Vector::NormType type) const
{
	double n = 0.0;
	switch (type) {
		case INF:
		case ONE:
			cerr << "\n Vector::norm: type not supported: " << type << endl;
			throw ExceptionT::kGeneralFail;
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
	/* dimension */
	if (!fArray)
		fArray = new dArrayT(source.Length());
	else
		fArray->Dimension(source.Length());

	/* clear values */
	*fArray = 0.0;
}

#endif /* __NOX__ */
