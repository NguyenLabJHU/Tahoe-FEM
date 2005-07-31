/* $Id: ElementMatrixT.h,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (03/17/1998)                                          */

#ifndef _ELEMENT_MATRIX_T_H_
#define _ELEMENT_MATRIX_T_H_

/* base class */
#include "dMatrixT.h"

class ElementMatrixT: public dMatrixT
{
public:

	/* format flags */
	enum FormatT {kNonSymmetric, kSymmetric, kSymmetricUpper, kDiagonal};

	/* constructors */
	ElementMatrixT(FormatT format);
	ElementMatrixT(int numrows, int numcols, FormatT format);
	ElementMatrixT(int squaredim, FormatT format);
	ElementMatrixT(const ElementMatrixT& source);
	
	/* format accessor */
	FormatT Format(void) const;
	
	/* (re-)setting the format */
	void SetFormat(FormatT format);
		
	/* assignment operators */
	ElementMatrixT& operator=(const dMatrixT& RHS);
	ElementMatrixT& operator=(const double value);
	
	/* take symmetric with values in upper only and copy to full */
	/* NOTE: error to call for kNonSymmetric matrices            */
	/* NOTE: not strictly const, but in essence                  */
	void CopySymmetric(void) const;

private:

	FormatT fFormat;

};

/* inlines */

/* format accessor */
inline ElementMatrixT::FormatT ElementMatrixT::Format(void) const { return fFormat; }
inline void ElementMatrixT::SetFormat(FormatT format)
{
	fFormat = format;
}

/* assigment operators */
inline ElementMatrixT& ElementMatrixT::operator=(const dMatrixT& RHS)
{
	dMatrixT::operator=(RHS);
	return *this;
}

inline ElementMatrixT& ElementMatrixT::operator=(const double value)
{
	dMatrixT::operator=(value);
	return *this;
}

#endif /* _ELEMENT_MATRIX_T_H_ */
