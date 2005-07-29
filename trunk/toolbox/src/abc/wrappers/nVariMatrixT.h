/* $Id: nVariMatrixT.h,v 1.9 2005-07-29 03:09:35 paklein Exp $ */
/* created: paklein (07/05/1998) */
#ifndef _N_VARI_MATRIX_T_H_
#define _N_VARI_MATRIX_T_H_

/* base class */
#include "VariBaseT.h"

/* direct members */
#include "nMatrixT.h"

namespace Tahoe {

/** wrapper for nMatrixT<>'s for dynamic re-sizing.
 * Manage changing memory using some headroom to cut down calls for
 * memory de/re-allocation */
template <class nTYPE>
class nVariMatrixT: public VariBaseT<nTYPE>
{
public:

	/** \name constructors */
	/*@{*/
	nVariMatrixT(void);
	nVariMatrixT(int headroom, nMatrixT<nTYPE>& ward);
	/*@}*/

	/** set the managed array. \note can only be set once */
	void SetWard(int headroom, nMatrixT<nTYPE>& ward);

	/** return true if the ward is already set */
	bool HasWard(void) const { return fWard != NULL; };
	
	/** \name set dimensions of the ward */
	/*@{*/
	void SetDimensions(int rows, int cols);
	void SetDimensions(int squaredim);
	/*@}*/

	/** \name dimension accessors of the ward */
	/*@{*/
	int Rows(void) const;
	int Cols(void) const;
	/*@}*/
	
	/** reference to the ward */
	const nMatrixT<nTYPE>& TheWard(void) const;
		
private:

	/** the managed array */
	nMatrixT<nTYPE>* fWard;
};

/*************************************************************************
 * Implementation
 *************************************************************************/

/* constructors */
template <class nTYPE>
nVariMatrixT<nTYPE>::nVariMatrixT(void): fWard(NULL) { }

template <class nTYPE>
nVariMatrixT<nTYPE>::nVariMatrixT(int headroom,
	nMatrixT<nTYPE>& ward):
	VariBaseT<nTYPE>(headroom),
	fWard(&ward)
{

}

/* set the managed array - can only be set ONCE */
template <class nTYPE>
void nVariMatrixT<nTYPE>::SetWard(int headroom, nMatrixT<nTYPE>& ward)
{
	/* inherited */
	this->SetHeadRoom(headroom);

	/* can only be called once */
	if (!fWard)
		fWard = &ward;
	else
		ExceptionT::GeneralFail("nVariMatrixT");
}
	
/* set length of the ward, fill extra space if specified */
template <class nTYPE>
inline void nVariMatrixT<nTYPE>::SetDimensions(int rows, int cols)
{
	/* ward must be set */
	if (!fWard) ExceptionT::GeneralFail("nVariMatrixT");

	/* update ArrayT data (don't copy old data) */
	SetAlias(*fWard, rows*cols, false);

	/* update rest */
	fWard->Set(rows, cols, fWard->Pointer());
}

template <class nTYPE>
inline void nVariMatrixT<nTYPE>::SetDimensions(int squaredim)
{
	/* bounce */
	SetDimensions(squaredim, squaredim);
}

/* dimensions accessors - of the ward */
template <class nTYPE>
inline int nVariMatrixT<nTYPE>::Rows(void) const
{
	/* ward must be set */
	if (!fWard) ExceptionT::GeneralFail("nVariMatrixT");

	return(fWard->Rows());
}

template <class nTYPE>
inline int nVariMatrixT<nTYPE>::Cols(void) const
{
	/* ward must be set */
	if (!fWard) ExceptionT::GeneralFail("nVariMatrixT");

	return(fWard->Cols());
}

/* reference to the ward */
template <class nTYPE>
const nMatrixT<nTYPE>& nVariMatrixT<nTYPE>::TheWard(void) const
{
	/* ward must be set */
	if (!fWard) ExceptionT::GeneralFail("nVariMatrixT");

	return(*fWard);
}

} // namespace Tahoe 
#endif /* _N_VARI_MATRIX_T_H_ */
