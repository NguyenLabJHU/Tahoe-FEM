/* $Id: nVariArray2DT.h,v 1.7 2002-11-25 07:04:16 paklein Exp $ */
/* created: paklein (04/18/1998) */
#ifndef _N_VARI_ARRAY2D_T_H_
#define _N_VARI_ARRAY2D_T_H_

/* base class */
#include "VariBaseT.h"

/* direct members */
#include "nArray2DT.h"

namespace Tahoe {

/** wrapper for nArray2DT<>'s. Adds dynamic re-sizing of the
 * major dimension, using some headroom to cut down calls for
 * memory de/re-allocation
 */
template <class nTYPE>
class nVariArray2DT: public VariBaseT<nTYPE>
{
public:

	/** \name constructors */
	/*@{*/
	nVariArray2DT(void);
	nVariArray2DT(int headroom, nArray2DT<nTYPE>& ward, int minordim);
	/*@}*/

	/** set the managed array - can only be set once */
	void SetWard(int headroom, nArray2DT<nTYPE>& ward, int minordim);
	
	/** \name set length of the ward
	 * Fills extra space and copy old data if specified */
	/*@{*/
	void Dimension(int majordim, int minordim);
	void SetMajorDimension(int majordim, const nTYPE& fill, bool copy_in);
	void SetMajorDimension(int majordim, bool copy_in);
	/*@}*/

	/** \name dimensions of the ward */
	/*@{*/
	int MajorDim(void) const;
	int MinorDim(void) const;
	/*@}*/
	
	/** reference to the ward */
	const nArray2DT<nTYPE>& TheWard(void) const;
		
private:

	/** \name the managed array */
	/*@{*/
	int fMinorDim;
	nArray2DT<nTYPE>* fWard;
	/*@}*/
};

/*************************************************************************
* Implementation
*************************************************************************/

/* constructors */
template <class nTYPE>
nVariArray2DT<nTYPE>::nVariArray2DT(void): 
	fWard(NULL),
	fMinorDim(0)
{ 

}

template <class nTYPE>
nVariArray2DT<nTYPE>::nVariArray2DT(int headroom,
	nArray2DT<nTYPE>& ward, int minordim): 
	fWard(NULL),
	fMinorDim(0)
{
	SetWard(headroom, ward, minordim);
}

/* set the managed array - can only be set ONCE */
template <class nTYPE>
void nVariArray2DT<nTYPE>::SetWard(int headroom, nArray2DT<nTYPE>& ward,
	int minordim)
{
	/* inherited */
	SetHeadRoom(headroom);

	/* can only be called once */
	if (!fWard)
	{
		fMinorDim = minordim;
		fWard     = &ward;
		if (fWard->MinorDim() > 0)
		{
			/* consistency check */
			if (fWard->MinorDim() != fMinorDim) throw ExceptionT::kSizeMismatch;
		}
		else
			/* set minor dimension */
			fWard->Set(0, fMinorDim, NULL);
	}
	else
		throw ExceptionT::kGeneralFail;
}
	
/* set length of the ward, fill extra space if specified */
template <class nTYPE>
inline void nVariArray2DT<nTYPE>::Dimension(int majordim, int minordim)
{
	/* set minor dimension */
	fMinorDim = minordim;

	/* set rest (don't copy old data) */
	SetMajorDimension(majordim, false);
}

template <class nTYPE>
inline void nVariArray2DT<nTYPE>::SetMajorDimension(int majordim, bool copy_in)
{
	/* ward must be set */
	if (!fWard) throw ExceptionT::kGeneralFail;

	/* update ArrayT data */
	SetAlias(*fWard, majordim*fMinorDim, copy_in);

	/* update rest */
	fWard->Set(majordim, fMinorDim, fWard->Pointer());
}

template <class nTYPE>
inline void nVariArray2DT<nTYPE>::SetMajorDimension(int majordim,
	const nTYPE& fill, bool copy_in)
{
	/* ward must be set */
	if (!fWard) throw ExceptionT::kGeneralFail;

	/* update ArrayT data */
	SetAlias(*fWard, majordim*fMinorDim, fill, copy_in);

	/* update rest */
	fWard->Set(majordim, fMinorDim, fWard->Pointer());
}

/* dimensions accessors - of the ward */
template <class nTYPE>
inline int nVariArray2DT<nTYPE>::MajorDim(void) const
{
	/* ward must be set */
	if (!fWard) throw ExceptionT::kGeneralFail;

	return(fWard->MajorDim());
}

template <class nTYPE>
inline int nVariArray2DT<nTYPE>::MinorDim(void) const
{
	/* ward must be set */
	if (!fWard) throw ExceptionT::kGeneralFail;

	return(fWard->MinorDim());
}

/* reference to the ward */
template <class nTYPE>
const nArray2DT<nTYPE>& nVariArray2DT<nTYPE>::TheWard(void) const
{
	/* ward must be set */
	if (!fWard) throw ExceptionT::kGeneralFail;

	return(*fWard);
}

} // namespace Tahoe 
#endif /* _N_VARI_ARRAY2D_T_H_ */
