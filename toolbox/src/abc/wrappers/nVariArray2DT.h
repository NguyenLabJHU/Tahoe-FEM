/* $Id: nVariArray2DT.h,v 1.1.1.1 2001-01-25 20:56:22 paklein Exp $ */
/* created: paklein (04/18/1998)                                          */
/* WRAPPER for nArray2DT<>'s to add dynamic re-sizing of the              */
/* major dimension, using some headroom to cut down calls for             */
/* memory de/re-allocation                                                */

#ifndef _N_VARI_ARRAY2D_T_H_
#define _N_VARI_ARRAY2D_T_H_

/* base class */
#include "VariBaseT.h"

/* direct members */
#include "nArray2DT.h"

template <class nTYPE>
class nVariArray2DT: public VariBaseT<nTYPE>
{
public:

	/* constructors */
	nVariArray2DT(void);
	nVariArray2DT(int headroom, nArray2DT<nTYPE>& ward, int minordim);

	/* set the managed array - can only be set once */
	void SetWard(int headroom, nArray2DT<nTYPE>& ward, int minordim);
	
	/* set length of the ward, fill extra space and copy old data if specified */
	void Dimension(int majordim, int minordim);
	void SetMajorDimension(int majordim, const nTYPE& fill, bool copy_in);
	void SetMajorDimension(int majordim, bool copy_in);

	/* dimensions accessors - of the ward */
	int MajorDim(void) const;
	int MinorDim(void) const;
	
	/* reference to the ward */
	const nArray2DT<nTYPE>& TheWard(void) const;
		
private:

	/* the managed array */
	int fMinorDim;
	nArray2DT<nTYPE>* fWard;
};

/*************************************************************************
* Implementation
*************************************************************************/

/* constructors */
template <class nTYPE>
nVariArray2DT<nTYPE>::nVariArray2DT(void): fWard(NULL) { }

template <class nTYPE>
nVariArray2DT<nTYPE>::nVariArray2DT(int headroom,
	nArray2DT<nTYPE>& ward, int minordim): fWard(NULL)
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
			if (fWard->MinorDim() != fMinorDim) throw eSizeMismatch;
		}
		else
			/* set minor dimension */
			fWard->Set(0, fMinorDim, NULL);
	}
	else
		throw eGeneralFail;
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
	if (!fWard) throw eGeneralFail;

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
	if (!fWard) throw eGeneralFail;

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
	if (!fWard) throw eGeneralFail;

	return(fWard->MajorDim());
}

template <class nTYPE>
inline int nVariArray2DT<nTYPE>::MinorDim(void) const
{
	/* ward must be set */
	if (!fWard) throw eGeneralFail;

	return(fWard->MinorDim());
}

/* reference to the ward */
template <class nTYPE>
const nArray2DT<nTYPE>& nVariArray2DT<nTYPE>::TheWard(void) const
{
	/* ward must be set */
	if (!fWard) throw eGeneralFail;

	return(*fWard);
}

#endif /* _N_VARI_ARRAY2D_T_H_ */
