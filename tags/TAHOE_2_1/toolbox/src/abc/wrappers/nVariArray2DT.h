/* $Id: nVariArray2DT.h,v 1.10 2004-05-12 17:21:01 paklein Exp $ */
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
	nVariArray2DT(nArray2DT<nTYPE>& ward);
	nVariArray2DT(int headroom, nArray2DT<nTYPE>& ward, int minordim);
	/*@}*/

	/** set the managed array - can only be set once */
	void SetWard(int headroom, nArray2DT<nTYPE>& ward, int minordim);
	
	/** return true if the ward is already set */
	bool HasWard(void) const { return fWard != NULL; };
	
	/** \name set length of the ward
	 * Fills extra space and copy old data if specified */
	/*@{*/
	void Dimension(int majordim, int minordim);
	void Dimension(const nArray2DT<nTYPE>& array2D) { Dimension(array2D.MajorDim(), array2D.MinorDim()); };
	void SetMajorDimension(int majordim, const nTYPE& fill, bool copy_in);
	void SetMajorDimension(int majordim, bool copy_in);

	/** exchange memory with the source array */
	void Swap(nArray2DT<nTYPE>& source);
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
	fMinorDim(0),
	fWard(NULL)
{ 

}

template <class nTYPE>
nVariArray2DT<nTYPE>::nVariArray2DT(nArray2DT<nTYPE>& ward):
	fMinorDim(0),
	fWard(NULL)
{
	SetWard(0, ward, 0);
}

template <class nTYPE>
nVariArray2DT<nTYPE>::nVariArray2DT(int headroom,
	nArray2DT<nTYPE>& ward, int minordim): 
	fMinorDim(0),
	fWard(NULL)
{
	SetWard(headroom, ward, minordim);
}

/* set the managed array - can only be set ONCE */
template <class nTYPE>
void nVariArray2DT<nTYPE>::SetWard(int headroom, nArray2DT<nTYPE>& ward,
	int minordim)
{
	const char caller[] = "nVariArray2DT<nTYPE>::SetWard";
	
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
			if (fWard->MinorDim() != fMinorDim) ExceptionT::SizeMismatch(caller);
		}
		else
			/* set minor dimension */
			fWard->Set(0, fMinorDim, NULL);
	}
	else
		ExceptionT::GeneralFail(caller, "ward already set");
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
	if (!fWard) ExceptionT::GeneralFail("nVariArray2DT<nTYPE>::SetMajorDimension", "ward not set");

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
	if (!fWard) ExceptionT::GeneralFail("nVariArray2DT<nTYPE>::SetMajorDimension", "ward not set");

	/* update ArrayT data */
	SetAlias(*fWard, majordim*fMinorDim, fill, copy_in);

	/* update rest */
	fWard->Set(majordim, fMinorDim, fWard->Pointer());
}

/* exchange memory with the source array */
template <class nTYPE>
void nVariArray2DT<nTYPE>::Swap(nArray2DT<nTYPE>& source)
{
	/* ward must be set */
	if (!fWard) ExceptionT::GeneralFail("nVariArray2DT<nTYPE>::Swap", "ward not set");

	/* current dimensions */
	int major_dim = source.MajorDim();
	int minor_dim = source.MinorDim();

	/* swap memory and dimension (yuk!) */
	VariBaseT<nTYPE>::Swap(source);
	nTYPE* mem;
	source.ReleasePointer(&mem);
	source.Set(source.Length()/fMinorDim, fMinorDim, mem);
	source.TakePointer(source.Length(), mem);
	
	/* set ward */
	fMinorDim = minor_dim;
	fWard->Free(); /* force reset */
	SetAlias(*fWard, major_dim*fMinorDim, true);
	fWard->Set(major_dim, fMinorDim, fWard->Pointer());
}

/* dimensions accessors - of the ward */
template <class nTYPE>
inline int nVariArray2DT<nTYPE>::MajorDim(void) const
{
	/* ward must be set */
	if (!fWard) ExceptionT::GeneralFail("nVariArray2DT<nTYPE>::MajorDim", "ward not set");

	return(fWard->MajorDim());
}

template <class nTYPE>
inline int nVariArray2DT<nTYPE>::MinorDim(void) const
{
	/* ward must be set */
	if (!fWard) ExceptionT::GeneralFail("nVariArray2DT<nTYPE>::MinorDim", "ward not set");

	return(fWard->MinorDim());
}

/* reference to the ward */
template <class nTYPE>
const nArray2DT<nTYPE>& nVariArray2DT<nTYPE>::TheWard(void) const
{
	/* ward must be set */
	if (!fWard) ExceptionT::GeneralFail("nVariArray2DT<nTYPE>::TheWard", "ward not set");

	return(*fWard);
}

} // namespace Tahoe 
#endif /* _N_VARI_ARRAY2D_T_H_ */
