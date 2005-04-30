/* $Id: VariArrayT.h,v 1.5 2005-04-30 21:14:41 paklein Exp $ */
/* created: paklein (04/18/1998) */
#ifndef _VARI_ARRAY_T_H_
#define _VARI_ARRAY_T_H_

/* base class */
#include "VariBaseT.h"

namespace Tahoe {

/** wrapper for ArrayT<>'s to add dynamic sizing.
 * Implemented efficient dynamic resizing of an ArrayT using
 * some headroom to cut down calls for memory de/re-allocation */
template <class TYPE>
class VariArrayT: public VariBaseT<TYPE>
{
public:

	/** \name constructors */
	/*@{*/
	VariArrayT(void);
	VariArrayT(int headroom, ArrayT<TYPE>& ward);
	/*@}*/

	/** set the managed array. \note can only be set once. */
	void SetWard(int headroom, ArrayT<TYPE>& ward);

	/** return true if the ward is already set */
	bool HasWard(void) const { return fWard != NULL; };
	
	/** \name set length of the ward
	 * Fill extra space and copy in old data if specified */
	/*@{*/	
	void SetLength(int length, bool copy_in);
	void SetLength(int length, const TYPE& fill, bool copy_in);
	/*@}*/	

	/** return the current length of the ward */
	int Length(void) const;

	/** free memory of self and ward */
	void Free(void);

private:

	/** the managed array */
	ArrayT<TYPE>* fWard;
};

/*************************************************************************
* Implementation
*************************************************************************/

/* constructors */
template <class TYPE>
VariArrayT<TYPE>::VariArrayT(void): fWard(NULL) { }

template <class TYPE>
VariArrayT<TYPE>::VariArrayT(int headroom, ArrayT<TYPE>& ward):
	VariBaseT<TYPE>(headroom),
	fWard(&ward)
{

}

/* set the managed array - can only be set ONCE */
template <class TYPE>
void VariArrayT<TYPE>::SetWard(int headroom, ArrayT<TYPE>& ward)
{
	this->SetHeadRoom(headroom);

	/* can only be called once */
	if (!fWard)
		fWard = &ward;
	else
		throw ExceptionT::kGeneralFail;
}
	
/* set length of the ward, fill extra space if specified */
template <class TYPE>
inline void VariArrayT<TYPE>::SetLength(int length, bool copy_in)
{
	/* ward must be set */
	if (!fWard) throw ExceptionT::kGeneralFail;

	/* use inherited function */
	SetAlias(*fWard, length, copy_in);
}

template <class TYPE>
inline void VariArrayT<TYPE>::SetLength(int length, const TYPE& fill, bool copy_in)
{
	/* ward must be set */
	if (!fWard) throw ExceptionT::kGeneralFail;

	/* use inherited function */
	SetAlias(*fWard, length, fill, copy_in);
}

/* return the current length of the ward */
template <class TYPE>
inline int VariArrayT<TYPE>::Length(void) const
{
	/* ward must be set */
	if (!fWard) throw ExceptionT::kGeneralFail;

	return(fWard->Length());
}

/* free memory of self and ward */
template <class TYPE>
inline void VariArrayT<TYPE>::Free(void)
{
	/* inherited */
	VariBaseT<TYPE>::Free();

	/* the ward */
	fWard->Set(0, NULL);
}

} // namespace Tahoe 
#endif /* _VARI_ARRAY_T_H_ */
