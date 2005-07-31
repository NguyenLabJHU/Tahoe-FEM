/* $Id: VariArrayT.h,v 1.1.1.1 2001-01-25 20:56:22 paklein Exp $ */
/* created: paklein (04/18/1998)                                          */
/* WRAPPER for ArrayT<>'s to add dynamic sizing using                     */
/* some headroom to cut down calls for memory                             */
/* de/re-allocation                                                       */

#ifndef _VARI_ARRAY_T_H_
#define _VARI_ARRAY_T_H_

/* base class */
#include "VariBaseT.h"

template <class TYPE>
class VariArrayT: public VariBaseT<TYPE>
{
public:

	/* constructors */
	VariArrayT(void);
	VariArrayT(int headroom, ArrayT<TYPE>& ward);

	/* set the managed array - can only be set ONCE */
	void SetWard(int headroom, ArrayT<TYPE>& ward);
	
	/* set length of the ward, fill extra space and copy in old
	 * data if specified */
	void SetLength(int length, bool copy_in);
	void SetLength(int length, const TYPE& fill, bool copy_in);

	/* return the current length of the ward */
	int Length(void) const;

	/* free memory of self and ward */
	void Free(void);

private:

	/* the managed array */
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
	SetHeadRoom(headroom);

	/* can only be called once */
	if (!fWard)
		fWard = &ward;
	else
		throw eGeneralFail;
}
	
/* set length of the ward, fill extra space if specified */
template <class TYPE>
inline void VariArrayT<TYPE>::SetLength(int length, bool copy_in)
{
	/* ward must be set */
	if (!fWard) throw eGeneralFail;

	/* use inherited function */
	SetAlias(*fWard, length, copy_in);
}

template <class TYPE>
inline void VariArrayT<TYPE>::SetLength(int length, const TYPE& fill, bool copy_in)
{
	/* ward must be set */
	if (!fWard) throw eGeneralFail;

	/* use inherited function */
	SetAlias(*fWard, length, fill, copy_in);
}

/* return the current length of the ward */
template <class TYPE>
inline int VariArrayT<TYPE>::Length(void) const
{
	/* ward must be set */
	if (!fWard) throw eGeneralFail;

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

#endif /* _VARI_ARRAY_T_H_ */
