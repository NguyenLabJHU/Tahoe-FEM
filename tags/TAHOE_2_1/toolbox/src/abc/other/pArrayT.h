/* $Id: pArrayT.h,v 1.7 2004-08-08 02:01:44 paklein Exp $ */
/* created: paklein (11/21/1996) */

#ifndef _P_ARRAY_T_H_
#define _P_ARRAY_T_H_

/* base class */
#include "ArrayT.h"

namespace Tahoe {

/* forward declarations */
template <class TYPEPtr> class ProxyTYPEPtr;

/** A class to help working with arrays of pointers. The data in the
 * array is of type (TYPE*). The pointers should be allocated with
 * operator new, not operator new []. The elements of the array are
 * initialized to NULL when the arrays are dimensioned. delete is
 * called for every element in the array in the class destructor. 
 * delete is also called on elements of the array before they are
 * over written with operator [].
 * Note: This template cannot be instantiated for a non-pointer type.
 * The memory management for the class calls delete for every
 * member of the array. */
template <class TYPEPtr>
class pArrayT: public ArrayT<TYPEPtr>
{
public:

	/** constructor zero length array */
	pArrayT(void);

	/** construct an array of the given length. All element of the
	 * array are initialized to NULL */
	explicit pArrayT(int length);

	/* destructor */
	~pArrayT(void);

	/** set the array size to the given length. No change occurs if the array
	 * is already the specified length. The previous contents of the array is
	 * not preserved. Use pArrayT::Resize to preserve the array contents. */
	void Dimension(int length);

	/** dimension the array to the new length keeping as much of the previous
	 * data as fits in the new space */
	void Resize(int new_length);

	/** \deprecated replaced by pArrayT::Dimension on 02/13/2002 */
	void Allocate(int length) { Dimension(length); };

	/** element accessor */
	ProxyTYPEPtr<TYPEPtr> operator[](int index);

	/** const element accessor */
	ProxyTYPEPtr<TYPEPtr> operator[](int index) const;

private:

	/** call delete for all members of the array */
	void DeleteAll(void);

	/** no shallow copies */
	void ShallowCopy(const pArrayT& RHS);

	/** no assigment operator */	 			  	
	void operator=(const pArrayT& RHS);
};

/* proxy - for element accessor */
template <class TYPEPtr>
class ProxyTYPEPtr
{
public:
	
	/* constructor */
	ProxyTYPEPtr(pArrayT<TYPEPtr>& parent, int dex);
	
	/* lvalue - no chaining */
	void operator=(TYPEPtr typeptr);
	
	/* rvalue - type conversion */
	operator TYPEPtr();

	/* rvalue - smart pointer */
	TYPEPtr operator->(); //CW wouldn't call functions with conversion

private:

	pArrayT<TYPEPtr>& fParent;
	int fDex;
};

/*************************************************************************
* Implementation
*************************************************************************/

/****** Proxy ******/

/* constructor */
template <class TYPEPtr>
inline ProxyTYPEPtr<TYPEPtr>::ProxyTYPEPtr(pArrayT<TYPEPtr>& parent, int dex):
	fParent(parent),
	fDex(dex)
{

}	

/* lvalue */
template <class TYPEPtr>
void ProxyTYPEPtr<TYPEPtr>::operator=(TYPEPtr typeptr)
{
	TYPEPtr* p = fParent.Pointer();
	
	/* free memory */
	if (p[fDex] != NULL) delete p[fDex];
	
	/* set value */
	p[fDex] = typeptr;
}
	
/* rvalue - type conversion */
template <class TYPEPtr>
inline ProxyTYPEPtr<TYPEPtr>::operator TYPEPtr()
{
	TYPEPtr* p = fParent.Pointer();
	return p[fDex];
}

/* rvalue - smart pointer */
template <class TYPEPtr>
inline TYPEPtr ProxyTYPEPtr<TYPEPtr>::operator->()
{
	return(ProxyTYPEPtr<TYPEPtr>::operator TYPEPtr());
}

/****** Proxy ******/

/* constructor */
template <class TYPEPtr>
pArrayT<TYPEPtr>::pArrayT(void) { }

template <class TYPEPtr>
pArrayT<TYPEPtr>::pArrayT(int length)
{
	Dimension(length); 
}

/* destructor */
template <class TYPEPtr>
pArrayT<TYPEPtr>::~pArrayT(void)
{
	/* delete all members of the list */
	DeleteAll();
}

/* allocate an array of the specified size.  Frees any existing
* memory */
template <class TYPEPtr>
void pArrayT<TYPEPtr>::Dimension(int length)
{
	/* reallocate if needed */
	if (fLength != length)
	{
		/* free existing array */
		if (fLength > 0) DeleteAll();

		/* allocate to new size */
		ArrayT<TYPEPtr>::Dimension(length);

		/* NULL all pointers */
		for (int i = 0; i < fLength; i++)
			fArray[i] = NULL;
	}
}

/* dimension the array to the new length keeping as much of the previous
 * data as fits in the new space */
template <class TYPEPtr>
void pArrayT<TYPEPtr>::Resize(int new_length)
{
	/* array is shrinking - free excess */
	for (int i = fLength; i > new_length; i--)
		delete fArray[i-1];

	/* inherited */
	ArrayT<TYPEPtr>::Resize(new_length, NULL);
}

/* element accessor */
template <class TYPEPtr>
inline ProxyTYPEPtr<TYPEPtr> pArrayT<TYPEPtr>::operator[](int index)
{
/* range checking */
#if __option (extended_errorcheck)
	if (index < 0 || index >= fLength) throw ExceptionT::kOutOfRange;
#endif

	return ProxyTYPEPtr<TYPEPtr>(*this,index);
}

template <class TYPEPtr>
inline ProxyTYPEPtr<TYPEPtr> pArrayT<TYPEPtr>::operator[](int index) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (index < 0 || index >= fLength) throw ExceptionT::kOutOfRange;
#endif

/* const_cast<> not supported */
#ifdef __SUNPRO_CC
	pArrayT<TYPEPtr>* const localthis = (pArrayT<TYPEPtr>* const) this;
	return( ProxyTYPEPtr<TYPEPtr>(*localthis,index) );
#else
	return( ProxyTYPEPtr<TYPEPtr>(*const_cast<pArrayT<TYPEPtr>*>(this),index) );
#endif
}

/*************************************************************************
* Private
*************************************************************************/

/* delete all members of the array */
template <class TYPEPtr>
void pArrayT<TYPEPtr>::DeleteAll(void)
{
	for (int i = 0; i < fLength; i++)
	{
		delete fArray[i];
		fArray[i] = NULL;
	}
}

} // namespace Tahoe 
#endif /* _P_ARRAY_T_H_ */
