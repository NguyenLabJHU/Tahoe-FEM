/* $Id: pArrayT.h,v 1.1.1.1 2001-01-25 20:56:23 paklein Exp $ */
/* created: paklein (11/21/1996)                                          */
/* This is the interface for a pointer array class.  The data in the      */
/* array is of type (TYPE*).  The members of the array is deleted when    */
/* the array is destructed.                                               */
/* Note: This template cannot be instantiated for a non-pointer type.     */
/* The memory management for the class calls delete for every             */
/* member of the array.                                                   */

#ifndef _P_ARRAY_T_H_
#define _P_ARRAY_T_H_

/* base class */
#include "ArrayT.h"

/* forward declarations */
template <class TYPEPtr> class ProxyTYPEPtr;

template <class TYPEPtr>
class pArrayT: public ArrayT<TYPEPtr>
{
public:

	/* constructors */
	pArrayT(void);
	pArrayT(int length);
	pArrayT(const pArrayT& source);

	/* destructor */
	~pArrayT(void);

	/* allocate an array of the specified size */
	void Allocate(int length);

	/* element accessor */
	ProxyTYPEPtr<TYPEPtr> operator[](int index);
	ProxyTYPEPtr<TYPEPtr> operator[](int index) const;

private:

	/* no shallow copies */
	void ShallowCopy(const pArrayT& RHS);

	/* no assigment operator */	 			  	
	void operator=(const pArrayT& RHS);

	/* no resizing */
	void Resize(void);

	/* delete all members of the array */
	void DeleteAll(void);
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
	Allocate(length); 
}

template <class TYPEPtr>
pArrayT<TYPEPtr>::pArrayT(const pArrayT& source)
{
	operator=(source); 
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
void pArrayT<TYPEPtr>::Allocate(int length)
{
	/* reallocate if needed */
	if (fLength != length)
	{
		/* free existing array */
		if (fLength > 0) DeleteAll();

		/* allocate to new size */
		ArrayT<TYPEPtr>::Allocate(length);

		/* NULL all pointers */
		for (int i = 0; i < fLength; i++)
			fArray[i] = NULL;
	}
}

/* element accessor */
template <class TYPEPtr>
inline ProxyTYPEPtr<TYPEPtr> pArrayT<TYPEPtr>::operator[](int index)
{
/* range checking */
#if __option (extended_errorcheck)
	if (index < 0 || index >= fLength) throw(eOutOfRange);
#endif

	return ProxyTYPEPtr<TYPEPtr>(*this,index);
}

template <class TYPEPtr>
inline ProxyTYPEPtr<TYPEPtr> pArrayT<TYPEPtr>::operator[](int index) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (index < 0 || index >= fLength) throw(eOutOfRange);
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

#endif /* _P_ARRAY_T_H_ */
