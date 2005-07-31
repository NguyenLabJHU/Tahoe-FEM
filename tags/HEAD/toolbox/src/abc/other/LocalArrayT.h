/* $Id: LocalArrayT.h,v 1.1.1.1 2001-01-25 20:56:25 paklein Exp $ */
/* created: paklein (07/10/1996)                                          */

#ifndef _LOCALARRAY_T_H_
#define _LOCALARRAY_T_H_

/* base class */
#include "dArrayT.h"

/* forward declarations */
class dArray2DT;

class LocalArrayT: public dArrayT
{	
public:

	/* array types */
	enum TypeT {kUnspecified,
	             kInitCoords,
                       kDisp,
                        kVel,
                        kAcc,
                 kCurrCoords,
                   kLastDisp,
                    kLastVel,
                   kLastAcc};		

	/* constructors */
	LocalArrayT(TypeT type);
	LocalArrayT(TypeT type, int numnodes, int minordim);
	LocalArrayT(const LocalArrayT& source);

	/* allocating */
	void Allocate(int numnodes, int minordim);
	void Set(int numnodes, int minordim, double*p);
	void SetType(TypeT type);
		
	/* accessors */
	TypeT Type(void) const;
	int NumberOfNodes(void) const;
	int MinorDim(void) const;
	
	/* element accessors */
	double& operator()(int majordim, int minordim) const;
	double* operator()(int minordim) const; //pointer to row

	/* assignment operator - does not set Type */
	LocalArrayT& operator=(const LocalArrayT& RHS);
	LocalArrayT& operator=(const double value);
	void Alias(const LocalArrayT& source);

	/* combining arrays - inserts all of source at start_node */
	void BlockCopyAt(const LocalArrayT& source, int start_node);
	
	/* return the vector with transposed indexing */
	void ReturnTranspose(nArrayT<double>& transpose) const;
	void FromTranspose(const nArrayT<double>& transpose);

	/* for registered arrays - preset source for SetLocal */
	int IsRegistered(void) const;
	void SetGlobal(const dArray2DT& global);
	void SetLocal(const ArrayT<int>& keys);

private:

	TypeT fType;
	int   fNumNodes;
	int   fMinorDim;
	
	/* source for SetLocal */
	const dArray2DT* fGlobal;
};

/* in-lines */

/* allocating */
inline void LocalArrayT::Allocate(int numnodes, int minordim)
{
	fNumNodes = numnodes;
	fMinorDim = minordim;
	
	/* call single argument allocate function */
	dArrayT::Allocate(fNumNodes*fMinorDim);
}

inline void LocalArrayT::Set(int numnodes, int minordim, double *p)
{
	fNumNodes = numnodes;
	fMinorDim = minordim;
	
	/* inherited */
	dArrayT::Set(fNumNodes*fMinorDim, p);
}

inline void LocalArrayT::SetType(TypeT type) { fType = type; }

/* accessors */
inline LocalArrayT::TypeT LocalArrayT::Type(void) const    { return fType;     }
inline int LocalArrayT::NumberOfNodes(void) const { return fNumNodes; }
inline int LocalArrayT::MinorDim(void) const      { return fMinorDim; }

/* element accessors */
inline double& LocalArrayT::operator()(int majordim, int minordim) const
{
#if __option (extended_errorcheck)
	if (majordim < 0 ||
	    majordim >= fNumNodes ||
	    minordim < 0 ||
	    minordim >= fMinorDim) throw eOutOfRange;
#endif

	return (fArray[minordim*fNumNodes + majordim]);
}

inline double* LocalArrayT::operator()(int majordim) const
{
#if __option (extended_errorcheck)
	if (majordim < 0 || majordim >= fNumNodes) throw eOutOfRange;
#endif

	return (fArray + majordim*fNumNodes);
}

/* assignment operator */
inline LocalArrayT& LocalArrayT::operator=(const double value)
{
	/* inherited */
	dArrayT::operator=(value);
	return *this;
}

/* for registered arrays - preset source for SetLocal */
inline int LocalArrayT::IsRegistered(void) const { return (fGlobal != NULL); }

/* make a shallow copy */
inline void LocalArrayT::Alias(const LocalArrayT& source)
{
	/* inherited */
	dArrayT::Alias(source);

	/* additional data */
	fType     = source.fType;
	fNumNodes = source.fNumNodes;
	fMinorDim = source.fMinorDim;
	fGlobal   = source.fGlobal;
}

#endif /* _LOCALARRAY_T_H_ */
