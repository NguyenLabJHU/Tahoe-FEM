/* $Id: LocalArrayT.h,v 1.14 2005-01-26 19:53:27 paklein Exp $ */
/* created: paklein (07/10/1996) */

#ifndef _LOCALARRAY_T_H_
#define _LOCALARRAY_T_H_

/* base class */
#include "dArrayT.h"

namespace Tahoe {

/* forward declarations */
class dArray2DT;

/** array class to facilitate working with subsets of a dArray2DT. The
 * principal function of this class is to gather data from specified rows
 * of a dArray2DT using LocalArrayT::SetLocal. The source dArray2DT must 
 * first be set using LocalArrayT::SetLocal. The data is transposed in the
 * process of being gathered, meaning that the nth row of a LocalArrayT
 * contains the data from the nth column of the source dArray2DT. The
 * data in the LocalArrayT is stored in row-major ordering. */
class LocalArrayT: public dArrayT
{	
public:

	/** array data types. Used by the node and element classes to 
	 * resolve the source for the array data. */
	enum TypeT {kUnspecified, /**< unspecified data type */
                       kDisp, /**< displacements */
                        kVel, /**< velocities */
                        kAcc, /**< accelerations */
                   kLastDisp, /**< displacements from the previous time step */
                    kLastVel, /**< velocities from the previous time step */
                    kLastAcc, /**< accelerations from the previous time step */
	             kInitCoords, /**< initial coordinates */
                 kCurrCoords  /**< current coordinates */ };

	/** \name constructors */
	/*@{*/
	/** default constructor. Constructs an array of zero length with type
	 * LocalArrayT::kUnspecified. */
	LocalArrayT(void);

	/** constructor. Constructs an array of zero length with the 
	 * specified type. */
	explicit LocalArrayT(TypeT type);

	/** constructor */
	LocalArrayT(TypeT type, int numnodes, int minordim);

	/** copy constructor */
	LocalArrayT(const LocalArrayT& source);
	/*@}*/

	/** dimension array */
	void Dimension(int numnodes, int minordim);

	/** \deprecated replaced by LocalArrayT::Dimension on 02/13/2002 */
	void Allocate(int numnodes, int minordim);
 
	/** \deprecated replaced by LocalArrayT::Alias */
	void Set(int numnodes, int minordim, double* p);

	/** create a shallow array */
	void Alias(int numnodes, int minordim, const double* p);
	
	/** copy data from an nArrayT.
	 * \param numnodes set LocalArrayT::NumberOfNodes
	 * \param minordim set LocalArrayT::MinorDim
	 * \param source data source array */
	void Copy(int numnodes, int minordim, const nArrayT<double>& source);

	/** set the array type */
	void SetType(TypeT type);
		
	/** array type */
	TypeT Type(void) const;

	/** major dimension of the array */
	int NumberOfNodes(void) const;

	/** minor dimension of the array */
	int MinorDim(void) const;
	
	/** element accessor */
	double& operator()(int majordim, int minordim);

	/** element accessor */
	const double& operator()(int majordim, int minordim) const;

	/** return a pointer to the specified row */
	double* operator()(int minordim);

	/** return a pointer to the specified row */
	const double* operator()(int minordim) const;

	/** assignment operator. This operator re-dimensions the array as
	 * needed and sets the source dArray2DT. The type of this array
	 * is not changed. */
	LocalArrayT& operator=(const LocalArrayT& RHS);

	/** assignment operator. Set all elements in the array to value */
	LocalArrayT& operator=(const double value);

	/** create an alias to the source array */
	void Alias(const LocalArrayT& source);

	/** return 1 if a source matrix has been set with LocalArrayT::SetGlobal,
	 * 0 otherwise. */
	int IsRegistered(void) const;

	/** set the source dArray2DT.
	 * \param global source for data collected with LocalArrayT::SetLocal.
	 *        The source must have the same MinorDim as this array. */
	void SetGlobal(const dArray2DT& global);

	/** gather values from the source dArray2DT. Values are collected from
	 * the rows of the source dArray2DT, being transposed in the process.
	 * \param keys list rows to gather. The length of keys must be the same
	 *        as the NumberOfNodes of this array. */
	void SetLocal(const ArrayT<int>& keys);
	
	/** return the vector with transposed indexing */
	void ReturnTranspose(nArrayT<double>& transpose) const;
	
	/** construct array with local ordering by transposing source */
	void FromTranspose(const nArrayT<double>& transpose);

	/** construct array with local ordering by transposing source */
	void FromTranspose(const double* transpose);

	/** scale and accumulate in local ordering by transposing source */
	void AddScaledTranspose(double scale, const nArrayT<double>& transpose);

	/** combining arrays - inserts all of source at start_node */
	void BlockCopyAt(const LocalArrayT& source, int start_node);

	/** collect subset 
	 * \param nodes indicies of nodes within the source to collect 
	 * \param source source array */
	void Collect(const ArrayT<int>& nodes, const LocalArrayT& source);

	/** compute the array average value
	 * \param avg returns with the average value: [minor dim] */
	void Average(dArrayT& avg) const;

private:

	/* parameters */
	TypeT fType;     /**< type designator for the data on the array */
	int   fNumNodes; /**< major dimension */
	int   fMinorDim; /**< minor dimension */
	
	/** source for LocalArrayT::SetLocal */
	const dArray2DT* fGlobal;
};

/* in-lines */

/* allocating */
inline void LocalArrayT::Dimension(int numnodes, int minordim)
{
	fNumNodes = numnodes;
	fMinorDim = minordim;
	
	/* call single argument allocate function */
	dArrayT::Dimension(fNumNodes*fMinorDim);
}

/* \deprecated replaced by LocalArrayT::Dimension on 02/13/2002 */
inline void LocalArrayT::Allocate(int numnodes, int minordim)
{
  Dimension(numnodes, minordim); 
}

inline void LocalArrayT::Set(int numnodes, int minordim, double *p)
{
	fNumNodes = numnodes;
	fMinorDim = minordim;
	
	/* inherited */
	dArrayT::Set(fNumNodes*fMinorDim, p);
}

inline void LocalArrayT::Alias(int numnodes, int minordim, const double *p)
{
	fNumNodes = numnodes;
	fMinorDim = minordim;
	
	/* inherited */
	dArrayT::Alias(fNumNodes*fMinorDim, p);
}

inline void LocalArrayT::SetType(TypeT type) { fType = type; }

/* accessors */
inline LocalArrayT::TypeT LocalArrayT::Type(void) const    { return fType;     }
inline int LocalArrayT::NumberOfNodes(void) const { return fNumNodes; }
inline int LocalArrayT::MinorDim(void) const      { return fMinorDim; }

/* element accessors */
inline double& LocalArrayT::operator()(int majordim, int minordim)
{
#if __option (extended_errorcheck)
	if (majordim < 0 ||
	    majordim >= fNumNodes ||
	    minordim < 0 ||
	    minordim >= fMinorDim) ExceptionT::OutOfRange("LocalArrayT::operator(,)");
#endif

	return (fArray[minordim*fNumNodes + majordim]);
}
inline const double& LocalArrayT::operator()(int majordim, int minordim) const
{
#if __option (extended_errorcheck)
	if (majordim < 0 ||
	    majordim >= fNumNodes ||
	    minordim < 0 ||
	    minordim >= fMinorDim) ExceptionT::OutOfRange("LocalArrayT::operator(,)");
#endif

	return (fArray[minordim*fNumNodes + majordim]);
}

inline double* LocalArrayT::operator()(int minordim)
{
#if __option (extended_errorcheck)
	if (minordim < 0 || minordim >= fMinorDim) ExceptionT::OutOfRange("LocalArrayT::operator()");
#endif

	return (fArray + minordim*fNumNodes);
}
inline const double* LocalArrayT::operator()(int minordim) const
{
#if __option (extended_errorcheck)
	if (minordim < 0 || minordim >= fMinorDim) ExceptionT::OutOfRange("LocalArrayT::operator()");
#endif

	return (fArray + minordim*fNumNodes);
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

inline void LocalArrayT::FromTranspose(const nArrayT<double>& transpose)
{
#if __option (extended_errorcheck)
	/* dimension check */
	if(fLength != transpose.Length()) ExceptionT::SizeMismatch("LocalArrayT::FromTranspose");
#endif

	/* call primitive function */
	FromTranspose(transpose.Pointer());
}

} /* namespace Tahoe */

#endif /* _LOCALARRAY_T_H_ */
