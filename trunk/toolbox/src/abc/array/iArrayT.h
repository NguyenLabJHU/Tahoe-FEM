/* $Id: iArrayT.h,v 1.11 2006-07-03 20:18:43 hspark Exp $ */
/* created: paklein (08/10/1996) */
#ifndef _IARRAY_T_H_
#define _IARRAY_T_H_

/* base class */
#include "nArrayT.h"

namespace Tahoe {

/** integer array class. Most functionality is inherited from nArrayT */
class iArrayT: public nArrayT<int>
{
public:

	/** \name constructors */
	/*@{*/
	iArrayT(void);
	explicit iArrayT(int length);
	iArrayT(const iArrayT& source);

	/** create an alias */
	iArrayT(int length, const int* p);
	/*@}*/
	
	/** \name assigment operators 
	 * Cannot provide operator with (int*) on the right-hand side since this is
	 * indistinguishable from an int. */
	/*@{*/
	iArrayT& operator=(const iArrayT& RHS); /**< assignment operator. Redimensions the array too match the source. */
	iArrayT& operator=(const int* pRHS); /**< assignment operator.  Copy as many values as fit */
	iArrayT& operator=(int value);  /**< set all elements in the array to value */
	/*@}*/

	/** flip specified values in the array
	 * \param from values to match
	 * \param to new value for matching entries
	 * \return number of values changed */
	int ChangeValue(int from, int to);

	/** count the instances of value in the array */
	int	Count(int value) const;
	
	/** return 1 of the value is present, 0 otherwise */
	int HasValue(int value) const;
	
	/** find the first occurence of a value in the array
	 * \param value value to match
	 * \param index returns with the index of the first occurence, if found, -1 otherwise
	 * \return 1 if value found, 0 otherwise */
	int HasValue(int value, int& index) const;

#if 0
	/** sort values in ascending order */
	void SortAscending(void);

	/** sort values in ascending order determined by the master array. Both *this
	 * and master are returned in ascending order. */
	void SortAscending(ArrayT<int>& master);

	/** sort values in ascending order determined by the master array. Both *this
	 * and master are returned in ascending order. */
	void SortAscending(ArrayT<double>& master);
#endif

	/** determine union of the given array. This function does allocate a map array
	 * the with a length the range of values in source. Determining the union requires
	 * 3 passes through the source array and 2 passes through the map. 
	 * \return a reference to *this. */
	iArrayT& Union(const nArrayT<int>& source);

	/** determine the union of the given arrays */
	iArrayT& Union(const ArrayT<const nArrayT<int>*>& source);
};

/* inlines */

/* assigment operators */
inline iArrayT& iArrayT::operator=(const iArrayT& RHS)
{
	nArrayT<int>::operator=(RHS);
	return *this;
}

inline iArrayT& iArrayT::operator=(const int* RHS)
{
	nArrayT<int>::operator=(RHS);
	return *this;
}

inline iArrayT& iArrayT::operator=(int value)
{
	nArrayT<int>::operator=(value);
	return *this;
}

#if 0
inline void iArrayT::SortAscending(void)
{
	/* inherited */
	nArrayT<int>::SortAscending();
}

inline void iArrayT::SortAscending(ArrayT<int>& master)
{
	/* inherited */
	nArrayT<int>::SortAscending(master);
}

inline void iArrayT::SortAscending(ArrayT<double>& master)
{
	/* utility */
	Tahoe::SortAscending(master, *this);
}
#endif

} /* namespace Tahoe */

#endif /* _IARRAY_T_H_ */
