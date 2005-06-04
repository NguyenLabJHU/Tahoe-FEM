/* $Id: AllGatherT.h,v 1.3 2005-06-04 16:59:42 paklein Exp $ */
#ifndef _GATHER_T_H_
#define _GATHER_T_H_

/* base class */
#include "MessageT.h"

/* direct members */
#include "iArrayT.h"

namespace Tahoe {

/** record for gathering values from all processes to all. Each
 * process sends the same data to all other processes. The data
 * from each process does not need to be the same length. */
class AllGatherT: public MessageT
{
public:

	/** constructor */
	AllGatherT(CommunicatorT& comm, int tag);

	/** \name initialize gather data
	 * \param my_data array of data sent from this processor to all. The
	 *        array sets the data type and size for the communication, so
	 *        it must be dimensioned. */
	/*@{*/
	void Initialize(const ArrayT<double>& my_data);
	void Initialize(const ArrayT<int>& my_data);
	void Initialize(int my_size);
	/*@}*/

	/** the size of the gathered array */
	int Total(void) const { return fTotal; };

	/** return true if all the amount of data sent from each processor
	 * is the same. */
	bool Equal(void) const { return fEqual; };

	/** the array of counts from every processor */
	const iArrayT& Counts(void) const { return fCounts; };

	/** \name perform the gather 
	 * Gather does not need to be called with the same data type used
	 * in call to AllGatherT::Initialize. When called with only one
	 * argument, the data from this partition must already be in the
	 * appropriate place within the destination array.
	 * \param my_data array of data sent from this processor to all
	 * \param gather destination for gathered data. Array must be at
	 *        least of length AllGatherT:Total.
	 */
	/*@{*/
	void AllGather(const nArrayT<double>& my_data, nArrayT<double>& gather);
	void AllGather(nArrayT<double>& gather);
	void AllGather(const nArrayT<int>& my_data, nArrayT<int>& gather);
	void AllGather(nArrayT<int>& gather);
	/*@}*/
		
private:

	/** true if all processes are same size */
	bool fEqual;

	/** counts from all */
	iArrayT fCounts;

	/** sum of values in AllGatherT::fCounts */
	int fTotal;

	/** displacements (in the receive buffer) if the data sizes are not equal.
	 * Displacement are calculated to pack the resulting data continuously
	 * in the gather buffer. */
	iArrayT fDisplacements;
};

/* inlines */
inline void AllGatherT::Initialize(const ArrayT<double>& my_data)
{
	Initialize(my_data.Length());
}

inline void AllGatherT::Initialize(const ArrayT<int>& my_data)
{
	Initialize(my_data.Length());
}

} /* namespace Tahoe */

#endif /* _GATHER_T_H_ */
