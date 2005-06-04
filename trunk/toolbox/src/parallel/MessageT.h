/* $Id: MessageT.h,v 1.3 2005-06-04 16:59:42 paklein Exp $ */
#ifndef _MESSAGE_T_H_
#define _MESSAGE_T_H_

namespace Tahoe {

/* forward declarations */
class CommunicatorT;
template <class TYPE> class ArrayT;

/** base class of records for parallel communication */
class MessageT
{
public:

	/** enumerator for value type */
	enum TypeT {
		Void,
		Integer,
		Double,
		String
	};

	/** constructor */
	MessageT(CommunicatorT& comm, int tag);

	/** destructor */
	virtual ~MessageT(void) {};

	/** message tag */
	int Tag(void) const { return fTag; };

protected:

	/** return true if all values are the same */
	bool Same(const ArrayT<int>& a) const;

protected:
 
 	/** the communicator */
 	CommunicatorT& fComm;
 	
 	/** data type for the message */
 	TypeT fType;

	/** message tag */
	int fTag;
};

} /* namespace Tahoe */

#endif /* _MESSAGE_T_H_ */
