/* $Id: MessageT.h,v 1.1 2002-12-05 08:25:19 paklein Exp $ */
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
	MessageT(CommunicatorT& comm);

protected:

	/** return true if all values are the same */
	bool Same(const ArrayT<int>& a) const;

protected:
 
 	/** the communicator */
 	CommunicatorT& fComm;
 	
 	/** data type for the message */
 	TypeT fType;
};

} /* namespace Tahoe */

#endif /* _MESSAGE_T_H_ */
