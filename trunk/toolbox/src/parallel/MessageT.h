/* $Id: MessageT.h,v 1.2 2003-01-27 06:42:48 paklein Exp $ */
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

	/** destructor */
	virtual ~MessageT(void) { };

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
