/* $Id: MessageT.cpp,v 1.1.2.1 2002-12-19 03:09:13 paklein Exp $ */
#include "MessageT.h"
#include "ArrayT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
const bool ArrayT<MessageT>::fByteCopy = false;
const bool ArrayT<MessageT*>::fByteCopy = true;
} /* namespace Tahoe */

/* constructor */
MessageT::MessageT(CommunicatorT& comm):
	fComm(comm),
	fType(Void)
{

}

/*************************************************************************
 * Protected
 *************************************************************************/

/* return true if all values are the same */
bool MessageT::Same(const ArrayT<int>& a) const
{
	if (a.Length() == 0)
		return true;
	else
	{
		int i = a[0];
		for (int j = 1; j < a.Length(); j++)
			if (a[j] != i)
				return false;
		return true;
	}
}
