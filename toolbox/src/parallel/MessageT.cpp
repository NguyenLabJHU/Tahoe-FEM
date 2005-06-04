/* $Id: MessageT.cpp,v 1.5 2005-06-04 16:59:42 paklein Exp $ */
#include "MessageT.h"
#include "ArrayT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<MessageT>::fByteCopy = false;
DEFINE_TEMPLATE_STATIC const bool ArrayT<MessageT*>::fByteCopy = true;
} /* namespace Tahoe */

/* constructor */
MessageT::MessageT(CommunicatorT& comm, int tag):
	fComm(comm),
	fType(Void),
	fTag(tag)
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
