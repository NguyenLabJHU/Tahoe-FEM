/* $Id: ExtractElement.h,v 1.1 2003-02-18 08:47:23 paklein Exp $ */
#ifndef _EXTRACT_ELEMENT_H_
#define _EXTRACT_ELEMENT_H_

#include "ExtractIOManager.h"

namespace Tahoe {

class ExtractElement: public ExtractIOManager
{
public:
	ExtractElement(ostream& message, istream& in, bool write);
  
protected:

	void Initialize(void);
	void TranslateVariables(void);
	
private:

	StringT fExtract_ID;
};

} /* namespace Tahoe */

#endif /* _EXTRACT_ELEMENT_H_ */
