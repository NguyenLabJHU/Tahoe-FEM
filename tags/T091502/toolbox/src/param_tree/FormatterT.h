/* $Id: FormatterT.h,v 1.1 2002-09-03 07:04:33 paklein Exp $ */
#ifndef _FORMATTER_T_H_
#define _FORMATTER_T_H_

/* forward declarations */
#include "ios_fwd_decl.h"
class ValueT;

namespace Tahoe {

/** base class for formatting output of ValueT's. There are
 * assumed to be two kinds of output for every format: value
 * and description. Output of the value writes the ValueT
 * with its current definition. Output of the description
 * writes the data description of the ValueT. */
class FormatterT
{
public:

	/** constructor */
	FormatterT(void);
	
	/** write the value */
	virtual void WriteValue(ostream& out, const ValueT& value) const = 0;
	
	/** write the data description */
	virtual void WriteDescription(ostream& out, const ValueT& value) const = 0;
};

} // namespace Tahoe 
#endif /* _FORMATTER_T_H_ */
