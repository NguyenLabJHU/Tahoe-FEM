/* $Id: FormatterT.h,v 1.2 2002-09-22 23:05:29 paklein Exp $ */
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
	
	/** \name writing values */
	/*@{*/
	/** initialize value file */
	virtual void InitValueFile(ostream& out, const StringT path) const = 0;

	/** close value file */
	virtual void CloseValueFile(ostream& out) const = 0;
	
	/** write the value */
	virtual void WriteValue(ostream& out, const ValueT& value) const = 0;
	/*@}*/

	/** \name writing data descriptions */
	/*@{*/
	/** initialize description file */
	virtual void InitDescriptionFile(ostream& out, const StringT path) const = 0;

	/** close description file */
	virtual void CloseDescriptionFile(ostream& out) const = 0;
	
	/** write the data description */
	virtual void WriteDescription(ostream& out, const ValueT& value) const = 0;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _FORMATTER_T_H_ */
