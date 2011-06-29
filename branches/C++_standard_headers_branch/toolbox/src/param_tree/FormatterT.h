/* $Id: FormatterT.h,v 1.5 2003-04-26 02:06:48 paklein Exp $ */
#ifndef _FORMATTER_T_H_
#define _FORMATTER_T_H_

/* forward declarations */
#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ParameterListT;

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

	/** destructor */
	virtual ~FormatterT(void) {};
	
	/** \name writing parameters
	 * All methods return true if successful. */	
	/*@{*/
	/** initialize parameter output stream. Stream must be ready to be written to. */
	virtual bool InitParameterFile(ostream& out) const = 0;

	/** close parameter file */
	virtual bool CloseParameterFile(ostream& out) const = 0;
	
	/** write the parameter list */
	virtual bool WriteParameterList(ostream& out, const ParameterListT& list) const = 0;
	/*@}*/

	/** \name writing data descriptions 
	 * All methods return true if successful. */
	/*@{*/
	/** initialize description output stream. Stream must be ready to be written to. */
	virtual bool InitDescriptionFile(ostream& out) const = 0;

	/** close description file */
	virtual bool CloseDescriptionFile(ostream& out) const = 0;
	
	/** write the data description */
	virtual bool WriteDescription(ostream& out, const ParameterListT& list) const = 0;
	/*@}*/

protected:

	/** \name tabbing */
	/*@{*/
	/** return the tab depth */
	int Depth(void) const { return fTabCount; };

	/** the tab string */
	const char* Tab(void) const { return fTabs; };
	
	/** increases the tab level by one and returns tab string. \note Function
	 * is not really const, but the tabbing info should not be used directly
	 * by sub-classes. */
	const char* TabOut(void) const;

	/** decreases the tab level by one and returns tab string. \note Function
	 * is not really const, but the tabbing info should not be used directly
	 * by sub-classes. */
	const char* TabIn(void) const;
	/*@}*/

private:

	/** \name tabbing */
	/*@{*/
	/** string used to generate current level of tabbing */
	char fTabs[11];

	/** current tab depth */
	int fTabCount;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _FORMATTER_T_H_ */
