/* $Id: DotLine_FormatterT.h,v 1.2 2003-05-04 22:59:53 paklein Exp $ */
#ifndef _DOT_LINE_FORMATTER_T_H_
#define _DOT_LINE_FORMATTER_T_H_

/* base class */
#include "FormatterT.h"

/* direct members */
#include "StringT.h"

namespace Tahoe {

/** formatter which writes the parameter list values in a nicely
 * formatted text style. This formatter does not implemented a 
 * parameter description. */
class DotLine_FormatterT: public FormatterT
{
public:

	/** constructor */
	DotLine_FormatterT(void);
	
	/** \name writing parameters
	 * All methods return true if successful. */	
	/*@{*/
	/** initialize parameter output stream. Stream must be ready to be written to. */
	virtual bool InitParameterFile(ostream& out) const;

	/** close parameter file */
	virtual bool CloseParameterFile(ostream& out) const;
	
	/** write the parameter list */
	virtual bool WriteParameterList(ostream& out, const ParameterListT& list) const;
	/*@}*/

	/** \name writing data descriptions 
	 * DotLine_FormatterT is only for formatting output of values and does
	 * not implement a parameter description */
	/*@{*/
	virtual bool InitDescriptionFile(ostream& out) const;
	virtual bool CloseDescriptionFile(ostream& out) const;
	virtual bool WriteDescription(ostream& out, const ParameterListT& list) const;
	/*@}*/

	/** set tabbing for each level of depth */
	void SetTabWidth(int tab_width) { fTabWidth = (tab_width > 0) ? tab_width : 0; };

private:

	/** return a row of dots which pads the given string */
	const StringT& Dots(const StringT& str) const;

	/** return a tabbing string */
	const StringT& Tab(void) const;

private:

	/** tree path */
	StringT fPath;
	
	/** \name dots */
	/*@{*/
	/** spaces for each level of depth, 0 by default */
	int fTabWidth;
	StringT fTab;
	
	/** dots string */
	StringT fDots;
	/*@}*/

};

} /* namespace Tahoe */

#endif /* _DOT_LINE_FORMATTER_T_H_ */
