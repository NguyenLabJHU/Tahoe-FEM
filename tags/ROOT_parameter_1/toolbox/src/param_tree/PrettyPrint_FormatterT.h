/* $Id: PrettyPrint_FormatterT.h,v 1.1 2003-04-26 02:09:18 paklein Exp $ */
#ifndef _PRETTY_PRINT_FORMATTER_T_H_
#define _PRETTY_PRINT_FORMATTER_T_H_

/* base class */
#include "FormatterT.h"

/* direct members */
#include "StringT.h"

namespace Tahoe {

/** formatter which writes the parameter list values in a nicely
 * formatted text style. This formatter does not implemented a 
 * parameter description. */
class PrettyPrint_FormatterT: public FormatterT
{
public:

	/** constructor */
	PrettyPrint_FormatterT(void);
	
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
	 * PrettyPrint_FormatterT is only for formatting output of values and does
	 * not implement a parameter description */
	/*@{*/
	virtual bool InitDescriptionFile(ostream& out) const;
	virtual bool CloseDescriptionFile(ostream& out) const;
	virtual bool WriteDescription(ostream& out, const ParameterListT& list) const;
	/*@}*/

private:

	/** return a row of dots which pads the given string */
	const StringT& Dots(const StringT& str) const;
	
private:

	/** \name tree paths */
	/*@{*/
	StringT fPath;
	StringT fLastPath;
	/*@}*/
	
	/** \name dots */
	/*@{*/
	StringT fDots;
	/*@}*/

};

} /* namespace Tahoe */

#endif /* _PRETTY_PRINT_FORMATTER_T_H_ */
