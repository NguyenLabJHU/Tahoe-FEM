/* $Id: expat_ParseT.h,v 1.4 2004-07-22 08:15:03 paklein Exp $ */
#ifndef _EXPAT_PARSE_T_H_
#define _EXPAT_PARSE_T_H_

/* class requires the expat library */
#ifdef __EXPAT__

#include "expat.h"

/* direct members */
#include "AutoArrayT.h"

namespace Tahoe {

/* forward declarations */
class StringT;
class ParameterListT;

/** interface to the expat parser library. More information about expat is
 * available from here: http://expat.sourceforge.net/. The parser uses a
 * makes use of a stack that is a static data member, so only a single
 * instance can be active at any time. */
class expat_ParseT
{
public:

	/** constructor */
	expat_ParseT(void);

	/** parse the given file placing values into the given parameter tree. If
	 * parsing fails, parser will throw ExceptionT::kBadInputValue. 
	 * \param file path to source file
	 * \param params destination for parsed information. All ParameterT's in
	 *        the parameter list will contain the data from the file as
	 *        strings. params will correspond to the root of the document
	 *        being parsed. */
	void Parse(const StringT& file, ParameterListT& params);

private:

	/** \name expat handlers */
	/*@{*/
	static void startElement(void *userData, const char *name, const char **atts);
	static void endElement(void *userData, const char *name);
	/*@}*/

private:

	/** the parameter list for the current call to Parse */
	static AutoArrayT<ParameterListT*> sListStack;

	/** root element of document being parsed */
	static ParameterListT* sRoot;
};

} /* namespace Tahoe */

#endif /* __EXPAT__ */

#endif /* _EXPAT_PARSE_T_H_ */
