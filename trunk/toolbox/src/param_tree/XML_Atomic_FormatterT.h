/* $Id: XML_Atomic_FormatterT.h,v 1.1 2002-11-16 20:50:21 paklein Exp $ */
#ifndef _XML_ATOMIC_FORMATTER_T_H_
#define _XML_ATOMIC_FORMATTER_T_H_

/* base class */
#include "FormatterT.h"

/* direct members */
#include "StringT.h";

namespace Tahoe {

/** formatter for XML for which all values are categorized into the
 * atomic types given by ValueT::TypeT */
class XML_Atomic_FormatterT: public FormatterT
{
public:

	/** constructor */
	XML_Atomic_FormatterT(void) {};
	
	/** DTD settings */
	void SetDTD(const StringT& dtd_path, const StringT& doc_root);
	
	/** \name writing parameters
	 * All methods return true if successful. */	
	/*@{*/
	/** initialize parameter output stream. Stream must be ready to be written to. */
	virtual bool InitParameterFile(ostream& out) const;

	/** close parameter file */
	virtual bool CloseParameterFile(ostream& out) const;
	
	/** write the value */
//	virtual bool WriteParameter(ostream& out, const ParameterT& parameter) const;
	/*@}*/

	/** \name writing data descriptions 
	 * All methods return true if successful. */
	/*@{*/
	/** initialize description output stream. Stream must be ready to be written to. */
	virtual bool InitDescriptionFile(ostream& out) const;

	/** close description file */
	virtual bool CloseDescriptionFile(ostream& out) const;
	
	/** write the data description */
	virtual bool WriteDescription(ostream& out, const ParameterT& parameter) const;
	/*@}*/

private:
	
	/** \name DTD parameters */
	/*@{*/
	StringT fDTD;
	StringT fDocumentRoot;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _XML_ATOMIC_FORMATTER_T_H_ */
