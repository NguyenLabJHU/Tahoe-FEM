/* $Id: XML_Atomic_FormatterT.h,v 1.5 2003-04-26 02:09:46 paklein Exp $ */
#ifndef _XML_ATOMIC_FORMATTER_T_H_
#define _XML_ATOMIC_FORMATTER_T_H_

/* base class */
#include "FormatterT.h"

/* direct members */
#include "StringT.h"

namespace Tahoe {

/* forward declarations */
class ParameterT;

/** formatter for XML for which all values are categorized into the
 * atomic types given by ValueT::TypeT */
class XML_Atomic_FormatterT: public FormatterT
{
public:

	/** constructor */
	XML_Atomic_FormatterT(void);
	
	/** DTD settings */
	void SetDTD(const StringT& doc_root, const StringT& dtd_path);
	
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
	 * All methods return true if successful. */
	/*@{*/
	/** initialize description output stream. Stream must be ready to be written to. */
	virtual bool InitDescriptionFile(ostream& out) const;

	/** close description file */
	virtual bool CloseDescriptionFile(ostream& out) const;

	/** write the data description. Since every ParameterListT is mapped to a
	 * <parameter_list> tag and every ParameterT is mapped to a <parameter> tag,
	 * the contents of list does not affect the data description file. */
	virtual bool WriteDescription(ostream& out, const ParameterListT& list) const;
	/*@}*/

private:

	/** write the value */
	void WriteParameter(ostream& out, const ParameterT& parameter) const;

private:
	
	/** \name DTD parameters */
	/*@{*/
	StringT fDTD;
	StringT fDocumentRoot;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _XML_ATOMIC_FORMATTER_T_H_ */
