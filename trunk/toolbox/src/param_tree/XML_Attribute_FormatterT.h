/* $Id: XML_Attribute_FormatterT.h,v 1.4 2003-04-26 02:09:46 paklein Exp $ */
#ifndef _XML_ATTRIBUTE_FORMATTER_T_H_
#define _XML_ATTRIBUTE_FORMATTER_T_H_

/* base class */
#include "FormatterT.h"

/* direct members */
#include "StringT.h"

namespace Tahoe {

/* forward declarations */
template <class TYPE> class BinaryTreeT;

/** formatter for XML for which parameter lists are mapped to
 * tags and single parameters are mapped to attributes */
class XML_Attribute_FormatterT: public FormatterT
{
public:

	/** constructor */
	XML_Attribute_FormatterT(void);
	
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

	/** write the data description */
	virtual bool WriteDescription(ostream& out, const ParameterListT& list) const;
	/*@}*/

private:

	/** write the data description. A list of tags is maintained in order to check
	 * for uniqueness in tags. ParameterListT with repeated names will not be processed.
	 * \param out output stream for description
	 * \param list parameter list being described
	 * \param tags list of tags needed to check for uniqueness of tags
	 * \return true if successful, false if problems occur, such as repeated tags */
	bool DoWriteDescription(ostream& out, const ParameterListT& list, BinaryTreeT<StringT>& tags) const;

	/** \name helpful functions for formatting */
	/*@{*/
	/** return the length of the longest parameter name */
	int ParameterWidth(const ParameterListT& list) const;

	/** return the length of the longest list name */
	int ListWidth(const ParameterListT& list) const;

	/** return the length of the longest reference name */
	int ReferenceWidth(const ParameterListT& list) const;
	/*@}*/

private:
	
	/** \name DTD parameters */
	/*@{*/
	StringT fDTD;
	StringT fDocumentRoot;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _XML_ATTRIBUTE_FORMATTER_T_H_ */
