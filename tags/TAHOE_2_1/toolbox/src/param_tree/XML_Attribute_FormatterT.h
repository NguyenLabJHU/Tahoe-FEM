/* $Id: XML_Attribute_FormatterT.h,v 1.8 2004-07-20 06:13:21 paklein Exp $ */
#ifndef _XML_ATTRIBUTE_FORMATTER_T_H_
#define _XML_ATTRIBUTE_FORMATTER_T_H_

/* base class */
#include "FormatterT.h"

/* direct members */
#include "BinaryTreeT.h"
#include "StringT.h"

namespace Tahoe {

/* forward declarations */
template <class TYPE> class BinaryTreeT;

/** formatter for XML for which parameter lists are mapped to
 * tags and single parameters are mapped to attributes */
class XML_Attribute_FormatterT: public FormatterT
{
public:

	/** enum for different document type defintions */
	enum DocTypeT {
  Undefined =-1,
		DTD = 0, /**< Document Type Definition */
		XSD = 1  /**< XML Schema */
	};

	/** constructor */
	XML_Attribute_FormatterT(DocTypeT doc_type);
	
	/** set the document description */
	void SetDocDescription(const StringT& doc_root, const StringT& description_path);
	
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

	/** \name accessors to counts 
	 * counts are reset by calls to XML_Attribute_FormatterT::InitDescriptionFile and
	 * accumulated during calls to XML_Attribute_FormatterT::WriteDescription */
	/*@{*/
	int ElementCount(void) const { return fElementCount; };
	int AttributeCount(void) const { return fAttributeCount; };
	int LimitCount(void) const { return fLimitCount; };
	/*@}*/

private:

	/** write the data description. A list of tags is maintained in order to check
	 * for uniqueness in tags. ParameterListT with repeated names will not be processed.
	 * \param out output stream for description
	 * \param list parameter list being described
	 * \param tags list of tags needed to check for uniqueness of tags
	 * \return true if successful, false if problems occur, such as repeated tags */
	bool DoWriteDTD(ostream& out, const ParameterListT& list, BinaryTreeT<StringT>& tags) const;

	/** write the XML schema data description. A list of tags is maintained in order to check
	 * for uniqueness in tags. ParameterListT with repeated names will not be processed.
	 * \param out output stream for description
	 * \param list parameter list being described
	 * \param tags list of tags needed to check for uniqueness of tags
	 * \return true if successful, false if problems occur, such as repeated tags */
	bool DoWriteXSD(ostream& out, const ParameterListT& list, BinaryTreeT<StringT>& tags) const;

	/** \name helpful functions for formatting */
	/*@{*/
	/** return the length of the longest parameter name */
	int ParameterWidth(const ParameterListT& list) const;

	/** return the length of the longest list name */
	int ListWidth(const ParameterListT& list) const;

	/** return the length of the longest reference name */
	//int ReferenceWidth(const ParameterListT& list) const;
	/*@}*/

private:
	
	/** \name document type parameters */
	/*@{*/
	DocTypeT fDocType;
	StringT  fPath;
	StringT  fDocumentRoot;
	/*@}*/

	/** \name counts
	 * Counts reset with call to XML_Attribute_FormatterT::InitDescriptionFile */
	/*@{*/
	int fElementCount;
	int fAttributeCount;
	int fLimitCount;
	/*@}*/
	
	/** collected list of tags which have already been described */
	BinaryTreeT<StringT> fTags;
};

} // namespace Tahoe 
#endif /* _XML_ATTRIBUTE_FORMATTER_T_H_ */
