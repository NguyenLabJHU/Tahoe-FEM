/* $Id: ParameterListT.h,v 1.4 2002-11-18 09:59:03 paklein Exp $ */
#ifndef _PARAMETER_LIST_T_H_
#define _PARAMETER_LIST_T_H_

/* direct members */
#include "ParameterT.h"
#include "AutoArrayT.h"

namespace Tahoe {

/** list of parameters.
 * A ParameterListT can contain three types of entries:
 *    -# plain ParameterT's
 *    -# nested lists of ParameterListT's reproduced in entirety
 *    -# references to ParameterListT's which are referred to by name and are not reproduced
 */
class ParameterListT
{
public:

	/** parameter occurrence */
	enum OccurrenceT {
		Once,       /**< exactly once */
		ZeroOrOnce, /**< zero or more times */
		OnePlus,    /**< one or more times */
		Any         /**< zero or any number of times */
	};

	/** constructor */
	ParameterListT(const StringT& name): fName(name) { };
	
	/** list name */
	const StringT& Name(void) const { return fName; };
	
	/** \name dimensions */
	/*@{*/
	/** number of parameters */
	int NumParameters(void) const { return fParameters.Length(); };

	/** number of nested parameter lists */
	int NumLists(void) const { return fParameterLists.Length(); };

	/** number of references to parameter lists */
	int NumReferences(void) const { return fReferences.Length(); };
	/*@}*/

#if 0	
	/** \name name space */
	/*@{*/
	/** (re-)set name space name */
	void SetNameSpace(const StringT& ns_name) { fNameSpace = ns_name; };

	/** return the name space name */
	const StringT& NameSpace(void) const { return fNameSpace; };
	/*@}*/
#endif

	/** \name adding items to the list */
	/*@{*/
	/** add a parameter. Returns true of there where no conflicts with
	 * existing parameters. The names of parameters cannot be repeated.
	 * By default, the ParameterListT::OccurrenceT is ParameterListT::Once. */
	bool AddParameter(const ParameterT& param, OccurrenceT occur = Once); 

	/** add a parameter list. Returns true of there where no conflicts with
	 * existing parameter lists. The names of parameter lists cannot be repeated.
	 * By default, the ParameterListT::OccurrenceT is ParameterListT::Once. */
	bool AddList(const ParameterListT& param_list, OccurrenceT occur = Once); 

	/** add a reference. Returns true of there where no conflicts with
	 * existing references. The names of reference cannot be repeated.
	 * By default, the ParameterListT::OccurrenceT is ParameterListT::Once. */
	bool AddReference(const StringT& ref, OccurrenceT occur = Once); 
	/*@}*/

	/** \name access to the list entries and occurrences */
	/*@{*/
	const ArrayT<ParameterT>&                  Parameters(void) const { return fParameters; };
	const ArrayT<ParameterListT::OccurrenceT>& ParameterOccurrences(void) const { return fParametersOccur; };
	const ArrayT<ParameterListT>&              Lists(void) const { return fParameterLists; };
	const ArrayT<ParameterListT::OccurrenceT>& ListOccurrences(void) const { return fParameterListsOccur; };
	const ArrayT<StringT>&                     References(void) const { return fReferences; };
	const ArrayT<ParameterListT::OccurrenceT>& ReferenceOccurrences(void) const { return fReferencesOccur; };
	/*@}*/

	/** \name description */
	/*@{*/
	void SetDescription(const StringT& description) { fDescription = description; };
	const StringT& Description(void) const { return fDescription; };
	/*@}*/

private:

	/** default constructor. Needed to allow making lists of lists */
	ParameterListT(void) {};

protected:

	/** parameters name space */
//	StringT fNameSpace;

	/** list name */
	StringT fName;

	/** description */
	StringT fDescription;

	/** \name simple parameters */
	/*@{*/
	AutoArrayT<ParameterT>  fParameters;
	AutoArrayT<OccurrenceT> fParametersOccur;
	/*@}*/

	/** \name nested parameters lists */
	/*@{*/
	AutoArrayT<ParameterListT> fParameterLists;
	AutoArrayT<OccurrenceT>    fParameterListsOccur;
	/*@}*/

	/** \name references to parameters lists */
	/*@{*/
	AutoArrayT<StringT>     fReferences;
	AutoArrayT<OccurrenceT> fReferencesOccur;
	/*@}*/
};

} /* namespace Tahoe */
#endif /* _PARAMETER_LIST_T_H_ */
