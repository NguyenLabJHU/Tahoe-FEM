/* $Id: ParameterListT.h,v 1.7 2003-04-26 02:07:36 paklein Exp $ */
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
		ZeroOrOnce, /**< zero or one time */
		OnePlus,    /**< one or more times */
		Any         /**< zero or any number of times */
	};

	/** constructor */
	ParameterListT(const StringT& name): fName(name) { };

	/** default constructor. Needed to allow making lists of lists */
	ParameterListT(void) {};
	
	/** \name list name */
	/*@{*/
	const StringT& Name(void) const { return fName; };
	void SetName(const StringT& name) { fName = name; };
	/*@}*/
	
	/** \name dimensions */
	/*@{*/
	/** number of parameters */
	int NumParameters(void) const { return fParameters.Length(); };

	/** number of nested parameter lists */
	int NumLists(void) const { return fParameterLists.Length(); };

	/** number of references to parameter lists */
	int NumReferences(void) const { return fReferences.Length(); };
	/*@}*/

	/** \name adding items to the list */
	/*@{*/
	/** add a parameter. Returns true of there where no conflicts with
	 * existing parameters. The names of parameters cannot be repeated.
	 * By default, the ParameterListT::OccurrenceT is ParameterListT::Once. */
	bool AddParameter(const ParameterT& param, OccurrenceT occur = Once); 

	bool AddParameter(int a, const StringT& name, OccurrenceT occur = Once);
	bool AddParameter(double x, const StringT& name, OccurrenceT occur = Once);
	bool AddParameter(const StringT& s, const StringT& name, OccurrenceT occur = Once);

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

	/** return the pointer to the given list or NULL if the list is not found */
	const ParameterListT* List(const StringT& name) const;

	/** return the non-const pointer to the given list or NULL if the list is not found */
	ParameterListT* List(const StringT& name);

	/** return the pointer to the given parameter or NULL if the list is not found */
	const ParameterT* Parameter(const StringT& name) const;

	/** return the non-const pointer to the given parameter or NULL if the list is not found */
	ParameterT* Parameter(const StringT& name);
	/*@}*/

	/** \name retrieving parameter values 
	 * Methods throw ExceptionT::kGeneralFail if the parameter is not found. */
	/*@{*/
	void GetParameter(const StringT& name, int& a) const;
	void GetParameter(const StringT& name, double& a) const;
	void GetParameter(const StringT& name, StringT& a) const;
	/*@}*/	

	/** \name description */
	/*@{*/
	void SetDescription(const StringT& description) { fDescription = description; };
	const StringT& Description(void) const { return fDescription; };
	/*@}*/
	
	/** create a validated parameter list. Take a raw list of parameters and a parameter 
	 * description and produce a validated parameter list. If the validated list cannot be 
	 * produced for any reason, the class throws a ExceptionT::kBadInputValue 
	 * \param source raw source list in which all parameters are stored as
	 *        strings as read from a file. 
	 * \param description parameter description list which is used to translate
	 *        values from the source to the appropriate data type, validating
	 *        against constraints and applying any unspecified default values. */
	void Validate(const ParameterListT& source, const ParameterListT& description);

private:

	/** clear all lists */
	void Clear(void);

protected:

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

inline bool ParameterListT::AddParameter(int a, const StringT& name, OccurrenceT occur)
{
	ParameterT parameter(a, name);
	return AddParameter(parameter, occur);
}
inline bool ParameterListT::AddParameter(double x, const StringT& name, OccurrenceT occur)
{
	ParameterT parameter(x, name);
	return AddParameter(parameter, occur);
}
inline bool ParameterListT::AddParameter(const StringT& s, const StringT& name, OccurrenceT occur)
{
	ParameterT parameter(s, name);
	return AddParameter(parameter, occur);
}

inline ParameterListT* ParameterListT::List(const StringT& name)
{
	/* const this */
	const ParameterListT* const this_ = (const ParameterListT* const) this;
	const ParameterListT* list = this_->List(name);
	return (ParameterListT*) list;
}

inline ParameterT* ParameterListT::Parameter(const StringT& name)
{
	/* const this */
	const ParameterListT* const this_ = (const ParameterListT* const) this;
	const ParameterT* parameter = this_->Parameter(name);
	return (ParameterT*) parameter;
}

} /* namespace Tahoe */

#endif /* _PARAMETER_LIST_T_H_ */
