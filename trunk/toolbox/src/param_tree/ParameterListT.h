/* $Id: ParameterListT.h,v 1.24 2005-04-14 01:17:46 paklein Exp $ */
#ifndef _PARAMETER_LIST_T_H_
#define _PARAMETER_LIST_T_H_

/* direct members */
#include "ParameterT.h"
#include "AutoArrayT.h"

namespace Tahoe {

/* forward declarations */
class ParameterInterfaceT;

/** list of parameters.
 * A ParameterListT can contain three types of entries:
 *    -# plain ParameterT's
 *    -# nested lists of ParameterListT's reproduced in entirety
 */
class ParameterListT
{
public:

	/** parameter and list occurrence */
	enum OccurrenceT {
		Undefined,  /**< undefined */
		Once,       /**< exactly once */
		ZeroOrOnce, /**< zero or one time */
		OnePlus,    /**< one or more times */
		Any         /**< zero or any number of times */
	};

	/** ordering of list items */
	enum ListOrderT {
		Sequence,  /**< ordered sequence of sublists */
		Choice     /**< choice of one the sublists */
	};

	/** constructor */
	ParameterListT(const char* name);

	/** default constructor. Needed to allow making lists of lists */
	ParameterListT(void);

	/** \name list name */
	/*@{*/
	const StringT& Name(void) const { return fName; };
	void SetName(const char* name) { fName = name; };
	/*@}*/
	
	/** \name set/get list attributes */
	/*@{*/
	/** the list order. Note that if the ParameterListT::ListOrderT is ParameterListT::Choice,
	 * the ParameterListT::OccurrenceT for each entry is ParameterListT::Once */
	ListOrderT ListOrder(void) const { return fListOrder; };

	/** set/change the list order */
	void SetListOrder(ListOrderT list_order);
	
	bool Inline(void) const { return fInline; };
	
	/** set/change inlining flag. ParameterListT::PlainList's with parameters cannot
	 * be converted to ParameterListT::Group's. */
	void SetInline(bool is_inline);

	/** flag which controls if duplicate names can be added with ParameterListT::AddList */
	bool DuplicateListNames(void) const { return fDuplicateListNames; };

	/** set flag which controls if duplicate names can be added with ParameterListT::AddList.
	 * If false, ParameterListT::AddList will return false if the name has already been added. */
	void SetDuplicateListNames(bool dup) { fDuplicateListNames = dup; };
	/*@}*/
	
	/** \name dimensions */
	/*@{*/
	/** number of parameters */
	int NumParameters(void) const { return fParameters.Length(); };

	/** number of nested parameter lists */
	int NumLists(void) const { return fParameterLists.Length(); };

	/** number of nested parameter lists with the given name */
	int NumLists(const char* name) const;
	/*@}*/

	/** \name adding items to the list */
	/*@{*/
	/** add a parameter. Returns true of there where no conflicts with
	 * existing parameters. The names of parameters cannot be repeated.
	 * By default, the ParameterListT::OccurrenceT is ParameterListT::Once. */
	bool AddParameter(const ParameterT& param, OccurrenceT occur = Once); 

	bool AddParameter(int a, const char* name, OccurrenceT occur = Once);
	bool AddParameter(double x, const char* name, OccurrenceT occur = Once);
	bool AddParameter(const char* s, const char* name, OccurrenceT occur = Once);
	bool AddParameter(bool b, const char* name, OccurrenceT occur = Once);
	bool AddParameter(ValueT::TypeT t, const char* name, OccurrenceT occur = Once);

	/** remove the first instance of a parameter. Returns true if the given parameter 
	 * was found and removed. */
	bool RemoveParameter(const char* name); 

	/** add a parameter list. Returns true if there where no conflicts with
	 * existing parameter lists. The names of parameter lists cannot be repeated.
	 * By default, the ParameterListT::OccurrenceT is ParameterListT::Once. */
	bool AddList(const ParameterListT& param_list, OccurrenceT occur = Once); 

	/** remove the first instance of a parameter list. Returns true if the given list 
	 * was found and removed. */
	bool RemoveList(const char* name); 
	/*@}*/

	/** \name access to the list entries and occurrences */
	/*@{*/
	const ArrayT<ParameterT>&                  Parameters(void) const { return fParameters; };
	const ArrayT<ParameterListT::OccurrenceT>& ParameterOccurrences(void) const { return fParametersOccur; };
	const ArrayT<ParameterListT>&              Lists(void) const { return fParameterLists; };
	const ArrayT<ParameterListT::OccurrenceT>& ListOccurrences(void) const { return fParameterListsOccur; };
	/*@}*/

	/** \name query access to parameters and sublists. Methods return NULL if the request
	 * cannot be completed. */
	/*@{*/		
	/** return the pointer to the given list. Returns a pointer to the nth instance of the
	 * given list or NULL if the list is not found or the instance is out of range. */
	const ParameterListT* List(const char* name, int instance = 0) const;

	/** return the non-const pointer to the given list. Returns a pointer to the nth 
	 * instance of the given list or NULL if the list is not found or the instance is 
	 * out of range. */
	ParameterListT* List(const char* name, int instance = 0);

	/** search for list by name. Returns a pointer to the nth list whose name contains
	 * the given search string or NULL if the list is not found or the instance is out 
	 * of range. */
	const ParameterListT* FindList(const char* search_name, int instance = 0) const;

	/** return the list associated a choice. Returns a pointer to the nth instance of the
	 * given list resulting from a choice or NULL if the list is not found or the instance 
	 * is out of range. 
	 * \param source ParameterInterfaceT which defined the choice
	 * \param choice_name name of the choice defined by ParameterInterfaceT::DefineInlineSub 
	 * \param instance occurrence within the list */
	const ParameterListT* ListChoice(const ParameterInterfaceT& source, const char* choice_name, int instance = 0) const;

	/** return the index to the given list. Returns index of the nth instance of the
	 * given list or -1 if the list is not found or the instance is out of range. */
	int ListIndex(const char* name, int instance = 0) const;

	/** return the pointer to the given parameter or NULL if the list is not found */
	const ParameterT* Parameter(const char* name) const;

	/** return the non-const pointer to the given parameter or NULL if the list is not found */
	ParameterT* Parameter(const char* name);

	/*@}*/

	/** access to a specific occurrence of a list. Method throws an exception of the
	 * request cannot be completed */
	const ParameterListT& GetList(const char* name, int instance = 0) const;

	/** access to a specific occurrence of a list. Method throws an exception of the
	 * request cannot be completed */
	const ParameterListT& GetList(int instance) const;

	/** \name retrieving parameter values 
	 * Methods throw ExceptionT::kGeneralFail if the parameter is not found. */
	/*@{*/
	void GetParameter(const char* name, int& a) const;
	void GetParameter(const char* name, double& a) const;
	void GetParameter(const char* name, StringT& a) const;
	void GetParameter(const char* name, bool& a) const;

	/** return the given parameter. Throws an exception of the parameter is not present */
	const ParameterT& GetParameter(const char* name) const;

	/** return the given parameter. Throws an exception of the parameter is not present */
	ParameterT& GetParameter(const char* name);

	/** return the list associated a choice. Throws an exception of the parameter is not present */
	const ParameterListT& GetListChoice(const ParameterInterfaceT& source, const char* choice_name, int instance = 0) const;

	/** search for parameter by name. Returns a pointer to the nth parameter whose name contains
	 * the given search string or NULL if the parameter is not found or the instance is out 
	 * of range. */
	const ParameterT* FindParameter(const char* search_name, int instance = 0) const;
	/*@}*/	

	/** \name description */
	/*@{*/
	void SetDescription(const char* description) { fDescription = description; };
	const StringT& Description(void) const { return fDescription; };
	/*@}*/

	/** clear all lists */
	void Clear(void);

protected:

	/** list name */
	StringT fName;
	
	/** list order */
	ListOrderT fListOrder;

	/** flag indicating if list is inline */
	bool fInline;
	
	/** flag indicating whether sublists may have duplicate names */
	bool fDuplicateListNames;

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
};

inline bool ParameterListT::AddParameter(int a, const char* name, OccurrenceT occur)
{
	ParameterT parameter(a, name);
	return AddParameter(parameter, occur);
}
inline bool ParameterListT::AddParameter(double x, const char* name, OccurrenceT occur)
{
	ParameterT parameter(x, name);
	return AddParameter(parameter, occur);
}
inline bool ParameterListT::AddParameter(const char* s, const char* name, OccurrenceT occur)
{
	ParameterT parameter(s, name);
	return AddParameter(parameter, occur);
}
inline bool ParameterListT::AddParameter(bool b, const char* name, OccurrenceT occur)
{
	ParameterT parameter(b, name);
	return AddParameter(parameter, occur);
}
inline bool ParameterListT::AddParameter(ValueT::TypeT t, const char* name, OccurrenceT occur)
{
	ParameterT parameter(t, name);
	return AddParameter(parameter, occur);
}

inline int ParameterListT::ListIndex(const char* name, int instance) const
{
	/* search list */
	int count = 0;
	for (int i = 0; i < fParameterLists.Length(); i++)
		if (fParameterLists[i].Name() == name)
			if (count++ == instance) 
				return i;

	/* fail */
	return -1;
}

inline const ParameterListT* ParameterListT::List(const char* name, int instance) const
{
	int dex = ListIndex(name, instance);
	return (dex > -1) ? fParameterLists.Pointer(dex) : NULL;
}

inline ParameterListT* ParameterListT::List(const char* name, int instance)
{
	/* const this */
	const ParameterListT* this_ = (const ParameterListT*) this;
	const ParameterListT* list = this_->List(name, instance);
	return (ParameterListT*) list;
}

inline const ParameterListT& ParameterListT::GetList(const char* name, int instance) const
{
	const ParameterListT* list = List(name, instance);
	if (!list) ExceptionT::GeneralFail("ParameterListT::GetList", "occurrence %d of list \"%s\" not found in \"%s\"",
		instance + 1, name, Name().Pointer());
	return *list;
}

/* access to a specific occurrence of a list. Method throws an exception of the 
 * request cannot be completed */
inline const ParameterListT& ParameterListT::GetList(int instance) const
{
	if (instance < 0 || instance >= fParameterLists.Length())
		ExceptionT::OutOfRange("ParameterListT::GetList", 
			"instance %d is out of range in \"%s\"",
			instance, Name().Pointer());
	return fParameterLists[instance];
}

inline ParameterT* ParameterListT::Parameter(const char* name)
{
	/* const this */
	const ParameterListT* this_ = (const ParameterListT*) this;
	const ParameterT* parameter = this_->Parameter(name);
	return (ParameterT*) parameter;
}

/* return the list associated a choice. Throws an exception of the parameter is not present */
inline const ParameterListT& ParameterListT::GetListChoice(const ParameterInterfaceT& source, const char* choice_name, int instance) const
{
	const ParameterListT* list = ListChoice(source, choice_name, instance);
	if (!list) ExceptionT::GeneralFail("ParameterListT::GetListChoice", "occurrence %d of choice \"%s\" not found in \"%s\"",
		instance + 1, choice_name, Name().Pointer());
	return *list;
}

} /* namespace Tahoe */

#endif /* _PARAMETER_LIST_T_H_ */
