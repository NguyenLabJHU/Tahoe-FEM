/* $Id: ParameterContainerT.h,v 1.9 2005-04-05 15:51:11 paklein Exp $ */
#ifndef _PARAMETER_CONTAINER_T_H_
#define _PARAMETER_CONTAINER_T_H_

/* base classes */
#include "ParameterInterfaceT.h"

/* direct members */
#include "ParameterListT.h"

namespace Tahoe {

/** standalone container of parameters. Class that both describes and contains a
 * user-definable list of parameters. Parameters are defined by the user through
 * the ParameterListT and then the class internally takes care of calls to the
 * ParameterInterfaceT interface. The only sublists which the container can
 * construct are those predefined for ParameterInterfaceT. All others will be
 * requested from the ParameterContainerT::fSubSource, if it has been set with
 * ParameterContainerT::SetSubSource. */
class ParameterContainerT: public ParameterInterfaceT
{
public:

	/** \name constructors */
	/*@{*/
	ParameterContainerT(const StringT& name);

	/** \name default constructor needed to define arrays */
	ParameterContainerT(void);
	/*@}*/

	/** \name adding items to the list */
	/*@{*/
	/** add a parameter. Returns true of there where no conflicts with
	 * existing parameters. The names of parameters cannot be repeated.
	 * By default, the ParameterListT::OccurrenceT is ParameterListT::Once. */
	bool AddParameter(const ParameterT& param, ParameterListT::OccurrenceT occur = ParameterListT::Once); 

	bool AddParameter(int a, const char* name, ParameterListT::OccurrenceT occur = ParameterListT::Once);
	bool AddParameter(double x, const char* name, ParameterListT::OccurrenceT occur = ParameterListT::Once);
	bool AddParameter(const char* s, const char* name, ParameterListT::OccurrenceT occur = ParameterListT::Once);
	bool AddParameter(bool b, const char* name, ParameterListT::OccurrenceT occur = ParameterListT::Once);
	bool AddParameter(ValueT::TypeT t, const char* name, ParameterListT::OccurrenceT occur = ParameterListT::Once);
	/*@}*/

	/** \name add a sublist */
	/*@{*/
	void AddSub(const StringT& name, 
		ParameterListT::OccurrenceT occur = ParameterListT::Once, 
		bool is_inline = false); 
	void AddSub(const SubListDescriptionT& sub);
	void AddSub(const ParameterContainerT& sub, ParameterListT::OccurrenceT occur = ParameterListT::Once,
		bool is_inline = false);
	/*@}*/

	/** set source for subs not defined by the container. The source must remain
	 * in scope at least as long as the container since the container may ask the
	 * source to construct any subs it does not define. Pass NULL to clear the
	 * sub source. */
	void SetSubSource(const ParameterInterfaceT* sub_source);

	/** a pointer to the ParameterInterfaceT of the given subordinate. If the container
	 * does not define the given sub, it will attempt to get it from ParameterContainerT::fSubSource,
	 * if the source has been defined. */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** return the description of the given inline subordinate parameter list.
     * If the container does not define the given sub, it will attempt to get 
     * it from ParameterContainerT::fSubSource, if the source has been defined. */
     virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** \name set/get list attributes */
	/*@{*/
	/** the order of subordinate lists */
	virtual ParameterListT::ListOrderT ListOrder(void) const;

	/** set/change the list order */
	void SetListOrder(ParameterListT::ListOrderT list_order);
	/*@}*/

	/** \name description */
	/*@{*/
	void SetDescription(const char* description) { fDescription = description; };
	const StringT& Description(void) const { return fDescription; };
	/*@}*/

protected:

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface. */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists
	 * \param sub_lists description of subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;
	/*@}*/

protected:

	/** list order */
	ParameterListT::ListOrderT fListOrder;

	/** flag indicating if list is inline */
	bool fInline;

	/** description */
	StringT fDescription;

	/** \name simple parameters */
	/*@{*/
	AutoArrayT<ParameterT> fParameters;
	AutoArrayT<ParameterListT::OccurrenceT> fParametersOccur;
	/*@}*/

	/** \name nested lists */
	/*@{*/
	/** sublists registered by ParameterContainerT::AddSub */
	SubListT fSubs;
	
	/** nested ParameterContainerT's registered ParameterContainerT::AddSub */
	AutoArrayT<ParameterContainerT> fContainers;
	AutoArrayT<ParameterListT::OccurrenceT> fContainersOccur;
	AutoArrayT<bool> fContainersInline;
	
	/** source for subs that are not defined by the container */
	const ParameterInterfaceT* fSubSource;
	/*@}*/

	/** \name maintain order of ParameterContainerT::fSubs and ParameterContainerT::fContainers
	 * as registered in the constructed list in ParameterContainerT::DefineSubs */
	/*@{*/
	AutoArrayT<char> fs_or_c;
	AutoArrayT<int> fs_or_c_index;	
	/*@}*/
};

/* inlines */

inline bool ParameterContainerT::AddParameter(int a, const char* name, ParameterListT::OccurrenceT occur)
{
	ParameterT parameter(a, name);
	return AddParameter(parameter, occur);
}
inline bool ParameterContainerT::AddParameter(double x, const char* name, ParameterListT::OccurrenceT occur)
{
	ParameterT parameter(x, name);
	return AddParameter(parameter, occur);
}
inline bool ParameterContainerT::AddParameter(const char* s, const char* name, ParameterListT::OccurrenceT occur)
{
	ParameterT parameter(s, name);
	return AddParameter(parameter, occur);
}
inline bool ParameterContainerT::AddParameter(bool b, const char* name, ParameterListT::OccurrenceT occur)
{
	ParameterT parameter(b, name);
	return AddParameter(parameter, occur);
}
inline bool ParameterContainerT::AddParameter(ValueT::TypeT t, const char* name, ParameterListT::OccurrenceT occur)
{
	ParameterT parameter(t, name);
	return AddParameter(parameter, occur);
}

inline ParameterListT::ListOrderT ParameterContainerT::ListOrder(void) const { return fListOrder; }
inline void ParameterContainerT::SetListOrder(ParameterListT::ListOrderT list_order) { fListOrder = list_order; }
//inline bool ParameterContainerT::Inline(void) const { return fInline; }
//inline void ParameterContainerT::SetInline(bool is_inline) { fInline = is_inline; }

} /* namespace Tahoe */

#endif /* _PARAMETER_CONTAINER_T_H_ */
