/* $Id: ParameterInterfaceT.h,v 1.5 2003-08-14 01:22:03 paklein Exp $ */
#ifndef _PARAMETER_INTERFACE_T_H_
#define _PARAMETER_INTERFACE_T_H_

/* direct members */
#include "ParameterListT.h"

namespace Tahoe {

/* forward declarations */
class StringT;
class SubListT;

/** abstract interface for classes which define and use parameters. There are
 * two types of parameters accessible through the interface:
 * -# parameters defined using ParameterInterfaceT::DefineParameters
 * -# subordinate parameters lists that are returned by ParameterInterfaceT::SubNames and may either be
 *    -# associated with a subordinate ParameterInterfaceT that must be returned by ParameterInterfaceT::NewSub
 *    -# defined as "inline". Inlined subordinate list do not contain parameters of the first kind and
 *       the interface that has defined the list as inline must define subordinates in the list with
 *       ParameterInterfaceT::DefineInlineSub.
 **/
class ParameterInterfaceT
{
public:

	/** constructor */
	ParameterInterfaceT(const StringT& name);

	/** destructor */
	virtual ~ParameterInterfaceT(void) {};

	/** \name identifier */
	/*@{*/
	const StringT& Name(void) const { return fName; };
	void SetName(const StringT& name);
	/*@}*/

	/** describe the parameters needed by the interface.
	 * \param list destination for the parameter descriptions. The list should have the
	 *        name corresponding to ParameterInterfaceT::Name. */
	virtual void DefineParameters(ParameterListT& list) const;

	/** extract validated parameters. Take a raw list of parameters and produce 
	 * a validated parameter list. If the validated list cannot be 
	 * produced for any reason, the class throws a ExceptionT::kBadInputValue 
	 * \param raw_list source list in which all parameters are stored as
	 *        strings, as read from a source file. 
	 * \param valid_list returns as a validated list witb values of the appropriate data 
	 *        type, validating against constraints and applying any unspecified default values. */
	virtual void ValidateParameterList(const ParameterListT& raw_list, ParameterListT& valid_list) const;

	/** accept parameter list.
	 * \param list input parameter list, which should be validated using ParameterInterfaceT::ValidateParameterList
	 *        to ensure the list conforms to the description defined by the interface. */
	virtual void TakeParameterList(const ParameterListT& list);

	/** \name subordinates that define parameters lists
	 * There are two types of subordinate lists, inlined and non-inlined. Non-inlined
	 * subordinates are returned with a call to ParameterInterfaceT::NewSub. Inlined
	 * subordinates are defined by the call to ParameterInterfaceT::DefineInlineSub. */
	/*@{*/
	/** the order of subordinate lists */
	virtual ParameterListT::ListOrderT ListOrder(void) const;
	
	/** information about subordinate parameter lists
	 * \param sub_lists description of subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list.
	 * Method will be called for each subordinate defined as inline by ParameterInterfaceT::SubNames
	 * or defined recursively by ParameterInterfaceT::DefineInlineSub. 
	 * \param sub name of the inlined subordinate list
	 * \param order defines whether list is a sequence or choice.
	 * \param sub_sub_lists description of contents of this sub list */
	virtual void DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
		SubListT& sub_sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate
	 * or NULL if the name is invalid. Responsibility for deleteting instantiations
	 * resides with the client who requested them. */
	virtual ParameterInterfaceT* NewSub(const StringT& list_name) const;
	/*@}*/

private:

	/** identifier */
	StringT fName;
};

/* forward declarations */
class SubListDescriptionT;

/** sub list descriptions */
class SubListT: public AutoArrayT<SubListDescriptionT>
{
public:

	/** constructor */
	SubListT(void) {};

	//TEMP
	~SubListT(void);

	/** \name add a sublist */
	/*@{*/
	void AddSub(const StringT& name, 
		ParameterListT::OccurrenceT occur = ParameterListT::Once, 
		bool is_inline = false); 
	void AddSub(const SubListDescriptionT& sub);
	/*@}*/
};

/** container for information about sublists */
class SubListDescriptionT
{
public:

	/** \name constructors */
	/*@{*/
	SubListDescriptionT(void) {};
	SubListDescriptionT(const StringT& name, ParameterListT::OccurrenceT occur = ParameterListT::Once, 
		bool is_inline = false);
	SubListDescriptionT(const SubListDescriptionT& source);
	/*@}*/

	//TEMP
	~SubListDescriptionT(void);

	/** \name read/write accessors */
	/*@{*/
	const StringT& Name(void) const { return fName; };
	const ParameterListT::OccurrenceT& Occurrence(void) const { return fOccurrence; };
	const bool& IsInline(void) const { return fIsInline; };
	/*@}*/

	/** assignment operator */
	SubListDescriptionT& operator=(const SubListDescriptionT& rhs);

private:

	/** \name list description */
	/*@{*/
	StringT fName;
	ParameterListT::OccurrenceT fOccurrence;
	bool fIsInline;
	/*@}*/
};

inline SubListDescriptionT::SubListDescriptionT(const StringT& name, 
	ParameterListT::OccurrenceT occur, bool is_inline):
	fName(name),
	fOccurrence(occur),
	fIsInline(is_inline)
{

}

/* add a sublist */
inline void SubListT::AddSub(const StringT& name, ParameterListT::OccurrenceT occur, bool is_inline)
{
	Append(SubListDescriptionT(name, occur, is_inline));
}

inline void SubListT::AddSub(const SubListDescriptionT& sub)
{
	Append(sub);
}

} /* namespace Tahoe */

#endif /* _PARAMETER_SOURCE_T_H_ */
