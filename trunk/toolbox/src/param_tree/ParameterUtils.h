/* $Id: ParameterUtils.h,v 1.4 2004-01-21 17:15:19 paklein Exp $ */
#ifndef _PARAMETER_UTILS_H_
#define _PARAMETER_UTILS_H_

/* direct members */
#include "ParameterInterfaceT.h"

namespace Tahoe {

/** named value */
template <class TYPE>
class NamedParameterT: public ParameterInterfaceT
{
public:

	/** constructor */
	NamedParameterT(const StringT& name);

	/** the value name */
	const StringT& ValueName(void) const { return fValueName; };

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list.
	 * \param list input parameter list, which should be validated using ParameterInterfaceT::ValidateParameterList
	 *        to ensure the list conforms to the description defined by the interface. */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** value */
	TYPE fValue;

private:

	/** name */
	StringT fValueName;
};

template <class TYPE>
NamedParameterT<TYPE>::NamedParameterT(const StringT& name):
	ParameterInterfaceT(name)
{

}

template <class TYPE>
void NamedParameterT<TYPE>::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);
	
	/* value name */
	list.AddParameter(fValueName, "name", ParameterListT::ZeroOrOnce);
	list.AddParameter(fValue, "value");
}

template <class TYPE>
void NamedParameterT<TYPE>::TakeParameterList(const ParameterListT& list)
{
	/* get name */
	const ParameterT* value_name = list.Parameter("name");
	if (value_name)
		fValueName = *value_name;

	/* the value */
	list.GetParameter("value", fValue);
}

/** an integer with default ParameterInterfaceT name "Integer" */
class IntegerParameterT: public NamedParameterT<int>
{
public:

	/** \name constructors */
	/*@{*/
	IntegerParameterT(void);
	IntegerParameterT(const StringT& name);
	/*@{*/
};

/** a double with default ParameterInterfaceT name "Double" */
class DoubleParameterT: public NamedParameterT<double>
{
public:

	/** \name constructors */
	/*@{*/
	DoubleParameterT(void);
	DoubleParameterT(const StringT& name);
	/*@{*/
};

/** a string with default ParameterInterfaceT name "String" */
class StringParameterT: public NamedParameterT<StringT>
{
public:

	/** \name constructors */
	/*@{*/
	StringParameterT(void);
	StringParameterT(const StringT& name);
	/*@{*/
};

/** base class for names lists where TYPE is a type derived from
 * ParameterInterfaceT */
template <class TYPE>
class NamedListT: public ParameterInterfaceT
{
public:

	/** constructors */
	NamedListT(const StringT& name);

	/** the list name */
	const StringT& ListName(void) const { return fListName; };

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);

	/** information about subordinate parameters */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& list_name) const;
	/*@}*/

private:

	/** list name */
	StringT fListName;
};

/* constructors */
template <class TYPE>
NamedListT<TYPE>::NamedListT(const StringT& name):
	ParameterInterfaceT(name)
{

}

/* describe the parameters needed by the interface */
template <class TYPE>
void NamedListT<TYPE>::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	/* add optional name */
	list.AddParameter(fListName, "name", ParameterListT::ZeroOrOnce);
}

/* information about subordinate parameters */
template <class TYPE>
void NamedListT<TYPE>::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);
	
	/* zero or more list entries */
	TYPE list_entry;
	sub_list.AddSub(list_entry.Name(), ParameterListT::Any);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
template <class TYPE>
ParameterInterfaceT* NamedListT<TYPE>::NewSub(const StringT& list_name) const
{
	TYPE list_entry;
	if (list_name == list_entry.Name())
		return new TYPE;
	else
		return ParameterInterfaceT::NewSub(list_name);
}

/* accept parameter list */
template <class TYPE>
void NamedListT<TYPE>::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* get name if defined */
	const ParameterT* list_name = list.Parameter("name");
	if (list_name) 
		fListName = *list_name;
}

/** list of integers. Defines a list of integers with default ParameterInterfaceT name 
 * "IntegerList" which contains zero or more "Integer" entries. */
class IntegerListT: public NamedListT<IntegerParameterT>
{
public:

	/** \name constructors */
	/*@{*/
	IntegerListT(const StringT& name);
	IntegerListT(void);
	/*@}*/
};

/** list of double's. Defines a list of double's with default ParameterInterfaceT name 
 * "DoubleList" which contains zero or more "Double" entries. */
class DoubleListT: public NamedListT<DoubleParameterT>
{
public:

	/** \name constructors */
	/*@{*/
	DoubleListT(const StringT& name);
	DoubleListT(void);
	/*@}*/
};

/** list of string's. Defines a list of string's with default ParameterInterfaceT name 
 * "StringList" which contains zero or more "String" entries. */
class StringListT: public NamedListT<StringParameterT>
{
public:

	/** \name constructors */
	/*@{*/
	StringListT(const StringT& name);
	StringListT(void);
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _PARAMETER_UTILS_H_ */
