/* $Id: ParameterUtils.h,v 1.10 2004-07-22 21:06:43 paklein Exp $ */
#ifndef _PARAMETER_UTILS_H_
#define _PARAMETER_UTILS_H_

/* direct members */
#include "ParameterInterfaceT.h"
#include "dArrayT.h"
#include "dMatrixT.h"

namespace Tahoe {

/** named value */
template <ParameterT::TypeT TYPE>
class NamedParameterT: public ParameterInterfaceT
{
public:

	/** constructors */
	NamedParameterT(const StringT& name);

	/** the value name */
	const StringT& ValueName(void) const { return fValueName; };

	/** implicit type conversion */
	operator const ValueT&() const { return fValue; };

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
	ValueT fValue;

private:

	/** name */
	StringT fValueName;
};

template <ParameterT::TypeT TYPE>
NamedParameterT<TYPE>::NamedParameterT(const StringT& name):
	ParameterInterfaceT(name),
	fValue(TYPE)
{

}

template <ParameterT::TypeT TYPE>
void NamedParameterT<TYPE>::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);
	
	/* value name */
	list.AddParameter(ParameterT::Word, "name", ParameterListT::ZeroOrOnce);
	list.AddParameter(TYPE, "value");
}

template <ParameterT::TypeT TYPE>
void NamedParameterT<TYPE>::TakeParameterList(const ParameterListT& list)
{
	/* get name */
	const ParameterT* value_name = list.Parameter("name");
	if (value_name)
		fValueName = *value_name;

	/* the value */
	fValue = list.GetParameter("value");
}

/** an integer with default ParameterInterfaceT name "Integer" */
class IntegerParameterT: public NamedParameterT<ParameterT::Integer>
{
public:

	/** \name constructors */
	/*@{*/
	IntegerParameterT(void);
	IntegerParameterT(const StringT& name);
	/*@{*/
};

/** a double with default ParameterInterfaceT name "Double" */
class DoubleParameterT: public NamedParameterT<ParameterT::Double>
{
public:

	/** \name constructors */
	/*@{*/
	DoubleParameterT(void);
	DoubleParameterT(const StringT& name);
	/*@{*/
};

/** a string (word) with default ParameterInterfaceT name "String" */
class StringParameterT: public NamedParameterT<ParameterT::Word>
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

	/** \name list length limits
	 * Limits can be set/re-set any time before NamedListT::DefineParameters is
	 * called. Pass -1 to clear length limits. \note: length limits are not
	 * fully supported yet. */
	/*@{*/
	void SetMinLength(int min_length) { fMinLength = min_length; };
	void SetMaxLength(int max_length) { fMaxLength = max_length; };
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);

	/** information about subordinate parameters */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	/*@}*/

private:

	/** list name */
	StringT fListName;
	
	/** \name list length limits */
	/*@{*/
	int fMinLength;
	int fMaxLength;
	/*@}*/
};

/* constructors */
template <class TYPE>
NamedListT<TYPE>::NamedListT(const StringT& name):
	ParameterInterfaceT(name),
	fMinLength(-1),
	fMaxLength(-1)
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
	if (fMinLength == 1)
		sub_list.AddSub(list_entry.Name(), ParameterListT::OnePlus);
	else
		sub_list.AddSub(list_entry.Name(), ParameterListT::Any);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
template <class TYPE>
ParameterInterfaceT* NamedListT<TYPE>::NewSub(const StringT& name) const
{
	TYPE list_entry;
	if (name == list_entry.Name())
		return new TYPE;
	else
		return ParameterInterfaceT::NewSub(name);
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

	/** extract string parameters to an array */
	static void Extract(const ParameterListT& list, ArrayT<StringT>& values);
};

/** vector */
class VectorParameterT: public ParameterInterfaceT
{
public:

	/** \name constructors */
	/*@{*/
	VectorParameterT(const StringT& name, int dim, char variable = 'v');
	VectorParameterT(int dim, char variable = 'v');
	
	/** construct extracting length from the name
	 * \param name_N name of the vector parameter list where N is an integer
	 *        defining the length of the list */
	VectorParameterT(const StringT& name_N, char variable = 'v');
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list.
	 * \param list input parameter list, which should be validated using ParameterInterfaceT::ValidateParameterList
	 *        to ensure the list conforms to the description defined by the interface. */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** implicit type conversion */
	operator const dArrayT&() const { return fVector; };

	/** extract parameters to a dArrayT */
	static void Extract(const ParameterListT& list, dArrayT& array, char variable = 'v');

protected:

	/** component names */
	char fVariable;

	/** values */
	dArrayT fVector;
};

/** matrix */
class MatrixParameterT: public ParameterInterfaceT
{
public:

	/** \name constructors */
	/*@{*/
	MatrixParameterT(const StringT& name, int row, int col, char variable = 'A');
	MatrixParameterT(int row, int col, char variable = 'A');

	/** construct extracting dimensions from the name
	 * \param name_NxM name of the matrix parameter list where N and M are the
	 *        integers defining the number of rows and columns, respectively. */
	MatrixParameterT(const StringT& name_NxM, char variable = 'A');
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list.
	 * \param list input parameter list, which should be validated using ParameterInterfaceT::ValidateParameterList
	 *        to ensure the list conforms to the description defined by the interface. */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** implicit type conversion */
	operator const dMatrixT&() const { return fMatrix; };

	/** extract parameters to a dMatrixT */
	static void Extract(const ParameterListT& list, dMatrixT& matrix, char variable = 'A');

protected:

	/** component names */
	char fVariable;

	/** copy values to symmetric positions during MatrixParameterT::TakeParameterList */
	bool fCopySymmetric;

	/** values */
	dMatrixT fMatrix;
};

} /* namespace Tahoe */

#endif /* _PARAMETER_UTILS_H_ */
