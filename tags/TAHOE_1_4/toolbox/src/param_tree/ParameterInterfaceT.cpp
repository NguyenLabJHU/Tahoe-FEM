/* $Id: ParameterInterfaceT.cpp,v 1.19 2004-07-12 21:49:59 paklein Exp $ */
#include "ParameterInterfaceT.h"
#include "ParameterListT.h"
#include "ParameterUtils.h"
#include "ParameterContainerT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<ParameterInterfaceT*>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<const ParameterInterfaceT*>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<SubListDescriptionT>::fByteCopy = false;
}

/* constructor */
ParameterInterfaceT::ParameterInterfaceT(const StringT& name)
{
	SetName(name);
}

void ParameterInterfaceT::SetName(const StringT& name)
{
	fName = name;
}

/* accept completed parameter list */
void ParameterInterfaceT::TakeParameterList(const ParameterListT& list) 
{
	const char caller[] = "ParameterInterfaceT::TakeParameterList";
	if (list.Name() != Name())
		ExceptionT::GeneralFail(caller, "list name \"%s\" must be \"%s\"",
			list.Name().Pointer(), Name().Pointer());
}

/* build complete parameter list description */
void ParameterInterfaceT::DefineParameters(ParameterListT& list) const
{
	const char caller[] = "ParameterInterfaceT::DefineParameters";
	if (list.Name() != Name())
		ExceptionT::GeneralFail(caller, "list name \"%s\" must be \"%s\"",
			list.Name().Pointer(), Name().Pointer());
}

/* validate the given parameter list */
void ParameterInterfaceT::ValidateParameterList(const ParameterListT& raw_list, ParameterListT& valid_list) const
{
	/* check name */
	const char caller[] = "ParameterInterfaceT::ValidateParameterList";
	if (raw_list.Name() != Name())
		ExceptionT::GeneralFail(caller, "raw list name \"%s\" must be \"%s\"",
			raw_list.Name().Pointer(), Name().Pointer());
	if (valid_list.Name() != Name())
		ExceptionT::GeneralFail(caller, "list name \"%s\" must be \"%s\"",
			raw_list.Name().Pointer(), Name().Pointer());

	/* collect description */
	ParameterListT description(Name());
	DefineParameters(description);

	/* parameter lists */
	const ArrayT<ParameterT>& raw_parameters = raw_list.Parameters();
	const ArrayT<ParameterT>& parameters = description.Parameters();
	const ArrayT<ParameterListT::OccurrenceT>& occurrences = description.ParameterOccurrences();

	/* look through parameters */
	for (int i = 0; i < parameters.Length(); i++)
	{
		const ParameterT& parameter = parameters[i];
		ParameterListT::OccurrenceT occurrence = occurrences[i];

		/* search entire list - parameters are not ordered */
		int count = 0;
		for (int j = 0; j < raw_parameters.Length(); j++) {

			const ParameterT& raw_parameter = raw_parameters[j];

			/* name match */
			if (raw_parameter.Name() == parameter.Name()) {
			
				/* check occurrence */
				if (count > 0 && (occurrence == ParameterListT::Once || occurrence == ParameterListT::ZeroOrOnce))
					ExceptionT::BadInputValue(caller, "parameter \"%s\" cannot appear in \"%s\" more than once",
						parameter.Name().Pointer(), Name().Pointer());
									
				/* get value */
				ParameterT new_parameter(parameter.Type(), parameter.Name());
				new_parameter.SetDescription(parameter.Description());
				if (raw_parameter.Type() == parameter.Type())
					new_parameter = raw_parameter;
				else if (raw_parameter.Type() == ParameterT::String || raw_parameter.Type() == ParameterT::Word)
				{
					/* convert to string */
					const StringT& value_str = raw_parameter;

					/* set from string */
					new_parameter.FromString(value_str);	
				}
				else
					ExceptionT::BadInputValue(caller, "source for \"%s\" must have type %d, %d, or %d: %d", 
						parameter.Name().Pointer(), ParameterT::String, ParameterT::Word, parameter.Type());

				/* increment count */
				count++;
				
				/* check limits */
				if (parameter.InBounds(new_parameter))
				{
					/* fix string-value pairs */
					if (parameter.Type() == ValueT::Enumeration)
						parameter.FixEnumeration(new_parameter);

					/* copy in limits */
					new_parameter.AddLimits(parameter.Limits());
				
					/* copy in default */
					const ValueT* default_value = parameter.Default();
					if (default_value)
						new_parameter.SetDefault(*default_value);
				
					/* add it */
					valid_list.AddParameter(new_parameter);
				}
				else 
				{
					/* check again, now verbose */
					parameter.InBounds(new_parameter, true);
					ExceptionT::BadInputValue(caller, "improper value for parameter \"%s\" in \"%s\"",
						parameter.Name().Pointer(), Name().Pointer());
				}
			}
		}
				
		/* look for default value */
		if (count == 0 && (occurrence == ParameterListT::Once || occurrence == ParameterListT::OnePlus)) 
		{
			const ValueT* default_value = parameter.Default();
			if (default_value) {
				ParameterT new_parameter(parameter.Type(), parameter.Name());
				new_parameter.SetDescription(parameter.Description());
				new_parameter = *default_value;

				/* check limits */
				if (parameter.InBounds(new_parameter))
				{
					/* copy in limits */
					new_parameter.AddLimits(parameter.Limits());

					/* copy in default */
					new_parameter.SetDefault(*default_value);
				
					/* add it */
					valid_list.AddParameter(new_parameter);
				}
				else 
				{	
					/* check again, now verbose */
					parameter.InBounds(new_parameter, true);
					ExceptionT::BadInputValue(caller, "improper value for parameter \"%s\" in \"%s\"",
						parameter.Name().Pointer(), Name().Pointer());
				}
			}
			else
				ExceptionT::BadInputValue(caller, 
					"required parameter \"%s\" is missing from \"%s\" and has no default",
					parameter.Name().Pointer(), Name().Pointer());
		}
	}
}

/* the order of subordinate lists */
ParameterListT::ListOrderT ParameterInterfaceT::ListOrder(void) const
{
	return ParameterListT::Sequence;
}

/* return the list of subordinate names */
void ParameterInterfaceT::DefineSubs(SubListT& sub_list) const
{
#pragma unused(sub_list)
//	sub_list.Dimension(0);
//NOTE: clearing the list causes classes which inherit ParameterInterfaceT more than once
//      to have an incomplete list of subs.
}

void ParameterInterfaceT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order,
	SubListT& sub_lists) const
{
	/* check for definition in NewSub */
	ParameterInterfaceT* inline_sub = NewSub(name);
	if (inline_sub) /* define inline sub */ {
	
		/* set sequence or choice */
		order = inline_sub->ListOrder();
		
		/* there should be no parameters */
		ParameterListT params(inline_sub->Name());
		inline_sub->DefineParameters(params);
		if (params.NumParameters() > 0)
			ExceptionT::GeneralFail("ParameterInterfaceT::DefineInlineSub", 
				"%d parameters not allowed in inline sub \"%s\"", params.NumParameters(), name.Pointer());
		
		/* get subs */
		inline_sub->DefineSubs(sub_lists);
		
		/* clean up */
		delete inline_sub;
	}
	else /* class must override DefineInlineSub to define the sub */ {
		order = ParameterListT::Sequence;
		sub_lists.Dimension(0);
	}
}

/* return a pointer to the ParameterInterfaceT */
ParameterInterfaceT* ParameterInterfaceT::NewSub(const StringT& name) const
{
	const char caller[] = "ParameterInterfaceT::NewSub";

	/* look for names ending in "_ID_list" */
	const char* find = strstr(name, "_ID_list");
	if (find && strlen(find) == 8) {
		StringListT* id_list = new StringListT(name);
		id_list->SetMinLength(1);
		return id_list;
	}

	if (name == "Integer")
		return new IntegerParameterT;
	else if (name == "IntegerList")
		return new IntegerListT;
	else if (name == "Double")
		return new DoubleParameterT;
	else if (name == "DoubleList")
		return new DoubleListT;
	else if (name == "String")
		return new StringParameterT;
	else if (name == "OrderedPair")
	{
		ParameterContainerT* pair = new ParameterContainerT("OrderedPair");
		ParameterT x(ParameterT::Double, "x");
		ParameterT y(ParameterT::Double, "y");
		pair->AddParameter(x);
		pair->AddParameter(y);
		return pair;
	}
	else if (strncmp("Vector_", name, 7) == 0) /* Vector_N */
		return new VectorParameterT(name);
	else if (strncmp("Matrix_", name, 7) == 0) /* Matrix_MxN */
		return new MatrixParameterT(name);
	else
		return NULL;
}

/* remove the first instance of the given sublist */
bool SubListT::RemoveSub(const char* name)
{
	/* scan */
	int index = -1;
	for (int i = 0; index == -1 && i < fLength; i++)
		if (fArray[i].Name() == name)
			index = i;

	/* remove */
	if (index != -1) {
		DeleteAt(index);
		return true;
	}	
	else
		return false;
}
