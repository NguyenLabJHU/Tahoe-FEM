/* $Id: ParameterTreeT.cpp,v 1.1.2.4 2003-05-03 09:06:52 paklein Exp $ */
#include "ParameterTreeT.h"
#include "ParameterInterfaceT.h"

using namespace Tahoe;

ParameterTreeT::ParameterTreeT(void) {}

ParameterTreeT::~ParameterTreeT(void)
{
	/* free parameter list data */
	for (int i = 0; i < fBranches.Length(); i++)
		delete fBranches[i];

	/* free instantiations */
	for (int i = 0; i < fDeleteMe.Length(); i++)
		delete fDeleteMe[i];
}

/* add a branch to the tree with the given root */
void ParameterTreeT::BuildDescription(ParameterInterfaceT& root)
{
	/* branch exists */
	if (fDictionary.HasKey(root.Name())) return;

	/* collect parameters into a new list */
	ParameterListT* root_list = new ParameterListT(root.Name());
	BuildBranch(root, *root_list);
	
	/* store branch */
	fBranches.Append(root_list);
}

/* add only parts of the branch contained in the guide parameter list */
void ParameterTreeT::BuildDescription(ParameterInterfaceT& root, const ParameterListT& guide)
{
	/* branch exists */
	if (fDictionary.HasKey(root.Name())) return;

	/* collect parameters into a new list */
	ParameterListT* root_list = new ParameterListT(root.Name());
	BuildBranch(root, guide, *root_list);
	
	/* store branch */
	fBranches.Append(root_list);
}

/* create a validated parameter list */
void ParameterTreeT::Validate(const ParameterInterfaceT& source, const ParameterListT& raw_list, 
	ParameterListT& valid_list)
{
	const char caller[] = "ParameterTreeT::Validate";
	const char* stage_name[] = {"parameters", "sub-lists"};
	int stage = 0;
	
	try { /* catch all exceptions */

	/* clear destination */
	valid_list.Clear();

	/* contents of subordinate lists */
	ArrayT<StringT> sub;
	ArrayT<ParameterListT::OccurrenceT> sub_occur;
	ArrayT<bool> sub_is_inline;

	/* validating self */
	if (source.Name() == raw_list.Name())
	{
		stage = 0;

		/* collect parameters */
		valid_list.SetName(raw_list.Name());
		source.ValidateParameters(raw_list, valid_list);
		stage = 1;
			
		/* add to the dictionary */
		fDictionary.Insert(source.Name(), &source);

		/* collect information on subordinate lists */
		source.SubNames(sub, sub_occur, sub_is_inline);

		/* list order */
		ParameterListT::ListOrderT order = source.ListOrder();
		if (order == ParameterListT::Sequence)
			ValidateSequence(source, raw_list.Lists(), valid_list, sub, sub_occur, sub_is_inline, true);
		else if (order == ParameterListT::Choice)
			ValidateChoice(source, raw_list.Lists(), valid_list, sub, sub_occur, sub_is_inline, true);
		else
			ExceptionT::GeneralFail(caller, "unrecognized list order %d", order);
	}
	else /* defining inlined list */
	{

		//TEMP
		ExceptionT::GeneralFail(caller, "why is there an inlined list here?????????");

#if 0
		/* list must be inlined */
		if (!raw_list.Inline())
			ExceptionT::GeneralFail(caller, "source \"%s\" cannot build list for non-inlined \"%s\"",
				source.Name().Pointer(), raw_list.Name().Pointer());
		
		/* collect information on inlined list */
		ParameterListT::ListOrderT order;
		source.DefineInlineSub(raw_list.Name(), order, sub, sub_occur, sub_is_inline);
#endif
	}
	
	} /* end try */

	/* rethrow */
	catch (ExceptionT::CodeT error) {
		ExceptionT::BadInputValue(caller, "exception %d validating %s in \"%s\"", 
			error, stage_name[stage], source.Name().Pointer());
	}
}

/* return the parameter list for the given branch or NULL if it does not exist */
const ParameterListT* ParameterTreeT::Branch(const StringT& name) const
{
	for (int i = 0; i < fBranches.Length(); i++)
		if (fBranches[i]->Name() == name)
			return fBranches[i];
		
	return NULL;
}

/*************************************************************************
 * Private
 *************************************************************************/

bool ParameterTreeT::ValidateChoice(
	const ParameterInterfaceT& source, 
	const ArrayT<ParameterListT>& raw_sub_list, 
	ParameterListT& valid_list,
	const ArrayT<StringT>& sub_name,
	const ArrayT<ParameterListT::OccurrenceT>& sub_occur,
	const ArrayT<bool>& sub_is_inline,
	bool throw_on_error)
{
	const char caller[] = "ParameterTreeT::ValidateChoice";

	/* can't match an empty list */
	if (raw_sub_list.Length() == 0 && sub_name.Length() > 0)
	{
		if (throw_on_error)
			ExceptionT::BadInputValue(caller, "empty source unable to match choice at position %d in \"%s\"",
				valid_list.NumLists()+1, valid_list.Name().Pointer());
		return false;
	}

	/* TEMP - just try to match first item */
	const ParameterListT& raw_sub_0 = raw_sub_list[0];
	const StringT& name_0 = raw_sub_0.Name();

	/* search possible sub names */
	for (int i = 0; i < sub_name.Length(); i++)
	{
		//TEMP
		if (sub_is_inline[i])
			ExceptionT::GeneralFail(caller, "inlined list \"%s\" not supported as choice in \"%s\"",
				sub_name[i].Pointer(), valid_list.Name().Pointer());
	
		/* check */
		if (sub_occur[i] != ParameterListT::Once)
			ExceptionT::GeneralFail(caller, "occurrence of items in choice must be Once");
	
		/* found match */
		if (sub_name[i] == name_0)
		{
			/* get parameter source */
			const ParameterInterfaceT* sub = NULL;
			if (fDictionary.HasKey(name_0))
				sub = fDictionary[name_0];
			else {
				sub = source.NewSub(name_0);
				fDeleteMe.Append(sub);
			}
		
			/* create validated sub-list */
			ParameterListT valid_sub_list;
			Validate(*sub, raw_sub_0, valid_sub_list);
			valid_list.AddList(valid_sub_list);

			return true;
		}
	}
	
	/* fail on passing through */
	if (throw_on_error)
		ExceptionT::BadInputValue(caller, "\"%s\" does not match choice at position %d in \"%s\"",
			name_0.Pointer(), valid_list.NumLists()+1, valid_list.Name().Pointer());
	return false;
}

bool ParameterTreeT::ValidateSequence(
	const ParameterInterfaceT& source, 
	const ArrayT<ParameterListT>& raw_sub_list, 
	ParameterListT& valid_list,
	const ArrayT<StringT>& sub_name,
	const ArrayT<ParameterListT::OccurrenceT>& sub_occur,
	const ArrayT<bool>& sub_is_inline,
	bool throw_on_error)
{
	const char caller[] = "ParameterTreeT::ValidateSequence";

	bool check_OK = true;
	int raw_dex = 0;
	int sub_dex = 0;
	for (sub_dex = 0; check_OK && sub_dex < sub_name.Length() && raw_dex < raw_sub_list.Length(); sub_dex++)
	{
		/* next list entry */
		bool is_inline = sub_is_inline[sub_dex];
		const StringT& name = sub_name[sub_dex];
		ParameterListT::OccurrenceT occur = sub_occur[sub_dex];
		
		/* inlined sub list */
		if (is_inline)
		{
			/* get sub-list description */
			ParameterListT::ListOrderT in_order;
			ArrayT<StringT> in_sub_name;
			ArrayT<ParameterListT::OccurrenceT> in_sub_occur;
			ArrayT<bool> in_sub_is_inline;
			source.DefineInlineSub(name, in_order, in_sub_name, in_sub_occur, in_sub_is_inline);

			switch (occur)
			{
				case ParameterListT::Once:
				case ParameterListT::ZeroOrOnce:
				{
					/* alias to tail of the parameter list */
					const ArrayT<ParameterListT> raw_sub_list_tail(raw_sub_list.Length() - raw_dex, raw_sub_list.Pointer(raw_dex));
					
					/* error handling */
					bool in_throw_on_error = (occur == ParameterListT::Once) ? throw_on_error : false;
				
					/* must have one */
					int last_list_count = valid_list.NumLists();
					if (in_order == ParameterListT::Sequence)
						check_OK = ValidateSequence(source, raw_sub_list_tail, valid_list, in_sub_name, in_sub_occur, in_sub_is_inline, in_throw_on_error);
					else if (in_order == ParameterListT::Choice)
						check_OK = ValidateChoice(source, raw_sub_list_tail, valid_list, in_sub_name, in_sub_occur, in_sub_is_inline, in_throw_on_error);
					else
						ExceptionT::GeneralFail(caller, "unrecognized list order %d", in_order);

					/* assumes the number of values added to valid list is the number removed from the raw list */
					raw_dex += valid_list.NumLists() - last_list_count; 
				
					/* check */
					if (!check_OK && in_throw_on_error)
						ExceptionT::BadInputValue(caller, "error with inlined item \"%s\" at position %d in \"%s\"", 
							name.Pointer(), valid_list.NumLists() + 1, source.Name().Pointer());
					break;
				}
				case ParameterListT::OnePlus:
				case ParameterListT::Any:
				{
					/* accept many */
					int count = 0;
					while (raw_dex < raw_sub_list.Length() && check_OK)
					{
						/* alias to tail of the parameter list */
						const ArrayT<ParameterListT> raw_sub_list_tail(raw_sub_list.Length() - raw_dex, raw_sub_list.Pointer(raw_dex));

						/* error handling */
						bool in_throw_on_error = (count == 0 && occur == ParameterListT::OnePlus) ? throw_on_error : false;
					
						/* search */
						int last_list_count = valid_list.NumLists();
						if (in_order == ParameterListT::Sequence)
							check_OK = ValidateSequence(source, raw_sub_list_tail, valid_list, in_sub_name, in_sub_occur, in_sub_is_inline, in_throw_on_error);
						else if (in_order == ParameterListT::Choice)
							check_OK = ValidateChoice(source, raw_sub_list_tail, valid_list, in_sub_name, in_sub_occur, in_sub_is_inline, in_throw_on_error);
						else
							ExceptionT::GeneralFail(caller, "unrecognized list order %d", in_order);

						if (check_OK) count++;

						/* assumes the number of values added to valid list is the number removed from the raw list */
						raw_dex += valid_list.NumLists() - last_list_count; 
					}
						
					/* must have at least one full pass */
					if (occur == ParameterListT::OnePlus && count < 1)
					{
						check_OK = false;
						if (throw_on_error)
							ExceptionT::BadInputValue(caller, "error with inlined item \"%s\" at position %d in \"%s\"", 
								name.Pointer(), valid_list.NumLists() + 1, source.Name().Pointer());
					}
					else /* no problem */
						check_OK = true;
					break;
				}
				default:
					ExceptionT::GeneralFail(caller, "unrecognized occurrence %d", occur);
			}		
		}
		else /* process non-inlined entries */
		{
			switch (occur)
			{
				case ParameterListT::Once:
				case ParameterListT::ZeroOrOnce:
				{
					/* look for one */
					int count = 0;
					if (raw_sub_list[raw_dex].Name() == name) 
					{
						count++;
					
						/* get parameter source */
						const ParameterInterfaceT* sub = NULL;
						if (fDictionary.HasKey(name))
							sub = fDictionary[name];
						else {
							sub = source.NewSub(name);
							fDeleteMe.Append(sub);
						}

						/* create validated sub-list */
						ParameterListT valid_sub_list;
						Validate(*sub, raw_sub_list[raw_dex], valid_sub_list);
						valid_list.AddList(valid_sub_list);
							
						/* next source item */
						raw_dex++;
					}

					/* check */
					if (occur == ParameterListT::Once && count != 1)
					{
						check_OK = false;						
						if (throw_on_error)
							ExceptionT::BadInputValue(caller, "item %d in \"%s\" must be \"%s\"",
								valid_list.NumLists() + 1, valid_list.Name().Pointer(), name.Pointer());
					}	
					break;
				}
				case ParameterListT::OnePlus:
				case ParameterListT::Any:
				{
					/* accept many */
					int count = 0;
					while (raw_dex < raw_sub_list.Length() && raw_sub_list[raw_dex].Name() == name)
					{
						count++;

						/* get parameter source */
						const ParameterInterfaceT* sub = NULL;
						if (fDictionary.HasKey(name))
							sub = fDictionary[name];
						else {
							sub = source.NewSub(name);
							fDeleteMe.Append(sub);
						}

						/* create validated sub-list */
						ParameterListT valid_sub_list;
						Validate(*sub, raw_sub_list[raw_dex], valid_sub_list);
						valid_list.AddList(valid_sub_list);
				
						/* next source item */
						raw_dex++;
					}
						
					/* must have at least one */
					if (occur == ParameterListT::OnePlus && count < 1)
					{
						check_OK = false;
						if (throw_on_error)
							ExceptionT::BadInputValue(caller, "list \"%s\" must contain at least one \"%s\" beginning at position %d",
								valid_list.Name().Pointer(), name.Pointer(), valid_list.NumLists() + 1);
					}
					break;
				}
				default:
					ExceptionT::GeneralFail(caller, "unrecognized occurrence %d", occur);
			}
		}
	}
			
	/* check any remaining */
	for (;check_OK && sub_dex < sub_name.Length(); sub_dex++)
		if (sub_occur[sub_dex] != ParameterListT::ZeroOrOnce && sub_occur[sub_dex] != ParameterListT::Any)
		{
			check_OK = false;
			if (throw_on_error)
				ExceptionT::GeneralFail(caller, "\"%s\" required at position %d in \"%s\"",
					sub_name[sub_dex].Pointer(), valid_list.NumLists()+1, valid_list.Name().Pointer());
		}

	return check_OK;
}

/* build the branch */
void ParameterTreeT::BuildBranch(ParameterInterfaceT& source, ParameterListT& params)
{
	const char caller[] = "ParameterTreeT::BuildBranch";
	
	/* contents of subordinate lists */
	ArrayT<StringT> sub_lists;
	ArrayT<ParameterListT::OccurrenceT> occur;
	ArrayT<bool> is_inline;

	/* defining self */
	if (source.Name() == params.Name()) 
	{	
		/* collect parameters */
		source.DefineParameters(params);
			
		/* add list to the dictionary */
		if (!fDictionary.Insert(source.Name(), &source))
			ExceptionT::GeneralFail(caller, "dictionary already contains \"%s\"",
				source.Name().Pointer());

		/* collect information on subordinate lists */
		source.SubNames(sub_lists, occur, is_inline);
	}
	else /* defining inlined list */
	{
		/* list must be inlined */
		if (!params.Inline())
			ExceptionT::GeneralFail(caller, "source \"%s\" cannot build list for non-inlined \"%s\"",
				source.Name().Pointer(), params.Name().Pointer());
		
		/* collect information on inlined list */
		ParameterListT::ListOrderT order;
		source.DefineInlineSub(params.Name(), order, sub_lists, occur, is_inline);
		params.SetListOrder(order);
	}

	/* define sublists */
	for (int i = 0; i < sub_lists.Length(); i++) {

		/* list is inline */
		if (is_inline[i])
		{
			/* source builds inline list */
			ParameterListT inline_params(sub_lists[i]);
			inline_params.SetInline(true);
			BuildBranch(source, inline_params);
			
			/* add list */
			if (!params.AddList(inline_params, occur[i]))
				ExceptionT::GeneralFail(caller, "could not add inline \"%s\" to list \"%s\"",
					inline_params.Name().Pointer(), params.Name().Pointer());
		}
		/* already in the tree */
		else if (fDictionary.HasKey(sub_lists[i]))
		{
			/* add a blank list as a placeholder */
			ParameterListT sub_params(sub_lists[i]);
			if (!params.AddList(sub_params, occur[i]))
				ExceptionT::GeneralFail(caller, "could not add sublist \"%s\" to list \"%s\"",
					sub_params.Name().Pointer(), params.Name().Pointer());
		}
		else /* new tree entry */
		{
			/* sublist */
			ParameterInterfaceT* sub = source.NewSub(sub_lists[i]);
			if (!sub)
				ExceptionT::GeneralFail(caller, "source \"%s\" did not return sublist \"%s\"",
					params.Name().Pointer(), sub_lists[i].Pointer());
					
			/* build the sublist */
			ParameterListT sub_params(sub->Name());
			BuildBranch(*sub, sub_params);
			
			/* add list */
			if (!params.AddList(sub_params, occur[i]))
				ExceptionT::GeneralFail(caller, "could not add sublist \"%s\" to list \"%s\"",
					sub_params.Name().Pointer(), params.Name().Pointer());
	
			/* add to list of items to delete */
			fDeleteMe.Append(sub);
		}
	}	
}

/* build the branch */
void ParameterTreeT::BuildBranch(ParameterInterfaceT& source, const ParameterListT& guide, 
	ParameterListT& params)
{
	class list_info {
		public:
			list_info(void): occur(ParameterListT::Undefined), parent(NULL) {};
		public:
			ParameterListT::OccurrenceT occur;
			list_info* parent;
	};

	const char caller[] = "ParameterTreeT::BuildBranch";

	/* collect parameters */
	source.DefineParameters(params);
	
	/* add list to the dictionary */
	if (!fDictionary.Insert(source.Name(), &source))
		ExceptionT::GeneralFail(caller, "dictionary already contains \"%s\"",
			source.Name().Pointer());

	/* sublist information */
	ArrayT<StringT> names;
	ArrayT<ParameterListT::OccurrenceT> occur;
	ArrayT<bool> is_inline;
	source.SubNames(names, occur, is_inline);
	
	/* construct sublists */
	const ArrayT<ParameterListT>& sub_lists = guide.Lists();
	for (int i = 0; i < sub_lists.Length(); i++) {

		/* sublist name */
		const StringT& sub_name = sub_lists[i].Name();

		/* fetch sublist */
		ParameterInterfaceT* sub = source.NewSub(sub_name);

		/* not found */
		if (!sub)
			ExceptionT::GeneralFail(caller, "source \"%s\" did not return sublist \"%s\"",
				params.Name().Pointer(), sub_name.Pointer());

		/* already in the tree */
		else if (fDictionary.HasKey(sub_name))
		{
			/* find occurrence */
			ParameterListT::OccurrenceT sub_occur = ParameterListT::Undefined;
			for (int j = 0; sub_occur == ParameterListT::Undefined && j < names.Length(); j++)
				if (names[j] == sub_name)
					sub_occur = occur[j];

			/* add a blank list as a placeholder */
			ParameterListT sub_params(sub_name);
			if (!params.AddList(sub_params, sub_occur))
				ExceptionT::GeneralFail(caller, "could not add sublist \"%s\" to list \"%s\"",
					sub_params.Name().Pointer(), params.Name().Pointer());
		}
		else /* new tree entry */
		{
			/* find occurrence */
			ParameterListT::OccurrenceT sub_occur = ParameterListT::Undefined;
			for (int j = 0; sub_occur == ParameterListT::Undefined && j < names.Length(); j++)
				if (names[j] == sub_name)
					sub_occur = occur[j];
					
			/* check if inlined list - could cache inlines lists in source */
			if (sub_occur == ParameterListT::Undefined)
			{
				for (int j = 0; sub_occur == ParameterListT::Undefined && j < is_inline.Length(); j++)
					if (is_inline[j])
					{
						/* get inlined list information */
						ParameterListT::ListOrderT order_j;
						ArrayT<StringT> names_j;
						ArrayT<ParameterListT::OccurrenceT> occur_j;
						ArrayT<bool> is_inline_j;
						source.DefineInlineSub(names[j], order_j, names_j, occur_j, is_inline_j);
						
						/* search */
						for (int k = 0; sub_occur == ParameterListT::Undefined && k < names_j.Length(); k++)
							if (names_j[k] == sub_name)
								sub_occur = occur[j]; /* gets occurrence of nested list */
					}					
			}
			
			/* not found */
			if (sub_occur == ParameterListT::Undefined)
				ExceptionT::GeneralFail(caller, "sublist \"%s\" not found in \"%s\" or in the first level of nested lists",
					sub_name.Pointer(), source.Name().Pointer());

			/* build the sublist */
			ParameterListT sub_params(sub->Name());
			BuildBranch(*sub, sub_lists[i], sub_params);

			/* add list */
			if (!params.AddList(sub_params, sub_occur))
				ExceptionT::GeneralFail(caller, "could not add sublist \"%s\" to list \"%s\"",
					sub_params.Name().Pointer(), params.Name().Pointer());
			
			/* add to list of items to delete */
			fDeleteMe.Append(sub);
		}
	}
}
