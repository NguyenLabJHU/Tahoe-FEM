/* $Id: ParameterTreeT.h,v 1.2.2.3 2003-05-03 17:22:09 paklein Exp $ */
#ifndef _PARAMETER_TREE_T_H_
#define _PARAMETER_TREE_T_H_

/* direct members */
#include "MapT.h"
#include "AutoArrayT.h"
#include "StringT.h"
#include "ParameterListT.h"

namespace Tahoe {

/* forward declarations */
class ParameterInterfaceT;

/** collect and manage parameter lists */
class ParameterTreeT
{
public:

	/** constructor */
	ParameterTreeT(void);

	/** destructor */
	~ParameterTreeT(void);

	/** add a branch to the tree with the given root */
	void BuildDescription(ParameterInterfaceT& root);

	/** add only parts of the branch contained in the guide parameter list */
	void BuildDescription(ParameterInterfaceT& root, const ParameterListT& guide);
//TEMP - needed ?

	/** create a validated parameter list. Take a raw list of parameters and produce 
	 * a validated parameter list. If the validated list cannot be 
	 * produced for any reason, the class throws a ExceptionT::kBadInputValue 
	 * \param source parameter interface associated with the list
	 * \param raw_list source list in which all parameters are stored as
	 *        strings, as read from a source file. 
	 * \param valid_list returns as a validated list witb values of the appropriate data 
	 *        type, validating against constraints and applying any unspecified default values. */
	void Validate(const ParameterInterfaceT& source, const ParameterListT& raw_list, 
		ParameterListT& valid_list);

	/** the branches in the tree */
	const ArrayT<ParameterListT*>& Branches(void) { return fBranches; };

	/** return the parameter list for the given branch or NULL if it does not exist */
	const ParameterListT* Branch(const StringT& name) const;

private:

	/** build the branch. Throws ExceptionT::GeneralFail if the source
	 * has already been added to the dictionary. */
	void BuildBranch(ParameterInterfaceT& source, ParameterListT& params);

	/** build only parts of the branch in the guide list. Throws ExceptionT::GeneralFail
	 * if the source has already been added to the dictionary. */
	void BuildBranch(ParameterInterfaceT& source, const ParameterListT& guide, 
		ParameterListT& params);

	/** \name sub list validation */
	/*@{*/
	bool ValidateSequence(const ParameterInterfaceT& source, 
		const ArrayT<ParameterListT>& raw_sub_list, 
		ParameterListT& valid_list,
		const ArrayT<StringT>& sub_name,
		const ArrayT<ParameterListT::OccurrenceT>& sub_occur,
		const ArrayT<bool>& sub_is_inline, bool throw_on_error);

	bool ValidateChoice(const ParameterInterfaceT& source, 
		const ArrayT<ParameterListT>& raw_sub_list, 
		ParameterListT& valid_list,
		const ArrayT<StringT>& sub_name,
		const ArrayT<ParameterListT::OccurrenceT>& sub_occur,
		const ArrayT<bool>& sub_is_inline, bool throw_on_error);
	/*@}*/

private:

	/** directionary of existing names */
	MapT<StringT, const ParameterInterfaceT*> fDictionary;

	/** branches of the tree */
	AutoArrayT<ParameterListT*> fBranches;
	
	/** list of pointers to delete */
	AutoArrayT<const ParameterInterfaceT*> fDeleteMe;
};

} /* namespace Tahoe */

#endif /* _PARAMETER_TREE_T_H_ */
