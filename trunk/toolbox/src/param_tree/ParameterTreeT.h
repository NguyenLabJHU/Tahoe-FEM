/* $Id: ParameterTreeT.h,v 1.4 2003-08-14 01:22:03 paklein Exp $ */
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
class SubListDescriptionT;

/** collect and manage parameter lists. There are two major components of this
 * interface:
 * - ParameterTreeT::BuildDescription adds a branch to the parameter beginning
 *   at the given ParameterInterfaceT. There branches are accessible through
 *   ParameterTreeT::Branches or ParameterTreeT::Branch.
 * - ParameterTreeT::Validate takes a raw parameter list and attempts to build
 *   a validated list beginning with the given ParameterInterfaceT.
 */
class ParameterTreeT
{
public:

	/** constructor */
	ParameterTreeT(void);

	/** destructor */
	~ParameterTreeT(void);

	/** add a branch to the tree with the given root */
	void BuildDescription(const ParameterInterfaceT& root);

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
	void BuildBranch(const ParameterInterfaceT& source, ParameterListT& params);

	/** \name sub list validation */
	/*@{*/
	bool ValidateSequence(const ParameterInterfaceT& source, 
		const ArrayT<ParameterListT>& raw_sub_list, 
		ParameterListT& valid_list,
		const AutoArrayT<SubListDescriptionT>& sub_list,
		bool throw_on_error);

	bool ValidateChoice(const ParameterInterfaceT& source, 
		const ArrayT<ParameterListT>& raw_sub_list, 
		ParameterListT& valid_list,
		const AutoArrayT<SubListDescriptionT>& sub_list,
		bool throw_on_error);
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
