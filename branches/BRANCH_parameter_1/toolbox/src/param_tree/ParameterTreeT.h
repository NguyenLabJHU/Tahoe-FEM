/* $Id: ParameterTreeT.h,v 1.2.2.1 2003-04-27 22:19:08 paklein Exp $ */
#ifndef _PARAMETER_TREE_T_H_
#define _PARAMETER_TREE_T_H_

/* direct members */
#include "MapT.h"
#include "AutoArrayT.h"
#include "StringT.h"

namespace Tahoe {

/* forward declarations */
class ParameterInterfaceT;
class ParameterListT;

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

private:

	/** directionary of existing names */
	MapT<StringT, ParameterInterfaceT*> fDictionary;

	/** branches of the tree */
	AutoArrayT<ParameterListT*> fBranches;
	
	/** list of pointers to delete */
	AutoArrayT<ParameterInterfaceT*> fDeleteMe;
};

} /* namespace Tahoe */

#endif /* _PARAMETER_TREE_T_H_ */
