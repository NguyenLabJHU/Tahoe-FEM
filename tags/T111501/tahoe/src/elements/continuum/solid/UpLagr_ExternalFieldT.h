/* $Id: UpLagr_ExternalFieldT.h,v 1.1 2001-07-14 01:16:40 paklein Exp $ */

#ifndef _UPDATED_LAGR_EX_FIELD_T_H_
#define _UPDATED_LAGR_EX_FIELD_T_H_

/* base class */
#include "UpdatedLagrangianT.h"

/* direct members */
#include "dRangeArrayT.h"

/** extension to the updated Lagrangian, finite strain solid
 * elements that reads external field variables from an
 * existing results database. the external data is interpolated
 * to the current time. \note the class currently
 * requires an ExodusII database for the external field
 * variables. */
class UpLagr_ExternalFieldT: public UpdatedLagrangianT
{
public:

	/** constructor */
	UpLagr_ExternalFieldT(FEManagerT& fe_manager);

	/** initialization */
	virtual void Initialize(void);

	/** external nodal field variables for the current element */
	const LocalArrayT& ExternalField(void) const;

	/** initialize current time increment. interpolate external field
	 * data from database to the current time step. */
	virtual void InitStep(void);

protected:

	/** increment current element. collects nodal variables of the
	 * external field for the current element. 
	 * \return false if the end of the element list has been hit. */
	virtual bool NextElement(void);	


private:

	/** assemble values into external field using the node map */
	void AssembleField(int col, double scale, const dArrayT& values);	

protected:

	/* external field variables */
	IOBaseT::FileTypeT fExternalFieldFormat; /**< file format of external file */
	StringT fExternalFieldFile; /**< (relative path) to external database */ 
	
	/** array of labels for the nodal field variables */
	ArrayT<StringT> fExternalFieldLabels;
	iArrayT fFieldVariableIndex; /**< index of field variable in database */

	/** list of times in the external results database */
	dRangeArrayT fTimeSteps;
	
	/** map of database node numbers */
	iArrayT fNodeMap;

	/** external field variables in local ordering */
	LocalArrayT fLocExternalField;
	
	/** external field variables for all nodes */
	dArray2DT fExternalField;
	
	/** vector used to read nodal values from the database */
	dArrayT fNodalValues;
};

/* inlines */
inline const LocalArrayT& UpLagr_ExternalFieldT::ExternalField(void) const
{
	return fLocExternalField;
}

#endif /* _UPDATED_LAGR_EX_FIELD_T_H_ */
