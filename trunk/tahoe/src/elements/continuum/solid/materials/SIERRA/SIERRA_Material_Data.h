/* $Id: SIERRA_Material_Data.h,v 1.2 2003-03-06 17:23:31 paklein Exp $ */
#ifndef _SIERRA_MAT_DATA_H_
#define _SIERRA_MAT_DATA_H_

#include "SIERRA_Material_Interface.h"

/* direct members */
#include "StringT.h"
#include "AutoArrayT.h"

namespace Tahoe {

/** data for a single Sierra material */
class SIERRA_Material_Data
{
public:

	/** constructor */
	SIERRA_Material_Data(const StringT& name, int XML_command_id, int modulus_flag);

	/** registration functions */
	/*@{*/
	void SetCheckFunction(Sierra_function_param_check check) { fCheckFunction = check; };
	
	/** function to do material computations */
	void SetCalcFunction(Sierra_function_material_calc calc) { fCalcFunction = calc; };
	
	/** function to do material initialization */
	void SetInitFunction(Sierra_function_material_init init) { fInitFunction = init; };

	/** register parser line ID */
	void AddXMLCommandID(int ID) { fXMLCommandID.AppendUnique(ID); };
	
	/** register input variable name */
	void AddInputVariable(const StringT& input_var);

	/** set the number of state variables */
	void SetNumStateVariables(int nsv) { fNumStateVars = nsv; }; 

	/** add a material property */
	void AddProperty(const StringT& name, double value);
	/*@}*/

	/** \name accessors */
	/*@{*/
	/** material name */
	const StringT& Name(void) const { return fName; };
	
	/** return the list of command IDs */
	const ArrayT<int>& XMLCommandID(void) const { return fXMLCommandID; };

	/** return the number of state variables */
	int NumStateVariables(void) { return fNumStateVars; }; 

	/** return the list of input variable names */
	const ArrayT<StringT>& InputVariables(void) { return fInputVariables; };

	/** return the list sized per variable */
	const ArrayT<int>& InputVariableSize(void) { return fInputVariableSize; };

	/** function to do checking/registration of parameters */
	Sierra_function_param_check CheckFunction(void) { return fCheckFunction; };
	
	/** function to do material computations */
	Sierra_function_material_calc CalcFunction(void) { return fCalcFunction; };
	
	/** function to do material initialization */
	Sierra_function_material_init InitFunction(void) { return fInitFunction; };

	/** array property names */
	const ArrayT<StringT>& PropertyNames(void) const { return fPropertyNames; };

	/** array property values */
	const ArrayT<double>&  PropertyValues(void) const { return fPropertyValues; };
	/*@}*/
	
private:

	/** material name */
	StringT fName;

	/** 0/1 = dont/do complete elastic constants */
	int fModulusFlag;

	/** number of state variables */
	int fNumStateVars;

	/** \name registered functions */
	/*@{*/
	/** function to do checking/registration of parameters */
	Sierra_function_param_check fCheckFunction;
	
	/** function to do material computations */
	Sierra_function_material_calc fCalcFunction;
	
	/** function to do material initialization */
	Sierra_function_material_init fInitFunction;
	/*@}*/
	
	/** list of XML command ids. The first command is the one defining the material
	 * itself */
	AutoArrayT<int> fXMLCommandID;
	 
	/** \name list of input variable names and offsets */
	/*@{*/
	AutoArrayT<StringT> fInputVariables;
	AutoArrayT<int>     fInputVariableSize;
	/*@}*/

	/** \name material properties */
	/*@{*/
	AutoArrayT<StringT> fPropertyNames;
	AutoArrayT<double>  fPropertyValues;
	/*@}*/
};

/* add a material property */
void SIERRA_Material_Data::AddProperty(const StringT& name, double value)
{
	if (fPropertyNames.AppendUnique(name))
		fPropertyValues.Append(value);
}

} /* namespace Tahoe */

#endif /* _SIERRA_MAT_DATA_H_ */
