/* $Id: SIERRA_Material_Data.h,v 1.7 2004-08-16 17:27:17 paklein Exp $ */
#ifndef _SIERRA_MAT_DATA_H_
#define _SIERRA_MAT_DATA_H_

#include "SIERRA_Material_Interface.h"

/* direct members */
#include "StringT.h"
#include "AutoArrayT.h"
#include "MapT.h"

namespace Tahoe {

/** data for a single Sierra material */
class SIERRA_Material_Data
{
public:

	/** enumeration of input variables */
	enum InputVariableT {
		kundefined,
		krot_strain_increment,
		krot_strain_inc,
		kvelocity_gradient,
		ktemperature_old,
		ktemperature_new
	};

	/** constructor */
	SIERRA_Material_Data(const StringT& name, int XML_command_id, int modulus_flag);

	/** registration functions */
	/*@{*/
	void SetCheckFunction(Sierra_function_param_check check) { fCheckFunction = check; };
	
	/** function to do material computations */
	void SetCalcFunction(Sierra_function_material_calc calc) { fCalcFunction = calc; };
	
	/** function to do material initialization */
	void SetInitFunction(Sierra_function_material_init init) { fInitFunction = init; };

	/** function to compute elastic moduli */
	void SetPCFunction(Sierra_pc_elastic_moduli_func func) { fPCFunction = func; };

	/** register parser line ID */
	void AddXMLCommandID(int ID) { fXMLCommandID.AppendUnique(ID); };
	
	/** register input variable name */
	void AddInputVariable(const StringT& input_var);

	/** set the number of state variables */
	void SetNumStateVariables(int nsv) { fNumStateVars = nsv; }; 
	/*@}*/

	/** \name accessors */
	/*@{*/
	/** unique material identifier */
	int ID(void) const { return fID; };
	
	/** material name */
	const StringT& Name(void) const { return fName; };
	
	/** return the list of command IDs */
	const ArrayT<int>& XMLCommandID(void) const { return fXMLCommandID; };

	/** return the number of state variables */
	int NumStateVariables(void) const { return fNumStateVars; }; 

	/** return the list of input variable names */
	const ArrayT<StringT>& InputVariables(void) const { return fInputVariables; };

	/** return the list sized per variable */
	const ArrayT<int>& InputVariableSize(void) const { return fInputVariableSize; };
	
	/** return the type of strain measure used by the model */
	InputVariableT StrainMeasure(void) const { return fStrainMeasure; };

	/** function to do checking/registration of parameters */
	Sierra_function_param_check CheckFunction(void) const { return fCheckFunction; };
	
	/** function to do material computations */
	Sierra_function_material_calc CalcFunction(void) const { return fCalcFunction; };
	
	/** function to do material initialization */
	Sierra_function_material_init InitFunction(void) const { return fInitFunction; };
	/*@}*/

	/** \name real constants for material properties */
	/*@{*/
	/** add a material property. Return the index of the value in the parameters array */
	int AddProperty(const StringT& name, double value);
	
	/** array property names */
	const ArrayT<StringT>& PropertyNames(void) const { return fPropertyNames; };

	/** array property values */
	const ArrayT<double>&  PropertyValues(void) const { return fPropertyValues; };

	/** return property by name */
	double Property(const StringT& name) const { return fPropertyMap[name]; }

	/** return true of the property has been defined */
	bool HasProperty(const StringT& name) const { return fPropertyMap.HasKey(name); }
	/*@}*/

	/** \name symbols */
	/*@{*/
	/** associate a symbol with a given name. Return true if the symbol or false if the name
	 * was already assigned a symbol. */
	bool AddSymbol(const StringT& name, const StringT& symbol);

	/** return symbol by name */
	const StringT& Symbol(const StringT& name) const { return fSymbolMap[name]; }
	/*@}*/	

	/** comparison function to use for finding real constant. Need to override
	 * the default StringT::operator> and StringT::operator< because the strings
	 * passed from Fortran do no have C/C++ line endings */
	static int Compare(const MapNodeT<StringT, double>& tree_node, 
	                   const MapNodeT<StringT, double>& test_node);
	
private:

	/** unique identifier for material data */
	int fID;

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
	
	/** function to compute elastic moduli */
	Sierra_pc_elastic_moduli_func fPCFunction;
	/*@}*/
	
	/** list of XML command ids. The first command is the one defining the material
	 * itself */
	AutoArrayT<int> fXMLCommandID;
	 
	/** \name list of input variable names and offsets */
	/*@{*/
	AutoArrayT<StringT> fInputVariables;
	AutoArrayT<int>     fInputVariableSize;
	
	/** deformation measure determined from calls to SIERRA_Material_Data::AddInputVariable */
	InputVariableT fStrainMeasure;
	/*@}*/

	/** \name real constants for material properties */
	/*@{*/
	AutoArrayT<StringT> fPropertyNames;
	AutoArrayT<double> fPropertyValues;
	MapT<StringT, double> fPropertyMap;
	/*@}*/

	/** symbols */
	MapT<StringT, StringT> fSymbolMap;

private:

	/** keep track of unique ID's */
	static int sNextID;
};

/* add a material property */
inline int SIERRA_Material_Data::AddProperty(const StringT& name, double value)
{
	if (fPropertyNames.AppendUnique(name))
	{
		fPropertyValues.Append(value);
		fPropertyMap.Insert(name, value);
		return fPropertyValues.Length() - 1;
	}
	else 
		return fPropertyNames.PositionOf(name);
}

/* add a symbol */
inline bool SIERRA_Material_Data::AddSymbol(const StringT& name, const StringT& symbol) {
	return fSymbolMap.Insert(name, symbol);
}

} /* namespace Tahoe */

#endif /* _SIERRA_MAT_DATA_H_ */
