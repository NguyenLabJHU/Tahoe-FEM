/* $Id: SIERRA_Material_Data.cpp,v 1.8 2004-08-08 02:02:57 paklein Exp $ */
#include "SIERRA_Material_Data.h"

using namespace Tahoe;

/* static data */
int SIERRA_Material_Data::sNextID = 0;

/* constructor */
SIERRA_Material_Data::SIERRA_Material_Data(const StringT& name, int XML_command_id, 
	int modulus_flag):
	fID(sNextID++),
	fName(name),
	fModulusFlag(modulus_flag),
	fNumStateVars(0),
	fCheckFunction(NULL),
	fCalcFunction(NULL),
	fInitFunction(NULL),
	fPCFunction(NULL),
	fStrainMeasure(kundefined)
{
	fXMLCommandID.Append(XML_command_id);
//	fPropertyMap.SetCompareFunction(SIERRA_Material_Data::Compare);
// NOTE: not needed now that all fortran strings should be translated to C strings
//       before calling SIERRA_Material_Data methods.
}

/* register input variable name */
void SIERRA_Material_Data::AddInputVariable(const StringT& input_var)
{
	const char caller[] = "SIERRA_Material_Data::AddInputVariable";

	int size = -1;
	if (input_var == "rot_strain_inc") {
		if (fStrainMeasure != kundefined)
			ExceptionT::GeneralFail(caller, "strain type is already set"); 
		fStrainMeasure = krot_strain_inc;
		size = 6;
	}
	else if (input_var == "rot_strain_increment") {
		if (fStrainMeasure != kundefined)
			ExceptionT::GeneralFail(caller, "strain type is already set"); 
		fStrainMeasure = krot_strain_increment;
		size = 6;
	}
	else if (input_var == "temperature_old")
		size = 1;
	else if (input_var == "temperature_new")
		size = 1;
	
	/* not found */
	if (size == -1)
		ExceptionT::GeneralFail(caller, "unrecognized variable: \"%s\"", input_var.Pointer());
	else if (fInputVariables.AppendUnique(input_var))
		fInputVariableSize.Append(size);
}

/* compare function that ignores the actual length of the test_value */
int SIERRA_Material_Data::Compare(
	const MapNodeT<StringT, double>& tree_node,
	const MapNodeT<StringT, double>& test_node)
{
	const StringT& tree_value = tree_node.Key();
	const StringT& test_value = test_node.Key();

	int len = tree_value.StringLength();
	return strncmp(tree_value, test_value, len);
};
