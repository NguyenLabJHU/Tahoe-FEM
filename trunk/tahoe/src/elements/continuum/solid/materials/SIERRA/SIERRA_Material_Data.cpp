/* $Id: SIERRA_Material_Data.cpp,v 1.2 2003-03-06 17:23:31 paklein Exp $ */
#include "SIERRA_Material_Data.h"

/* constructor */
SIERRA_Material_Data::SIERRA_Material_Data(const StringT& name, int XML_command_id, 
	int modulus_flag):
	fName(name),
	fModulusFlag(modulus_flag),
	fNumStateVars(0),
	fCheckFunction(NULL),
	fCalcFunction(NULL),
	fInitFunction(NULL)
{
	fXMLCommandID.Append(XML_command_id);
}

/* register input variable name */
void SIERRA_Material_Data::AddInputVariable(const StringT& input_var)
{
	int size = -1;
	if (input_var == "rot_strain_inc")
		size = 6;
	//others?
	
	/* not found */
	if (size == -1)
		ExceptionT::GeneralFail("SIERRA_Material_Data::AddInputVariable",
			"unrecognized variable: \"%s\"", input_var.Pointer());
	else if (fInputVariables.AppendUnique(input_var))
		fInputVariableSize.Append(size);
}
