/* $Id: SIERRA_Material_Interface.cpp,v 1.1 2003-03-06 17:23:31 paklein Exp $ */
#include "SIERRA_Material_Interface.h"
#include "SIERRA_Material_DB.h"
#include "SIERRA_Material_Data.h"

/* retrieve a named value */
void get_real_constant(double* destination, const double* all_values, 
	const char* value_name)
{
	int index = SIERRA_Material_DB::RealIndex(value_name);
	*destination = all_values[index];
}

/* retrieve the index for the value for the given material */
void get_var_index(int* index, int* num_workset_elem, const char* variable_name, 
	const char* material_name)
{
	const char caller[] = "get_var_index";

	SIERRA_Material_Data* mat = SIERRA_Material_DB::Material(material_name);
	
	/* input variables */
	const AutoArrayT<StringT>& vars = mat->InputVariables();
	const AutoArrayT<int>& sizes = mat->InputVariableSize();

	if (vars.Length() != 1)
		ExceptionT::GeneralFail(caller, "assuming only 1 variable: %d", vars.Length());
	//NOTE: until we also get information about the number of integration
	//      points per element, we cannot figure out the offset into the
	//      variable array to the start of the given variable. Data is assumed
	//      to be stored across the workset as: 
	//		{{{a}_1, {a}_2, ..., {a}_nwsn}, {{b}_1, {b}_2, ..., {b}_nwsn}, ...}
	//
	//      where {a}_i = the variables for all integration points of element i
#pragma unused(num_workset_elem)
	
	if (vars[0] != variable_name)
		ExceptionT::GeneralFail(caller, "variable not found: %s", variable_name);
	else
		*index = sizes[0];
}

/* register the material model */
void register_material(int* XML_command_id, Sierra_function_param_check check_func, 
	int* modulus_flag, const char* material_name)
{
	/* initialize material record */
	SIERRA_Material_DB::InitMaterial(material_name, *XML_command_id, *modulus_flag);

	/* set check function */
	SIERRA_Material_Data* mat = SIERRA_Material_DB::Material(material_name);
	mat->SetCheckFunction(check_func);
}

/* register function to do material computations */
void register_process_func(Sierra_function_material_calc calc_func, const char* material_name)
{
	SIERRA_Material_Data* mat = SIERRA_Material_DB::Material(material_name);
	mat->SetCalcFunction(calc_func);
}

/* register function to do material initialization */
void register_init_func(Sierra_function_material_init init_func, const char* material_name)
{
	SIERRA_Material_Data* mat = SIERRA_Material_DB::Material(material_name);
	mat->SetInitFunction(init_func);
}

/* register the number of state variables */
void register_num_state_vars(int* nsv, const char* material_name)
{
	SIERRA_Material_Data* mat = SIERRA_Material_DB::Material(material_name);
	mat->SetNumStateVariables(*nsv);
}

/* register the data */
void register_input_var(const char* variable_name, const char* material_name)
{
	SIERRA_Material_Data* mat = SIERRA_Material_DB::Material(material_name);
	mat->AddInputVariable(variable_name);
}

/* register the XML commands that specify material parameters */
void register_parser_line(int* XML_command_id, const char* material_name)
{
	SIERRA_Material_Data* mat = SIERRA_Material_DB::Material(material_name);
	mat->AddXMLCommandID(*XML_command_id);
}
