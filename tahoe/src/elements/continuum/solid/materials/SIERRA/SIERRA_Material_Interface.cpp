/* $Id: SIERRA_Material_Interface.cpp,v 1.9 2004-08-08 02:02:57 paklein Exp $ */
#include "SIERRA_Material_Interface.h"
#include "SIERRA_Material_DB.h"
#include "SIERRA_Material_Data.h"

using namespace Tahoe;

const int kStringBufferSize = 255;
static char StringBuffer[kStringBufferSize];

/* register a named value */
void FORTRAN_NAME(register_real_constant)(double* value, const int* mat_vals, 
	const char* value_name, int value_name_len)
{
	/* fetch material */
	SIERRA_Material_Data* mat = SIERRA_Material_DB::Material(*mat_vals);
	
	/* add value - no check on success */
	f2c_string(value_name, value_name_len, StringBuffer, kStringBufferSize);
	mat->AddProperty(StringBuffer, *value);
}

void FORTRAN_NAME(get_real_constant)(double* destination, const int* mat_vals, 
	const char* value_name, int value_name_len)
{
	/* fetch material */
	SIERRA_Material_Data* mat = SIERRA_Material_DB::Material(*mat_vals);
	
	/* retrieve value */
	f2c_string(value_name, value_name_len, StringBuffer, kStringBufferSize);
	*destination = mat->Property(StringBuffer);
}

/* retrieve the index for the value for the given material */
void FORTRAN_NAME(get_var_index)(int* index, int* num_workset_elem, const char* variable_name, 
	const char* material_name, int variable_name_len, int material_name_len)
{
#pragma unused(variable_name_len)

	const char caller[] = "get_var_index";

	//NOTE: until we also get information about the number of integration
	//      points per element, we cannot figure out the offset into the
	//      variable array to the start of the given variable. Data is assumed
	//      to be stored across the workset as: 
	//		{{{a}_1, {a}_2, ..., {a}_nwsn}, {{b}_1, {b}_2, ..., {b}_nwsn}, ...}
	//
	//      where {a}_i = the variables for all integration points of element i
	if (*num_workset_elem != 1)
		ExceptionT::GeneralFail(caller, "expecting 1 workset element: %d", *num_workset_elem);

	/* fetch material */
	f2c_string(material_name, material_name_len, StringBuffer, kStringBufferSize);
	SIERRA_Material_Data* mat = SIERRA_Material_DB::Material(StringBuffer);
	
	/* input variables */
	const AutoArrayT<StringT>& vars = mat->InputVariables();
	const AutoArrayT<int>& sizes = mat->InputVariableSize();

	/* compute offset */
	*index = 1; /* FORTRAN numbering! */
	for (int i = 0; i < vars.Length(); i++)
		if (strncmp(vars[i], variable_name, vars[i].StringLength()) == 0)
			return; /* found */
		else
			*index += sizes[i];

	/* fail */
	ExceptionT::GeneralFail(caller, "variable not found: %s", StringBuffer);
}

/* register the material model */
void FORTRAN_NAME(register_material)(int* XML_command_id, Sierra_function_param_check check_func, 
	int* modulus_flag, const char* material_name, int material_name_len)
{
	f2c_string(material_name, material_name_len, StringBuffer, kStringBufferSize);

	/* initialize material record */
	SIERRA_Material_DB::InitMaterial(StringBuffer, *XML_command_id, *modulus_flag);

	/* set check function */
	SIERRA_Material_Data* mat = SIERRA_Material_DB::Material(StringBuffer);
	mat->SetCheckFunction(check_func);
}

/* register function to do material computations */
void FORTRAN_NAME(register_process_func)(Sierra_function_material_calc calc_func, 
	const char* material_name, int material_name_len)
{
	f2c_string(material_name, material_name_len, StringBuffer, kStringBufferSize);
	SIERRA_Material_Data* mat = SIERRA_Material_DB::Material(StringBuffer);
	mat->SetCalcFunction(calc_func);
}

/* register function to do material initialization */
void FORTRAN_NAME(register_init_func)(Sierra_function_material_init init_func, 
	const char* material_name, int material_name_len)
{
	f2c_string(material_name, material_name_len, StringBuffer, kStringBufferSize);
	SIERRA_Material_Data* mat = SIERRA_Material_DB::Material(StringBuffer);
	mat->SetInitFunction(init_func);
}

/* register function to compute tangent moduli */
void FORTRAN_NAME(register_pc_elastic_moduli_func)(Sierra_pc_elastic_moduli_func pc_func, 
	const char* func_name, int func_name_len)
{
	f2c_string(func_name, func_name_len, StringBuffer, kStringBufferSize);
	SIERRA_Material_Data* mat = SIERRA_Material_DB::Material(StringBuffer);
	mat->SetPCFunction(pc_func);
}

/* register function evaluation */
void FORTRAN_NAME(register_func_eval)(const int* matvals, const char* func_name, int func_name_len)
{
#pragma unused(matvals)
#pragma unused(func_name)
#pragma unused(func_name_len)
//NOTE: nothing to do. Function information is retrieved during the call to
//     apub_fortran_fctn_eval.
}

/* register the number of state variables */
void FORTRAN_NAME(register_num_state_vars)(int* nsv, 
	const char* material_name, int material_name_len)
{
	f2c_string(material_name, material_name_len, StringBuffer, kStringBufferSize);
	SIERRA_Material_Data* mat = SIERRA_Material_DB::Material(StringBuffer);
	mat->SetNumStateVariables(*nsv);
}

/* register the data */
void FORTRAN_NAME(register_input_var)(const char* variable_name, const char* material_name,
	int variable_name_len, int material_name_len)
{
	f2c_string(material_name, material_name_len, StringBuffer, kStringBufferSize);
	SIERRA_Material_Data* mat = SIERRA_Material_DB::Material(StringBuffer);
	f2c_string(variable_name, variable_name_len, StringBuffer, kStringBufferSize);
	mat->AddInputVariable(StringBuffer);
}

/* register the XML commands that specify material parameters */
void FORTRAN_NAME(register_parser_line)(int* XML_command_id, const char* material_name, int material_name_len)
{
	f2c_string(material_name, material_name_len, StringBuffer, kStringBufferSize);
	SIERRA_Material_Data* mat = SIERRA_Material_DB::Material(StringBuffer);
	mat->AddXMLCommandID(*XML_command_id);
}

/* evaluate function */
void FORTRAN_NAME(apub_fortran_fctn_eval)(const int* matvals, const double* arg, double* out,
	const char* func_name, int func_name_len)
{
	/* fetch material */
	SIERRA_Material_Data* mat = SIERRA_Material_DB::Material(*matvals);

	/* name passed in is key to function's symbol */
	f2c_string(func_name, func_name_len, StringBuffer, kStringBufferSize);

	/* get symbol */
	const StringT& function_symbol = mat->Symbol(StringBuffer);

	/* evaluate function */
	*out = SIERRA_Material_DB::Evaluate(function_symbol, *arg);
}

/* error reporting */
void FORTRAN_NAME(report_error)(int* code, const char* error_string, int error_string_len)
{
	f2c_string(error_string, error_string_len, StringBuffer, kStringBufferSize);
	ExceptionT::GeneralFail("SIERRA_Material_Interface::report_error", "code %d: \"%s\"",
		*code, StringBuffer);
}

/* convert a fortran character array into a C string */
void f2c_string(const char* f_string, int f_string_len, char* buffer, int buffer_len)
{
	/* check buffer size */
	if (f_string_len >= buffer_len - 1)
		ExceptionT::GeneralFail("f2c_string", "insufficient buffer size: %d >= %d",
			f_string_len, buffer_len - 1);

	/* copy string */
	memcpy(buffer, f_string, f_string_len*sizeof(char));
	
	/* terminate */
	buffer[f_string_len] = '\0';
}
