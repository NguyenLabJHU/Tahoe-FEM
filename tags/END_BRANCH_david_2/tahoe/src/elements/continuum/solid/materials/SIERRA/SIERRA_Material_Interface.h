/* $Id: SIERRA_Material_Interface.h,v 1.8 2004-08-08 02:02:57 paklein Exp $ */
#ifndef __SIERRA_MAT_INTERFACE_H__
#define __SIERRA_MAT_INTERFACE_H__

/* Fortran to C/C++ name mangling */
#include "fortran_names.h"

#ifdef __cplusplus 
extern "C" {
#endif

/** \name retrieving parameter information */
/*@{*/ 
/** register a named value */
extern void FORTRAN_NAME(register_real_constant)(double* value, const int* mat_vals, 
	const char* value_name, int value_name_len);

/** retrieve a named value */
extern void FORTRAN_NAME(get_real_constant)(double* destination, const int* mat_vals, 
	const char* value_name, int value_name_len);

/** retrieve the index for the value for the given material. The numbering uses
 * Fortran conventions, so the first index is 1. */
extern void FORTRAN_NAME(get_var_index)(int* index, int* num_workset_elem, const char* variable_name, 
	const char* material_name, int variable_name_len, int material_name_len);
/*@}*/ 

/** \name Sierra callback function types */
/*@{*/
/** function to do checking/registration of parameters. This function
 * is registered with a call to register_material.
 * \param parameters array of material parameters
 */
typedef void (*Sierra_function_param_check)(int* mat_vals);

/** function to do material computations
 * \param nelem number of elements in the workset
 * \param dt time step
 * \param vars_input array of driving input variables, i.e., strain, temperature, etc.
 *        The array is dimension [nelem] x [ivars_size]    
 * \param ivars_size number of input variables per stress point
 * \param stress_old (rotated) stress from the previous time increment
 * \param stress_new destination for current stress
 * \param isize_state number of state variables per stress point
 * \param state_old state variables from the previous time increment
 * \param state_new destination for updated state variables
 * \param matvals array of material parameters
 */
typedef void (*Sierra_function_material_calc)(int* nelem, double* dt,
	double* vars_input, int* ivars_size, 
	double* stress_old, double* stress_new, 
	int* isize_state, double* state_old, double* state_new, 
	int* matvals);

/** function to do material initialization
 * \param nelem number of elements in the workset
 * \param dt time step
 * \param vars_input array of driving input variables, i.e., strain, temperature, etc.
 *        The array is dimension [nelem] x [ivars_size]    
 * \param ivars_size number of input variables per stress point
 * \param isize_state number of state variables
 * \param state_old state variables from the previous time increment
 * \param state_new destination for updated state variables
 * \param matvals array of material parameters
 */
typedef void (*Sierra_function_material_init)(int* nelem, double* dt, 
	double* vars_input, int* ivars_size,
	int* isize_state, double* state_old, double* state_new, 
	int* matvals);

/** function to compute tangent moduli (not used) */
typedef void (*Sierra_pc_elastic_moduli_func)(void);
/*@}*/

/** \name registration functions */
/*@{*/
/** register the material model.
 * \param XML_command_id XML command id for the material parameter block
 * \param function to do checking/registration of parameters
 * \param modulus_flag  0/1 = dont/do complete elastic constants
 * \param material_name name for material model
*/
extern void FORTRAN_NAME(register_material)(int* XML_command_id, Sierra_function_param_check check_func, 
	int* modulus_flag, const char* material_name, int material_name_len);

/** register function to do material computations */
extern void FORTRAN_NAME(register_process_func)(Sierra_function_material_calc calc_func, const char* material_name, int material_name_len);

/** register function to do material initialization */
extern void FORTRAN_NAME(register_init_func)(Sierra_function_material_init init_func, const char* material_name, int material_name_len);

/** register the number of state variables */
extern void FORTRAN_NAME(register_num_state_vars)(int* nsv, const char* material_name, int material_name_len);

/** register the data that the material model needs from the element in
 * order to do its computations */
extern void FORTRAN_NAME(register_input_var)(const char* variable_name, const char* material_name, int variable_name_len, int material_name_len);

/** register the XML commands that specify material parameters */
extern void FORTRAN_NAME(register_parser_line)(int* XML_command_id, const char* material_name, int material_name_len);

/** register function to compute tangent moduli */
extern void FORTRAN_NAME(register_pc_elastic_moduli_func)(Sierra_pc_elastic_moduli_func pc_func, const char* func_name, int func_name_len);

/** register function evaluation */
extern void FORTRAN_NAME(register_func_eval)(const int* matvals, const char* func_name, int func_name_len);
/*@}*/

/** evaluate function */
extern void FORTRAN_NAME(apub_fortran_fctn_eval)(const int* matvals, const double* arg, double* out,
	const char* func_name, int func_name_len);

/** error reporting */
extern void FORTRAN_NAME(report_error)(int* code, const char* error_string, int error_string_len);

/** convert a fortran character array into a C string */
extern void f2c_string(const char* f_string, int f_string_len, char* buffer, int buffer_len);

#ifdef __cplusplus 
}
#endif

#endif /* __SIERRA_MAT_INTERFACE_H__ */
