/* $Id: SIERRA_HypoElastic.c,v 1.9 2003-12-28 08:21:22 paklein Exp $ */
#include "SIERRA_Material_Interface.h"
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* function prototypes */
void SIERRA_HypoElastic_reg(void);

/* parameter check function */
void SIERRA_HypoElastic_check(int* matvals);

/* state variable init function */
void SIERRA_HypoElastic_init(int* nelem, double* dt, int* nsv, 
	double* state_old, double* state_new, int* matvals, int* ncd);

/* calculation function */
void SIERRA_HypoElastic_calc(int* nelem, double* dt,
	double* vars_input, int* ivars_size, double* stress_old, double* stress_new, 
	int* nsv, double* state_old, double* state_new, int* matvals, int* ncd);

/*******/

void SIERRA_HypoElastic_reg(void)
{
	const char model_name[] = "HYPOELASTIC";
	const char var_name[] = "rot_strain_increment";
	int num_state = 1;
	int iflag = 1;
	int material_code = 9999;

	/* register the material model */
	FORTRAN_NAME(register_material)(&material_code, SIERRA_HypoElastic_check, &iflag, model_name, strlen(model_name));

	/* register function to do material computations */
	FORTRAN_NAME(register_process_func)(SIERRA_HypoElastic_calc, model_name, strlen(model_name));

	/* register function to do material initialization */
	FORTRAN_NAME(register_init_func)(SIERRA_HypoElastic_init, model_name, strlen(model_name));

	/* register the number of state variables that the model needs */
	FORTRAN_NAME(register_num_state_vars)(&num_state, model_name, strlen(model_name));

	/* register the data that the material model needs from the element in
	 * order to do its computations */
	FORTRAN_NAME(register_input_var)(var_name, model_name, strlen(var_name), strlen(model_name));
}

void SIERRA_HypoElastic_check(int* matvals)
{
	/* fetch material properties */
	double bulk_modulus, two_mu;
	FORTRAN_NAME(get_real_constant)(&bulk_modulus, matvals, "BULK_MODULUS", strlen("BULK_MODULUS"));
	FORTRAN_NAME(get_real_constant)(&two_mu, matvals,"TWO_MU", strlen("TWO_MU"));
	
	if (bulk_modulus < 0.0 || two_mu < 0.0) {
		printf("{kappa, 2 mu} = {%g, %g}\n", bulk_modulus, two_mu);
		abort();
	}
}

/* function to do material initialization */
void SIERRA_HypoElastic_init(int* nelem, double* dt, int* nsv, 
	double* state_old, double* state_new, int* matvals, int* ncd)
{
#pragma unused(dt)
#pragma unused(matvals)
#pragma unused(ncd)

	int i, j;
	for (j = 0; j < *nelem; j++)
	{
		for (i = 0; i < *nsv; i++)
		{
			state_old[i] = 0.0;
			state_new[i] = 0.0;
		}
		state_old += *nsv;
		state_new += *nsv;
	}
}

/* function to do material computations */
void SIERRA_HypoElastic_calc(int* nelem, double* dt,
	double* vars_input, int* ivars_size, double* stress_old, double* stress_new, 
	int* nsv, double* state_old, double* state_new, int* matvals, int* ncd)
{
#pragma unused(dt)
#pragma unused(ncd)

	const char model_name[] = "HYPOELASTIC";
	int i, istrain;
	double bulk_modulus = 0.0;
	double two_mu = 0.0;
	double* dstrain = NULL;
	double kappa_minus_2_mu_by_3, kappa_plus_4_mu_by_3;
	
	/* fetch material properties */
	FORTRAN_NAME(get_real_constant)(&bulk_modulus, matvals, "BULK_MODULUS", strlen("BULK_MODULUS"));
	FORTRAN_NAME(get_real_constant)(&two_mu, matvals, "TWO_MU", strlen("TWO_MU"));
	kappa_minus_2_mu_by_3 = bulk_modulus - two_mu/3.0;
	kappa_plus_4_mu_by_3 = bulk_modulus + 2.0*two_mu/3.0;

	/* offset to strain variable - starts at 1 by Fortran numbering convention */
	FORTRAN_NAME(get_var_index)(&istrain, nelem, "rot_strain_increment", model_name,
		strlen("rot_strain_increment"), strlen(model_name));
	dstrain = vars_input + istrain - 1;
	
	/* loop over stress points */
	for (i = 0; i < *nelem; i++)
	{
		/* compute updated stress */
		stress_new[0] = stress_old[0] + dstrain[0]*kappa_plus_4_mu_by_3 + (dstrain[1] + dstrain[2])*kappa_minus_2_mu_by_3;
		stress_new[1] = stress_old[1] + dstrain[1]*kappa_plus_4_mu_by_3 + (dstrain[0] + dstrain[2])*kappa_minus_2_mu_by_3;
		stress_new[2] = stress_old[2] + dstrain[2]*kappa_plus_4_mu_by_3 + (dstrain[0] + dstrain[1])*kappa_minus_2_mu_by_3;

		stress_new[3] = stress_old[3] + two_mu*dstrain[3];
		stress_new[4] = stress_old[4] + two_mu*dstrain[4];
		stress_new[5] = stress_old[5] + two_mu*dstrain[5];
		
		/* compute updated strain energy density */
		state_new[0] = state_old[0] 
			+ stress_new[0]*dstrain[0] + stress_new[1]*dstrain[1] + stress_new[2]*dstrain[2]
       + 2.0*(stress_new[3]*dstrain[3] + stress_new[4]*dstrain[4] + stress_new[5]*dstrain[5]);

		/* next stress point */
		dstrain += *ivars_size;
		stress_old += 6;
		stress_new += 6;
		state_old += *nsv;
		state_new += *nsv;
	}
}

#ifdef __cplusplus
} /* extern "C" */
#endif
