/* $Id: ParadynPairT.h,v 1.3 2002-12-05 07:07:45 paklein Exp $ */
#ifndef _PARADYN_PAIR_T_H_
#define _PARADYN_PAIR_T_H_

/* base class */
#include "PairPropertyT.h"
#include "ParadynT.h"

/* direct members */
#include "StringT.h"
#include "dArray2DT.h"

namespace Tahoe {

/* forward declarations */
class dArrayT;

/** pair interaction constructed for Paradyn (EAM) format. The
 * potential parameters are read from the given file, from which
 * a cubic spline is calculated which is then used to evaluate
 * the potential and its derivatives. */
class ParadynPairT: public PairPropertyT, protected ParadynT
{
public:

	/** constructor. Reads parameters from file and computes the
	 * coefficients of a cubic spline through the evenly spaced
	 * values of the potential read from the file. */
	ParadynPairT(const StringT& param_file);

	/** write properties to output */
	virtual void Write(ostream& out) const;

	/** \name return interaction functions */
	/*@{*/
	/** return a pointer to the energy function */
	virtual EnergyFunction getEnergyFunction(void);

	/** return a pointer to the force function */
	virtual ForceFunction getForceFunction(void);

	/** return a pointer to the stiffness function */
	virtual StiffnessFunction getStiffnessFunction(void);

	/** return Paradyn-style coefficients table.
	 * returns false if no table is available. */
	virtual bool getParadynTable(const double** coeff, double& dr, 
		int& row_size, int& num_rows) const;
	/*@}*/

	/** the coefficients array */
	const dArray2DT& Coefficients(void) const { return fCoefficients; };

private:

	/** \name interaction functions */
	/*@{*/
	static double Energy(double r_ab, double* data_a, double* data_b);
	static double Force(double r_ab, double* data_a, double* data_b);
	static double Stiffness(double r_ab, double* data_a, double* data_b);
	/*@}*/

private:

	/** path to source file */
	StringT fParams;
	
	/** description from parameters file */
	StringT fDescription;
	
	/** \name Paradyn file parameters */
	/*@{*/
	int fAtomicNumber;
	double fLatticeParameter;
	StringT fStructure;	
	
	/** cut-off distance */
	double f_cut;
	
	/** 1/dr */
	double f_1bydr;
	/*@}*/

	/** coefficients of interpolant */
	dArray2DT fCoefficients;
	
	/** \name static parameters 
	 * There parameters are use during evaluation of the static interaction 
	 * functions. These are copied to static when a function pointer is requested */
	/*@{*/
	static int     s_nr;
	static double  s_1bydr;
	static double* s_coeff;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _HARMONIC_PAIR_T_H_ */
