/* $Id: ParadynPairT.h,v 1.6 2003-10-28 23:31:52 paklein Exp $ */
#ifndef _PARADYN_PAIR_T_H_
#define _PARADYN_PAIR_T_H_

/* base class */
#include "PairPropertyT.h"

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
class ParadynPairT: public PairPropertyT
{
public:

	/** constructor. Reads parameters from file and computes the
	 * coefficients of a cubic spline through the evenly spaced
	 * values of the potential read from the file. */
	ParadynPairT(const StringT& param_file);
	ParadynPairT(void);

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

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/** \name interaction functions */
	/*@{*/
	static double Energy(double r_ab, double* data_a, double* data_b);
	static double Force(double r_ab, double* data_a, double* data_b);
	static double Stiffness(double r_ab, double* data_a, double* data_b);
	/*@}*/

	/** compute the coefficients. Translated from the Paradyn routine interpolate_pair.F
	 * by Steve Plimpton, SNL-NM.
	 * \param f function values at evenly spaced intervals 
	 * \param dx intervals between data points
	 * \param coeff destination of coefficient. Allocated during the function call. */
	static void ComputeCoefficients(const ArrayT<double>& f, double dx, dArray2DT& coeff);

	/** initialize potential from stream */
	void ReadParameters(const StringT& param_file);

private:

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
	 * There parameters are used during evaluation of the static interaction 
	 * functions. These are copied to static when a function pointer is requested */
	/*@{*/
	static int     s_nr;
	static double  s_1bydr;
	static double* s_coeff;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _PARADYN_PAIR_T_H_ */



