/* $Id: ParadynEAMT.h,v 1.4 2004-04-09 02:03:05 hspark Exp $ */
#ifndef _PARADYN_EAM_T_H_
#define _PARADYN_EAM_T_H_

/* direct members */
#include "StringT.h"
#include "dArray2DT.h"

/* base class */
#include "EAMPropertyT.h"

namespace Tahoe {

/* forward declarations */
class dArrayT;

class ParadynEAMT: public EAMPropertyT
{
public:

	/** constructor. Reads parameters from file and computes the
	 * coefficients of a cubic spline through the evenly spaced
	 * values of the potential,electron density read from the 
         * file. */
	ParadynEAMT(const StringT& param_file);

	/** write properties to output */
	virtual void Write(ostream& out) const;

	/** \name return interaction functions */
	/*@{*/
	/** return a pointer to the energy function */
	virtual PairEnergyFunction getPairEnergy(void);
	virtual EmbedEnergyFunction getEmbedEnergy(void);
	virtual EDEnergyFunction getElecDensEnergy(void);

	/** return a pointer to the force function */
	virtual PairForceFunction getPairForce(void);
	virtual EmbedForceFunction getEmbedForce(void);
	virtual EDForceFunction getElecDensForce(void);

	/** return a pointer to the stiffness function */
	virtual PairStiffnessFunction getPairStiffness(void);
	virtual EmbedStiffnessFunction getEmbedStiffness(void);
	virtual EDStiffnessFunction getElecDensStiffness(void);

	/** return Paradyn-style coefficients table.
	 * returns false if no table is available. */
	virtual bool getParadynTable(const double** coeff, double& dr, 
		int& row_size, int& num_rows) const;
	/*@}*/

	/** the coefficients array */
	const dArray2DT& PairCoefficients(void) const { return fPairCoeff; };
	const dArray2DT& EmbedCoefficients(void) const { return fEmbedCoeff; };
	const dArray2DT& ElectronDensityCoefficients(void) const { return fElectronDensityCoeff; };

	/** add accessor function for lattice parameter */
	virtual double GetLatticeParameter(void) const { return fLatticeParameter; };

private:

	/** \name interaction functions */
	/*@{*/
	static double PairEnergy(double r_ab, double* data_a, double* data_b);
	static double EmbeddingEnergy(double rho_ab, double* data_a, double* data_b);
	static double ElecDensEnergy(double r_ab, double* data_a, double* data_b);

	static double PairForce(double r_ab, double* data_a, double* data_b);
	static double EmbeddingForce(double rho_ab, double* data_a, double* data_b);
	static double ElecDensForce(double r_ab, double* data_a, double* data_b);

	static double PairStiffness(double r_ab, double* data_a, double* data_b);
	static double EmbeddingStiffness(double rho_ab, double* data_a, double* data_b);
	static double ElecDensStiffness(double r_ab, double* data_a, double* data_b);
	/*@}*/

	/** \name auxiliary functions */
	static double EnergyAux(double r_ab,int n, double inc, double* coeff);
	static double ForceAux(double r_ab,int n, double inc, double* coeff);
	static double StiffnessAux(double r_ab,int n, double inc, double* coeff);


	/** compute the coefficients. Translated from the Paradyn routine interpolate.F
	 * by Steve Plimpton, SNL-NM.
	 * \param f function values at evenly spaced intervals 
	 * \param dx intervals between data points
	 * \param coeff destination of coefficient. Allocated during the function call. */
	static void ComputeCoefficients(const ArrayT<double>& f, double dx, dArray2DT& coeff);

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
	double f_inc;
	/** 1/rho */
	double rho_inc;
	/*@}*/

	/** \name coefficients of interpolant */
	/*@{*/
	dArray2DT fPairCoeff;
	dArray2DT fEmbedCoeff;
	dArray2DT fElectronDensityCoeff;
	/*@}*/

	/** \name static parameters 
	 * There parameters are used during evaluation of the static interaction 
	 * functions. These are copied to static when a function pointer is requested */
	/*@{*/
	static int     s_nr;
	static double  s_f_inc;
	static double* s_Paircoeff;


	static int     s_np;
	static double  s_e_inc;
	static double* s_Embcoeff;


	static double* s_ElecDenscoeff;

	/*@}*/
};

} /* namespace Tahoe */

#endif /* _PARADYN_EAM_T_H_ */
