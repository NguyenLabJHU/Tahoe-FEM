#ifndef _PARADYN_T_H_
#define _PARADYN_T_H_

/* direct members */
#include "StringT.h"
#include "dArray2DT.h"

/* base class */
#include "ParticlePropertyT.h"

namespace Tahoe {

/* forward declarations */
class dArrayT;

class ParadynT: public ParticlePropertyT
{
public:

	/** constructor. Reads parameters from file and computes the
	 * coefficients of a cubic spline through the evenly spaced
	 * values of the potential,electron density read from the 
         * file. */
	ParadynT(const StringT& param_file);

	/** write properties to output */
	virtual void Write(ostream& out) const;

	/** \name EAM 
	 * Since the interaction functions will be called many times during
	 * a calculation, we do not want these functions to be virtual. Instead,
	 * we define some virtual functions that return pointers to functions
	 * that do the evaluation of the pair interactions. In this way, a single
	 * virtual function call is needed to return a pointer to a static
	 * member function, which is then called many times. */
	/*@{*/

	typedef double (*PairEnergyFunction)(double r_ab, double* data_a, double* data_b);
	typedef double (*EmbedEnergyFunction)(double rho_ab, double* data_a, double* data_b);
	typedef double (*EDEnergyFunction)(double r_ab, double* data_a, double* data_b);

	typedef double (*PairForceFunction)(double r_ab, double* data_a, double* data_b);
	typedef double (*EmbedForceFunction)(double rho_ab, double* data_a, double* data_b);
	typedef double (*EDForceFunction)(double r_ab, double* data_a, double* data_b);

	typedef double (*PairStiffnessFunction)(double r_ab, double* data_a, double* data_b);
	typedef double (*EmbedStiffnessFunction)(double rho_ab, double* data_a, double* data_b);
	typedef double (*EDStiffnessFunction)(double r_ab, double* data_a, double* data_b);
	

	/** \name return interaction functions */
	/*@{*/
	/** return a pointer to the energy function */
	ParadynT::PairEnergyFunction getPairEnergy(void);
	ParadynT::EmbedEnergyFunction getEmbedEnergy(void);
	ParadynT::EDEnergyFunction getElecDensEnergy(void);

	/** return a pointer to the force function */
	ParadynT::PairForceFunction getPairForce(void);
        ParadynT::EmbedForceFunction getEmbedForce(void);
	ParadynT::EDForceFunction getElecDensForce(void);

	/** return a pointer to the stiffness function */
	ParadynT::PairEnergyFunction getPairStiffness(void);
	ParadynT::EmbedEnergyFunction getEmbedStiffness(void);
	ParadynT::EDEnergyFunction getElecDensStiffness(void);

	/** return Paradyn-style coefficients table.
	 * returns false if no table is available. */
	virtual bool getParadynTable(const double** coeff, double& dr, 
		int& row_size, int& num_rows) const;
	/*@}*/

	/** the coefficients array */
	const dArray2DT& PairCoefficients(void) const { return fPairCoeff; };
	const dArray2DT& EmbedCoefficients(void) const { return fEmbedCoeff; };
	const dArray2DT& ElectronDensityCoefficients(void) const { return fElectronDensityCoeff; };

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

	/** coefficients of interpolant */
	dArray2DT fPairCoeff;
	dArray2DT fEmbedCoeff;
	dArray2DT fElectronDensityCoeff;
	
	
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

#endif /* _PARADYN_T_H_ */






