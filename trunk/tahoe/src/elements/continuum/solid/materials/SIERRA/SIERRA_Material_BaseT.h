/* $Id: SIERRA_Material_BaseT.h,v 1.5 2004-07-15 08:27:31 paklein Exp $ */
#ifndef _SIERRA_MAT_BASE_T_H_
#define _SIERRA_MAT_BASE_T_H_

/* base classes */
#include "FSIsotropicMatT.h"

/* direct members */
#include "StringT.h"
#include "iArrayT.h"
#include "dArray2DT.h"

namespace Tahoe {

/* forward declarations */
class SpectralDecompT;
class SIERRA_Material_Data;
class ParameterListT;

/** base class for wrappers around Sierra material models. Sub-classes \e must 
 * overload two methods: SIERRA_Material_BaseT::Register_SIERRA_Material, which
 * should register all material parameters with SIERRA_Material_DB, and
 * SIERRA_Material_BaseT::SetOutputVariables, which defines the labels and indecies 
 * for values in the state variable array that will be written with material output. */
class SIERRA_Material_BaseT: public FSIsotropicMatT
{
public:

	/** constructor */
	SIERRA_Material_BaseT(ifstreamT& in, const FSMatSupportT& support);

	/** destructor */
	~SIERRA_Material_BaseT(void);

	/** required parameter flags */
	virtual bool Need_F_last(void) const { return true; };

	/** form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/** \name initialization */
	/*@{*/
	/** initialization. Called immediately after constructor to allow
	 * class specific initializations. */
	virtual void Initialize(void);

	/** return true if model needs ContinuumMaterialT::PointInitialize
	 * to be called for every integration point of every element as
	 * part of the model initialization. \return false by default. */
	virtual bool NeedsPointInitialization(void) const;
	
	/** model initialization. Called per integration point for every
	 * element using the model. Deformation variables are available
	 * during this call. */
	virtual void PointInitialize(void);
	/*@}*/

	/** update internal variables. Called once per element for all
	 * elements using the model, hence no deformation variables are
	 * available during this call. */
	virtual void UpdateHistory(void);

	/** restore internal variables to their state at the beginning of
	 * the current time increment. Called once per element for all
	 * elements using the model, hence no deformation variables are
	 * available during this call. */
	virtual void ResetHistory(void);

	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. The value is not guaranteed to
	 * persist during intervening calls to any other non-const
	 * accessor. \return 1/3 of the trace of the three-dimensional
	 * stress tensor, regardless of the dimensionality of the
	 * problem. */
	virtual double Pressure(void) const;
	/*@}*/

	/** \name material description */
	/*@{*/
	/** material tangent moduli */
	virtual const dMatrixT& C_IJKL(void);

	/** 2nd Piola-Kirchhoff stress */
	virtual const dSymMatrixT& S_IJ(void);
	/*@}*/

	/** returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);

	/** \name material output */
	/*@{*/
	/* returns the number of variables computed for nodal extrapolation
	 * during for element output, i.e., some of the state variables. */
	virtual int NumOutputVariables(void) const;
	
	/** return the labels for the output variables */
	virtual void OutputLabels(ArrayT<StringT>& labels) const;

	/** compute the output variables */
	virtual void ComputeOutput(dArrayT& output);
	/*@}*/

protected:

	/** \name required by sub-classes */
	/*@{*/
	/** call the SIERRA registration function */
	virtual void Register_SIERRA_Material(void) const = 0;

	/** set material output */
	virtual void SetOutputVariables(iArrayT& variable_index,
		ArrayT<StringT>& output_labels) const = 0;
	/*@}*/

	/** \name load/store element data */
	/*@{*/
	void Load(ElementCardT& element, int ip);
	void Store(ElementCardT& element, int ip);
	/*@}*/

private:

	/** \name conversion functions */
	/*@{*/
	void SIERRA_to_dSymMatrixT(const double* pA, dSymMatrixT& B) const;
	void dSymMatrixT_to_SIERRA(const dSymMatrixT& A, double* pB) const;
	/*@}*/

	/** compute strains, rotated stresses, etc. */
	void Set_Calc_Arguments(void);

	/** \name parameter handling */
	/*@{*/
	/** read input parameters. Read until encountering a line beginning with
	 * "end" and ending with the name of the ParameterListT. Nested "begin" is 
	 * are processed recursively. Parameters must follow the pattern
	 *    value_name = value
	 */
	void Read_SIERRA_Input(ifstreamT& in, ParameterListT& param_list) const;

	SIERRA_Material_Data* Process_SIERRA_Input(ParameterListT& param_list);
	/*@}*/

protected:

	/** \name Sierra_function_material_calc arguments */
	/*@{*/
	nArrayT<double> fdstran;     /**< rotated strain increment */
	nArrayT<double> fstress_old; /**< (rotated) stress from the previous time increment */
	nArrayT<double> fstress_new; /**< destination for updated stress */
	nArrayT<double> fstate_old;  /**< state variables from the previous time increment */
	nArrayT<double> fstate_new;  /**< destination for updated state variables */
	nArrayT<double> fmatvals;    /**< array of material parameters */
	/*@}*/

private:

	/** dump debug info */
	int fDebug;

	/** material name */
	StringT fMaterialName;
	
	/** model name */
	StringT fMaterialModelName;

	/** tangent type */
	GlobalT::SystemTypeT fTangentType;
	
	/** material data card */
	SIERRA_Material_Data* fSIERRA_Material_Data;

	/* work space */
	dMatrixT    fModulus;            // return value
	dSymMatrixT fStress;             // return value
	double fPressure; /**< pressure for the most recent calculation of the stress */
	
	/** \name polar decomposition work space */
	/*@{*/
	SpectralDecompT* fDecomp;
	dMatrixT fF_rel;
	dMatrixT fA_nsd;
	dSymMatrixT fU1, fU2, fU1U2;
	/*@}*/

	/** \name material output data */
	/*@{*/
	iArrayT fOutputIndex;
	ArrayT<StringT> fOutputLabels;
	/*@}*/

	/** calc function argument array storage */
	/*@{*/
	nArrayT<double> fArgsArray;

	/** storage block size (per ip) */
	int fBlockSize;
	/*@}*/

	/** number of SIERRA_Material_BaseT instances.
	 * SIERRA_Material_BaseT::sSIERRA_Material_DB is constructed when the
	 * first one is instantiated and is deleted when the last one is freed. */
	static int sSIERRA_Material_count;
};

/* inlines */
inline GlobalT::SystemTypeT SIERRA_Material_BaseT::TangentType(void) const
{
	return fTangentType;
}

} /* namespace Tahoe */

#endif /* _SIERRA_MAT_BASE_T_H_ */
