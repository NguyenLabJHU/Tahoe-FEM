/* $Id: SolidMaterialT.h,v 1.1.1.1.2.5 2001-07-02 21:54:37 paklein Exp $ */
/* created: paklein (11/20/1996)                                          */
/* Defines the interface for elastic continuum materials.                 */

#ifndef _STRUCTURAL_MATERIALT_H_
#define _STRUCTURAL_MATERIALT_H_

#include "GlobalT.h"

/* base class */
#include "ContinuumMaterialT.h"

/* direct members */
#include "dMatrixT.h"

/* forward declarations */
class ifstreamT;
class ElementBaseT;
class dMatrixT;
class ThermalDilatationT;
class LoadTime;
class dSymMatrixT;
class LocalArrayT;
class ElasticT;

/** base class for constitutive models for solids */
class SolidMaterialT: public ContinuumMaterialT
{
public:

	/** constructor */
	SolidMaterialT(ifstreamT& in, const ContinuumElementT& element);

	/** destructor */
	~SolidMaterialT(void);

	/** write parameters */
	virtual void Print(ostream& out) const;

	/* spatial description */
	virtual const dMatrixT& c_ijkl(void) = 0;  /**< spatial tangent modulus */
	virtual const dSymMatrixT& s_ij(void) = 0; /**< Cauchy stress */

	/* material description */
	virtual const dMatrixT& C_IJKL(void) = 0;  /**< material tangent moduli */
	virtual const dSymMatrixT& S_IJ(void) = 0; /**< 2nd Piola-Kirchhoff stress */

	/** strain energy density */
	virtual double StrainEnergyDensity(void) = 0;

	/** acoustical tensor.
	 * \param normal wave propagation direction
	 * \return acoustical tensor */
	virtual const dSymMatrixT& AcousticalTensor(const dArrayT& normal) = 0;

	/** acoustic wave speeds.
	 * \param normal wave propagation direction
	 * \param speeds the computed acoustic wave speeds */
	void WaveSpeeds(const dArrayT& normal, dArrayT& speeds);

	/* required parameter flags - all false by default */
	virtual bool NeedDisp(void) const     { return false; };
	virtual bool NeedLastDisp(void) const { return false; };
	virtual bool NeedVel(void) const      { return false; };

	/* returns true if the material has internal forces in the unloaded
	 * configuration, ie thermal strains */
	virtual int HasInternalStrain(void) const;

	/* thermal accessors */
	int ThermalLTfNumber(void) const;
	void SetThermalLTfPtr(const LoadTime* LTfPtr);
	
	double ThermalElongation(void) const; //percentage
		//wrapper functions exist to complete construction of fThermal
		//w/o requiring the header file. For efficiency, actual communication
		//with fThermal is done directly, and therefore requires the header
	 	
	/** \return mass density */
	double Density(void) const;
	
	/** Rayleigh damping. \return mass proportional damping coefficient */
	double MassDamping(void) const;

	/** Rayleigh damping. \return stiffness proportional damping coefficient */
	double StiffnessDamping(void) const;

	/** test for localization. check for bifurcation using current
	 * Cauchy stress and the spatial tangent moduli.
	 * \param normal orientation of the localization if localized
	 * \return 1 if the determinant of the acoustical tensor is negative
	 * or 0 if the determinant is positive. */
	virtual int IsLocalized(dArrayT& normal);

protected:

	/* thermal */
	ThermalDilatationT*	fThermal;

	/* mass density */
	double fDensity;

private:	

	double fMassDamp;
	double fStiffDamp;
	
};

/* returns the density */
inline double SolidMaterialT::Density(void) const { return fDensity; }

/* access to Rayleigh damping parameters */
inline double SolidMaterialT::MassDamping(void) const { return fMassDamp; }
inline double SolidMaterialT::StiffnessDamping(void) const { return fStiffDamp;}

#endif /* _STRUCTURAL_MATERIALT_H_ */
