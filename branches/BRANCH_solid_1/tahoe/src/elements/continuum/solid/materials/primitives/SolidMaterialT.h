/* $Id: SolidMaterialT.h,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
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

class SolidMaterialT: public ContinuumMaterialT
{
public:

	/* constructor */
	SolidMaterialT(ifstreamT& in, const ElasticT& element);

	/* destructor */
	~SolidMaterialT(void);

	/* print parameters */
	virtual void Print(ostream& out) const;

	/* spatial description */
	virtual const dMatrixT& c_ijkl(void) = 0;	// spatial tangent moduli
	virtual const dSymMatrixT& s_ij(void) = 0;	// Cauchy stress

	/* material description */
	virtual const dMatrixT& C_IJKL(void) = 0;	// material tangent moduli
	virtual const dSymMatrixT& S_IJ(void) = 0;	// PK2 stress

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void) = 0;

	/* return the acoustical tensor and wave speeds */
	virtual const dSymMatrixT& AcousticalTensor(const dArrayT& normal) = 0;
	void WaveSpeeds(const dArrayT& normal, dArrayT& speeds);

	/* required parameter flags */
	virtual bool NeedDisp(void) const;
	virtual bool NeedLastDisp(void) const;
	virtual bool NeedVel(void) const;

	/* returns true if the material has internal forces in the unloaded
	 * configuration, ie thermal strains */
	virtual int HasInternalStrain(void) const;

	/* thermal accessors */
	int ThermalLTfNumber(void) const;
	void SetThermalLTfPtr(const LoadTime* LTfPtr);
	
	//double& PlanarDilatationFactor(void);
	//DEV - this is messy
	
	double ThermalElongation(void) const; //percentage
		//wrapper functions exist to complete construction of fThermal
		//w/o requiring the header file. For efficiency, actual communication
		//with fThermal is done directly, and therefore requires the header
	 	
	/* returns the density */
	double Density(void) const;
	
	/* access to Rayleigh damping parameters */
	double MassDamping(void) const;
	double StiffnessDamping(void) const;

	/* returns 1 if the strain localization conditions if satisfied,
	 * .ie if the acoustic tensor has zero (or negative eigenvalues),
	 * for the current conditions (current integration point and strain
	 * state). If localization is detected, the normal (current config)
	 * to the surface is returned in normal */
	virtual int IsLocalized(dArrayT& normal);

protected:

//DEV
//	dMatrixT fModuli;	
	double   fDensity;

	/* thermal */
	ThermalDilatationT*	fThermal;

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
