/* $Id: SolidMaterialT.h,v 1.22 2005-02-10 17:59:56 raregue Exp $ */
/* created: paklein (11/20/1996) */
#ifndef _STRUCTURAL_MATERIALT_H_
#define _STRUCTURAL_MATERIALT_H_

#include "GlobalT.h"

/* base class */
#include "ContinuumMaterialT.h"

/* direct members */
#include "dMatrixT.h"
#include "ThermalDilatationT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class ElementBaseT;
class dMatrixT;
class ThermalDilatationT;
class ScheduleT;
class dSymMatrixT;
class LocalArrayT;
class SolidElementT;

/** base class for constitutive models for solids */
class SolidMaterialT: public ContinuumMaterialT
{
public:

	/** \name 2D constrain options */
	enum ConstraintT {
		kNoConstraint = 0, /**< no constraint, material is 3D */
		kPlaneStress = 1, /**< plane stress */
		kPlaneStrain = 2  /**< plane strain */};

	/** constructor */
	SolidMaterialT(void);

	/** destructor */
	~SolidMaterialT(void);

	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void) = 0;
	
	/** spatial elastic modulus */
	virtual const dMatrixT& ce_ijkl(void) = 0;

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void) = 0;

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. The value is not guaranteed to
	 * persist during intervening calls to any other non-const
	 * accessor. \return 1/3 of the trace of the three-dimensional
	 * stress tensor, regardless of the dimensionality of the
	 * problem. */
	virtual double Pressure(void) const = 0;
	/*@}*/

	/** \name material description */
	/*@{*/
	/** material tangent moduli */
	virtual const dMatrixT& C_IJKL(void) = 0;

	/** 2nd Piola-Kirchhoff stress */
	virtual const dSymMatrixT& S_IJ(void) = 0;
	/*@}*/

	/** 2D constrain options or kNoConstraint::kNoConstraint if the material
	 * is not 2D */
	ConstraintT Constraint(void) const { return fConstraint; };

	/** \name queries */
	/*@{*/
	/** return true if the material can produce localization */
	virtual bool HasLocalization(void) const { return false; };

	/** return true if the material generates heat. The returns false unless 
	 * overridden. */
	virtual bool HasIncrementalHeat(void) const { return false; };

	virtual bool NeedDisp(void) const     { return false; };
	virtual bool NeedLastDisp(void) const { return false; };
	virtual bool NeedVel(void) const      { return false; };
	
	/** return true if the density varies with position */
	virtual bool HasChangingDensity(void) const { return false; };
	/*@}*/

	/** incremental heat generation (energy/volume). The value should be the amount of
	 * heat associated with the updated stress calculated with the most recent call to 
	 * SolidMaterialT::s_ij or SolidMaterialT::S_IJ */
	virtual double IncrementalHeat(void);

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
	
	/** return the strain in the material. The definition of strain will be
	 * dependent on the subclass */
	virtual void Strain(dSymMatrixT& strain) = 0;

	/** returns true if the material has internal forces in the unloaded
	 * configuration, i.e. thermal strains */
	int HasThermalStrain(void) const;

	/** returns the schedule number for the imposed thermal strain */
	int ThermalStrainSchedule(void) const;

	/** set the schedule for the prescribed temperature */
	void SetThermalSchedule(const ScheduleT* LTfPtr);
	
	/** return the thermal expansion rate as a percentage */
	double ThermalElongation(void) const;
	 	
	/** \return mass density */
	virtual double Density(void);

	/** test for localization. check for bifurcation using current
	 * Cauchy stress and the spatial tangent moduli.
	 * \param normals orientation of the localization if localized
	 * \return true if the determinant of the acoustical tensor is negative
	 * or false if the determinant is positive. */
	//virtual int IsLocalized(dArrayT& normal);
	virtual bool IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list. This function also checks if thermal strain are 
	 * being imposed and if the material supports thermal strain, using 
	 * SolidMaterialT::SupportsThermalStrain. */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** \name internal variables */
	/*@{*/
	virtual const iArrayT& InternalDOF(void) const;
	virtual const dArrayT& InternalStressVars(void);
	virtual const dArrayT& InternalStrainVars(void);
	/*@}*/
	
	/** return true if material implementation supports imposed thermal
	 * strains. */
	virtual bool SupportsThermalStrain(void) const { return false; };

protected:

	/* thermal */
	ThermalDilatationT*	fThermal;

	/* mass density */
	double fDensity;

	/** 2D constrain option */
	ConstraintT fConstraint;

	/** coefficient of thermal expansion in strain/unit temperature */
	double fCTE;
};

/* incremental heat generation */
inline double SolidMaterialT::IncrementalHeat(void) { return 0.0; }

/* returns the density */
inline double SolidMaterialT::Density(void) { return fDensity; }

/* imposed thermal strains */
inline int SolidMaterialT::HasThermalStrain(void) const { return fThermal->IsActive(); }
inline int SolidMaterialT::ThermalStrainSchedule(void) const { return fThermal->ScheduleNum(); }
inline void SolidMaterialT::SetThermalSchedule(const ScheduleT* LTfPtr) { fThermal->SetSchedule(LTfPtr); }
inline double SolidMaterialT::ThermalElongation(void) const { return fThermal->PercentElongation(); }


} // namespace Tahoe 
#endif /* _STRUCTURAL_MATERIALT_H_ */
