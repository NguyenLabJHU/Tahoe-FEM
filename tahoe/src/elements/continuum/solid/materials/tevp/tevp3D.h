/* $Id: tevp3D.h,v 1.10 2003-11-21 22:46:58 paklein Exp $ */
/* created:  Harold Park (06/25/2001) */
#ifndef _TEVP_3D_H_
#define _TEVP_3D_H_

/* base classes */
#include "FSSolidMatT.h"
#include "IsotropicT.h"
#include "iArrayT.h"

namespace Tahoe {

/* forward declarations */
class ElementCardT;

/** Thermoelasto-viscoplastic material used to generate shear bands */
class tevp3D: public FSSolidMatT, public IsotropicT
{
 public:
  /* constructor */
  tevp3D(ifstreamT& in, const FSMatSupportT& support);
  
  /* materials initialization */
  virtual bool NeedsPointInitialization(void) const { return true; }
  virtual void PointInitialize(void);

  /* required parameter flags */
  virtual bool NeedVel(void) const;

  /* update internal variables */
  virtual void UpdateHistory(void);

  /* reset internal variables to last converged solution */
  virtual void ResetHistory(void);

  /* print parameters */
  virtual void Print(ostream& out) const;
  virtual void PrintName(ostream& out) const;

	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. See SolidMaterialT::Pressure
	 * for more information. */
	virtual double Pressure(void) const { return fStress.Trace()/3.0; };
	/*@}*/

  /* material description - not implemented */
  virtual const dMatrixT& C_IJKL(void);  // material tangent moduli
  virtual const dSymMatrixT& S_IJ(void);  // PK2 stress

  /* strain energy density */
  virtual double StrainEnergyDensity(void);

  /* returns the number of variables computed for nodal extrapolation
   * during for element output, ie internal variables */
  virtual int NumOutputVariables(void) const;
  virtual void OutputLabels(ArrayT<StringT>& labels) const;
  virtual void ComputeOutput(dArrayT& output);

 private:
  /* computational functions */

  /* deformation gradient, rate of deformation, spin */
  void ComputeF(void);
  void ComputeD(void);
  dMatrixT& ComputeSpin(void);
  void ComputeEbtotCtconXxii(void);
  void ComputePP(void); 
  double ComputeEcc(void);
  void ComputeDmat(void);   // Original elastic coefficient tensor
  dArrayT& ComputeEP_tan(void);  // Modulus correction
  void ComputeSmlp(void);
  enum LoadingStatusT {kIsPlastic = 0, kIsElastic = 1, kReset = 3};

  void AllocateElement(ElementCardT& element); // If element/IP goes plastic

  /* Enumerated data types/definitions */
  enum InternalVariablesT {kTemp = 0,   // Temperature
                             kSb = 1,   // Effective Stress
                             kEb = 2};  // Effective Strain
  enum ModelT {kTevp = 0,          // Thermo-elasto-viscoplastic
               kFluid = 1,
               kCrack = 2};        // Fluid model

  /* Output values/internal variable functions below - these functions
   * should ONLY be called AFTER the stress and modulus have been computed */
  double ComputeFluidTemperature(void);
  double ComputeViscoTemperature(void);
  double ComputeEffectiveStress(void);
  double ComputeFluidEffectiveStrain(void);
  double ComputeViscoEffectiveStrain(void);
  void CheckCriticalCriteria(ElementCardT& element, int ip);
  int CheckIfPlastic(ElementCardT& element, int ip);
  /* load element data for the specified integration point */
  void LoadData(const ElementCardT& element, int ip); 
  /* element level data */
  void Update(ElementCardT& element);
  void Reset(ElementCardT& element);
  /* These take the internal variable array and expand it to a matrix, or
   * vice versa */
  dMatrixT& ArrayToMatrix(const dArrayT& StressArray);
  dArrayT& MatrixToArray(const dSymMatrixT& StressMatrix);  
  dSymMatrixT& Return3DStress(const dMatrixT& StressMatrix);
  dSymMatrixT& ArrayToSymMatrix3D(const dArrayT& StressArray);

 protected:
  /* return values */
  dSymMatrixT fStress;
  dMatrixT fModulus;
  double fStrainEnergyDensity;   // How do I define this for this material?

  /* element level internal variables */
  dArrayT fInternal;             // Internal variables
  dArrayT fTempKirchoff;      // Store the Kirchoff stress from the previous
                            // timestep (S11, S22, S33, S32, S31, S21)
  dArrayT fTempCauchy;      // Store the Cauchy stress from the previous
                            // timestep (S11, S22, S33, S32, S31, S21)
 private:

	double fDt;           // Timestep

  /* work space */
  dMatrixT fFtot;              // Total deformation gradient (3D)
  dSymMatrixT fDtot;         // Symmetric version of rate of deformation
  dMatrixT fGradV;             // Velocity gradient (3D)
  const LocalArrayT* fLocVel;  // Nodal velocities
//  const LocalArrayT& fLocDisp; // Nodal displacements (not necessary)
  dMatrixT fF_temp;            // Deformation gradient to work with
  dMatrixT fSpin;              // Spin tensor
  int fCriticalStrain;         // Checks if critical strain criteria is met
  double fEbtot;               // Incremental effective strain
  double fXxii;                // Used in computing tangent modulus
  double fCtcon;               // Used in computing tangent modulus
  dArrayT fPP;
  dMatrixT fDmat;              // Original elastic coefficient tensor
  dArrayT fEP_tan;             // Plastic correction to Jaumann stress rate
  double fEcc;                 
  dMatrixT fStressMatrix;      // Expand stress array to a matrix
  dArrayT fStressArray;        // Flatten stress matrix to an array
  dSymMatrixT fStill3D;        // Still 3D version of stress
  dSymMatrixT fStress3D;
  dArrayT fSmlp;
  dSymMatrixT fSymStress3D; 
  double fJ;                   // Jacobian of deformation gradient
  double fVisc;                // Original viscosity

  /* output variables/internal variables */
  double fTemperature;         // Temperature
  double fSb;                  // Effective stress
  double fEb;                  // Effective strain
  
  /* Global material constants - class scope variables */
  double Temp_0, El_E, El_V, El_K, El_G, Sb0, Rho0, Eb0, Eb0tot, BigN, Smm;
  double Alpha_T, Delta, Theta, Kappa, Cp, Chi, Ccc, Pcp;
  double Epsilon_1, Epsilon_2, Epsilon_rate, Gamma_d, Mu_d, SigCr;
  double Xi;
};

} // namespace Tahoe 
#endif /* _TEVP_3D_H_ */
                                



