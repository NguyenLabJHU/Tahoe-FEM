/* $Id: tevp2D.h,v 1.3 2001-05-01 19:03:10 hspark Exp $ */
/* Thermoelasto-viscoplastic material used to generate shear bands */
/* Created:  Harold Park (04/04/2001) */
/* Last Updated:  Harold Park (05/01/2001) */

#ifndef _TEVP_2D_H_
#define _TEVP_2D_H_

/* base classes */
#include "FDStructMatT.h"
#include "IsotropicT.h"
#include "iArrayT.h"
#include "Material2DT.h"

/* forward declarations */
class ShapeFunctionT;
class ElementCardT;

class tevp2D: public FDStructMatT, public IsotropicT, public Material2DT
{
 public:
  /* constructor */
  tevp2D(ifstreamT& in, const ElasticT& element);
  
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

  /* spatial description - this IS implemented */
  virtual const dMatrixT& c_ijkl(void);  // spatial tangent moduli
  virtual const dSymMatrixT& s_ij(void);  // Cauchy stress

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
  
  /* accessor functions to be inlined - these should return the value from
   * the previous timestep */
  double GetTimeStep(void) const;   // Returns the current timestep
  double GetTime(void) const;       // Returns the total time elapsed
  double GetEbtot(void) const;      
  double GetCtcon(void) const;
  double GetXxii(void) const;
  dArrayT GetPP(void) const;
  dMatrixT GetDmat(void) const;     // Return the elastic modulus tensor
  dArrayT GetSmlp(void) const;

 private:
  /* computational functions */

  /* deformation gradient, rate of deformation, spin */
  void ComputeGradients(void);
  dMatrixT& CauchyToKirchoff(dMatrixT temp_stress);
  double ComputeEbtot(void);     
  // Computes the incremental effective strain
  double ComputeXxii(void);
  double ComputeCtcon(void);
  dArrayT& ComputePP(void); 
  double ComputeEcc(void);
  dMatrixT& ComputeDmat(void);   // Original elastic coefficient tensor
  dArrayT& ComputeEP_tan(void);  // Modulus correction
  dArrayT& ComputeSmlp(void);
  enum LoadingStatusT {kIsPlastic = 0, kIsElastic = 1, kReset = 3};
  // Should LoadingStatusT be protected?
  void AllocateElement(ElementCardT& element); // If element/IP goes plastic

  /* Enumerated data types/definitions */
  enum InternalVariablesT {kTemp = 0,   // Temperature
			     kSb = 1,   // Effective Stress
                             kEb = 2};  // Effective Strain
  enum ModelT {kTevp = 0,          // Thermo-elasto-viscoplastic
               kFluid = 1};        // Fluid model
  enum StessComponentsT {kSig11 = 0,
			 kSig12 = 1,
			 kSig22 = 2,   // fTempStress stored like this...
			 kSig33 = 3};

  /* Output values/internal variable functions below - these functions
   * should ONLY be called AFTER the stress and modulus have been computed */
  double ComputeTemperature(ElementCardT& element, int ip);
  double ComputeEffectiveStress(void);
  double ComputeEffectiveStrain(ElementCardT& element, int ip);
  int CheckCriticalStrain(ElementCardT& element, int ip);
  int CheckIfPlastic(ElementCardT& element, int ip);
  /* load element data for the specified integration point */
  void LoadData(const ElementCardT& element, int ip); 
  /* element level data */
  void Update(ElementCardT& element);
  void Reset(ElementCardT& element);
  /* These take the internal variable array and expand it to a matrix, or
   * vice versa */
  dMatrixT ArrayToMatrix(dArrayT StressArray);
  dArrayT MatrixToArray(dSymMatrixT StressMatrix);  
  dSymMatrixT Return3DStress(dMatrixT StressMatrix);

 protected:
  /* return values */
  dSymMatrixT fStress;
  dMatrixT fModulus;
  double fStrainEnergyDensity;   // How do I define this for this material?
  
/* element level internal variables */
  dArrayT fInternal;             // Internal variables
  dArrayT fTempStress;      // Store the Kirchoff stress from the previous
                            // timestep (Sig11, Sig12=Sig21, Sig22, Sig33)

 private:

  const double& fDt;           // Timestep
  const double& fTime;         // Total time elapsed

  /* work space */
  dMatrixT fFtot_2D;           // Deformation gradient 2D
  dMatrixT fFtot;              // Total deformation gradient (3D)
  dMatrixT fDtot;              // Rate of deformation (3D)
  dMatrixT fFinv_2D;           // Inverse of F needed to compute D
  dMatrixT fGradV_2D;          // Velocity gradient (2D)
  dMatrixT fGradV;             // Velocity gradient (3D)
  dMatrixT fKirchoff;          // Kirchoff stress tensor
  const LocalArrayT& fLocVel;  // Nodal velocities
  const LocalArrayT& fLocDisp; // Nodal displacements (not necessary)
  dMatrixT fF_temp;            // Deformation gradient to work with
  double fSpin;                // Spin scalar - only Spin(1,2) non-zero...
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
  dSymMatrixT fStress3D;       // 3D version of stress -> fStress
  dArrayT fSmlp;

  /* output variables/internal variables */
  double fTemperature;         // Temperature
  double fSb;                  // Effective stress
  double fEb;                  // Effective strain
  
  /* Global material constants - class scope variables */
  double Temp_0, El_E, El_V, El_K, El_G, Sb0, Rho0, Eb0, Eb0tot, BigN, Smm;
  double Alpha_T, Delta, Theta, Kappa, Cp, Chi, Ccc, Pcp;
  double Epsilon_1, Epsilon_2, Epsilon_rate, Gamma_d, Mu_d, SigCr;
  double Xi;

  /* shape functions */
  const ShapeFunctionT& fShapes;   // Needed to compute velocity gradient

};

/* inline functions / accessor functions */
inline double tevp2D::GetTimeStep(void) const { return fDt; }
inline double tevp2D::GetTime(void) const { return fTime; }
inline double tevp2D::GetEbtot(void) const { return fEbtot; }
inline double tevp2D::GetCtcon(void) const { return fCtcon; }
inline double tevp2D::GetXxii(void) const { return fXxii; }
inline dArrayT tevp2D::GetPP(void) const { return fPP; }
inline dMatrixT tevp2D::GetDmat(void) const { return fDmat; }
inline dArrayT tevp2D::GetSmlp(void) const { return fSmlp; }

#endif /* _TEVP_2D_H_ */
				





