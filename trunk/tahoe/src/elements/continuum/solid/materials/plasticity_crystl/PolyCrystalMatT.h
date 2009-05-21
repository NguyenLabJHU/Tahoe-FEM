/* $Id: PolyCrystalMatT.h,v 1.14 2009-05-21 22:30:27 tdnguye Exp $ */
#ifndef _POLY_CRYSTAL_MAT_T_H_
#define _POLY_CRYSTAL_MAT_T_H_

#include "FDHookeanMatT.h"

#include "NLCSolverWrapperPtr.h"

#include "GlobalT.h"
#include "ArrayT.h"
#include "dArrayT.h"
#include "Array2DT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "LocalArrayT.h"
#include "ifstreamT.h"

namespace Tahoe {

class SlipGeometry;
class LatticeOrient;
class CrystalElasticity;
class SlipKinetics;
class SlipHardening;
class NLCSolver;

class SolidElementT;

class PolyCrystalMatT : public FDHookeanMatT
{
 public:
  // constructor
  PolyCrystalMatT(void);

  // destructor
  virtual ~PolyCrystalMatT();

	/** material has history variables */
	virtual bool HasHistory(void) const { return true; };

	/** returns true. Derived classes override ContinuumMaterialT::NeedsPointInitialization */
	virtual bool NeedsPointInitialization(void) const;

	/** model initialization. Called per integration point for every
	 * element using the model. Deformation variables are available
	 * during this call. Uses PolyCrystalMatT::NumVariablesPerElement to
	 * allocate the element storage, and then calls PolyCrystalMatT::InitializeVariables
	 * to initialize the state variable space. */
	void PointInitialize(void);

  // required parameter flag
  virtual bool NeedLastDisp() const;

  /* required parameter flags */
  virtual bool Need_F_last(void) const { return true; };

  // some methods to set/initialize member data
  virtual void SetSlipKinetics() = 0;
  virtual void SetSlipHardening() = 0;
  virtual int  NumberOfUnknowns() const = 0;

  // methods invoked from nonlinear constitutive solver
  virtual void FormRHS(const dArrayT& array, dArrayT& rhs) = 0;
  virtual void FormLHS(const dArrayT& array, dMatrixT& lhs) = 0;

  // solve for crystal state / restores/saves solution during subincrementation
  virtual void IterateOnCrystalState(bool& stateConverged, int subIncr) = 0;
  virtual void RestoreSavedSolution() = 0;
  virtual void SaveCurrentSolution() = 0;

  // general accesors
  const double& TimeStep() const;
  ifstreamT& Input_x();
  const int NumGrain() const;
  const int NumSlip() const;
  const dArrayT& MaterialProperties() const;

  // accesors to kinetic equation and slip hardening models
  SlipKinetics& GetSlipKinetics() const;
  SlipHardening& GetSlipHardening() const;

  // some needed accesors (by slip hardening classes mainly)
  const dArrayT& GetResolvedShearStress() const;
  const dArrayT& GetIncrSlipShearStrain() const;

  int Size(void) { return FSMatSupport().Size(); }
  int Rank(void) { return FSMatSupport().Rank(); }

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

 protected:
  /* set (material) tangent modulus */
  virtual void SetModulus(dMatrixT& modulus);

  // subincrementation procedure to compute crystal state
  void SolveCrystalState();

  // function to compute 3D deformations regardless of dimensionality of the
  // problem. For 2D, the out-of-plane direction is x3 and the deformation
  // is assumed to be plane strain
  void Compute_Ftot_3D(dMatrixT& F_3D) const;
  void Compute_Ftot_3D(dMatrixT& F_3D, int ip) const;
  void Compute_Ftot_last_3D(dMatrixT& F_3D) const;
  void Compute_Ftot_last_3D(dMatrixT& F_3D, int ip) const;

  // 4th order tensor transformation: Co_ijkl = F_iI F_jJ F_kK f_lL Ci_IJKL
  void FFFFC_3D(dMatrixT& Co, dMatrixT& Ci, const dMatrixT& F);

 private:

	/** returns the number of variables stored per element */
	virtual int  NumVariablesPerElement(void) = 0;

	/** called by PolyCrystalMatT::PointInitialize */
	virtual void InitializeCrystalVariables(ElementCardT& element) = 0;

  /** return true if material implementation supports imposed thermal
    * strains. This material does not support multiplicative thermal
    * strains. FDHookeanMatT has been updated, but this class needs
    * another look. */
  virtual bool SupportsThermalStrain(void) const { return false; };

  // slip system geometry
  void SetSlipSystems();

  // read lattice orientation data, construct array fEuler
  void SetLatticeOrientation();

  // crystal elasticity
  void SetCrystalElasticity();

  // solver for nonlinear constitutive equations
  void SetConstitutiveSolver();

 protected:

	// time step
	double fdt;

	// references to displacements 
	const LocalArrayT* fLocLastDisp;
	const LocalArrayT* fLocDisp;

	// stream for crystal input data
	ifstreamT fInput;

  // number crystals at each IP
  int fNumGrain;

  // number of slip systems
  int fNumSlip;

  // data to control iterations on state
  int fMaxIterState;
  double fTolerState;

  // control code for supporting classes
  int fCrystalType;
  int fODFCode;
  int fElastCode;
  int fSolverCode;
  int fKinEqnCode;
  int fHardCode;

  // steps to output texture
  int fODFOutInc;

  // iteration counter for local Newton (NLCSolver) and state
  int fIterCount;
  int fIterState;

  // pointers to supporting classes
  SlipGeometry*      fSlipGeometry;
  LatticeOrient*     fLatticeOrient;
  CrystalElasticity* fElasticity;
  SlipKinetics*      fKinetics;
  SlipHardening*     fHardening;
  NLCSolver*         fSolver;

  // handle to NLCSolver
  NLCSolverWrapperPtr fSolverPtr;

  // material properties
  dArrayT fMatProp;

  // total deformation gradients
  dMatrixT fFtot_n;
  dMatrixT fFtot;
  dMatrixT fFt;

  // incremental slip shearing rate / resolve shear stres
  dArrayT fDGamma;
  dArrayT fGamma; //integrated plastic slip
  dArrayT fTau;

  // Schmidt tensor in crystal/sample coords
  ArrayT<dMatrixT> fZc;
  ArrayT<dMatrixT> fZ;

  // slip tensors in crystal/sample coords
  ArrayT<dArrayT> fSlipSc;   // slip direction
  ArrayT<dArrayT> fSlipMc;   // slip plane normal
  ArrayT<dArrayT> fSlipS;
  ArrayT<dArrayT> fSlipM;

  // array for Euler angles at integration point
  ArrayT<dArrayT> fangles;

  // huge temporary array to hold all euler angles
  ArrayT<Array2DT<dArrayT> > fEuler; 

  // rotation matrix from Euler angles
  dMatrixT fRotMat;

  // crystal and aggregate Cauchy stress
  dSymMatrixT fs_ij;
  dSymMatrixT fsavg_ij;

  // crystal and aggregate Moduli
  dMatrixT fc_ijkl;
  dMatrixT fcavg_ijkl;
};

// general needed accesors 
inline const double& PolyCrystalMatT::TimeStep() const { return fdt; }

inline ifstreamT& PolyCrystalMatT::Input_x() { return fInput; }
inline const int PolyCrystalMatT::NumGrain() const { return fNumGrain; }
inline const int PolyCrystalMatT::NumSlip() const { return fNumSlip; }
inline const dArrayT& PolyCrystalMatT::MaterialProperties() const { return fMatProp; }

inline SlipKinetics& PolyCrystalMatT::GetSlipKinetics() const { return *fKinetics; }
inline SlipHardening& PolyCrystalMatT::GetSlipHardening() const { return *fHardening; }

inline const dArrayT& PolyCrystalMatT::GetResolvedShearStress() const { return fTau; }
inline const dArrayT& PolyCrystalMatT::GetIncrSlipShearStrain() const { return fDGamma; }

} // namespace Tahoe 
#endif /* _POLY_CRYSTAL_MAT_T_H_ */
