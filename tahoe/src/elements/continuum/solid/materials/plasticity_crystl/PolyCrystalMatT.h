/*
  File: PolyCrystalMatT.h
*/

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
#include "ifstreamT.h" 
#include "LocalArrayT.h"

class SlipGeometry;
class LatticeOrient;
class CrystalElasticity;
class SlipKinetics;
class SlipHardening;
class NLCSolver;

class ElasticT;

class PolyCrystalMatT : public FDHookeanMatT
{
 public:
  // constructor
  PolyCrystalMatT(ifstreamT& in, const FiniteStrainT& element);

  // destructor
  virtual ~PolyCrystalMatT();

  // allocate space/initialize crystal arrays (all)
  virtual bool NeedsInitialization() const;
  virtual void Initialize();

  // required parameter flag
  virtual bool NeedLastDisp() const;

  	/** required parameter flags */
	virtual bool Need_F_last(void) const { return true; };

  // some methods to set/initialize member data
  virtual void SetSlipKinetics() = 0;
  virtual void SetSlipHardening() = 0;
  virtual void InitializeCrystalVariables() = 0;
  virtual int  NumVariablesPerElement() = 0;

  // methods invoked from nonlinear constitutive solver
  virtual void FormRHS(const dArrayT& dgamma, dArrayT& rhs) = 0;
  virtual void FormLHS(const dArrayT& dgamma, dMatrixT& lhs) = 0;

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
  virtual const dArrayT& GetResolvedShearStress() const = 0;
  virtual const dArrayT& GetIncrSlipShearStrain() const = 0;

  // print data read
  virtual void Print(ostream& out) const;

 protected:

	/* set (material) tangent modulus */
	virtual void SetModulus(dMatrixT& modulus);

  // print name
  virtual void PrintName(ostream& out) const;

  // allocate all elements at once
  void AllocateElements();

 private:

 	/** return true if material implementation supports imposed thermal
	 * strains. This material does not support multiplicative thermal
	 * strains. FDHookeanMatT has been updated, but this class needs
	 * another look. */
	virtual bool SupportsThermalStrain(void) const { return false; };

  // slip system geometry
  void SetSlipSystems();

  // read lattice orientation data, construct array fEuler
  void SetLatticeOrientation() ;

  // crystal elasticity
  void SetCrystalElasticity();

  // solver for nonlinear constitutive equations
  void SetConstitutiveSolver();

 protected:
  // current time & time step
  const double& ftime;
  double fdt;

  // status flag
  const GlobalT::StateT& fStatus;

  // references to displacements 
  const LocalArrayT& fLocLastDisp;
  const LocalArrayT& fLocDisp;

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

  // iteration counter for NLCSolver (DGamma) and state
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

  // Schmidt tensor in crystal coords
  ArrayT<dMatrixT> fZc;

  // array for Euler angles at integration point
  ArrayT<dArrayT> fangles;

  // huge temporary array to hold all euler angles
  ArrayT<Array2DT<dArrayT> > fEuler; 

  // aggregate Cauchy stress
  dSymMatrixT fsavg_ij;

  // aggregate Moduli
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

#endif /* _POLY_CRYSTAL_MAT_T_H_ */
