/* $Id: EVPFDBaseT.h,v 1.8 2003-01-29 07:35:06 paklein Exp $ */
#ifndef _EVP_FD_BASE_T_H_
#define _EVP_FD_BASE_T_H_

#include "IsotropicT.h"
#include "FDHookeanMatT.h"

#include "NLCSolverWrapperPtr.h"
#include "KineticEqnBase.h"

#include "GlobalT.h"
#include "ArrayT.h"
#include "dArrayT.h"
#include "Array2DT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "ifstreamT.h" 

namespace Tahoe {

class NLCSolver;
class SolidElementT;

class EVPFDBaseT : public FDHookeanMatT, public IsotropicT
{
 public:
  // constructor
  EVPFDBaseT(ifstreamT& in, const FSMatSupportT& support);

  // destructor
  virtual ~EVPFDBaseT();

  // allocate space/initialize crystal arrays (all)
  virtual void Initialize();

	/** returns true. Derived classes override ContinuumMaterialT::NeedsPointInitialization */
	virtual bool NeedsPointInitialization() const;

	/** model initialization. Called per integration point for every
	 * element using the model. Deformation variables are available
	 * during this call. Uses EVPFDBaseT::NumVariablesPerElement to
	 * allocate the element storage, and then calls EVPFDBaseT::InitializeVariables
	 * to initialize the state variable space. */
	virtual void PointInitialize(void);

	/** apply initialize current time step */
	virtual void InitStep(void);

  // required parameter flag
  virtual bool NeedLastDisp() const;
  
  	/** required parameter flags */
	virtual bool Need_F_last(void) const { return true; };

  // PVFs defined in derived classes 
  virtual void SetKineticEquation() = 0;

  // PVFs invoked from NLCSolver and defined in derived classes
  virtual void FormRHS(const dArrayT& variab, dArrayT& rhs) = 0;
  virtual void FormLHS(const dArrayT& variab, dMatrixT& lhs) = 0;

  // general accessors
  ifstreamT& Input();
  double Temperature();

  // print data read
  virtual void Print(ostream& out) const;

 protected:

	/* set (material) tangent modulus */
	virtual void SetModulus(dMatrixT& modulus);

  // print name
  virtual void PrintName(ostream& out) const;

  // size of system of nonlinear equations
  virtual int GetNumberOfEqns();

	// function to compute 3D deformation regardless of dimensionality of the
	// problem. For 2D, the out-of-plane direction is x3 and the deformation
	// is assumed to be plane strain
	void Compute_Ftot_3D(dMatrixT& F_3D) const;	
	void Compute_Ftot_last_3D(dMatrixT& F_3D) const;	

 private:

	/** return the number variables stored per element */
	virtual int NumVariablesPerElement() = 0;

	/** initialize the state variables for the given element */
	virtual void InitializeVariables(ElementCardT& element) = 0;

 	/** return true if material implementation supports imposed thermal
	 * strains. This material does not support multiplicative thermal
	 * strains. FDHookeanMatT has been updated, but this class needs
	 * another look. */
	virtual bool SupportsThermalStrain(void) const { return false; };

  // solver for nonlinear constitutive equations
  void SetConstitutiveSolver();

 protected:

	// current time & time step
	// const double& ftime;
	double fdt;

	// temperature
	double fTheta;

  // pointers to displacements
  const LocalArrayT* fLocDisp;
  const LocalArrayT* fLocLastDisp;

  // stream for input data
  ifstreamT fInput;

  // Lame's elastic constants
  double fmu;
  double flambda;
  double fbulk;

  // control code for supporting classes
  int fSolverCode;
  int fKinEqnCode;

  // iteration counter of NLCSolver
  int fIterCount;

  // pointers to supporting classes
  KineticEqnBase* fKineticEqn;
  NLCSolver*  fSolver;

  // handle to NLCSolver
  NLCSolverWrapperPtr fSolverPtr;

  // 3D total deformation gradient
  dMatrixT fFtot;

  // Cauchy stress
  dSymMatrixT fs_ij;

  // Tangent Moduli
  dMatrixT fc_ijkl;
};

inline ifstreamT& EVPFDBaseT::Input() { return fInput; }
inline double EVPFDBaseT::Temperature() { return fTheta; }

} // namespace Tahoe 
#endif /* _EVP_FD_BASE_T_H_ */
