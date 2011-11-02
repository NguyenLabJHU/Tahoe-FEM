
#if !defined(_FSDEMatViscoT_)
#define _FSDEMatViscoT_

#include <cassert>

#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "FSSolidMatT.h"
#include "FSDEMatSupportViscoT.h"
#include "SpectralDecompT.h"

namespace Tahoe {

  class FSDEMatViscoT: public FSSolidMatT
  {

  public:

    // constructors
    FSDEMatViscoT();

    // set parameters
    void DefineParameters(ParameterListT& list) const;
    void TakeParameterList(const ParameterListT& list);

    // information about subordinate parameter lists
    virtual void DefineSubs(SubListT& sub_list) const;

    // Interface required by Tahoe
    double StrainEnergyDensity();

    // material mechanical tangent modulus
    virtual const dMatrixT& C_IJKL();

    // material electromechanical tangent modulus
    virtual const dMatrixT& E_IJK();

    // material electric tangent modulus
    virtual const dMatrixT& B_IJ();

    // Second Piola-Kirchhoff stress
    virtual const dSymMatrixT& S_IJ();

    // electric displacement
    virtual const dArrayT& D_I();

    // electric field
    virtual const dArrayT& E_I();

    // spatial mechanical tangent modulus
    virtual const dMatrixT& c_ijkl();

    // Cauchy stress
    virtual const dSymMatrixT& s_ij();

    // pressure associated with the last computed stress
    double Pressure() const;
	
    // accessors and mutators for material constants
    void SetElectricPermittivity(double epsilon);
    double GetElectricPermittivity() const;

    void SetFSDEMatSupportVisco(const FSDEMatSupportViscoT* support);

	// FUNCTIONS BELOW COPIED FROM RGViscoelasticityT.h
	/** return true if the material has history variables */
	virtual bool HasHistory(void) const { return true; };
	
	/*Initialize history variable*/
	virtual bool NeedsPointInitialization(void) const {return true;}; 
	virtual void PointInitialize(void);              

	/* update/reset internal variables */
	virtual void UpdateHistory(void); // element at a time
	virtual void ResetHistory(void);  // element at a time
	/* apply pre-conditions at the current time step */
	virtual void InitStep(void){ FSSolidMatT::InitStep(); };
	
	/* form of tangent matrix (symmetric by default) */
//	virtual GlobalT::SystemTypeT TangentTypeDE(void) const;

	/*Returns eigenvalues of viscous deformation gradient
	Assumes that current values of Cv and Cvn have been loaded using Load(ElementCardT& element, int ip)*/
	const dArrayT& Compute_Eigs_vDE(const int process_id);
	const dArrayT& Compute_Eigs_vnDE(const int process_id);
	
	void LoadDE(ElementCardT& element, int ip);
	void StoreDE(ElementCardT& element, int ip);

	/* Dimension internal state variables*/
	/*derived class must call RGViscoelaticity::SetStateVariables(fNumProcess)
	  to dimension internal state variable arrays if fNumProcess > 1 (default value)*/
	void SetStateVariablesDE (const int numprocess);

  protected:

    const FSDEMatSupportViscoT* fFSDEMatSupportVisco;

	// FUNCTIONS BELOW COPIED FROM RGViscoelasticityT.h
	/* construct symmetric rank-4 mixed-direction tensor (6.1.44) */
  	void MixedRank4_2DDE(const dArrayT& a, const dArrayT& b, dMatrixT& rank4_ab) const;
  	void MixedRank4_3DDE(const dArrayT& a, const dArrayT& b, dMatrixT& rank4_ab) const;

  private:

    void Initialize();

    const dMatrixT RightCauchyGreenDeformation();
    const dArrayT ElectricField();
    const dArrayT ElectricField(int ip);

    // data

  public:

    static const char* Name;

  protected:

	// FUNCTIONS BELOW COPIED FROM RGViscoelasticityT.h
	/*internal state variables. Dimension numprocess<nsd x nsd>*/
	ArrayT<dSymMatrixT> fC_v;
	ArrayT<dSymMatrixT> fC_vn;
	
	/* number of nonequilibrium processes*/
	/* must be set in derived classes before TakeParameterList is called*/
	/* default value is 1*/
	int fNumProcess;
	
	/*number of state variables*/
	int fnstatev;
	
	/* internal state variables array*/
	dArrayT fstatev;

  private:

    double fElectricPermittivity;
    double fEnergyDensity;
    double fMu;
    double fNrig;
    double fLambda;
	double fKappa;

    dArrayT fElectricField;
    dArrayT fElectricDisplacement;
	dArrayT fParams;
	
    dSymMatrixT fStress;
    dMatrixT fTangentMechanicalElec;
    dMatrixT fTangentMechanical;
    dMatrixT fTangentElectromechanical;
    dMatrixT fTangentElectrical;

	// FUNCTIONS BELOW COPIED FROM RGViscoelasticityT.h
	/* spectral operations */
	SpectralDecompT fSpectralDecompRef;

  };

} // namespace Tahoe

#include "FSDEMatViscoT.i.h"

#endif // _FSDEMatViscoT_
