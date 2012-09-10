
#if !defined(_FSDEMatQ1P0ViscoT_)
#define _FSDEMatQ1P0ViscoT_

#include <cassert>

#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "FSSolidMatT.h"
#include "FSDEMatSupportQ1P0ViscoT.h"
#include "SpectralDecompT.h"	// From RGViscoelasticityT.h
#include "PotentialT.h"	// These from RGSplitT2.h
#include "C1FunctionT.h"

namespace Tahoe {

  class FSDEMatQ1P0ViscoT: public FSSolidMatT
  {

  public:

    // constructors
    FSDEMatQ1P0ViscoT();

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

	// NEW:  non-equilibrium spatial tangent modulus, conversion to material neq modulus
	const dMatrixT& c_ijkl_neq();
	const dMatrixT& C_IJKL_NEQ();

    // Cauchy stress
    virtual const dSymMatrixT& s_ij();

	// NEW:  non-equilibrium spatial stress, conversion to material neq stress
	const dSymMatrixT& s_ij_neq();
	const dSymMatrixT& S_IJ_NEQ();

	// Q1P0 STUFF
	virtual const dMatrixT& b_ij();
	virtual const dArrayT& d_i();
	virtual const dMatrixT& e_ijk();

    // pressure associated with the last computed stress
    double Pressure() const;
	
    // accessors and mutators for material constants
    void SetElectricPermittivity(double epsilon);
    double GetElectricPermittivity() const;

    void SetFSDEMatSupportQ1P0Visco(const FSDEMatSupportQ1P0ViscoT* support);

	/* FUNCTIONS BELOW COPIED FROM RGViscoelasticityT.h */
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
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* Returns eigenvalues of viscous deformation gradient
	Assumes that current values of Cv and Cvn have been loaded using Load(ElementCardT& element, int ip) */
	const dArrayT& Compute_Eigs_v(const int process_id);
	const dArrayT& Compute_Eigs_vn(const int process_id);
	
	void Load(ElementCardT& element, int ip);
	void Store(ElementCardT& element, int ip);

	/* Dimension internal state variables*/
	/* derived class must call RGViscoelaticity::SetStateVariables(fNumProcess)
	  to dimension internal state variable arrays if fNumProcess > 1 (default value) */
	void SetStateVariables (const int numprocess);

	/* FUNCTIONS BELOW COPIED FROM RGSplitT2.h */
	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** compute mechanical strains */
//	virtual const dMatrixT& MechanicalDeformation(void);
	const dMatrixT& MechanicalDeformation(void);
	
	/** compute thermal strains */
//	virtual const dMatrixT& ThermalDeformation_Inverse(void);
	const dMatrixT& ThermalDeformation_Inverse(void);

  protected:

    const FSDEMatSupportQ1P0ViscoT* fFSDEMatSupportQ1P0Visco;

	/* FUNCTIONS BELOW COPIED FROM RGViscoelasticityT.h */
	/* construct symmetric rank-4 mixed-direction tensor (6.1.44) */
  	void MixedRank4_2D(const dArrayT& a, const dArrayT& b, dMatrixT& rank4_ab) const;
  	void MixedRank4_3D(const dArrayT& a, const dArrayT& b, dMatrixT& rank4_ab) const;

  private:

    void Initialize();

    const dMatrixT RightCauchyGreenDeformation();
    const dArrayT ElectricField();
    const dArrayT ElectricField(int ip);

	/* FUNCTIONS BELOW COPIED FROM RGSplitT2.h */
	virtual void Compute_Calg(const dArrayT& tau_dev, const dSymMatrixT& dtau_dev, const double& tau_m, 
						const double& dtau_m, dMatrixT& Calg, const int type);
	virtual void ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
	                   dArrayT& eigenstress, dSymMatrixT& eigenmodulus, const int process_num);  

    // data

  public:

    static const char* Name;

  protected:


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
    dMatrixT fTangentElectromechanicalSpatial;       
    dMatrixT fTangentElectrical;

	/* FUNCTIONS BELOW COPIED FROM RGViscoelasticityT.h */
	/* spectral operations */
	SpectralDecompT fSpectralDecompRef;

	/* FUNCTIONS BELOW COPIED FROM RGViscoelasticityT.h */
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

	/* FUNCTIONS BELOW COPIED FROM RGSplitT2.h */
	/* return values */
	dMatrixT fModulusNEQ;
	dSymMatrixT fStressNEQ;
	
	/* spectral operations */
	SpectralDecompT fSpectralDecompSpat;

	/*mechanical strains*/
	dMatrixT fF_M;
	
	/*thermal strains*/
	dMatrixT fF_T_inv;
	
	/*work space*/
	dSymMatrixT fb;
	dSymMatrixT fbe;
	dSymMatrixT fb_tr;
	dMatrixT fF3D;
	dSymMatrixT fInverse;
	
	dArrayT     fEigs;
	dArrayT     fEigs_e;
	dArrayT     fEigs_tr;
	dArrayT     fEigs_dev;

	dArrayT	    ftau_EQ;
	dArrayT     ftau_NEQ;

	dSymMatrixT fStress3D;
	dSymMatrixT fDtauDe_EQ;
	dSymMatrixT fDtauDe_NEQ;
	dMatrixT fCalg;

	dMatrixT    fModulus3D;
	dMatrixT    fModMat;
  	dMatrixT    fiKAB;
	
	/*potential*/
	ArrayT<PotentialT*> fPot;
  	/*viscosities*/
	ArrayT<C1FunctionT*> fVisc_s;
	ArrayT<C1FunctionT*> fVisc_b;	
	double fietaS;
	double fietaB;


  };

} // namespace Tahoe

//#include "FSDEMatQ1P0ViscoT.i.h"

#endif // _FSDEMatQ1P0ViscoT_
