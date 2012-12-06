
#if !defined(_FSDEMatQ1P0SurfaceT_)
#define _FSDEMatQ1P0SurfaceT_

#include <cassert>

#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "FSSolidMatT.h"
#include "FSDEMatSupportQ1P0SurfaceT.h"

namespace Tahoe {

/* From FCC3D_Surf */
class FCCLatticeT_Q1P0Surf;
class PairPropertyT;
class BondLatticeT;

  class FSDEMatQ1P0SurfaceT: public FSSolidMatT
  {

  public:

    // constructors
    FSDEMatQ1P0SurfaceT();

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

	// Q1P0Surface STUFF
	virtual const dMatrixT& b_ij();
	virtual const dArrayT& d_i();
	virtual const dMatrixT& e_ijk();

    // pressure associated with the last computed stress
    double Pressure() const;
	
    // accessors and mutators for material constants
    void SetElectricPermittivity(double epsilon);
    double GetElectricPermittivity() const;

    void SetFSDEMatSupportQ1P0Surface(const FSDEMatSupportQ1P0SurfaceT* support);

	/* FCC3D_Surf */
	/** destructor */
	~FSDEMatQ1P0SurfaceT(void);

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** \name Cauchy-Born parameters */
	/*@{*/
	/** return a reference to the bond lattice */
	const BondLatticeT& BondLattice(void) const;

	/** reference volume */
	double CellVolume(void) const { return fAtomicVolume; };

	/** nearest neighbor distance */
	double NearestNeighbor(void) const { return fNearestNeighbor; };
	/*@}*/

	/** thickness of surface layer to subtract off of bulk */
	double SurfaceThickness(void) const { return fSurfaceThickness; };

	/** compute the symetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli);	
	
	/** symmetric 2nd Piola-Kirchhoff stress */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2);

  protected:

    const FSDEMatSupportQ1P0SurfaceT* fFSDEMatSupportQ1P0Surface;

	/* FCC3D_Surf */
	/** return the equi-axed stretch at which the stress is zero. This method
	 * assumes the material is isotropic when subject to equi-axed stretch. */
	double ZeroStressStretch(void);	

  private:

    void Initialize();

    const dMatrixT RightCauchyGreenDeformation();
    const dArrayT ElectricField();
    const dArrayT ElectricField(int ip);

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

	/* FCC3D_Surf */
	/** nearest neighbor distance */
	double fNearestNeighbor;

	/** surface layer thickness */
	double fSurfaceThickness;

	/** bond information */
	FCCLatticeT_Q1P0Surf* fFCCLattice_Q1P0Surf;

	/** pair interaction potential */
	PairPropertyT* fPairProperty;

	/** \name work space */
	/*@{*/
	dMatrixT fBondTensor4;
	dArrayT  fBondTensor2;
	/*@}*/

	/** atomic volume */
	double fAtomicVolume;

	/** atomic area for surface cauchy-born */
	double fAtomicArea;

	/** dummy full bond density array */
	/* THIS IS THE REFERENCE VOLUME/AREA */
	dArrayT fFullDensity;
	
	/** flag to indicate whether stress calculation for output should include
	 * the full bond density */
	bool fFullDensityForStressOutput;

  };

} // namespace Tahoe

#endif // _FSDEMatQ1P0SurfaceT_
