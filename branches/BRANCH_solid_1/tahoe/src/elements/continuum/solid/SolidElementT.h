/* $Id: SolidElementT.h,v 1.4.2.2 2001-06-29 01:21:14 paklein Exp $ */
/* created: paklein (05/28/1996)                                          */

#ifndef _ELASTIC_T_H_
#define _ELASTIC_T_H_

/* base class */
#include "ContinuumElementT.h"

/* direct members */
#include "dArray2DT.h"
#include "LocalArrayT.h"
#include "dSymMatrixT.h"

/* forward declarations */
class ShapeFunctionT;
class StructuralMaterialT;
class StringT;

/** base class of elements for deformation of solids */
class SolidElementT: public ContinuumElementT
{
public:
	
	/** list/index of nodal outputs */
	enum NodalOutputCodeT {
		iNodalCoord = 0, /**< (reference) coordinates */
 	     iNodalDisp = 1, /**< displacements */
       iNodalStress = 2, /**< extrapolated stresses */
         iPrincipal = 3, /**< extrapolated principal stresses */
     iEnergyDensity = 4, /**< extrapolated strain energy density */
        iWaveSpeeds = 5, /**< extrapolated local wave speeds */
      iMaterialData = 6};/**< extrapolated  model output */

	/** list/index of element outputs */
	enum ElementOutputCodeT {
	      iCentroid = 0, /**< (reference) centroid coordinates */
		      iMass = 1, /**< integrated element mass */
	  iStrainEnergy = 2, /**< integrated strain energy */
	 iKineticEnergy = 3, /**< integrated kinetic energy */
    iLinearMomentum = 4, /**< integrated linear momentum */
          iIPStress = 5, /**< integration point stresses */
    iIPMaterialData = 6};/**< integration point material model output */
      
	/** constructor */
	SolidElementT(FEManagerT& fe_manager);

	/* accessors */
	const LocalArrayT& LastDisplacements(void) const;
	const LocalArrayT& Velocities(void) const;
	const LocalArrayT& Accelerations(void) const;
	
	/** initialization. called immediately after constructor */
	virtual void Initialize(void);

	/** set the controller */
	virtual void SetController(eControllerT* controller);

	/** form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* solution calls */
	virtual void AddNodalForce(int node, dArrayT& force);
	virtual void AddLinearMomentum(dArrayT& momentum);
	
	/* returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void);

	/* compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);
	
protected:

	/** construct list of materials from the input stream */
	virtual void ReadMaterialData(ifstreamT& in);

	/* print element group data */
	virtual void PrintControlData(ostream& out) const;
	virtual void EchoOutputCodes(ifstreamT& in, ostream& out);

	/* initialization functions */
	virtual void SetLocalArrays(void);
	virtual void SetShape(void);

	/* form shape functions and derivatives */
	virtual void SetGlobalShape(void);

	/* construct the effective mass matrix */
	virtual void LHSDriver(void);
	void ElementLHSDriver(void);

	/* form the residual force vector */
	virtual void RHSDriver(void);
	void ElementRHSDriver(void);

	/* increment current element */
	virtual bool NextElement(void);	
	
	/* form the element stiffness matrix */
	virtual void FormStiffness(double constK);
	
	/* compute the effective acceleration and velocities based
	 * on the algorithmic flags formXx and the given constants
	 * constXx.
	 *
	 *		acc_eff  = constMa acc  + constCv a vel
	 *      vel_eff  = 0
	 *      disp_eff = constKd disp + constCv b vel
	 *
	 * where a and b are the Rayleigh damping coefficients.  No
	 * effective velocity since it's accounted for in the effective
	 * a and d.
	 *
	 * Note: In the process, the function collects the required
	 *       local arrays.
	 */
	virtual void ComputeEffectiveDVA(int formBody,
		int formMa, double constMa, int formCv, double constCv,
		int formKd, double constKd);

	/* body force */
	void FormRayleighMassDamping(double constM);

	/* damping force */
	virtual void FormCv(double constC);
	
	/* internal force */
	virtual void FormKd(double constK);

	/* return a pointer to a new material list */
	virtual MaterialListT* NewMaterialList(int size) const;

	/* driver for calculating output values */
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	                           const iArrayT& e_codes, dArray2DT& e_values);

	/** indicies of elements in the list of material needs */
	enum MaterialNeedsT {kNeedDisp = 0,
	                     kNeedVel  = 1,
	                 KNeedLastDisp = 2};

private:

	/* construct output labels array */
	virtual void SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;
	virtual void SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;
	virtual void GenerateOutputLabels(const iArrayT& n_counts, ArrayT<StringT>& n_labels, 
		const iArrayT& e_counts, ArrayT<StringT>& e_labels) const;

protected:

	/* control data */
	int	fMassType;	
	int	fStrainDispOpt;

	/* propagation direction for wave speeds */
	dArrayT fNormal;
	
	/* arrays with local ordering */
	LocalArrayT fLocLastDisp; /**< last converged displacements */
	LocalArrayT fLocVel;      /**< nodal velocities */
	LocalArrayT fLocAcc;      /**< nodal accelerations */

	/* run time */
	StructuralMaterialT*  fCurrMaterial;
	ArrayT<ArrayT<bool> > fMaterialNeeds;

	/* work space */
	dMatrixT    fD;      /* constitutive matrix        */
	dMatrixT    fB;      /* strain-displacement matrix */
	dSymMatrixT fStress; /* stress vector              */	

	/* parameters */
	static const int NumNodalOutputCodes;
	static const int NumElementOutputCodes;
};

/* accessors */
inline const LocalArrayT& SolidElementT::LastDisplacements(void) const { return fLocLastDisp; }
inline const LocalArrayT& SolidElementT::Velocities(void) const { return fLocVel; }
inline const LocalArrayT& SolidElementT::Accelerations(void) const { return fLocAcc; }

#endif /* _ELASTIC_T_H_ */
