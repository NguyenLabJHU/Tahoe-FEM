/* $Id: SolidElementT.h,v 1.1.1.1 2001-01-29 08:20:39 paklein Exp $ */
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

class SolidElementT: public ContinuumElementT
{
public:
	
	enum OutputCodeT {iNodalCoord = 0, // (reference) nodal coordinates
                       iNodalDisp = 1, // nodal displacements
                     iNodalStress = 2, // nodal stresses
                       iPrincipal = 3, // principal stresses
                   iEnergyDensity = 4, // nodal strain energy density
                      iWaveSpeeds = 5, // wave speeds
                    iMaterialData = 6};// material model output

	/* constructor */
	SolidElementT(FEManagerT& fe_manager);

	/* accessors */
	const LocalArrayT& LastDisplacements(void) const;
	const LocalArrayT& Velocities(void) const;
	const LocalArrayT& Accelerations(void) const;
	
	/* data initialization */
	virtual void Initialize(void);

	/* set the controller */
	virtual void SetController(eControllerT* controller);

	/* form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* solution calls */
	virtual void AddNodalForce(int node, dArrayT& force);
	virtual void AddLinearMomentum(dArrayT& momentum);
	
	/* returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void);

	/* compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);

protected:

	/* print element group data */
	virtual void PrintControlData(ostream& out) const;
	virtual void EchoOutputCodes(ifstreamT& in, ostream& out);

	/* initialization functions */
	virtual void SetLocalArrays(void);
	virtual void SetShape(void);

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

	/* driver for nodal value calculations */
	virtual void ComputeNodalValues(const iArrayT& codes);

private:

	/* construct output labels array */
	virtual void SetOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts);
	virtual void GenerateOutputLabels(const iArrayT& counts,
		ArrayT<StringT>& labels) const;

protected:

	/* control data */
	int	fMassType;	
	int	fStrainDispOpt;

	/* propagation direction for wave speeds */
	dArrayT fNormal;
	
	/* arrays with local ordering */
	LocalArrayT fLocLastDisp; // last converged disp's local ordering
	LocalArrayT fLocVel;      // velocities with local ordering
	LocalArrayT fLocAcc;      // accelerations with local ordering
// NOTE: compute approximate velocities and accelerations for the
//       for quasi-static case????

	/* run time */
	StructuralMaterialT* fCurrMaterial;

	/* work space */
	dMatrixT    fD;      /* constitutive matrix        */
	dMatrixT    fB;      /* strain-displacement matrix */
	dSymMatrixT fStress; /* stress vector              */	

	/* parameters */
	static const int NumOutputCodes;
};

/* accessors */
inline const LocalArrayT& SolidElementT::LastDisplacements(void) const { return fLocLastDisp; }
inline const LocalArrayT& SolidElementT::Velocities(void) const { return fLocVel; }
inline const LocalArrayT& SolidElementT::Accelerations(void) const { return fLocAcc; }

#endif /* _ELASTIC_T_H_ */
