/* $Id: StaggeredMultiScaleT.h,v 1.11 2003-02-03 04:40:20 paklein Exp $ */ 
#ifndef _STAGGERED_MULTISCALE_T_H_ 
#define _STAGGERED_MULTISCALE_T_H_ 

/* base classes */
#include "ElementBaseT.h"

/* direct members */
#include "LocalArrayT.h"
#include "GeometryT.h"

/* base multiscale classes */
#include "FEA.h"
#include "VMS.h"
#include "FEA_FormatT.h"

namespace Tahoe {

/* forward declarations */
class ShapeFunctionT;

/** StaggeredMultiScaleT: This class contains methods pertaining to kinematics of
 * a dual field formulation. These include deformation gradients Fe and Fp
 * and gradients such as Grad(ue) and Grad(up) as examples.
 * Sandia National Laboratory and the University of Michigan **/

class StaggeredMultiScaleT: public ElementBaseT 
{
 //----- Class Methods -----------
	
 public:

	/** constructor */
	StaggeredMultiScaleT(const ElementSupportT& support, const FieldT& coarse, 
		const FieldT& fine);

	/** destructor */
	~StaggeredMultiScaleT(void);

	/** data initialization */
	virtual void Initialize(void); 

	/** return true if the element contributes to the solution of the
	 * given group. ElementBaseT::InGroup returns true if group is the
	 * same as the group of the FieldT passed in to ElementBaseT::ElementBaseT. */
	virtual bool InGroup(int group) const;

	/** close current time increment. Called if the integration over the
	 * current time increment was successful. */
	virtual void CloseStep(void);

	/** collecting element group equation numbers. See ElementBaseT::Equations
	 * for more information */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/** return a const reference to the run state flag */
	virtual GlobalT::SystemTypeT TangentType(void) const;
	
	/** accumulate the residual force on the specified node
	 * \param node test node
	 * \param force array into which to assemble to the residual force */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);
	
	/** returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void);

	/** write element group parameters to out */
	virtual void PrintControlData(ostream& out) const;

	/** register element for output */
	virtual void RegisterOutput(void);

	/** write element output */
	virtual void WriteOutput(void);

	/** compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);
	/*@}*/

	/** \name restart functions */
	/*@{*/
	/** write restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::ReadRestart implementation. */
	virtual void WriteRestart(ostream& out) const;

	/** read restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::WriteRestart implementation. */
	virtual void ReadRestart(istream& in);
	/*@}*/

protected:

	/** \name drivers called by ElementBaseT::FormRHS and ElementBaseT::FormLHS */
	/*@{*/
	/** form group contribution to the stiffness matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT);

	/** form group contribution to the residual */
	virtual void RHSDriver(void);
	/*@}*/

	void Select_Equations ( const int &iCoarseScale, const int &iFineScale );

private:

	//----- Class Data -----------
	
public:	
protected:
private:

	/** Data at time steps n and n+1 used by both Coarse and Fine */
	//VMS_VariableT n,np1; // <-- keep local scope in elmt loop for now 

	/** Gradients with respect to reference coodinates */
	FEA_dMatrixT fGRAD_ua, fGRAD_ua_n, fGRAD_ub, fGRAD_ub_n, fSigma;

	/** \name  values read from input in the constructor */
	/*@{*/
	/** element geometry */
	GeometryT::CodeT fGeometryCode;

	/** number of integration points */
	int	fNumIP;
	/*@}*/

	/** \name element displacements in local ordering */
	/*@{*/
	LocalArrayT ua;     /**< fine scale displacement */
	LocalArrayT ua_n; 	/**< fine scale displacement from previous increment */
	LocalArrayT del_ua; /**< the Newton-R update i.e. del_ua = ua - ua_n (ua subcript n+1 implied) */
	LocalArrayT ub;     /**< coarse scale displacement */
	LocalArrayT ub_n; 	/**< coarse scale displacement from previous increment */
	LocalArrayT del_ub; /**< the Newton-R update i.e. del_ub = ub - ub_n (ub subcript n+1 implied) */
	dArrayT			del_ua_vec;  	/** need in vector for i.e. { {ua1_1,ua1_2},{ua2_1,ua2_2}, ... } */
	dArrayT			del_ub_vec;		/** need in vector for i.e. { {ub1_1,ub1_2},{ub2_1,ub2_2}, ... } */
	/*@}*/
	
	/** \name shape functions wrt to current coordinates */
	/*@{*/
	/** shape functions and derivatives. The derivatives are wrt to the 
	 * coordinates in StaggeredMultiScaleT::fCurrCoords, which are the
	 * current coordinates */
	ShapeFunctionT* fShapes;
	
	FEA_ShapeFunctionT fFEA_Shapes;

	/** reference coordinates */
	LocalArrayT fInitCoords;     

	/** current coordinates */
	LocalArrayT fCurrCoords;
	/*@}*/

	/** the BLACK BOXS */
	/** VMF is acronym for "Variational Multi-Field". Class VMF_Virtual_WorkT is simply
	 * linearizaiton of the Virtual Work Equaiton, using a decomposiiton u = ua + ub.  This
	 * class can be used for any multi-field formulation where u is decomposed as stated, 
	 * (to include Variational Multi-Scale (VMS) formulation). */

	/* Data Storage */
	dMatrixT 	fKa_I, 	fKb_I;  
	dMatrixT 	fKa_II, fKb_II; 
	dArrayT 	fFint_I;
	dArrayT		fFint_II;

	/* Multi-Field Element Formulators */
	CoarseScaleT* fEquation_I;	
	FineScaleT* 	fEquation_II;

	/* Multi-Field Materials */
	VMF_MaterialT* fFineMaterial;
	VMF_MaterialT* fCoarseMaterial;

	/* Conversion methods: puts data in FEA format (very little cost in perspective) */
	FEA_FormatT Convert;

	/** the coarse scale field */
	const FieldT& fCoarse;
	
	/** the fine scale field */
	const FieldT& fFine;	

	/** equations per element for the fine scale. The coarse scale equations are
	 * in ElementBaseT::fEqnos and are handled by ElementBaseT. */
	iArray2DT fEqnos_fine;

	/** \name state variable storage *
	 * State variables are handled ABAQUS-style. For every iteration, the state 
	 * variables from the previous increment are passed to the element, which 
	 * updates the values in place. Each row in the array is the state variable
	 * storage for all integration points for an element */
	/*@{*/
	dArray2DT fdState_new;
	dArray2DT fdState;

	iArray2DT fiState_new;
	iArray2DT fiState;
	/*@}*/
	
	/** \name output */
	/*@{*/
	/** output ID */
	int fOutputID;
	
	/** integration point stresses. Calculated and stored during 
	 * StaggeredMultiScaleT::RHSDriver */
	dArray2DT fIPStress;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _STAGGERED_MULTISCALE_T_H_ */



