/* $Id: AugLagSphereT.h,v 1.9 2004-01-05 07:12:40 paklein Exp $ */
/* created: paklein (03/24/1999) */
#ifndef _AUGLAG_SPHERE_T_H_
#define _AUGLAG_SPHERE_T_H_

/* base classes */
#include "PenaltySphereT.h"
#include "DOFElementT.h"

namespace Tahoe {

/* forward declarations */
class XDOF_ManagerT;
class FieldT;

class AugLagSphereT: public PenaltySphereT, public DOFElementT
{
public:

	/* constructor */
	AugLagSphereT(FEManagerT& fe_manager, XDOF_ManagerT* XDOF_nodes, 
		const FieldT& field, const dArray2DT& coords, const dArray2DT& disp);

	/* initialize data */
	virtual void Initialize(void);
	virtual void SetEquationNumbers(void);

	/* append element equations numbers to the list */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);
	virtual void Connectivities(AutoArrayT<const iArray2DT*>& connects_1,
		AutoArrayT<const RaggedArray2DT<int>*>& connects_2,
		AutoArrayT<const iArray2DT*>& equivalent_nodes) const;

	/* restarts */
	virtual void ReadRestart(istream& in);
	virtual void WriteRestart(ostream& out) const;

	/* finalize step */
	virtual void CloseStep(void);

	/** tangent
	 * \param sys_type "maximum" tangent type needed by the solver. The GlobalT::SystemTypeT
	 *        enum is ordered by generality. The solver should indicate the most general
	 *        system type that is actually needed. */
	virtual void ApplyLHS(GlobalT::SystemTypeT sys_type);

	/* returns the array for the DOF tags needed for the current config */
	virtual void SetDOFTags(void);
	virtual iArrayT& DOFTags(int tag_set);

	/* generate nodal connectivities - does nothing here */
	virtual void GenerateElementData(void);

	/* return the contact elements */
	virtual const iArray2DT& DOFConnects(int tag_set) const;

	/* restore the DOF values to the last converged solution */
	virtual void ResetDOF(dArray2DT& DOF, int tag_set) const;

	/* returns 1 if group needs to reconfigure DOF's, else 0 */
	virtual int Reconfigure(void);

	/** restore any state data to the previous converged state */
	virtual void ResetState(void) { };

	/** return the equation group to which the generate degrees of
	 * freedom belong. */
	virtual int Group(void) const;	
	
private:

	/* accumulate the contact force vector fContactForce */
	virtual void ComputeContactForce(double kforce);

private:

	/** XDOF manager */
	XDOF_ManagerT* fXDOF_Nodes;

	/** the field */
	const FieldT& fField;
	
	/* contact equation sets (shallow copy of contact node equs) */
	iArray2DT fContactEqnos2D;
	iArray2DT fContactTags;
	
	/* contact DOF tags and DOF's */
	iArrayT fContactDOFtags;
	dArrayT fLastDOF;
};

} // namespace Tahoe 
#endif /* _AUGLAG_SPHERE_T_H_ */
