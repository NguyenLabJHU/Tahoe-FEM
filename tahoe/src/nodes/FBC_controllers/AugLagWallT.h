/* $Id: AugLagWallT.h,v 1.1 2001-09-11 23:43:48 paklein Exp $ */

#ifndef _AUGLAG_WALL_T_H_
#define _AUGLAG_WALL_T_H_

/* base classes */
#include "PenaltyWallT.h"
#include "DOFElementT.h"

/* forward declarations */
class XDOF_ManagerT;

/** rigid barrier with augmented Lagrangian enforcement of
 * non-interpenetration */
class AugLagWallT: public PenaltyWallT, public DOFElementT
{
public:

	/* constructor */
	AugLagWallT(FEManagerT& fe_manager, XDOF_ManagerT* XDOF_nodes, const iArray2DT& eqnos,
		const dArray2DT& coords, const dArray2DT* vels);

	/* initialize data */
	virtual void Initialize(void);
	virtual void SetEquationNumbers(void);

	/* append element equations numbers to the list */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);
	virtual void Connectivities(AutoArrayT<const iArray2DT*>& connects_1,
		AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;

	/* restarts */
	virtual void ReadRestart(istream& in);
	virtual void WriteRestart(ostream& out) const;

	/* finalize step */
	virtual void CloseStep(void);

	/* tangent term */
	virtual void ApplyLHS(void);

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
	
private:

	/* accumulate the contact force vector fContactForce */
	virtual void ComputeContactForce(double kforce);

private:

	/* nodemanager(s) */
	XDOF_ManagerT* fXDOF_Nodes;
	
	/* contact equation sets (shallow copy of contact node equs) */
	iArray2DT fContactEqnos2D;
	iArray2DT fContactTags;
	
	/* contact DOF tags and DOF's */
	iArrayT fContactDOFtags;
	dArrayT fLastDOF;
};

#endif /* _AUGLAG_WALL_T_H_ */
