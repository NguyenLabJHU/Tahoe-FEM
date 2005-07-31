/* $Id: AugLagContact2DT.h,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (05/31/1998)                                          */

#ifndef _AUGLAG_CONTACT2D_T_H_
#define _AUGLAG_CONTACT2D_T_H_

/* base classes */
#include "Contact2DT.h"
#include "DOFElementT.h"

/* direct members */
#include "AutoArrayT.h"

/* forward declarations */
class iGridManager2DT;
class XDOF_ManagerT;

class AugLagContact2DT: public Contact2DT, public DOFElementT
{
public:

	/* constructor */
	AugLagContact2DT(FEManagerT& fe_manager, XDOF_ManagerT* XDOF_nodes);

	/* allocates space and reads connectivity data */
	virtual void Initialize(void);

	/* append element equations numbers to the list */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);
	
	/* returns the array for the DOF tags needed for the current config */
	virtual iArrayT& SetDOFTags(void);
	virtual const iArrayT& DOFTags(void) const;

	/* generate nodal connectivities */
	virtual void GenerateElementData(void);
	// NOTE: since the sequence of setting global equation
	//       number is controlled externally, responsibility
	//       for calling the element group to (self-) configure
	//       is also left to calls from the outside. otherwise
	//       it's tough to say whether data requested by the group
	//       is current.

	/* return the contact elements */
	virtual const iArray2DT& DOFConnects(void) const;

	/* restore the DOF values to the last converged solution */
	virtual void ResetDOF(dArray2DT& DOF) const;

	/* returns 1 if group needs to reconfigure DOF's, else 0 */
	virtual int Reconfigure(void);

	/* element level reconfiguration for the current solution */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

	/* append connectivities */
	virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
		AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;

	/* restart functions */
	virtual void ReadRestart(istream& in);
	virtual void WriteRestart(ostream& out) const;
	//TEMP - restarts have not been tested. these functions
	//       throw exceptions
		 	
protected:

	/* print element group data */
	virtual void PrintControlData(ostream& out) const;
		 	
	/* construct the effective mass matrix */
	virtual void LHSDriver(void);

	/* construct the residual force vector */
	virtual void RHSDriver(void);
	
protected:

	double fr; // regularization parameter

	/* extended interaction data */
	iArray2DT fXDOFConnectivities;
	iArray2DT fXDOFEqnos;
	
	/* nodemanager with external DOF's */
	XDOF_ManagerT* fXDOF_Nodes;
	
	/* contact DOF tags */
	iArrayT fContactDOFtags; // VARIABLE

	/* constraint history */
	iArrayT fLastActiveMap; // VARIABLE
	dArrayT fLastDOF;       // VARIABLE

private:

	/* dynamic work space managers for element arrays */
	nVariArray2DT<int> fXDOFConnectivities_man;		
	nVariArray2DT<int> fXDOFEqnos_man;

};

#endif /* _AUGLAG_CONTACT2D_T_H_ */
