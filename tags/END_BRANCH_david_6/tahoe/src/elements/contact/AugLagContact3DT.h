/* $Id: AugLagContact3DT.h,v 1.3 2004-07-15 08:26:08 paklein Exp $ */
#ifndef _AUGLAG_CONTACT3D_T_H_
#define _AUGLAG_CONTACT3D_T_H_

/* base classes */
#include "Contact3DT.h"
#include "DOFElementT.h"

/* direct members */
#include "AutoArrayT.h"

namespace Tahoe {

/** contact enforcement in 3D using an augmented Lagrangian formulation.
 * Formulation by J. Heegaard and A. Curnier, IJNME \b 36, 569-593, 1993. */
class AugLagContact3DT: public Contact3DT, public DOFElementT
{
public:

	/** constructor */
	AugLagContact3DT(const ElementSupportT& support);

	/** \name implementation of the DOFElementT interface */
	/*@{*/
	/* append element equations numbers to the list */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);
	
	/* returns the array for the DOF tags needed for the current config */
	virtual void SetDOFTags(void);
	virtual iArrayT& DOFTags(int tag_set);

	/* generate nodal connectivities */
	virtual void GenerateElementData(void);

	/* return the contact elements */
	virtual const iArray2DT& DOFConnects(int tag_set) const;

	/* restore the DOF values to the last converged solution */
	virtual void ResetDOF(dArray2DT& DOF, int tag_set) const;

	/* returns 1 if group needs to reconfigure DOF's, else 0 */
	virtual int Reconfigure(void);

	/** restore any state data to the previous converged state */
	virtual void ResetState(void) { };

	/* the group */
	virtual int Group(void) const { return Contact3DT::Group(); }
	/*@}*/

	/* element level reconfiguration for the current solution */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

	/* append connectivities */
	virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
		AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;

	/** \name restart functions
	 * 	\note restarts have not been tested. these functions throw 
	 * ExceptionT::kGeneralFail. */
	/*@{*/
	virtual void ReadRestart(istream& in);
	virtual void WriteRestart(ostream& out) const;
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	
protected:

	/** step in setting contact configuration. Intercepted here so that
	 * the last contact configuration can be stored */
	virtual bool SetActiveInteractions(void);

	/* construct the effective mass matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

	/* construct the residual force vector */
	virtual void RHSDriver(void);
	
protected:

	/** regularization parameter */
	double fr;

	/* extended interaction data */
	iArray2DT fXDOFConnectivities;
	iArray2DT fXDOFEqnos;
	
	/* contact DOF tags */
	iArrayT fContactDOFtags; // VARIABLE

	/* constraint history */
	iArrayT fLastActiveMap; // VARIABLE
	dArrayT fLastDOF;       // VARIABLE

private:

	/** \name dynamic work space managers for element arrays */
	/*@{*/
	nVariArray2DT<int> fXDOFConnectivities_man;		
	nVariArray2DT<int> fXDOFEqnos_man;
	/*@}*/

	/** \name element coords and displacements */
	/*@{*/
	dArray2DT fElCoord;
	dArray2DT fElRefCoord;
	dArray2DT fElDisp;
	/*@}*/
	
	/** \name work space */
	/*@{*/
	dMatrixT fdc_du;
	dMatrixT fdn_du;
	dMatrixT fM1;
	dMatrixT fM2;
	dArrayT  fV1;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _AUGLAG_CONTACT3D_T_H_ */
