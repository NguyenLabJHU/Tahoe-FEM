/* $Id: SWDiamondT.h,v 1.8 2004-07-15 08:30:17 paklein Exp $ */
/* created: paklein (03/19/1997) */

#ifndef _SWDIAMOND_T_H_
#define _SWDIAMOND_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "LocalArrayT.h"
#include "RodMaterialT.h"

/* templates */
#include "pArrayT.h"


namespace Tahoe {

class SWDiamondT: public ElementBaseT
{
public:

	/* constructor */
	SWDiamondT(const ElementSupportT& support, const FieldT& field);

	/* form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* NOT implemented. Returns an zero force vector */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);

	/* returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void);

	/* append element equations numbers to the list */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);
	
	/* writing output */
	virtual void RegisterOutput(void);
	virtual void WriteOutput(void);

	/* compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);
		//TEMP: for now, does nothing
	 			  	
protected: /* for derived classes only */
	 	
	/* called by FormRHS and FormLHS */
	virtual void LHSDriver(GlobalT::SystemTypeT);
	virtual void RHSDriver(void);

	/* print element group data */
	virtual void PrintControlData(ostream& out) const;
		
	/* element data */
	virtual void ReadMaterialData(ifstreamT& in);	
	virtual void WriteMaterialData(ostream& out) const;
	virtual void EchoConnectivityData(ifstreamT& in, ostream& out);

	/* call AFTER 2 and 3 body node lists are set */
	virtual void ConfigureElementData(void);

	/* element list increment */
	virtual bool Next2Body(void);
	virtual bool Next3Body(void);
	
	/* element calculations */
	double Energy3Body(void);
	void Force3Body(void);
	void Stiffness3Body(void);

	double Energy2Body(void);
	void Force2Body(void);
	void Stiffness2Body(void);

	/* print connectivity element data */
	void PrintConnectivityData(ostream& out);

private: /* potential functions and their derivatives */

	/* 3 body potentials */
	double U3body(double r1, double r2, double c12) const;

		/* 1st derivs */
	double Dr1U3body(double r1, double r2, double c12) const;
	double Dr2U3body(double r1, double r2, double c12) const;
	double Dc12U3body(double r1, double r2, double c12) const;

		/* 2nd derivs */
	double DDr1U3body(double r1, double r2, double c12) const;
	double DDr2U3body(double r1, double r2, double c12) const;
	double DDc12U3body(double r1, double r2, double c12) const;
		
		/* mixed derivs */
	double Dr1Dr2U3body(double r1, double r2, double c12) const;
	double Dr1Dc12U3body(double r1, double r2, double c12) const;
	double Dr2Dc12U3body(double r1, double r2, double c12) const;

	/* compute entire Hessian at once */
	void ComputeHessian_3Body(dMatrixT& hessian,
		double r1, double r2, double c12);

	/* 2 body potential and derivatives */
	double U2body(double r) const;
	double DU2body(double r) const;
	double DDU2body(double r) const;

protected:

	/* unit scaling */
	double	feps;

	/* 2 body potential */
	double	fA;
	double	fdelta;
	 	
	/* 3 body potential */
	double	fgamma;
	double	flambda;
	
	double	frcut;		/* cut-off distance */
	double	fa;			/* lattice spacing parameter */
	
	/* derived values */
	double	fB;

	/* I/O ID */
	int       fOutputID;
	iArrayT   fNodesUsed;
	iArray2DT fOutputConnects;

	/* 3 Body */
	ElementMatrixT& fK_3Body;
	dArrayT&        fF_3Body;

	LocalArrayT	fLocX_3Body;
	LocalArrayT fLocd_3Body;

	AutoArrayT<ElementCardT>& List_3Body;
	
	iArray2DT fNodes_3Body;
	iArray2DT fEqnos_3Body;

	/* 2 Body */
	ElementMatrixT fK_2Body;
	dArrayT        fF_2Body;

	LocalArrayT	fLocX_2Body;
	LocalArrayT fLocd_2Body;

	AutoArrayT<ElementCardT> List_2Body;

	iArray2DT	fNodes_2Body;		
	iArray2DT	fEqnos_2Body;			

	/* work space */
	dMatrixT	fHessian_3Body; //Hessian of 3-body potential wrt {r1,r2,c12}

};

} // namespace Tahoe 
#endif /* _SWDIAMOND_T_H_ */
