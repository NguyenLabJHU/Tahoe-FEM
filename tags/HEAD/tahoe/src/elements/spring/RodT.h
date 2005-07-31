/* $Id: RodT.h,v 1.1.1.1 2001-01-29 08:20:34 paklein Exp $ */
/* created: paklein (10/22/1996)                                          */

#ifndef _ROD_T_H_
#define _ROD_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "LocalArrayT.h"
#include "RodMaterialT.h"

/* templates */
#include "pArrayT.h"

class RodT: public ElementBaseT
{
public:

	/* constructor */
	RodT(FEManagerT& fe_manager);
	
	/* initialization */
	virtual void Initialize(void);

	/* form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* NOT implemented. Returns an zero force vector */
	virtual void AddNodalForce(int node, dArrayT& force);
			
	/* returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void);
	
	/* writing output */
	virtual void RegisterOutput(void);
	virtual void WriteOutput(IOBaseT::OutputModeT mode);

	/* compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);
		//TEMP: for now, does nothing

	/* Element type parameters */
	static const int kRodTndof = 2; /* number of degrees of freedom per node */
	static const int  kRodTnsd = 2; /* number of spatial dimensions */
	 			  	
protected: /* for derived classes only */
	 	
	/* called by FormRHS and FormLHS */
	virtual void LHSDriver(void);
	virtual void RHSDriver(void);

	/* increment current element */
	virtual bool NextElement(void);
		
	/* element data */
	virtual void ReadMaterialData(ifstreamT& in);	
	virtual void WriteMaterialData(ostream& out) const;
	
	/* element calculations */
	double ElementEnergy(void);
	void ElementForce(double constKd);
	void ElementStiffness(double constK);

protected:

	/* material data */
	pArrayT<RodMaterialT*> fMaterialsList; 	
	RodMaterialT*	       fCurrMaterial;

	LocalArrayT	fLocInitCoords; /* refcoords with local ordering */
	LocalArrayT fLocDisp;       /* displacements with local ordering */

};

#endif /* _ROD_T_H_ */
