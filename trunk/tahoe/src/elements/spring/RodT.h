/* $Id: RodT.h,v 1.4 2002-01-06 06:58:37 cbhovey Exp $ */
/* created: paklein (10/22/1996) */

#ifndef _ROD_T_H_
#define _ROD_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "LocalArrayT.h"
#include "RodMaterialT.h"

/* templates */
#include "pArrayT.h"

/** \note the RodT class doesn't provide complete support for the
 * different time integration schemes implemented using the
 * controller classes. need to add something like the
 * ComputeEffectiveDVA functions from the continuum element
 * classes to provide contributions to the global equations
 * which are consistent with the time integration algorithm.
 * PAK (05/30/1999) */
class RodT: public ElementBaseT
{
public:

	/** spring potential types */
	enum PotentialCodeT {
        kQuad = 1, /**< linear spring */
       kLJ612 = 2 /**< Lennard-Jones 6-12 potential */
		};

	/** constructor */
	RodT(FEManagerT& fe_manager);
	
	/** initialization */
	virtual void Initialize(void);

	/** form of tangent matrix */
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

	/* Element type parameters */
	static const int kRodTndof; /* number of degrees of freedom per node */
	static const int  kRodTnsd; /* number of spatial dimensions */
	 			  	
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

	/** return true if connectivities are changing */
	virtual bool ChangingGeometry(void) const { return false; };

protected:

	/** reference ID for sending output */
	int fOutputID;

	/* material data */
	pArrayT<RodMaterialT*> fMaterialsList; 	
	RodMaterialT*	       fCurrMaterial;

	LocalArrayT	fLocInitCoords; /* refcoords with local ordering */
	LocalArrayT fLocDisp;       /* displacements with local ordering */
};

#endif /* _ROD_T_H_ */
