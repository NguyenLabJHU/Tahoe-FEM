/* $Id: ParticleT.h,v 1.2 2002-10-20 22:48:27 paklein Exp $ */
#ifndef _ROD_T_H_
#define _ROD_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "LocalArrayT.h"
#include "C1FunctionT.h"

/* templates */
#include "pArrayT.h"

/** base class for particle types */
class ParticleT: public ElementBaseT
{
public:

	/** constructor */
	ParticleT(const ElementSupportT& support, const FieldT& field);
	
	/** initialization. Completely overrides ElementBaseT::Initialize */
	virtual void Initialize(void);

	/** form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* NOT implemented. Returns an zero force vector */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);
			
	/* returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void);
	
	/* writing output */
	virtual void RegisterOutput(void);
	virtual void WriteOutput(IOBaseT::OutputModeT mode);
//NOTE: for parallel calculations - write what you've got.

	/* compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);
	 			  	
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

	/** local to global tag map. Neighborlists constructed using
	 * group-local numbering */
	iArrayT GlobalTag;

private:

	/** \name work space */
	/*@{*/
	/** constant matrix needed to compute the stiffness */
	dMatrixT fOneOne;

	/** pair vector */
	dArrayT fBond;

	/** current coordinates for one pair bond */
//	dArray2DT fPairCoords;
	/*@}*/
};

#endif /* _ROD_T_H_ */
