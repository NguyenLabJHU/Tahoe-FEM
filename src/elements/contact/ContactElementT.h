/* $Id: ContactElementT.h,v 1.8 2001-06-28 23:57:14 rjones Exp $ */

#ifndef _CONTACT_ELEMENT_T_H_
#define _CONTACT_ELEMENT_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "pArrayT.h"
#include "LocalArrayT.h"
#include "dArray2DT.h"
#include "nVariArray2DT.h"
#include "ContactSurfaceT.h"
#include "ContactSearchT.h"


class ContactElementT: public ElementBaseT
{
public:

	/* constructor */
	ContactElementT(FEManagerT& fe_manager);

	/* destructor */
	virtual ~ContactElementT(void);

	/* form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* element level reconfiguration for the current solution */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

	/* initialization after constructor */
	virtual void Initialize(void);

	/* solution calls */
	virtual void AddNodalForce(int node, dArrayT& force); //not implemented

	/* Returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void); // not implemented
	
	/* writing output */
	virtual void RegisterOutput(void);
	virtual void WriteOutput(IOBaseT::OutputModeT mode);

	/* compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);  // not implemented

        /* append element equations numbers to the list */
        virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
                AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/* append connectivities */
	virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
		AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;
	virtual void ConnectsX(AutoArrayT<const iArray2DT*>& connects) const;
		// returns no (NULL) geometry connectivies
	 	
        /* surface specification modes */
        enum SearchParametersT { kGapTol = 0,
				kXiTol ,
				kSearchNumParameters};
        enum EnforcementParametersT { kPass = 0, // these are enf. specific!!
				kPenalty ,
				kConsistentTangent ,
				kSmithFerranteA,
				kSmithFerranteB,
				kEnfNumParameters};
	

protected:
	/* contact surfaces */
	ArrayT<ContactSurfaceT> fSurfaces; 

        /* search interaction parameters, symmetric matrix */
        nMatrixT<dArrayT> fSearchParameters ;

	// this will have kPass, kPenalty etc
	nMatrixT<dArrayT> fEnforcementParameters ;

	/* print element group data */
	virtual void PrintControlData(ostream& out) const;
	
	/* initialization steps */
	virtual void EchoConnectivityData(ifstreamT& in, ostream& out);
#if 0
	virtual void SetWorkSpace(void);
#endif

	/* generate contact element data  */
	bool SetContactConfiguration(void);

	/* update contact element data  */
	bool UpdateContactConfiguration(void);

	/* search pointer */
	ContactSearchT* fContactSearch ;

	/* link surfaces in ConnectsU - for graph */
	iArray2DT fSurfaceLinks;
	
private:
        /* surface specification modes */
        enum SurfaceSpecModeT { kSideSets = 1};
};

#endif /* _CONTACT_ELEMENT_T_H_ */
