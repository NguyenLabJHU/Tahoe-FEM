/* $Id: ContactT.h,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (12/11/1997)                                          */
/* Cohesive surface element group - non-directional cohesive              */
/* energy function, ie. fails in compression.                             */

#ifndef _CONTACT_T_H_
#define _CONTACT_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "pArrayT.h"
#include "LocalArrayT.h"
#include "dArray2DT.h"
#include "nVariArray2DT.h"

class ContactT: public ElementBaseT
{
public:

	/* constructor */
	ContactT(FEManagerT& fe_manager, int numfacetnodes);

	/* destructor */
	virtual ~ContactT(void);

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

	/* append connectivities */
	virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
		AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;
	virtual void ConnectsX(AutoArrayT<const iArray2DT*>& connects) const;
		// returns no (NULL) geometry connectivies
	 	
protected:

	/* surface specification modes */
	enum SurfaceSpecModeT {kNodesOnFacet = 0,
                               kSideSets = 1,
                           kBodyBoundary = 2};

	/* striker node specification */
	enum StrikerSpecModeT {kListStrikers = 0,
                           kSurfaceNodes = 1,
                            kAllStrikers = 2,
                          kContactBodies = 3};

	/* print element group data */
	virtual void PrintControlData(ostream& out) const;
	
	/* initialization steps */
	virtual void EchoConnectivityData(ifstreamT& in, ostream& out);
	virtual void SetWorkSpace(void);

	/* generate contact element data - return true if configuration has
	 * changed since the last call */
	bool SetContactConfiguration(void);

	/* steps in setting contact configuration */
	virtual bool SetActiveInteractions(void) = 0; // "internal" data
	virtual void SetConnectivities(void) = 0; // "external" data - interface to FEManager

	/* surface input functions */
	void InputNodesOnFacet(ifstreamT& in, iArray2DT& facets);
	void InputSideSets(ifstreamT& in, ostream& out, iArray2DT& facets);
	void InputBodyBoundary(ifstreamT& in, ArrayT<iArray2DT>& surfaces,
		int& surface);

	/* generate striker list from surfaces */
	void StrikersFromSurfaces(void);

private:

	/* read strikers */
	void ReadStrikers(ifstreamT& in, ostream& out);

protected:

	/* control parameters */
	int fNumFacetNodes;

	/* tags on contact surfaces */
	ArrayT<iArray2DT> fSurfaces; // lists of facets

	/* database info */
	iArrayT	  fStrikerTags; // should be variable
	dArray2DT fStrikerCoords; // should be variable
		// only used for search grid

	/* by-striker data */
	// could also map the strikers onto the variable sized data
	// for the active strikers.
	iArrayT fActiveMap; // map to active strikers
	AutoArrayT<int> fActiveStrikers; // global numbers of active strikers
	AutoArrayT<int> fHitSurface;     // contact surface for each active striker
	AutoArrayT<int> fHitFacets;      // facet of contact surface

	/* link surfaces in ConnectsU - for graph */
	iArray2DT fSurfaceLinks;
	
private:

	/* dynamic work space managers for element arrays */
	nVariArray2DT<int> fConnectivities_man;		
	nVariArray2DT<int> fEqnos_man;
	
//NOTE: need to generate overall list of facets (all same "geometry"?)
//      to return as ConnectsX for domain decompositition, and what to
//      due with strikers? For now do not return decomposition connects
};

#endif /* _CONTACT_T_H_ */
