/* $Id: ContactSurfaceT.h,v 1.14 2001-09-14 00:27:16 rjones Exp $ */


#ifndef _CONTACT_SURFACE_T_H_
#define _CONTACT_SURFACE_T_H_

/* base class */
#include "SurfaceT.h"

/* direct members */
#include "ArrayT.h"
#include "SurfaceT.h"
#include "ContactNodeT.h"
#include "nMatrixT.h"

/* forward declarations */
class ofstreamT;
class FEManagerT;
//class nMatrixT<dArrayT>;

/* 
a ContactSurface will only have one opposing face per
node and be considered "smooth" i.e. the full boundary 
surface of a cube will be made up of 6 surfaces
*/

class ContactSurfaceT : public SurfaceT
{
  public:
	/* constructor */
	ContactSurfaceT(void);

	/* destructor */
	~ContactSurfaceT(void);

	/* allocate contact node array */
	void Initialize(const NodeManagerT* node_manager);

	/* set contact status */
	void SetContactStatus(nMatrixT<dArrayT>& enforcement_parameters);
	void UpdateContactStatus(nMatrixT<dArrayT>& enforcement_parameters);

	/* potential connectivities based on growing/sliding contact */
	void SetPotentialConnectivity(int num_multipliers);

	/* access functions */
	inline ArrayT<ContactNodeT*>& ContactNodes(void) 
		{return fContactNodes;}
	inline RaggedArray2DT<int>& Connectivities(void)
		{return fConnectivities;}
	inline RaggedArray2DT<int>& EqNums(void)
		{return fEqNums;}  // this can NOT be const
	bool IsInConnectivity
		(int primary_local_node, int secondary_global_node) const;

	void PrintContactArea(ostream& out) const;
	void PrintGap(ostream& out) const;
	void PrintGap(ofstream& out) const;
	void PrintNormals(ofstream& out) const;
	void PrintStatus(ostream& out) const;

	inline void InitializeMultiplierMap(void)
		{fMultiplierMap = -1;}
	void DetermineMultiplierExtent(void);
	void TagMultiplierMap(const ArrayT<FaceT*>&  faces);
	inline iArrayT&  MultiplierTags(void) 
		{return fMultiplierTags;} // this can NOT be const	
	inline const iArrayT&  MultiplierTags(void) const
		{return fMultiplierTags;} 
	inline const iArrayT&  MultiplierMap(void) const
		{return fMultiplierMap;}	
	void AllocateMultiplierTags(dArray2DT& multiplier_values);
	void ResetMultipliers(dArray2DT& multiplier_values);
	void MultiplierTags(iArrayT& local_nodes, iArrayT& multiplier_tags);
	iArray2DT& RealGhostNodePairs(void);


  protected:
        /* nodal arrays */
	ArrayT <ContactNodeT*>  fContactNodes ; 

	int fNumPotentialContactNodes;

	/* potential connectivities for the time step */
	RaggedArray2DT<int> fConnectivities;

	/* space for associated equation numbers */
	RaggedArray2DT<int> fEqNums;
#if 0
	/* for frictional slip */
	ArrayT <ContactNodeT*>  fPreviousContactNodes;
#endif
	/* Multiplier Data, which is variable size */
	/* global multiplier "node" tags for active nodes */
	iArrayT fMultiplierTags; 
	/* hash for local node to active nodes */
	iArrayT fMultiplierMap; 
	/* multiplier history */
	iArrayT fLastMultiplierMap; 
	dArray2DT fLastMultiplierValues; 
	iArray2DT fRealGhostNodePairs;


};

#endif /* _CONTACT_SURFACE_T_H_ */
