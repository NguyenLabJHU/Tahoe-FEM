/* $Id: ContactSurfaceT.h,v 1.29 2003-12-20 01:22:14 rjones Exp $ */
#ifndef _CONTACT_SURFACE_T_H_
#define _CONTACT_SURFACE_T_H_

/* base class */
#include "SurfaceT.h"

/* direct members */
#include "ArrayT.h"
#include "Array2DT.h"
#include "dArray2DT.h"
#include "nMatrixT.h"

namespace Tahoe {

/* forward declarations */
class ofstreamT;
class ContactNodeT;

/** 
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
	void Initialize(const ElementSupportT& support, int num_multipliers);

	/* potential connectivities based on growing/sliding contact */
	void SetPotentialConnectivity(void);

	/* potential connectivities for multipliers */
	void SetMultiplierConnectivity(void);


	/* access functions */
	inline ArrayT<ContactNodeT*>& ContactNodes(void) 
		{return fContactNodes;}
	inline RaggedArray2DT<int>& Connectivities(void) { return fConnectivities; }
	inline const RaggedArray2DT<int>& Connectivities(void) const { return fConnectivities; }
	inline RaggedArray2DT<int>& EqNums(void)
		{return fEqNums;}  // this can NOT be const
	bool IsInConnectivity
		(int primary_local_node, int secondary_global_node) const;

	/* output */
	void CollectOutput(iArrayT& OutputFlags, dArray2DT& values) const;
	void PrintContactArea(ostream& out) const;
	void PrintGaps(ostream& out) const;
	void PrintGaps(ofstream& out) const;
	void PrintNormals(ostream& out) const;
	void PrintNormals(ofstream& out) const;
	void PrintStatus(ostream& out) const;
	void PrintMultipliers(ostream& out) const;
	void PrintMultipliers(ofstream& out) const;


	void InitializeMultiplierMap(void);
	void DetermineMultiplierExtent(void);
	void TagMultiplierMap(const ArrayT<FaceT*>&  faces);
	inline iArrayT&  MultiplierTags(void) 
		{return fMultiplierTags;} // this can NOT be const	
	inline const iArrayT&  MultiplierTags(void) const
		{return fMultiplierTags;} 
	inline const iArrayT&  MultiplierMap(void) const
		{return fMultiplierMap;}	
	void AllocateMultiplierTags(void);
	void ResetMultipliers(dArray2DT& multiplier_values) const;
	void MultiplierTags
		(const iArrayT& local_nodes, iArrayT& multiplier_tags) const;
	void MultiplierValues
		(const iArrayT& local_nodes, nArray2DT<double*>& multiplier_values) const;
	void MultiplierValues
		(const iArrayT& local_nodes, dArray2DT& multiplier_values) const;
	const iArray2DT& DisplacementMultiplierNodePairs(void);
	inline bool HasMultiplier(int i) 
		{return fMultiplierMap[i] > -1;} 
	inline void AliasMultipliers(const dArray2DT& multipliers)
		{fMultiplierValues.Alias(multipliers);}
	inline double& Multiplier(int tag, int i) 
		{return fMultiplierValues(fMultiplierMap[tag],i);} 

  protected:
	/* contact nodes */
	ArrayT <ContactNodeT*>  fContactNodes ; 

	int fNumMultipliers;
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
	dArray2DT fMultiplierValues; 
	/* global multiplier "node" tags for active nodes */
	iArrayT fMultiplierTags; 
	/* hash for local u node to active p nodes */
	iArrayT fMultiplierMap; 
	/* multiplier history */
	iArrayT fLastMultiplierMap; 
	dArray2DT fLastMultiplierValues; 
	/* u-x pairs for ContactElementT::DOFConnects */
	iArray2DT fDisplacementMultiplierNodePairs;


};

} // namespace Tahoe 
#endif /* _CONTACT_SURFACE_T_H_ */
