/* $Id: ParticleT.h,v 1.5 2002-11-22 01:49:45 paklein Exp $ */
#ifndef _PARTICLE_T_H_
#define _PARTICLE_T_H_

/* base class */
#include "ElementBaseT.h"

namespace Tahoe {

/** forward declarations */
class iGridManagerT;
class PotentialT;

/** base class for particle types */
class ParticleT: public ElementBaseT
{
public:

	/** constructor */
	ParticleT(const ElementSupportT& support, const FieldT& field);

	/** destructor */
	~ParticleT(void);
	
	/** initialization */
	virtual void Initialize(void);

	/** form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* NOT implemented. Returns an zero force vector */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);
			
	/* returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void) { return 0.0; };
	
	/* writing output */
	virtual void RegisterOutput(void);
	virtual void WriteOutput(void);

	/* compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);
	 			  	
protected: /* for derived classes only */

	/** echo element connectivity data. Reads parameters that define
	 * which nodes belong to this ParticleT group. */
	virtual void EchoConnectivityData(ifstreamT& in, ostream& out);
	
	/** generate labels for output data */
	virtual void GenerateOutputLabels(ArrayT<StringT>& labels) const;

	/** return true if connectivities are changing */
	virtual bool ChangingGeometry(void) const;

	/** generate neighborlist
	 * \param particle_tags global tags for which to determine neighhors
	 * \param distance distance over which to search for neighbors
	 * \param double_list if true the neighbor lists will contain two references for
	 *        every neighbor interaction. For tags A and B, this means B will appear
	 *        in the neighbor list for A and vice versa. if true, neighbor lists will
	 *        only contain neighbors for which A > B.
	 * \param neighbors list of neighbors for every tag in particle_tags as rows in the
	 *        2D array. The list for each tag begins with the tag itself. Therefore, all
	 *        lists are at least length 1. The numbering of the tags in the neighbor list
	 *        is in terms of the tags provided in particle_tags. */
	void GenerateNeighborList(const iArrayT& particle_tags, double distance, bool double_list,
		RaggedArray2DT<int>& neighbors);

protected:

	/** reference ID for sending output */
	int fOutputID;

	/** \name particle attributes */
	/*@{*/
	dArrayT fMass;
	iArrayT fType;
	/*@}*/
	
	/** \name cached calculated values */
	/*@{*/
	dArrayT   fEnergy;
	dArray2DT fForce;
	/*@{*/
	
	/** \name group running averages.
	 * Values are averages over {n1, n2,...,nN} steps */
	/*@{*/
	iArrayT fAverages;
	dArrayT fKE_avg;
	dArrayT fPE_avg;
	/*@}*/

	/** \name local to global tag map.
	 * Used for things like neighbor lists */
	/*@{*/
	/** list of node set IDs define particles in this group. If ALL
	 * nodes are particles, this will be length 1 with fID[0] = "ALL" */
	ArrayT<StringT> fID;
	
	/** map of global = fGlobalTag[local], which is compact in the local tags */
	iArrayT fGlobalTag;

	/** connectivities used to define the output set. Just an alias to the
	 * ParticleT::fGlobalTag. */
	iArray2DT fPointConnectivities;
	/*@}*/
	
	/** search grid */
	iGridManagerT* fGrid;
};

} /* namespace Tahoe */

#endif /* _PARTICLE_T_H_ */
