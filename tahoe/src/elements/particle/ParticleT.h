/* $Id: ParticleT.h,v 1.8 2002-12-04 06:32:40 paklein Exp $ */
#ifndef _PARTICLE_T_H_
#define _PARTICLE_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "nVariArray2DT.h"
#include "InverseMapT.h"

namespace Tahoe {

/** forward declarations */
class iGridManagerT;

/** base class for particle types */
class ParticleT: public ElementBaseT
{
public:

	/** enum for particle property types */
	enum PropertyT {
        kHarmonicPair = 0, /**< harmonic pair potential */
    kLennardJonesPair = 1, /**< Jennard-Jones 6/12 pair potential */
         kParadynPair = 2  /**< pair potential in Paradyn (EAM) format */
	};
	
	/** stream extraction operator */
	friend istream& operator>>(istream& in, ParticleT::PropertyT& property);	

	/** constructor */
	ParticleT(const ElementSupportT& support, const FieldT& field);

	/** destructor */
	~ParticleT(void);
	
	/** initialization */
	virtual void Initialize(void);

	/** form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/** NOT implemented. Returns an zero force vector */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);
			
	/** returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void) { return 0.0; };
	
	/** resgiter for writing output. Uses the uses sub-class implementations
	 * of ParticleT::GenerateOutputLabels to register the particle group for
	 * output. Sub-classes also need to implemented the WriteOutput method. */
	virtual void RegisterOutput(void);

	/** write output. ParticleT::WriteOutput only writes search grid statistics.
	 * Sub-classes are responsible for writing data for each particle, given the
	 * variables names returned by ParticleT::GenerateOutputLabels. */
	virtual void WriteOutput(void);

	/* compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);

	/** trigger reconfiguration */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

	/** \name restart functions */
	/*@{*/
	/** write restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::ReadRestart implementation. */
	virtual void WriteRestart(ostream& out) const;

	/** read restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::WriteRestart implementation. */
	virtual void ReadRestart(istream& in);
	/*@}*/
	 			  	
protected: /* for derived classes only */

	/** echo element connectivity data. Reads parameters that define
	 * which nodes belong to this ParticleT group. */
	virtual void EchoConnectivityData(ifstreamT& in, ostream& out);
	
	/** generate labels for output data */
	virtual void GenerateOutputLabels(ArrayT<StringT>& labels) const = 0;

	/** return true if connectivities are changing */
	virtual bool ChangingGeometry(void) const;

	/** set neighborlists and any other system configuration information
	 * based on the current information. */
	virtual void SetConfiguration(void) = 0;

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

	/** construct the list of properties from the given input stream */
	virtual void EchoProperties(ifstreamT& in, ofstreamT& out) = 0;

	/** assemble particle mass matrix into LHS of global equation system
	 * \param mass mass associated with each particle type */
	void AssembleParticleMass(const dArrayT& mass);

protected:

	/** reference ID for sending output */
	int fOutputID;

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

	/** number of steps between reseting neighbor lists */
	int fReNeighborIncr;

	/** \name map from global tag to local tag */
	InverseMapT fGlobalToLocal;

	/** \name particle properties */
	/*@{*/
	/** number of types. Read during ParticleT::EchoConnectivityData. */
	int fNumTypes;

	/** particle type for global tag */
	iArrayT fType;
	
	/** map of particle types to properties: {type_a, type_b} -> property number */
	nMatrixT<int> fPropertiesMap;
	/*@}*/

	/** \name cached calculated values */
	/*@{*/
//	dArrayT   fEnergy;
	dArray2DT fForce;
	/*@{*/
	
	/** \name group running averages.
	 * Values are averages over {n1, n2,...,nN} steps */
	/*@{*/
//	iArrayT fAverages;
//	dArrayT fKE_avg;
//	dArrayT fPE_avg;
	/*@}*/
	
	/** the search grid */
	iGridManagerT* fGrid;

private:

	/** count between resetting neighbor lists */
	int fReNeighborCounter;
};

} /* namespace Tahoe */

#endif /* _PARTICLE_T_H_ */
