/* $Id: ParticleT.h,v 1.34.6.2 2005-05-31 06:14:36 paklein Exp $ */
#ifndef _PARTICLE_T_H_
#define _PARTICLE_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "VariArrayT.h"
#include "nVariArray2DT.h"
#include "InverseMapT.h"
#include "dArray2DT.h"

namespace Tahoe {

/** forward declarations */
class iGridManagerT;
class CommManagerT;
class ParticlePropertyT;
class dSPMatrixT; //TEMP
class InverseMapT;
class RandomNumberT;
class ThermostatBaseT;

/** base class for particle types */
class ParticleT: public ElementBaseT
{
public:

	/** constructor */
	ParticleT(const ElementSupportT& support);

	/** destructor */
	~ParticleT(void);

	/** (re-)set the system configuration */
	virtual GlobalT::InitStatusT UpdateConfiguration(void);

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

	/** \name restart functions */
	/*@{*/
	/** write restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::ReadRestart implementation. */
	virtual void WriteRestart(ofstreamT& out) const;

	/** read restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::WriteRestart implementation. */
	virtual void ReadRestart(ifstreamT& in);
	/*@}*/

	/** define the particles to skip. This is a list of nodes though "owned" 
	 * by this processor and appearing in the list of particles, should be skipped 
	 * in the calculation of energy, force, etc. This method must be called
	 * whenever the list changes as it is used to reset some internal data. An
	 * empty list should be passed to clear any previous lists. The list may
	 * contain duplicates. Note this method does not trigger recalculation of
	 * the neighborlists. This can be triggered explicitly with a call to
	 * ParticleT::SetConfiguration */
	void SetSkipParticles(const iArrayT& skip);

	/** set neighborlists and any other system configuration information
	 * based on the current information. ParticleT::SetConfiguration
	 * simply stores a list of coordinates in the current configuration
	 * which is used to determine if SetConfiguration needs to be called
	 * again. */
	virtual void SetConfiguration(void);

	/** compute the part of the stiffness matrix */
	virtual void FormStiffness(const InverseMapT& col_to_col_eq_row_map,
		const iArray2DT& col_eq, dSPMatrixT& stiffness) = 0;

	/** contribution to the nodal residual forces. Return the contribution of this element
	 * group to the residual for the given solver group. ParticleT::InternalForce 
	 * returns the internal force calculated with the latest call to ElementBaseT::FormRHS. */
	virtual const dArray2DT& InternalForce(int group);

	/** read/write access to the properties map */
	nMatrixT<int>& PropertiesMap(void) { return fPropertiesMap; };

	/** add the element group's contribution to the lumped (scalar) mass of the given nodes */
	virtual void LumpedMass(const iArrayT& nodes, dArrayT& mass) const;

	/** modify the thickness of the image atom layer. The image atom layer is used to enforce
	 * periodic boundary conditions. Values less that ParticleT::fNeighborDistance are not allowed */
	void SetPeriodicSkin(double skin) { fPeriodicSkin = (skin > fNeighborDistance) ? skin : fNeighborDistance; };

#ifdef __NO_RTTI__
	/** \name fixes for environments without working RTTI */
	/*@{*/
	/** cast this to ParticleT */
	virtual ParticleT* dynamic_cast_ParticleT(void) { return this; };
	/*@}*/
#endif

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list. Properties information is extracted from the
	 * list with the call to ParticleT::ExtractProperties, a purely virtual
	 * method that should be overridden by derived types specifying their own
	 * properties information. */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** atomic kinetic energies */
	void AtomicKineticEnergies(dArrayT& ke);

protected: /* for derived classes only */

	/** \name initialization methods */
	/*@{*/
	/** extract the properties information from the parameter list
	 * \name list source for parameters 
	 * \name type_names list of type names
	 * \name properties passed in empty and should return with a list of pointers to properties
	 *       information. The pointers should be dynamically allocated and will be freed by
	 *       ParticleT on destruction.
	 * \name properties_map map dimension (number of types) x (number of types) w */
	virtual void ExtractProperties(const ParameterListT& list, const ArrayT<StringT>& type_names,
		ArrayT<ParticlePropertyT*>& properties, nMatrixT<int>& properties_map) = 0;
	/*@}*/

	/** generate labels for output data */
	virtual void GenerateOutputLabels(ArrayT<StringT>& labels) const = 0;

	/** return true if connectivities are changing */
	virtual bool ChangingGeometry(void) const;

	/** generate neighborlist
	 * \param particle_tags global tags for which to determine neighhors. If NULL, find
	 *        neighbors for all nodes.
	 * \param distance distance over which to search for neighbors
	 * \param neighbors list of neighbors for every tag in particle_tags as rows in the
	 *        2D array. The list for each tag begins with the tag itself. Therefore, all
	 *        lists are at least length 1. The numbering of the tags in the neighbor list
	 *        is in terms of the tags provided in particle_tags.
	 * \param double_list if true the neighbor lists will contain two references for
	 *        every neighbor interaction. For tags A and B, this means B will appear
	 *        in the neighbor list for A and vice versa. If false, neighbor lists will
	 *        only contain neighbors for A satisfying B > A.
	 * \param full_list if double_list is false, passing full_list true will add B to
	 *        the neighbor list of A regardless of which tag is higher in order to
	 *        produce a full list of interactions for A. If double_list is true,
	 *        this flag has no effect. The list generated with full_list true differs
	 *        from the one generated by double_list under two circumstances:
	 *     -# if B < A, the A-B interaction would normally appear in the neighbor list of B. Tags
	 *        marked to be skipped, using ParticleT::SetSkipParticles, do not have neighbor lists,
	 *        so if B is marked to be skipped, the A-B interaction must be included in the
	 *        neighborlist of A.
	 *     -# when there is no assembly of partial forces across processor boundaries, all
	 *        bonds across processor boundaries need to be reproduced on both processors.
	 */
	void GenerateNeighborList(const ArrayT<int>* particle_tags, double distance, 
		RaggedArray2DT<int>& neighbors, bool double_list, bool full_list);

	/** assemble particle mass matrix into LHS of global equation system
	 * \param mass mass associated with each particle type */
	void AssembleParticleMass(const dArrayT& mass);

	/** return the maximum distance between the current coordinates and those
	 * stored in ParticleT::fReNeighborCoords. The maximum is over local atoms
	 * only. */
	double MaxDisplacement(void) const;
	
	/** Apply thermostatting/damping constraints to forces */
	void ApplyDamping(const RaggedArray2DT<int>& fNeighbors);

	/** Read in damping and thermostatting paramters. Called from
	 *  Initialize */
	virtual void SetDamping(const ParameterListT& list);

	/** return a new pair property or NULL if the name is invalid */
	ThermostatBaseT* New_Thermostat(const StringT& name, bool throw_on_fail) const;

protected:

	/** reference ID for sending output */
	int fOutputID;

	/** \name local to global tag map.
	 * Used for things like neighbor lists */
	/*@{*/
	/** list of node set IDs define particles in this group. If ALL
	 * nodes are particles, this will be length 1 with fID[0] = "ALL" */
	ArrayT<StringT> fID;
	
	/** connectivities used to define the output set. Just an alias to the
	 * ParticleT::fGlobalTag. */
	iArray2DT fPointConnectivities;
	/*@}*/

	/** the neighboring cut-off distance */
	double fNeighborDistance;

	/** distance to nearest neighbors */
	double fLatticeParameter;

	/** characteristic distance to nearest neighbor */
	double fNearestNeighborDistance;	

	/** maximum distance an atom can displace before the neighbor lists are reset */
	double fReNeighborDisp;

	/** maximum number of steps between reseting neighbor lists */
	int fReNeighborIncr;

	/** \name periodic boundary conditions */
	/*@{*/
	/** thickess of the image atoms layer used to enforce periodic boundary conditions. This
	 * initialized to ParticleT::fNeighborDistance and can be modified with ParticleT::SetPeriodicSkin */
	double fPeriodicSkin;

	/** periodic lower and upper bounds for each coordinate direction */
	dArray2DT fPeriodicBounds;
	
	/** schedule for stretching of periodic boundaries. The schedules should
	 * have a value of 1 at time = 0.0 and increase or decrease from there
	 * to define the stretch in the periodic BC's as a function of time. */
	ArrayT<const ScheduleT*> fStretchSchedule;
	double fStretchTime;
	/*@}*/

	/** \name particle properties */
	/*@{*/
	/** number of types. Read during ParticleT::EchoConnectivityData. */
//	int fNumTypes;

	/** names of each particle type */
	ArrayT<StringT> fTypeNames;

	/** particle type for global tag */
	AutoArrayT<int> fType;
	
	/** message id for exchanging the ParticleT::fType array */
	int fTypeMessageID;
	
	/** map of particle types to properties: {type_a, type_b} -> property number */
	nMatrixT<int> fPropertiesMap;

	/** particle properties list */
	ArrayT<ParticlePropertyT*> fParticleProperties;
	
	/** list of active particles. The pointer will be NULL if all nodes in the
	 * CommManagerT::PartitionNodes list are active. */
	AutoArrayT<int>* fActiveParticles;
	/*@}*/

	/** \name cached calculated values */
	/*@{*/
//	dArrayT   fEnergy;
	dArray2DT fForce;
	nVariArray2DT<double> fForce_man;
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

	/** coordinates of partition nodes at the last re-neighbor. The list of coordinares
	 * is collected during ParticleT::SetConfiguration */
	dArray2DT fReNeighborCoords;
	double fDmax;  /**< maximum distance between the current
	                    coordinates and the coordinates in ParticleT::fReNeighborCoords.
	                    This value is computed during ParticleT::RelaxSystem. */

	/* Damping, thermostatting variables */
//	bool QisDamped;
//	RandomNumberT* fRandom;
	ArrayT<ThermostatBaseT*> fThermostats;

	/** \name workspace for ParticlePairT::RHSDriver. Used to accumulate the force for
	 * a single row of ParticlePairT::fNeighbors. */
	/*@{*/
	dArrayT fForce_list;
	VariArrayT<double> fForce_list_man;

	/** constant matrix needed to compute the stiffness */
	dMatrixT fOneOne;
	/*@}*/
	
	/* calculate slip vector and strain tensor */
	void Calc_Slip_and_Strain(dArray2DT &s_values, RaggedArray2DT<int> &RefNearestNeighbors, const int &kEulerLagr) ;

	int Combination(int n, int k);

	/* compute the centrosymmetry parameter */
	void Calc_CSP(const RaggedArray2DT<int>& NearestNeighbors, dArrayT& csp);
	
	// compute the coordination number
	void Calc_CN(const RaggedArray2DT<int>& NearestNeighbors, iArrayT& cnarray);
    
	void SetRefNN(RaggedArray2DT<int> &NearestNeighbors,RaggedArray2DT<int> &RefNearestNeighbors);

private:

	/** count between resetting neighbor lists */
	int fReNeighborCounter;
	bool fhas_periodic;
	bool fhas_moving_periodic;
	dArrayT fPeriodicLengths;
};

} /* namespace Tahoe */

#endif /* _PARTICLE_T_H_ */
