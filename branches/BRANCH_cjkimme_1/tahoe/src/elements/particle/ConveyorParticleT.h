/* $Id: ConveyorParticleT.h,v 1.1.2.1 2003-09-18 21:03:36 cjkimme Exp $ */
#ifndef _CONVEYOR_PARTICLE_T_H_
#define _CONVEYOR_PARTICLE_T_H_

/* base class */
#include "ParticlePairT.h"

/* direct members */


namespace Tahoe {


/** class to interact with AtomicConveyorT KBC_Controller for steady-state crack
 *  propagation. This class primarily exists to modify interactions above
 *  and below the pre-crack plane as well as to handle the ramped damping left
 *  and right BCs as they change when material is cut and pasted from left to
 *  right. */
class ConveyorParticleT: public ParticlePairT
{
public:

	/** constructor */
	ConveyorParticleT(const ElementSupportT& support, const FieldT& field);

	/** destructor */
	~ConveyorParticleT(void);
	
	/** initialization */
	virtual void Initialize(void);

	/** \name restart functions */
	/*@{*/
	/** write restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::ReadRestart implementation. */
	virtual void WriteRestart(ostream& out) const;

	/** read restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::WriteRestart implementation. */
	virtual void ReadRestart(istream& in);
	/*@}*/
	
	/** turn off interactions between atoms above and below crack plane */
	void CreateNoninteractingAtoms(iArrayT& topAtoms, iArrayT& bottomAtoms);

protected: /* for derived classes only */

	

protected:

	/** \name Backups of ParticleT types */
	/*@{*/
	/** number of types. Read during ParticleT::EchoConnectivityData. */
	int fNumTypes_Actual;

	/** particle type for global tag -- fType will hold types that may have been
	 *  modified by the presence of the precrack. */
	iArrayT fType_Actual;
	
	/** map of particle types to properties: {type_a, type_b} -> property number */
	nMatrixT<int> fPropertiesMap_Actual;
	/*@}*/
	
	/** particle properties list */
	ArrayT<ParticlePropertyT*> fParticleProperties;
	
};

} /* namespace Tahoe */

#endif /* _CONVEYOR_PARTICLE_T_H_ */
