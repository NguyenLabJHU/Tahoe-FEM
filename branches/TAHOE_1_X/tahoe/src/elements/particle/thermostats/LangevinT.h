/* $Id: LangevinT.h,v 1.4 2003-10-30 17:15:21 paklein Exp $ */
#ifndef _LANGEVIN_T_H_
#define _LANGEVIN_T_H_

#include "ios_fwd_decl.h"

/* base class */
#include "ThermostatBaseT.h"

/* direct members */
#include "iArrayT.h"
#include "RandomNumberT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** insert witty comment here */
class LangevinT: public ThermostatBaseT
{
public:
	
	/** constructor */
	LangevinT(ifstreamT& in, const int& nsd, const double& dt);
	LangevinT(void);

	/** destructor */
	virtual ~LangevinT(void) {};
	
	/** augment/overwrite forces with new ones */
	virtual void ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, AutoArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties);
	
	/** write properties to output */
	virtual void Write(ostream& out) const;
	
	/** write restart information */
	virtual void WriteRestart(ostream& out) const;
	
	/** read restart information */
	virtual void ReadRestart(istream& in);
	
	void SetRandNumGenerator(RandomNumberT* frand);
	
protected:

	/** \name properties */
	/*@{*/
	double fAmp;
	/*@}*/
	
	RandomNumberT* fRandom;
};

inline void LangevinT::SetRandNumGenerator(RandomNumberT* frand)
{
	fRandom = frand;
}


} /* namespace Tahoe */

#endif /* _LANGEVIN_T_H_ */
