/* $Id: ThermostatBaseT.h,v 1.7.20.1 2005-06-05 06:23:55 paklein Exp $ */
#ifndef _THERMOSTAT_BASE_T_H_
#define _THERMOSTAT_BASE_T_H_

/* base class */
#include "ParameterInterfaceT.h"

#include "ios_fwd_decl.h"

/* direct members */
#include "iArrayT.h"
#include "dArrayT.h"
#include "RandomNumberT.h"
#include "RaggedArray2DT.h"
#include "AutoArrayT.h"
#include "ScheduleT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class dArray2DT;
class ParticlePropertyT;
class BasicSupportT;

/** base class for thermostatting and damping */
class ThermostatBaseT: public ParameterInterfaceT
{
public:

	enum ThermostatT {
					kFree = 0, /**< you figure it out */
    	  	  	  kDamped = 1, /**< velocity-dependent damping */
    			kLangevin = 2, /**< Langevin (stochastic) thermostat */
	  	  	  kNoseHoover = 3, /**< Nose-Hoover determistic thermostat */
	  	 kGaussIsokinetic = 4, /**< Gaussian Isokinetic deterministic thermostat */
	  	   kRampedDamping = 5  /**< Damping coefficient varies over space  */
	};
   
	enum ControlTypeT {
   					kNodes = 0, /**< apply to some or all nodesets */
    			   kRegion = 1 /**< apply to nodes in region of space */
	}; 

	/** constructor */
	ThermostatBaseT(const BasicSupportT& support);

	/** destructor */
	virtual ~ThermostatBaseT(void) {};
	
	/** damping coefficient */
	double Beta(void) const { return fBeta; };
	
	/** temperature */
	virtual double Temperature(void) const { return fTemperature; }
	
	/** nodes that are thermostatted */
	iArrayT& NodeList(void); 
	
	/** write restart information */
	virtual void WriteRestart(ostream& out) const;
	
	/** read restart information */
	virtual void ReadRestart(istream& in);
	
	/** augment/overwrite forces with new ones */
	virtual void ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, ArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** \name particle picking methods */
	/*@{*/
	/** initialize nodes to thermostat by NodeSet */		
	virtual void InitNodeSets(const ParameterListT& pick_nodes);
	
	/** initialize nodes to thermostat by spatial coords */
	virtual void InitRegion(const ParameterListT& pick_region);
	/*@}*/

	/** generate a node set based on spatial region */
	void NodesInRegion(const dArray2DT& coords, const ArrayT<int>* partition_nodes);
	
protected:

	/** \name host code support */
	const BasicSupportT& fSupport;

	/** \name properties */
	/*@{*/
	double fBeta;	
	double fTemperature;
//	double fTimeStep;
	/*@}*/
	
	/** \name nodes that are thermostatted */
	/*@{*/
	bool fAllNodes;
	iArrayT fNodes;
	/*@}*/
	
	/** Number of spatial dimensions */
//	int fSD;
		
	/** \name Region paramters */
	/*@{*/
	/** Bounding box */
	dArrayT fxmin, fxmax;
	/** Re-check for particles in region every nIncs timesteps */
	int nIncs; 
	/*@}*/
	
	const ScheduleT* fTemperatureSchedule;
	double fTemperatureScale;
};

inline iArrayT& ThermostatBaseT::NodeList(void)
{
	return fNodes;
}

//inline void ThermostatBaseT::SetRandNumGenerator(RandomNumberT* frand)
//{
//	fRandom = frand;
//}

} /* namespace Tahoe */

#endif /* _THERMOSTAT_BASE_T_H_ */
