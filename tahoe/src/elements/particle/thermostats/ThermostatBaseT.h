/* $Id: ThermostatBaseT.h,v 1.1 2003-04-16 18:15:54 cjkimme Exp $ */
#ifndef _THERMOSTAT_BASE_T_H_
#define _THERMOSTAT_BASE_T_H_

#include "ios_fwd_decl.h"

/* direct members */
#include "iArrayT.h"
#include "dArrayT.h"
#include "RandomNumberT.h"
#include "RaggedArray2DT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class dArray2DT;

/** base class for thermostatting and damping */
class ThermostatBaseT
{
public:

	enum ThermostatT {
					kFree = 0, /**< you figure it out */
    	  	  	  kDamped = 1, /**< velocity-dependent damping */
    			kLangevin = 2, /**< Langevin (stochastic) thermostat */
	  	  	  kNoseHoover = 3,
	  kGaussianIsokinetic = 4
	};
   
	enum ControlTypeT {
   					kNodes = 0, /**< apply to some or all nodesets */
    			   kRegion = 1 /**< apply to nodes in region of space */
	}; 
	
	/** stream extraction operators */
	friend istream& operator>>(istream& in, ThermostatBaseT::ThermostatT& property);	

	/** constructor */
	ThermostatBaseT(ifstreamT& in, int nsd, double dt);

	/** destructor */
	virtual ~ThermostatBaseT(void) {};
	
	/** damping coefficient */
	double Beta(void) const { return fBeta; };
	
	/** temperature */
	virtual double Temperature(void) const { return fTemperature; }
	
	/** nodes that are thermostatted */
	iArrayT& NodeList(void); 
	
	/** write properties to output */
	virtual void Write(ostream& out) const;
	
	/** write restart information */
	virtual void WriteRestart(ostream& out) const;
	
	/** read restart information */
	virtual void ReadRestart(istream& in);
	
	/** augment/overwrite forces with new ones */
	virtual void ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
					dArray2DT& forces);
	
	void NodesInRegion(const dArray2DT& coords,	
					const ArrayT<int>* partition_nodes);

	
protected:

	/** \name properties */
	/*@{*/
	double fBeta;	
	double fTemperature;
	double fTimeStep;
	/*@}*/
	
	/** Nodes that are thermostatted */
	iArrayT fNodes;
	
	/** True if stochastic force is returned */
//	bool QLangevin;
	
	/** Number of spatial dimensions */
	int fSD;
	
//	RandomNumberT* fRandom;
	
	/** \name Region paramters */
	/*@{*/
	/** Bounding box */
	dArrayT fxmin, fxmax;
	/** Re-check for particles in region every nIncs timesteps */
	int nIncs; 
	/*@}*/
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
