/* $Id: ParaDynPrescribedT.h,v 1.1.2.1 2003-09-18 21:03:44 cjkimme Exp $ */
#ifndef _PARADYN_PRESCRIBED_T_H_
#define _PARADYN_PRESCRIBED_T_H_

/* base class */
#include "KBC_ControllerT.h"

/* direct members */
#include "ScheduleT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "BasicFieldT.h"

namespace Tahoe {

/* forward declarations */
class BasicFieldT;
class RandomNumberT;

/** Prescribes displacements and velocities read from ParaDyn format files.
 */
class ParaDynPrescribedT: public KBC_ControllerT
{
public:	

	/** constructor */
	ParaDynPrescribedT(NodeManagerT& node_manager, BasicFieldT& field);

	/** destructor */
 	~ParaDynPrescribedT(void);

	/** initialize data. Must be called immediately after construction */
	virtual void Initialize(ifstreamT& in);

	/** write class parameters */
	void WriteParameters(ostream& out) const;

	/** do at start of timestep */
	virtual void InitStep(void);

	/** Initialize to appropriate temperature */
	virtual void InitialCondition(void);
	
	virtual bool IsICController(void) { return false; }
	
	static dArray2DT* SendCoordinates(void);
	
	static dArray2DT* SendVelocities(void);

protected:
	
	void SetBCCards(void);
	/*@}*/

	/** Read data from a Paradyn dump file */	
	void ReadArray(StringT& file_name, dArray2DT& info, int numGarbageLines);
	
	/** Construct an input file name */
	void MakeFileName(StringT& fileName, StringT& rootName, StringT& suffixName, int index);

protected:

	/** the field */
	BasicFieldT& fField;
	
	/** Throw a bone to the KBC Cards */
	ScheduleT fDummySchedule;
	
	/** thermostatted nodes */
	iArrayT fNodes;

	/** files */
	StringT disp_file_root, disp_file_suffix, vel_file_root, vel_file_suffix;
	
	/** ranges for file names */
	int n_start, n_end, n_index;
	
	/** prescribed DOFs */
	dArray2DT coords, vels;

	static dArray2DT* coord_ptr;
	static dArray2DT* vel_ptr;
};

} // namespace Tahoe 
#endif /*  _PARADYN_PRESCRIBED_T_H_ */
