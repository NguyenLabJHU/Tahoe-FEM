/* $Id: From2Dto3DT.h,v 1.1 2003-05-28 23:13:07 cjkimme Exp $ */
/* created: paklein (06/23/1999) */

#ifndef _FROM_2D_TO_3D_T_H_
#define _FROM_2D_TO_3D_T_H_

/* base class */
#include "SurfacePotentialT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** Wrapper to convert a 2D cohesive model to one in 3D. It is assumed that
 *  the normal response of the model is indepedent of dimensionality while 
 *  the tangential response in 3D is isotropic in the tangent plane. This class
 *  is in principle extendable to 2D models it currently does not construct, but
 *  it is not designed to handle well those 2D models whose output depends on
 *  dimensionality, e.g. one whose output is a position-like vector that would
 *  need to be up-converted to 3D to have unambiguous meaning. */
class From2Dto3DT: public SurfacePotentialT
{
public:

	/** constructors */
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	From2Dto3DT(ifstreamT& in, int code, const double& time_step);
#endif
	/** constructor for use in SIERRA */
	From2Dto3DT(dArrayT& params);
	
	~From2Dto3DT(void);

	/** return the number of state variables needed by the model */
	int NumStateVariables(void) const;
	
	/** initialize the state variable array. By default, initialization
	 * involves only setting the array to zero. */
	virtual void InitStateVariables(ArrayT<double>& state);

	/** dissipated energy */
	virtual double FractureEnergy(const ArrayT<double>& state);

	/** incremental heat. The amount of energy per unit undeformed area
	 * released as heat over the current increment */
	virtual double IncrementalHeat(const dArrayT& jump, const ArrayT<double>& state);

	/** potential energy */
	virtual double Potential(const dArrayT& jump_u, const ArrayT<double>& state);
	
	/** surface traction. Internal variables are integrated over the current
	 * time step. */	
	virtual const dArrayT& Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate);

	/** tangent stiffness */
	virtual const dMatrixT& Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma);

	/** surface status */
	virtual StatusT Status(const dArrayT& jump_u, const ArrayT<double>& state);

	/** write model name to output */
	virtual void PrintName(ostream& out) const;

	/** write model parameters */
	virtual void Print(ostream& out) const;
	
	virtual int NumOutputVariables(void) const;
	
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	
	virtual void ComputeOutput(const dArrayT& jump, const ArrayT<double>& state, 
		dArrayT& output);
	
private:

	SurfacePotentialT* f2DModel;

};

} // namespace Tahoe 
#endif /* _XU_NEEDLE_3D_T_H_ */
