/* $Id: CSEAnisoT.h,v 1.35 2004-07-15 08:25:57 paklein Exp $ */
/* created: paklein (11/19/1997) */
#ifndef _CSE_ANISO_T_H_
#define _CSE_ANISO_T_H_

/* base class */
#include "CSEBaseT.h"

/* direct members */
#include "pArrayT.h"
#include "RaggedArray2DT.h"
#include "dArray2DT.h"
#include "Array2DT.h"
#include "LocalArrayT.h"

namespace Tahoe {

/* forward declarations */
class SurfacePotentialT;
#ifndef _FRACTURE_INTERFACE_LIBRARY_
class TiedPotentialBaseT;
#endif

/** Cohesive surface elements with vector argument cohesive relations. */
class CSEAnisoT: public CSEBaseT
{
public:

	/** constructors */
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	CSEAnisoT(const ElementSupportT& support);
#else
	CSEAnisoT(ElementSupportT& support, bool rotate);
#endif
	
	/* destructor */
	~CSEAnisoT(void);

	/* form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/** prepare for a sequence of time steps */
	virtual void InitialCondition(void);

	/** close current time increment */
	virtual void CloseStep(void);

#ifndef _FRACTURE_INTERFACE_LIBRARY_
	/** write restart data to the output stream. */
	virtual void WriteRestart(ostream& out) const;

	/** read restart data to the output stream. */
	virtual void ReadRestart(istream& in);

	/** state variable array */
	RaggedArray2DT<double>& StateVariables(void) { return fStateVariables; };
#else

  	/* send restart array */
	virtual void WriteRestart(double* outgoingData) const;
	
	/* receive restart array */
	virtual void ReadRestart(double* incomingData);
#endif	

#ifdef _FRACTURE_INTERFACE_LIBRARY_	
	/* Initialize fields passed in from the outside */
	virtual void InitStep(void);
#endif

	/** compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);

	/** set the active elements.
	 * \param array of status flags for all elements in the group */
	virtual void SetStatus(const ArrayT<StatusT>& status);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* tangent matrix and force vector */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);
	virtual void RHSDriver(void);

	/* nodal value calculations */
	virtual void SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;
	virtual void SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;

	/* compute output values */
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
		const iArrayT& e_codes, dArray2DT& e_values);

	/* construct output labels array */
	virtual void GenerateOutputLabels(const iArrayT& n_codes, ArrayT<StringT>& n_labels,
		const iArrayT& e_codes, ArrayT<StringT>& e_labels) const;

	/* write all current element information to the stream */
	virtual void CurrElementInfo(ostream& out) const;

	/* operations with pseudo rank 3 (list in j) matrices */
	void u_i__Q_ijk(const dArrayT& u, const ArrayT<dMatrixT>& Q,
		dMatrixT& Qu);

	void Q_ijk__u_j(const ArrayT<dMatrixT>& Q, const dArrayT& u,
		dMatrixT& Qu);
		
	/* Fake output to send to TiedNodesT */
	void ComputeFreeNodesForOutput(void);

	/* get output from surrounding continuum */
	void StoreBulkOutput(void);

	/* Compute surface values from bulk values */	
	void SurfaceValuesFromBulk(const ElementCardT& element, iArrayT& ndIndices,
			dArray2DT& elementVals, LocalArrayT& nodal_values);
		
	/* Interpolate bulk quantities to integration points, rotating to local frame if desired */	
	void FromNodesToIPs(bool rotate, dArrayT& localFrameIP, LocalArrayT& nodal_values);

	/* Query tied cohesive model to see if nodes' status will change */	
	void UntieOrRetieNodes(int elNum, int nnd, const TiedPotentialBaseT* tiedpot, 
							ArrayT<double>& state, dArrayT& localFrameIP);

	/* initialize thermal source terms*/
	void InitializeTemperature(const FieldT* temperature);

protected:

	/* shape function - if (fRotate) then wrt current config */
	bool fRotate;
	SurfaceShapeT* fCurrShapes;

	/* cohesive surface potentials */
	iArrayT fNumStateVariables;
	pArrayT<SurfacePotentialT*> fSurfPots;
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	ArrayT<TiedPotentialBaseT*> fTiedPots;
//	TiedPotentialBaseT* tiedpot;
	bool qRetieNodes;

	/** \name state variable storage arrays. 
	 * arrays have dimensions: [nel] x [nip * nvar] */
	/*@{*/
	RaggedArray2DT<double> fStateVariables;
	RaggedArray2DT<double> fStateVariables_last;
	/*@}*/

	const GlobalT::StateT& fRunState;

#else
	/*In SIERRA, we're assuming only one cohesive law per element
	 *block, but keeping the array 2D so we can loop over its
	 *minorDim of 1. Most importantly, SIERRA controls and allocates
	 *state variables, so these should just point to the data they
	 *give us. 
	 */
	dArray2DT fStateVariables;
	dArray2DT fStateVariables_last;
#endif

	/** incremental heat sources for each element block */
	ArrayT<dArray2DT> fIncrementalHeat;
	
	/* coordinate transformation */
	dMatrixT fQ;     // t'_i = Q_ji t_j, where t' is in the local frame
	dArrayT  fdelta; // gap vector (local frame)
	dArrayT  fT;     // traction vector (global frame)
	dMatrixT fddU;	 // surface stiffness (local frame)
	
	ArrayT<dMatrixT> fdQ; // list representation of rank 3 of dQ_ij/du_k

	double fIPArea; /**< reference area associated with the current integration point */
	
	/* work space (for tangent) */
	dMatrixT fnsd_nee_1;
	dMatrixT fnsd_nee_2;

	/* variables for calculating nodal info */
	bool fCalcNodalInfo;
	int fNodalInfoCode, iTiedFlagIndex;
	dArray2DT fNodalQuantities;
	iArrayT iBulkGroups, otherInds;
	
	/* if nodes are tied, keep track of free nodes per element */
	dArray2DT freeNodeQ, freeNodeQ_last;
};

} // namespace Tahoe 
#endif /* _CSE_ANISO_T_H_ */
