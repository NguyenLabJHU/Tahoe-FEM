/* $Id: FEManagerT_THK.h,v 1.14 2005-04-06 15:42:00 paklein Exp $ */
#ifndef _FE_MANAGER_THK_H_
#define _FE_MANAGER_THK_H_

/* element configuration header */
#include "ElementsConfig.h"
#include "DevelopmentElementsConfig.h"
#if defined(BRIDGING_ELEMENT) && defined(BRIDGING_ELEMENT_DEV)

/* base class */
#include "FEManagerT_bridging.h"
#include "dMatrixT.h"

namespace Tahoe {

/* forward declarations */
class ParticlePairT;

/** FEManagerT to support the time history kernel (THK) formulation */
class FEManagerT_THK: public FEManagerT_bridging
{
public:

	/** constructor */
	FEManagerT_THK(const StringT& input, ofstreamT& output, CommunicatorT& comm,
		const ArrayT<StringT>& argv, TaskT task);

	/** return array containing atom numbers of boundary and ghost atoms - 2D version **/
	const iArrayT& InterpolationNodes2D(void);
	
	/** return array containing atom numbers of boundary and ghost atoms - 3D veresion **/
	const iArrayT& InterpolationNodes3D(void);
	
	/** predictor routine for FEM solution interpolated to MD boundary atoms.
		predictor and corrector combined because of constant acceleration assumption.  **/
	void BAPredictAndCorrect(double timestep, dArray2DT& badisp, dArray2DT& bavel, dArray2DT& baacc);
	
	/** calculate THK force on boundary atoms for 2D disp/force formulation **/
	const dArray2DT& THKForce(const StringT& bridging_field, const dArray2DT& badisp);

	/** calculate THK disp for ghost atoms for 3D disp/disp formulation **/
	const dArray2DT& THKDisp(const StringT& bridging_field, const dArray2DT& badisp);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	// accessor for the ghostoffmapping
	const nMatrixT<int> GetGhostMap(void) { return fghostoffmap;};

private:

	/** 2D Bridging Scale Initialization */
	void Initialize2D(void);

	/** 3D Bridging Scale Initialization */
	void Initialize3D(void);

	/** compute theta tables for 2D disp/force formulation */
	void ComputeThetaTables2D(const StringT& data_file);
	                
	/** compute theta tables for 3D disp/disp formulation - ADD 3 more data_files as input */
	void ComputeThetaTables3D(const StringT& data_file);

private:

	/** \name input parameters */
	/*@{*/
	int fNcrit;
	double fLatticeParameter;
	StringT fThetaFile, fGhostMapFile;
	ArrayT<StringT> fTHKNodes;
	/*@}*/

	int fN_times, fNumstep_crit, fNeighbors;
	iArray2DT fTop, fBottom;
	dArray2DT fTHKforce, fGaussdisp, fTHKdisp, fInitdisp;
	iArray2DT fTop20, fTop21, fTop30, fTop31, fBottom20, fBottom21, fBottom30, fBottom31;
	iArrayT fTopatoms, fBottomatoms, fTopatoms2, fBottomatoms2, fToprow, fToprow2, fBottomrow, fBottomrow2, fShift;
	iArrayT fTA0, fTA1, fBA0, fBA1, fBottomrow0, fBottomrow1, fToprow0, fToprow1;

	/** interpolation points */
	iArrayT fInterpolationNodes;
	
	/** theta values: [neighbor] x [time x n_theta_values] */
	/** 2D Theta */
	ArrayT<dArray2DT> fThetaTable, fThetaTableT, fThetaTableB;
	
	/** 3D Thetas */
	ArrayT<dArray2DT> fTheta11t, fTheta12t, fTheta21t, fTheta22t;
	ArrayT<dArray2DT> fTheta11b, fTheta12b, fTheta21b, fTheta22b;
	
	/** displacement history: [n_boundary_atoms] x [time x ndof] */
	ArrayT<dArray2DT> fHistoryTablet, fHistoryTableb, fHistoryTablet1, fHistoryTableb1;
	
	// DEF added these
	// ghostoffmap matrix
	nMatrixT<int> fghostoffmap;
	iArrayT fthk_bound_lengths;
	iArray2DT fthk_boundary_atoms;
	int fmax_thk_bound_length;
	int fnumsets;
	/** displacement history: [boundary_n] x [[n_boundary_atoms] x [time x ndof]] */
	ArrayT< ArrayT<dArray2DT> > fHistoryTable;
	/** boundary neighbors: [boundary_n] x [[n_boundary_atoms] x [n_neighbors]] */
	ArrayT<iArray2DT> fbound_neighbor_atoms;
};

} /* namespace Tahoe */

#endif /* BRIDGING_ELEMENT && BRIDGING_ELEMENT_DEV */
#endif /* _FE_MANAGER_THK_H_ */
