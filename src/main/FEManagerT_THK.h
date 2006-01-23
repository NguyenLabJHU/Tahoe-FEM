/* $Id: FEManagerT_THK.h,v 1.19 2006-01-23 23:17:04 d-farrell2 Exp $ */

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

	/** return array containing atom numbers of boundary and ghost atoms **/
	const iArrayT& InterpolationNodes(void);
	
	/** predictor routine for FEM solution interpolated to MD boundary atoms.
		predictor and corrector combined because of constant acceleration assumption.  **/
	void BAPredictAndCorrect(double timestep, dArray2DT& badisp, dArray2DT& bavel, dArray2DT& baacc);
	
	/** calculate THK displacement for ghost atoms for 2/3D disp formulation **/
	const dArray2DT& THKDisp(const StringT& bridging_field, const dArray2DT& badisp);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;
	
	/** set up new subordinate parameter list */
	ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	
	/** 2D/3D MD/THK and BSM THK Initialization */
	void InitializeTHK(bool ignore_continuum);
	
	/** accessor for the ghostoffmapping */
	const nMatrixT<int> GetGhostMap(void) { return fghostoffmap;};

	/*@}*/

	

private:
	
	/** perform neighbor search for THK boundary atoms, 2D */
	void DoNeighSearch2D(void);
	
	/** perform neighbor search for THK boundary atoms, 3D */
	void DoNeighSearch3D(void);
	
	/** find the ghost atom properties map */
	void DoGhostMap(void);
	
	/** compute theta tables for 2D/3D disp/disp or disp/force formulation (doesn't matter, its all the same) */
	void ComputeThetaTables(void);

private:

	/** \name input parameters */
	/*@{*/
	int fNcrit;
	double fTcut, fLatticeParameter, fSearchParameter;
	StringT fThetaFile, fGhostMapFile;
	ArrayT<StringT> fTHKNodes, fTHKGhostNodes;
	/*@}*/

	int fN_times, fNumstep_crit, fNeighbors;
	dArray2DT fTHKforce, fTHKdisp;
	iArrayT fShift;

	/** interpolation points */
	iArrayT fInterpolationNodes;
	
// DEF added these:
	
	// ghostoffmap matrix
	nMatrixT<int> fghostoffmap;
	iArrayT fthk_bound_lengths;
	iArray2DT fthk_boundary_atoms;
	int fmax_thk_bound_length;
	int fnumsets;
	
	int ftotal_b_atoms, ftotal_g_atoms; // total number of boundary atoms, ghost atoms
	
	// array of THK BC plane normals
	ArrayT<dArrayT> fTHK_normals; 
	
	/** atoms in each node set (ghost atoms): [boundary_n] x [n_boundary_atoms]  (no repeats) */
	ArrayT<iArrayT> fghost_set_atoms;
	
	/** atoms in each node set (real atoms): [boundary_n] x [n_boundary_atoms]  (has repeats) */
	ArrayT<iArrayT> fbound_set_atoms;
	
	/** displacement history: [boundary_n] x [[n_boundary_atoms] x [time x ndof]] (has repeats) */
	ArrayT< ArrayT<dArray2DT> > fHistoryTable;
	
	/** THK values: [boundary_n] x [[neighbor] x [time x n_theta_values]] (has repeats)*/
	ArrayT< ArrayT<dArray2DT> > fThetaTable_array;
	
	/** boundary neighbors: [boundary_n] x [[n_boundary_atoms] x [n_neighbors]] (has repeats)*/
	ArrayT<iArray2DT> fbound_neighbor_atoms;
	
	/** file containing the fourier coefficients for the sine series used to calculate Theta matrices */
	// size: [n_sets] - 1 entry per set
	ArrayT<StringT> fThetaFile_array;
	
	// parameter which allows the scaling of Tcrit and the Bn's for the sine series (defaults to 1)
	// This allows fourier sine coeffs calculated for k,m = 1 to be used for any k, m it is sqrt(k/m) for the desired values
	double fOmega_sys;
	
};

} /* namespace Tahoe */

#endif /* BRIDGING_ELEMENT && BRIDGING_ELEMENT_DEV */
#endif /* _FE_MANAGER_THK_H_ */
