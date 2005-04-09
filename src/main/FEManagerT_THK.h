/* $Id: FEManagerT_THK.h,v 1.15 2005-04-09 18:27:33 d-farrell2 Exp $ */

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
	const dArray2DT& THKForce2D(const StringT& bridging_field, const dArray2DT& badisp); // was THKForce

	/** calculate THK force for ghost atoms for 3D disp/disp formulation **/
	const dArray2DT& THKForce3D(const StringT& bridging_field, const dArray2DT& badisp); // was THKDisp

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
	/*@}*/

	// accessor for the ghostoffmapping
	const nMatrixT<int> GetGhostMap(void) { return fghostoffmap;};

private:

	/** 2D Bridging Scale Initialization */
	void Initialize2D(void);

	/** 3D Bridging Scale Initialization */
	void Initialize3D(void);

	/** compute theta tables for 2D disp/force formulation */
	void ComputeThetaTables2D(void);
	                
	/** compute theta tables for 3D disp/disp formulation */
	void ComputeThetaTables3D(void);

private:

	/** \name input parameters */
	/*@{*/
	int fNcrit;
	double fTcrit, fLatticeParameter;
	StringT fThetaFile, fGhostMapFile;
	ArrayT<StringT> fTHKNodes;
	/*@}*/

	int fN_times, fNumstep_crit, fNeighbors;
	dArray2DT fTHKforce;
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
	
	int ftotal_b_atoms; // total number of boundary atoms
	
	// array of THK BC plane normals
	ArrayT<dArrayT> fTHK_normals; 
	
	/** atoms in each node set: [boundary_n] x [n_boundary_atoms]  (has repeats) */
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
	
	// information for special atoms (corners and edges)
	StringT fSpecAtomFile;
	int fnum_spec_atoms; // number of special atoms
	ArrayT<iArrayT> fSpecAtomInfo; // special atom information array size: [# spec atoms] x [# sets atom is in + 2]
	iArrayT fSpecAtomID; // array of the atom number of the special atoms
};

} /* namespace Tahoe */

#endif /* BRIDGING_ELEMENT && BRIDGING_ELEMENT_DEV */
#endif /* _FE_MANAGER_THK_H_ */
