/* $Id: FEManagerT_THK.h,v 1.5 2003-07-11 16:45:19 hspark Exp $ */
#ifndef _FE_MANAGER_THK_H_
#define _FE_MANAGER_THK_H_

/* element configuration header */
#include "ElementsConfig.h"
#ifdef BRIDGING_ELEMENT

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
	FEManagerT_THK(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
		ifstreamT& bridging_input);

	/** initialize members */
	virtual void Initialize(InitCodeT init = kFull);

	/** \name solution steps. 
	 * See FEManagerT for more information */
	/*@{*/
	/** initialize the current time increment for all groups */
	virtual ExceptionT::CodeT InitStep(void);

	/** close the current time increment for all groups */
	virtual ExceptionT::CodeT CloseStep(void);
	/*@}*/
	
	/** return array containing atom numbers of boundary and ghost atoms **/
	const iArrayT& InterpolationNodes(void);
	
	/** predictor routine for FEM solution interpolated to MD boundary atoms.
		predictor and corrector combined because of constant acceleration assumption.  **/
	void BAPredictAndCorrect(double timestep, dArray2DT& badisp, dArray2DT& bavel, dArray2DT& baacc);
	
	/** calculate THK force on boundary atoms **/
	const dArray2DT& THKForce(const dArray2DT& badisp);

	/** impose gaussian wave initial displacements */
	const dArray2DT& GaussianWave(void);

private:

	/** compute theta tables */
	void ComputeThetaTables(const StringT& data_file);
	                
private:

	dMatrixT fTheta;
	int fNcrit, fN_times, fNumstep_crit;
	iArray2DT fTop, fBottom;
	dArray2DT fTHKforce, fGaussdisp;
	iArrayT fTopatoms, fBottomatoms, fToprow, fBottomrow, fShift;

	/** interpolation points */
	iArrayT fInterpolationNodes;
	
	/** theta values: [neighbor] x [time x n_theta_values] */
	ArrayT<dArray2DT> fThetaTable;
	
	/** displacement history: [n_boundary_atoms] x [time x ndof] */
	ArrayT<dArray2DT> fHistoryTablet, fHistoryTableb;
};

} /* namespace Tahoe */

#endif /* BRIDGING_ELEMENT */
#endif /* _FE_MANAGER_THK_H_ */
