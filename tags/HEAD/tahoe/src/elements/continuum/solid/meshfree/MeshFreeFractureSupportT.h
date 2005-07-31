/* $Id: MeshFreeFractureSupportT.h,v 1.1.1.1 2001-01-29 08:20:39 paklein Exp $ */
/* created: paklein (02/15/2000)                                          */

#ifndef _MESHFREE_FRACTURE_T_H_
#define _MESHFREE_FRACTURE_T_H_

/* base class */
#include "MeshFreeElementSupportT.h"

/* direct members */
#include "dArray2DT.h"
#include "nVariArray2DT.h"
#include "dSymMatrixT.h"

/* forward declarations */
class StructuralMaterialT;
class FrontT;
class SamplingSurfaceT;

class MeshFreeFractureSupportT: public MeshFreeElementSupportT
{
public:

	enum FractureCriterionT {kNoCriterion = 0,
	                       kMaxHoopStress = 1,
	                         kMaxTraction = 2,
	                            kAcoustic = 3};
	friend istream& operator>>(istream& in,
		MeshFreeFractureSupportT::FractureCriterionT& criterion);

	/* constructor */
	MeshFreeFractureSupportT(ifstreamT& in);

	/* destructor */
	~MeshFreeFractureSupportT(void);

	/* cutting facets */
	int NumFacetNodes(void) const;
	const dArray2DT& Facets(void) const;

	/* new facets data */
	const ArrayT<int>& ResetFacets(void) const;
	const dArray2DT& InitTractions(void) const;

	/* write output */
	void WriteOutput(ostream& out);

	/* initialize/finalize time increment */
	void InitStep(void) { return; }
	void CloseStep(void);
	void ResetStep(void); // restore last converged state
	
	/* returns true if the crack growth is possible */
	bool HasActiveCracks(void) const;
	
	/* fracture criterion */
	FractureCriterionT FractureCriterion(void) const;
	
protected:

	/* initialization */
	void InitSupport(ifstreamT& in, ostream& out, AutoArrayT<ElementCardT>& elem_cards,
		const iArrayT& surface_nodes, int numDOF, int max_node_num,
		const StringT& model_file, IOBaseT::FileTypeT format);

	/* returns true for crack growth */
	bool CheckGrowth(StructuralMaterialT& material, LocalArrayT& disp,
		bool verbose);

private:

	/* sampling surface status code */
	enum SurfaceStatusT {kON, kMarked, kOFF};

	/* steps in InitSupport() */
	void InitCuttingFacetsAndFronts(ifstreamT& in, ostream& out);
	void InitSamplingSurfaces(ifstreamT& in, ostream& out);

	/* initial active cracks from stream data */
	void InitializeFronts(ifstreamT& in, ostream& out);

	/* steps in checking growth */
	bool CheckFronts(StructuralMaterialT& material, LocalArrayT& disp, bool verbose);
	bool CheckSurfaces(StructuralMaterialT& material, LocalArrayT& disp, bool verbose);

	/* initialize the cutting facet database */
	void InitFacetDatabase(int num_facet_nodes);
	
	/* return the specified metric and returns the associated traction
	 * in the local frame defined by the transformation Q. n is the
	 * surface normal in the current configuration expressed in the
	 * global frame. Call only after configuring the meshfree field
	 * at the current point. Return value has sign convention that
	 * "more positive" is closer to failed */
	double ComputeCriterion(StructuralMaterialT& material, const dMatrixT& Q,
		const dArrayT& n, FractureCriterionT criterion, double critical_value,
		dArrayT& t_local);

	//TEMP
	void SetStreamPrefs(ostream& out);

private:

	/* crack database */
	int fNumFacetNodes;
	dArray2DT fFacets;
	nVariArray2DT<double> fFacetman;
	AutoArrayT<int> fResetFacets;
	dArray2DT fInitTractions; // tractions at initiation in local frame
	nVariArray2DT<double> fInitTractionMan;

	/* crack fronts */
	double fs_i;  // fraction of fs_u to trigger insertion
	double fda;   // crack extension increment
	double fda_s; // fraction of fda along which stress is sampled
	double fcone; // max angle to check
	int    fn_s;  // number of sampling points (around strain ahead)
	ArrayT<FrontT*> fFrontList;

	/* sampling surfaces */
	ArrayT<SamplingSurfaceT*> fSamplingSurfaces;

	/* failure stress */
	FractureCriterionT fCriterion;
	double fs_u; // scalar surface initiation criterion
	             // need if there are active front or sampling surfaces

	/* failure criterion work space */
	dSymMatrixT fhoop;
	dArrayT ftmp_nsd;
};

/* inlines */
inline int MeshFreeFractureSupportT::NumFacetNodes(void) const
{
	return fNumFacetNodes;
}

inline const dArray2DT& MeshFreeFractureSupportT::Facets(void) const
{
	return fFacets;
}

inline const ArrayT<int>& MeshFreeFractureSupportT::ResetFacets(void) const
{
	return fResetFacets;
}

inline const dArray2DT& MeshFreeFractureSupportT::InitTractions(void) const
{
	return fInitTractions;
}

/* returns true if the crack growth is possible */
inline bool MeshFreeFractureSupportT::HasActiveCracks(void) const
{
	if (fFrontList.Length() == 0 &&
	    fSamplingSurfaces.Length() == 0) return false;
	else return true;
}

/* fracture criterion */
inline MeshFreeFractureSupportT::FractureCriterionT
MeshFreeFractureSupportT::FractureCriterion(void) const
{
	return fCriterion;
}

#endif /* _MESHFREE_FRACTURE_T_H_ */
