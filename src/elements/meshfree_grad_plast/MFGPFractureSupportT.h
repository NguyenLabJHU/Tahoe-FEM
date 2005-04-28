/* $Id: MFGPFractureSupportT.h,v 1.2 2005-04-28 01:22:11 kyonten Exp $ */
#ifndef _MFGP_FRACTURE_T_H_
#define _MFGP_FRACTURE_T_H_

/* base class */
#include "MFGPElementSupportT.h"

/* direct members */
#include "dArray2DT.h"
#include "nVariArray2DT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */
class SolidMaterialT;
class FrontT;
class SamplingSurfaceT;

/** support for meshfree calculations including representation of cracks using
 * cutting surfaces */
class MFGPFractureSupportT: public MFGPElementSupportT
{
public:

	enum FractureCriterionT {kNoCriterion = 0,
	                       kMaxHoopStress = 1,
	                         kMaxTraction = 2,
	                            kAcoustic = 3};

	/** constructor */
	MFGPFractureSupportT(void);

	/** destructor */
	virtual ~MFGPFractureSupportT(void);

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

	/** check for extension of active crack fronts. Material properties are evaluated
	 * at the sampling points using the constitutive model and the displacement array.
	 * This call collects the list of facets returned by MeshFreeFractureSupportT::ResetFacets
	 * \param material pointer to the bulk constitutive model. If passes as NULL, no check
	 *        is performed, but the reset facet list is cleared.
	 * \param disp pointer to the nodal displacement data. Can by passed as NULL \e only if
	 *        the material pointer is also NULL.
	 * \param verbose pass as true to write debugging data to cout
	 * \return true of new facets have been inserted */
	bool CheckGrowth(SolidMaterialT* material, LocalArrayT* disp,
		bool verbose);
	
	/** initialization of meshless information. This method must be called once after 
	 * a call to MeshFreeElementSupportT::TakeParameterList */
	virtual void InitSupport(ostream& out, AutoArrayT<ElementCardT>& elem_cards, 
		const iArrayT& surface_nodes, int numDOF, int max_node_num, ModelManagerT* model);

	/** initialization of meshless information. This method must be called once after 
	 * a call to MeshFreeElementSupportT::TakeParameterList */
	virtual void InitSupport(int numDOF_displ, int numDOF_plast);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** translate integer to FractureCriterionT */
	static FractureCriterionT int2FractureCriterionT(int i);

private:

	/* sampling surface status code */
	enum SurfaceStatusT {kON, kMarked, kOFF};

	/* steps in InitSupport() */
	void InitCuttingFacetsAndFronts(ostream& out);
	void InitSamplingSurfaces(ifstreamT& in, ostream& out);

	/* initial active cracks from stream data */
	void InitializeFronts(ifstreamT& in, ostream& out);

	/* steps in checking growth */
	bool CheckFronts(SolidMaterialT& material, LocalArrayT& disp, bool verbose);
	bool CheckSurfaces(SolidMaterialT& material, LocalArrayT& disp, bool verbose);

	/* initialize the cutting facet database */
	void InitFacetDatabase(int num_facet_nodes);
	
	/* return the specified metric and returns the associated traction
	 * in the local frame defined by the transformation Q. n is the
	 * surface normal in the current configuration expressed in the
	 * global frame. Call only after configuring the meshfree field
	 * at the current point. Return value has sign convention that
	 * "more positive" is closer to failed */
	double ComputeCriterion(SolidMaterialT& material, const dMatrixT& Q,
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
inline int MFGPFractureSupportT::NumFacetNodes(void) const
{
	return fNumFacetNodes;
}

inline const dArray2DT& MFGPFractureSupportT::Facets(void) const
{
	return fFacets;
}

inline const ArrayT<int>& MFGPFractureSupportT::ResetFacets(void) const
{
	return fResetFacets;
}

inline const dArray2DT& MFGPFractureSupportT::InitTractions(void) const
{
	return fInitTractions;
}

/* returns true if the crack growth is possible */
inline bool MFGPFractureSupportT::HasActiveCracks(void) const
{
	if (fFrontList.Length() == 0 &&
	    fSamplingSurfaces.Length() == 0) return false;
	else return true;
}

/* fracture criterion */
inline MFGPFractureSupportT::FractureCriterionT
MFGPFractureSupportT::FractureCriterion(void) const
{
	return fCriterion;
}

} // namespace Tahoe 
#endif /* _MFGP_FRACTURE_T_H_ */
