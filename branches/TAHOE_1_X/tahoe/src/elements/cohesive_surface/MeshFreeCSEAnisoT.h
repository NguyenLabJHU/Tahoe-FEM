/* $Id: MeshFreeCSEAnisoT.h,v 1.9 2004-01-05 07:34:30 paklein Exp $ */
/* created: paklein (06/08/2000) */

#ifndef _MF_CSE_ANISO_T_H_
#define _MF_CSE_ANISO_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "GeometryT.h"
#include "LocalArrayT.h"
#include "LocalArrayGroupT.h"
#include "pArrayT.h"
#include "nArrayGroupT.h"
#include "nMatrixGroupT.h"
#include "VariRaggedArray2DT.h"
#include "dArray2DT.h"
#include "nVariArray2DT.h"

namespace Tahoe {

/* forward declaration */
class MeshFreeFractureSupportT;
class MeshFreeSurfaceShapeT;
class SurfacePotentialT;

/** cohesive surfaces in a meshfree domain. Surfaces are
* added adaptively during the MeshFreeCSEAnisoT::RelaxSystem
* through communication with a MeshFreeFractureSupportT */
class MeshFreeCSEAnisoT: public ElementBaseT
{
public:

	/* constructor */
	MeshFreeCSEAnisoT(const ElementSupportT& support, const FieldT& field);

	/* destructor */
	~MeshFreeCSEAnisoT(void);
	
	/* form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* initialize class data */
	virtual void Initialize(void);

	/* start of new time sequence */
	virtual void InitialCondition(void);

	/* finalize time increment */
	virtual void CloseStep(void);

	/* resets to the last converged solution */
	virtual GlobalT::RelaxCodeT ResetStep(void);

	/* element level reconfiguration for the current solution */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

	/* solution calls */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);

	/* returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void); //not implemented

	/* append element equations numbers to the list */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/* writing output */
	virtual void RegisterOutput(void);
	virtual void WriteOutput(void);

	/* compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);

	/* appends group connectivities to the array (X -> geometry, U -> field) */
	virtual void ConnectsX(AutoArrayT<const iArray2DT*>& connects) const;
	virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	             AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;

	/* returns 1 if DOF's are interpolants of the nodal values */
	virtual int InterpolantDOFs(void) const;

	/* element status flags */
	enum StatusFlagT {kOFF = 0,
                      kON = 1,
                  kMarked = 2};

protected:

	/* print element group data */
	virtual void PrintControlData(ostream& out) const;

	/* element data */
	virtual void EchoConnectivityData(ifstreamT& in, ostream& out);

	/* tangent matrix and force vector */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);
	virtual void RHSDriver(void);

	/* write all current element information to the stream */
	virtual void CurrElementInfo(ostream& out) const;

private:

	/* initialize facets in the reset list */
	void InitializeNewFacets(void);

	/* set element work space dimensions */
	void SetNumberOfNodes(int nnd);

	/* operations with pseudo rank 3 (list in j) matrices */
	void u_i__Q_ijk(const dArrayT& u, const ArrayT<dMatrixT>& Q,
		dArrayT& nee_vec, dMatrixT& Qu);

	void Q_ijk__u_j(const ArrayT<dMatrixT>& Q, const dArrayT& u,
		dMatrixT& Qu);

protected:

	/* parameters */
	GeometryT::CodeT fGeometryCode;
	int fNumIntPts;
	int fOutputArea;
	int fMFElementGroup;

	/* meshfree domain support */
	MeshFreeFractureSupportT* fMFFractureSupport;

	/* meshfree surface support */
	MeshFreeSurfaceShapeT* fMFSurfaceShape;

	/* cohesive surface potentials */
	SurfacePotentialT* fSurfacePotential;

	/* local arrays */
	LocalArrayT fLocDisp;

	/* work space */
	double   fFractureArea;
	dMatrixT fNEEmat;
	dArrayT  fNEEvec;	

	/* coordinate transformation */
	dMatrixT fQ;     // t'_i = Q_ji t_j, where t' is in the local frame
	dArrayT  fdelta; // gap vector (local frame)
	dArrayT  fT;     // traction vector (global frame)
	dMatrixT fddU_l; // surface stiffness (local frame)
	dMatrixT fddU_g; // surface stiffness (global frame)
	ArrayT<dMatrixT> fdQ; // [nsd] : [nsd] x [nee]
		// list representation of rank 3 dQ_ij/du_k

	/* work space (for tangent) */
	dMatrixT fnsd_nee_1, fnsd_nee_2; // [nsd] x [nee]
	
	/* global equation numbers */
	VariRaggedArray2DT<int> fElemEqnosEX;

	/* facet active flags */
	AutoArrayT<StatusFlagT> fActiveFlag;
	
	/* dynamic work space managers */
	LocalArrayGroupT fLocGroup;
	nArrayGroupT<double>  fNEEArray;
	nMatrixGroupT<double> fNEEMatrix;
	nMatrixGroupT<double> fMatrixManager; // [nsd x [nee]

	/* needed to initialize decohesion laws */
	dArrayT fInitTraction;

	/* cohesive law variable storage */
	dArray2DT fd_Storage;      // [nel] x [nip] : [n_state]
	dArray2DT fd_Storage_last; // [nel] x [nip] : [n_state]
	nVariArray2DT<double> fd_Storage_man;      // [nel] x [nip] : [n_state]
	nVariArray2DT<double> fd_Storage_last_man; // [nel] x [nip] : [n_state]
};

} // namespace Tahoe 
#endif /* _MF_CSE_ANISO_T_H_ */
