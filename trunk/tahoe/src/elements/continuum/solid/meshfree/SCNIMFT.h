/* $Id */
#ifndef _SCNIMF_T_H_
#define _SCNIMF_T_H_

/* base class */
#include "ElementBaseT.h"

/* direct members */
#include "VariArrayT.h"
#include "LocalArrayT.h"
#include "MeshFreeNodalShapeFunctionT.h"
#include "SSMatSupportT.h"
#include "MaterialListT.h"
#include "nVariArray2DT.h"
#include "InverseMapT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "RaggedArray2DT.h"
#include "GeometryBaseT.h"
#include "MaterialListT.h"
#include "ScheduleT.h"

namespace Tahoe {

/** forward declarations */
class iGridManagerT;
class CommManagerT;
class dSPMatrixT; //TEMP
class InverseMapT;
class ifstreamT;
class ofstreamT;

#ifdef __QHULL__
class CompGeomT;
#endif

/** base class for particle types */
class SCNIMFT: public ElementBaseT
{
public:

	/** constructor */
	SCNIMFT(const ElementSupportT& support, const FieldT& field);
	SCNIMFT(const ElementSupportT& support);

	/** destructor */
	~SCNIMFT(void);
	
	/** initialization */
	virtual void Initialize(void);

	/** form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/** NOT implemented. Returns an zero force vector */
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);
			
	/** returns the energy as defined by the derived class types */
	virtual double InternalEnergy(void) { return 0.0; };
	
	/** resgiter for writing output. Uses the sub-class implementations
	 * of ParticleT::GenerateOutputLabels to register the particle group for
	 * output. Sub-classes also need to implemented the WriteOutput method. */
	virtual void RegisterOutput(void);

	/** write output. ParticleT::WriteOutput only writes search grid statistics.
	 * Sub-classes are responsible for writing data for each particle, given the
	 * variables names returned by ParticleT::GenerateOutputLabels. */
	virtual void WriteOutput(void);

	/* compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);

	/** trigger reconfiguration */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

	/** \name restart functions */
	/*@{*/
	/** write restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::ReadRestart implementation. */
	virtual void WriteRestart(ostream& out) const;

	/** read restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::WriteRestart implementation. */
	virtual void ReadRestart(istream& in);
	/*@}*/

	/** */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);
	
	/** */
	virtual void RHSDriver(void);
	
	/** Generate local equation numbers */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
						AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/** \name types needed for the Voronoi diagram calculation */
	/*@{*/
#ifndef __QHULL__
	/** Basic structure -- hullMap[i][j] gives index of jth vertex of structure i */
    typedef ArrayT<iArrayT> ConvexHullMap;

    /** Voronoi diagram facet structure is an array of convex hulls  */
    typedef ArrayT< ConvexHullMap > VoronoiDiagramMap;
#else	
	/** Basic structure -- hullMap[i][j] gives index of jth vertex of structure i */
    typedef CompGeomT::ConvexHullMap ConvexHullMap;

    /** Voronoi diagram facet structure is an array of convex hulls  */
    typedef CompGeomT::VoronoiDiagramMap VoronoiDiagramMap;
#endif
	/*@}*/

protected: /* for derived classes only */

	/** echo element connectivity data. Reads parameters that define
	 * which nodes belong to this ParticleT group. */
	virtual void EchoConnectivityData(ifstreamT& in, ostream& out);
	
	void ReadMaterialData(ifstreamT& in);
	
	void WriteMaterialData(ostream& out) const;
	
	MaterialListT* NewMaterialList(int nsd, int size);
	
	/** generate labels for output data */
	virtual void GenerateOutputLabels(ArrayT<StringT>& labels);

	/** return true if connectivities are changing */
	virtual bool ChangingGeometry(void) const;

	/** assemble particle mass matrix into LHS of global equation system
	 * \param mass mass associated with each particle type */
	void AssembleParticleMass(const dArrayT& mass);

	/** translate internal storage of bVector to Strain-Displacement matrix */	
	void bVectorToMatrix(double *bVector, dMatrixT& BJ);
	
	/** compute B matrices for strain smoothing/nodal integration */
	virtual void ComputeBMatrices(void);
	
	/** write out Voronoi diagram data */
	void VoronoiDiagramToFile(ofstreamT& vout);
	
	/** read in Voronoi diagram data */
	void VoronoiDiagramFromFile(ifstreamT& vin);
	
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
		SubListT& sub_sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& list_name) const;
	/*@}*/

protected:

	/** reference ID for sending output */
	int fOutputID;
	
	/** connectivities used to define the output set. Just an alias to the
	 * ParticleT::fGlobalTag. */
	iArray2DT fPointConnectivities;
	/*@}*/

	/** \name cached calculated values */
	/*@{*/
	dArray2DT fForce;
	nVariArray2DT<double> fForce_man;
	/*@{*/
	
	/** spatial dimensionality */
	int fSD;

	/** indices of nodes */
	iArrayT fNodes;
	/** the coordinates of the nodes */
	dArray2DT fDeloneVertices;
	
	/** midpoints of each of these are centroids of Voronoi facets */
	iArray2DT fDeloneEdges;
	
	/** Number of Delone edges that are not on the body bounday */
	int nInteriorDeloneEdges;
	
	/** additional edges with nodes as endpoints for boundary integration */
	//iArray2DT fBoundaryDeloneEdges;
	
	/** centroids of the facets corresponding to those edges */
	dArray2DT fBoundaryDeloneCentroids;
	
	/** additional edges associated only with one node */
	iArrayT fNonDeloneEdges;
	
	/** centroids of the facets for those edges */
	dArray2DT fNonDeloneCentroids;
	dArray2DT fNonDeloneNormals;
	
	/** dual of the Delone Edges -- area of Voronoi Facets */
	dArrayT fDualAreas;
	
	/** areas of boundary facets */
	dArrayT fIntegrationWeights;

	GeometryBaseT* fFakeGeometry;

	/** list of materials */
	MaterialListT* fMaterialList;
	SSMatSupportT* fSSMatSupport;
	
	/** workspace for strain smoothing */
	ArrayT<dArray2DT> bVectors;
	
	/* output control */
	iArrayT	fNodalOutputCodes;
	iArrayT	fElementOutputCodes;
	  	
	/* body force vector */
	const ScheduleT* fBodySchedule; /**< body force schedule */
	dArrayT fBody; /**< body force vector   */
	
	/** \name arrays with local ordering */
	/*@{*/
	LocalArrayT fLocInitCoords;   /**< initial coords with local ordering */
	//LocalArrayT fLocDisp;	      /**< displacements with local ordering  */ 
	/*@}*/
	
	/** shape functions */
	MeshFreeNodalShapeFunctionT* fNodalShapes;
	
	/** underlying Element connectivities. Needed only for MLS stuff right now */
	ArrayT<const iArray2DT*> fElementConnectivities;

	/** equation numbers */
	RaggedArray2DT<int> fEqnos;


#ifdef __QHULL__	
	CompGeomT* fVoronoi;
#else
	void* fVoronoi;
#endif

	ConvexHullMap fVoronoiCells;
	VoronoiDiagramMap fVoronoiFacetIndices;
	ArrayT<dArrayT> fVoronoiFacetAreas;
	ArrayT<dArray2DT> fVoronoiFacetNormals;
	dArrayT fVoronoiCellVolumes;
	dArray2DT fVoronoiVertices;

};

} /* namespace Tahoe */

#endif /* _SCNIMF_T_H_ */


