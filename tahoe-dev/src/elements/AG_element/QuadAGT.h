/* $Id: QuadAGT.h,v 1.1 2006-04-30 04:09:15 regueiro Exp $ */
/* created: ebrahimi (02/01/2006) */
 #ifndef _QUAD_AG_T_H_
 #define _QUAD_AG_T_H_

 /* base class */
 #include "GeometryBaseT.h"
 
 //#include "ParameterListT.h"

 namespace Tahoe {

 /** 2D quadrilaterialAG parent domain */
 class QuadAGT: public GeometryBaseT
 {
 public:

     /** constructors */
     QuadAGT(int numnodes);
     
     /** destructor */
	 ~QuadAGT(void);

     /** return the geometry code */
     virtual GeometryT::CodeT Geometry(void) const { return kQuadrilateralAG; };

     /** evaluate the shape functions. See
      * GeometryBaseT::EvaluateShapeFunctions for documentation */
     virtual void EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na) const;

     /** evaluate the shape functions and gradients. See
      * GeometryBaseT::EvaluateShapeFunctions for documentation */
     virtual void EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na,
         dArray2DT& DNa) const;

     /** evaluate the shape functions and gradients. See
      * GeometryBaseT::SetLocalShape for documentation */
     virtual void SetLocalShape(dArray2DT& Na, ArrayT<dArray2DT>& Na_x,
         dArrayT& weights) const;

     /** compute gradients of the "bubble" modes */
     virtual void BubbleModeGradients(ArrayT<dArray2DT>& Na_x) const;

     /** set the values of the nodal extrapolation matrix */
     virtual void SetExtrapolation(dMatrixT& extrap) const;

     /** integration point gradient matrix */
     virtual void IPGradientTransform(int ip, dMatrixT& transform) const;

     /** return the local node numbers for each facet of the element
      * numbered to produce at outward normal in the order: vertex
      * nodes, mid-edge nodes, mid-face nodes */
     virtual void NodesOnFacet(int facet, iArrayT& facetnodes) const;
     virtual void NumNodesOnFacets(iArrayT& num_nodes) const;

     /** return the local node numbers for each edge of element where
      * the row number corresponds with the canonical numbering of the
      * edges in the element. Edges are only defined for 3D domain geometries. */
     virtual void NodesOnEdges(iArray2DT& nodes_on_edges) const;

     /** returns the nodes on each facet needed to determine neighbors
      * across facets */
     virtual void NeighborNodeMap(iArray2DT& facetnodes) const;

     /** return geometry and number of nodes on each facet */
     virtual void FacetGeometry(ArrayT<CodeT>& facet_geom,
         iArrayT& facet_nodes) const;

     /** return true if the given point is within the domain defined by
      * the list of coordinates.
      * \param coords list of coordinates defining the domain
      * \param point test point coordinates */
     virtual bool PointInDomain(const LocalArrayT& coords, const dArrayT& point) const;

     /** return the integration point whose domain contains the given point in the
      * parent domain coordinates */
     virtual int IPDomain(int nip, const dArrayT& coords) const;

     /** \name nodal subdomains, see GeometryBaseT for more information */
     /*@{*/
     /** subdomain geometry */
     virtual GeometryT::CodeT NodalSubDomainGeometry(void) const;

     /** number of nodes defining the nodal subdomain */
     virtual int NodalSubDomainNumPoints(void) const;

     /** compute the coordinates of the points defining the nodal subdomain */
     virtual void NodalSubDomainCoordinates(const LocalArrayT& coords, int node,
         LocalArrayT& subdomain_coords) const;
     /*@}*/
     
	/** sign function */
	int Sign(double a) const;
	
	/** calculating the AG kernel */
	double Gamma(double rs) const;
	
	/** calculating Derivative of the AG kernel */
	double GammaPrim(double rs) const;
	
protected:

	/** m, n */
	//double m_param;
	//double n_param;
	
	/** AGType */
	//int AGType; 
	
	//double Gamma0, Gamma_r, Gamma_s, GammaPrim_r, GammaPrim_s;
	
	/** pointer to AG parameter list */
	//ParameterListT* fQuadAG;
	
 };

 } /* namespace Tahoe */

 #endif /* _QUAD_AG_T_H_ */

