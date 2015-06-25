#include "FSDielectricElastomerQ1P0SurfaceT.h"
#include "FSDEMatQ1P02DT.h"
#include "FSDEMatSupportQ1P02DT.h"
#include "ShapeFunctionT.h"

/* From TLCBSurfaceT */
#include "ModelManagerT.h"
#include "MaterialListT.h"
#include "eIntegratorT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "InverseMapT.h"
#include "RowAutoFill2DT.h"
#include "OutputSetT.h"

#include "FEManagerT.h"
#include "TimeManagerT.h"

#include <iostream>
#include <algorithm>

using namespace std;

namespace Tahoe {

/* vector functions */
inline static void CrossProduct(const double* A, const double* B, double* AxB) {   
      AxB[0] = A[1]*B[2] - A[2]*B[1];
      AxB[1] = A[2]*B[0] - A[0]*B[2];
	AxB[2] = A[0]*B[1] - A[1]*B[0];
};

inline static double Dot(const double* A, const double* B){ 
      return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
};

inline static void Vector(const double* start, const double* end, double* v) {
      v[0] = end[0] - start[0];
      v[1] = end[1] - start[1];
      v[2] = end[2] - start[2];
};

inline static void Scale(double* v, double scale) {
      v[0] *= scale;
      v[1] *= scale;
      v[2] *= scale;
};

inline static void Sum(const double* A, const double* B, double* AB) {
      AB[0] = A[0] + B[0];
	AB[1] = A[1] + B[1];
	AB[2] = A[2] + B[2];
};

FSDielectricElastomerQ1P0SurfaceT::FSDielectricElastomerQ1P0SurfaceT(const ElementSupportT& support) :
    FSDielectricElastomerQ1P02DT(support),
    fLocCurrCoords(LocalArrayT::kCurrCoords), fSurfTension(0), fSurfaceCBSupport(NULL)
{
      SetName("dielectric_elastomer_Q1P0Elastocapillary");
}


/* Destructor */
FSDielectricElastomerQ1P0SurfaceT::~FSDielectricElastomerQ1P0SurfaceT()
{
      //cout << "jhf" << endl;
      /* TLCBSurface Stuff */
	/* free surface models */
	delete fSurfaceCBSupport;

      //cout << "dffgd" << endl;
}


// ------------------------------------------------specify parameters needed by the interface--------------------------------------------------------- //
void FSDielectricElastomerQ1P0SurfaceT::DefineParameters(ParameterListT& list) const
{
      /* Inherited from the most similar element */
      FSDielectricElastomerQ1P02DT::DefineParameters(list);

      /* Taking the output_surface */
      ParameterT output_surface(ParameterT::Boolean, "output_surface");
      output_surface.SetDefault(false);
      list.AddParameter(output_surface);
}

  /* ----------------------------- information about subordinate parameter lists ------------------------------------------------------------------------*/
void FSDielectricElastomerQ1P0SurfaceT::DefineSubs(SubListT& sub_list) const
{

      /* inherited from the most similar element */
      FSDielectricElastomerQ1P02DT::DefineSubs(sub_list);

      /* list of passivated surfaces - side set list */
      sub_list.AddSub("passivated_surface_ID_list", ParameterListT::ZeroOrOnce);

}


  // ---------------------------------------------------- accept parameter list ------------------------------------------------------------------------ //
void FSDielectricElastomerQ1P0SurfaceT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "FSDielectricElastomerQ1P0SurfaceT::TakeParameterList";

      /* inherited from the most similar element for bulk element */
      FSDielectricElastomerQ1P02DT::TakeParameterList(list);

      /* Shape Functions and related values */
      const ShapeFunctionT& shape = ShapeFunction();
      int nsd = shape.NumSD();                          // # of spatial dimensions in problem
      int nfs = shape.NumFacets();                      // # of total possible surface facets?
      int nsi = shape.FacetShapeFunction(0).NumIP();    // # IPs per surface face (2 x 1 = 2 for 1D surface)
      int nfn = shape.FacetShapeFunction(0).NumNodes(); // # nodes on each surface face?
      int nen = NumElementNodes();                      // # nodes on each element
      int nme = nen * nsd;                              // # of mechanical degree of freedoms per element

      /* mechanical due to surface */
      fAmm_mat2.Dimension(nme, nme);

      /* support for the surface model */
      fF_Surf_List.Dimension(nsi);

      /* Need to actually place values into fF_Surf_List when testing (identity) */
      for (int i = 0; i < fF_Surf_List.Length(); i++)
            fF_Surf_List[i].Dimension(nsd);

      /* DUMMY INITIALIZE fF_Surf_List - SPECIFY DEFORMATION GRADIENT (copied from TotalLagrangianCBSurfaceT) */
      fF_Surf_List[0].Identity();

      /* Back to normal flow */
      /* Support for surface models */
      fSurfaceCBSupport = new FSMatSupportT(nsd, nsi);
      fSurfaceCBSupport->SetContinuumElement(this);
      fSurfaceCBSupport->SetDeformationGradient(&fF_Surf_List);



      /* Do we need to redefine this in "canonical" normal order? */
      /* 2D */
      double normals_dat[4*2] = {
            1.0, 0.0,
            -1.0, 0.0,
            0.0, 1.0,
            0.0,-1.0,
      };
      dArray2DT normals(4, 2, normals_dat);

      /* Allocating the normal vectors */
      fNormal.Dimension(nfs);
      for (int i = 0; i < nfs; i++)
      {
            fNormal[i].Dimension(nsd);
            fNormal[i] = normals(i);
      }

      /* find parameter list for the bulk material */
      int num_blocks = list.NumLists("large_strain_element_block");
      if (num_blocks > 1)
      ExceptionT::GeneralFail(caller, "expecting only 1 not %d element blocks", num_blocks);
      const ParameterListT& block = list.GetList("large_strain_element_block");
      const ParameterListT& mat_list = block.GetListChoice(*this, "large_strain_material_choice");
      const ArrayT<ParameterListT>& mat = mat_list.Lists();
      const ParameterListT& bulk_params = mat[0]; /* getting bulk parameters */

      fSurfTension = bulk_params.GetParameter("gamma"); /* Reading the value of gamma (surface tension) */


      fT_0 = bulk_params.GetParameter("t_0");

      bool output_surface = list.GetParameter("output_surface");

      /* collect surface element information */
      ArrayT<StringT> block_ID;
      ElementBlockIDs(block_ID);
      ModelManagerT& model_manager = ElementSupport().ModelManager();
      model_manager.BoundingElements(block_ID, fSurfaceElements, fSurfaceElementNeighbors);

      /* determine normal type of each face */
      dMatrixT Q(nsd);
      dMatrixT jacobian(nsd, nsd-1);
      LocalArrayT face_coords(LocalArrayT::kCurrCoords, nfn, nsd); // kCurrCoords: current, kInitCoords: initial
      iArrayT face_nodes(nfn), face_nodes_index(nfn);
      ElementSupport().RegisterCoordinates(face_coords);
      fSurfaceElementFacesType = fSurfaceElementNeighbors;
      fSurfaceElementFacesType = -1;

      /* nodes per surface type */
      RowAutoFill2DT<int> nodes_on_surfaces(nfs, 25);

      for (int i = 0; i < fSurfaceElements.Length(); i++) /* loop over the surface elements for one element fSurfaceElements.Length() = 1 then i = 0*/
      {
            for (int j = 0; j < fSurfaceElementNeighbors.MinorDim(); j++) /* loop over faces */
            {
                  if (fSurfaceElementNeighbors(i,j) == -1) /* no neighbor => surface for one element (i,j) = (0,0), (0,1), (0,2) and (0,3)*/
                  {
                        /* face parent domain */
                        const ParentDomainT& surf_shape = shape.FacetShapeFunction(j);

                        /* collect coordinates of face nodes */
                        ElementCardT& element_card = ElementCard(fSurfaceElements[i]);
                        shape.NodesOnFacet(j, face_nodes_index);  // fni = 4 nodes of surface face
                        face_nodes.Collect(face_nodes_index, element_card.NodesX());
                        face_coords.SetLocal(face_nodes);

                        /* face normal (using 1st integration point) */
                        surf_shape.DomainJacobian(face_coords, 0, jacobian);
                        surf_shape.SurfaceJacobian(jacobian, Q);  // Last column of Q is normal vector to surface face

                        /* match to face normal types, i.e. normals_dat */
                        int normal_type = -1;
                        for (int k = 0; normal_type == -1 && k < fNormal.Length(); k++)
                        {
                              if ((Q.DotCol(nsd-1, fNormal[k]) - 1.0) < -kSmall)
                                    normal_type = -1;
                              else
                                    normal_type = k;
                        }
                        /* no match */
                        if (normal_type == -1)
                              ExceptionT::GeneralFail(caller, "could not classify normal on face %d of element %d", j+1, fSurfaceElements[i]+1);

                        /* store */
                        fSurfaceElementFacesType(i,j) = normal_type;

                        /* collect face nodes */
                        if (output_surface) nodes_on_surfaces.AppendUnique(normal_type, face_nodes);
                  }
            }
      }

      /* copy nodes on face types data */
      if (output_surface) {
            fSurfaceNodes.Dimension(nodes_on_surfaces.MajorDim());
            for (int i = 0; i < fSurfaceNodes.Length(); i++) {
                  fSurfaceNodes[i].Dimension(nodes_on_surfaces.MinorDim(i));
                  fSurfaceNodes[i].Copy(nodes_on_surfaces(i));
            }
      }

      /* process passivated surfaces */
      const ParameterListT* passivated_surfaces = list.List("passivated_surface_ID_list");
      if (passivated_surfaces) {
            ArrayT<StringT> ss_ID;
            StringListT::Extract(*passivated_surfaces, ss_ID);
            if (ss_ID.Length() > 0) {

                  /* need map into the surface element list */
                  InverseMapT surf_elem_map;
                  surf_elem_map.SetOutOfRange(InverseMapT::Throw);
                  surf_elem_map.SetMap(fSurfaceElements);

                  /* model manager */
                  ModelManagerT& model = ElementSupport().ModelManager();

                  /* loop over side set ID's */
                  for (int i = 0; i < ss_ID.Length(); i++) {
                        const StringT& id = ss_ID[i];

                          /* side set parameters */
                        iArray2DT sides = model.SideSet(id);
                        const StringT& block_id = model.SideSetGroupID(id);

                          if (sides.MajorDim() > 0) {

                              /* convert to element numbering within the group */
                              iArrayT elems(sides.MajorDim());
                              sides.ColumnCopy(0, elems);
                              BlockToGroupElementNumbers(elems, block_id);
                              sides.SetColumn(0, elems);

                              /* mark passivated faces */
                              for (int j = 0; j < sides.MajorDim(); j++) {
                                    int elem = sides(j,0);
                                    int s_elem = surf_elem_map.Map(elem);
                                    int side = sides(j,1);

                                    /* mark as non-surface (by giving negative neighbor id) */
                                    fSurfaceElementNeighbors(s_elem,side) = -2;
                              }
                        }
                  }
            }
      }

      // cout << "\033[1;31mTakeParameterList looks ok\033[0m\n";
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* ------------------------------------------------ calculate the LHS of residual, or element stiffness matrix ------------------------------------------------*/
void FSDielectricElastomerQ1P0SurfaceT::FormStiffness(double constK)
{
	//cout.setf(ios::fixed,ios::floatfield);
	//cout.precision(2);
	//cout << fTimeInit << endl;



	//cout << "Can you see this?" << endl;
	const char caller[] = "FSDielectricElastomerQ1P0SurfaceT::FormStiffness";

	/* Inherited */
	FSDielectricElastomerQ1P02DT::FormStiffness(constK);
	dMatrixT::SymmetryFlagT format = (fLHS.Format()
			== ElementMatrixT::kNonSymmetric)
            ? dMatrixT::kWhole
            : dMatrixT::kUpperOnly;

      const ShapeFunctionT& shape = ShapeFunction();
      int nsd = shape.NumSD();                          // # of spatial dimensions in problem
      int nfs = shape.NumFacets();                      // # of total possible element faces
      int nsi = shape.FacetShapeFunction(0).NumIP();    // # IPs per surface face
      int nfn = shape.FacetShapeFunction(0).NumNodes(); // # nodes on each surface face
      int nen = NumElementNodes();                      // # nodes in bulk element
      int nme = nen * nsd;

      iArrayT counter(nen);

      /* loop over surface elements */
      dMatrixT jacobian(nsd, nsd-1);
      iArrayT face_nodes(nfn), face_nodes_index(nfn);
      LocalArrayT face_coords(LocalArrayT::kCurrCoords, nfn, nsd);
      ElementSupport().RegisterCoordinates(face_coords);

      dMatrixT K1;
      dMatrixT K2;
      dMatrixT K_Total;
      dArrayT fB;
      dMatrixT K3(nen);

      K1.Dimension(nen, nen);
      K2.Dimension(nen, nen);
      fB.Dimension(nen);
      K_Total.Dimension(nen, nen);

 	 double CurrTime = ElementSupport().Time(); // Obtaining current time

 	 double fNewSurfTension = min(fSurfTension, (CurrTime*fSurfTension)/fT_0); // Ramping up the surface tension

 	 //fNewSurfTension = fSurfTension;
 	 //cout << fNewSurfTension << endl;
 	 //cout << fSurfTension << endl;

 	 // Should be loop over number of elements on surface!! //
 	 for (int i = 0; i < fSurfaceElements.Length(); i++)
 	 {

 		 /* bulk element information */
 		 int element = fSurfaceElements[i];

 		 if (element == CurrElementNumber()) // Is the current element a surface element?
 		 {
 			 const ElementCardT& element_card = ElementCard(element);
 			 fLocInitCoords.SetLocal(element_card.NodesX()); /* reference coordinates over bulk element (collects first x coords and then y coords) i.e.  [x1 x2 x3 x4 y1 y2 y3 y4] */
 			 fLocDisp.SetLocal(element_card.NodesU()); /* displacements over bulk element */

 			 //for (int k = 0; k < fLocInitCoords.Length(); k ++){
 			 //	cout << fLocInitCoords[k] << endl;
 			 //}
 			 fB = 0.0;
 			 K1 = 0.0;
 			 K2 = 0.0;
 			 K3 = 0.0;
 			 K_Total = 0.0;
 			 fAmm_mat2 = 0.0;


 			 for (int j = 0; j < fSurfaceElementNeighbors.MinorDim(); j++) /* loop over faces */
 			 {
 				 if (fSurfaceElementNeighbors(i,j) == -1) /* no neighbor => surface */
 				 {
 					 /* face parent domain */
 					 const ParentDomainT& surf_shape = shape.FacetShapeFunction(j);

 					 /* collect coordinates of face nodes */
 					 ElementCardT& element_card = ElementCard(fSurfaceElements[i]);
 					 shape.NodesOnFacet(j, face_nodes_index);  // fni = 4 nodes of surface face
 					 face_nodes.Collect(face_nodes_index, element_card.NodesX());
 					 face_coords.SetLocal(face_nodes);

 					 K1 = 0.0; K2 = 0.0; fB = 0.0;

 					 K1(0,0) = 1.0; K1(1,1) = 1.0;
 					 K1(2,2) = 1.0; K1(3,3) = 1.0;
 					 K1(0,2) = -1.0; K1(2,0) = -1.0;
 					 K1(3,1) = -1.0; K1(1,3) = -1.0;

 					 double x_1 = face_coords[0];
 					 double x_2 = face_coords[1];
 					 double y_1 = face_coords[2];
 					 double y_2 = face_coords[3];

 					 /* cout << nfn;
                  	cout << "i = " << i << " , j = " << j << endl;
                  	cout << "fSurfaceElementFacesType("<< i << ", " << j << ") = " << fSurfaceElementFacesType(i,j) << endl;
                   	cout << "fSurfaceElements = " << fSurfaceElements[i] << endl;
                   	cout << "Coordinates:" << endl;
                  	cout << "(x1, y1) = " << "(" << x_1 << ", " << y_1 << ")" << endl;
                   	cout << "(x2, y2) = " << "(" << x_2 << ", " << y_2 << ")" << endl;
                   	cout << endl; */

 					 /* For 2D cubic element: nen = 4 */
 					 fB[0] = (x_1 - x_2);
 					 fB[1] = (y_1 - y_2);
 					 fB[2] = (x_2 - x_1);
 					 fB[3] = (y_2 - y_1);

 					 K2.Outer(fB, fB); /* Multiplying B to BT */

 					 /* Length of the surface */
 					 double L_e = sqrt((x_1 - x_2)*(x_1 - x_2) + (y_1 - y_2)*(y_1 - y_2));

 					 double coeff1 =  fNewSurfTension/L_e;
 					 double coeff2 = -fNewSurfTension/(L_e*L_e*L_e);

 					 K_Total = 0.0;

 					 K1 *= coeff1;
 					 K2 *= coeff2;
 					 K_Total += K1;
 					 K_Total += K2;



 					 /* Constructing fAmm_mat */
 					 int normaltype = fSurfaceElementFacesType(i, j);
 					 counter = CanonicalNodes(normaltype);

 					 for (int n = 0; n < nen; n++) {
 					 	 fAmm_mat2(counter[n], counter[0]) = fAmm_mat2(counter[n], counter[0]) + K_Total(n ,0);
 					 	 fAmm_mat2(counter[n], counter[1]) = fAmm_mat2(counter[n], counter[1]) + K_Total(n ,1);
 					 	 fAmm_mat2(counter[n], counter[2]) = fAmm_mat2(counter[n], counter[2]) + K_Total(n ,2);
 					 	 fAmm_mat2(counter[n], counter[3]) = fAmm_mat2(counter[n], counter[3]) + K_Total(n ,3);
 					 }



 					 /*if (normaltype == 3) {
 						 for (int n = 0; n < nen; n++)
 							 for (int m = 0; m < nen; m++)
 								 fAmm_mat2(counter[n], counter[m]) = fAmm_mat2(counter[n], counter[m]) + K_Total(n ,m);
 					 } else if (normaltype == 0) {
 						 for (int n = 0; n < nen; n++)
 							 for (int m = 0; m < nen; m++)
 								 fAmm_mat2(counter[n], counter[m]) = fAmm_mat2(counter[n], counter[m]) + K_Total(n ,m);
 					 } else if (normaltype == 2) {
 						 for (int n = 0; n < nen; n++)
 							 for (int m = 0; m < nen; m++)
 								 fAmm_mat2(counter[n], counter[m]) = fAmm_mat2(counter[n], counter[m]) + K_Total(n ,m);
 					 } else {
 						 for (int n = 0; n < nen; n++)
 							 for (int m = 0; m < nen; m++)
 								 fAmm_mat2(counter[n], counter[m]) = fAmm_mat2(counter[n], counter[m]) + K_Total(n ,m);
 					 } */

 					 int order = fIntegrator->Order();
 					 //cout << constK << endl;
 					 if (order == 2)
 						 fAmm_mat2 *= constK;


 					 /* End of constructing fAmm_mat */

 				 } /* End of if */

 			 } /* End of surface edge loop */

 			 fLHS.AddBlock(0, 0, fAmm_mat2);

 			 // Saving the fLHS matrix
 			  /*ofstream myLHS;
            	myLHS.open("fAmm_mat2.txt");
            	for (int i = 0; i < fAmm_mat2.Rows(); i++)
            	{
            		for (int j = 0; j < fAmm_mat2.Cols(); j++)
            		{
            			// myLHS << "fLHS(" << i << "," << j << ")= " << fLHS(i, j); // List the values of fLHS(i,j)
            			if (fAmm_mat2(i, j) == 0)
            				myLHS << "0.00000" << " ";
            			else
            				myLHS << fAmm_mat2(i, j) << " "; // See the matrix form of fLHS for a quick look
            		}
            	myLHS << endl;
            	}
            	myLHS.close()*/
 		 }
 	 } /* End of element loop */

 	 // cout << "\033[1;31mFormStiffness looks good!\033[0m" << endl;

} /* End of Fucntion FormStiffness */

/* ------------------------------------------------ Compute RHS, or residual of element equations ------------------------------------------------*/
void FSDielectricElastomerQ1P0SurfaceT::FormKd(double constK)
{
      const char caller[] = "FSDielectricElastomerQ1P0SurfaceT::FormKd";

      /* Inherited */
      FSDielectricElastomerQ1P02DT::FormKd(constK);

      /* dimensions */
      const ShapeFunctionT& shape = ShapeFunction();
      int nsd = shape.NumSD();                          // # of spatial dimensions in problem
      int nfs = shape.NumFacets();                      // # of total possible element faces
      int nsi = shape.FacetShapeFunction(0).NumIP();    // # IPs per surface face
      int nfn = shape.FacetShapeFunction(0).NumNodes(); // # nodes on each surface face
      int nen = NumElementNodes();                      // # nodes in bulk element
      int nme = nen * nsd;

      iArrayT counter(nen);

      /* loop over surface elements */
      dMatrixT jacobian(nsd, nsd-1);
      LocalArrayT face_coords(LocalArrayT::kCurrCoords, nfn, nsd); // kCurrCoords = current coordinates
      iArrayT face_nodes(nfn), face_nodes_index(nfn);
      ElementSupport().RegisterCoordinates(face_coords);

      dArrayT fD;
      dArrayT R_Total;
      dArrayT R;

      fD.Dimension(nen);
      R_Total.Dimension(nme);
      R.Dimension((nsd + 1) * nen);

  	  double CurrTime = ElementSupport().Time();

  	  double fNewSurfTension = min(fSurfTension, (CurrTime*fSurfTension)/fT_0);

  	  //cout << fNewSurfTension << endl;

  	 //fNewSurfTension = fSurfTension;

      for (int i = 0; i < fSurfaceElements.Length(); i++)
      {
            /* bulk element information */

    	  	int element = fSurfaceElements[i];
            if (element == CurrElementNumber()) // Is the current element a surface element?
            {
            	const ElementCardT& element_card = ElementCard(element);
            	fLocInitCoords.SetLocal(element_card.NodesX()); /* reference coordinates over bulk element */
            	fLocDisp.SetLocal(element_card.NodesU()); /* displacements over bulk element */

            	fD = 0.0;
            	R_Total = 0.0;
            	R = 0.0;

            	for (int j = 0; j < fSurfaceElementNeighbors.MinorDim(); j++) /* loop over faces */
            	{
            		if (fSurfaceElementNeighbors(i,j) == -1) /* no neighbor => surface */
            		{
                        	/* face parent domain */
                        	const ParentDomainT& surf_shape = shape.FacetShapeFunction(j);

                        	/* collect coordinates of face nodes */
                        	ElementCardT& element_card = ElementCard(fSurfaceElements[i]);
                        	shape.NodesOnFacet(j, face_nodes_index);  // fni = 4 nodes of surface face
                        	face_nodes.Collect(face_nodes_index, element_card.NodesX());
                        	face_coords.SetLocal(face_nodes);

                        	double x_1 = face_coords[0];
                        	double x_2 = face_coords[1];
                        	double y_1 = face_coords[2];
                        	double y_2 = face_coords[3];

                        	fD = 0.0;
                        	// Length of the surface
                        	double L_e = sqrt((x_1 - x_2)*(x_1 - x_2) + (y_1 - y_2)*(y_1 - y_2));

                        	double coeff3 = -fNewSurfTension/(L_e);

                        	// cout << "coeff3= " << coeff3 << endl;

                        	fD[0] = coeff3*(x_1 - x_2);
                        	fD[1] = coeff3*(y_1 - y_2);
                        	fD[2] = coeff3*(x_2 - x_1);
                        	fD[3] = coeff3*(y_2 - y_1);

                        	//                        D *= coeff3;
                        	//R_Total = 0.0;

                        	// cout << D[0] << D[1] << D[2] << D[3] << endl;

                        	int normaltype = fSurfaceElementFacesType(i,j);
                        	counter = CanonicalNodes(normaltype);

                        	R_Total[counter[0]] = R_Total[counter[0]] + fD[0];
                        	R_Total[counter[1]] = R_Total[counter[1]] + fD[1];
                        	R_Total[counter[2]] = R_Total[counter[2]] + fD[2];
                        	R_Total[counter[3]] = R_Total[counter[3]] + fD[3];



                        	/*if (normaltype == 3)
                        		for (int m = 0; m < nen; m++)
                        			R_Total[counter[m]] = R_Total[counter[m]] + fD[m];
                        	else if (normaltype == 0)
                        		for (int m = 0; m < nen; m++)
                        			R_Total[counter[m]] = R_Total[counter[m]] + fD[m];
                        	else if (normaltype == 2)
                        		for (int m = 0; m < nen; m++)
                        			R_Total[counter[m]] = R_Total[counter[m]] + fD[m];
                        	else
                        		for (int m = 0; m < nen; m++)
                        			R_Total[counter[m]] = R_Total[counter[m]] + fD[m];*/


            		} /* End of if */

            	}  /* End of surface edge loop */

            	R.CopyIn(0, R_Total);
            	fRHS += R;


            	/* Saving RHS to a file */
            	/* ofstream myRHS;
            	myRHS.open("fRHS_with_surface.txt");
            	for (int i = 0; i < fRHS.Length(); i++)
            	{
                  	  myRHS << "fRHS(" << i << ") = " << fRHS[i] << endl;
            	}
            	myRHS.close(); */
            	}
      	  } /* End of element loop */
      	  // cout << "\033[1;31mFormKd looks good\033[0m" << endl;

}
/***********************************************************************
 * Protected
 ***********************************************************************/
iArrayT FSDielectricElastomerQ1P0SurfaceT::CanonicalNodes(const int normaltype)
{
	/* Return nodes for canonical (psi, eta) element based on normal type */
	int nen = NumElementNodes();
	iArrayT counter(nen);
	if (normaltype == 3) {
		counter[0] = 0;
		counter[1] = 1;
		counter[2] = 2;
		counter[3] = 3;
	} else if (normaltype == 0) {
		counter[0] = 2;
		counter[1] = 3;
		counter[2] = 4;
		counter[3] = 5;
	} else if (normaltype == 2) {
		counter[0] = 4;
		counter[1] = 5;
		counter[2] = 6;
		counter[3] = 7;
	} else {
		counter[0] = 6;
		counter[1] = 7;
		counter[2] = 0;
		counter[3] = 1;
	}


	return counter;
}

/***********************************************************************
 * Private
 ***********************************************************************/

  // ------------------------------------------------ extrapolate from integration points and compute output nodal/element values ------------------------ //
void FSDielectricElastomerQ1P0SurfaceT::ComputeOutput(const iArrayT& n_codes,
      dArray2DT& n_values, const iArrayT& e_codes, dArray2DT& e_values)
{
      // cout << "\033[1;31mCan you see this?\033[0m" << endl;
      /* inherited - bulk contribution */
      FSDielectricElastomerQ1P02DT::ComputeOutput(n_codes, n_values, e_codes, e_values);

}

} // namespace Tahoe
