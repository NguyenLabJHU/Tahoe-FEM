/* $Id: UpLagFiberCompT.cpp,v 1.6 2006-10-31 15:54:10 rjones Exp $ */
/* created: paklein (07/03/1996) */
#include "UpLagFiberCompT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "toolboxConstants.h"
#include "ParameterContainerT.h"
#include "ModelManagerT.h"
#include "SolidMaterialT.h"
#include "ShapeFunctionT.h"
#include "Traction_CardT.h"

#include "FSFiberMatT.h"
#include "FSFiberMatSupportT.h"
#include "FSFiberMatListT.h"

#include "VariLocalArrayT.h"

using namespace Tahoe;

/* constructor */
UpLagFiberCompT::UpLagFiberCompT(const ElementSupportT& support):
	UpdatedLagrangianT(support),
	fFiberSupport(NULL)
{
	SetName("uplag_fiber_comp_planar");
}

UpLagFiberCompT::~UpLagFiberCompT(void) {
	delete fFiberSupport;
}

/* information about subordinate parameter lists */
void UpLagFiberCompT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidElementT::DefineSubs(sub_list);	

	/* element block/material specification */
	sub_list.AddSub("fiber_comp_element_block", ParameterListT::OnePlus);
	sub_list.AddSub("fiber_orientations", ParameterListT::OnePlus);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* UpLagFiberCompT::NewSub(const StringT& name) const
{
	/* inherited */
	if (name == "fiber_comp_element_block")
	{
		ParameterContainerT* block = new ParameterContainerT(name);
		
		/* list of element block ID's (defined by ElementBaseT) */
		block->AddSub("block_ID_list", ParameterListT::Once);
	
		/* choice of materials lists */
		block->AddSub("fiber_comp_material", ParameterListT::Once);
	
		/* set this as source of subs */
		block->SetSubSource(this);
		
		return block;
	}
	else if (name == "fiber_orientations") /* fiber orientations */
	{
    ParameterContainerT* choice = new ParameterContainerT(name);
    choice->SetListOrder(ParameterListT::Choice);
		choice->SetSubSource(this);

		{
		ParameterContainerT ss_spec("side_set");
		ss_spec.AddParameter(ParameterT::Word, "side_set_ID");
		ParameterT coord_sys(ParameterT::Enumeration, "coordinate_system");
		coord_sys.AddEnumeration("global", Traction_CardT::kCartesian);
		coord_sys.AddEnumeration( "local", Traction_CardT::kLocal);
		coord_sys.SetDefault(Traction_CardT::kCartesian);
		ss_spec.AddParameter(coord_sys);
		ss_spec.AddSub("DoubleList", ParameterListT::OnePlus); 		

		choice->AddSub(ss_spec);
		}

		{
		ParameterContainerT ellipsoid("ellipsoid");

		ParameterT block_ID(fID, "block_ID");
		block_ID.SetDefault("all");
		ellipsoid.AddParameter(block_ID);

		LimitT lower(0.0, LimitT::Lower);
		ParameterT Rx(ParameterT::Double, "Rx");
		Rx.AddLimit(lower);
		Rx.SetDefault(1.0);
		ellipsoid.AddParameter(Rx);
		ParameterT Ry(ParameterT::Double, "Ry");
		Ry.AddLimit(lower);
		Ry.SetDefault(1.0);
		ellipsoid.AddParameter(Ry);
		ParameterT Rz(ParameterT::Double, "Rz");
		Rz.AddLimit(lower);
		Rz.SetDefault(1.0);
		ellipsoid.AddParameter(Rz);
		ParameterT Cx(ParameterT::Double, "Cx");
		Cx.SetDefault(0.0);
		ellipsoid.AddParameter(Cx);
		ParameterT Cy(ParameterT::Double, "Cy");
		Cy.SetDefault(0.0);
		ellipsoid.AddParameter(Cy);
		ParameterT Cz(ParameterT::Double, "Cz");
		Cz.SetDefault(0.0);
		ellipsoid.AddParameter(Cz);

		ParameterT normal(ParameterT::Enumeration, "projection_normal");
		normal.AddEnumeration("x", 1);
		normal.AddEnumeration("y", 2);
		normal.AddEnumeration("z", 3);
		normal.SetDefault(3);
		ellipsoid.AddParameter(normal);

		ParameterT coord_sys(ParameterT::Enumeration, "coordinate_system");
		coord_sys.AddEnumeration("cartesian", UpLagFiberCompT::kCartesian);
		coord_sys.AddEnumeration("polar", UpLagFiberCompT::kPolar);
		coord_sys.SetDefault(UpLagFiberCompT::kCartesian);
		ellipsoid.AddParameter(coord_sys);
		ellipsoid.AddSub("DoubleList", ParameterListT::OnePlus); 		

		ellipsoid.SetDescription("((x-Cx)/Rx)^2 + ((y-Cy)/Ry)^2 + ((z-Cz)/Rz)^2 = 1");

		choice->AddSub(ellipsoid);
		}

		return choice;
	}
	else /* inherited */
		return SolidElementT::NewSub(name);
}

/* describe the parameters needed by the interface */
void UpLagFiberCompT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	UpdatedLagrangianT::DefineParameters(list);
}

/* accept parameter list */
void UpLagFiberCompT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	UpdatedLagrangianT::TakeParameterList(list);
	/*store fibers in element list*/
	int num_elem = NumElements();
	fFiber_list.Dimension(num_elem);
	ReadFiberVec(list);
	for (int i = 0; i < NumElements(); i++)
	{
//		cout << "\nelement: "<<i
//			 << "\t"<<fFiber_list[i];
		int num_fibers = fFiber_list[i].MajorDim();
		if (num_fibers == 0)
			ExceptionT::GeneralFail("UpLagFiberCompT::TakeParameterList",
			 "Fiber orientations not specified for element %d", i);
	}
	
}

/* extract the list of material parameters */
void UpLagFiberCompT::CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const
{
	const char caller[] = "UpLagFiberCompT::CollectMaterialInfo";
	
	/* initialize */
	mat_params.Clear();

	/* set materials list name */
	mat_params.SetName("fiber_comp_material");
	
	/* collected material parameters */
	int num_blocks = all_params.NumLists("fiber_comp_element_block");
	for (int i = 0; i < num_blocks; i++) {

		/* block information */	
		const ParameterListT& block = all_params.GetList("fiber_comp_element_block", i);
		
		/* collect material parameters */
		const ParameterListT& mat_list = block.GetList(mat_params.Name());
		const ArrayT<ParameterListT>& mat = mat_list.Lists();
		mat_params.AddList(mat[0]);
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* construct a new material support and return a pointer */
MaterialSupportT* UpLagFiberCompT::NewMaterialSupport(MaterialSupportT* p) const
{
	/* allocate */
	if (NumDOF() != NumSD())
	    ExceptionT::GeneralFail("UpLagFiberCompT::NewMaterialSupport", "ndof != nsd not supported");

	if (!p) p = new FSFiberMatSupportT(NumDOF(), NumIP());

	/* inherited initializations */
	FiniteStrainT::NewMaterialSupport(p);
	
	/* set fiber orientation vectors */
	FSFiberMatSupportT* ps = TB_DYNAMIC_CAST(FSFiberMatSupportT*, p);
	if (ps) {
		ps->SetFibers(&fFiber_list);
	}

	return p;
}

/* construct materials manager and read data */
MaterialListT* UpLagFiberCompT::NewMaterialList(const StringT& name, int size)
{
	/* resolve number of spatial dimensions */
	/* no match */
	if (name != "fiber_comp_material")
		return NULL;

	if (size > 0)
	{
		/* material support */
		if (!fFiberSupport) {
			fFiberSupport = TB_DYNAMIC_CAST(FSFiberMatSupportT*, NewMaterialSupport());
			if (!fFiberSupport) ExceptionT::GeneralFail("UpLagFiberCompT::NewMaterialList");
		}
		/* allocate */
		return new FSFiberMatListT(size, *fFiberSupport);
	}
	else
		return new FSFiberMatListT;
}

/* read in fiber orientation information information */
void UpLagFiberCompT::ReadFiberVec(const ParameterListT& list)
{
	const char caller[] = "UpLagFiberCompT::TakeFiberVec";

	int num_sets = list.NumLists("fiber_orientations");
	if (num_sets > 0)
	{
		for (int i = 0; i < num_sets; i++)
		{
			const ParameterListT& fibers = list.GetListChoice(*this, "fiber_orientations",i);

			if (fibers.Name() == "side_set") 
				ReadSideSetVec(fibers);
			else if (fibers.Name() == "ellipsoid") 
				ReadAnalyticVec(fibers);
			else
				ExceptionT::GeneralFail(caller, "invalid surface specification");
		}
	}
		else
			ExceptionT::GeneralFail(caller, "no fibers defined");
}

void UpLagFiberCompT::ReadSideSetVec(const ParameterListT& fibers)
{
	const char caller[] = "UpLagFiberCompT::ReadSideSetVec";

		int nsd = NumSD();
		/* model manager */
		ModelManagerT& model = ElementSupport().ModelManager();
	
		/* temp space */
		StringT block_ID;    /*block id of sideset*/
		iArray2DT localsides; /*numside x 1 (elem number) + numfacetnodes (node numbers)*/ 
		Traction_CardT::CoordSystemT coord_sys;  /*global cartesian or local element coords*/
		dArray2DT values;     /*p_vec*/

		/* nodes on element facets */
		iArrayT num_facet_nodes;
		fShapes->NumNodesOnFacets(num_facet_nodes);
		
		/* coordinates of facet nodes  register with initial coordinates */
		LocalArrayT coords(LocalArrayT::kInitCoords);
		VariLocalArrayT coord_man(25, coords, nsd);
		ElementSupport().RegisterCoordinates(coords);

		iArrayT facet_nodes_loc;  /*facet nodes in local element numbering*/
		iArrayT facet_nodes_glob;  /*facet nodes in global numbering*/

		/*jacobian of surface mapping*/
		dMatrixT jacobian(nsd, nsd - 1);
		/*rotation tensor*/
		dMatrixT Q(nsd);
		dMatrixT Qbar(nsd);
		Qbar.Identity();

			
		/* side set */
		const StringT& ss_ID = fibers.GetParameter("side_set_ID");
		/*reads in element and facet numbers of sideset SS_ID*/
		localsides = model.SideSet(ss_ID);
		/*number of sides in set*/
		int num_sides = localsides.MajorDim();
		if (num_sides > 0)
		{
				/*block ID of set set*/
				block_ID = model.SideSetGroupID(ss_ID);
				coord_sys = Traction_CardT::int2CoordSystemT(fibers.GetParameter("coordinate_system"));

				/* switch to elements numbering within the group */
				iArray2DT& side_set = localsides;            /*sides info for set i*/
				/* side a: elem #, facet #*/
				iArrayT elems(num_sides);
				/*copy element numbers of sideset into iArrayT elems*/
				side_set.ColumnCopy(0, elems);
				/*convert from local block numbering to global group numbering of elements*/
				BlockToGroupElementNumbers(elems, block_ID);
				/*copy group element numbering in side_set array*/
				side_set.SetColumn(0, elems);

				/* all facets in set must have the same number of nodes */
				int num_nodes = num_facet_nodes[side_set(0,1)];
				for (int f = 0; f < num_sides; f++)
					if (num_facet_nodes[side_set(f,1)] != num_nodes)
						ExceptionT::BadInputValue(caller, "faces side set \"%s\" have different numbers of nodes",
							ss_ID.Pointer());

				/* read in fiber orientation vectors*/
				dArray2DT& p_vec = values;
				int num_fibers = fibers.NumLists("DoubleList");
				p_vec.Dimension(num_fibers, nsd);
				if (num_fibers < 2)
					ExceptionT::GeneralFail(caller, "expecting at least two fiber orientations");
				
				for (int f = 0; f < num_fibers; f++) 
				{
					const ParameterListT& P = fibers.GetList("DoubleList", f);
					int dim = P.NumLists("Double");
					if (dim != nsd)
						ExceptionT::GeneralFail(caller, "expecting orientation vector length %d not %d",
							nsd, dim);
							
					double* p = p_vec(f); 
					/* same for all face nodes */
					for (int k = 0; k < nsd; k++)
						p[k] = P.GetList("Double", k).GetParameter("value");
				}

				if (coord_sys == Traction_CardT::kCartesian)
				{	
					/*store fiber orientation vector in element list*/
					for (int j = 0; j < num_sides; j++)
					{					
						int elem = side_set(j,0);
						fFiber_list[elem] = p_vec;
					}	
				}
				else if (coord_sys == Traction_CardT::kLocal)
				{
					for (int j = 0; j < num_sides; j++)
					{
						int elem = side_set(j,0);
						/*dimension and initialize*/
						dArray2DT& P_vec = fFiber_list[elem];
						P_vec.Dimension(num_fibers, nsd);
						P_vec = 0.0;
						
						/*Rotate from loc parent coord to global cartesian coord*
						* and store fiber orientation vectors in element list */
						/* get facet local node numbers */
						int facet = side_set(j,1);
						fShapes->NodesOnFacet(facet, facet_nodes_loc);
						int nnd = facet_nodes_loc.Length();
						coord_man.SetNumberOfNodes(nnd);

						facet_nodes_glob.Dimension(nnd);
						facet_nodes_glob.Collect(facet_nodes_loc,fElementCards[elem].NodesX());
						/*get global coordinates of facet nodes*/
						coords.SetLocal(facet_nodes_glob);
						
						/* surface shape functions */
						const ParentDomainT& surf_shape = ShapeFunction().FacetShapeFunction(facet);
						int nip = surf_shape.NumIP();
						double scale = 1.0/nip;
						for (int k = 0; k < num_fibers; k++)
						{						
							const double* p_loc = p_vec(k);
							double* p_glb = P_vec(k);
							
							/*average over element ips*/
							for (int l = 0; l < nip; l++)
							{
								surf_shape.DomainJacobian(coords, l, jacobian);
								double detj = surf_shape.SurfaceJacobian(jacobian, Q);
								/*rotate parent coords such that xi aligns with global x*/
								int sense = (Q[8] > .01*kSmall) ? 1 : -1;

								/*default*/
//								Qbar = Q;
								if (nnd == 4 || nnd == 8 || nnd ==9) /*case quad*/
								{
									const double x3 = coords(2,0);
									const double y3 = coords(2,1);
									const double x1 = coords(0,0);
									const double y1 = coords(0,1);
									if ((x3 - x1)*sense >.01*kSmall && (y3 -y1)>0.01*kSmall)
									{
										/*Qhat = {{0,1},{-1,0}}; Qbar = Q.Qhat*/
										Qbar[0] = Q[0];
										Qbar[1] = Q[1];
										Qbar[2] = Q[2];
										Qbar[3] = Q[3];
										Qbar[4] = Q[4];
										Qbar[5] = Q[5];
										Qbar[6] = Q[6];
										Qbar[7] = Q[7];
										Qbar[8] = Q[8];										
									}
									if ((x1 - x3)*sense > 0.01*kSmall && (y3-y1) >0.01*kSmall)
									{
										/*Qhat = {{0,1},{-1,0}}; Qbar = Q.Qhat*/
										Qbar[0] = -Q[3];
										Qbar[1] = -Q[4];
										Qbar[2] = -Q[5];
										Qbar[3] = Q[0];
										Qbar[4] = Q[1];
										Qbar[5] = Q[2];
										Qbar[6] = Q[6];
										Qbar[7] = Q[7];
										Qbar[8] = Q[8];										
									}
									else if ((x1 - x3)*sense > 0.01*kSmall && (y1 - y3) > 0.01*kSmall)
									{
										/*Qhat = {{-1,0},{0,-1}}; Qbar = Q.Qhat*/
										Qbar[0] = -Q[0];
										Qbar[1] = -Q[1];
										Qbar[2] = -Q[2];
										Qbar[3] = -Q[3];
										Qbar[4] = -Q[4];
										Qbar[5] = -Q[5];
										Qbar[6] = Q[6];
										Qbar[7] = Q[7];
										Qbar[8] = Q[8];
									}
									else if ((x3 - x1)*sense > 0.01*kSmall && (y1 - y3) > 0.01*kSmall)
									{
										/*Qhat = {{0,-1},{1,0}}; Qbar = Q.Qhat*/
										Qbar[0] = Q[3];
										Qbar[1] = Q[4];
										Qbar[2] = Q[5];
										Qbar[3] = -Q[0];
										Qbar[4] = -Q[1];
										Qbar[5] = -Q[2];
										Qbar[6] = Q[6];
										Qbar[7] = Q[7];
										Qbar[8] = Q[8];
									}
								}
								Qbar.Multx(p_loc, p_glb, scale, dMatrixT::kAccumulate);
							}
						}
					}
				}
		}
		else
				ExceptionT::GeneralFail(caller, "empty side set: num sides = 0");
}

void UpLagFiberCompT::ReadAnalyticVec(const ParameterListT& fibers)
{
	const char caller[] = "UpLagFiberCompT::ReadAnalyticVec";

	int nsd = NumSD();
	const dArray2DT& coordinates = ElementSupport().InitialCoordinates();

	StringT block_ID = fibers.GetParameter("block_ID");

	/*ellipsoid parameters : ((x-Cx)/Rx)^2 + ((y-Cy)/Ry)^2 + ((z-Cz)/Rz)^2 = 1 */
	double Cx = fibers.GetParameter("Cx");
	double Cy = fibers.GetParameter("Cy");
	double Cz = fibers.GetParameter("Cz");
	dArrayT C(nsd);
	C[0] = Cx; C[1] = Cy; C[2] = Cz;
	double Rx = fibers.GetParameter("Rx");
	double Ry = fibers.GetParameter("Ry");
	double Rz = fibers.GetParameter("Rz");
	dArrayT R(nsd);
	R[0] = Rx; R[1] = Ry; R[2] = Rz;
	/* projection plane */
	int inormal = fibers.GetParameter("projection_normal");
 	int i1 =0, i2 = 1, i3 = 2;
	if      (inormal ==1) { i1 =1; i2 = 2; i3 = 0;}
	else if (inormal ==2) { i1 =2; i2 = 0; i3 = 1;}
	/* cartesian or polar */
	int coor_sys = fibers.GetParameter("coordinate_system");

	/* read fiber orientation vectors*/
	dArray2DT p_vec;
	int num_fibers = fibers.NumLists("DoubleList");
	p_vec.Dimension(num_fibers, nsd-1);
	if (num_fibers < 2)
		ExceptionT::GeneralFail(caller, "expecting at least two fiber orientations");
	for (int f = 0; f < num_fibers; f++) 
	{
		const ParameterListT& P = fibers.GetList("DoubleList", f);
		int dim = P.NumLists("Double");
		if (dim != nsd-1)
			ExceptionT::GeneralFail(caller, "expecting orientation vector length %d not %d",
			nsd, dim);
							
		double* p = p_vec(f); 
		/* same for all face nodes */
		for (int k = 0; k < dim; k++)
						p[k] = P.GetList("Double", k).GetParameter("value");
	}

	/* loop over elements in the block */
	dArrayT xc(nsd),n(nsd),pn(nsd-1),pp(nsd-1);
	int block_dex = 0;
	int block_count = 0;
	const ElementBlockDataT* block_data = fBlockData.Pointer(block_dex);
	Top();
	while (NextElement())
	{
    /* reset block info (skip empty) */
    while (block_count == block_data->Dimension()) {
      block_data = fBlockData.Pointer(++block_dex);
      block_count = 0;
    }
    block_count++;

		if (block_ID == block_data->ID() || block_ID == "all" ) {

			/*dimension and initialize*/
			dArray2DT& P_vec = fFiber_list[CurrElementNumber()];
			P_vec.Dimension(num_fibers, nsd);
			P_vec = 0.0;

			/* calculate element centroid */
			iArrayT nodes = CurrentElement().NodesX();
			int nen = NumElementNodes();
			xc = 0.0;
			for (int i = 0; i < nen; i++)
			{
				for (int j = 0; j < coordinates.MinorDim(); j++)
					xc [j] += coordinates(nodes[i],j);
			}
			xc /= nen;

			/* position to normal map */
			for (int i = 0; i < nsd; i++) n[i] = (xc[i]-C[i])/(R[i]*R[i]);
			n /= n.Magnitude();

			if (coor_sys == UpLagFiberCompT::kPolar) 
			{
				/* projected normal = e_r */
				pn[0] = (xc[i1]-C[i1])/(R[i1]*R[i1]);
				pn[1] = (xc[i2]-C[i2])/(R[i2]*R[i2]);
				pn /= pn.Magnitude();
			}

			/* project s.t. in-plane direction unchanged and vector is unit */
			for (int k = 0; k < num_fibers; k++)
			{
				const double* p = p_vec(k); /* in plane direction */
				dArrayT q;  q.Alias(nsd,P_vec(k)); /* on surface direction */

				if (coor_sys == UpLagFiberCompT::kPolar) 
				{
					/* rotate in-plane vector to polar system */
					pp[0] = p[0]*pn[0] - p[1]*pn[1]; 
					pp[1] = p[0]*pn[1] + p[1]*pn[0]; 
					q[i1] =  pp[0]*n[i3];
					q[i2] =  pp[1]*n[i3];
					q[i3] = -pp[0]*n[i1]-pp[1]*n[i2];
				}
				else
				{
					q[i1] =  p[0]*n[i3];
					q[i2] =  p[1]*n[i3];
					q[i3] = -p[0]*n[i1]-p[1]*n[i2];
				}
				q /= q.Magnitude();
//#define MY_DEBUG 1
#ifdef MY_DEBUG
				cout << "block ID: " << block_data->ID() << "; ID :" << block_ID << "\n";
				cout << "q: " << q[0] << "  " << q[1] << "  " << q[2] << " " << xc[0] << "  " << xc[1] << "  " << xc[2] << "\n";
				dArrayT q2(nsd);
				if(k==0) { q2 = q; }
				else{
					cout << "# dot: " << q[0]*q2[0]+q[1]*q2[1]+q[2]*q2[2] << "\n";
				}
#endif
			}
#ifdef MY_DEBUG
			cout << "n: " << n[0] << "  " << n[1] << "  " << n[2] << " " << xc[0] << "  " << xc[1] << "  " << xc[2] << "\n";
#endif
		}
	}
}
