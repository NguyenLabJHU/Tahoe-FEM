/* $Id: UpLagFiberCompT.cpp,v 1.1 2006-08-03 01:10:40 thao Exp $ */
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
	
		/* choice of materials lists (inline) */
		block->AddSub("fiber_comp_material", ParameterListT::Once);
	
		/* set this as source of subs */
		block->SetSubSource(this);
		
		return block;
	}
	/* body force */
	else if (name == "fiber_orientations") /* fiber orientations */
	{
		ParameterContainerT* fiber_orient = new ParameterContainerT(name);

		fiber_orient->AddParameter(ParameterT::Word, "side_set_ID");
		
		ParameterT coord_sys(ParameterT::Enumeration, "coordinate_system");
		coord_sys.AddEnumeration("global", Traction_CardT::kCartesian);
		coord_sys.AddEnumeration( "local", Traction_CardT::kLocal);
		coord_sys.SetDefault(Traction_CardT::kCartesian);
		fiber_orient->AddParameter(coord_sys);

		fiber_orient->AddSub("DoubleList", ParameterListT::OnePlus); 		
		
		return fiber_orient;
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

		int nsd = NumSD();
		/* model manager */
		ModelManagerT& model = ElementSupport().ModelManager();
	
		/* temp space */
		ArrayT<StringT> block_ID(num_sets);    /*block id of sideset*/
	    ArrayT<iArray2DT> localsides(num_sets); /*numside x 1 (elem number) + numfacetnodes (node numbers)*/ 
	    ArrayT<Traction_CardT::CoordSystemT> coord_sys(num_sets);  /*global cartesian or local element coords*/
	    ArrayT<dArray2DT> values(num_sets);     /*p_vec*/

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

	    /* loop over natural BC's */
	    int tot_num_sides = 0;
	    for (int i = 0; i < num_sets; i++) 
	   	{
	    	const ParameterListT& fibers = list.GetList("fiber_orientations", i);
			
	    	/* side set */
	    	const StringT& ss_ID = fibers.GetParameter("side_set_ID");
			/*reads in element and facet numbers of sideset SS_ID*/
			localsides[i] = model.SideSet(ss_ID);
			/*number of sides in set*/
			int num_sides = localsides[i].MajorDim();
			tot_num_sides += num_sides;
			if (num_sides > 0)
			{
				/*block ID of set set*/
				block_ID[i] = model.SideSetGroupID(ss_ID);
				coord_sys[i] = Traction_CardT::int2CoordSystemT(fibers.GetParameter("coordinate_system"));

				/* switch to elements numbering within the group */
				iArray2DT& side_set = localsides[i];            /*sides info for set i*/
																/* side a: elem #, facet #*/
				iArrayT elems(num_sides);
				/*copy element numbers of sideset into iArrayT elems*/
				side_set.ColumnCopy(0, elems);
				/*convert from local block numbering to global group numbering of elements*/
				BlockToGroupElementNumbers(elems, block_ID[i]);
				/*copy group element numbering in side_set array*/
				side_set.SetColumn(0, elems);

				/* all facets in set must have the same number of nodes */
				int num_nodes = num_facet_nodes[side_set(0,1)];
				for (int f = 0; f < num_sides; f++)
					if (num_facet_nodes[side_set(f,1)] != num_nodes)
						ExceptionT::BadInputValue(caller, "faces side set \"%s\" have different numbers of nodes",
							ss_ID.Pointer());

				/* read in fiber orientation vectors*/
				dArray2DT& p_vec = values[i];
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

				if (coord_sys[i] == Traction_CardT::kCartesian)
				{	
					/*store fiber orientation vector in element list*/
					for (int j = 0; j < num_sides; j++)
					{					
						int elem = side_set(j,0);
						fFiber_list[elem] = p_vec;
					}	
				}
				else if (coord_sys[i] == Traction_CardT::kLocal)
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
//								cout << "\nloc_node: \n"<<facet_nodes_loc;
//								cout << "\nglob_node: \n"<<facet_nodes_glob;
//								cout << "\nelem: "<<elem<<"\t facet: "<<facet<<"\t ip: "<<l;
//								cout <<"\nQ: "<<Q;
//								cout << "\njacobian: "<< jacobian;
								Q.Multx(p_loc, p_glb, scale, dMatrixT::kAccumulate);
							}
						}
					}
				}
			}
			else
				ExceptionT::GeneralFail(caller, "empty side set: num sides = 0");
		}
	}
}
