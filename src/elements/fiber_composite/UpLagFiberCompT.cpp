/* $Id: UpLagFiberCompT.cpp,v 1.16 2013-02-01 19:44:45 tahoe.kziegler Exp $ */
/* created: paklein (07/03/1996) */
#include "UpLagFiberCompT.h"

#include <cmath>
#include <iostream>
#include <iomanip>

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

//#define DEBUG
/* vector functions */
const double Pi = acos(-1.0);

inline static void CrossProduct(const double* A, const double* B, double* AxB)
{ AxB[0] = A[1]*B[2] - A[2]*B[1];
  AxB[1] = A[2]*B[0] - A[0]*B[2];
  AxB[2] = A[0]*B[1] - A[1]*B[0];
};

inline static void DotProduct(const double*A,const double* B, double c)
{ 
c=A[0]*B[0]+A[1]*B[1]+A[2]*B[2];
};

using namespace Tahoe;
/* constructor */
UpLagFiberCompT::UpLagFiberCompT(const ElementSupportT& support):
	SimoQ1P0(support),
//	UpdatedLagrangianT(support),
	fFiberSupport(NULL)
{
	SetName("uplag_fiber_comp_planar");
	
	int uplag = 4;
	
	cout << "\nuplag: "<< uplag <<endl;
}

UpLagFiberCompT::~UpLagFiberCompT(void) {
	delete fFiberSupport;
	
	int uplag5 = 4;
	
	cout << "\nuplag5: "<< uplag5 <<endl;
	
}

/* information about subordinate parameter lists */
void UpLagFiberCompT::DefineSubs(SubListT& sub_list) const
{
	
	int uplag1 = 4;
	
	cout << "\nuplag1: "<< uplag1 <<endl;
	/* inherited */
	SolidElementT::DefineSubs(sub_list);	

	/* element block/material specification */
	sub_list.AddSub("fiber_comp_element_block", ParameterListT::OnePlus);
	sub_list.AddSub("fiber_orientations", ParameterListT::OnePlus);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* UpLagFiberCompT::NewSub(const StringT& name) const
{
	
	int uplag4 = 4;
	
	cout << "\nuplag4: "<< uplag4 <<endl;
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

		{
		ParameterContainerT ellipse("revolved_ellipse");

		ParameterT block_ID(fID, "block_ID");
		block_ID.SetDefault("all");
		ellipse.AddParameter(block_ID);

		LimitT lower(0.0, LimitT::Lower);
		ParameterT Rx(ParameterT::Double, "major_axis");
		Rx.AddLimit(lower);
		Rx.SetDefault(1.0);
		ellipse.AddParameter(Rx);
		ParameterT Ry(ParameterT::Double, "minor_axis");
		Ry.AddLimit(lower);
		Ry.SetDefault(1.0);
		ellipse.AddParameter(Ry);
		ParameterT Cx(ParameterT::Double, "center_x");
		Cx.SetDefault(0.0);
		ellipse.AddParameter(Cx);
		ParameterT Cy(ParameterT::Double, "center_y");
		Cy.SetDefault(0.0);
		ellipse.AddParameter(Cy);

		ParameterT phi(ParameterT::Double, "rotation_angle_degrees");
		phi.SetDefault(0.0);
		ellipse.AddParameter(phi);
		
		ellipse.SetDescription("[(x-Cx)cos(phi) - (y-Cy)sin(phi)]^2/Rx^2 + [(x-Cx)sin(phi) + (y-Cy)cos(phi)]^2/Ry^2  = 1");

		choice->AddSub(ellipse);
		}

		{
		ParameterContainerT brick("brick");
			
		ParameterT block_ID(fID, "block_ID");
		block_ID.SetDefault("all");
		brick.AddParameter(block_ID);

		brick.SetDescription("Specify normal vector N and fiber orientation vectors in DoubleList");
		ParameterT Nx(ParameterT::Double, "Nx");
		Nx.SetDefault(0.0);
		brick.AddParameter(Nx);
		ParameterT Ny(ParameterT::Double, "Ny");
		Ny.SetDefault(0.0);
		brick.AddParameter(Ny);
		ParameterT Nz(ParameterT::Double, "Nz");
		Nz.SetDefault(1.0);
		brick.AddParameter(Nz);
			
		brick.AddSub("DoubleList", ParameterListT::OnePlus); 		
			
		choice->AddSub(brick);
		}
		
		{
		ParameterContainerT datainputfile("Input_File");
		ParameterT datafile(ParameterT::String,"data_input_file_root");
		datafile.SetDescription("Contains the fiber orientations in the lab coordinate system");
		datainputfile.AddParameter(datafile);
		choice->AddSub(datainputfile);	
		}
		
	 return choice;
	}
	else /* inherited */
		return SolidElementT::NewSub(name);
	
	int uplag6 = 4;
	cout << "\nuplag6: "<< uplag6 <<endl;
	
}

/* describe the parameters needed by the interface */
void UpLagFiberCompT::DefineParameters(ParameterListT& list) const
{
	
	int uplag2 = 4;
	
	cout << "\nuplag2: "<< uplag2 <<endl;
	
	/* inherited */
	SimoQ1P0::DefineParameters(list);
//	UpdatedLagrangianT::DefineParameters(list);
}

/* accept parameter list */
void UpLagFiberCompT::TakeParameterList(const ParameterListT& list)
{
	
	int uplag3 = 4;
	
	cout << "\nuplag3: "<< uplag3 <<endl;
	/* inherited */
	SimoQ1P0::TakeParameterList(list);

	/*store fibers in element list*/
	int num_elem = NumElements();
	fFiber_list.Dimension(num_elem);
		
	ReadFiberVec(list);

	for (int i = 0; i < NumElements(); i++)
	{
//		cout << "\nelement: "<<i
//			 << "\t"<<fFiber_list[i];
		int num_fibers = fFiber_list[i].MajorDim();		/*fFiber_list contains at least one vector for fiber orientation and one vector for the out of plane normal*/
		num_fibers--;
		if (num_fibers <1)
			ExceptionT::GeneralFail("UpLagFiberCompT::TakeParameterList",
			 "Fiber orientations not specified for element %d", i);
	}
	
}

/* extract the list of material parameters */
void UpLagFiberCompT::CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const
{
	
	int uplag7 = 4;
	
	cout << "\nuplag7: "<< uplag7 <<endl;
	
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

		fFiberSupport->SetElementCards(&fElementCards);

		/* allocate */
		return new FSFiberMatListT(size, *fFiberSupport);
	}
	else
		return new FSFiberMatListT;
}

/* read in fiber orientation information information */
void UpLagFiberCompT::ReadFiberVec(const ParameterListT& list)
{
	const char caller[] = "UpLagFiberCompT::ReadFiberVec";

//	cout << "\nreading: "<<endl;
	int num_sets = list.NumLists("fiber_orientations");
	
	//cout << "num_sets:  " << num_sets << endl;

	if (num_sets > 0)
			{
		for (int i = 0; i < num_sets; i++)
		{
			const ParameterListT& fibers = list.GetListChoice(*this, "fiber_orientations",i);
			
			
			if (fibers.Name() == "side_set") 
				ReadSideSetVec(fibers);
			else if (fibers.Name() == "ellipsoid") 
				ReadAnalyticVec(fibers);
			else if (fibers.Name() == "revolved_ellipse") 
				ReadAxi(fibers);
			else if (fibers.Name() == "brick") 
				BrickVec(fibers);
			else if (fibers.Name() == "Input_File")
				Readfiberfile(fibers);
		else
				ExceptionT::GeneralFail(caller, "invalid surface specification");
		}
	}
		else
			ExceptionT::GeneralFail(caller, "no fibers defined");
}

void UpLagFiberCompT::Readfiberfile(const ParameterListT& fibers)
/*In this function fFiber_list is not used in the standard way in which it is defined
 but it works if you only have a plane with one fiber family- 
 this will be updated shortly to assure it is used in the conventional way*/
{
const char caller[] = "UpLagFiberCompT::Readfiberfile";

fUserFile = fibers.GetParameter("data_input_file_root");
fDataInput.open(fUserFile);
int nsd = NumSD();

if (fDataInput.is_open())
  {
	cout << "fiber file successfully open"<<"\n";
		/* the file should contain as many lines as elements*/
		int num_elem = NumElements();
		double fx1, fx2, fy1, fy2, fz1, fz2, mod1, mod2, mod3, c;
		double n[3] = {0.0,0.0,0.0};	// element normal;
		double f[3] = {0.0,0.0,0.0};	// fiber direction
		double g[3] = {0.0,0.0,0.0};	// second vector in the fiber plane
		int elem_num;
	  cout<< num_elem<<endl;
		for (int k=0;k<num_elem; k++) //change to while file is still good;
		{
			fDataInput >> elem_num;
			elem_num--;
			//cout << "elem number"<< elem_num<<"\n";
			dArray2DT& P_vec = fFiber_list[elem_num];
			P_vec.Dimension(3,nsd);
			/*fiber plane defined by fx1 and fx2*/
			fDataInput >> fx1;
			fDataInput >> fy1;
			fDataInput >> fz1;
			fDataInput >> fx2;
			fDataInput >> fy2;
			fDataInput >> fz2;
			// calculate the norm of the input vector (should be 1 but in case)
			mod1 = sqrt(fx1*fx1+fy1*fy1+fz1*fz1);
			mod2 = sqrt(fx2*fx2+fy2*fy2+fz2*fz2);
			
			// normalize the input vector
			f[0]=fx1/mod1; f[1]=fy1/mod1;f[2]=fz1/mod1;
			g[0]=fx2/mod2; g[1]=fy2/mod2;g[2]=fz2/mod2;
			// need to check if the vector are orthonormal.
			// we will make them orthonormal using a Schmidt orthonormalization process
			// f2=f2-dot(f1,f2)f1;
			DotProduct(f, g, c);
			g[0]=g[0]-c*f[0]; // g_new.f=g.f-(g.f)(f.f)=0
			g[1]=g[1]-c*f[1];
			g[2]=g[2]-c*f[2];
			// normalize the new vector g
			mod2=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);	
			g[0]=g[0]/mod2;
			g[1]=g[1]/mod2;
			g[2]=g[2]/mod2;

			// now f and g should be orthonormal		
			P_vec(0,0) = f[0];
			P_vec(0,1) = f[1];
			P_vec(0,2) = f[2];
			
			P_vec(1,0) = g[0];
			P_vec(1,1) = g[1];
			P_vec(1,2) = g[2];
			

			/*normal to fiber plane = n*/
			CrossProduct(f,g,n);
			P_vec(2,0) = n[0];
			P_vec(2,1) = n[1];
			P_vec(2,2) = n[2];
			
			}
		}

	  else
			{
			cout << "Error opening file fiber direction file \n";
	}
	//for (int i=0;i<NumElements(); i++) {
//			dArray2DT& P_vec = fFiber_list[i];
//			cout << P_vec(0,0)<<" "<<P_vec(0,1)<< " "<<P_vec(0,2)<<endl;
//			cout << P_vec(1,0)<<" "<<P_vec(1,1)<< " "<<P_vec(1,2)<<endl;
//			cout << P_vec(2,0)<<" "<<P_vec(2,1)<< " "<<P_vec(2,2)<<endl;
//		cout << "----------------------"<<endl;
//		
//	}
	 
fDataInput.close(); 
	cout<< "fiber vector file closed"<<endl;
 }


void UpLagFiberCompT::ReadSideSetVec(const ParameterListT& fibers)
{
	const char caller[] = "UpLagFiberCompT::ReadSideSetVec";

		int nsd = NumSD();
		/* model manager */
		ModelManagerT& model = ElementSupport().ModelManager();
	
		/* temp space */
		StringT block_ID;    /*block id of sideset*/
		iArray2DT localsides; /*numside x 1 (elem number) + 1 (facet #1)*/ 
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
		dMatrixT Q_avg(nsd);
		Q_avg = 0.0;
			
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
//			cout <<"\nblock_ID: "<< block_ID;
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
			if (num_fibers ==0)
					ExceptionT::GeneralFail(caller, "expecting at least one fiber");
			p_vec.Dimension(num_fibers, nsd);
			
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
	
			/*Calculate rotation matrix.  Assume all facets have the same normal*/				
			int elem = side_set(0,0);
			/*Rotate from loc parent coord to global cartesian coord*
				* and store fiber orientation vectors in element list */
			/* get facet local node numbers */
			int facet = side_set(0,1);
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
			/*average over element ips*/
			for (int l = 0; l < nip; l++)
			{
				surf_shape.DomainJacobian(coords, l, jacobian);
				double detj = surf_shape.SurfaceJacobian(jacobian, Q);
				Q_avg.AddScaled(scale, Q);
			}

			/*Assigns fiber directions based on rotation matrix and coordinate system*/
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
				if (block_ID == block_data->ID() || block_ID == "all" ) 
				{
					/*dimension and initialize*/
					dArray2DT& P_vec = fFiber_list[CurrElementNumber()];
					P_vec.Dimension(num_fibers+1, nsd);
					P_vec = 0.0;
					if (coord_sys == Traction_CardT::kCartesian)
					{
						for (int k = 0; k < num_fibers; k++) /*stores fiber vectors*/
							for (int l = 0; l< nsd; l++)
								P_vec(k,l)= p_vec(k,l);		
					}
					else if (coord_sys == Traction_CardT::kLocal)
					{
						for (int k = 0; k < num_fibers; k++)
						{						
							const double* p_loc = p_vec(k);
							double* p_glb = P_vec(k);
							Q_avg.Multx(p_loc, p_glb);  /*rotates vector to local coordinate system*/
						}	
					}
									
					for (int l = 0; l< nsd; l++)				/*stores normal vector*/
						P_vec(num_fibers,l) = Q_avg(nsd-1,l);
					
				}
			}
	}
	else
		ExceptionT::GeneralFail(caller, "empty side set: num sides = 0");
}

/*applies for coordinate basis (e_r, e_phi, e_theta) where e_r represents the surface normal and e_theta is the direction of revolution.  */
/*for the revolved surface, z=0, e_theta coincides with e_z of the lab coordinate.  The fiber plane is defined by e_\phi and e_theta*/
void UpLagFiberCompT::ReadAxi(const ParameterListT& fibers)
{
	const char caller[] = "UpLagFiberCompAxiT::ReadAxi";

	int nsd = 3;
	const dArray2DT& coordinates = ElementSupport().InitialCoordinates();

	StringT block_ID = fibers.GetParameter("block_ID");

	/*ellipsoid parameters : [(x-Cx)cos(phi) - (y-Cy)sin(phi)]^2/Rx^2 + [(x-Cx)sin(phi) + (y-Cy)cos(phi)]^2/Ry^2  = 1 */
	double Cx = fibers.GetParameter("center_x");
	double Cy = fibers.GetParameter("center_y");

	double A = fibers.GetParameter("major_axis");
	double B = fibers.GetParameter("minor_axis");
	
	double angle = fibers.GetParameter("rotation_angle_degrees");
	angle *= Pi/180.;
	
	/* loop over elements in the block */
	double n_r[3] = {0.0,0.0,0.0};	// element normal;
	double t_phi[3]= {0.0,0.0,0.0}; // tangent vector;
	double t_theta[3] ={0.0,0.0,1.0}; //tangent vector;

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

		if (block_ID == block_data->ID() || block_ID == "all" ) 
		{

			/*dimension and initialize*/
			dArray2DT& P_vec = fFiber_list[CurrElementNumber()];
			P_vec.Dimension(3, nsd);
			P_vec = 0.0;

			/* calculate element centroid */
			iArrayT nodes = CurrentElement().NodesX();
			int nen = NumElementNodes();
			double xc = 0.0;
			double yc = 0.0;
			double zc = 0.0;
			for (int i = 0; i < nen; i++)
			{
				xc += coordinates(nodes[i],0);
				yc += coordinates(nodes[i],1);
				zc += coordinates(nodes[i],2);
			}
			xc /= nen; 
			yc /= nen; 
			zc /= nen;
			
			/*calculate rotation about e_theta*/
			double theta = atan2(-zc, xc);
			/*rotate centroid back to plane z=0*/
			double xs = xc*cos(theta) - zc*sin(theta);
			double ys = yc;
//#ifdef DEBUG
			double zs = xc*sin(theta) + zc*cos(theta);
			if(zs > 1e-6)
				cout << "\nzs: "<<zs;
//#endif DEBUG			
			/*calculate rotation angle phi, polar coordinate, of point (xs, ys)*/
			double c1 = (xs-Cx)*sin(angle) + (ys-Cy)*cos(angle);
			c1 /= B;
			double c2 = (xs-Cx)*cos(angle) - (ys-Cy)*sin(angle);
			c2 /= A;
			double phi = atan2(c1, c2);
			
	
			/*calculate meridional and circumferential directions associated with centroid*/
			double dxs = B*cos(phi)*sin(angle) - A*sin(phi)*cos(angle);
			double dys = B*cos(phi)*cos(angle) + A*sin(phi)*sin(angle);
			double ds = sqrt(dxs*dxs + dys*dys);
			double sgn_xs = xs/sqrt(xs*xs);
			
			t_phi[0] = dxs/ds*cos(theta);
			t_phi[1] = dys/ds;
			t_phi[2] = -dxs/ds*sin(theta);
			
			/*calculate e_theta*/
			t_theta[0] = sin(theta);
			t_theta[1] = 0.0;
			t_theta[2] = cos(theta);
			
			/*calculate normal at element centroid*/				
			CrossProduct(t_phi,t_theta, n_r);
			
			/*fiber plane defined by p1= e_phi, p2 = e_theta*/
			P_vec(0,0) = t_phi[0];
			P_vec(0,1) = t_phi[1];
			P_vec(0,2) = t_phi[2];
			
			P_vec(1,0) = t_theta[0];
			P_vec(1,1) = t_theta[1];
			P_vec(1,2) = t_theta[2];

			/*normal to fiber plane = e_r*/
			P_vec(2,0) = n_r[0];
			P_vec(2,1) = n_r[1];
			P_vec(2,2) = n_r[2];

#ifdef DEBUG
//		cout << "\nfiber vectors: "<<P_vec;
#endif DEBUG			
		}	
	}
	
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
			nsd-1, dim);
							
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
			P_vec.Dimension(num_fibers+1, nsd);
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

//					cout << "\nelem: " << CurrElementNumber();
//					cout << "\ncentroid: "<< xc[0] << "  " << xc[1] << "  " << xc[2];
			/* project s.t. in-plane direction unchanged and vector is unit */
/*			int elem = CurrElementNumber();
			if (elem == 1791 || elem == 2333 || elem == 2674 || elem == 2659 || elem == 2612 || elem == 1771)
			{
					cout << "\nelem: "<<CurrElementNumber();
					cout << "\ncentroid: "<<xc;
					cout << "\nnormal: " << n[0] << "  " << n[1] << "  " << n[2] << "\n";
			}*/
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
				
			}
			//temp: orthogonalize the system:
			//q2 = cross(n,q1)
			const double* p = p_vec(1); /* in plane direction */
			dArrayT q2;  q2.Alias(nsd,P_vec(1)); /* on surface direction */
			dArrayT q1;  q1.Alias(nsd,P_vec(0)); /* on surface direction */
			q2[i1] = n[i2]*q1[i3]-n[i3]*q1[i2];
			q2[i2] = n[i3]*q1[i1]-n[i1]*q1[i3];
			q2[i3] = n[i1]*q1[i2]-n[i2]*q1[i1];
/*
				if (elem == 1791 || elem == 2333 || elem == 2674 || elem == 2659 || elem == 2612 || elem == 1771)
				{
						cout << "\nfiber " << 1;
						cout << "\nq1: "<<q1;
						cout << "\nq2: "<<q2;
				}
*/			
			for (int l=0;l<nsd;l++)
				P_vec(num_fibers,l) = n[l];  /*stores surface normal*/

		}
	}
}


void UpLagFiberCompT::BrickVec(const ParameterListT& fibers)
{
	const char caller[] = "UpLagFiberCompT::ReadAnalyticVec";

	int nsd = NumSD();
	
	cout << "nsd" << nsd <<endl;

	
	const dArray2DT& coordinates = ElementSupport().InitialCoordinates();
	StringT block_ID = fibers.GetParameter("block_ID");

	cout << "block_ID" << block_ID <<endl;

	
	/* read fiber orientation vectors*/
	dArray2DT p_vec;
	
	cout << "p_vec" << p_vec <<endl;

	
	int num_fibers = fibers.NumLists("DoubleList");
	
	cout << "num_fibers" << num_fibers <<endl;

	
	p_vec.Dimension(num_fibers, nsd);
	
	cout << "p_vec2" << p_vec <<endl;

	
	for (int f = 0; f < num_fibers; f++) 
	{
		const ParameterListT& P = fibers.GetList("DoubleList", f);
		
		int dim = P.NumLists("Double");
		
		cout << "dim" << dim <<endl;

		if (dim != nsd)
			ExceptionT::GeneralFail(caller, "expecting orientation vector length %d not %d",
			nsd, dim);
							
		double* p = p_vec(f);
		
		cout << "f:  " << f <<endl;
		cout << "p:  " << p[0] <<endl;
		cout << "p1:  " << p[1] <<endl;
		cout << "p2:  " << p[2] <<endl;

		
		/* same for all face nodes */
		for (int k = 0; k < dim; k++)
						p[k] = P.GetList("Double", k).GetParameter("value");
		
		cout << "f:  " << f <<endl;
		cout << "p:  " << p[0] <<endl;
		cout << "p1:  " << p[1] <<endl;
		cout << "p2:  " << p[2] <<endl;
	}

	cout << "pvec_new:  " << p_vec <<endl;

	
	/* loop over elements in the block */
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
		
		cout << "block_data:  " << block_data <<endl;
		
    }
    block_count++;
		
		cout << "block_count:  " << block_count <<endl;


		if (block_ID == block_data->ID() || block_ID == "all" ) {
			
			cout << "block_ID" << block_ID <<endl;

			cout << "fFiber_List1" << fFiber_list[0] <<endl;


			/*dimension and initialize*/
			dArray2DT& P_vec = fFiber_list[CurrElementNumber()];
						
			P_vec.Dimension(num_fibers+1, nsd);
			P_vec = 0.0;
			
			cout << "fFiber_List2" << fFiber_list[0] <<endl;


			/* project s.t. in-plane direction unchanged and vector is unit */
			for (int k = 0; k < num_fibers; k++)
			{
				
				cout << "k:  " << k <<endl;

				
				P_vec(k,0) = p_vec(k,0);		
				P_vec(k,1) = p_vec(k,1);		
				P_vec(k,2) = p_vec(k,2);


			}
			
			P_vec(num_fibers,0) = fibers.GetParameter("Nx");
			P_vec(num_fibers,1) = fibers.GetParameter("Ny");
			P_vec(num_fibers,2) = fibers.GetParameter("Nz");
			
			cout << "P_vector:  " << P_vec <<endl;

			
					
			cout << "P1:  " << P_vec(0,0) <<endl;
			cout << "P2:  " << P_vec(0,1) <<endl;
			cout << "P3:  " << P_vec(0,2) <<endl;
			cout << "P4:  " << P_vec(1,0) <<endl;
			cout << "P5:  " << P_vec(1,1) <<endl;
			cout << "P6:  " << P_vec(1,2) <<endl;
			cout << "P7:  " << P_vec(2,0) <<endl;
			cout << "P8:  " << P_vec(2,1) <<endl;
			cout << "P9:  " << P_vec(2,2) <<endl;
			
			cout << "fFiber_List3" << fFiber_list[0] <<endl;

			
			
		}
	}
	
}

