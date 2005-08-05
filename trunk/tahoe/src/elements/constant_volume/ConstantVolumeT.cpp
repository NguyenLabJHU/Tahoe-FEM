/* $Id: ConstantVolumeT.cpp,v 1.1 2005-08-05 01:00:10 alindblad Exp $ */
#include "ConstantVolumeT.h"

#include "ModelManagerT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"
#include "SurfaceShapeT.h"
#include "iArrayT.h"
#include "iNodeT.h"
#include "eIntegratorT.h"
#include "iGridManagerT.h"
#include "OutputSetT.h"
#include "ScheduleT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "C1FunctionT.h"

using namespace Tahoe;

/* constructor */
ConstantVolumeT::ConstantVolumeT(const ElementSupportT& support):
	ElementBaseT(support),
	fGrid(NULL),
	fNEE_vec_man(0, true),
	fNEE_mat_man(0, true),
	fFace2_man(0, true),
	fGrad_d_man(0, fGrad_d)
{
	SetName("constant_volume");

	/* register dynamically resized arrays */
	fNEE_vec_man.Register(fRHS);
	fNEE_vec_man.Register(fNEE_vec);

	fNEE_mat_man.Register(fLHS);
	fNEE_mat_man.Register(fNEE_mat);

	fFace2_man.Register(fIPCoords2);
	fFace2_man.Register(fIPNorm2);	
}

/* destructor */
ConstantVolumeT::~ConstantVolumeT(void)
{
	for (int i = 0; i < fShapes.Length(); i++)
		delete fShapes[i];
	for (int i = 0; i < fCurrShapes.Length(); i++)
		delete fCurrShapes[i];
	delete fGrid;
}

/* element level reconfiguration for the current solution */
GlobalT::RelaxCodeT ConstantVolumeT::RelaxSystem(void)
{
    WriteCallLocation("RelaxSystem");
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

	/* generate interaction element data */
	bool changed = SetConfiguration();

	/* minimal test of new-ness */
	if (!changed)
		return relax;
	else {

		/* write statistics */
		ostream& out = ElementSupport().Output();
		out << "\n Constant volume: group " << ElementSupport().ElementGroupNumber(this) + 1 << '\n';
		out << " Time                           = " << ElementSupport().Time() << '\n';
	       	//out << " Active face pairs              = " << fSurface1.Length() << '\n';

		return GlobalT::MaxPrecedence(relax, GlobalT::kReEQ);
	}
}

/* solution calls */
void ConstantVolumeT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
#pragma unused(field)
#pragma unused(node)
#pragma unused(force)
//not implemented
    WriteCallLocation("AddNodalForce");
}

/* Returns the energy as defined by the derived class types */
double ConstantVolumeT::InternalEnergy(void)
{
    WriteCallLocation("InternalEnergy");
//not implemented
	return 0.0;
}

/* writing output - nothing to write */
void ConstantVolumeT::RegisterOutput(void)
{ 
    WriteCallLocation("RegisterOutput");
  /*// 	/* output labels */
  //	const char* displ_labels[3] = {"D_X", "D_Y", "D_Z"};
  //	const char* tract_labels[3] = {"t_X", "t_Y", "t_Z"};
  //	if (NumDOF() > 3) ExceptionT::SizeMismatch("ConstantVolumeT::RegisterOutput", "unsupported ndof > %d", NumDOF());
  //	ArrayT<StringT> n_labels(2*NumDOF());
  //	const char** label = displ_labels;
  //	for (int i = 0; i < n_labels.Length(); i++)
  //	{
  //		n_labels[i] = *label++;
  //		if (i == NumDOF() - 1)
  //			label = tract_labels;
  //	}

  //	/* register each surface */
  //	for (int i = 0; i < fOutputID.Length(); i++)
  //	{
  //		/* set output specifier */
  //		OutputSetT output_set(fShapes[i]->GeometryCode(), fSurfaces[i], n_labels);
		
  //		/* register and get output ID */
  //		fOutputID[i] = ElementSupport().RegisterOutput(output_set);
  //	}*/
}

void ConstantVolumeT::WriteOutput(void)
{
  /*not implemented: nothing to output*/
    WriteCallLocation("WriteOutput");
}

/* compute specified output parameter and send for smoothing */
void ConstantVolumeT::SendOutput(int kincode)
{
#pragma unused(kincode)
//not implemented: tractions/forces
    WriteCallLocation("SendOutput");
}

/* appends group connectivities to the array */
void ConstantVolumeT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
    WriteCallLocation("ConnectsU");
	/* inherited */
	ElementBaseT::ConnectsU(connects_1, connects_2);
	
	/* append face-pair connectivities */
	connects_2.Append(&fFaceConnectivities);
}

/* returns no (NULL) geometry connectivies */
void ConstantVolumeT::ConnectsX(AutoArrayT<const iArray2DT*>& connects) const
{
#pragma unused (connects)

    WriteCallLocation("ConnectsX");
}

/* collecting element group equation numbers */
void ConstantVolumeT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
    WriteCallLocation("Equations");
	/* inherited */
	ElementBaseT::Equations(eq_1, eq_2);
	
	/* collect equations */
	Field().SetLocalEqnos(fFaceConnectivities, fFaceEquations);

	/* add equations to the array */
	eq_2.Append(&fFaceEquations);
}

/* describe the parameters needed by the interface */
void ConstantVolumeT::DefineParameters(ParameterListT& list) const
{
    WriteCallLocation("DefineParameters");
	/* inherited */
	ElementBaseT::DefineParameters(list);
}

/* information about subordinate parameter lists */
void ConstantVolumeT::DefineSubs(SubListT& sub_list) const
{
    WriteCallLocation("DefineSubs");
	/* inherited */
	ElementBaseT::DefineSubs(sub_list);

	/* surfac */
	sub_list.AddSub("volume_surface", ParameterListT::Once);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ConstantVolumeT::NewSub(const StringT& name) const
{
    WriteCallLocation("NewSub");
	/* try C1 function */
	C1FunctionT* function = C1FunctionT::New(name);
	if (function) 
		return function;

	/* other subs */
	if (name == "volume_surface")
	{
		ParameterContainerT* volume_surface = new ParameterContainerT(name);
		volume_surface->SetListOrder(ParameterListT::Choice);
				
		/* surface from side set  */
		ParameterContainerT surface_side_set("surface_side_set");
		surface_side_set.AddParameter(ParameterT::Word, "side_set_ID");
		volume_surface->AddSub(surface_side_set);	
	
		return volume_surface;	
	}
	else /* inherited */
		return ElementBaseT::NewSub(name);
}

/* accept parameter list */
void ConstantVolumeT::TakeParameterList(const ParameterListT& list)
{
    WriteCallLocation("TakeParameterList");
	const char caller[] = "ConstantVolumeT::TakeParameterList";

	/* inherited */
	ElementBaseT::TakeParameterList(list);

	/* dimensioning */
	fFace2_man.SetMinorDimension(NumSD());

	/* construct the adhesive surfaces */
	ExtractSurfaces(list);
	
	/* set up work space */
	SetWorkSpace();
	
	/* set initial configuration */
	SetConfiguration();	

	/* write statistics */
	ostream& out = ElementSupport().Output();
	out << "\n Constant volume: group " << ElementSupport().ElementGroupNumber(this) + 1 << '\n';
	out << " Time                           = " << ElementSupport().Time() << '\n';
	//out << " Active face pairs              = " << fSurface1.Length() << '\n';
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* form group contribution to the stiffness matrix */
void ConstantVolumeT::LHSDriver(GlobalT::SystemTypeT)
{
  /* NEED TO IMPLEMENT */
    WriteCallLocation("LHSDriver"); 
}

/* form group contribution to the residual */
void ConstantVolumeT::RHSDriver(void)
{
  /* NEED TO IMPLEMENT */
    WriteCallLocation("RHSDriver");
}

/** construct the adhesive surfaces */
void ConstantVolumeT::ExtractSurfaces(const ParameterListT& list)
{
        WriteCallLocation("ExtractSurfaces");
	const char caller[] = "ConstantVolumeT::ExtractSurfaces";

	/* get surfaces */
	int num_surfaces = list.NumLists("constant_volume");

	/* read surfaces */
	AutoArrayT<iArray2DT*> surfaces;
	AutoArrayT<GeometryT::CodeT> geom;
	for (int i = 0; i < num_surfaces; i++)
	{
		const ParameterListT& surface_spec = list.GetListChoice(*this, "constant_volume", i);
	
		ArrayT<GeometryT::CodeT> new_geom;
		ArrayT<iArray2DT> new_surfaces;
	
		if (surface_spec.Name() == "surface_side_set") 
		{
			new_geom.Dimension(1);
			new_surfaces.Dimension(1);
			InputSideSets(surface_spec, new_geom[0], new_surfaces[0]);
		}
		else
			ExceptionT::GeneralFail(caller, "unrecognized constant volume surface \"%s\"",
				surface_spec.Name().Pointer());

		/* collect */
		geom.Append(new_geom);
		for (int j = 0; j < new_surfaces.Length(); j++)
		{
			iArray2DT* surf = new iArray2DT;
			surf->Swap(new_surfaces[j]);
			surfaces.Append(surf);
		}

		/* next */
		num_surfaces += new_surfaces.Length() - 1;
		i += new_surfaces.Length() - 1;
	}
	
	/* grab surface data */
	fSurfaces.Dimension(num_surfaces);
	for (int i = 0; i < fSurfaces.Length(); i++)
	{
		/* facets data */
		fSurfaces[i].Swap(*surfaces[i]);

		/* delete temp space */
		delete surfaces[i];
	}

	/* output stream */
	ofstreamT& out = ElementSupport().Output();
	bool print_input = ElementSupport().PrintInput();
	
	/* echo data */
	out << " Volume surfaces:\n";
	out << setw(kIntWidth) << "surface"
	    << setw(kIntWidth) << "facets"
	    << setw(kIntWidth) << "size" << '\n';
	for (int j = 0; j < fSurfaces.Length(); j++)
	{		
	  	iArray2DT& surface = fSurfaces[j];

	  	out << setw(kIntWidth) << j+1
	  	    << setw(kIntWidth) << surface.MajorDim()
	  	    << setw(kIntWidth) << surface.MinorDim() << "\n\n";
	  	
	 	 /* set offset for output */
	 	 if (print_input) {
	 	 	surface++;
	  		surface.WriteNumbered(out);
	 	 	surface--;
	  		out << '\n';
		}	
	}
	
	/* count non-empty */
	int surface_count = 0;
	for (int j = 0; j < fSurfaces.Length(); j++)
	  	if (fSurfaces[j].MajorDim() > 0) 
	  		surface_count++;
	
	/* remove empty surfaces */
	if (surface_count != fSurfaces.Length())
	{
		out << " Found empty surfaces:\n\n";
		ArrayT<iArray2DT> tmp_surfaces(surface_count);
		surface_count = 0;
		for (int i = 0; i < fSurfaces.Length(); i++)
		{
	  		iArray2DT& surface = fSurfaces[i];
			if (surface.MajorDim() == 0)
				out << " removing surface " << i+1 << '\n';
			else
				tmp_surfaces[surface_count++].Swap(surface);
		}
		
		/* exchange */
		fSurfaces.Swap(tmp_surfaces);
	}

	/* other per surface data */
	fLocInitCoords.Dimension(fSurfaces.Length());
	fLocCurrCoords.Dimension(fSurfaces.Length());
	fShapes.Dimension(fSurfaces.Length());
	fShapes = NULL;
	fCurrShapes.Dimension(fSurfaces.Length());
	fCurrShapes = NULL;
	fFaceForce.Dimension(fSurfaces.Length());
	fCurrentFaceArea.Dimension(fSurfaces.Length());
	for (int i = 0; i < fSurfaces.Length(); i++)
	{
		/* local array for initial coordinates */
		fLocInitCoords[i].SetType(LocalArrayT::kInitCoords);
		fLocInitCoords[i].Dimension(fSurfaces[i].MinorDim(), NumSD());
		ElementSupport().RegisterCoordinates(fLocInitCoords[i]);

		/* local array for current coordinates */
		fLocCurrCoords[i].SetType(LocalArrayT::kCurrCoords);
		fLocCurrCoords[i].Dimension(fSurfaces[i].MinorDim(), NumSD());
		ElementSupport().RegisterCoordinates(fLocCurrCoords[i]);
	
		int nen = 2*fLocInitCoords[i].NumberOfNodes();
		/* surface shape functions over undeformed configuration */
		fShapes[i] = new SurfaceShapeT(geom[i], NumIP(geom[i]), 
			nen, nen/2, NumSD(), fLocInitCoords[i]);

		/* initialize */
		fShapes[i]->Initialize();

		/* surface shape functions over current configuration */
		fCurrShapes[i] = new SurfaceShapeT(*fShapes[i], fLocCurrCoords[i]), 			

		/* initialize */
		fCurrShapes[i]->Initialize();
		
		/* space to store face forces */
		fFaceForce[i].Dimension(fSurfaces[i].MajorDim(), NumDOF());
		fFaceForce[i] = 0.0;

		/* current area */
		fCurrentFaceArea[i].Dimension(fSurfaces[i].MajorDim());
		fCurrentFaceArea[i] = 0.0;
	}

	fOutputID.Dimension(fSurfaces.Length());
	fOutputID = -1;
}

void ConstantVolumeT::SetWorkSpace(void)
{
    WriteCallLocation("SetWorkSpace");
	/* count total number of faces */
	int num_faces = 0;
	for (int i = 0; i < fSurfaces.Length(); i++)
		num_faces += fSurfaces[i].MajorDim();
		
	/* set face index array */
	fFaceIndex.Dimension(num_faces, kFaceIndexDim);
	int* p = fFaceIndex.Pointer();
	for (int i = 0; i < fSurfaces.Length(); i++)
	{
		const iArray2DT& surface = fSurfaces[i];
		for (int j = 0; j < surface.MajorDim(); j++)
		{
			*p++ = i; /* surface number */		
			*p++ = j; /* local index on surface */
		}
	}
	
	/* centroid coordinates */
	fFaceCentroids.Dimension(num_faces, NumSD());
	fFaceCentroids = 0.0;
}

/* generate element data - return true if configuration has
* changed since the last call */
bool ConstantVolumeT::SetConfiguration(void)
{
    WriteCallLocation("SetConfiguration");
  //	/* set all shape functions to the first ip */
  //	for (int i = 0; i < fCurrShapes.Length(); i++)
  //		fCurrShapes[i]->SetIP(0);

  //	/* compute the face "centroids" a normal */
  //	dArray2DT normals(fFaceIndex.MajorDim(), NumSD());
  //	dMatrixT Q(NumSD());
  //	int n_index = NumSD() - 1;
  //	dArrayT normal;
  //	iArrayT face_nodes;
  //	dArrayT centroid;
  //	const ElementSupportT& support = ElementSupport();
  //	for (int i = 0; i < fFaceIndex.MajorDim(); i++)
  //	{
  //		/* current facet information */
  //		int surface_index = fFaceIndex(i, kSurface);
  //		const iArray2DT& surface = fSurfaces[surface_index];
  //		SurfaceShapeT& shape = *fCurrShapes[surface_index];
  //		LocalArrayT& coords = fLocCurrCoords[surface_index];
  //		int local_index = fFaceIndex(i, kLocalIndex);
  //		surface.RowAlias(local_index, face_nodes);
  //		fFaceCentroids.RowAlias(i, centroid);
  //		
  //		/* collect current coordinates */
  //		coords.SetLocal(face_nodes);
  //	
  //		/* compute average coordinates */
  //		coords.Average(centroid);
  //		
  //		/* local surface axes */
  //		normals.RowAlias(i, normal);		
  //		shape.Jacobian(Q);
  //		Q.CopyColumn(n_index, normal);
  //	}
  //
  //	/* reset the search grids */
  //	if (!fGrid) fGrid = new iGridManagerT(kAvgCellNodes, -1, fFaceCentroids, NULL);
  //	fGrid->Reset();
  //
  //	/* store old configuration */
  //	ArrayT<int> s1_last(fSurface1);
  //	ArrayT<int> s2_last(fSurface2);
  //	fSurface1.Dimension(0);
  //	fSurface2.Dimension(0);
  //	
  //	/* count nodes in each face pair */
  //	iArrayT node_counts(fSurface1.Length());
  //	for (int i = 0; i < node_counts.Length(); i++)
  //	{
  //		int s1 = fFaceIndex(fSurface1[i], kSurface);
  //		int s2 = fFaceIndex(fSurface2[i], kSurface);
  //	
  //		node_counts[i]  = fSurfaces[s1].MinorDim();
  //		node_counts[i] += fSurfaces[s2].MinorDim();
  //	}

  //	/* generate connectivities */
  //	iArrayT elem, surf1, surf2;
  //	fFaceConnectivities.Configure(node_counts);
  //	for (int i = 0; i < fFaceConnectivities.MajorDim(); i++)
  //	{
  //		/* surface index */
  //		int s1 = fFaceIndex(fSurface1[i], kSurface);
  //		int s2 = fFaceIndex(fSurface2[i], kSurface);
  //	
  //		/* local face index */
  //		int i1 = fFaceIndex(fSurface1[i], kLocalIndex);
  //		int i2 = fFaceIndex(fSurface2[i], kLocalIndex);
  //	
  //		/* set aliases */
  //		fSurfaces[s1].RowAlias(i1, surf1);
  //		fSurfaces[s2].RowAlias(i2, surf2);
  //		fFaceConnectivities.RowAlias(i, elem);
  //		
  //		/* copy in */
  //		elem.CopyIn(0, surf1);
  //		elem.CopyIn(surf1.Length(), surf2);
  //	}
  //	
  //	/* allocate equations array */
  //	fFaceEquations.Configure(node_counts, NumDOF());
  //	fFaceEquations = -1;

  //	/* true if changed */
  //	return (s1_last != fSurface1) || (s2_last != fSurface2);

  /* TEMP - until I figure out what this is doing */
  return true;
}

/***********************************************************************
 * Private
 ***********************************************************************/

void ConstantVolumeT::InputSideSets(const ParameterListT& list, GeometryT::CodeT& geom, iArray2DT& facets)
{
	const char caller[] = "ConstantVolumeT::InputSideSets";

	WriteCallLocation("InputSideSets");
	/* extract side set ID */
	StringT ss_ID;
	ss_ID = list.GetParameter("side_set_ID");

	/* read side set faces */
	ModelManagerT& model = ElementSupport().ModelManager();
	ArrayT<GeometryT::CodeT> facet_geom;
	iArrayT facet_nodes;
	model.SideSet(ss_ID, facet_geom, facet_nodes, facets);
	geom = facet_geom[0];
}

/* return the number of integration points to use for the given face geometry */
int ConstantVolumeT::NumIP(GeometryT::CodeT code) const
{
    WriteCallLocation("NumIP");
	switch (code) 
	{
		case GeometryT::kLine:
			return 2;
		case GeometryT::kQuadrilateral:
		case GeometryT::kTriangle:
			return 4;
		default:
			ExceptionT::GeneralFail("ConstantVolumeT::NumIP", "unrecognized geometry %d", code);
	}
	return 0;
}

/* takes a string indicating the function calling ConstantVolumeT::WriteCallLocation
** and writes the current call location to cout
** FOR DEBUGGING AND CODE ACCLIMATION PURPOSES ONLY*/
void ConstantVolumeT::WriteCallLocation( char* loc ) const {
    cout << "*****************************************************\n";
    cout << "\tInside of ConstantVolumeT::" << loc << endl;
    cout << "*****************************************************\n";
}
