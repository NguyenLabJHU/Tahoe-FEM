/* $Id: AdhesionT.cpp,v 1.1.2.4 2002-10-19 17:49:39 paklein Exp $ */
#include "AdhesionT.h"

#include "ModelManagerT.h"
#include "fstreamT.h"
#include "SurfaceShapeT.h"
#include "iArrayT.h"
#include "iNodeT.h"
#include "eControllerT.h"
#include "iGridManagerT.h"

/* interaction functions */
#include "LennardJones612.h"
#include "SmithFerrante.h"

using namespace Tahoe;

const int kAvgCellNodes = 10;

/* constructor */
AdhesionT::AdhesionT(const ElementSupportT& support, const FieldT& field):
	ElementBaseT(support, field),
	fGrid(NULL),
	fCutOff(0.0),
	fAdhesion(NULL),
	fNEE_vec_man(0),
	fNEE_mat_man(0),
	fFace2_man(0, NumSD()),
	fGrad_d_man(0, fGrad_d)
{
	/* register dynamically resized arrays */
	fNEE_vec_man.Register(fRHS);
	fNEE_vec_man.Register(fNEE_vec);

	fNEE_mat_man.Register(fLHS);
	fNEE_mat_man.Register(fNEE_mat);

	fFace2_man.Register(fIPCoords2);
	fFace2_man.Register(fIPNorm2);
	
	/* read cut-off distance */
	ifstreamT& in = ElementSupport().Input();
	in >> fCutOff;
	if (fCutOff < kSmall) throw ExceptionT::kBadInputValue;

	/* set adhesion function */
	int code;
	in >> code;
	switch (code)
	{
		case C1FunctionT::kLennardJones:
		{	
			double A;
			in >> A;
			fAdhesion = new LennardJones612(A);
			break;
		}	
		case C1FunctionT::kSmithFerrante:
		{
			double A, B;
			in >> A >> B;
			fAdhesion = new SmithFerrante(A,B,0.0);
			break;
		}
		default:
			cout << "\n AdhesionT::Initialize: unrecognized function: " << code << endl;
			throw ExceptionT::kBadInputValue;	
	}
}

/* destructor */
AdhesionT::~AdhesionT(void)
{
	for (int i = 0; i < fShapes.Length(); i++)
		delete fShapes[i];
	for (int i = 0; i < fCurrShapes.Length(); i++)
		delete fCurrShapes[i];
	delete fGrid;
	delete fAdhesion;
}

/* element level reconfiguration for the current solution */
GlobalT::RelaxCodeT AdhesionT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::RelaxSystem();

	/* generate interaction element data */
	bool changed = SetConfiguration();

	/* minimal test of new-ness */
	if (!changed)
		return relax;
	else
		return GlobalT::MaxPrecedence(relax, GlobalT::kReEQ);
}

/* initialization after constructor */
void AdhesionT::Initialize(void)
{
	/* inherited */
	ElementBaseT::Initialize();

	/* set up work space */
	SetWorkSpace();
	
	/* set initial configuration */
	SetConfiguration();	
}

/* solution calls */
void AdhesionT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
#pragma unused(field)
#pragma unused(node)
#pragma unused(force)
//not implemented
}

/* Returns the energy as defined by the derived class types */
double AdhesionT::InternalEnergy(void)
{
//not implemented
	return 0.0;
}

/* writing output - nothing to write */
void AdhesionT::RegisterOutput(void) { }
void AdhesionT::WriteOutput(IOBaseT::OutputModeT mode)
{
#pragma unused(mode)

	/* statistics */
	ostream& out = ElementSupport().Output();
	out << "\n Surface adhesion: group " << ElementSupport().ElementGroupNumber(this) + 1 << '\n';
	out << " Time                           = " << ElementSupport().Time() << '\n';
	out << " Active face pairs              = " << fSurface1.Length() << '\n';
	if (fSurface1.Length() > 0)
	{
		out << setw(kIntWidth) << "surf 1";
		out << setw(kIntWidth) << "facet";
		out << setw(kIntWidth) << "surf 2";
		out << setw(kIntWidth) << "facet" << '\n';
		for (int i = 0; i < fSurface1.Length(); i++)
		{
			/* surface index */
			int s1 = fFaceIndex(fSurface1[i], kSurface);
			int s2 = fFaceIndex(fSurface2[i], kSurface);
	
			/* local face index */
			int i1 = fFaceIndex(fSurface1[i], kLocalIndex);
			int i2 = fFaceIndex(fSurface2[i], kLocalIndex);
		
			out << setw(kIntWidth) << s1+1
			    << setw(kIntWidth) << i1+1
			    << setw(kIntWidth) << s2+1
			    << setw(kIntWidth) << i2+1 << '\n';
		}
		out << endl;
	}
}

/* compute specified output parameter and send for smoothing */
void AdhesionT::SendOutput(int kincode)
{
#pragma unused(kincode)
//not implemented: tractions/forces
}

/* appends group connectivities to the array */
void AdhesionT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
	/* inherited */
	ElementBaseT::ConnectsU(connects_1, connects_2);
	
	/* append face-pair connectivities */
	connects_2.Append(&fFaceConnectivities);
}

/* returns no (NULL) geometry connectivies */
void AdhesionT::ConnectsX(AutoArrayT<const iArray2DT*>& connects) const
{
#pragma unused (connects)
}

/* collecting element group equation numbers */
void AdhesionT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	/* inherited */
	ElementBaseT::Equations(eq_1, eq_2);
	
	/* collect equations */
	Field().SetLocalEqnos(fFaceConnectivities, fFaceEquations);

	/* add equations to the array */
	eq_2.Append(&fFaceEquations);
}

/***********************************************************************
* Protected
***********************************************************************/

/* form group contribution to the stiffness matrix */
void AdhesionT::LHSDriver(void)
{
	/* time-stepping parameters */
	double constK = 0.0;
	int     formK = fController->FormK(constK);
	if (!formK) return;
	
	/* dimensions */
	int nsd = NumSD();

	/* work space */
	iArrayT nodes1, nodes2, equations;
	dArrayT ipx1, ipx2, v_12(nsd);
	dArrayT n1, n2, Na;
	dMatrixT Q1(nsd), Q2(nsd), shNaMat;
	AutoArrayT<double> j2w2_list, jump;

	/* loop over active face pairs */
	for (int i = 0; i < fSurface1.Length(); i++)
	{
		/* surface index */
		int s1 = fFaceIndex(fSurface1[i], kSurface);
		int s2 = fFaceIndex(fSurface2[i], kSurface);
	
		/* local face index */
		int i1 = fFaceIndex(fSurface1[i], kLocalIndex);
		int i2 = fFaceIndex(fSurface2[i], kLocalIndex);
		
		/* face node numbers */
		fSurfaces[s1].RowAlias(i1, nodes1);
		fSurfaces[s2].RowAlias(i2, nodes2);

		/* local coordinate arrays */
		fLocInitCoords[s1].SetLocal(nodes1);
		fLocInitCoords[s2].SetLocal(nodes2);		
		fLocCurrCoords[s1].SetLocal(nodes1);
		fLocCurrCoords[s2].SetLocal(nodes2);		
	
		/* surface shape functions */
		SurfaceShapeT& shape1 = *fShapes[s1];
		SurfaceShapeT& shape2 = *fShapes[s2];
		SurfaceShapeT& curr_shape1 = *fCurrShapes[s1];
		SurfaceShapeT& curr_shape2 = *fCurrShapes[s2];

		/* resize working arrays for the pair */
		fNEE_vec_man.Dimension(fFaceEquations.MinorDim(i), false);
		fNEE_mat_man.Dimension(fFaceEquations.MinorDim(i), fFaceEquations.MinorDim(i));
		fLHS = 0.0;
		jump.Dimension(fFaceConnectivities.MinorDim(i));
		fGrad_d_man.SetDimensions(nsd, nsd*jump.Length());

		/* resize working arrays for face 2 */
		int nip2 = shape2.NumIP();
		j2w2_list.Dimension(nip2);
		fFace2_man.SetMajorDimension(nip2, false);

		/* double-loop over integration points */
		shape1.TopIP();
		while (shape1.NextIP())
		{
			int ip1 = shape1.CurrIP();

			/* integration point coordinates */
			ipx1 = curr_shape1.IPCoords();
			double j1w1 = shape1.Jacobian()*shape1.IPWeight();
			curr_shape1.Jacobian(Q1);
			Q1.ColumnAlias(nsd-1, n1);
			
			/* face 1 shape functions */
			Na.Set(shape1.NumFacetNodes(), jump.Pointer());
			shape1.Shapes(Na);
			Na *= -1.0;

			shape2.TopIP();
			while (shape2.NextIP())
			{
				/* data for face 2 */
				int ip2 = shape2.CurrIP();
				fIPCoords2.RowAlias(ip2, ipx2);
				fIPNorm2.RowAlias(ip2, n2);
				double& j2w2 = j2w2_list[ip2];
			
				/* calculate once and store */
				if (ip1 == 0)
				{
					ipx2 = curr_shape2.IPCoords();
					j2w2 = shape2.Jacobian()*shape2.IPWeight();
					curr_shape2.Jacobian(Q2);
					Q2.CopyColumn(nsd-1, n2);
				}
			
				/* gap vector from face 1 to face 2 */
				v_12.DiffOf(ipx2, ipx1);
				double d = v_12.Magnitude();
			
				if (fabs(d) > kSmall &&
				    dArrayT::Dot(v_12, n1) > 0.0 &&
				    dArrayT::Dot(v_12, n2) < 0.0)
				{
					double k  = j1w1*j2w2*constK;			
					double k2 = k*(fAdhesion->DFunction(d))/d;
					double k1 = (k*fAdhesion->DDFunction(d) - k2)/(d*d);
					
					/* face 2 shape functions */
					Na.Set(shape2.NumFacetNodes(), jump.Pointer(shape1.NumFacetNodes()));
					shape2.Shapes(Na);
					
					/* form d_jump/du */
					shNaMat.Set(1, jump.Length(), jump.Pointer());
					fGrad_d.Expand(shNaMat, nsd, dMatrixT::kOverwrite);

					/* accumulate */
					fGrad_d.MultTx(v_12, fNEE_vec);
					fNEE_mat.Outer(fNEE_vec, fNEE_vec);
					fLHS.AddScaled(k1, fNEE_mat);
				
					fNEE_mat.MultATB(fGrad_d, fGrad_d);
					fLHS.AddScaled(k2, fNEE_mat);
				}
			}
		}
		
		/* assemble */
		fFaceEquations.RowAlias(i, equations);
		ElementSupport().AssembleLHS(Group(), fLHS, equations);
	}
}

/* form group contribution to the residual */
void AdhesionT::RHSDriver(void)
{
	/* time-stepping parameters */
	double constKd = 0.0;
	int     formKd = fController->FormKd(constKd);
	if (!formKd) return;
	
	/* dimensions */
	int nsd = NumSD();

	/* work space */
	iArrayT nodes1, nodes2, equations;
	dArrayT ipx1, ipx2, v_12(nsd);
	dArrayT n1, n2, Na;
	dMatrixT Q1(nsd), Q2(nsd), shNaMat;
	AutoArrayT<double> j2w2_list, jump;

	/* loop over active face pairs */
	for (int i = 0; i < fSurface1.Length(); i++)
	{
		/* surface index */
		int s1 = fFaceIndex(fSurface1[i], kSurface);
		int s2 = fFaceIndex(fSurface2[i], kSurface);

		/* local face index */
		int i1 = fFaceIndex(fSurface1[i], kLocalIndex);
		int i2 = fFaceIndex(fSurface2[i], kLocalIndex);

		/* face node numbers */
		fSurfaces[s1].RowAlias(i1, nodes1);
		fSurfaces[s2].RowAlias(i2, nodes2);

		/* local coordinate arrays */
		fLocInitCoords[s1].SetLocal(nodes1);
		fLocInitCoords[s2].SetLocal(nodes2);		
		fLocCurrCoords[s1].SetLocal(nodes1);
		fLocCurrCoords[s2].SetLocal(nodes2);		
	
		/* surface shape functions */
		SurfaceShapeT& shape1 = *fShapes[s1];
		SurfaceShapeT& shape2 = *fShapes[s2];
		SurfaceShapeT& curr_shape1 = *fCurrShapes[s1];
		SurfaceShapeT& curr_shape2 = *fCurrShapes[s2];

		/* resize working arrays for the pair */
		fNEE_vec_man.Dimension(fFaceEquations.MinorDim(i), false);
		fRHS = 0.0;
		jump.Dimension(fFaceConnectivities.MinorDim(i));
		fGrad_d_man.SetDimensions(nsd, nsd*jump.Length());

		/* resize working arrays for face 2 */
		int nip2 = shape2.NumIP();
		j2w2_list.Dimension(nip2);
		fFace2_man.SetMajorDimension(nip2, false);

		/* double-loop over integration points */
		shape1.TopIP();
		while (shape1.NextIP())
		{
			int ip1 = shape1.CurrIP();

			/* integration point coordinates */
			ipx1 = curr_shape1.IPCoords();
			double j1w1 = shape1.Jacobian()*shape1.IPWeight();
			curr_shape1.Jacobian(Q1);
			Q1.ColumnAlias(nsd-1, n1);
			
			/* face 1 shape functions */
			Na.Set(shape1.NumFacetNodes(), jump.Pointer());
			shape1.Shapes(Na);
			Na *= -1.0;

			shape2.TopIP();
			while (shape2.NextIP())
			{
				/* data for face 2 */
				int ip2 = shape2.CurrIP();
				fIPCoords2.RowAlias(ip2, ipx2);
				fIPNorm2.RowAlias(ip2, n2);
				double& j2w2 = j2w2_list[ip2];
			
				/* calculate once and store */
				if (ip1 == 0)
				{
					ipx2 = curr_shape2.IPCoords();
					j2w2 = shape2.Jacobian()*shape2.IPWeight();
					curr_shape2.Jacobian(Q2);
					Q2.CopyColumn(nsd-1, n2);
				}
			
				/* gap vector from face 1 to face 2 */
				v_12.DiffOf(ipx2, ipx1);
				double d = v_12.Magnitude();

				if (fabs(d) > kSmall &&
				    dArrayT::Dot(v_12, n1) > 0.0 &&
				    dArrayT::Dot(v_12, n2) < 0.0)
				{
					/* adhesive force */
					double dphi =-j1w1*j2w2*constKd*(fAdhesion->DFunction(d));
					
					/* face 2 shape functions */
					Na.Set(shape2.NumFacetNodes(), jump.Pointer(shape1.NumFacetNodes()));
					shape2.Shapes(Na);

					/* form d_jump/du */
					shNaMat.Set(1, jump.Length(), jump.Pointer());
					fGrad_d.Expand(shNaMat, nsd, dMatrixT::kOverwrite);
				
					/* accumulate */
					fGrad_d.MultTx(v_12, fNEE_vec);
					fRHS.AddScaled(dphi/d, fNEE_vec);
				}
			}
		}
		
		/* assemble */
		fFaceEquations.RowAlias(i, equations);
		ElementSupport().AssembleRHS(Group(), fRHS, equations);
	}
}

/* print element group data */
void AdhesionT::PrintControlData(ostream& out) const
{
	/* inherited */
	ElementBaseT::PrintControlData(out);

	/* adhesive force information */
	out << " Adhesive interaction. . . . . . . . . . . . . . =\n";
	fAdhesion->PrintName(out);
	out << " Parameters:\n";
	fAdhesion->Print(out);
	out << " Interaction cut-off distance. . . . . . . . . . = " << fCutOff << endl;
}

/* echo bodies and striker nodes */
void AdhesionT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
	int num_surfaces;
	in >> num_surfaces;
	out << " Number of surfaces. . . . . . . . . . . . . . . = "
	    << num_surfaces << '\n';
	if (num_surfaces < 1) throw ExceptionT::kBadInputValue;

	/* read surfaces */
	AutoArrayT<iArray2DT*> surfaces;
	AutoArrayT<GeometryT::CodeT> geom;
	for (int i = 0; i < num_surfaces; i++)
	{
		ArrayT<GeometryT::CodeT> new_geom;
		ArrayT<iArray2DT> new_surfaces;
	
		int spec_mode;
		in >> spec_mode;
		switch (spec_mode)
		{
			case kNodesOnFacet:
			{
				new_geom.Dimension(1);
				new_surfaces.Dimension(1);
				InputNodesOnFacet(in, new_geom[0], new_surfaces[0]);
				break;
			}
			case kSideSets:
			{
				new_geom.Dimension(1);
				new_surfaces.Dimension(1);
				InputSideSets(in, new_geom[0], new_surfaces[0]);
				break;
			}
			case kBodyBoundary:
			{
				/* may resize the surfaces array */
				InputBodyBoundary(in, new_geom, new_surfaces);
				num_surfaces = fSurfaces.Length();
				break;
			}
			default:
				cout << "\n AdhesionT::EchoConnectivityData: unknown surface specification\n";
				cout <<   "     mode " << spec_mode << " for surface " << i+1 << '\n';
				throw ExceptionT::kBadInputValue;
		}
		
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
	
	/* set per surface data */
	fSurfaces.Dimension(num_surfaces);
	fLocInitCoords.Dimension(fSurfaces.Length());
	fLocCurrCoords.Dimension(fSurfaces.Length());
	fShapes.Dimension(fSurfaces.Length());
	fShapes = NULL;
	fCurrShapes.Dimension(fSurfaces.Length());
	fCurrShapes = NULL;
	for (int i = 0; i < fSurfaces.Length(); i++)
	{
		/* facets data */
		fSurfaces[i].Swap(*surfaces[i]);

		/* local array for initial coordinates */
		fLocInitCoords[i].SetType(LocalArrayT::kInitCoords);
		fLocInitCoords[i].Dimension(fSurfaces[i].MinorDim(), NumSD());
		ElementSupport().RegisterCoordinates(fLocInitCoords[i]);

		/* local array for current coordinates */
		fLocCurrCoords[i].SetType(LocalArrayT::kCurrCoords);
		fLocCurrCoords[i].Dimension(fSurfaces[i].MinorDim(), NumSD());
		ElementSupport().RegisterCoordinates(fLocCurrCoords[i]);
	
		/* surface shape functions over undeformed configuration */
		fShapes[i] = new SurfaceShapeT(geom[i], NumIP(geom[i]), 
			2*fLocInitCoords[i].NumberOfNodes(), NumSD(), fLocInitCoords[i]);

		/* initialize */
		fShapes[i]->Initialize();

		/* surface shape functions over current configuration */
		fCurrShapes[i] = new SurfaceShapeT(*fShapes[i], fLocCurrCoords[i]), 			

		/* initialize */
		fCurrShapes[i]->Initialize();
	}
	
	/* delete temp space */
	for (int i = 0; i < surfaces.Length(); i++)
		delete surfaces[i];

	/* echo data and correct numbering offset */
	out << " Adhesive surfaces:\n";
	out << setw(kIntWidth) << "surface"
	    << setw(kIntWidth) << "facets"
	    << setw(kIntWidth) << "size" << '\n';
	int surface_count = 0;
	for (int j = 0; j < fSurfaces.Length(); j++)
	{		
	  	iArray2DT& surface = fSurfaces[j];

	  	out << setw(kIntWidth) << j+1
	  	    << setw(kIntWidth) << surface.MajorDim()
	  	    << setw(kIntWidth) << surface.MinorDim() << "\n\n";
	  	
	  	/* set offset for output */
	  	if (ElementSupport().PrintInput())
	  	{
	  		surface++;
	  		surface.WriteNumbered(out);
	  		surface--;
	  		out << '\n';
	  	}
	  	
	  	/* count non-empty */
	  	if (surface.MajorDim() > 0) surface_count++;
	}	
	
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
}

void AdhesionT::SetWorkSpace(void)
{
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
	
	/* construct grid */
	if (!fGrid)
		fGrid = new iGridManagerT(kAvgCellNodes, -1, fFaceCentroids, NULL);
}

/* generate element data - return true if configuration has
* changed since the last call */
bool AdhesionT::SetConfiguration(void)
{
	/* set all shape functions to the first ip */
	for (int i = 0; i < fCurrShapes.Length(); i++)
		fCurrShapes[i]->SetIP(0);

	/* compute the face "centroids" a normal */
	dArray2DT normals(fFaceIndex.MajorDim(), NumSD());
	dMatrixT Q(NumSD());
	int n_index = NumSD() - 1;
	dArrayT normal;
	iArrayT face_nodes;
	dArrayT centroid;
	const ElementSupportT& support = ElementSupport();
	for (int i = 0; i < fFaceIndex.MajorDim(); i++)
	{
		/* current facet information */
		int surface_index = fFaceIndex(i, kSurface);
		const iArray2DT& surface = fSurfaces[surface_index];
		SurfaceShapeT& shape = *fCurrShapes[surface_index];
		LocalArrayT& coords = fLocCurrCoords[surface_index];
		int local_index = fFaceIndex(i, kLocalIndex);
		surface.RowAlias(local_index, face_nodes);
		fFaceCentroids.RowAlias(i, centroid);
		
		/* collect current coordinates */
		coords.SetLocal(face_nodes);
	
		/* compute average coordinates */
		coords.Average(centroid);
		
		/* local surface axes */
		normals.RowAlias(i, normal);		
		shape.Jacobian(Q);
		Q.CopyColumn(n_index, normal);
	}
	
	/* reset the search grids */
	fGrid->Reset();

	/* store old configuration */
	ArrayT<int> s1_last(fSurface1);
	ArrayT<int> s2_last(fSurface2);
	fSurface1.Dimension(0);
	fSurface2.Dimension(0);

	/* search for interacting faces */
	dArrayT vec_ij(NumSD());
	for (int i = 0; i < fFaceIndex.MajorDim(); i++)
	{
		int i_surface = fFaceIndex(i, kSurface);
	
		/* get potential interactions */
		const AutoArrayT<iNodeT>& hits = fGrid->HitsInRegion(fFaceCentroids(i), 2.0*fCutOff);
		for (int jj = 0; jj < hits.Length(); jj++)
		{
			int j = hits[jj].Tag();
			int j_surface = fFaceIndex(j, kSurface);
			
			/* filter on index info first */
			if (i_surface < j_surface || /* different surfaces */
				(i_surface == j_surface && i < j)) /* same surface */
			{
				/* vector from centroids of i to j */
				vec_ij.DiffOf(fFaceCentroids(j), fFaceCentroids(i));
			
				/* surfaces "facing" each other */
				if (normals.DotRow(i,vec_ij) > 0.0 &&
					normals.DotRow(j,vec_ij) < 0.0)
				{
					fSurface1.Append(i);
					fSurface2.Append(j);
				}
			}
		}
	}
	
	/* count nodes in each face pair */
	iArrayT node_counts(fSurface1.Length());
	for (int i = 0; i < node_counts.Length(); i++)
	{
		int s1 = fFaceIndex(fSurface1[i], kSurface);
		int s2 = fFaceIndex(fSurface2[i], kSurface);
	
		node_counts[i]  = fSurfaces[s1].MinorDim();
		node_counts[i] += fSurfaces[s2].MinorDim();
	}

	/* generate connectivities */
	iArrayT elem, surf1, surf2;
	fFaceConnectivities.Configure(node_counts);
	for (int i = 0; i < fFaceConnectivities.MajorDim(); i++)
	{
		/* surface index */
		int s1 = fFaceIndex(fSurface1[i], kSurface);
		int s2 = fFaceIndex(fSurface2[i], kSurface);
	
		/* local face index */
		int i1 = fFaceIndex(fSurface1[i], kLocalIndex);
		int i2 = fFaceIndex(fSurface2[i], kLocalIndex);
	
		/* set aliases */
		fSurfaces[s1].RowAlias(i1, surf1);
		fSurfaces[s2].RowAlias(i2, surf2);
		fFaceConnectivities.RowAlias(i, elem);
		
		/* copy in */
		elem.CopyIn(0, surf1);
		elem.CopyIn(surf1.Length(), surf2);
	}
	
	/* allocate equations array */
	fFaceEquations.Configure(node_counts, NumDOF());
	fFaceEquations = -1;

	/* true if changed */
	return (s1_last != fSurface1) || (s2_last != fSurface2);
}

/***********************************************************************
* Private
***********************************************************************/

/* surface input functions */
void AdhesionT::InputNodesOnFacet(ifstreamT& in, GeometryT::CodeT& geom, iArray2DT& facets)
{
	int num_facets = -1;
	int num_facet_nodes = -1;
	in >> geom >> num_facets >> num_facet_nodes;
	if (num_facets < 0 || num_facet_nodes < 0) throw ExceptionT::kBadInputValue;
	
	/* dimension */
	facets.Dimension(num_facets, num_facet_nodes);
	
	/* read */
	in >> facets;
	
	/* correct numbering */
	facets--;
}

void AdhesionT::InputSideSets(ifstreamT& in, GeometryT::CodeT& geom, iArray2DT& facets)
{
	/* read data from parameter file */
	ArrayT<StringT> ss_ID;
	bool multidatabasesets = false; /* change to positive and the parameter file format changes */
	ModelManagerT& model = ElementSupport().Model();
	model.SideSetList(in, ss_ID, multidatabasesets);

	if (ss_ID.Length () != 1) {
		cout << "\n AdhesionT::InputSideSets: Model Manager read more than one side set, not programmed for this." << endl;
		throw ExceptionT::kBadInputValue;
	}

	/* read side set faces */
	ArrayT<GeometryT::CodeT> facet_geom;
	iArrayT facet_nodes;
	model.SideSet(ss_ID[0], facet_geom, facet_nodes, facets);
	geom = facet_geom[0];
}

void AdhesionT::InputBodyBoundary(ifstreamT& in, ArrayT<GeometryT::CodeT>& geom,
	ArrayT<iArray2DT>& surfaces)
{
	/* gather element group info */
	int num_blocks = -1;
	in >> num_blocks;
	if (num_blocks < 1) throw ExceptionT::kBadInputValue;
	ArrayT<StringT> IDs(num_blocks);
	for (int i = 0; i < IDs.Length(); i++)
		in >> IDs[i];

	/* get sets of facet */
	GeometryT::CodeT geometry;
	iArrayT surface_nodes;
	ElementSupport().Model().SurfaceFacets(IDs, geometry, surfaces, surface_nodes);
	
	/* face geometries */
	geom.Dimension(surfaces.Length());
	geom = geometry;
}

/* return the number of integration points to use for the given face geometry */
int AdhesionT::NumIP(GeometryT::CodeT code) const
{
	switch (code) 
	{
		case GeometryT::kLine:
			return 2;
		case GeometryT::kQuadrilateral:
		case GeometryT::kTriangle:
			return 4;
		default:
		{
			cout << "\n AdhesionT::NumIP: unrecognized geometry: " << code << endl;
			throw ExceptionT::kGeneralFail;
		}
	}
	return 0;
}
