/* $Id: AdhesionT.cpp,v 1.21 2005-07-20 06:53:45 paklein Exp $ */
#include "AdhesionT.h"

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

/* interaction functions */
#include "LennardJones612.h"
#include "SmithFerrante.h"
#include "ModSmithFerrante.h"

using namespace Tahoe;

const int kAvgCellNodes = 10;

/* constructor */
AdhesionT::AdhesionT(const ElementSupportT& support):
	ElementBaseT(support),
	fGrid(NULL),
	fCutOff(0.0),
	fPenalizePenetration(0),
	fAllowSameSurface(0),
	fAdhesion(NULL),
	fNEE_vec_man(0, true),
	fNEE_mat_man(0, true),
	fFace2_man(0, true),
	fGrad_d_man(0, fGrad_d)
{
	SetName("adhesion");

	/* register dynamically resized arrays */
	fNEE_vec_man.Register(fRHS);
	fNEE_vec_man.Register(fNEE_vec);

	fNEE_mat_man.Register(fLHS);
	fNEE_mat_man.Register(fNEE_mat);

	fFace2_man.Register(fIPCoords2);
	fFace2_man.Register(fIPNorm2);	
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
	else {

		/* write statistics */
		ostream& out = ElementSupport().Output();
		out << "\n Surface adhesion: group " << ElementSupport().ElementGroupNumber(this) + 1 << '\n';
		out << " Time                           = " << ElementSupport().Time() << '\n';
		out << " Active face pairs              = " << fSurface1.Length() << '\n';

		return GlobalT::MaxPrecedence(relax, GlobalT::kReEQ);
	}
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
void AdhesionT::RegisterOutput(void)
{ 
	/* output labels */
	const char* displ_labels[3] = {"D_X", "D_Y", "D_Z"};
	const char* tract_labels[3] = {"t_X", "t_Y", "t_Z"};
	if (NumDOF() > 3) ExceptionT::SizeMismatch("AdhesionT::RegisterOutput", "unsupported ndof > %d", NumDOF());
	ArrayT<StringT> n_labels(2*NumDOF());
	const char** label = displ_labels;
	for (int i = 0; i < n_labels.Length(); i++)
	{
		n_labels[i] = *label++;
		if (i == NumDOF() - 1)
			label = tract_labels;
	}

	/* register each surface */
	for (int i = 0; i < fOutputID.Length(); i++)
	{
		/* set output specifier */
		OutputSetT output_set(fShapes[i]->GeometryCode(), fSurfaces[i], n_labels);
		
		/* register and get output ID */
		fOutputID[i] = ElementSupport().RegisterOutput(output_set);
	}
}

void AdhesionT::WriteOutput(void)
{
	/* statistics */
	ostream& out = ElementSupport().Output();
	out << "\n Surface adhesion: group " << ElementSupport().ElementGroupNumber(this) + 1 << '\n';
	out << " Time                           = " << ElementSupport().Time() << '\n';
	out << " Active face pairs              = " << fSurface1.Length() << '\n';
	if (fSurface1.Length() > 0 && ElementSupport().Logging() == GlobalT::kVerbose)
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
	
	if (ElementSupport().Logging() != GlobalT::kSilent) {
		out << " Search grid statistics:\n";
		fGrid->WriteStatistics(out);
		out.flush();
	}
	
	/* loop over surfaces */
	int ndof = NumDOF();
	const int nout = 2*ndof;
	for (int i = 0; i < fOutputID.Length();i++)
	{
		const iArray2DT& surface = fSurfaces[i];
		int nfn = surface.MinorDim();
	
		/* reset averaging workspace */
		ElementSupport().ResetAverage(nout);
	
		/* work space */
		dArray2DT nodal_all(nfn, nout);
		dArray2DT disp(nfn, ndof), force(nfn, ndof);

		/* local array */
		LocalArrayT loc_disp(LocalArrayT::kDisp, nfn, ndof);
		Field().RegisterLocal(loc_disp);

		/* loop over faces */
		iArrayT face_nodes, face_num(nfn);
		for (int j = 0; j < surface.MajorDim(); j++)
		{
			/* face nodes */
			surface.RowAlias(j, face_nodes);
		
			/* displacements */
			loc_disp.SetLocal(face_nodes);
			loc_disp.ReturnTranspose(disp);
			
			/* nodal tractions - same over whole face */
			double area = fCurrentFaceArea[i][j];
			if (area > kSmall) {
				face_num = j;
				force.RowCollect(face_num, fFaceForce[i]);
				force /= area;
			} else force = 0.0;

			/* copy in the cols (in sequence of output) */
			int colcount = 0;
			nodal_all.BlockColumnCopyAt(disp , colcount); colcount += disp.MinorDim();
			nodal_all.BlockColumnCopyAt(force, colcount); 

			/* accumulate */
			ElementSupport().AssembleAverage(face_nodes, nodal_all);
		}

		/* get nodally averaged values */
		dArray2DT n_values(surface.MinorDim(), nout);
		ElementSupport().OutputUsedAverage(n_values);

		/* send to output */
		dArray2DT e_values;
		ElementSupport().WriteOutput(fOutputID[i], n_values, e_values);
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

/* describe the parameters needed by the interface */
void AdhesionT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ElementBaseT::DefineParameters(list);

	ParameterT penalize_penetration(ParameterT::Boolean, "penalize_penetration");
	penalize_penetration.SetDefault(false);
	list.AddParameter(penalize_penetration);

	ParameterT self_interaction(ParameterT::Boolean, "self_interaction");
	self_interaction.SetDefault(false);
	list.AddParameter(self_interaction);

	ParameterT cut_off(ParameterT::Double, "cut_off");
	cut_off.AddLimit(kSmall, LimitT::Lower);
	list.AddParameter(cut_off);
}

/* information about subordinate parameter lists */
void AdhesionT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ElementBaseT::DefineSubs(sub_list);

	/* interaction potential choice */
	sub_list.AddSub("adhesion_potential", ParameterListT::Once);

	/* surfaces */
	sub_list.AddSub("adhesion_surface", ParameterListT::OnePlus);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* AdhesionT::NewSub(const StringT& name) const
{
	/* try C1 function */
	C1FunctionT* function = C1FunctionT::New(name);
	if (function) 
		return function;

	/* other subs */
	if (name == "adhesion_potential") {
	
		ParameterContainerT* adhesion_potential = new ParameterContainerT(name);
		adhesion_potential->SetListOrder(ParameterListT::Choice);
		adhesion_potential->SetSubSource(this);
	
		/* choice of parameters */
		adhesion_potential->AddSub("Lennard-Jones_6-12");
		adhesion_potential->AddSub("Smith-Ferrante");
		adhesion_potential->AddSub("modified_Smith-Ferrante");
	
		return adhesion_potential;
	}
	else if (name == "adhesion_surface")
	{
		ParameterContainerT* adhesion_surface = new ParameterContainerT(name);
		adhesion_surface->SetListOrder(ParameterListT::Choice);
		
		/* schedule function to scale the surface interaction */
		adhesion_surface->AddParameter(ParameterT::Integer, "interaction_scaling_schedule", ParameterListT::ZeroOrOnce);
		
		/* surface from side set  */
		ParameterContainerT surface_side_set("surface_side_set");
		surface_side_set.AddParameter(ParameterT::Word, "side_set_ID");
		adhesion_surface->AddSub(surface_side_set);
	
		/* surfaces from body block boundaries */
		ParameterContainerT block_boundaries("block_boundaries");
		block_boundaries.AddSub("block_ID_list");
		adhesion_surface->AddSub(block_boundaries);
	
		return adhesion_surface;	
	}
	else /* inherited */
		return ElementBaseT::NewSub(name);
}

/* accept parameter list */
void AdhesionT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "AdhesionT::TakeParameterList";

	/* inherited */
	ElementBaseT::TakeParameterList(list);

	/* dimensioning */
	fFace2_man.SetMinorDimension(NumSD());

	/* parameters */
	fPenalizePenetration = list.GetParameter("penalize_penetration");
	fAllowSameSurface= list.GetParameter("self_interaction");
	fCutOff = list.GetParameter("cut_off");

	/* construct potential */
	const ParameterListT& adhesion_potential = list.GetListChoice(*this, "adhesion_potential");
	fAdhesion = C1FunctionT::New(adhesion_potential.Name());
	if (!fAdhesion) ExceptionT::GeneralFail(caller, "could not construct \"%s\"", adhesion_potential.Name().Pointer());
	fAdhesion->TakeParameterList(adhesion_potential);

	/* construct the adhesive surfaces */
	ExtractSurfaces(list);
	
	/* set up work space */
	SetWorkSpace();
	
	/* set initial configuration */
	SetConfiguration();	

	/* write statistics */
	ostream& out = ElementSupport().Output();
	out << "\n Surface adhesion: group " << ElementSupport().ElementGroupNumber(this) + 1 << '\n';
	out << " Time                           = " << ElementSupport().Time() << '\n';
	out << " Active face pairs              = " << fSurface1.Length() << '\n';
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* form group contribution to the stiffness matrix */
void AdhesionT::LHSDriver(GlobalT::SystemTypeT)
{
	/* time-stepping parameters */
	double constK = 0.0;
	int     formK = fIntegrator->FormK(constK);
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
	
	     if (fAllowSameSurface || s1 != s2 ) {

		/* local face index */
		int i1 = fFaceIndex(fSurface1[i], kLocalIndex);
		int i2 = fFaceIndex(fSurface2[i], kLocalIndex);

		/* interaction scaling function */
		double scale = 0.5*(((fScaling[s1]) ? fScaling[s1]->Value(): 1.0) + 
                            ((fScaling[s2]) ? fScaling[s2]->Value(): 1.0));
		
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
			
				if (  (fPenalizePenetration || 
							fabs(d) > kSmall) &&
				    dArrayT::Dot(v_12, n1) > 0.0 &&
				    dArrayT::Dot(v_12, n2) < 0.0)
				{
					double k  = scale*j1w1*j2w2*constK;			
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
}

/* form group contribution to the residual */
void AdhesionT::RHSDriver(void)
{
	/* time-stepping parameters */
	double constKd = 0.0;
	int     formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;
	
	/* dimensions */
	int nsd = NumSD();

	/* work space */
	iArrayT nodes1, nodes2, equations;
	dArrayT ipx1, ipx2, v_12(nsd);
	dArrayT n1, n2, Na;
	dMatrixT Q1(nsd), Q2(nsd), shNaMat;
	AutoArrayT<double> j2w2_list, jump;

	/* init cached values */
	for (int i = 0; i < fSurfaces.Length(); i++)
	{
		fFaceForce[i] = 0.0;
		fCurrentFaceArea[i] = 0.0;
	}
	
	/* current area */
	double a1, a2;

	/* loop over active face pairs */
	for (int i = 0; i < fSurface1.Length(); i++)
	{
	     /* surface index */
	     int s1 = fFaceIndex(fSurface1[i], kSurface);
	     int s2 = fFaceIndex(fSurface2[i], kSurface);

	     if (fAllowSameSurface || s1 != s2 ) {

		/* local face index */
		int i1 = fFaceIndex(fSurface1[i], kLocalIndex);
		int i2 = fFaceIndex(fSurface2[i], kLocalIndex);

		/* interaction scaling function */
		double scale = 0.5*(((fScaling[s1]) ? fScaling[s1]->Value(): 1.0) + 
                           ((fScaling[s2]) ? fScaling[s2]->Value(): 1.0));

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
		a1 = a2 = 0.0;
		shape1.TopIP();
		while (shape1.NextIP())
		{
			int ip1 = shape1.CurrIP();

			/* integration point coordinates */
			ipx1 = curr_shape1.IPCoords();
			double j1w1 = shape1.Jacobian()*shape1.IPWeight();
			a1 += curr_shape1.Jacobian(Q1)*shape1.IPWeight();
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
					a2 += curr_shape2.Jacobian(Q2)*shape2.IPWeight();
					Q2.CopyColumn(nsd-1, n2);
				}
			
				/* gap vector from face 1 to face 2 */
				v_12.DiffOf(ipx2, ipx1);
				double d = v_12.Magnitude();

				if ( (fPenalizePenetration ||
							fabs(d) > kSmall) &&
				    dArrayT::Dot(v_12, n1) > 0.0 &&
				    dArrayT::Dot(v_12, n2) < 0.0)
				{
					/* adhesive force */
					double dphi =-scale*j1w1*j2w2*constKd*(fAdhesion->DFunction(d));
					
					/* face 2 shape functions */
					Na.Set(shape2.NumFacetNodes(), jump.Pointer(shape1.NumFacetNodes()));
					shape2.Shapes(Na);

					/* form d_jump/du */
					shNaMat.Set(1, jump.Length(), jump.Pointer());
					fGrad_d.Expand(shNaMat, nsd, dMatrixT::kOverwrite);
				
					/* accumulate */
					fGrad_d.MultTx(v_12, fNEE_vec);
					fRHS.AddScaled(dphi/d, fNEE_vec);
					
					/* integrate force on face */
					double f =-dphi/constKd/d;
					fFaceForce[s1].AddToRowScaled(i1, f, v_12);
					fFaceForce[s2].AddToRowScaled(i2,-f, v_12);
				}
			}
		}
		
		/* store current area */
		fCurrentFaceArea[s1][i1] = a1;
		fCurrentFaceArea[s2][i2] = a2;
		
		/* assemble */
		fFaceEquations.RowAlias(i, equations);
		ElementSupport().AssembleRHS(Group(), fRHS, equations);
	     }
	}
}

/** construct the adhesive surfaces */
void AdhesionT::ExtractSurfaces(const ParameterListT& list)
{
	const char caller[] = "AdhesionT::ExtractSurfaces";

	/* get surfaces */
	int num_surfaces = list.NumLists("adhesion_surface");

	/* read surfaces */
	AutoArrayT<int> schedules;
	AutoArrayT<iArray2DT*> surfaces;
	AutoArrayT<GeometryT::CodeT> geom;
	for (int i = 0; i < num_surfaces; i++)
	{
		const ParameterListT& surface_spec = list.GetListChoice(*this, "adhesion_surface", i);
	
		ArrayT<GeometryT::CodeT> new_geom;
		ArrayT<iArray2DT> new_surfaces;
	
		if (surface_spec.Name() == "surface_side_set") 
		{
			new_geom.Dimension(1);
			new_surfaces.Dimension(1);
			InputSideSets(surface_spec, new_geom[0], new_surfaces[0]);
		}
		else if (surface_spec.Name() == "block_boundaries")
		{		
			/* may resize the surfaces array */
			InputBodyBoundary(surface_spec, new_geom, new_surfaces);
			num_surfaces += new_surfaces.Length() - 1;
		}
		else
			ExceptionT::GeneralFail(caller, "unrecognized adhesion surface \"%s\"",
				surface_spec.Name().Pointer());

		/* collect */
		geom.Append(new_geom);
		for (int j = 0; j < new_surfaces.Length(); j++)
		{
			iArray2DT* surf = new iArray2DT;
			surf->Swap(new_surfaces[j]);
			surfaces.Append(surf);
		}
		
		/* scaling schedule */
		const ParameterListT& surface = list.GetList("adhesion_surface", i);
		const ParameterT* scaling = surface.Parameter("interaction_scaling_schedule");
		if (scaling) {
			int schedule = *scaling;
			schedule--;
			schedules.Append(schedule);	
		}
		else
			schedules.Append(-1);

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
	out << " Adhesive surfaces:\n";
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
	fScaling.Dimension(fSurfaces.Length());
	fScaling = NULL;
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
		
		/* get scaling function */
		if (schedules[i] > -1) fScaling[i] = ElementSupport().Schedule(schedules[i]);
	}

	fOutputID.Dimension(fSurfaces.Length());
	fOutputID = -1;
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
	fFaceCentroids = 0.0;
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
	if (!fGrid) fGrid = new iGridManagerT(kAvgCellNodes, -1, fFaceCentroids, NULL);
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
		const AutoArrayT<iNodeT>& hits = fGrid->HitsInRegion(fFaceCentroids(i), fCutOff);		
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
				double dist = vec_ij.Magnitude();
			
				/* surfaces "facing" each other */
				if (dist < fCutOff &&
				    normals.DotRow(i,vec_ij) > 0.0 &&
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

void AdhesionT::InputSideSets(const ParameterListT& list, GeometryT::CodeT& geom, iArray2DT& facets)
{
	const char caller[] = "AdhesionT::InputSideSets";

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

void AdhesionT::InputBodyBoundary(const ParameterListT& list, ArrayT<GeometryT::CodeT>& geom,
	ArrayT<iArray2DT>& surfaces)
{
	/* collect block ID's */
	ArrayT<StringT> IDs;
	StringListT::Extract(list.GetList("block_ID_list"), IDs);

	/* get sets of facet */
	GeometryT::CodeT geometry;
	iArrayT surface_nodes;
	ElementSupport().ModelManager().SurfaceFacets(IDs, geometry, surfaces, surface_nodes);
	
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
			ExceptionT::GeneralFail("AdhesionT::NumIP", "unrecognized geometry %d", code);
	}
	return 0;
}
