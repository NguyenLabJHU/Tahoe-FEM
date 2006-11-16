/* $Id: PressureBCT.cpp,v 1.3 2006-11-16 16:29:42 r-jones Exp $ */
#include "PressureBCT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include <ctype.h>

#include "ofstreamT.h"
#include "GlobalT.h"
#include "FieldT.h"
#include "ModelManagerT.h"
#include "CommManagerT.h"

#include "ScheduleT.h"
#include "eIntegratorT.h"
#include "IOBaseT.h"
#include "OutputSetT.h"
#include "CommunicatorT.h"
#include "ElementBaseT.h"
#include "StringT.h"

#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "FieldSupportT.h"

#include "RaggedArray2DT.h"

#define LOCAL_DEBUG

// 2 DO:
// * add penalty to enforce volume constraint
// * add multipler to enforce volume constraint
// look at: tahoe/src/elements/constant_volume

using namespace Tahoe;

const double Pi = acos(-1.0);

const double permutation[3][3][3] = 
{ 0, 0, 0, // 1
  0, 0, 1,
  0,-1, 0,
  0, 0,-1, // 2
  0, 0, 0, 
  1, 0, 0,
  0, 1, 0, // 3
 -1, 0, 0,
  0, 0, 0};

const int iperm[3][3] = 
{ -1, 2, 1,
   2,-1, 0,
   1, 0,-1};

const double psign[3][3] = 
{ 0, 1,-1,
 -1, 0, 1,
  1,-1, 0};

/* vector functions */
inline static void CrossProduct(const double* A, const double* B, double* AxB)
{ AxB[0] = A[1]*B[2] - A[2]*B[1];
  AxB[1] = A[2]*B[0] - A[0]*B[2];
  AxB[2] = A[0]*B[1] - A[1]*B[0];
};

inline static double Dot(const double* A, const double* B)
{ return A[0]*B[0] + A[1]*B[1] + A[2]*B[2]; };

inline static void Vector(const double* start, const double* end, double* v)
{
	v[0] = end[0] - start[0];
	v[1] = end[1] - start[1];
	v[2] = end[2] - start[2];
};


PressureBCT::PressureBCT(void):
	fScheduleScale(1.0),
	fPenalty(1.0),
	fControlType(kPressureControl),
	fndir(2),
	fnsd(3)
{
	SetName("pressure_bc");
}

/* destructor */
PressureBCT::~PressureBCT(void) 
{ 
}

/* form of tangent matrix */
GlobalT::SystemTypeT PressureBCT::TangentType(void) const
{
	return GlobalT::kNonSymmetric;
}

// Note: the prescribed bc acts over faces of domain elements there are
//       no new connectivities
void PressureBCT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	// for volume control all faces affect each other 
	if (fControlType == kVolumeControl) {
		// equation numnbers
		if (true) {
			Field().SetLocalEqnos(fConnectivities, fEquationNumbers);
		}
		else { // multiplier
		}
		eq_1.Append(&fEquationNumbers);
	}

	/* inherited */
//	FBC_ControllerT::Equations(eq_1, eq_2);
}

void PressureBCT::Connectivities(AutoArrayT<const iArray2DT*>& connects_1,
  AutoArrayT<const RaggedArray2DT<int>*>& connects_2,
  AutoArrayT<const iArray2DT*>& equivalent_nodes) const
{
	if (fControlType == kVolumeControl) {
		connects_1.Append(&fConnectivities);
	}

	/* inherited */
//	FBC_ControllerT::Connectivities(connects_1,connects_2,equivalent_nodes);
}

int PressureBCT::Connectivities(iArray2DT& connectivities)
{
	// connectivities : one vector of all nodes 
	int size = 0;
	if (fControlType == kVolumeControl) {
		// connectivities : one vector of all nodes 
		for (int i = 0; i < fSurfaces.Length(); i++) {
			size += fSurfaces[i].Length();
		}
		connectivities.Dimension(1,size);
		int ii=0;
		for (int i = 0; i < fSurfaces.Length(); i++) {
			const iArray2DT& faces = fSurfaces[i];
			for (int j = 0; j < faces.MajorDim(); j++)
			{		
				for (int k = 0; k < faces.MinorDim(); k++)
				{		
					connectivities(0,ii++) = faces(j,k);
				}
			}
		}
	}
	return size;
}

void PressureBCT::InitialCondition(void)
{
	// dimension 
	fReaction.Dimension(fnsd);

	const dArray2DT& Coords = FieldSupport().InitialCoordinates();
	const iArray2DT& Eqnos  = Field().Equations();
	dArray2DT coord;

	/* compute initial volume */
	double volume,area;
	fVolume0 = 0.0;
	for (int i = 0; i < fSurfaces.Length(); i++)
	{		
		iArray2DT& faces = fSurfaces[i];
		int nnodes = faces.MinorDim();
		DomainIntegrationT& domain = *(fDomains[i]);
		coord.Dimension(nnodes,fnsd);

		for (int j = 0; j < faces.MajorDim(); j++)
		{		
			const int* pface = faces(j);
			coord.RowCollect(pface, Coords);
			ComputeVolume(domain,coord,volume,area);
			fVolume0 += volume;
		}
	}
}

void PressureBCT::ReadRestart(istream& in)
{
	/* inherited */
	FBC_ControllerT::ReadRestart(in);
}

void PressureBCT::WriteRestart(ostream& out) const
{
	/* inherited */
	FBC_ControllerT::WriteRestart(out);

}

/* compute the nodal contribution to the residual force vector */
void PressureBCT::ApplyRHS(void)
{
	double constK = 0.0;
	int formK = fIntegrator->FormKd(constK);
	if (!formK) return;

	const dArray2DT& Coords = FieldSupport().CurrentCoordinates();
	const iArray2DT& Eqnos  = Field().Equations();
	dArray2DT coord;
	dArray2DT force;
	iArray2DT eqns;


	/* compute volume */
	fVolume = 0.0;
	fArea = 0.0;
	double volume, area;
	for (int i = 0; i < fSurfaces.Length(); i++)
	{		
		iArray2DT& faces = fSurfaces[i];
		int nnodes = faces.MinorDim();
		DomainIntegrationT& domain = *(fDomains[i]);
		coord.Dimension(nnodes,fnsd);
		eqns.Dimension(nnodes,fnsd);

		for (int j = 0; j < faces.MajorDim(); j++)
		{		
			const int* pface = faces(j);
			coord.RowCollect(pface, Coords);
			ComputeVolume(domain,coord,volume,area);
			fVolume += volume;
			fArea += area;
		}
	}
	if (fArea < 0) fArea = -fArea;

	// target pressure or volume
	if (fControlType == kVolumeControl) {
		double target_volume_change = fScheduleScale*(fSchedule->Value());
		target_volume_change *= constK;
		// penalty
		fPressure = fPenalty * (fVolume - fVolume0 - target_volume_change);
	}
	else {
		fPressure = fScheduleScale*(fSchedule->Value());
		fPressure *= constK;
	}

	/* compute forces */
	fReaction = 0.0;
	for (int i = 0; i < fSurfaces.Length(); i++)
	{		
		iArray2DT& faces = fSurfaces[i];
		int nnodes = faces.MinorDim();
		DomainIntegrationT& domain = *(fDomains[i]);
		coord.Dimension(nnodes,fnsd);
		force.Dimension(nnodes,fnsd);
		eqns.Dimension(nnodes,fnsd);

		for (int j = 0; j < faces.MajorDim(); j++)
		{		
			const int* pface = faces(j);
			coord.RowCollect(pface, Coords);
			ComputeForce(domain,coord,force);
			force *= fPressure;
			eqns.RowCollect(pface, Eqnos);
			FieldSupport().AssembleRHS(fGroup, force, eqns);
			for (int k = 0; k < nnodes; k++) {
				for (int l = 0; l < fnsd; l++) { fReaction[l] += force(k,l); }
			}
		}
	}
}


/* tangent term */
void PressureBCT::ApplyLHS(GlobalT::SystemTypeT sys_type)
{

	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;

	const dArray2DT& Coords = FieldSupport().CurrentCoordinates();
	const iArray2DT& Eqnos  = Field().Equations();
	dArray2DT coord;
	ElementMatrixT stiffness(ElementMatrixT::kNonSymmetric);
	iArray2DT eqns;

	/* compute stiffness */
	for (int i = 0; i < fSurfaces.Length(); i++)
	{		
	  iArray2DT& faces = fSurfaces[i];
		int nnodes = faces.MinorDim();
		DomainIntegrationT& domain = *(fDomains[i]);
		coord.Dimension(nnodes,fnsd);
		stiffness.Dimension(nnodes*fnsd);
		eqns.Dimension(nnodes,fnsd);

		for (int j = 0; j < faces.MajorDim(); j++)
		{		
			const int* pface = faces(j);
			coord.RowCollect(pface, Coords);
			ComputeStiffness(domain,coord,stiffness);
			stiffness *= -fPressure ;
			eqns.RowCollect(pface, Eqnos);
			FieldSupport().AssembleLHS(fGroup, stiffness, eqns);
		}
	}

	if (fControlType == kVolumeControl) {
	iArray2DT eqns2;
	dArray2DT force;
	dArray2DT delV;
	for (int i = 0; i < fSurfaces.Length(); i++)
	{		
	  iArray2DT& faces = fSurfaces[i];
		int nnodes = faces.MinorDim();
		DomainIntegrationT& domain = *(fDomains[i]);
		coord.Dimension(nnodes,fnsd);
		force.Dimension(nnodes,fnsd);
		delV.Dimension(nnodes,fnsd);
		stiffness.Dimension(nnodes*fnsd);
		eqns.Dimension(nnodes,fnsd);
		eqns2.Dimension(nnodes,fnsd);

		for (int j = 0; j < faces.MajorDim(); j++)
		{		
			const int* pface = faces(j);
			eqns.RowCollect(pface, Eqnos);
			coord.RowCollect(pface, Coords);
			ComputeForce(domain,coord,force);
			force *= -fPenalty;
			for (int ii = 0; ii < fSurfaces.Length(); ii++)
			{
	  		iArray2DT& faces2 = fSurfaces[ii];
				for (int jj = 0; jj < faces2.MajorDim(); jj++)
				{		
					const int* pface = faces2(jj);
					eqns2.RowCollect(pface, Eqnos);
					coord.RowCollect(pface, Coords);
					ComputeVolumeStiffness(domain,coord,delV);
					stiffness.Outer(force,delV);
					FieldSupport().AssembleLHS(fGroup, stiffness, eqns, eqns2);
				}
			}
		}
	}
	}
}

/* apply kinematic boundary conditions */
void PressureBCT::InitStep(void)
{
	/* the time step */
	//double time_step = FieldSupport().TimeStep();
}

/* finalize step */
void PressureBCT::CloseStep(void)
{
}

/* reset to the last known solution */
void PressureBCT::Reset(void)
{
}

/* update constrain forces */
GlobalT::RelaxCodeT PressureBCT::RelaxSystem(void)
{
	GlobalT::RelaxCodeT relax = FBC_ControllerT::RelaxSystem();
	
	/* re-center */
	GlobalT::RelaxCodeT my_relax = GlobalT::kNoRelax;

	/* return */
	return GlobalT::MaxPrecedence(relax, my_relax);
}

/* register data for output */
void PressureBCT::RegisterOutput(void)
{
}

/* writing results */
void PressureBCT::WriteOutput(ostream& out) const
{
	//int d_width = out.precision() + kDoubleExtra;

	/* mp support */
	//const CommunicatorT& comm = FieldSupport().Communicator();

	int step = FieldSupport().StepNumber();
	double time = FieldSupport().Time();

	out << "\n P r e s s u r e  R e g i o n   D a t a :\n\n";
	out << step << " " << time 
              << " pressure: " << fPressure 
							<< " volume: " << fVolume - fVolume0; 
#ifdef LOCAL_DEBUG
	out << " Volume: " << fVolume 
	    << " Volume0: " <<  fVolume0 ;
	out << " reaction: " << fReaction[fndir] 
	    << " area: " <<  fArea 
	    << " r/a: " << fReaction[fndir]/fArea;
#endif
	out << "\n";
	
}

/* describe the parameters needed by the interface */
void PressureBCT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FBC_ControllerT::DefineParameters(list);

	list.AddParameter(ParameterT::Integer, "schedule");

	ParameterT sscale(ParameterT::Double, "schedule_scale");
	sscale.SetDefault(1.0);
	list.AddParameter(sscale);

	ParameterT control(ParameterT::Enumeration, "control");
	control.AddEnumeration("pressure", kPressureControl);
	control.AddEnumeration("volume", kVolumeControl);
	control.SetDefault(kPressureControl);
	list.AddParameter(control);

	ParameterT penalty(ParameterT::Double, "penalty");
	penalty.SetDefault(1.0);
	penalty.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(penalty);

	ParameterT normal(ParameterT::Enumeration, "normal");
	normal.AddEnumeration("x", 0);
	normal.AddEnumeration("y", 1);
	normal.AddEnumeration("z", 2);
	normal.SetDefault(2);
	list.AddParameter(normal);

}

/* information about subordinate parameter lists */
void PressureBCT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FBC_ControllerT::DefineSubs(sub_list);

	/* surface */
	sub_list.AddSub("surfaces",ParameterListT::OnePlus);	
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* PressureBCT::NewSub(const StringT& name) const
{
	if (name == "surfaces")
	{
		ParameterContainerT* surface = new ParameterContainerT(name);
		surface->SetListOrder(ParameterListT::Choice);

		/* surface from side set  */
		ParameterContainerT surface_side_set("surface_side_set");
		surface_side_set.AddParameter(ParameterT::Word, "side_set_ID");
		surface->AddSub(surface_side_set);

		/* surfaces from body boundary */
		ParameterContainerT body_boundary("body_boundary");
		body_boundary.AddParameter(ParameterT::Integer, "body_element_group");
		surface->AddSub(body_boundary);

		return surface;
	}
	else /* inherited */
		return FBC_ControllerT::NewSub(name);
}

/* accept parameter list */
void PressureBCT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "PressureBCT::TakeParameterList";

	/* inherited */
	FBC_ControllerT::TakeParameterList(list);


	int schedule = list.GetParameter("schedule");
	schedule--;
	fSchedule = FieldSupport().Schedule(schedule);

	fScheduleScale = list.GetParameter("schedule_scale");

	fControlType =  list.GetParameter("control");

  fPenalty = list.GetParameter("penalty");

	fndir = list.GetParameter("normal");

	/* dimension */
	fnsd = FieldSupport().NumSD();
	if (fnsd != 3) 
	  ExceptionT::GeneralFail(caller, " only valid for 3D problems");


	/* get parameters */
	const ParameterListT* surfaces = list.List("surfaces");
	if (surfaces) {
		/* get surfaces */
		int num_surfaces = list.NumLists("surfaces");
		fSurfaces.Dimension(num_surfaces);
		fDomains.Dimension(num_surfaces);
		GeometryT::CodeT geometry = GeometryT::kQuadrilateral;
		for (int i = 0; i < fSurfaces.Length(); i++)
		{
			const ParameterListT& surface_spec 
				= list.GetListChoice(*this, "surfaces", i);

			if (surface_spec.Name() == "surface_side_set")
				InputSideSets(surface_spec, geometry,fSurfaces[i]);
			else if (surface_spec.Name() == "body_boundary")
			{
				ExceptionT::GeneralFail(caller, " body boundary method not implemented");
/*
				// may resize the surfaces array 
				InputBodyBoundary(surface_spec, fSurfaces, i);
				num_surfaces = fSurfaces.Length();
*/
			}
			else
				ExceptionT::GeneralFail(caller, "unrecognized surface \"%s\"",
					surface_spec.Name().Pointer());

			int nnodes = fSurfaces[i].MinorDim();
			int nip = nnodes;
			fDomains[i] = 
				new DomainIntegrationT(geometry, nip, nnodes);
			(fDomains[i])->Initialize();
		}
	}

	/* output stream */
	ofstreamT& out = FieldSupport().Output();
 	bool print_input = FieldSupport().PrintInput();

	/* echo data  */
	out << " Pressure surfaces:\n";
	out << setw(kIntWidth) << "surface"
	    << setw(kIntWidth) << "facets"
	    << setw(kIntWidth) << "size" << '\n';
	for (int j = 0; j < fSurfaces.Length(); j++)
	{
			iArray2DT& surface = fSurfaces[j];

			out << setw(kIntWidth) << j+1
			    << setw(kIntWidth) << surface.MajorDim()
			    << setw(kIntWidth) << surface.MinorDim() << "\n\n";

#ifdef LOCAL_DEBUG
		/* verbose */
		if (true) {
			surface++;
			surface.WriteNumbered(out);
			surface--;
			out << '\n';
		}
#endif
	}

	// set-up connectivities
	int size = Connectivities(fConnectivities);
	fEquationNumbers.Dimension(1,fnsd*size);
}

void PressureBCT::InputSideSets(const ParameterListT& list, 
		GeometryT::CodeT geom_code, iArray2DT& facets)
{
	/* extract side set ID */
	StringT ss_ID;
	ss_ID = list.GetParameter("side_set_ID");

	/* read side set faces */
	ModelManagerT& model = FieldSupport().ModelManager();
	ArrayT<GeometryT::CodeT> facet_geom;
	iArrayT facet_nodes;
	model.SideSet(ss_ID, facet_geom, facet_nodes, facets);
	geom_code = facet_geom[0];
	//cout << "GEOMETRY CODE: " << geom_code << "\n";
}

void PressureBCT:: ComputeVolume(DomainIntegrationT& domain, dArray2DT& coord, 
		double& volume, double& area)
{
	volume = 0.0;
	area = 0.0;
	domain.TopIP();
	const double* wgs = domain.IPWeights();
	// quadrature
	while (domain.NextIP())
	{
		/* length nnodes */
		const double* N = domain.IPShape();
		const double* T1 = domain.IPDShape(0);
		const double* T2 = domain.IPDShape(1);
		double t1[3] = {0.0,0.0,0.0};
		double t2[3] = {0.0,0.0,-1.0};
		double  n[3] = {0.0,0.0,0.0};
		for (int j = 0; j < fnsd; j++)
		{		
			t1[j] = coord.DotColumn(j,T1);
			t2[j] = coord.DotColumn(j,T2);
		}
		CrossProduct(t1,t2,n);
		double x_ndir = coord.DotColumn(fndir,N);

		const double wg = wgs[domain.CurrIP()];
		// enclosed volume has opposite outward unit normal from the surface
		volume -= n[fndir]*x_ndir;
		area   += n[fndir];
	}
}

void PressureBCT:: ComputeForce(DomainIntegrationT& domain, dArray2DT& coord, 
		dArray2DT& force)
{
	force = 0.0;
	domain.TopIP();
	const double* wgs = domain.IPWeights();
	// quadrature
	while (domain.NextIP())
	{
		/* length nnodes */
		const double* T1 = domain.IPDShape(0);
		const double* T2 = domain.IPDShape(1);
		double t1[3] = {0.0,0.0,0.0};
		double t2[3] = {0.0,0.0,-1.0};
		double  n[3] = {0.0,0.0,0.0};
		for (int j = 0; j < fnsd; j++)
		{		
			t1[j] = coord.DotColumn(j,T1);
			t2[j] = coord.DotColumn(j,T2);
		}
		CrossProduct(t1,t2,n);

		const double* S = domain.IPShape();
		const double wg = wgs[domain.CurrIP()];
		for (int k = 0; k < force.MajorDim(); k++) {
			for (int j = 0; j < force.MinorDim(); j++) {		
				force(k,j) -= S[k]*n[j]*wg; 
			}
		}
	}
}

void PressureBCT:: ComputeStiffness(DomainIntegrationT& domain, dArray2DT& coord, 
		ElementMatrixT& stiffness)
{
	stiffness = 0.0;
	domain.TopIP();
	int nnodes =  domain.ParentDomain().NumNodes();
	const double* wgs = domain.IPWeights();
	while (domain.NextIP())
	{		
		/* length nnodes */
		const double* T1 = domain.IPDShape(0);
		const double* T2 = domain.IPDShape(1);
		double t1[3] = {0.0,0.0,0.0};
		double t2[3] = {0.0,0.0,-1.0};
		for (int j = 0; j < fnsd; j++)
		{		
			t1[j] = coord.DotColumn(j,T1);
			t2[j] = coord.DotColumn(j,T2);
		}

		const double* S = domain.IPShape();
		const double wg = wgs[domain.CurrIP()];
		int row = 0, col = 0, k =0;
		for (int I = 0; I < nnodes; I++) { 
			for (int i = 0; i < fnsd; i++) {
				col = 0;
				for (int J = 0; J < nnodes; J++) {		
					for (int j = 0; j < fnsd; j++) {
						if ((k = iperm[i][j]) > -1 ) 
							stiffness(row,col) += S[I]*
								psign[i][j]*( t1[k]*T2[J] - t2[k]*T1[J] )*wg; 
						col++;
					}
				}
				row++;
			}
		}
	}
}

void PressureBCT:: ComputeVolumeStiffness(DomainIntegrationT& domain, dArray2DT& coord, 
		dArray2DT& delV)
{
	delV = 0.0;
	domain.TopIP();
	int nnodes =  domain.ParentDomain().NumNodes();
	const double* wgs = domain.IPWeights();
	// quadrature
	while (domain.NextIP())
	{
		/* length nnodes */
		const double* N = domain.IPShape();
		const double* T1 = domain.IPDShape(0);
		const double* T2 = domain.IPDShape(1);
		double t1[3] = {0.0,0.0,0.0};
		double t2[3] = {0.0,0.0,-1.0};
		double  n[3] = {0.0,0.0,0.0};
		for (int j = 0; j < fnsd; j++)
		{		
			t1[j] = coord.DotColumn(j,T1);
			t2[j] = coord.DotColumn(j,T2);
		}
		CrossProduct(t1,t2,n);

		double x_ndir = coord.DotColumn(fndir,N);

		const double* S = domain.IPShape();
		const double wg = wgs[domain.CurrIP()];
		int k;
		for (int J = 0; J < nnodes; J++) { 
			// (e_3 . n da) (e_3 . N)
			delV(J,fndir) -= n[fndir]*S[J]*wg; 
			// (e_3 . x)  (e_3 .  ee (t1 T2 - t2 T1))
			for (int j = 0; j < fnsd; j++) {
				if ((k = iperm[fndir][j]) > -1 ) 
					delV(J,j) += x_ndir *
						psign[fndir][j]*( t1[k]*T2[J] - t2[k]*T1[J] )*wg; 
			}
		}
	}
}

// functions needed for the single pressure multiplier --------------------
#if 0

void PressureBCT::SetDOFTags(void) { fMultiplierTags.Dimension(1); }

iArrayT& PressureBCT::DOFTags(int tag_set) { return fMultiplierTags; }

void PressureBCT::GenerateElementData(void)
{
	fMultiplierTags(1,2);
//  fMultiplierTags.SetColumn(0, fContactNodes);
  fMultiplierTags.SetColumn(1, fMultiplierTags);

}

const iArray2DT& PressureBCT::DOFConnects(int tag_set) const
{
	return fMultiplierConnects;
}

void PressureBCT::ResetDOF(dArray2DT& DOF, int tag_set) const
{
	dArrayT constraints;
	constraints.Alias(DOF);
//	constraints = fLastDOF;
}

int PressureBCT::Reconfigure(void) { return 0; }

int PressureBCT::Group(void) const { return Field().Group(); };
#endif

