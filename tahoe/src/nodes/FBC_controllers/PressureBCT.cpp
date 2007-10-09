/* $Id: PressureBCT.cpp,v 1.5 2007-10-09 23:17:47 r-jones Exp $ */
// created : rjones 2006
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
#include "XDOF_ManagerT.h"

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

//#define LOCAL_DEBUG

// 2 DO:

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
	fPenalty(0.0),
	fControlType(kPressureControl),
	fUseMultipliers(false),
	fndir(2),
	fnsd(3)
// NOTE: should initialize all member data
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
		if (fUseMultipliers) {
			int nuconn = fConnectivities.Length() -1;
			int ndof_u = Field().NumDOF();
			iArray2DT tmp_conn;
			tmp_conn.Alias(1,nuconn,fConnectivities.Pointer());
			iArray2DT tmp_eqnum;
			tmp_eqnum.Alias(1,nuconn*ndof_u,fEquationNumbers.Pointer());
			Field().SetLocalEqnos(tmp_conn, tmp_eqnum);
			const iArray2DT& xeqs = FieldSupport().XDOF_Manager().XDOF_Eqnos(this, 0);
			fEquationNumbers(0,nuconn*ndof_u) = xeqs(0,0);
#ifdef LOCAL_DEBUG
			cout << "x dof:  " << xeqs(0,0) << "\n";
			cout << "x node: " << fMultiplierConnects(0,0) << "\n";
			cout << "x tag:  " << fMultiplierConnects(0,1) << "\n";
#endif
		}
		else {
			Field().SetLocalEqnos(fConnectivities, fEquationNumbers);
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

void PressureBCT::SetConnectivities(void)
{
	// displacement connectivities : one vector of all nodes 
	int size = 0;
	if (fControlType == kVolumeControl) {
		// connectivities : one vector of all nodes (plus mutliplier tag)
		size = fFaces.Length();
		if (fUseMultipliers) { size++;}
		fConnectivities.Dimension(1,size);
		int ii=0;
		for (int j = 0; j < fFaces.MajorDim(); j++)
		{		
			for (int k = 0; k < fFaces.MinorDim(); k++)
			{		
				fConnectivities(0,ii++) = fFaces(j,k);
			}
		}
		if (fUseMultipliers) {fConnectivities(0,ii) = fMultiplierTags[0];}
		// equation numbers
		int ndof_u = Field().NumDOF();
		size = ndof_u*ii;
		if (fUseMultipliers) { size++;}
		fEquationNumbers.Dimension(1,size);
	}
}

void PressureBCT::InitialCondition(void)
{
	// dimension 
	fReaction.Dimension(fnsd);

	const dArray2DT& Coords = FieldSupport().InitialCoordinates();
	const iArray2DT& Eqnos  = Field().Equations();

	/* compute initial volume */
	double volume,area;
	fVolume0 = 0.0;
	for (int j = 0; j < fFaces.MajorDim(); j++)
	{		
			const int* pface = fFaces(j);
			fcoord.RowCollect(pface, Coords);
			ComputeVolume(fcoord,volume,area);
			fVolume0 += volume;
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

	/* compute volume */
	const dArray2DT& Coords = FieldSupport().CurrentCoordinates();
	fVolume = 0.0;
	fArea = 0.0;
	double volume, area;
	for (int j = 0; j < fFaces.MajorDim(); j++)
	{		
		const int* pface = fFaces(j);
		fcoord.RowCollect(pface, Coords);
		ComputeVolume(fcoord,volume,area);
		fVolume += volume;
		fArea += area;
	}
	if (fArea < 0) fArea = -fArea;

	// target pressure or volume
	if (fControlType == kVolumeControl) {
		double target_volume_change = fScheduleScale*(fSchedule->Value());
		target_volume_change *= constK;
		// penalty
		double constraint = fVolume - fVolume0 - target_volume_change;
		fPressure = fPenalty * constraint;
		// add multiplier
		if (fUseMultipliers) {
			const dArray2DT& multipliers 
				= FieldSupport().XDOF_Manager().XDOF(this, 0);
			fPressure += multipliers(0,0);
			// constraint
			const iArray2DT& xeqns 
				= FieldSupport().XDOF_Manager().XDOF_Eqnos(this, 0);
			dArray2DT xforce;
			xforce.Dimension(1,1);
			xforce(0,0) = constraint;
			FieldSupport().AssembleRHS(fGroup, xforce, xeqns);
		}
	}
	else {
		fPressure = fScheduleScale*(fSchedule->Value());
		fPressure *= constK;
	}

	/* compute forces */
	const iArray2DT& Eqnos  = Field().Equations();
	fReaction = 0.0;
	for (int j = 0; j < fFaces.MajorDim(); j++)
	{		
			const int* pface = fFaces(j);
			fcoord.RowCollect(pface, Coords);
			ComputeForce(fcoord,fforce);
			fforce *= fPressure;
			feqns.RowCollect(pface, Eqnos);
			FieldSupport().AssembleRHS(fGroup, fforce, feqns);
			// compute reactions
			for (int k = 0; k < fnnodes; k++) {
				for (int l = 0; l < fnsd; l++) { fReaction[l] += fforce(k,l); }
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
	// NOTE : global matrix looks symmetric
	ElementMatrixT stiffness(ElementMatrixT::kNonSymmetric);
	// can use SetFormat after construction
	stiffness.Dimension(fnnodes*fnsd);

	/* compute stiffness: int N p Delta( n da )*/
	for (int j = 0; j < fFaces.MajorDim(); j++)
	{		
			const int* pface = fFaces(j);
			fcoord.RowCollect(pface, Coords);
			ComputeStiffness(fcoord,stiffness);
			stiffness *= -fPressure ;
			feqns.RowCollect(pface, Eqnos);
			FieldSupport().AssembleLHS(fGroup, stiffness, feqns);
	}

	// stiffness due to volume constraint: int N pen Delta V da
	// where V = int z da
if (fControlType == kVolumeControl) {
	bool xdof_not_assembled = false;
	if ( fUseMultipliers ) {
		xdof_not_assembled = true;
	}
	for (int j = 0; j < fFaces.MajorDim(); j++)
	{
		const int* pface = fFaces(j);
		feqns.RowCollect(pface, Eqnos);
		fcoord.RowCollect(pface, Coords);
		ComputeForce(fcoord,fforce);
		if ( fUseMultipliers ) {
			fdelP = fforce;
		}
		fforce *= -fPenalty;
		for (int jj = 0; jj < fFaces.MajorDim(); jj++)
		{		
			const int* pface = fFaces(jj);
			feqns2.RowCollect(pface, Eqnos);
			fcoord.RowCollect(pface, Coords);
			ComputeVolumeStiffness(fcoord, fdelV);
			if ( fUseMultipliers &&  xdof_not_assembled  ) {
				ElementMatrixT Vstiffness(ElementMatrixT::kNonSymmetric);
				Vstiffness.Alias(1,fnnodes*fnsd,fdelV.Pointer());
				Vstiffness *= -1;
				const iArray2DT& xeqns 
					= FieldSupport().XDOF_Manager().XDOF_Eqnos(this, 0);
				FieldSupport().AssembleLHS(fGroup, Vstiffness, xeqns, feqns2);
			}
			stiffness.Outer(fforce,fdelV);
			FieldSupport().AssembleLHS(fGroup, stiffness, feqns, feqns2);
		}
		xdof_not_assembled = false; // assemble delV once
		if ( fUseMultipliers ) {
			ElementMatrixT Pstiffness(ElementMatrixT::kNonSymmetric);
			Pstiffness.Alias(fnnodes*fnsd,1,fdelP.Pointer());
			Pstiffness *= -1;
			const iArray2DT& xeqns 
				= FieldSupport().XDOF_Manager().XDOF_Eqnos(this, 0);
			FieldSupport().AssembleLHS(fGroup, Pstiffness, feqns, xeqns);
		}
	}

}
}

/* apply kinematic boundary conditions */
void PressureBCT::InitStep(void)
{
}

/* finalize step */
void PressureBCT::CloseStep(void)
{
  /* store last converged DOF array */
  if (fUseMultipliers) {
    dArrayT xdof;
    xdof.Alias(FieldSupport().XDOF_Manager().XDOF(this, 0));
    fLastDOF = xdof;
  }
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
#if 0
  ArrayT<StringT> n_labels(fnsd);
  n_labels[0] = "normal_X";
  n_labels[1] = "normal_Y";
  n_labels[2] = "normal_Z";

  OutputSetT output_set(GeometryT::kPoint, fFlattenedNodeSets, n_labels);
  fOutputID = field_support.RegisterOutput(output_set);
#endif
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
	if (fUseMultipliers) {
		const dArray2DT& multipliers 
			= FieldSupport().XDOF_Manager().XDOF(this, 0);
		out << " multiplier: " << multipliers(0,0);
	}
#endif
	out << "\n";

#if 0
  /* send output */
  dArray2DT e_values;
  field_support.WriteOutput(fOutputID, n_values, e_values);
#endif
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

	ParameterT use_multiplier(ParameterT::Enumeration, "use_multiplier");
	use_multiplier.AddEnumeration("false", 0);
	use_multiplier.AddEnumeration("true", 1);
	use_multiplier.SetDefault(0);
	list.AddParameter(use_multiplier);

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

	/* surface : a collection of side sets */
	sub_list.AddSub("side_set_ID_list");// str list must end in "_list"
}

/* accept parameter list */
void PressureBCT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "PressureBCT::TakeParameterList";

	/* inherited */
	FBC_ControllerT::TakeParameterList(list);

	/* basic parameters */
	int schedule = list.GetParameter("schedule");
	schedule--;
	fSchedule = FieldSupport().Schedule(schedule);

	fScheduleScale = list.GetParameter("schedule_scale");

	fControlType =  list.GetParameter("control");

	if (fControlType == kVolumeControl) {
	fUseMultipliers = list.GetParameter("use_multiplier");

	fPenalty = list.GetParameter("penalty");

	fndir = list.GetParameter("normal");
	}

	/* dimension */
	fnsd = FieldSupport().NumSD();
	if (fnsd != 3) 
	  ExceptionT::GeneralFail(caller, " only valid for 3D problems");

	/* surface : collection of side sets/faces */
	StringListT::Extract(list.GetList("side_set_ID_list"),  fssetIDs);
	ModelManagerT& model = FieldSupport().ModelManager();
	if (fssetIDs.Length()) {
		// size
		int nfaces = 0;
		for (int i = 0; i < fssetIDs.Length(); i++)
		{
			nfaces += model.SideSetLength(fssetIDs[i]);
		}
		// fill
		GeometryT::CodeT geometry;
		iArrayT fFaces_alias;
		int size = 0;	
		for (int i = 0; i < fssetIDs.Length(); i++)
		{
			ArrayT<GeometryT::CodeT> face_geom;
			iArrayT face_nodes;
			iArray2DT faces;
			model.SideSet(fssetIDs[i], face_geom, face_nodes, faces);
#if 0
			if (model.IsSideSetLocal(fssetIDs[i]) ) {
					ExceptionT::GeneralFail(caller, " all faces must globally numbered");
			}
// void 	SideSetLocalToGlobal (const StringT &element_ID, const iArray2DT
// &local, iArray2DT &global)
#endif
			if (i == 0) {
				geometry = face_geom[0];
				fnnodes = faces.MinorDim();
				fFaces.Dimension(nfaces, fnnodes);
				fFaces_alias.Alias(fFaces.Length(),fFaces.Pointer());
				int nip = fnnodes;
				fDomain = new DomainIntegrationT(geometry, nip, fnnodes);
				fDomain->Initialize();
			}
			else {
				if (face_geom[0] != geometry) { 
					ExceptionT::GeneralFail(caller, " all faces must be same type");
				}
			}
			// copy
			iArrayT faces_alias;
			faces_alias.Alias(faces.Length(),faces.Pointer());
			fFaces_alias.CopyIn(size,faces_alias);
			size += faces.Length();
		}
	}
	else {
		ExceptionT::GeneralFail(caller, " no surfaces defined");
	}

#ifdef LOCAL_DEBUG
	/* output stream */
	ofstreamT& out = FieldSupport().Output();
 	bool print_input = FieldSupport().PrintInput();

	/* echo data  */
	out << " Pressure surface: " << fFaces.MajorDim() << " faces\n";
	out << setw(kIntWidth) << "surface"
	    << setw(kIntWidth) << "faces"
	    << setw(kIntWidth) << "size" << '\n';
	/* verbose */
	fFaces++;
	fFaces.WriteNumbered(out);
	fFaces--;
	out << '\n';
#endif

	// attach single pressure multiplier to first node
	if (fUseMultipliers) {
		fMultiplierNodes.Dimension(1);
		fMultiplierNodes[0] = fFaces(0,0);
		/* register with node manager - sets initial fContactDOFtags */
		iArrayT set_dims(1);
		set_dims = 1; // one pressure multiplier
		FieldSupport().XDOF_Manager().XDOF_Register(this, set_dims);
	}

	// set-up connectivities
	SetConnectivities();

	// set-up workspace
	fcoord.Dimension(fnnodes,fnsd);
	fforce.Dimension(fnnodes,fnsd);
	feqns.Dimension(fnnodes,fnsd);
	if (fControlType == kVolumeControl) {
	feqns2.Dimension(fnnodes,fnsd);
	fdelV.Dimension(fnnodes,fnsd);
	if (fUseMultipliers) {
	fdelP.Dimension(fnnodes,fnsd);
	}
	}
}

void PressureBCT:: ComputeVolume(dArray2DT& coord, double& volume, double& area)
{
	DomainIntegrationT domain = *fDomain;
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
#ifdef LOCAL_DEBUG
		double  x[3] = {0.0,0.0,0.0};
		for (int j = 0; j < fnsd; j++) x[j] = coord.DotColumn(j,N);
		cout << "x: " << x[0] << " " << x[1] << " " << x[2] << " "
		     << "n: " << n[0] << " " << n[1] << " " << n[2] << "\n";
#endif
	}
}

void PressureBCT:: ComputeForce(dArray2DT& coord, dArray2DT& force)
{
	DomainIntegrationT domain = *fDomain;
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

void PressureBCT:: ComputeStiffness(dArray2DT& coord, ElementMatrixT& stiffness)
{
	DomainIntegrationT domain = *fDomain;
	stiffness = 0.0;
	domain.TopIP();
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
		for (int I = 0; I < fnnodes; I++) { 
			for (int i = 0; i < fnsd; i++) {
				col = 0;
				for (int J = 0; J < fnnodes; J++) {		
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

void PressureBCT:: ComputeVolumeStiffness(dArray2DT& coord, dArray2DT& delV)
{
	DomainIntegrationT domain = *fDomain;
	delV = 0.0;
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

		const double* S = domain.IPShape();
		const double wg = wgs[domain.CurrIP()];
		int k;
		for (int J = 0; J < fnnodes; J++) { 
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

// functions for the single pressure multiplier of the volume constraint

void PressureBCT::SetDOFTags(void) { fMultiplierTags.Dimension(1); }

iArrayT& PressureBCT::DOFTags(int tag_set) { return fMultiplierTags; }

void PressureBCT::GenerateElementData(void)
{
	fMultiplierConnects.Dimension(1,2);
  fMultiplierConnects.SetColumn(0, fMultiplierNodes);
  fMultiplierConnects.SetColumn(1, fMultiplierTags);

}

const iArray2DT& PressureBCT::DOFConnects(int tag_set) const
{ return fMultiplierConnects; }

void PressureBCT::ResetDOF(dArray2DT& DOF, int tag_set) const
{
	dArrayT constraints;
	constraints.Alias(DOF);
	constraints = fLastDOF;
}

int PressureBCT::Reconfigure(void) { return 0; }

int PressureBCT::Group(void) const { return Field().Group(); };

