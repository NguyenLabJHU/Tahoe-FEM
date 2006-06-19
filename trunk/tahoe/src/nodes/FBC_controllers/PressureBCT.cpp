/* $Id: PressureBCT.cpp,v 1.1 2006-06-19 15:25:34 r-jones Exp $ */
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

#define LOCAL_DEBUG

using namespace Tahoe;

const double Pi = acos(-1.0);

const double permutation[3][3][3] = { 0, 0, 0, // 1
	                              0, 0, 1,
	                              0,-1, 0,
                                      0, 0,-1, // 2
	                              0, 0, 0, 
	                              1, 0, 0,
                                      0, 1, 0, // 3
	                             -1, 0, 0,
	                              0, 0, 0};
				      
const int iperm[3][3] = { -1, 2, 1,
	                   2,-1, 0,
	                   1, 0,-1};

const double psign[3][3] = { 0, 1,-1,
	                    -1, 0, 1,
	                     1,-1, 0};

/* vector functions */
inline static void CrossProduct(const double* A, const double* B, double* AxB)
{   AxB[0] = A[1]*B[2] - A[2]*B[1];
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
	fpscale(1.0)
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

// Note: since this bc acts over faces of domain elements there are
//       no new connectivities
void PressureBCT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	/* inherited */
	FBC_ControllerT::Equations(eq_1, eq_2);
}

void PressureBCT::InitialCondition(void)
{
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

	int nsd = FieldSupport().NumSD();
	const dArray2DT& Coords = FieldSupport().CurrentCoordinates();
	//const dArray2DT& Coords = FieldSupport().InitialCoordinates();
	const iArray2DT& Eqnos  = Field().Equations();
	dArray2DT coord;
	dArray2DT force;
	iArray2DT eqns;

	double pval = fpscale*(fSchedule->Value());
        pval *= constK;

	/* compute forces */
	for (int i = 0; i < fSurfaces.Length(); i++)
	{		
	  	iArray2DT& faces = fSurfaces[i];
		int nnodes = faces.MinorDim();
		DomainIntegrationT& domain = *(fDomains[i]);
		coord.Dimension(nnodes,nsd);
		force.Dimension(nnodes,nsd);
		eqns.Dimension(nnodes,nsd);

		for (int j = 0; j < faces.MajorDim(); j++)
		{		
			const int* pface = faces(j);
			coord.RowCollect(pface, Coords);
			Compute_Force(domain,coord,force);
			force *= pval;
			eqns.RowCollect(pface, Eqnos);
			FieldSupport().AssembleRHS(fGroup, force, eqns);
		}
	}
}


/* tangent term */
void PressureBCT::ApplyLHS(GlobalT::SystemTypeT sys_type)
{

	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;

	int nsd = FieldSupport().NumSD();
	const dArray2DT& Coords = FieldSupport().CurrentCoordinates();
	const iArray2DT& Eqnos  = Field().Equations();
	dArray2DT coord;
	ElementMatrixT stiffness(ElementMatrixT::kNonSymmetric);
	iArray2DT eqns;

	double pval = fpscale*(fSchedule->Value());
        pval *= constK;

	/* compute stiffness */
	for (int i = 0; i < fSurfaces.Length(); i++)
	{		
	  	iArray2DT& faces = fSurfaces[i];
		int nnodes = faces.MinorDim();
		DomainIntegrationT& domain = *(fDomains[i]);
		coord.Dimension(nnodes,nsd);
		stiffness.Dimension(nnodes*nsd);
		eqns.Dimension(nnodes,nsd);

		for (int j = 0; j < faces.MajorDim(); j++)
		{		
			const int* pface = faces(j);
			coord.RowCollect(pface, Coords);
			Compute_Stiffness(domain,coord,stiffness);
			stiffness *= -pval ;
			eqns.RowCollect(pface, Eqnos);
			FieldSupport().AssembleLHS(fGroup, stiffness, eqns);
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

	out << "\n P r e s s u r e  R e g i o n   D a t a :\n\n";
	
}

/* describe the parameters needed by the interface */
void PressureBCT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FBC_ControllerT::DefineParameters(list);

	ParameterT pscale(ParameterT::Double, "pressure_scale");
	pscale.SetDefault(1.0);
	pscale.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(pscale);
    
        list.AddParameter(ParameterT::Integer, "schedule");
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

	fpscale = list.GetParameter("pressure_scale");

	int schedule = list.GetParameter("schedule");
	schedule--;
	fSchedule = FieldSupport().Schedule(schedule);


	/* dimension */
	int nsd = FieldSupport().NumSD();
	if (nsd != 3) 
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
	                        InputSideSets(surface_spec, 
					geometry,fSurfaces[i]);
/*
	                else if (surface_spec.Name() == "body_boundary")
	                {
	                        // may resize the surfaces array 
	                        InputBodyBoundary(surface_spec, fSurfaces, i);
	                        num_surfaces = fSurfaces.Length();
	                }
*/
	                else
	                        ExceptionT::GeneralFail(caller, 
				"unrecognized surface \"%s\"",
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
	out << " Contact surfaces:\n";
	out << setw(kIntWidth) << "surface"
	    << setw(kIntWidth) << "facets"
	    << setw(kIntWidth) << "size" << '\n';
	for (int j = 0; j < fSurfaces.Length(); j++)
	{		
	  	iArray2DT& surface = fSurfaces[j];

	  	out << setw(kIntWidth) << j+1
	  	    << setw(kIntWidth) << surface.MajorDim()
	  	    << setw(kIntWidth) << surface.MinorDim() << "\n\n";
  	
		/* verbose */
		if (true) {
			surface++;
			surface.WriteNumbered(out);
			surface--;
			out << '\n';
	  	}
	}


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

void PressureBCT:: Compute_Force(DomainIntegrationT& domain, dArray2DT& coord, 
		dArray2DT& force)
{
	force = 0.0;
	domain.TopIP();
	int nsd = domain.NumSD()+1;
	const double* wgs = domain.IPWeights();
	while (domain.NextIP())
	{		
		/* length nnodes */
		const double* T1 = domain.IPDShape(0);
		const double* T2 = domain.IPDShape(1);
		double t1[3] = {0.0,0.0,0.0};
		double t2[3] = {0.0,0.0,-1.0};
		double  n[3] = {0.0,0.0,0.0};
		for (int j = 0; j < nsd; j++)
		{		
			t1[j] = coord.DotColumn(j,T1);
			t2[j] = coord.DotColumn(j,T2);
		}
		CrossProduct(t1,t2,n);

		const double* S = domain.IPShape();
		const double wg = wgs[domain.CurrIP()];
		for (int k = 0; k < force.MajorDim(); k++) {
		for (int j = 0; j < force.MinorDim(); j++) {		
			force(k,j) -= S[k]*n[j]*wg; }}
	}
}

void PressureBCT:: Compute_Stiffness(DomainIntegrationT& domain, dArray2DT& coord, 
		ElementMatrixT& stiffness)
{
	stiffness = 0.0;
	domain.TopIP();
	int nsd = domain.NumSD()+1;
	int nnodes =  domain.ParentDomain().NumNodes();
	const double* wgs = domain.IPWeights();
	while (domain.NextIP())
	{		
		/* length nnodes */
		const double* T1 = domain.IPDShape(0);
		const double* T2 = domain.IPDShape(1);
		double t1[3] = {0.0,0.0,0.0};
		double t2[3] = {0.0,0.0,-1.0};
		for (int j = 0; j < nsd; j++)
		{		
			t1[j] = coord.DotColumn(j,T1);
			t2[j] = coord.DotColumn(j,T2);
		}

		const double* S = domain.IPShape();
		const double wg = wgs[domain.CurrIP()];
		int row = 0, col = 0, k =0;
		for (int I = 0; I < nnodes; I++) { 
		for (int i = 0; i < nsd; i++) {
		col = 0;
		for (int J = 0; J < nnodes; J++) {		
		for (int j = 0; j < nsd; j++) {
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
