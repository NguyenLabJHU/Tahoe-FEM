#include "BridgingScaleT.h"

#include <iostream.h>
#include <iomanip.h>

#include "fstreamT.h"
//#include "FEManagerT.h"
#include "ModelManagerT.h"
#include "StructuralMaterialT.h"
#include "ShapeFunctionT.h"
#include "DomainIntegrationT.h"
#include "eControllerT.h"
#include "Traction_CardT.h"
#include "iAutoArrayT.h"
#include "OutputSetT.h"
#include "ScheduleT.h"

//TEMP: all this for general traction BC implementation?
#include "VariArrayT.h"
#include "nVariArray2DT.h"
#include "VariLocalArrayT.h"

/* services */
#include "EdgeFinderT.h"
#include "GraphT.h"

/* materials lists */
#include "MaterialListT.h"
#include "Material2DT.h"

/* constructor */

using namespace Tahoe;

BridgingScaleT::BridgingScaleT(const ElementSupportT& support, 
	const FieldT& field):
	ElementBaseT(support, field),
	fMaterialList(NULL),
	fBodySchedule(NULL),
	fBody(NumDOF()),
	fTractionBCSet(0),
	fShapes(NULL),
	fLocInitCoords(LocalArrayT::kInitCoords),
	fLocDisp(LocalArrayT::kDisp),
	fDOFvec(NumDOF())
{
	ifstreamT&  in = ElementSupport().Input();
	ostream&    out = ElementSupport().Output();
		
	/* control parameters */
	in >> fGeometryCode; //TEMP - should actually come from the geometry database
	in >> fNumIP;
}

/* destructor */
BridgingScaleT::~BridgingScaleT(void)
{	
	delete fShapes;
	delete fMaterialList;
}

/* accessors */
/* interpolate the nodal field values to the current integration point */
void BridgingScaleT::IP_Interpolate(const LocalArrayT& nodal_u, dArrayT& ip_u) const
{
    /* computed by shape functions */
    ShapeFunction().InterpolateU(nodal_u, ip_u);
}

void BridgingScaleT::IP_Interpolate(const LocalArrayT& nodal_u, dArrayT& ip_u, int ip) const
{
    /* computed by shape functions */
    ShapeFunction().InterpolateU(nodal_u, ip_u, ip);
}

/* allocates space and reads connectivity data */
void BridgingScaleT::Initialize(void)
{
	/* inherited */
	ElementBaseT::Initialize();
	
	/* allocate work space */
	fNEEvec.Allocate(NumElementNodes()*NumDOF());

	/* initialize local arrays */
	SetLocalArrays();

	/* construct shape functions */
	SetShape();

	/* streams */
	ifstreamT& in = ElementSupport().Input();
	ostream&  out = ElementSupport().Output();

	/* output print specifications */
	EchoOutputCodes(in, out);

	/* echo material properties */
	ReadMaterialData(in);	
	WriteMaterialData(out);

	/* get form of tangent */
	GlobalT::SystemTypeT type = TangentType();
	
	/* set form of element stiffness matrix */
	if (type == GlobalT::kSymmetric)
		fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);
	else if (type == GlobalT::kNonSymmetric)
		fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
	else if (type == GlobalT::kDiagonal)
		fLHS.SetFormat(ElementMatrixT::kDiagonal);
}

void BridgingScaleT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	/* inherited */
	ElementBaseT::Equations(eq_1, eq_2);

	/* mark traction BC data as old */
	fTractionBCSet = 0;
}

/* form of tangent matrix */
GlobalT::SystemTypeT BridgingScaleT::TangentType(void) const
{
	/* initialize to lowest precedence */
	GlobalT::SystemTypeT type = GlobalT::kDiagonal;

	for (int i = 0; i < fMaterialList->Length(); i++)
	{
		GlobalT::SystemTypeT e_type = (*fMaterialList)[i]->TangentType();
	
		/* using type precedence */
		type = (e_type > type) ? e_type : type;
	}
	
	return type;
}

/* initialize/finalize step */
void BridgingScaleT::InitStep(void)
{
	/* inherited */
	ElementBaseT::InitStep();

	/* set material variables */
	fMaterialList->InitStep();
}

/* initialize/finalize step */
void BridgingScaleT::CloseStep(void)
{
	/* inherited */
	ElementBaseT::CloseStep();

	/* set material variables */
	fMaterialList->CloseStep();

	/* update element level internal variables */
	if (fMaterialList->HasHistoryMaterials())
	{
		Top();
		while (NextElement())
		{
		        /* compute bridging scale related displacements */
		        ComputeError();
			ComputeFineScaleU();

			ElementCardT& element = CurrentElement();
			if (element.IsAllocated())
			{
				ContinuumMaterialT* pmat = (*fMaterialList)[element.MaterialNumber()];

				/* material update function */
				pmat->UpdateHistory();
			}
		}
	}
}

/* resets to the last converged solution */
void BridgingScaleT::ResetStep(void)
{
	/* inherited */
	ElementBaseT::ResetStep();

	/* update material internal variables */
	if (fMaterialList->HasHistoryMaterials())
	{
		Top();
		while (NextElement())
		{
			ElementCardT& element = CurrentElement();		
			if (element.IsAllocated())
			{
				ContinuumMaterialT* pmat = (*fMaterialList)[element.MaterialNumber()];

				/* material reset function */
				pmat->ResetHistory();
			}
		}
	}
}

/* writing output */
void BridgingScaleT::RegisterOutput(void)
{
//NOTE: could loop over each output mode and register
//      it with the output separately. for now just register
//      "kAtInc"
	
	/* nodal output */
	iArrayT n_counts;
	SetNodalOutputCodes(IOBaseT::kAtInc, fNodalOutputCodes, n_counts);

	/* element output */
	iArrayT e_counts;
	SetElementOutputCodes(IOBaseT::kAtInc, fElementOutputCodes, e_counts);
	ArrayT<StringT> block_ID(fBlockData.Length());
	for (int i = 0; i < block_ID.Length(); i++)
		block_ID[i] = fBlockData[i].ID();

	/* collect variable labels */
	ArrayT<StringT> n_labels(n_counts.Sum());
	ArrayT<StringT> e_labels(e_counts.Sum());
	GenerateOutputLabels(n_counts, n_labels, e_counts, e_labels);

	/* set output specifier */
	StringT set_ID;
	set_ID.Append(ElementSupport().ElementGroupNumber(this) + 1);
	OutputSetT output_set(set_ID, fGeometryCode, block_ID, fConnectivities,
		n_labels, e_labels, false);
		
	/* register and get output ID */
	fOutputID = ElementSupport().RegisterOutput(output_set);
}

//NOTE - this function is/was identical to CSEBaseT::WriteOutput
void BridgingScaleT::WriteOutput(IOBaseT::OutputModeT mode)
{
//TEMP - not handling general output modes yet
	if (mode != IOBaseT::kAtInc)
	{
		cout << "\n BridgingScaleT::WriteOutput: only handling \"at increment\"\n"
		     <<   "     print mode. SKIPPING." << endl;
		return;
	}

	/* map output flags to count of values */
	iArrayT n_counts;
	SetNodalOutputCodes(mode, fNodalOutputCodes, n_counts);
	iArrayT e_counts;
	SetElementOutputCodes(mode, fElementOutputCodes, e_counts);

	/* calculate output values */
	dArray2DT n_values;
	dArray2DT e_values;
	ComputeOutput(n_counts, n_values, e_counts, e_values);

	/* send to output */
	ElementSupport().WriteOutput(fOutputID, n_values, e_values);
}


/***********************************************************************
* Protected
***********************************************************************/

namespace Tahoe {

/* stream extraction operator */
istream& operator>>(istream& in, BridgingScaleT::MassTypeT& type)
{
	int i_type;
	in >> i_type;
	switch (i_type)
	{
		case BridgingScaleT::kNoMass:
			type = BridgingScaleT::kNoMass;
			break;
		case BridgingScaleT::kConsistentMass:
			type = BridgingScaleT::kConsistentMass;
			break;
		case BridgingScaleT::kLumpedMass:
			type = BridgingScaleT::kLumpedMass;
			break;
		default:
			cout << "\n BridgingScaleT::MassTypeT: unknown type: "
			<< i_type<< endl;
			throw eBadInputValue;	
	}
	return in;
}

}

/* initialize local arrays */
void BridgingScaleT::SetLocalArrays(void)
{
	/* dimension */
	fLocInitCoords.Allocate(NumElementNodes(), NumSD());
	fLocDisp.Allocate(NumElementNodes(), NumDOF());

	/* set source */
	ElementSupport().RegisterCoordinates(fLocInitCoords);
	Field().RegisterLocal(fLocDisp);	
}

/* form the residual force vector */
void BridgingScaleT::RHSDriver(void)
{
	// should call derived class function here	
}

/* form global shape function derivatives */
void BridgingScaleT::SetGlobalShape(void)
{
	/* fetch (initial) coordinates */
	SetLocalX(fLocInitCoords);
	
	/* compute shape function derivatives */
	fShapes->SetDerivatives();
}

/* print element group data */
void BridgingScaleT::PrintControlData(ostream& out) const
{
	/* inherited */
	ElementBaseT::PrintControlData(out);

	out << " Element geometry code . . . . . . . . . . . . . = " << fGeometryCode << '\n';
	out << "    eq." << GeometryT::kPoint         << ", point\n";
	out << "    eq." << GeometryT::kLine          << ", line\n";
	out << "    eq." << GeometryT::kQuadrilateral << ", quadrilateral\n";
	out << "    eq." << GeometryT::kTriangle	  << ", triangle\n";
	out << "    eq." << GeometryT::kHexahedron	  << ", hexahedron\n";
	out << "    eq." << GeometryT::kTetrahedron   << ", tetrahedron\n";
	out << " Number of integration points. . . . . . . . . . = " << fNumIP    << '\n';
}

void BridgingScaleT::ReadMaterialData(ifstreamT& in)
{
	/* construct material list */
	int size;
	in >> size;
	fMaterialList = NewMaterialList(size);
	if (!fMaterialList) throw eOutOfMemory;

	/* read */
	fMaterialList->ReadMaterialData(in);
	
	/* check range */
	for (int i = 0; i < fBlockData.Length(); i++)
		if (fBlockData[i].MaterialID() < 0 ||
		    fBlockData[i].MaterialID() >= size)
		{
			cout << "\n BridgingScaleT::ReadMaterialData: material number "
			     << fBlockData[i].MaterialID() + 1 << '\n';
			cout<<    "     for element block " << i + 1 << " is out of range" << endl;
			throw eBadInputValue;
		}
}

/* use in conjunction with ReadMaterialData */
void BridgingScaleT::WriteMaterialData(ostream& out) const
{
	fMaterialList->WriteMaterialData(out);

	/* flush buffer */
	out.flush();
}

/* write all current element information to the stream */
void BridgingScaleT::CurrElementInfo(ostream& out) const
{
	/* inherited */
	ElementBaseT::CurrElementInfo(out);
	dArray2DT temp;
	temp.Allocate(fLocInitCoords.NumberOfNodes(), fLocInitCoords.MinorDim());
	
	out <<   " initial coords:\n";
	temp.Allocate(fLocInitCoords.NumberOfNodes(), fLocInitCoords.MinorDim());
	fLocInitCoords.ReturnTranspose(temp);
	temp.WriteNumbered(out);

	out <<   " displacements:\n";
	temp.Allocate(fLocDisp.NumberOfNodes(), fLocDisp.MinorDim());
	fLocDisp.ReturnTranspose(temp);
	temp.WriteNumbered(out);
}

/***********************************************************************
* Private
***********************************************************************/

/* return the default number of element nodes */
int BridgingScaleT::DefaultNumElemNodes(void) const
{
	switch (fGeometryCode)
	{
		case GeometryT::kLine:
			return 2;
		case GeometryT::kQuadrilateral:
			return 4;
		case GeometryT::kTriangle:
			return 3;
		case GeometryT::kHexahedron:
			return 8;
		case GeometryT::kTetrahedron:
			return 4;
		case GeometryT::kPentahedron:
			return 6;
		default:
			cout << "\n BridgingScaleT::DefaultNumElemNodes: unknown geometry code: "
			     << fGeometryCode << endl;
			return 0;
	}
}
//NOTE: needed because ExodusII does not store ANY information about
//      empty element groups, which causes trouble for parallel execution
//      when a partition contains no element from a group.

void BridgingScaleT::ComputeError(void)
{
  /* compute the error caused by projecting the "exact" (MD) solution onto a
   * finite dimensional basis set (FEM basis) */



}

void BridgingScaleT::ComputeFineScaleU(void)
{
  /* compute the fine scale displacement, ie the "exact" (MD) solution minus
   * error interpolated over a given finite element */



}
