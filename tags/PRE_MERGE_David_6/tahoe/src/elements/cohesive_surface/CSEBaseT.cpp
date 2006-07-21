/* $Id: CSEBaseT.cpp,v 1.39 2005-07-12 05:50:11 paklein Exp $ */
/* created: paklein (11/19/1997) */
#include "CSEBaseT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include "ifstreamT.h"
#include "SurfaceShapeT.h"
#include "iAutoArrayT.h"
#include "OutputSetT.h"
#include "ElementSupportT.h"
#include "ModelManagerT.h"
#include "ParameterContainerT.h"
#include "CommunicatorT.h"

using namespace Tahoe;

/* initialize static data */
const int CSEBaseT::NumNodalOutputCodes = 5;
static const char* NodalOutputNames[5] = {
	"coordinates",
	"displacements",
	"displacement_jump",
	"traction",
	"material_output"};

const int CSEBaseT::NumElementOutputCodes = 3;
static const char* ElementOutputNames[3] = {
	"centroid",
	"cohesive_energy",
	"avg_traction"};

CSEBaseT::AndOrT CSEBaseT::int2AndOrT(int i) {
	if (i == kAND) return kAND;
	else if (i == kOR) return kOR;
	else ExceptionT::GeneralFail("CSEBaseT::int2AndOrT", "%d out of range", i);
	return kAND; /* dummy */
}
CSEBaseT::OpT CSEBaseT::int2OpT(int i) {
	if (i == kEqual) return kEqual;
	else if (i == kLess) return kLess;
	else if (i == kGreater) return kGreater;
	else ExceptionT::GeneralFail("CSEBaseT::int2OpT", "%d out of range", i);
	return kEqual; /* dummy */
}
CSEBaseT::CoordinateT CSEBaseT::int2CoordinateT(int i) {
	if (i == kX) return kX;
	else if (i == kY) return kY;
	else if (i == kZ) return kZ;
	else ExceptionT::GeneralFail("CSEBaseT::int2CoordinateT", "%d out of range", i);
	return kX; /* dummy */
}

#ifndef _FRACTURE_INTERFACE_LIBRARY_
/* constructor */
CSEBaseT::CSEBaseT(const ElementSupportT& support):
	ElementBaseT(support),
	fGeometryCode(GeometryT::kNone),
	fNumIntPts(-1),
	fAxisymmetric(false),
	fCloseSurfaces(false),
	fOutputArea(false),
	fLocInitCoords1(LocalArrayT::kInitCoords),
	fLocCurrCoords(LocalArrayT::kCurrCoords),
	fFractureArea(0.0),
	fShapes(NULL),
	fOutputGlobalTractions(false),
	fpc_AndOr(kAND)
{
	SetName("CSE_base");	
}
#else
/* constructor */
CSEBaseT::CSEBaseT(ElementSupportT& support):
	ElementBaseT(support),
	fGeometryCode(GeometryT::kNone),
	fNumIntPts(-1),
	fAxisymmetric(false),
	fCloseSurfaces(false),
	fOutputArea(false),
	fLocInitCoords1(LocalArrayT::kInitCoords),
	fLocCurrCoords(LocalArrayT::kCurrCoords),
	fFractureArea(0.0),
	fShapes(NULL),
	fOutputGlobalTractions(false),
	fpc_AndOr(kAND)
{
	SetName("CSE_base");

	int i_code = ElementSupport().ReturnInputInt(ElementSupportT::kGeometryCode);
	fGeometryCode = GeometryT::int2CodeT(i_code);
	fNumIntPts =  ElementSupport().ReturnInputInt(ElementSupportT::kNumIntPts);
	fCloseSurfaces =  ElementSupport().ReturnInputInt(ElementSupportT::kCloseSurface);
	fOutputArea =  ElementSupport().ReturnInputInt(ElementSupportT::kOutputArea);

	/* checks */
	if (NumSD() == 2 && fGeometryCode != GeometryT::kLine)
	{
		throw ExceptionT::kBadInputValue;
	}
	else if (NumSD() == 3 &&
	         fGeometryCode != GeometryT::kQuadrilateral &&
	         fGeometryCode != GeometryT::kTriangle)
	{
		throw ExceptionT::kBadInputValue;
	}
	
	if (fCloseSurfaces != 0 &&
	    fCloseSurfaces != 1) throw ExceptionT::kBadInputValue;
	if (fOutputArea != 0 &&
	    fOutputArea != 1) throw ExceptionT::kBadInputValue;
}
#endif // _FRACTURE_INTERFACE_LIBRARY_

/* destructor */
CSEBaseT::~CSEBaseT(void)
{
	delete fShapes;
	fShapes = NULL;
}

/* initial condition/restart functions (per time sequence) */
void CSEBaseT::InitialCondition(void)
{
	const char caller[] = "CSEBaseT::InitialCondition";

	/* inherited */
	ElementBaseT::InitialCondition();

	/* initialize element status flags */
	int nel = NumElements();
	for (int i = 0; i < nel; i++)
		fElementCards[i].Flag() = ElementCardT::kON;

	/* create pre-crack */
	if (fpc_coordinate.Length() > 0)
	{
		/* node map of facet 1 */
		iArrayT facet1;
		(fShapes->NodesOnFacets()).RowAlias(0, facet1);

		/* loop over elements */
		dArrayT ip_X(NumSD());
		Top();
		while (NextElement())
		{
			/* current element */
			ElementCardT& element = CurrentElement();
		
			/* get ref geometry (1st facet only) */
			fNodes1.Collect(facet1, element.NodesX());
			fLocInitCoords1.SetLocal(fNodes1);
		
			/* loop over integration points */
			bool all_ip_off = true;
			fShapes->TopIP();
			while (all_ip_off && fShapes->NextIP())
			{
				/* integration point coordinates */
				fShapes->Interpolate(fLocInitCoords1, ip_X);

				/* test the integration point */
				bool ip_off = (fpc_AndOr == kAND) ? true : false;
				for (int i = 0; i < fpc_coordinate.Length(); i++)
				{
					/* evaluate test */
					bool off = true;
					double value = fpc_value[i];
					double coord = ip_X[fpc_coordinate[i]];
					if (fpc_op[i] == kEqual)
						off = (fabs(coord - value) < kSmall);
					else if (fpc_op[i] == kLess)
						off = (coord < value);
					else if (fpc_op[i] == kGreater)
						off = (coord > value);
					else
						ExceptionT::GeneralFail(caller, "unrecognized operator %d", fpc_op[i]);
				
					/* accumulate test */
					if (fpc_AndOr == kAND)
						ip_off = (off && ip_off);
					else /* or */
						ip_off = (off || ip_off);					
				}
				
				/* accumulate integration point */
				all_ip_off = (all_ip_off && ip_off);
			}
			
			/* turn element off */
			if (all_ip_off) element.Flag() = ElementCardT::kOFF;
		}
	}
}

#ifdef _FRACTURE_INTERFACE_LIBRARY_	
	/* Initialize fields passed in from the outside */
void CSEBaseT::InitStep(void) {};
#endif

/* finalize time increment */
void CSEBaseT::CloseStep(void)
{
	/* inherited */
	ElementBaseT::CloseStep();

	/* deactivate marked elements */
	int nel = NumElements();
	for (int i = 0; i < nel; i++)
	{
		ElementCardT::StatusT& flag = fElementCards[i].Flag();
		flag = (flag == ElementCardT::kMarked) ? ElementCardT::kOFF : flag;
	}
}

/* resets to the last converged solution */
GlobalT::RelaxCodeT CSEBaseT::ResetStep(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = ElementBaseT::ResetStep();

	/* unset marks */
	int nel = NumElements();
	for (int i = 0; i < nel; i++)
	{
		ElementCardT::StatusT& flag = fElementCards[i].Flag();
		flag = (flag == ElementCardT::kMarked) ? ElementCardT::kON : flag;
	}
	
	return relax;
}

#ifndef _FRACTURE_INTERFACE_LIBRARY_
/* solution calls */
void CSEBaseT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
//TEMP
#pragma unused(field)
#pragma unused(node)
#pragma unused(force)
	//not implemented
}
#endif

/* returns the energy as defined by the derived class types */
double CSEBaseT::InternalEnergy(void) { return 0.0; } //not implemented

/* writing output */
void CSEBaseT::RegisterOutput(void)
{
//NOTE: could loop over each output mode and register
//      it with the output separately. for now just register
//      "kAtInc"

	/* "deformed" geometry */
	GeometryT::CodeT geo_code;
	switch (fGeometryCode)
	{
		case GeometryT::kLine:		
			geo_code = GeometryT::kQuadrilateral;
			break;

		case GeometryT::kQuadrilateral:
			geo_code = GeometryT::kHexahedron;
			break;

		case GeometryT::kTriangle:
			geo_code = GeometryT::kPentahedron;
			break;

		default:	
			ExceptionT::GeneralFail("CSEBaseT::RegisterOutput",
				"could not translate geometry code %d to a pseudo-geometry for the volume",
				fGeometryCode);
	}	

	/* nodal output */
	iArrayT n_counts;
	SetNodalOutputCodes(IOBaseT::kAtInc, fNodalOutputCodes, n_counts);

	/* element output */
	iArrayT e_counts;
	SetElementOutputCodes(IOBaseT::kAtInc, fElementOutputCodes, e_counts);
	
	/* collect variable labels */
	ArrayT<StringT> n_labels(n_counts.Sum());
	ArrayT<StringT> e_labels(e_counts.Sum());
	GenerateOutputLabels(n_counts, n_labels, e_counts, e_labels);

#ifndef _FRACTURE_INTERFACE_LIBRARY_

	/* collect output connectivities */
	ModelManagerT& model = ElementSupport().ModelManager();
	ArrayT<const iArray2DT*> output_connects;		
	model.ElementGroupPointers(fOutputBlockID, output_connects);

	/* set output specifier */
	OutputSetT output_set(geo_code, fOutputBlockID, output_connects, n_labels, e_labels, false);

	/* register and get output ID */
	fOutputID = ElementSupport().RegisterOutput(output_set);
#else
	ElementSupport().RegisterOutput(n_labels,e_labels);
#endif
}

//NOTE - this function is identical to ContinuumElementT::WriteOutput
void CSEBaseT::WriteOutput(void)
{
	/* fracture area */
	if (fOutputArea)
	{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
		farea_out << setw(kDoubleWidth) << ElementSupport().Time();
		farea_out << setw(kDoubleWidth) << fFractureArea << endl;
#endif
	}

	/* regular output */
	IOBaseT::OutputModeT mode = IOBaseT::kAtInc;

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

/* compute specified output parameter and send for smoothing */
void CSEBaseT::SendOutput(int kincode)
{
	/* output flags */
	iArrayT flags(fNodalOutputCodes.Length());

	/* set flags to get desired output */
	flags = IOBaseT::kAtNever;
	switch (kincode)
	{
		case NodalDisp:
		    flags[NodalDisp] = NumDOF();
			break;
		case NodalDispJump:
		    flags[NodalDispJump] = 1;
			break;
		case NodalTraction:
		    flags[NodalTraction] = 1;
			break;
		case MaterialData:
		    flags[MaterialData] = 1;
			break;
		default:
			ExceptionT::BadInputValue("CSEBaseT::SendKinematic", "invalid output code: %d", kincode);
	}

	/* number of output values */
	iArrayT n_counts;
	SetNodalOutputCodes(IOBaseT::kAtInc, flags, n_counts);

	/* reset averaging workspace */
	ElementSupport().ResetAverage(n_counts.Sum());

	/* set flags for no element output */
	iArrayT e_counts(fElementOutputCodes.Length());
	e_counts = 0;

	/* generate nodal values */
	dArray2DT e_values, n_values;
	ComputeOutput(n_counts, n_values, e_counts, e_values);
}

/* resolve the output variable label into the output code and offset within the output. */
void CSEBaseT::ResolveOutputVariable(const StringT& variable, int& code, int& offset)
{
	/* search output labels */
	code = -1;
	offset = -1;
	iArrayT e_counts(NumElementOutputCodes);
	e_counts = 0;
	iArrayT n_codes(NumNodalOutputCodes);
	for (int i = 0; code == -1 && i < NumNodalOutputCodes; i++)
	{
		ArrayT<StringT> n_labels, e_labels;
		n_codes = 0;
		n_codes[i] = 1;
		
		iArrayT n_counts;
		SetNodalOutputCodes(IOBaseT::kAtInc, n_codes, n_counts);
		GenerateOutputLabels(n_counts, n_labels, e_counts, e_labels);
		
		for (int j = 0; offset == -1 && j < n_labels.Length(); j++)
			if (n_labels[j] == variable) /* found */ {
				code = i;
				offset = j;
			}
	}
	
	/* inherited */
	if (code == -1 || offset == -1)
		ElementBaseT::ResolveOutputVariable(variable, code, offset);
}

/* describe the parameters needed by the interface */
void CSEBaseT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ElementBaseT::DefineParameters(list);

	ParameterT axisymmetric(fAxisymmetric, "axisymmetric");
	axisymmetric.SetDefault(fAxisymmetric);
	list.AddParameter(axisymmetric);

	ParameterT close_surfaces(fCloseSurfaces, "close_surfaces");
	close_surfaces.SetDefault(fCloseSurfaces);
	list.AddParameter(close_surfaces);
	
	ParameterT output_area(fOutputArea, "output_area");
	output_area.SetDefault(fOutputArea);
	list.AddParameter(output_area);
}

/* information about subordinate parameter lists */
void CSEBaseT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ElementBaseT::DefineSubs(sub_list);

	/* geometry and integration rule (inline) */
	sub_list.AddSub("surface_geometry", ParameterListT::Once, true);

	/* define pre-crack region */
	sub_list.AddSub("pre_crack", ParameterListT::ZeroOrOnce);

	/* nodal output codes (optional) */
	sub_list.AddSub("surface_element_nodal_output", ParameterListT::ZeroOrOnce);
	sub_list.AddSub("surface_element_element_output", ParameterListT::ZeroOrOnce);
}

/* return the description of the given inline subordinate parameter list */
void CSEBaseT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order,
	SubListT& sub_lists) const
{
	if (name == "surface_geometry")
	{
		/* choice */
		order = ParameterListT::Choice;

		/* element geometries */
		sub_lists.AddSub(GeometryT::ToString(GeometryT::kLine));		
		sub_lists.AddSub(GeometryT::ToString(GeometryT::kQuadrilateral));
		sub_lists.AddSub(GeometryT::ToString(GeometryT::kTriangle));
	}
	else /* inherited */
		ElementBaseT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT */
ParameterInterfaceT* CSEBaseT::NewSub(const StringT& name) const
{
	/* look for geometry */
	ParameterInterfaceT* geom = GeometryT::New(name);
	if (geom)
		return geom;

	if (name == "surface_element_nodal_output")
	{
		ParameterContainerT* node_output = new ParameterContainerT(name);
		
		/* all false by default */
		for (int i = 0; i < NumNodalOutputCodes; i++) {
			ParameterT output(ParameterT::Integer, NodalOutputNames[i]);
			output.SetDefault(1);
			node_output->AddParameter(output, ParameterListT::ZeroOrOnce);
		}

		return node_output;
	}
	else if (name == "surface_element_element_output")
	{
		ParameterContainerT* element_output = new ParameterContainerT(name);
		
		/* all false by default */
		for (int i = 0; i < NumElementOutputCodes; i++) {
			ParameterT output(ParameterT::Integer, ElementOutputNames[i]);
			output.SetDefault(1);
			element_output->AddParameter(output, ParameterListT::ZeroOrOnce);
		}

		return element_output;	
	}
	else if (name == "pre_crack")
	{
		ParameterContainerT* pre_crack = new ParameterContainerT(name);
		pre_crack->SetDescription("rules defining the pre-crack");
		pre_crack->SetSubSource(this);
		
		/* how to apply the rules */
		ParameterT and_or(ParameterT::Enumeration, "and_or");
		and_or.AddEnumeration("AND", kAND);
		and_or.AddEnumeration( "OR", kOR);
		and_or.SetDefault(fpc_AndOr);
		pre_crack->AddParameter(and_or);
	
		/* must have at least one rule */
		pre_crack->AddSub("pre_crack_rule", ParameterListT::OnePlus);
	
		return pre_crack;
	}
	else if (name == "pre_crack_rule")
	{
		ParameterContainerT* pre_crack_rule = new ParameterContainerT(name);
	
		ParameterT direction(ParameterT::Enumeration, "coordinate");
		direction.AddEnumeration("x", 0);
		direction.AddEnumeration("y", 1);
		direction.AddEnumeration("z", 2);
		direction.SetDefault(0);
		pre_crack_rule->AddParameter(direction);

		ParameterT op(ParameterT::Enumeration, "op");
		op.AddEnumeration(  "equal", 0);
		op.AddEnumeration(   "less", 1);
		op.AddEnumeration("greater", 2);
		op.SetDefault(0);
		pre_crack_rule->AddParameter(op);	
	
		pre_crack_rule->AddParameter(ParameterT::Double, "value");

		return pre_crack_rule;
	}
	else /* inherited */
		return ElementBaseT::NewSub(name);
}

/* accept parameter list */
void CSEBaseT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "CSEBaseT::TakeParameterList";

	/* axisymmetry */
	fAxisymmetric = list.GetParameter("axisymmetric");
	if (fAxisymmetric && NumSD() != 2) /* check */
		ExceptionT::GeneralFail(caller, "expecting 2 dimensions not %d with axisymmetry", NumSD());

	/* element geometry - need to set geometry before calling inherited method */
	const ParameterListT& geom = list.GetListChoice(*this, "surface_geometry");
	fGeometryCode = GeometryT::string2CodeT(geom.Name());
	fNumIntPts = geom.GetParameter("num_ip");

	/* inherited */
	ElementBaseT::TakeParameterList(list);

	/* take parameters */
	fCloseSurfaces = list.GetParameter("close_surfaces");
	fOutputArea = list.GetParameter("output_area");

	/* pre-crack */
	const ParameterListT* pre_crack = list.List("pre_crack");
	if (pre_crack) {
		fpc_AndOr = int2AndOrT(pre_crack->GetParameter("and_or"));
		int num_rules = pre_crack->NumLists("pre_crack_rule");
		fpc_coordinate.Dimension(num_rules);
		fpc_op.Dimension(num_rules);
		fpc_value.Dimension(num_rules);
		for (int i = 0; i < num_rules; i++) {
			const ParameterListT& pre_crack_rule = pre_crack->GetList("pre_crack_rule", i);
			fpc_coordinate[i] = int2CoordinateT(pre_crack_rule.GetParameter("coordinate"));
			fpc_op[i] = int2OpT(pre_crack_rule.GetParameter("op"));
			fpc_value[i] = pre_crack_rule.GetParameter("value");
		}
	}

	/* nodal output codes */
	fNodalOutputCodes.Dimension(NumNodalOutputCodes);
	fNodalOutputCodes = 0;
	const ParameterListT* nodal_output = list.List("surface_element_nodal_output");
	if (nodal_output)
		for (int i = 0; i < NumNodalOutputCodes; i++)
		{
			/* look for entry */
			const ParameterT* nodal_value = nodal_output->Parameter(NodalOutputNames[i]);
			if (nodal_value) {
				int do_write = *nodal_value;
				if (do_write == 1)
					fNodalOutputCodes[i] = 1;
				else if (i == NodalTraction && do_write == 2) {
					fNodalOutputCodes[i] = 1;
					fOutputGlobalTractions = true;
				}
			}
		}

	/* element output codes */
	fElementOutputCodes.Dimension(NumElementOutputCodes);
	fElementOutputCodes = 0;
	const ParameterListT* element_output = list.List("surface_element_element_output");
	if (element_output)
		for (int i = 0; i < NumElementOutputCodes; i++)
		{
			/* look for entry */
			const ParameterT* element_value = element_output->Parameter(ElementOutputNames[i]);
			if (element_value) {
				int do_write = *element_value;
				if (do_write == 1)
					fElementOutputCodes[i] = 1;
			}
		}

	/* dimensions */
	int num_facet_nodes = NumFacetNodes();

	/* initialize local arrays */
	fLocInitCoords1.Dimension(num_facet_nodes, NumSD());
	fLocCurrCoords.Dimension(NumElementNodes(), NumSD());
	ElementSupport().RegisterCoordinates(fLocInitCoords1);
	ElementSupport().RegisterCoordinates(fLocCurrCoords);

	/* construct surface shape functions */
	fShapes = new SurfaceShapeT(fGeometryCode, fNumIntPts, NumElementNodes(), 
		num_facet_nodes, NumDOF(), fLocInitCoords1);
	if (!fShapes) throw ExceptionT::kOutOfMemory;
	fShapes->Initialize();

	/* work space */
	fNodes1.Dimension(num_facet_nodes);
	int nee = NumElementNodes()*NumDOF();
	fNEEvec.Dimension(nee);
	fNEEmat.Dimension(nee);

	/* close surfaces */
	if (fCloseSurfaces) CloseSurfaces();

#ifndef _FRACTURE_INTERFACE_LIBRARY_
	/* output stream */
	if (fOutputArea == 1)
	{
		/* generate file name */
		StringT name = ElementSupport().InputFile();
		name.Root();
		name.Append(".grp", ElementSupport().ElementGroupNumber(this) + 1);
		name.Append(".fracture");
		
		/* open stream */
		farea_out.open(name);
	}
#endif
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* define the elements blocks for the element group */
void CSEBaseT::CollectBlockInfo(const ParameterListT& list, ArrayT<StringT>& block_ID,  
	ArrayT<int>& mat_index) const
{
	const char caller[] = "CSEBaseT::CollectBlockInfo";

	/* inherited */
	ElementBaseT::CollectBlockInfo(list, block_ID, mat_index);

	/* quick exit */
	if (block_ID.Length() == 0) {
		ElementSupport().Communicator().Max(0); /* complete Max below */
		return;
	}

	/* write output over the original connectivities */
	CSEBaseT* non_const_this = (CSEBaseT*) this;
	non_const_this->fOutputBlockID = block_ID;

	/* geometry information */
	ModelManagerT& model = ElementSupport().ModelManager();
	int nel, nen = 0;
	for (int i = 0; nen == 0 && i < block_ID.Length(); i++)
		model.ElementGroupDimensions (block_ID[i], nel, nen);

	/* group may be empty for multi-processor calculations */
	nen = ElementSupport().Communicator().Max(nen);

	/* check for higher order elements */
	int nsd = NumSD();
	if ((nsd == 2 && nen != 4 && nen != 6) || 
	    (nsd == 3 && nen != 8 && nen != 16))
	{
#ifndef _FRACTURE_INTERFACE_LIBRARY_	
		/* message */
		ostream& out = ElementSupport().Output();
		cout << "\n " << caller << ": detected higher order elements\n";
		out  << "\n " << caller << ": detected higher order elements\n";
#endif

		/* nen: 8 -> 6 */
		int map_2D[] = {0, 1, 2, 3, 4, 6}; 

		/* nen: 20 -> 16 */
		int map_3D[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}; 

		/* node map */
		iArrayT map((nsd == 2) ? 6 : 16, (nsd == 2) ? map_2D : map_3D); 

		/* loop over connectivity blocks */
		for (int b = 0; b < block_ID.Length(); b++)
		{
			/* send new connectivities to model manager */
			StringT& new_id = block_ID[b];
			StringT old_id = new_id;
			new_id.Append(b+1, 3);
			
			/* see if new_id is already present */
			int new_dex = model.ElementGroupIndex(new_id);
			if (new_dex == ModelManagerT::kNotFound)
			{
#ifndef _FRACTURE_INTERFACE_LIBRARY_	
				/* message */
		     	cout << "     translating element block ID " << old_id << endl;	     	
		     	out  << "     translating element block ID " << old_id << endl;
#endif
				/* translate */
				const iArray2DT& source = model.ElementGroup(fOutputBlockID[b]);
				iArray2DT dest(source.MajorDim(), map.Length());
				for (int i = 0; i < dest.MajorDim(); i++) {
					int* a = dest(i);
					const int* b = source(i);
					for (int j = 0; j < map.Length(); j++)
						*a++ = b[map[j]];	
				}

				/* send new connectivities to model manager */
				if (!model.RegisterElementGroup (new_id, dest, GeometryT::kNone, true))
					ExceptionT::GeneralFail(caller, "could not register element block ID \"%s\"",
						new_id.Pointer());
			}

#ifndef _FRACTURE_INTERFACE_LIBRARY_	
			/* message */
			cout << "     block ID " << old_id << " replaced by ID " << new_id << endl;		
			out  << "     block ID " << old_id << " replaced by ID " << new_id << endl;
#endif
		}
	}
}

void CSEBaseT::SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const
{
	/* initialize */
	counts.Dimension(flags.Length());
	counts = 0;
	
	if (flags[NodalCoord] == mode)
		counts[NodalCoord] = NumSD();
	if (flags[NodalDisp] == mode)
		counts[NodalDisp] = NumDOF();
	if (flags[NodalDispJump] == mode)
		counts[NodalDispJump] = 1;
	if (flags[NodalTraction] == mode)
		counts[NodalTraction] = 1;
}

void CSEBaseT::SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
	/* initialize */
	counts.Dimension(flags.Length());
	counts = 0;

	if (flags[Centroid] == mode)
		counts[Centroid] = NumSD();
	if (flags[CohesiveEnergy] == mode)
		counts[CohesiveEnergy] = 1;
	if (flags[Traction] == mode)
		counts[Traction] = 1;
}

/* construct output labels array */
void CSEBaseT::GenerateOutputLabels(const iArrayT& n_codes, ArrayT<StringT>& n_labels,
	const iArrayT& e_codes, ArrayT<StringT>& e_labels) const
{
	/* allocate nodal output labels */
	n_labels.Dimension(n_codes.Sum());

	int count = 0;
	if (n_codes[NodalDisp])
	{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
		/* labels from the field */
		const ArrayT<StringT>& labels = Field().Labels();
		for (int i = 0; i < labels.Length(); i++)
			n_labels[count++] = labels[i];
#else
		const char* labels[] = {"D_1", "D_2", "D_3"};
		for (int i = 0; i < NumSD(); i++)
			n_labels[count++] = labels[i];
#endif
	}

	if (n_codes[NodalCoord])
	{
		const char* xlabels[] = {"x1", "x2", "x3"};
		for (int i = 0; i < NumSD(); i++)
			n_labels[count++] = xlabels[i];
	}

	if (n_codes[NodalDispJump]) n_labels[count++] = "jump";
	if (n_codes[NodalTraction]) n_labels[count++] = "Tmag";
	
	/* allocate nodal output labels */
	e_labels.Dimension(e_codes.Sum());
	count = 0;
	if (e_codes[Centroid])
	{
		const char* xlabels[] = {"xc_1", "xc_2", "xc_3"};
		for (int i = 0; i < NumSD(); i++)
			e_labels[count++] = xlabels[i];
	}
	if (e_codes[CohesiveEnergy]) e_labels[count++] = "phi";
	if (e_codes[Traction]) e_labels[count++] = "Tmag";
}

/* write all current element information to the stream */
void CSEBaseT::CurrElementInfo(ostream& out) const
{
#ifndef _FRACTURE_INTERFACE_LIBRARY_
	/* inherited */
	ElementBaseT::CurrElementInfo(out);
	dArray2DT temp;

	out <<   " initial coords:\n";
	temp.Dimension(fLocInitCoords1.NumberOfNodes(), fLocInitCoords1.MinorDim());
	fLocInitCoords1.ReturnTranspose(temp);
	temp.WriteNumbered(out);

	out <<   " current coords:\n";
	temp.Dimension(fLocCurrCoords.NumberOfNodes(), fLocCurrCoords.MinorDim());
	fLocCurrCoords.ReturnTranspose(temp);
	temp.WriteNumbered(out);
#else
#pragma unused(out)
#endif
}

/***********************************************************************
* Private
***********************************************************************/

/* close surfaces to zero gap */
void CSEBaseT::CloseSurfaces(void) const
{
	/* get coordinates */
	const dArray2DT& init_coords = ElementSupport().InitialCoordinates();

	/* local nodes numbers on each facet */
	const iArray2DT& facetnodes = fShapes->NodesOnFacets();
		
	/* collapse elements */
	int nel = NumElements();
	for (int i = 0; i < nel; i++)
	{			
		const int* pfacet1 = facetnodes(0);
		const int* pfacet2 = facetnodes(1);
		const iArrayT& elemnodes = fElementCards[i].NodesX();
		const int* nodes = elemnodes.Pointer();

		for (int j = 0; j < facetnodes.MinorDim(); j++)
		{
			/* facet coordinates */		
			double* px1 = const_cast<double*>(init_coords(nodes[*pfacet1++]));
			double* px2 = const_cast<double*>(init_coords(nodes[*pfacet2++]));
				
			for (int k = 0; k < NumSD(); k++)
			{
				double x_mid = 0.5*(*px1 + *px2);
				*px1++ = x_mid;
				*px2++ = x_mid;
			}
		}		
	}
}

/* return the default number of element nodes */
int CSEBaseT::DefaultNumElemNodes(void) const
{
	/* return number of element nodes given the facet geometry */
	switch (fGeometryCode)
	{
		case GeometryT::kLine:
			return 4;
		case GeometryT::kQuadrilateral:
			return 8;
		case GeometryT::kTriangle:
			return 6;
		default:
#ifndef _FRACTURE_INTERFACE_LIBRARY_
			cout << "\n CSEBaseT::DefaultNumElemNodes: unknown geometry code: "
			     << fGeometryCode << endl;
#endif
			return 0;
	}
}
//NOTE: needed because ExodusII does not store ANY information about
//      empty element groups, which causes trouble for parallel execution
//      when a partition contains no element from a group.
