/* $Id: ACME_Contact3DT.cpp,v 1.6 2004-07-15 08:26:08 paklein Exp $ */
/* created: paklein (10/15/2000) */

#include "ACME_Contact3DT.h"
#include "ElementSupportT.h"

/* library support options */
#ifdef __ACME__


/* constructor */

using namespace Tahoe;

ACME_Contact3DT::ACME_Contact3DT(const ElementSupportT& support, const FieldT& field):
	PenaltyContact3DT(support, field),
	fContactSearch(NULL)
{

}

/* destructor */
ACME_Contact3DT::~ACME_Contact3DT(void)
{
	delete fContactSearch;
}

/* initialization after constructor */
void ACME_Contact3DT::Initialize(void)
{
	/* inherited */
	PenaltyContact3DT::Initialize();

	/* write ACME version to output */
	ostream& out = ElementSupport().Output();
	out << " ACME version number . . . . . . . . . . . . . . = "
	    << fContactSearch->Version() << '\n';
}

/* restart functions */
void ACME_Contact3DT::ReadRestart(istream& in)
{
#pragma unused(in)
	cout << "\n ACME_Contact3DT::ReadRestart: not implemented" << endl;
	throw ExceptionT::kGeneralFail;
}

void ACME_Contact3DT::WriteRestart(ostream& out) const
{
#pragma unused(out)
	cout << "\n AugLagContact2DT::WriteRestart: not implemented" << endl;
	throw ExceptionT::kGeneralFail;
}

/***********************************************************************
* Protected
***********************************************************************/

/* element data - accepts surfaces as side sets only */
void ACME_Contact3DT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
	int num_surfaces;
	in >> num_surfaces;
	out << " Number of contact surfaces. . . . . . . . . . . = "
	    << num_surfaces << '\n';
	if (num_surfaces < 1) throw ExceptionT::kBadInputValue;

	/* read contact surfaces */
	fSurfaces.Dimension(num_surfaces);
	for (int i = 0; i < fSurfaces.Length(); i++)
		InputSideSets(in, out, fSurfaces[i]);

	/* convert quad faces */
//NOTE: remove to keep tri's	
	out << "\n ACME_Contact3DT::EchoConnectivityData: subdividing\n"
	    <<   "     4-nodes facets on surfaces into (2) triangular facets" << endl;
	for (int ii = 0; ii < fSurfaces.Length(); ii++)
		ConvertQuadToTri(fSurfaces[ii]);
//NOTE: remove to keep tri's

	/* echo data and correct numbering offset */
	out << " Contact surfaces:\n";
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
	  	if (fFEManager.PrintInput())
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
		out << " Found empty contact surfaces:\n\n";
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
	
	/* collect striker nodes from surfaces */
	StrikersFromSurfaces();

	/* echo */
	if (ElementSupport().PrintInput())
	{
		out << "\n Striker nodes:\n";
		fStrikerTags++;
		out << fStrikerTags.wrap(8) << '\n';
		fStrikerTags--;	
	}
	
	/* allocate striker coords */
	fStrikerCoords.Dimension(fStrikerTags.Length(), fNumSD);
}

/* steps in setting contact configuration */
bool ACME_Contact3DT::SetActiveInteractions(void)
{
	int last_num_active = fActiveStrikers.Length();

	/* collect contact node coordinates */
	fStrikerCoords.RowCollect(fStrikerTags, fNodes->CurrentCoordinates());

	/* set ACME - only 1 node block */
	ContactSearch::ContactNode_Configuration config =
		ContactSearch::CURRENT_CONFIG;
	ContactSearch::ContactErrorCode error =
		fContactSearch->Set_Node_Block_Configuration(config, 1,
			fStrikerCoords.Pointer());
	if (error != ContactSearch::NO_ERROR)
	{
		cout << "\n ACME_Contact3DT::SetActiveInteractions: Set_Node_Block_Configuration\n"
		     <<   "     returned error code: " << error << endl;
		throw ExceptionT::kGeneralFail;	
	}

	/* search for interactions */
	error = fContactSearch->Static_Search_1_Configuration();
	if (error != ContactSearch::NO_ERROR)
	{
		cout << "\n ACME_Contact3DT::SetActiveInteractions: Static_Search_1_Configuration\n"
		     <<   "     returned error code: " << error << endl;
		throw ExceptionT::kGeneralFail;	
	}

	/* number of interactions */
	int num_interactions, data_size;
	fContactSearch->Size_NodeFace_Interactions(num_interactions, data_size);
	fnode_block_ids.Resize(num_interactions, false);
	//fnode_indices_in_block.Resize(num_interactions, false);
	//fface_block_ids.Resize(num_interactions, false);
	//fface_indices_in_block.Resize(num_interactions, false);
	fface_proc.Resize(num_interactions, false);
	fdata.Resize(data_size*num_interactions, false); //not used for anything yet

	/* write directly into class work space */
	fActiveStrikers.Resize(num_interactions, false);
	fHitSurface.Resize(num_interactions, false);
	fHitFacets.Resize(num_interactions, false);
	
	/* collect interactions */
	fContactSearch->Get_NodeFace_Interactions(
		fnode_block_ids.Pointer(),
		fActiveStrikers.Pointer(),
		fHitSurface.Pointer(),
		fHitFacets.Pointer(),
		fface_proc.Pointer(),
		fdata.Pointer());
		
	/* translate configuration */
	for (int i = 0; i < num_interactions; i++)
	{
		/* just 1 block for now */
		if (fnode_block_ids[i] != 1)
		{
			cout << "\n ACME_Contact3DT::SetActiveInteractions: unexpected node\n"
			     <<   "     block ID: " << fnode_block_ids[i] << endl;
			throw ExceptionT::kGeneralFail;
		}
		
		/* striker */
		fActiveStrikers[i] = fStrikerTags[fActiveStrikers[i] - 1]; //OFFSET
		
		/* surface */
		fHitSurface[i]--; //OFFSET
		
		/* facet on surface */
		fHitFacets[i]--; //OFFSET
	}

	/* simple check of changing configuration */
	if (last_num_active == 0 && fActiveStrikers.Length() == 0)
		return false;
	else
		return true;
}

/***********************************************************************
* Protected
***********************************************************************/

/* constuct ContactSearch object */
void ACME_Contact3DT::SetWorkSpace(void)
{
	/* inherited */
	PenaltyContact3DT::SetWorkSpace();

	/* ACME parameters */
	int Dimensionality = 3; //TEMP - must be 3D
	int Number_of_States = 1; //TEMP - must be 1
	int Number_of_Entity_Keys = fSurfaces.Length();
	int Number_of_Node_Blocks = 1; //TEMP - must be 1

	Node_Block_Types.Dimension(Number_of_Node_Blocks);
	Node_Block_Types[0] = ContactSearch::NODE; //TEMP must be NODE

	Number_Nodes_in_Blocks.Dimension(Number_of_Node_Blocks);
	Number_Nodes_in_Blocks[0] = fStrikerTags.Length(); //TEMP - all NODE in 1
	const int* Node_Global_IDs = fStrikerTags.Pointer(); //all facet nodes as strikers
	int  Number_of_Face_Blocks = fSurfaces.Length();

	Face_Block_Types.Dimension(Number_of_Face_Blocks);
//NOTE: set this based on the number of nodes	
	Face_Block_Types = ContactSearch::TRIFACEL3; //this end has tri facets only
//NOTE: set this based on the number of nodes	

	Number_Faces_in_Blocks.Dimension(Number_of_Face_Blocks);
	for (int i = 0; i < Number_Faces_in_Blocks.Length(); i++)
		Number_Faces_in_Blocks[i] = fSurfaces[i].MajorDim();		

	GenerateACMEConnectivities(Connectivity);

	int  Number_of_Nodal_Comm_Partners = 0; //TEMP - must be 0 for now
	Nodal_Comm_Proc_IDs.Dimension(1); //TEMP - unused
	Number_Nodes_to_Partner.Dimension(1); //TEMP - unused
	Communication_Nodes.Dimension(1); //TEMP - unused
		 	
	/* set up ACME search */
	ContactSearch::ContactErrorCode error;
	fContactSearch = new ContactSearch(Dimensionality,
		Number_of_States,
		Number_of_Entity_Keys,
		Number_of_Node_Blocks,
		Node_Block_Types.Pointer(),
		Number_Nodes_in_Blocks.Pointer(),
		Node_Global_IDs,
		Number_of_Face_Blocks,
		Face_Block_Types.Pointer(),
		Number_Faces_in_Blocks.Pointer(),
		Connectivity.Pointer(),
		Number_of_Nodal_Comm_Partners,
		Nodal_Comm_Proc_IDs.Pointer(),
		Number_Nodes_to_Partner.Pointer(),
		Communication_Nodes.Pointer(),
		mpi_communicator,
		error);
	if (!fContactSearch) throw ExceptionT::kOutOfMemory;
	if (error != ContactSearch::NO_ERROR)
	{
		cout << "\n ACME_Contact3DT::InitContactSearch: ContactSearch constructor\n"
		     <<   "     returned error code: " << error << endl;
		throw ExceptionT::kGeneralFail;	
	}

	/* set search data */
	bool self_contact = false;
	bool two_pass = false;
	double normal_tolerance = 0.1;
	double tangent_tolerance = 0.0;
	fSearchData.Dimension(3*Number_of_Entity_Keys*Number_of_Entity_Keys);
	double* pdata = fSearchData.Pointer();
	for (int i3 = 0; i3 < Number_of_Entity_Keys; i3++)
		for (int i2 = 0; i2 < Number_of_Entity_Keys; i2++)
		{
			/* interaction flag */
			if (two_pass)
			{
				if (i2 == i3)
					*pdata++ = (self_contact) ? 1 : 0;
				else
					*pdata++ = 1;
			}
			else /* one-pass */
			{
				if (self_contact)
					*pdata++ = (i2 >= i3) ? 1 : 0;
				else
					*pdata++ = (i2 > i3) ? 1 : 0;
			}

			/* normal search tolerance */
			*pdata++ = normal_tolerance;

			/* tangential search tolerance */
			*pdata++ = tangent_tolerance;
		}
		
	/* check size */
	error = fContactSearch->Check_Search_Data_Size(3, Number_of_Entity_Keys);
	if (error != ContactSearch::NO_ERROR)
	{
		cout << "\n ACME_Contact3DT::InitContactSearch: Check_Search_Data_Size\n"
		     <<   "     returned error code: " << error << endl;
		throw ExceptionT::kGeneralFail;	
	}
	
	/* set data */
	fContactSearch->Set_Search_Data(fSearchData.Pointer());
	
	/* configure options */
	double dummy = 1.0;
	fContactSearch->Set_Search_Option(ContactSearch::MULTIPLE_INTERACTIONS,
		ContactSearch::INACTIVE, &dummy);
	fContactSearch->Set_Search_Option(ContactSearch::NORMAL_SMOOTHING,
		ContactSearch::INACTIVE, &dummy);
}

/* translate native facet data to ACME format */
void ACME_Contact3DT::GenerateACMEConnectivities(iArrayT& connectivities)
{
	/* dimension */
	int tot_num_facets = 0;
	for (int i = 0; i < fSurfaces.Length(); i++)
		tot_num_facets += fSurfaces[i].MajorDim();
	connectivities.Dimension(tot_num_facets*fNumFacetNodes);

	/* node number -> striker tag map */
	int max, shift;
	fStrikerTags.MinMax(shift, max);
	int range = max - shift + 1;
	iArrayT map(range);
	map = -1;

	/* set map */
	for (int j = 0; j < fStrikerTags.Length(); j++)
	{
		int& index = map[fStrikerTags[j] - shift];
		if (index != -1)
		{
			cout << "\n ACME_Contact3DT::GenerateACMEConnectivities: striker\n"
			     <<   "     tag " << fStrikerTags[j]
			     << " appears more than once in the list of strikers." << endl;
			throw ExceptionT::kGeneralFail;
		}
		else
			index = j;
	}

	/* set 1D connectivities array */
	int* pconnect = connectivities.Pointer();
	for (int k = 0; k < fSurfaces.Length(); k++)
	{
		int* psurf = fSurfaces[k].Pointer();
		int dim = fSurfaces[k].Length();
		for (int i = 0; i < dim; i++)
			*pconnect++ = map[*psurf++ - shift] + 1; // FORTRAN numbering
	}
	
	/* check */
	int min;
	connectivities.MinMax(min, max);
	if (min != 1 || max != fStrikerTags.Length())
	{
		cout << "\n ACME_Contact3DT::GenerateACMEConnectivities: error\n"
		     <<   "     generating ACME connectivies {min, max} = {" << min
		     << ", " << max << "}" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

#endif /* __ACME__ */
