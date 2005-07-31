/* $Id: ModelFileT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (12/15/1999)                                          */

#include "ModelFileT.h"

#include <iostream.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

#include "fstreamT.h"
#include "iArray2DT.h"
#include "ExodusT.h"

/* parameters */
const char  sComment = '#';
const char* sVersion = "1.0";

const char* keyword[] = {"version",
"title",
"dimensions",
"nodes",
"elements",
"nodesets",
"sidesets",
"set"};

enum WordIndex {kversion = 0,
ktitle,
kdimensions,
knodes,
kelements,
knodesets,
ksidesets,
kset};

/* constructor */
ModelFileT::ModelFileT(void):
	fMode(kClosed),
	fNumNodes(0),
	fDimension(0)
{

}

/* destructor */
ModelFileT::~ModelFileT(void) { if (fMode != kClosed) Close(); }

/* translate */
ModelFileT::StatusT ModelFileT::Translate(const ExodusT& exo_file)
{
	/* output file name */
	StringT new_file;
	new_file.Root(exo_file.Filename());
	new_file.Append(".geom");
	OpenWrite(new_file, true);
	
	/* coordinates */
	int nnd = exo_file.NumNodes();
	int nsd = exo_file.NumDimensions();
	dArray2DT coords(nnd, nsd);
	exo_file.ReadCoordinates(coords);
	if (PutCoordinates(coords) != kOK) return kFail;
	coords.Free();
	
	/* element blocks */
	int nblk = exo_file.NumElementBlocks();
	iArrayT blk_id(nblk);
	exo_file.ElementBlockID(blk_id);
	for (int i = 0; i < nblk; i++)
	{
		/* read */
		int nel, nen;
		exo_file.ReadElementBlockDims(blk_id[i], nel, nen);
		GeometryT::CodeT geometry_code;
		iArray2DT connects(nel, nen);
		exo_file.ReadConnectivities(blk_id[i], geometry_code, connects);
	
		/* add */
		if (PutElementSet(blk_id[i], connects) != kOK) return kFail;
	}
	
	/* node sets */
	int nns = exo_file.NumNodeSets();
	iArrayT ns_id(nns);
	exo_file.NodeSetID(ns_id);
	for (int j = 0; j < nns; j++)
	{
		/* read */
		int nsn = exo_file.NumNodesInSet(ns_id[j]);
		iArrayT nodes(nsn);
		exo_file.ReadNodeSet(ns_id[j], nodes);
	
		/* add */
		if (PutNodeSet(ns_id[j], nodes) != kOK) return kFail;
	}

	/* side sets */
	int nss = exo_file.NumSideSets();
	iArrayT ss_id(nss);
	exo_file.SideSetID(ss_id);
	for (int k = 0; k < nss; k++)
	{
		/* read */
		int ssd = exo_file.NumSidesInSet(ss_id[k]);
		iArray2DT sides(ssd, 2);
		int block_ID;
		exo_file.ReadSideSet(ss_id[k], block_ID, sides);
		
		/* add */
		if (PutSideSet(ss_id[k], block_ID, sides) != kOK) return kFail;
	}

	/* close and write */
	Close();
	return kOK;
}

/* open file */
ModelFileT::StatusT ModelFileT::OpenRead(const StringT& file_name)
{
	/* no file open */
	if (fMode != kClosed) return kFail;

	/* see if file exists */	
	ifstreamT tmp(sComment, file_name);
	if (tmp.is_open() && CheckVersion(tmp) == kOK)
	{
		fFileName = file_name;
		fMode = kRead;
		
		/* read file information */
		StatusT status = GetInformation();
		if (status != kOK)
			cout << "\n ModelFileT::OpenRead: GetInformation: error"
			     << endl;
		
		return status;
	}
	else
		return kFail;
}

ModelFileT::StatusT ModelFileT::OpenWrite(const StringT& file_name,
	bool extern_file)
{
	/* no file open */
	if (fMode != kClosed)
		return kFail;
	else
	{
		fFileName = file_name;
		fMode = kWrite;
		fExternFile = extern_file;
		return kOK;
	}
}

/* close */
void ModelFileT::Close(void)
{
	if (fMode == kWrite)
	{
		/* write data to file */
		WriteFile(fExternFile);

		/* set dimensions */
		fNumNodes = 0;
		fDimension = 0;

		/* free all memory */
		fCoordinates.Free();
		fElementID.Free();
		fNodeSetID.Free();
		fSideSetID.Free();
		for (int i = 0; i < fElementSets.Length(); i++)
			delete fElementSets[i];
		for (int j = 0; j < fNodeSets.Length(); j++)
			delete fNodeSets[j];
		for (int k = 0; k < fSideSets.Length(); k++)
			delete fSideSets[k];
	}

	fMode = kClosed;
}

/* title */
ModelFileT::StatusT ModelFileT::PutTitle(const StringT& title)
{
	if (fMode != kWrite)
		return kFail;
	else
	{
		fTitle = title;
		return kOK;
	}
}

ModelFileT::StatusT ModelFileT::GetTitle(StringT& title) const
{
	if (fMode != kRead)
		return kFail;
	else
	{
		ifstreamT in(sComment, fFileName);
		if (AdvanceStream(in, keyword[ktitle]) == kOK)
		{
			/* advance */
			in.next_char();
			title.GetLineFromStream(in);
	
			return in.good() ? kOK : kFail;
		}
		else
			return kFail;
	}
}

/* coordinates */
ModelFileT::StatusT ModelFileT::PutCoordinates(const dArray2DT& coords)
{
	if (fMode != kWrite || fCoordinates.MajorDim() != 0)
		return kFail;
	else
	{
		/* store dimensions */
		fNumNodes  = coords.MajorDim();
		fDimension = coords.MinorDim();

		/* write coordinates */
		if (fExternFile)
		{
			StringT extern_file_name(fFileName);
			extern_file_name.Append(".nd");
					
			ofstream ex_out;
			OpenStream(ex_out, extern_file_name);
			ex_out << fNumNodes  << '\n';
			ex_out << fDimension << '\n';
			coords.WriteNumbered(ex_out);
		}
		/* store coordinates */
		else
			fCoordinates = coords;
		
		return kOK; // assume OK
	}
}

ModelFileT::StatusT ModelFileT::GetDimensions(int& num_nodes,
	int& dimension) const
{
	if (fMode != kRead)
		return kFail;
	else
	{
		num_nodes = fNumNodes;
		dimension = fDimension;
		return kOK;
	}
}

ModelFileT::StatusT ModelFileT::GetCoordinates(dArray2DT& coords) const
{
	if (fMode != kRead) return kFail;
	
	ifstreamT in(sComment, fFileName);
	if (AdvanceStream(in, keyword[knodes]) == kOK)
	{
		ifstreamT in2;
		ifstreamT& src = OpenExternal(in, in2, cout, false,
			"ModelFileT::GetCoordinates: file not found");

		int num_nodes, dimension;
		src >> num_nodes >> dimension;
		coords.Allocate(num_nodes, dimension);
		if (coords.Length() > 0) coords.ReadNumbered(src);
		
		return in2.good() ? kOK : kFail;
	}
	else
		return kFail;
}

/* element sets */
ModelFileT::StatusT ModelFileT::PutElementSet(int ID, const iArray2DT& set)
{
	if (fMode != kWrite) return kFail;

	/* ID must be unique */
	int dex;
	if (fElementID.ColumnHasValue(0, ID, dex))
		return kFail;
	else
	{
		/* store set info */
		int dim = fElementID.MajorDim();
		if (dim == 0)
			fElementID.Allocate(1, 3);
		else
			fElementID.Resize(dim + 1);
		
		fElementID(dim, 0) = ID;
		fElementID(dim, 1) = set.MajorDim();
		fElementID(dim, 2) = set.MinorDim();
			
		fElementSets.Resize(dim + 1);

		/* write set */
		if (fExternFile)
		{
			fElementSets[dim] = NULL;

			StringT extern_file_name(fFileName);
			extern_file_name.Append(".es", dim);
			
			ofstream ex_out;
			OpenStream(ex_out, extern_file_name);
			ex_out << set.MajorDim() << '\n';
			ex_out << set.MinorDim() << '\n';
			set.WriteNumbered(ex_out);
		
			return kOK; // assume OK
		}
		else
		{
			fElementSets[dim] = new iArray2DT(set);
			return !fElementSets[dim] ? kFail : kOK;
		}
	}
}

ModelFileT::StatusT ModelFileT::GetElementSetID(iArrayT& ID) const
{
	if (fMode != kRead)
		return kFail;
	else
	{
		ID.Allocate(fElementID.MajorDim());
		fElementID.ColumnCopy(0, ID);
		return kOK;
	}
}

ModelFileT::StatusT ModelFileT::GetElementSetDimensions(int ID,
	int& num_elements, int& dimension) const
{
	int dex;
	if (fMode != kRead || !fElementID.ColumnHasValue(0, ID, dex))
		return kFail;
	else
	{
		num_elements = fElementID(dex, 1);
		dimension    = fElementID(dex, 2);
		return kOK;
	}
}

ModelFileT::StatusT ModelFileT::GetElementSet(int ID, iArray2DT& set) const
{
	int dex;
	if (fMode != kRead || !fElementID.ColumnHasValue(0, ID, dex))
		return kFail;
	else
	{
		/* advance stream to element set */
		ifstreamT in(sComment, fFileName);
		if (AdvanceStreamToSubsection(in, keyword[kelements],
			keyword[kset], dex) == kOK)
		{
			ifstreamT in2;
			ifstreamT& src = OpenExternal(in, in2, cout, false,
				"ModelFileT::GetElementSet: file not found");

			int num_elements;
			int num_element_nodes;
			src >> num_elements >> num_element_nodes;

			/* check */
			int nel, nen;
			GetElementSetDimensions(ID, nel, nen);
			if (nel != num_elements || nen != num_element_nodes)
			{
				cout << "\n ModelFileT::GetElementSet: internal dimension error for ID: "
				     << ID << endl;
				return kFail;
			}
			
			set.Allocate(num_elements, num_element_nodes);
			if (set.Length() > 0) set.ReadNumbered(src);
		
			return src.good() ? kOK : kFail;
		}
		else
			return kFail;
	}
}

/* node sets */
ModelFileT::StatusT ModelFileT::PutNodeSet(int ID, const iArrayT& set)
{
	if (fMode != kWrite) return kFail;

	/* ID must be unique */
	int dex;
	if (fNodeSetID.ColumnHasValue(0, ID, dex))
		return kFail;
	else
	{
		/* store set info */
		int dim = fNodeSetID.MajorDim();
		if (dim == 0)
			fNodeSetID.Allocate(1, 2);
		else
			fNodeSetID.Resize(dim + 1);
			
		fNodeSetID(dim, 0) = ID;
		fNodeSetID(dim, 1) = set.Length();
	
		fNodeSets.Resize(dim + 1);

		/* write set */
		if (fExternFile)
		{
			fNodeSets[dim] = NULL;
		
			StringT extern_file_name(fFileName);
			extern_file_name.Append(".ns", dim);

			ofstream ex_out;
			OpenStream(ex_out, extern_file_name);
			ex_out << set.Length() << '\n';
			ex_out << set.wrap(10) << '\n';
			
			return kOK; // assume OK	
		}
		/* store set */
		else
		{
			fNodeSets[dim] = new iArrayT(set);
			return !fNodeSets[dim] ? kFail : kOK;
		}
	}
}

ModelFileT::StatusT ModelFileT::GetNodeSetID(iArrayT& ID) const
{
	if (fMode != kRead)
		return kFail;
	else
	{
		int num_sets = fNodeSetID.MajorDim();
		ID.Allocate(num_sets);
		if (num_sets > 0) fNodeSetID.ColumnCopy(0, ID);
		return kOK;
	}
}

ModelFileT::StatusT ModelFileT::GetNodeSetDimensions(int ID,
	int& num_nodes) const
{
	int dex;
	if (fMode != kRead || !fNodeSetID.ColumnHasValue(0, ID, dex))
		return kFail;
	else
	{
		num_nodes = fNodeSetID(dex, 1);
		return kOK;
	}
}

ModelFileT::StatusT ModelFileT::GetNodeSet(int ID, iArrayT& set) const
{
	int dex;
	if (fMode != kRead || !fNodeSetID.ColumnHasValue(0, ID, dex))
{
cout << "\n ModelFileT::GetNodeSet: wrong mode or ID not found" << endl;
		return kFail;
}
	else
	{
		/* advance stream to node set */
		ifstreamT in(sComment, fFileName);
		if (AdvanceStreamToSubsection(in, keyword[knodesets],
			keyword[kset], dex) == kOK)
		{
			ifstreamT in2;
			ifstreamT& src = OpenExternal(in, in2, cout, false,
				"ModelFileT::GetNodeSet: file not found");

			int num_nodes;
			src >> num_nodes;
			
			/* check */
			int nnd;
			GetNodeSetDimensions(ID, nnd);
			if (nnd != num_nodes)
			{
				cout << "\n ModelFileT::GetNodeSet: internal dimension error for ID: "
				     << ID << endl;
				return kFail;
			}
			
			set.Allocate(num_nodes);
			if (set.Length() > 0) src >> set;
		
			return src.good() ? kOK : kFail;
		}
		else
			return kFail;
	}
}

ModelFileT::StatusT ModelFileT::GetNodeSets(const iArrayT& ID, iArrayT& set) const
{
	/* get total number of nodes */
	int num_nodes = 0;
	for (int i = 0; i < ID.Length(); i++)
	{
		/* read set size */
		int count;
		if (GetNodeSetDimensions(ID[i], count) != kOK) return kFail;
					
		/* add up */
		num_nodes += count;
	}		

	/* allocate */
	set.Allocate(num_nodes);
	iArrayT tmp_set;
	int count = 0;
	for (int j = 0; j < ID.Length(); j++)
	{
		/* read node set */
		if (GetNodeSet(ID[j], tmp_set) != kOK) return kFail;
					
		/* copy in */
		set.CopyPart(count, tmp_set, 0, tmp_set.Length());
		count += tmp_set.Length();
	}		

	/* falls through if successful */
	return kOK;
}

/* side sets */
ModelFileT::StatusT ModelFileT::PutSideSet(int ID, int element_set_ID,
	const iArray2DT& set)
{
	if (fMode != kWrite) return kFail;

	/* ID must be unique */
	int dex;
	if (fSideSetID.ColumnHasValue(0, ID, dex))
		return kFail;
	else
	{
		/* store set info */
		int dim = fSideSetID.MajorDim();
		if (dim == 0)
			fSideSetID.Allocate(1, 3);
		else
			fSideSetID.Resize(dim + 1);

		fSideSetID(dim, 0) = ID;
		fSideSetID(dim, 1) = element_set_ID;
		fSideSetID(dim, 2) = set.MajorDim();
		
		fSideSets.Resize(dim + 1);

		/* write set */
		if (fExternFile)
		{
			fSideSets[dim] = NULL;
		
			StringT extern_file_name(fFileName);
			extern_file_name.Append(".ss", dim);

			ofstream ex_out;
			OpenStream(ex_out, extern_file_name);
			ex_out << set.MajorDim() << '\n';
			ex_out << set << '\n';
			
			return kOK; // assume OK
		}
		/* store set */
		else
		{
			fSideSets[dim] = new iArray2DT(set);
			return !fSideSets[dim] ? kFail : kOK;
		}
	}
}

ModelFileT::StatusT ModelFileT::GetSideSetID(iArrayT& ID) const
{
	if (fMode != kRead)
		return kFail;
	else
	{
		int num_sets = fSideSetID.MajorDim();
		ID.Allocate(num_sets);
		if (num_sets > 0) fSideSetID.ColumnCopy(0, ID);
		return kOK;
	}
}

ModelFileT::StatusT ModelFileT::GetSideSetDimensions(int ID,
	int& num_sides) const
{
	int dex;
	if (fMode != kRead || !fSideSetID.ColumnHasValue(0, ID, dex))
		return kFail;
	else
	{
		num_sides = fSideSetID(dex, 2);
		return kOK;
	}
}

ModelFileT::StatusT ModelFileT::GetSideSet(int ID, int& element_set_ID,
	iArray2DT& set) const
{
	int dex;
	if (fMode != kRead || !fSideSetID.ColumnHasValue(0, ID, dex))
		return kFail;
	else
	{
		/* advance stream to node set */
		ifstreamT in(sComment, fFileName);
		if (AdvanceStreamToSubsection(in, keyword[ksidesets],
			keyword[kset], dex) == kOK)
		{
			element_set_ID = fSideSetID(dex, 1);
	
			ifstreamT in2;
			ifstreamT& src = OpenExternal(in, in2, cout, false,
				"ModelFileT::GetSideSet: file not found");

			int num_sides;
			src >> num_sides;
			
			/* check */
			int ns;
			GetSideSetDimensions(ID, ns);
			if (ns != num_sides)
			{
				cout << "\n ModelFileT::GetSideSet: internal dimension error for ID: "
				     << ID << endl;
				return kFail;
			}
			
			set.Allocate(num_sides, 2);
			if (set.Length() > 0) src >> set;
		
			return src.good() ? kOK : kFail;
		}
		else
			return kFail;
	}
}

/**************************************************************************
* Private
**************************************************************************/

/* return 1 if version is current */
ModelFileT::StatusT ModelFileT::CheckVersion(ifstreamT& in) const
{
	StatusT status = AdvanceStream(in, keyword[kversion]);
	if (status == kOK)
	{
		StringT version;
		in >> version;
		return version == sVersion ? kOK : kFail;
	}
	else
		return kFail;
}

/* advance to line after next occurence of key */
ModelFileT::StatusT ModelFileT::AdvanceStream(istream& in,
	const char* key) const
{
	int found = 0;
	char word[255], line[255];
	while (!found && in.good())
	{
		in >> word;
		if (word[0] == '*')
		{
			ToLower(word);
			if (strcmp(word + 1, key) == 0)
				found = 1;
		}
		else
			in.getline(line, 254);
	}
	
	return found ? kOK : kFail;
}

ModelFileT::StatusT ModelFileT::AdvanceStreamToSubsection(istream& in,
	const char* section, const char* subsection, int index) const
{
	if (AdvanceStream(in, section) == kOK)
	{
		int i = 0;
		StatusT status = AdvanceStream(in, subsection);
		while (status == kOK && i < index)
		{
			i++;
			status = AdvanceStream(in, subsection);
		}
		
		return status;
	}
	else
		return kFail;
}	

/* get set information from file */
ModelFileT::StatusT ModelFileT::GetInformation(void)
{
	ifstreamT in(sComment, fFileName);
	StatusT status = AdvanceStream(in, keyword[kdimensions]);
	if (status == kOK)
	{
		in >> fNumNodes;
		in >> fDimension;
		
		int num_elem_sets;
		in >> num_elem_sets;
		fElementID.Allocate(num_elem_sets, 3);
		if (fElementID.Length() > 0) in >> fElementID;

		int num_node_sets;
		in >> num_node_sets;
		fNodeSetID.Allocate(num_node_sets, 2);
		if (fNodeSetID.Length() > 0) in >> fNodeSetID;

		int num_side_sets;
		in >> num_side_sets;
		fSideSetID.Allocate(num_side_sets, 3);
		if (fSideSetID.Length() > 0) in >> fSideSetID;
		
		return in.good() ? kOK : kFail;
	}
	else
		return kFail;
}

/* write data to file */
void ModelFileT::WriteFile(bool extern_file) const
{
	/* open and format stream */
	ofstream out;
	OpenStream(out, fFileName);

	out << "*" << keyword[kversion] << '\n';
	out << sVersion << '\n';

	out << "*" << keyword[ktitle] << '\n';
	out << fTitle << '\n';

	out << "*" << keyword[kdimensions] << '\n';
	out << fNumNodes << '\n';
	out << fDimension << '\n';
	out << fElementID.MajorDim() << " # number of element sets\n";
	out << fElementID << '\n';
	out << fNodeSetID.MajorDim() << " # number of node sets\n";
	out << fNodeSetID << '\n';
	out << fSideSetID.MajorDim() << " # number of side sets\n";
	out << fSideSetID << '\n';

	out << "*" << keyword[knodesets] << '\n';
	for (int i = 0; i < fNodeSets.Length(); i++)
	{
		out << "*" << keyword[kset] << '\n';
		
		const iArrayT& set = *(fNodeSets[i]);
		if (extern_file)
		{
			StringT extern_file_name(fFileName);
			extern_file_name.Append(".ns", i);
			out << extern_file_name << '\n';			
		}		
		else
		{
			out << set.Length() << '\n';
			out << set.wrap(10) << '\n';
		}
	}

	out << "*" << keyword[ksidesets] << '\n';
	for (int j = 0; j < fSideSets.Length(); j++)
	{
		out << "*" << keyword[kset] << '\n';
		
		const iArray2DT& set = *(fSideSets[j]);
		if (extern_file)
		{
			StringT extern_file_name(fFileName);
			extern_file_name.Append(".ss", j);
			out << extern_file_name << '\n';			
		}		
		else
		{
			out << set.MajorDim() << '\n';
			out << set << '\n';
		}
	}

	out << "*" << keyword[kelements] << '\n';
	for (int k = 0; k < fElementSets.Length(); k++)
	{
		out << "*" << keyword[kset] << '\n';
		
		const iArray2DT& set = *(fElementSets[k]);
		if (extern_file)
		{
			StringT extern_file_name(fFileName);
			extern_file_name.Append(".es", k);
			out << extern_file_name << '\n';
		}		
		else
		{
			out << set.MajorDim() << '\n';
			out << set.MinorDim() << '\n';
			set.WriteNumbered(out);
		}
	}

	out << "*" << keyword[knodes] << '\n';
	if (extern_file)
	{
		StringT extern_file_name(fFileName);
		extern_file_name.Append(".nd");
		out << extern_file_name << '\n';
	}		
	else
	{
		out << fCoordinates.MajorDim() << '\n';
		out << fCoordinates.MinorDim() << '\n';
		fCoordinates.WriteNumbered(out);
	}
}

ifstreamT& ModelFileT::OpenExternal(ifstreamT& in,  ifstreamT& in2,
	ostream& out, bool verbose, const char* fail) const
{
	/* check for external file */
	char nextchar = in.next_char();
	if (isdigit(nextchar))
		return in;
	else
	{
		/* open external file */
		StringT file;
		in >> file;
		if (verbose) out << " external file: " << file << '\n';
			
		/* open stream */
		file.ToNativePathName();
		in2.open(file);
		if (!in2.is_open())
		{
			if (fail) cout << fail << ": " << file << endl;
			throw eBadInputValue;
		}

		/* set comments */
		if (in.skip_comments()) in2.set_marker(in.comment_marker());

		return in2;
	}
}

/* open output file */
ostream& ModelFileT::OpenStream(ofstream& out, const StringT& file_name) const
{
	//out.close(); // error to close more than once???
	out.open(file_name);
	
	/* format */
	out.precision(DBL_DIG);
	out.setf(ios::showpoint);
	out.setf(ios::right, ios::adjustfield);
	out.setf(ios::scientific, ios::floatfield);
	
	return out;
}

/* convert string to lower case */
void ModelFileT::ToLower(char* str) const
{
	int count = 0;
	while (*str != '\0' && ++count < 255)
	{
	    *str = tolower(*str);
		str++;
	}
}	
