/* $Id: PartitionT.cpp,v 1.12 2004-06-17 06:38:27 paklein Exp $ */
/* created: paklein (11/16/1999) */
#include "PartitionT.h"

#include "ifstreamT.h"
#include "GraphT.h"
#include "AutoArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "RaggedArray2DT.h"
#include "StringT.h"

using namespace Tahoe;

/* parameters */
const int kHeadRoom = 20;                // percent
const char* sPartitionTVersion = "v1.1"; // version marker

/* constructor */
PartitionT::PartitionT(void):
	fNumPartitions(0),
	fID(-1),
	fScope(kUnSet),
	fDecompType(kUndefined)
{

}

/* returns true if version if current */
bool PartitionT::CheckVersion(const StringT& version)
{
	return version == sPartitionTVersion;
}

/* nodes assigned to the partition  - internal + border */
void PartitionT::PartitionNodes(iArrayT& nodes, NumberScopeT scope) const
{
	/* allocate */
	nodes.Dimension(fNodes_i.Length() + fNodes_b.Length());

	/* partition nodes in local numbering */
	nodes.CopyPart(0, fNodes_i, 0, fNodes_i.Length());
	nodes.CopyPart(fNodes_i.Length(), fNodes_b, 0, fNodes_b.Length());

	/* map to global numbers */
	if (scope == kGlobal) SetNodeScope(kGlobal, nodes);
}

/* communication nodes */
const iArrayT* PartitionT::NodesIn(int ID) const
{
	int ID_dex;
	if (fCommID.HasValue(ID, ID_dex))
		return fNodes_in.Pointer(ID_dex);
	else
		return NULL;
}

const iArrayT* PartitionT::NodesOut(int ID) const
{
	int ID_dex;
	if (fCommID.HasValue(ID, ID_dex))
		return fNodes_out.Pointer(ID_dex);
	else
		return NULL;
}

/* set node info */
void PartitionT::Set(int num_parts, int id, const iArrayT& part_map, 
	const GraphT& graph)
{
	/* total number of partitions */
	fNumPartitions = num_parts;
	if (fNumPartitions < 1) throw ExceptionT::kGeneralFail;

	/* set ID */
	fID = id;
	if (fID < 0 || fID >= fNumPartitions) throw ExceptionT::kOutOfRange;
		
	/* numbering is global */
	fScope = kGlobal;

	/* resolve internal/boundary nodes */
	ClassifyNodes(part_map, graph);
	
	/* set send information */
	SetReceive(part_map);
	
	/* clear inverse maps */
	fInvNodeMap.Free();
	for (int i = 0; i < fInvElementMap.Length(); i++)
		fInvElementMap.Free();
}

void PartitionT::Set(int num_parts, int id, const iArrayT& part_map, const ArrayT<const iArray2DT*>& connects_1,
	const ArrayT<const RaggedArray2DT<int>*>& connects_2)
{
	/* total number of partitions */
	fNumPartitions = num_parts;
	if (fNumPartitions < 1) throw ExceptionT::kGeneralFail;

	/* set ID */
	fID = id;
	if (fID < 0 || fID >= fNumPartitions) throw ExceptionT::kOutOfRange;
		
	/* numbering is global */
	fScope = kGlobal;

	/* resolve internal/boundary nodes */
	ClassifyNodes(part_map, connects_1, connects_2);
	
	/* set send information */
	SetReceive(part_map);
	
	/* clear inverse maps */
	fInvNodeMap.Free();
	for (int i = 0; i < fInvElementMap.Length(); i++)
		fInvElementMap.Free();
}

/* store outgoing data (count maps onto fCommID) */
void PartitionT::SetOutgoing(const ArrayT<iArrayT>& nodes_out)
{
	/* checks */
	if (nodes_out.Length() != fCommID.Length())
	{
		cout << "\n PartitionT::SetIncoming: expecting dimension of outgoing ID's ("
		     << nodes_out.Length() << ")\n";
		cout <<   "     to be the same as the Comm ID list (" << fCommID.Length()
		     << ")" << endl;
		throw ExceptionT::kSizeMismatch;
	}
	
	/* only at global scope for now */
	if (fScope != kGlobal)
	{
		cout << "\n PartitionT::SetOutgoing: number scope must be global" << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* store */
	fNodes_out.Dimension(nodes_out.Length());
	for (int i = 0; i < fNodes_out.Length(); i++)
		fNodes_out[i] = nodes_out[i];
}

void PartitionT::SetScope(NumberScopeT scope)
{
	/* quick exit */
	if (scope == fScope || scope == kUnSet) return;

	/* convert nodal data */
	int shift;
	iArrayT node_map;
	SetNodeMap(scope, node_map, shift);

	MapValues(node_map, shift, fNodes_i);
	MapValues(node_map, shift, fNodes_b);		
	MapValues(node_map, shift, fNodes_e);

	for (int k = 0; k < fNodes_in.Length(); k++)
		MapValues(node_map, shift, fNodes_in[k]);

	for (int j = 0; j < fNodes_out.Length(); j++)
		MapValues(node_map, shift, fNodes_out[j]);

	/* convert element data */
	for (int i = 0; i < fElementMap.Length(); i++)
		if (fElementMap[i].Length() > 0)
		{
			int shift;
			iArrayT element_map;
			SetElementMap(scope, fElementBlockID[i], element_map, shift);
		
			MapValues(element_map, shift, fElements_i[i]);
			MapValues(element_map, shift, fElements_b[i]);
		}

	/* set scope */
	fScope = scope;
}

void PartitionT::InitElementBlocks(const ArrayT<StringT>& blockID)
{
	/* copy */
	fElementBlockID = blockID;
	
	/* allocate memory */
	int n = fElementBlockID.Length();
	fElements_i.Dimension(n);
	fElements_b.Dimension(n);
	fElementMap.Dimension(n);
	fElementMapShift.Dimension(n);
	fInvElementMap.Dimension(n);
}

/* collect internal and border elements */
void PartitionT::SetElements(const StringT& blockID, const iArray2DT& connects)
{
	/* only at global scope for now */
	if (fScope != kGlobal)
	{
		cout << "\n PartitionT::SetElements: number scope must be global"
		     << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* resolve ID */
	int dex = ElementBlockIndex(blockID, "SetElements");

	/* node range */
	int min, max;
	connects.MinMax(min, max);

	/* mark nodes as internal or border */
	int range = max - min + 1;
	ArrayT<StatusT> node_map(range);
	node_map = kExternal;
	MapStatus(kInternal, fNodes_i, node_map, min);
	MapStatus(  kBorder, fNodes_b, node_map, min);
	
	//TEMP
	MapStatus(  kBorder, fNodes_e, node_map, min); // temp fix 1
	
	AutoArrayT<int> elements_i(20);
	AutoArrayT<int> elements_b(20);
	
	/* sort elements into internal/external */
	int nel = connects.MajorDim();
	int nen = connects.MinorDim();
	const int* pelem = connects.Pointer();
	for (int i = 0; i < nel; i++)
	{
		int n_i = 0;
		int n_b = 0;
		int n_e = 0;
		for (int j = 0; j < nen; j++)
		{
			int dex = *pelem++ - min;
			if (dex >= 0 && dex <= max)
			{
				StatusT status = node_map[dex];
				if (status == kInternal)
					n_i++;
				else if (status == kBorder) // counting border or external nodes
					n_b++;
				else
					n_e++;
			}
			else
				n_e++;
		}
	
		//if (n_i > 0) // does not reconstruct "overlap" region
		if (n_e == 0) // temp fix 1
		{
			if (n_b > 0)
				elements_b.Append(i);
			else
				elements_i.Append(i);			
		}
// if connectsX where labelled by connectsU
#if 0
		if (n_i != 0 || n_b != 0)
		{
			if (n_e > 0)
				elements_b.Append(i);
			else
				elements_i.Append(i);			
		}
#endif

	}

	/* copy in */
	int n_i = elements_i.Length();
	fElements_i[dex].Dimension(n_i);
	elements_i.CopyInto(fElements_i[dex]);

	int n_b = elements_b.Length();
	fElements_b[dex].Dimension(n_b);
	elements_b.CopyInto(fElements_b[dex]);
	
	/* set map */
	fElementMap[dex].Dimension(n_i + n_b);
	fElementMap[dex].CopyPart(0, fElements_i[dex], 0, n_i);
	fElementMap[dex].CopyPart(n_i, fElements_b[dex], 0, n_b);
}

/* check cross-references - returns 1 if OK */
int PartitionT::CrossCheck(const PartitionT& that) const
{
	int i1, i2;
	int q1, q2;

	/* receiving */
	q1 = fCommID.HasValue(that.fID, i1);
	q2 = (that.fCommID).HasValue(fID, i2);
	if (q1 != q2)
		return 0;
	else if (q1)
	{
		int n1 = fNodes_in[i1].Length();
		int n2 = (that.fNodes_out[i2]).Length();
		if (n1 != n2)
			return 0;
		else
		{
			const int* p1 = fNodes_in[i1].Pointer();
			const int* p2 = (that.fNodes_out[i2]).Pointer();
			for (int i = 0; i < n1; i++)
				if (fNodeMap[*p1++] != (that.fNodeMap)[*p2++]) return 0;
				// need to check global/local here
		}
	}

	/* sending */
	q1 = fCommID.HasValue(that.fID, i1);
	q2 = (that.fCommID).HasValue(fID, i2);
	if (q1 != q2)
		return 0;
	else if (q1)
	{
		int n1 = fNodes_out[i1].Length();
		int n2 = (that.fNodes_in[i2]).Length();
		if (n1 != n2)
			return 0;
		else
		{
			const int* p1 = fNodes_out[i1].Pointer();
			const int* p2 = (that.fNodes_in[i2]).Pointer();
			for (int i = 0; i < n1; i++)
				if (fNodeMap[*p1++] != (that.fNodeMap)[*p2++]) return 0;
		}
	}

	/* true on fall through */
	return 1;
}

namespace Tahoe {

/* I/O */
ostream& operator<<(ostream& out, const PartitionT& partition)
{
	out << sPartitionTVersion << '\n';
	out << setw(6) << partition.fNumPartitions << " # number of parts\n";
	out << setw(6) << partition.fID            << " # partition number\n";
	out << setw(6) << partition.fScope         << " # numbering scope\n";
	out << setw(6) << partition.fDecompType    << " # decomposition type\n";
	
	// nodal information
	out << "# internal nodes:\n";
	out << (partition.fNodes_i).Length() << '\n';
	out << (partition.fNodes_i).wrap_tight(10) << '\n'; // internal nodes
	
	out << "# border nodes:\n";
	out << (partition.fNodes_b).Length() << '\n';
	out << (partition.fNodes_b).wrap_tight(10) << '\n'; // border nodes	
	
	out << "# external nodes:\n";
	out << (partition.fNodes_e).Length()  << '\n';
	out << (partition.fNodes_e).wrap_tight(10) << '\n'; // external nodes
	
	// receive/send information
	out << "# comm ID list:\n";
	out << (partition.fCommID).Length() << '\n';
	out << (partition.fCommID).wrap_tight(10) << '\n'; // ID's of communicating partitions

	out << "# incoming nodes:\n";
	for (int i = 0; i < (partition.fNodes_in).Length(); i++)
	{
		out << (partition.fNodes_in)[i].Length() << '\n';
		out << (partition.fNodes_in)[i].wrap_tight(10) << '\n';
	}

	out << "# outgoing nodes:\n";
	for (int j = 0; j < (partition.fNodes_out).Length(); j++)
	{
		out << (partition.fNodes_out)[j].Length() << '\n';
		out << (partition.fNodes_out)[j].wrap_tight(10) << '\n';
	}
	
	// element information
	out << "# number of element blocks:\n";
	out << (partition.fElements_i).Length() << '\n';

	out << "# element block ID's:\n";
	int wrap = 0;
	for (int i = 0; i < partition.fElementBlockID.Length(); i++)
	{
		if (wrap++ == 10) {
			out << '\n';
			wrap = 0;
		}
		out << partition.fElementBlockID[i] << " ";
	}
	out << '\n';

	out << "# internal elements (by block):\n";
	for (int k = 0; k < (partition.fElements_i).Length(); k++)
	{
		out << (partition.fElements_i)[k].Length() << '\n';
		out << (partition.fElements_i)[k].wrap_tight(10) << '\n'; // internal elements
	}
	
	out << "# external elements (by block):\n";
	for (int l = 0; l < (partition.fElements_b).Length(); l++)
	{
		out << (partition.fElements_b)[l].Length() << '\n';
		out << (partition.fElements_b)[l].wrap_tight(10) << '\n'; // internal elements
	}

	// global node map
	out << "# node maps:\n";
	out << (partition.fNodeMap).Length() << '\n';
	out << (partition.fNodeMap).wrap_tight(10) << '\n'; // global[local]

	// block global element numbering map
	out << "# element maps (by block):\n";
	for (int m = 0; m < (partition.fElementMap).Length(); m++)
	{
		out << (partition.fElementMap)[m].Length() << '\n';
		out << (partition.fElementMap)[m].wrap_tight(10) << '\n';
	}

	return out;
}

}

/* operator support */
ifstreamT& PartitionT::Read(ifstreamT& in)
{
	/* set comment marker */
	char old_marker = in.comment_marker();
	in.set_marker(CommentMarker());

	/* check version */
	StringT version;
	in >> version;
	if (version != sPartitionTVersion)
		ExceptionT::BadInputValue("operator>>PartitionT", 
			"file version %s is not the current %s", version.Pointer(), sPartitionTVersion);

	in >> fNumPartitions; // number of parts
	in >> fID;            // partition number
	in >> fScope;         // numbering scope
	in >> fDecompType;    // decomposition type

	int length;
	
	// nodal information
	in >> length;
	fNodes_i.Dimension(length);
	in >> fNodes_i; // internal nodes	
	
	in >> length;
	fNodes_b.Dimension(length);
	in >> fNodes_b; // border nodes	

	in >> length;
	fNodes_e.Dimension(length);
	in >> fNodes_e; // external nodes
	
	// receive/send information
	in >> length;
	fCommID.Dimension(length);
	in >> fCommID; // ID's of communicating partitions

	fNodes_in.Dimension(length);
	for (int i = 0; i < length; i++)
	{
		int dim;
		in >> dim;
		fNodes_in[i].Dimension(dim);
		in >> fNodes_in[i];
	}

	fNodes_out.Dimension(length);
	for (int j = 0; j < length; j++)
	{
		int dim;
		in >> dim;
		fNodes_out[j].Dimension(dim);
		in >> fNodes_out[j];
	}
	
	// element information
	in >> length;
	fElementBlockID.Dimension(length);
	for (int i = 0; i < fElementBlockID.Length(); i++)
		in >> fElementBlockID[i];
	InitElementBlocks(fElementBlockID);	

	for (int k = 0; k < fElements_i.Length(); k++)
	{
		int dim;
		in >> dim;
		fElements_i[k].Dimension(dim);
		in >> fElements_i[k]; // internal elements
	}

	for (int l = 0; l < fElements_b.Length(); l++)
	{
		int dim;
		in >> dim;
		fElements_b[l].Dimension(dim);
		in >> fElements_b[l]; // internal elements
	}
	
	// global node map
	in >> length;
	fNodeMap.Dimension(length);
	in >> fNodeMap; // global[local]
	if (length != (fNodes_i.Length() +
	               fNodes_b.Length() +
                   fNodes_e.Length())) throw ExceptionT::kBadInputValue;

	// block global element numbering map
	for (int m = 0; m < fElementMap.Length(); m++)
	{
		int dim;
		in >> dim;
		fElementMap[m].Dimension(dim);
		in >> fElementMap[m];
	}

	/* restore comment marker */
	if (in.skip_comments())
		in.set_marker(old_marker);
	else
		in.clear_marker();
	
	return in;
}

/* resolve element block ID to index */
int PartitionT::ElementBlockIndex(const StringT& blockID, const char* caller) const
{
	int dex = -1;
	for (int i = 0; dex == -1 && i < fElementBlockID.Length(); i++)
		if (fElementBlockID[i] == blockID)
			dex = i;

	if (dex == -1)
	{
		const char* this_routine = "ElementBlockIndex";
		const char* str = (caller != NULL) ? caller : this_routine;
		cout << "\n PartitionT::" << str << ": block ID not found: "
		     << blockID << endl;
		throw ExceptionT::kGeneralFail;
	}
	return dex;
}

/* number transformations */
void PartitionT::MapValues(const iArrayT& map, int shift, ArrayT<int>& values) const
{
	int  n_map = map.Length();
	int  nn = values.Length();
	int* pn = values.Pointer();
	for (int i = 0; i < nn; i++)
	{
		//TEMP ?
		int value = *pn - shift;
		if (value < 0 || value >= n_map)
		{
			cout << "\n PartitionT::MapValues: value " << value
			     << " at position " << i << " is out of\n"
			     <<   "     range {0," << n_map << "}" << endl;
			throw ExceptionT::kOutOfRange;
		}
	
		*pn = map[value];
		pn++;
	}
}

const iArrayT& PartitionT::InverseNodeMap(int& index_shift) const
{
	/* construct inverse node map */
	if (fInvNodeMap.Length() == 0)
	{
		/* cast away const-ness */
		PartitionT* tmp = (PartitionT*) this;
		MakeInverseMap(tmp->fNodeMap, tmp->fInvNodeMap, tmp->fNodeMapShift);	
	}
	index_shift = fNodeMapShift;
	return fInvNodeMap;
}

/* returns indeces of global nodes that lie within the partition */
void PartitionT::ReturnPartitionNodes(const iArrayT& global_nodes,
	iArrayT& partition_indices) const
{
	/* make inverse map */
	iArrayT inv_map;
	int shift;
	MakeInverseMap(fNodeMap, inv_map, shift);

	AutoArrayT<int> tmp(20);
	for (int i = 0; i < global_nodes.Length(); i++)
	{
		int snode = global_nodes[i] - shift;
		if (snode > -1 && snode < inv_map.Length())
			if (inv_map[snode] > -1)
				tmp.Append(i);

	}

	/* copy to return value */
	partition_indices.Dimension(tmp.Length());
	tmp.CopyInto(partition_indices);
}

/* returns indeces of (block) global elements that lie within
* the partition */
void PartitionT::ReturnPartitionElements(const StringT& blockID,
	const iArrayT& global_elements, iArrayT& partition_indices) const
{
	/* make inverse map */
	iArrayT inv_map;
	int shift;
	MakeInverseMap(ElementMap(blockID), inv_map, shift);

	AutoArrayT<int> tmp(20);
	for (int i = 0; i < global_elements.Length(); i++)
	{
		int selement = global_elements[i] - shift;
		if (selement > -1 && selement < inv_map.Length())
			if (inv_map[selement] > -1)
				tmp.Append(i);
	}

	/* copy to return value */
	partition_indices.Dimension(tmp.Length());
	tmp.CopyInto(partition_indices);
}

/* return the node to processor map */
void PartitionT::ReturnProcessorMap(ArrayT<int>& n2p) const
{
	/* dimension/initialize */
	n2p.Dimension(fNodeMap);
	n2p = fID;
	
	/* label nodes from elsewhere */
	for (int i = 0; i < fCommID.Length(); i++)
	{
		int proc = fCommID[i];
		const iArrayT& nodes = fNodes_in[i];
		for (int j = 0; j < nodes.Length(); j++)
			n2p[nodes[j]] = proc;
	}
}

/* mapping functions (assumes scope is currently the opposite) */
void PartitionT::SetNodeScope(NumberScopeT scope, ArrayT<int>& nodes) const
{
	/* quick exit */
	if (nodes.Length() == 0) return;

	/* select map */
	int shift;
	iArrayT map;
	SetNodeMap(scope, map, shift);

	/* apply map */
	MapValues(map, shift, nodes);
}

void PartitionT::SetElementScope(NumberScopeT scope, const StringT& blockID, ArrayT<int>& elements) const
{
	/* quick exit */
	if (elements.Length() == 0) return;

	/* select map */
	int shift;
	iArrayT map;
	SetElementMap(scope, blockID, map, shift);
	
	/* apply map */
	MapValues(map, shift, elements);
}

namespace Tahoe {

/* input operator for scope */
istream& operator>>(istream& in, PartitionT::NumberScopeT& scope)
{
	int i_scope;
	in >> i_scope;
	switch (i_scope)
	{
		case PartitionT::kUnSet:
			scope = PartitionT::kUnSet;
			break;
		case PartitionT::kLocal:
			scope = PartitionT::kLocal;
			break;
		case PartitionT::kGlobal:
			scope = PartitionT::kGlobal;
			break;
		default:
			ExceptionT::BadInputValue("operator>>PartitionT::NumberScopeT", 
				"unrecognized scope: %d", i_scope);	
	}
	return in;
}

/* input operator for scope */
istream& operator>>(istream& in, PartitionT::DecompTypeT& t)
{
	int i_t;
	in >> i_t;
	switch (i_t)
	{
		case PartitionT::kUndefined:
			t = PartitionT::kUndefined;
			break;
		case PartitionT::kGraph:
			t = PartitionT::kGraph;
			break;
		case PartitionT::kAtom:
			t = PartitionT::kAtom;
			break;
		case PartitionT::kSpatial:
			t = PartitionT::kSpatial;
			break;
		default:
			ExceptionT::BadInputValue("operator>>PartitionT::DecompTypeT", 
				"unrecognized scope: %d", i_t);	
	}
	return in;
}

} // namespace Tahoe

/************************************************************************
* Private
************************************************************************/

/* make inverse map (filled with -1) */
void PartitionT::MakeInverseMap(const iArrayT& map, iArrayT& inv_map,
	int& shift) const
{
  if (map.Length() == 0)
	{
	  shift = 0;
	  inv_map.Dimension(0);
	}
  else
	{
	/* range */
	int max;
	map.MinMax(shift, max);
	int range = max - shift + 1;
	
	/* dimension */
	inv_map.Dimension(range);
	inv_map = -1;

	/* make map */
	int dim = map.Length();
	for (int i = 0; i < dim; i++)
		inv_map[map[i] - shift] = i;
	}
}

/* set node info */
void PartitionT::ClassifyNodes(const iArrayT& part_map,
	const GraphT& graph)
{
	/* work space */
	AutoArrayT<int> nodes_i(kHeadRoom);
	AutoArrayT<int> nodes_b(kHeadRoom);
	AutoArrayT<int> nodes_e(kHeadRoom);
	AutoArrayT<int> commID(kHeadRoom);

	/* resolve internal/boundary nodes */
	int nnd = part_map.Length();
	iArrayT mixed(nnd);
	mixed = 0;
	int* pmix = mixed.Pointer();
	for (int i = 0; i < nnd; i++)
	{
		int part_i = part_map[i];
		if (part_i == fID)
		{
			int degree = graph.Degree(i);
			const int *pedge = graph.Edges(i);
			for (int j = 0; j < degree; j++)
			{
				int part_j = part_map[*pedge];
				if (part_j != part_i)
				{
					*pmix = 1;
					nodes_e.AppendUnique(*pedge);
					commID.AppendUnique(part_j);
				}
				pedge++;
			}

			if (*pmix)
				nodes_b.Append(i);
			else
				nodes_i.Append(i);
		}
			
		pmix++;
	}
	
	/* store */
	fNodes_i.Dimension(nodes_i.Length());
	nodes_i.CopyInto(fNodes_i);

	fNodes_b.Dimension(nodes_b.Length());
	nodes_b.CopyInto(fNodes_b);	

	fNodes_e.Dimension(nodes_e.Length());
	nodes_e.CopyInto(fNodes_e);

	fCommID.Dimension(commID.Length());
	commID.CopyInto(fCommID);

#if 0
//NOTE: the other version of ClassifyNodes() inherently sorts the node
//      lists in ascending order and this step allows the two approaches
//      to be verified against each other.
	cout << "\n PartitionT::ClassifyNodes: NOTE: sorting node lists in ascending order" << endl;
	fNodes_i.SortAscending();
	fNodes_b.SortAscending();
	fNodes_e.SortAscending();
	fCommID.SortAscending();
#endif
	
	/* generate node map (just number sequentially through _i, _b, _e) */
	fNodeMap.Dimension(fNodes_i.Length() +
	                  fNodes_b.Length() +
	                  fNodes_e.Length()); // sets sequence for local node
	                                      // numbers - DO NOT CHANGE
	
	/* copy in (with no resequencing) */
	fNodeMap.CopyPart(0, fNodes_i, 0, fNodes_i.Length());
	fNodeMap.CopyPart(fNodes_i.Length(), fNodes_b, 0, fNodes_b.Length());
	fNodeMap.CopyPart(fNodes_i.Length() + fNodes_b.Length(), fNodes_e, 0, fNodes_e.Length());
}

void PartitionT::ClassifyNodes(const iArrayT& part_map, const ArrayT<const iArray2DT*>& connects_1,
	const ArrayT<const RaggedArray2DT<int>*>& connects_2)
{
	/* node classification */
	int nnd = part_map.Length();
	ArrayT<StatusT> status(nnd);
	status = kUnset;
	
	/* mark all partition nodes as internal */
	for (int i = 0; i < nnd; i++)
		if (part_map[i] == fID)
			status[i] = kInternal;

	/* run through connectivities in group 1 */
	for (int j = 0; j < connects_1.Length(); j++)
	{
		/* set dimensions */
		const iArray2DT& connects = *(connects_1[j]);
		int nel = connects.MajorDim();
		int nen = connects.MinorDim();
	
		/* scan elements */
		const int* pel = connects.Pointer();
		for (int i = 0; i < nel; i++)
		{
			/* look for mixed elements */
			bool has_internal = false;
			bool has_external = false;
			for (int a = 0; a < nen && (!has_internal || !has_external); a++)
			{
				int node = pel[a];
				if (status[node] == kInternal || status[node] == kBorder)
					has_internal = true;
				else
					has_external = true;
			}
		
			/* re-mark nodes */
			if (has_internal && has_external)
			{
				for (int a = 0; a < nen; a++)
				{
					int node = pel[a];
					if (status[node] == kInternal)
						status[node] = kBorder;
					else if (status[node] == kUnset)
						status[node] = kExternal;
				}
			}
		
			/* next */
			pel += nen;
		}
	}

	/* run through connectivities in group 1 */
	for (int k = 0; k < connects_2.Length(); k++)
	{
		/* set dimensions */
		const RaggedArray2DT<int>& connects = *(connects_2[k]);
		int nel = connects.MajorDim();
	
		/* scan elements */
		for (int i = 0; i < nel; i++)
		{
			/* current element */
			int nen = connects.MinorDim(i);
			const int* pel = connects(i);

			/* look for mixed elements */
			bool has_internal = false;
			bool has_external = false;
			for (int a = 0; a < nen && (!has_internal || !has_external); a++)
			{
				int node = pel[a];
				if (status[node] == kInternal || status[node] == kBorder)
					has_internal = true;
				else
					has_external = true;
			}
		
			/* re-mark nodes */
			if (has_internal && has_external)
			{
				for (int a = 0; a < nen; a++)
				{
					int node = pel[a];
					if (status[node] == kInternal)
						status[node] = kBorder;
					else if (status[node] == kUnset)
						status[node] = kExternal;
				}
			}
		}
	}

	/* count-em up */
	iArrayT comm_list(fNumPartitions);
	comm_list = 0;
	int n_i, n_b, n_e;
	n_i = n_b = n_e = 0;
	StatusT* pstat = status.Pointer();
	for (int l = 0; l < nnd; l++)
	{
		if (*pstat == kInternal)
			n_i++;
		else if (*pstat == kBorder)
			n_b++;
		else if (*pstat == kExternal)
		{
			n_e++;
			comm_list[part_map[l]] = 1;
		}
		pstat++;
	}

	/* store (sorted) comm list */
	fCommID.Dimension(comm_list.Count(1));
	int dex = 0;
	for (int ll = 0; ll < fNumPartitions; ll++)
		if (comm_list[ll] == 1)
			fCommID[dex++] = ll;

	/* allocate node lists */
	fNodes_i.Dimension(n_i);
	fNodes_b.Dimension(n_b);
	fNodes_e.Dimension(n_e);

	/* sort-em out */
	n_i = n_b = n_e = 0;
	pstat = status.Pointer();
	for (int m = 0; m < nnd; m++)
	{
		if (*pstat == kInternal)
			fNodes_i[n_i++] = m;
		else if (*pstat == kBorder)
			fNodes_b[n_b++] = m;
		else if (*pstat == kExternal)
			fNodes_e[n_e++] = m;
		pstat++;
	}
	
	/* generate node map (just number sequentially through _i, _b, _e) */
	fNodeMap.Dimension(fNodes_i.Length() +
	                  fNodes_b.Length() +
	                  fNodes_e.Length()); // sets sequence for local node
	                                      // numbers - DO NOT CHANGE
	
	/* copy in (with no resequencing) */
	fNodeMap.CopyPart(0, fNodes_i, 0, fNodes_i.Length());
	fNodeMap.CopyPart(fNodes_i.Length(), fNodes_b, 0, fNodes_b.Length());
	fNodeMap.CopyPart(fNodes_i.Length() + fNodes_b.Length(), fNodes_e, 0, fNodes_e.Length());
}

/* map status of (in range) parts into status_map */
void PartitionT::MapStatus(StatusT status, const iArrayT& part,
	ArrayT<StatusT>& status_map, int offset)
{
	int  range = status_map.Length();
	int length = part.Length();
	const int* ppart = part.Pointer();
	for (int i = 0; i < length; i++)
	{
		int dex = *ppart++ - offset;
		if (dex >= 0 && dex < range) status_map[dex] = status;
	}
}

/* set send nodes/partition information */
void PartitionT::SetReceive(const iArrayT& part_map)
{
	/* (partition - min) -> index in fCommID */
	int min;
	iArrayT ID_map;
	MakeInverseMap(fCommID, ID_map, min);

	/* count nodes to each ID */
	iArrayT counts(ID_map.Length());
	counts = 0;
	int n_e = fNodes_e.Length();
	int* pexternal = fNodes_e.Pointer();
	for (int j = 0; j < n_e; j++)
	{
		int dex = part_map[*pexternal++] - min;
		counts[dex]++;
	}
	
	/* allocate communication map */
	fNodes_in.Dimension(fCommID.Length());
	for (int k = 0; k < fCommID.Length(); k++)
		fNodes_in[k].Dimension(counts[fCommID[k] - min]);

	/* store communication map */
	counts = 0;
	pexternal = fNodes_e.Pointer();
	for (int i = 0; i < n_e; i++)
	{
		int dex = part_map[*pexternal] - min;
		iArrayT& row = fNodes_in[ID_map[dex]];
		int& count = counts[dex];
		row[count++] = *pexternal++;
	}
}

/* set numbering maps */
void PartitionT::SetNodeMap(NumberScopeT scope, iArrayT& map, int& shift) const
{
	if (scope == kLocal)
	{
		/* non-const this */
		PartitionT* this_tmp = (PartitionT*) this;
	
		/* construct inverse node map */
		if (fInvNodeMap.Length() == 0)
			MakeInverseMap(fNodeMap,
			               this_tmp->fInvNodeMap,
			               this_tmp->fNodeMapShift);

		shift = fNodeMapShift;
		map.Alias(fInvNodeMap);
	}
	else
	{
		shift = 0;
		map.Alias(fNodeMap);
	}
}

void PartitionT::SetElementMap(NumberScopeT scope, const StringT& blockID, iArrayT& map,
	int& shift) const
{
	int dex = ElementBlockIndex(blockID, "SetElementMap");
	if (scope == kLocal)
	{
		/* construct inverse element map */
		if (fInvElementMap[dex].Length() == 0)
			MakeInverseMap(fElementMap[dex], 
				const_cast<iArrayT&>(fInvElementMap[dex]),
				const_cast<int&>(fElementMapShift[dex]));
				
		shift = fElementMapShift[dex];
		map.Alias(fInvElementMap[dex]);
	}
	else
	{
		shift = 0;
		map.Alias(fElementMap[dex]);
	}
}
