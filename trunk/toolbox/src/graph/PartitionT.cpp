/* $Id: PartitionT.cpp,v 1.1.1.1 2001-01-25 20:56:27 paklein Exp $ */
/* created: paklein (11/16/1999)                                          */
/* graph partition information (following NEMESIS data model)             */

#include "PartitionT.h"

#include "fstreamT.h"
#include "GraphT.h"
#include "AutoArrayT.h"
#include "iArray2DT.h"
#include "StringT.h"

/* parameters */
const int kHeadRoom = 20;             // percent
const char* sPartitionTVersion = "v1.0"; // version marker

/* constructor */
PartitionT::PartitionT(void):
	fNumPartitions(0),
	fID(-1),
	fScope(kUnSet)
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
	nodes.Allocate(fNodes_i.Length() + fNodes_b.Length());

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
void PartitionT::Set(int num_parts, int id, const iArrayT& part_map, const GraphT& graph)
{
	/* total number of partitions */
	fNumPartitions = num_parts;
	if (fNumPartitions < 1) throw eGeneralFail;

	/* set ID */
	fID = id;
	if (fID < 0 || fID >= fNumPartitions) throw eOutOfRange;
		
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
		throw eSizeMismatch;
	}
	
	/* only at global scope for now */
	if (fScope != kGlobal)
	{
		cout << "\n PartitionT::SetOutgoing: number scope must be global" << endl;
		throw eGeneralFail;
	}

	/* store */
	fNodes_out.Allocate(nodes_out.Length());
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

void PartitionT::InitElementBlocks(const iArrayT& blockID)
{
	/* copy */
	fElementBlockID = blockID;
	
	/* allocate memory */
	int n = fElementBlockID.Length();
	fElements_i.Allocate(n);
	fElements_b.Allocate(n);
	fElementMap.Allocate(n);
	fElementMapShift.Allocate(n);
	fInvElementMap.Allocate(n);
}

/* collect internal and border elements */
void PartitionT::SetElements(int blockID, const iArray2DT& connects)
{
	/* only at global scope for now */
	if (fScope != kGlobal)
	{
		cout << "\n PartitionT::SetElements: number scope must be global"
		     << endl;
		throw eGeneralFail;
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
	
	AutoArrayT<int> elements_i(20, true);
	AutoArrayT<int> elements_b(20, true);
	
	/* sort elements into internal/external */
	int nel = connects.MajorDim();
	int nen = connects.MinorDim();
	int* pelem = connects.Pointer();
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
	fElements_i[dex].Allocate(n_i);
	elements_i.CopyInto(fElements_i[dex]);

	int n_b = elements_b.Length();
	fElements_b[dex].Allocate(n_b);
	elements_b.CopyInto(fElements_b[dex]);
	
	/* set map */
	fElementMap[dex].Allocate(n_i + n_b);
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
			int* p1 = fNodes_in[i1].Pointer();
			int* p2 = (that.fNodes_out[i2]).Pointer();
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
			int* p1 = fNodes_out[i1].Pointer();
			int* p2 = (that.fNodes_in[i2]).Pointer();
			for (int i = 0; i < n1; i++)
				if (fNodeMap[*p1++] != (that.fNodeMap)[*p2++]) return 0;
		}
	}

	/* true on fall through */
	return 1;
}

/* I/O */
ostream& operator<<(ostream& out, const PartitionT& partition)
{
	out << sPartitionTVersion << '\n';
	out << partition.fNumPartitions << '\n'; // number of parts
	out << partition.fID            << '\n'; // partition number
	out << partition.fScope         << '\n'; // numbering scope
	
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
	out << (partition.fElementBlockID).wrap_tight(10) << '\n';

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

/* operator support */
ifstreamT& PartitionT::Read(ifstreamT& in)
{
	/* check version */
	StringT version;
	in >> version;
	if (version != sPartitionTVersion)
	{
		cout << "\n operator>>PartitionT&: file version " << version
		     << " is not current: " << sPartitionTVersion << endl;
		throw eGeneralFail;
	}

	int length;

	in >> fNumPartitions; // number of parts
	in >> fID;            // partition number
	in >> fScope;         // numbering scope
	
	// nodal information
	in >> length;
	fNodes_i.Allocate(length);
	in >> fNodes_i; // internal nodes	
	
	in >> length;
	fNodes_b.Allocate(length);
	in >> fNodes_b; // border nodes	

	in >> length;
	fNodes_e.Allocate(length);
	in >> fNodes_e; // external nodes
	
	// receive/send information
	in >> length;
	fCommID.Allocate(length);
	in >> fCommID; // ID's of communicating partitions

	fNodes_in.Allocate(length);
	for (int i = 0; i < length; i++)
	{
		int dim;
		in >> dim;
		fNodes_in[i].Allocate(dim);
		in >> fNodes_in[i];
	}

	fNodes_out.Allocate(length);
	for (int j = 0; j < length; j++)
	{
		int dim;
		in >> dim;
		fNodes_out[j].Allocate(dim);
		in >> fNodes_out[j];
	}
	
	// element information
	in >> length;
	fElementBlockID.Allocate(length);
	in >> fElementBlockID;
	InitElementBlocks(fElementBlockID);	

	for (int k = 0; k < fElements_i.Length(); k++)
	{
		int dim;
		in >> dim;
		fElements_i[k].Allocate(dim);
		in >> fElements_i[k]; // internal elements
	}

	for (int l = 0; l < fElements_b.Length(); l++)
	{
		int dim;
		in >> dim;
		fElements_b[l].Allocate(dim);
		in >> fElements_b[l]; // internal elements
	}
	
	// global node map
	in >> length;
	fNodeMap.Allocate(length);
	in >> fNodeMap; // global[local]
	if (length != (fNodes_i.Length() +
	               fNodes_b.Length() +
fNodes_e.Length())) throw eBadInputValue;

	// block global element numbering map
	for (int m = 0; m < fElementMap.Length(); m++)
	{
		int dim;
		in >> dim;
		fElementMap[m].Allocate(dim);
		in >> fElementMap[m];
	}
	
	return in;
}

/* resolve element block ID to index */
int PartitionT::ElementBlockIndex(int blockID, const char* caller) const
{
	int dex = 0;
	if (!fElementBlockID.HasValue(blockID, dex))
	{
		const char* this_routine = "ElementBlockIndex";
		const char* str = (caller != NULL) ? caller : this_routine;
		cout << "\n PartitionT::" << str << ": block ID not found: "
		     << blockID << endl;
		throw eGeneralFail;
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
			throw eOutOfRange;
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

	AutoArrayT<int> tmp(20, true);
	for (int i = 0; i < global_nodes.Length(); i++)
	{
		int snode = global_nodes[i] - shift;
		if (snode > -1 && snode < inv_map.Length())
			if (inv_map[snode] > -1)
				tmp.Append(i);

	}

	/* copy to return value */
	partition_indices.Allocate(tmp.Length());
	tmp.CopyInto(partition_indices);
}

/* returns indeces of (block) global elements that lie within
* the partition */
void PartitionT::ReturnPartitionElements(int blockID,
	const iArrayT& global_elements, iArrayT& partition_indices) const
{
	/* make inverse map */
	iArrayT inv_map;
	int shift;
	MakeInverseMap(ElementMap(blockID), inv_map, shift);

	AutoArrayT<int> tmp(20, true);
	for (int i = 0; i < global_elements.Length(); i++)
	{
		int selement = global_elements[i] - shift;
		if (selement > -1 && selement < inv_map.Length())
			if (inv_map[selement] > -1)
				tmp.Append(i);
	}

	/* copy to return value */
	partition_indices.Allocate(tmp.Length());
	tmp.CopyInto(partition_indices);
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

void PartitionT::SetElementScope(NumberScopeT scope, int blockID, ArrayT<int>& elements) const
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
			cout << "\n operator>>PartitionT::NumberScopeT: unknown value: " << i_scope << endl;
			throw eBadInputValue;
	}

	return in;
}

/************************************************************************
* Private
************************************************************************/

/* make inverse map (filled with -1) */
void PartitionT::MakeInverseMap(const iArrayT& map, iArrayT& inv_map,
	int& shift) const
{
	/* range */
	int max;
	map.MinMax(shift, max);
	int range = max - shift + 1;
	
	/* dimension */
	inv_map.Allocate(range);
	inv_map = -1;

	/* make map */
	int dim = map.Length();
	for (int i = 0; i < dim; i++)
		inv_map[map[i] - shift] = i;	
}

/* set node info */
void PartitionT::ClassifyNodes(const iArrayT& part_map,
	const GraphT& graph)
{
	/* work space */
	AutoArrayT<int> nodes_i(kHeadRoom, true);
	AutoArrayT<int> nodes_b(kHeadRoom, true);
	AutoArrayT<int> nodes_e(kHeadRoom, true);
	AutoArrayT<int> commID(kHeadRoom, true);

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
			int *pedge = graph.Edges(i);
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
	fNodes_i.Allocate(nodes_i.Length());
	nodes_i.CopyInto(fNodes_i);

	fNodes_b.Allocate(nodes_b.Length());
	nodes_b.CopyInto(fNodes_b);	

	fNodes_e.Allocate(nodes_e.Length());
	nodes_e.CopyInto(fNodes_e);	

	fCommID.Allocate(commID.Length());
	commID.CopyInto(fCommID);
	
	/* generate node map (just number sequentially through _i, _b, _e) */
	fNodeMap.Allocate(fNodes_i.Length() +
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
	int* ppart = part.Pointer();
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
	fNodes_in.Allocate(fCommID.Length());
	for (int k = 0; k < fCommID.Length(); k++)
		fNodes_in[k].Allocate(counts[fCommID[k] - min]);

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

void PartitionT::SetElementMap(NumberScopeT scope, int blockID, iArrayT& map,
	int& shift) const
{
	int dex = ElementBlockIndex(blockID, "SetElementMap");
	if (scope == kLocal)
	{
		/* construct inverse element map */
		if (fInvElementMap[dex].Length() == 0)
			MakeInverseMap(fElementMap[dex], fInvElementMap[dex],
				fElementMapShift[dex]);
				
		shift = fElementMapShift[dex];
		map.Alias(fInvElementMap[dex]);
	}
	else
	{
		shift = 0;
		map.Alias(fElementMap[dex]);
	}
}
