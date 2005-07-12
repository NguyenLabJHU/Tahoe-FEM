/* $Id: CommManagerT.cpp,v 1.17.2.11 2005-06-27 16:53:54 paklein Exp $ */
#include "CommManagerT.h"
#include "CommunicatorT.h"
#include "ModelManagerT.h"
#include "NodeManagerT.h"
#include "PartitionT.h"
#include "InverseMapT.h"
#include "FieldT.h"
#include "SpatialGridT.h"
#include <float.h>
#include "ifstreamT.h"
#include "ofstreamT.h"

/* message types */
#include "AllGatherT.h"
#include "PointToPointT.h"
#include "CartesianShiftT.h"

using namespace Tahoe;

/* debugging */
#define CommManagerT_DEBUG_SPATIAL 1
#include "TextOutputT.h"

CommManagerT::CommManagerT(CommunicatorT& comm, ModelManagerT& model_manager):
	fComm(comm),
	fModelManager(model_manager),
	fSize(fComm.Size()),
	fRank(fComm.Rank()),
	fSkin(0.0),
	fPartition(NULL),
	fNodeManager(NULL),
	fNumRealNodes(0),
	fd_send_buffer_man(fd_send_buffer),
	fd_recv_buffer_man(fd_recv_buffer),
	fi_send_buffer_man(fi_send_buffer),
	fi_recv_buffer_man(fi_recv_buffer),
	fCartesianShift(NULL),
	fNodalAttributes(0, false),
	fAttributeMessageID(kNULLMessageID)
{
	/* dimension the attributes array group */
	fNodalAttributes.Dimension(fModelManager.NumNodes(), false);
}

CommManagerT::~CommManagerT(void)
{
	/* free any remaining communications */
	for (int i = 0; i < fCommunications.Length(); i++)
		delete fCommunications[i];

	/* free any remaining ghost communications */
	for (int i = 0; i < fGhostCommunications.Length(); i++)
		delete fGhostCommunications[i];

	delete fCartesianShift;
	delete fPartition;
}

/* set or clear partition information */
void CommManagerT::ReadPartition(const StringT& part_file)
{
	const char caller[] = "CommManagerT::ReadPartition";

	int token = 1;

	/* check file */
	fPartFile = part_file;
	ifstreamT part_in('#', fPartFile);
	if (!part_in.is_open()) {
		cout << "\n " << caller << ": could not open file: " << fPartFile << endl;
		token = 0;	
	}
	if (fComm.Sum(token) != Size()) ExceptionT::GeneralFail(caller, "error reading parition files");
	
	/* read partition information */
	fPartition = new PartitionT;
	part_in >> (*fPartition);
	if (fPartition->ID() != Rank()) {
		cout << "\n " << caller << "partition ID " << fPartition->ID() << " does not match process rank " << Rank() << endl;
		token = 0;
	}
	if (fComm.Sum(token) != Size()) ExceptionT::GeneralFail(caller, "partition file error");

	/* node map */
	fNodeMap.Alias(fPartition->NodeMap());
		
	/* nodes owned by this processor */
	fProcessor.Dimension(fNodeMap.Length());
	fPartition->ReturnProcessorMap(fProcessor);
	CollectPartitionNodes(fProcessor, fComm.Rank(), fPartitionNodes, fExternalNodes);

	/* clear the inverse map */
	fPartitionNodes_inv.ClearMap();
}

/* set or clear node manager information */
void CommManagerT::SetNodeManager(NodeManagerT* node_manager)
{
	fNodeManager = node_manager;

	/* dimension space for periodic BC */
	int nsd = fModelManager.NumDimensions();
	fIsPeriodic.Dimension(nsd);
	fIsPeriodic = false;
	fPeriodicBoundaries.Dimension(nsd, 2);
	fPeriodicBoundaries = 0.0;
	fPeriodicLength.Dimension(nsd);
	fPeriodicLength = 0.0;
}

/* set boundaries */
void CommManagerT::SetPeriodicBoundaries(int i, double x_i_min, double x_i_max)
{
	if (x_i_min > x_i_max) ExceptionT::GeneralFail();
	fIsPeriodic[i] = true;
	fPeriodicBoundaries(i,0) = x_i_min;
	fPeriodicBoundaries(i,1) = x_i_max;
	fPeriodicLength[i] = x_i_max - x_i_min;
}

/* unset boundaries */
void CommManagerT::ClearPeriodicBoundaries(int i) { fIsPeriodic[i] = false; }

/* perform actions needed the first time CommManagerT::Configure is called. */
void CommManagerT::Initialize(void)
{
	const char caller[] = "CommManagerT::Initialize";

	/* check */
	if (!fPartition && Size() > 1)
		ExceptionT::GeneralFail(caller, "partition not set");

	/* serial */
	if (Size() == 1) {
		fNumRealNodes = fModelManager.NumNodes();
		fProcessor.Dimension(fNumRealNodes);
		fProcessor = fComm.Rank();
		return;
	}
	
	/* graph decomposition - nothing to do? */
	if (fPartition->DecompType() == PartitionT::kGraph)
		return;
	else if (fPartition->DecompType() == PartitionT::kIndex)
	{
		/* change the numbering scope of partition data to global */
		fPartition->SetScope(PartitionT::kGlobal);
	
		/* number of nodes owned by this partition */
		int npn = fPartition->NumPartitionNodes();
		int nsd = fModelManager.NumDimensions();

		/* total number of nodes */
		AllGatherT all_gather(fComm, 0);
		all_gather.Initialize(npn*nsd);
		fNumRealNodes = all_gather.Total()/nsd;

		/* collect global reference coordinate list */
		const dArray2DT& reference_coords = fModelManager.Coordinates();
		if (reference_coords.MajorDim() != fNumRealNodes)
		{
			/* collect */
			dArray2DT coords_all(fNumRealNodes, nsd);		
			all_gather.AllGather(reference_coords, coords_all);

			/* reset model manager */
			fModelManager.UpdateNodes(coords_all, true);
		}

		/* dimension the attributes array group */
		fNodalAttributes.Dimension(fModelManager.NumNodes(), true);

		/* translate element sets to global numbering */
		const ArrayT<StringT>& el_id = fModelManager.ElementGroupIDs();
		for (int i = 0; i < el_id.Length(); i++)
		{
			/* copy existing connectivities */
			iArray2DT elems = fModelManager.ElementGroup(el_id[i]);

			/* map scope of nodes in connectivities */
			fPartition->SetNodeScope(PartitionT::kGlobal, elems);
			
			/* re-set the connectivities */
			fModelManager.UpdateElementGroup(el_id[i], elems, true);
		}

		/* translate node sets to global numbering and all gather */
		const ArrayT<StringT>& nd_id = fModelManager.NodeSetIDs();
		for (int i = 0; i < nd_id.Length(); i++)
		{
			/* copy node set */
			iArrayT partial_nodes = fModelManager.NodeSet(nd_id[i]);

			/* map scope of nodes in the node set */
			fPartition->SetNodeScope(PartitionT::kGlobal, partial_nodes);
			
			/* collect all */
			AllGatherT gather_node_set(fComm, 0);
			gather_node_set.Initialize(partial_nodes.Length());
			iArrayT all_nodes(gather_node_set.Total());
			gather_node_set.AllGather(partial_nodes, all_nodes);
			
			/* re-set the node set */
			fModelManager.UpdateNodeSet(nd_id[i], all_nodes, true);		
		}

		/* translate side sets to global numbering */	
		const ArrayT<StringT>& ss_id = fModelManager.SideSetIDs();
		//TEMP
		if (ss_id.Length() > 0) cout << "\n CommManagerT::FirstConfigure: skipping size sets" << endl;

		/* set node-to-processor map - could also get from n/(part size) */
		fProcessor.Dimension(fModelManager.NumNodes());
		const iArrayT& counts = all_gather.Counts();
		int proc = 0;
		int* p = fProcessor.Pointer();
		for (int i = 0; i < counts.Length(); i++)
		{
			int n_i = counts[i]/nsd;
			for (int j = 0; j < n_i; j++)
				*p++ = proc;
			proc++;
		}
		
		/* clear the node numbering map - each processor knows about all nodes */
		fNodeMap.Dimension(0);
		
		/* collect partition nodes (general algorithm) */
		CollectPartitionNodes(fProcessor, fComm.Rank(), fPartitionNodes, fExternalNodes);
		
		/* clear the inverse map */
		fPartitionNodes_inv.ClearMap();
		
		/* (potentially) all are adjacent to external nodes */
		//fBorderNodes.Alias(fPartitionNodes);
		
		/* all the rest are external */
		fExternalNodes.Dimension(fNumRealNodes - fPartitionNodes.Length());
		int lower, upper;
		lower = upper = fNumRealNodes;
		if (fPartitionNodes.Length() > 0) {
			lower = fPartitionNodes.First();
			upper = fPartitionNodes.Last();
		}
		int dex = 0;
		for (int i = 0; i < lower; i++)
			fExternalNodes[dex++] = i;
		for (int i = upper+1; i < fNumRealNodes; i++)
			fExternalNodes[dex++] = i;
	}
	else if (fPartition->DecompType() == PartitionT::kSpatial)
	{
		int nsd = fModelManager.NumDimensions();

		/* determine the bounds and neighbors of this processor - (reference coordinates) */
		fBounds.Dimension(nsd, 2);         /* [nsd] x {low, high} */
		fAdjacentCommID.Dimension(nsd, 2); /* [nsd] x {low, high} */
		GetProcessorBounds(fModelManager.Coordinates(), fBounds, fAdjacentCommID);

		/* swap pattern {+x, +y, +z,..., -x, -y, -z,...} */
		fSwap.Dimension(2*nsd, 2);
		for (int i = 0; i < nsd; i++)
		{
			/* forward */
			fSwap(i,0) = i; /* direction {x, y, z,...} */
			fSwap(i,1) = 1; /* sign +/- */
	
			/* back */
			fSwap(i+nsd,0) = i; /* direction {x, y, z,...} */
			fSwap(i+nsd,1) = 0; /* sign +/- */
		}

#ifdef CommManagerT_DEBUG_SPATIAL
cout << "\n processor bounds = \n" << fBounds << '\n';
cout << "\n adjacent = \n" << fAdjacentCommID << '\n';
cout << "\n swaps = \n" << fSwap << '\n';
#endif /* CommManagerT_DEBUG_SPATIAL */
		
		/* dimension the swapping node lists */
		fSendNodes.Dimension(2*nsd);
		fRecvNodes.Dimension(2*nsd);

		/* initialize shift communications */
		fCartesianShift = new CartesianShiftT(fComm, 0, fAdjacentCommID, fSwap, fSendNodes, fRecvNodes);
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized decomposition type %d", fPartition->DecompType());
}

/* Update system configuration, including enforcement of periodic boundary conditions */
void CommManagerT::UpdateConfiguration(void)
{
	const char caller[] = "CommManagerT::UpdateConfiguration";

#ifdef CommManagerT_DEBUG_SPATIAL
/* debugging */
cout << "\n " << caller << ": IN\n";

/* partition nodes */
cout << "\npartition nodes: " << fPartitionNodes.Length() << '\n';
if (1 && fNodeMap.Length() > 0) {
	iArrayT tmp(fPartitionNodes.Length());
	tmp.Collect(fPartitionNodes, fNodeMap);
	tmp++;
	cout << tmp.wrap(10) << '\n';
}
cout << endl;
#endif /* CommManagerT_DEBUG_SPATIAL */

	PartitionT::DecompTypeT decomp_type = (fPartition) ? fPartition->DecompType() : PartitionT::kUndefined;
	if (decomp_type == PartitionT::kIndex || decomp_type == PartitionT::kUndefined) /* serial or index decomposition */
	{
		/* update image nodes */
		EnforcePeriodicBoundaries();

		/* update attributes */
		int nnd = fNodeManager->NumNodes();
		fNodalAttributes.Dimension(nnd, true);
		if (fAttributeMessageID == kNULLMessageID) fAttributeMessageID = Init_AllGather(MessageT::Integer, 1);
		iArray2DT attr_wrapper;
		for (int i = 0; i < fNodalAttributes.NumRegistered(); i++) {
			attr_wrapper.Alias(nnd, 1, fNodalAttributes.Array(i).Pointer());
			AllGather(fAttributeMessageID, attr_wrapper);
		}
	}
	else if (decomp_type == PartitionT::kSpatial) /* spatial decomposition */
	{
//TEMP
for (int i = 0; i < fIsPeriodic.Length(); i++)
	if (fIsPeriodic[i])
		ExceptionT::GeneralFail(caller, "periodic boundary conditions not supported");
//TEMP

		/* work space */
		iArray2DT i_values;
		nVariArray2DT<int> i_values_man;
		dArray2DT new_init_coords;
		nVariArray2DT<double> new_init_coords_man;
		dArray2DT new_curr_coords;
		nVariArray2DT<double> new_curr_coords_man;

		/* determine the bounds and neighbors of this processor - (reference coordinates) */		
		int npn = fPartitionNodes.Length(); /* current number of nodes owned by this partition */
		dArray2DT curr_coords(npn, fModelManager.NumDimensions(), fNodeManager->CurrentCoordinates().Pointer()); /* processor nodes are 0...(npn-1) */
		GetProcessorBounds(curr_coords, fBounds, fAdjacentCommID);

		/* wrap in try block */
		int token = 1;
		try {
			/* init data needed for reconfiguring across processors */
			InitConfigure(i_values, i_values_man, new_init_coords, new_init_coords_man, new_curr_coords, new_curr_coords_man);
	
			/* distribute nodes over the spatial grid */
			Distribute(i_values, i_values_man, new_init_coords, new_init_coords_man, new_curr_coords, new_curr_coords_man);

#ifdef CommManagerT_DEBUG_SPATIAL
cout << "\ncaller: number of partition nodes: " << new_init_coords.MajorDim() << '\n';
#endif
			
			/* exchange border nodes */
			SetExchange(i_values, i_values_man, new_init_coords, new_init_coords_man, new_curr_coords, new_curr_coords_man);
	
			/* finalize configuration */
			CloseConfigure(i_values, new_init_coords);

#ifdef CommManagerT_DEBUG_SPATIAL
/* debugging */
cout << "\ncaller: OUT\n";

/* coordinates */
cout << "\nreference coordinates:" << fNodeManager->NumNodes() << '\n';
//TextOutputT::WriteNodeValues(cout, fNodeMap, fNodeManager->InitialCoordinates());
//cout << "\ncurrent coordinates:\n";
//TextOutputT::WriteNodeValues(cout, fNodeMap, fNodeManager->CurrentCoordinates());

/* node sets */
const ArrayT<StringT>& ns_ID = fModelManager.NodeSetIDs();
for (int i = 0; i < ns_ID.Length(); i++) {
	const iArrayT& ns = fModelManager.NodeSet(ns_ID[i]);
	cout << "\n node set: " << ns_ID[i] << '\n';
	for (int j = 0; j < ns.Length(); j++)
		cout << setw(kIntWidth) << ns[j]+1
		     << setw(kIntWidth) << fNodeMap[ns[j]]+1 << '\n';
	cout << '\n';
}

/* partition nodes */
cout << "\npartition nodes: " << fPartitionNodes.Length() << '\n';
if (1) {
	iArrayT tmp(fPartitionNodes.Length());
	tmp.Collect(fPartitionNodes, fNodeMap);
	tmp++;
	cout << tmp.wrap(10) << '\n';
}
cout << endl;
#endif /* CommManagerT_DEBUG_SPATIAL */
		}		
		catch (ExceptionT::CodeT error) { token = 0; }
		
		/* synch and check */
		if (fComm.Sum(token) != Size()) ExceptionT::GeneralFail(caller);
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized decomposition type %d", decomp_type);
}

/* return true if the list of partition nodes may be changing */
bool CommManagerT::PartitionNodesChanging(void) const
{
	if (fPartition && fPartition->DecompType() == PartitionT::kSpatial)
		return true;
	else
		return false;
}

/* inverse of fPartitionNodes list */
const InverseMapT* CommManagerT::PartitionNodes_inv(void) const
{
	const ArrayT<int>* nodes = PartitionNodes();
	if (!nodes)
		return NULL;
	else
	{
		/* cast away const-ness */
		CommManagerT* non_const_this = (CommManagerT*) this;
	
		/* set inverse map */
		if (nodes->Length() > 0 && fPartitionNodes_inv.Length() == 0) {
			nArrayT<int> nd;
			nd.Alias(*nodes);
			non_const_this->fPartitionNodes_inv.SetMap(nd);
		}
		
		return &fPartitionNodes_inv;
	}
}

/* set up */
int CommManagerT::Init_AllGather(MessageT::TypeT t, int num_vals)
{
	const char caller[] = "CommManagerT::Init_AllGather";

	/* not multiprocessor - return dummy id */
	if (!fPartition) return -1;

	/* store the number of values */
	fNumValues.Append(num_vals);

	int message_tag = fCommunications.Length();
	PartitionT::DecompTypeT decomp_type = fPartition->DecompType();
	if (decomp_type == PartitionT::kGraph) /* non-blocking send-receive */
	{
		/* new point-to-point */
		PointToPointT* p2p = new PointToPointT(fComm, message_tag, *fPartition);

		/* allocate buffers */
		p2p->Initialize(t, num_vals);
	
		/* store */
		fCommunications.Append(p2p);
		
		/* no ghost nodes allowed */
		fGhostCommunications.Append(NULL);
	}
	else if (decomp_type == PartitionT::kIndex) /* all gather */
	{
		/* new all gather */
		AllGatherT* all_gather = new AllGatherT(fComm, message_tag);

		/* set message size */
		all_gather->Initialize(fPartition->NumPartitionNodes()*num_vals);
		
		/* store */
		fCommunications.Append(all_gather);

		/* new all gather for gh */
		AllGatherT* ghost_all_gather = new AllGatherT(fComm, message_tag);

		/* set message size */
		ghost_all_gather->Initialize(fPBCNodes.Length()*num_vals);
		
		/* store */
		fGhostCommunications.Append(ghost_all_gather);
	}
	else if (decomp_type == PartitionT::kSpatial) /* shifts */
	{
		/* dummies */
		fCommunications.Append(NULL);
		fGhostCommunications.Append(NULL);
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized decomp type: %d", decomp_type);

	/* ID is just the array position */
	return message_tag;
}

/* clear the persistent communication */
void CommManagerT::Clear_AllGather(int id)
{
	/* ignore dummy id (not multiprocessor) */
	if (id == -1) return;

	/* check */
	if (id < 0 || id >= fCommunications.Length())
		ExceptionT::OutOfRange("CommManagerT::Clear_AllGather");

	/* clear the number of values */
	fNumValues[id] = 0;

	/* free memory and purge from array */
	delete fCommunications[id];
	fCommunications[id] = NULL;

	/* ghost nodes - free memory and purge from array */
	delete fGhostCommunications[id];
	fGhostCommunications[id] = NULL;
}

/* perform the all gather */
void CommManagerT::AllGather(int id, nArray2DT<double>& values)
{
	const char caller[] = "CommManagerT::AllGather";

	/* not multiprocessor - ghost nodes only */
	if (!fPartition) {

		/* copy values to ghost nodes */
		int ngh = fPBCNodes.Length();
		for (int i = 0; i < ngh; i++)
			values.SetRow(fPBCNodes_ghost[i], values(fPBCNodes[i]));

		/* nothing else */
		return;
	}

	PartitionT::DecompTypeT decomp_type = fPartition->DecompType();
	if (decomp_type == PartitionT::kGraph) /* non-blocking send-receive */
	{
		/* retrieve pointer */
		PointToPointT* p2p = TB_DYNAMIC_CAST(PointToPointT*, fCommunications[id]);
		if (!p2p) ExceptionT::GeneralFail(caller);

		/* do it */
		p2p->AllGather(values);
	}
	else if (decomp_type == PartitionT::kIndex) /* all gather */
	{
		/* retrieve pointer */
		AllGatherT* all_gather = TB_DYNAMIC_CAST(AllGatherT*, fCommunications[id]);
		if (!all_gather) ExceptionT::GeneralFail(caller);
		
		/* exchange */
		dArray2DT exchange;
		exchange.Alias(fNumRealNodes, values.MinorDim(), values(0));
		all_gather->AllGather(exchange);

		/* ghost nodes */
		if (fModelManager.NumNodes() > fNumRealNodes) {
		
			/* copy values for local ghost nodes */
			int ngh = fPBCNodes.Length();
			for (int i = 0; i < ngh; i++)
				values.SetRow(fPBCNodes_ghost[i], values(fPBCNodes[i]));

			/* retrieve communications pointer */
			AllGatherT* all_gather = TB_DYNAMIC_CAST(AllGatherT*, fGhostCommunications[id]);
			if (!all_gather) ExceptionT::GeneralFail(caller);

			/* exchange */
			exchange.Alias(fModelManager.NumNodes() - fNumRealNodes, values.MinorDim(), values(fNumRealNodes));
			all_gather->AllGather(exchange);
		}		
	}
	else if (decomp_type == PartitionT::kSpatial) /* shifts */
	{
		/* check */
		if (!fCartesianShift) ExceptionT::GeneralFail(caller, "shift not set");

		/* exchange */
		fCartesianShift->AllGather(values);
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized decomp type: %d", decomp_type);
}

void CommManagerT::AllGather(int id, nArray2DT<int>& values)
{
	const char caller[] = "CommManagerT::AllGather";

	/* not multiprocessor - ghost nodes only */
	if (!fPartition) {

		/* copy values to ghost nodes */
		int ngh = fPBCNodes.Length();
		for (int i = 0; i < ngh; i++)
			values.SetRow(fPBCNodes_ghost[i], values(fPBCNodes[i]));

		/* nothing else to do */
		return;
	}

	PartitionT::DecompTypeT decomp_type = fPartition->DecompType();
	if (decomp_type == PartitionT::kGraph) /* non-blocking send-receive */
	{
		/* retrieve pointer */
		PointToPointT* p2p = TB_DYNAMIC_CAST(PointToPointT*, fCommunications[id]);
		if (!p2p) ExceptionT::GeneralFail(caller);

		/* do it */
		p2p->AllGather(values);
	}
	else if (decomp_type == PartitionT::kIndex) /* all gather */
	{
		/* retrieve communications pointer */
		AllGatherT* all_gather = TB_DYNAMIC_CAST(AllGatherT*, fCommunications[id]);
		if (!all_gather) ExceptionT::GeneralFail(caller);
		
		/* exchange */
		iArray2DT exchange;
		exchange.Alias(fNumRealNodes, values.MinorDim(), values(0));
		all_gather->AllGather(exchange);
		
		/* ghost nodes */
		if (fModelManager.NumNodes() > fNumRealNodes) {
		
			/* copy values to ghost nodes */
			int ngh = fPBCNodes.Length();
			for (int i = 0; i < ngh; i++)
				values.SetRow(fPBCNodes_ghost[i], values(fPBCNodes[i]));

			/* retrieve communications pointer */
			AllGatherT* all_gather = TB_DYNAMIC_CAST(AllGatherT*, fGhostCommunications[id]);
			if (!all_gather) ExceptionT::GeneralFail(caller);

			/* exchange */
			exchange.Alias(fModelManager.NumNodes() - fNumRealNodes, values.MinorDim(), values(fNumRealNodes));
			all_gather->AllGather(exchange);
		}		
	}
	else if (decomp_type == PartitionT::kSpatial) /* shifts */
	{
		/* check */
		if (!fCartesianShift) ExceptionT::GeneralFail(caller, "shift not set");

		/* exchange */
		fCartesianShift->AllGather(values);
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized decomp type: %d", decomp_type);
}

/* register an array of nodal attributes */
int CommManagerT::RegisterNodalAttribute(nArrayT<int>& array) {
	return fNodalAttributes.Register(array);
}

/* send output data for writing */
void CommManagerT::WriteOutput(int print_step)
{
	const char caller[] = "CommManagerT::WriteOutput";

	/* write changing part file */
	if (PartitionNodesChanging() && fPartition) {

		//TEMP
		if (DecompType() != PartitionT::kSpatial)
			ExceptionT::GeneralFail(caller, "changing partition not supported with decomp type %d", DecompType());

		/* open stream */
		StringT file = fPartFile;
		file.Append(".ps", print_step, 4);
		ofstreamT out(file);
		if (!out.is_open())
			ExceptionT::GeneralFail(caller, "could not open file \"%s\"", file.Pointer());

		/* write partition information */
		out << *fPartition << '\n';
	}
}

/*************************************************************************
 * Private
 *************************************************************************/

/* collect partition nodes */
void CommManagerT::CollectPartitionNodes(const ArrayT<int>& n2p_map, int part, 
	AutoArrayT<int>& part_nodes, AutoArrayT<int>& external_nodes) const
{
	int count = 0;
	const int* p = n2p_map.Pointer();
	int len = n2p_map.Length();
	for (int i = 0; i < len; i++)
		if (*p++ == part)
			count++;
			
	part_nodes.Dimension(count);
	external_nodes.Dimension(len - count);
	int* pn = part_nodes.Pointer();
	int* en = external_nodes.Pointer();
	p = n2p_map.Pointer();
	count = 0;
	for (int i = 0; i < len; i++)
		if (*p++ == part)
			*pn++ = i;
		else
			*en++ = i;
}

/* enforce the periodic boundary conditions */
void CommManagerT::EnforcePeriodicBoundaries(void)
{
	const char caller[] = "CommManagerT::EnforcePeriodicBoundaries";
	if (!fNodeManager) ExceptionT::GeneralFail(caller, "node manager node set");

	/* only implemented for index decomposition (or serial) */
	if (fPartition && fPartition->DecompType() != PartitionT::kIndex)
		ExceptionT::GeneralFail(caller, "wrong decomposition type %d", fPartition->DecompType());
	
	/* reference coordinates */
	const dArray2DT& reference_coords = fModelManager.Coordinates();

	/* the coordinate update field */
	dArray2DT* field = fNodeManager->CoordinateUpdate();
	if (!field) ExceptionT::GeneralFail(caller, "no coordinate update array");
	dArray2DT& displacement = *field;

	/* nodes owned by this partition */
	const ArrayT<int>* partition_nodes = PartitionNodes();
	int nnd = (partition_nodes) ? partition_nodes->Length() : reference_coords.MajorDim();

	/* current number of "ghost" nodes */
	int ngn = fComm.Sum(fPBCNodes.Length());

	/* width of the communication layer */
	double skin = fSkin;

	/* reset ghost nodes */
	fPBCNodes.Dimension(0);
	fPBCNodes_ghost.Dimension(0);
	fPBCNodes_face.Dimension(0);

	/* loop over directions */
	int x_lower = 1; /* ...000001 */
	int x_upper = 2; /* ...000010 */
	bool has_periodic = false;
	for (int i = 0; i < fIsPeriodic.Length(); i++)
	{
		if (fIsPeriodic[i])
		{
			has_periodic = true;

			/* coordinate limits */
			double x_min = fPeriodicBoundaries(i,0);
			double x_max = fPeriodicBoundaries(i,1);
			double x_len = x_max - x_min;
			
			/* existing number of ghosts nodes */
			int ngh = fPBCNodes.Length();

			/* nodes owned by this partition */
			for (int j = 0; j < nnd; j++)
			{
				/* current coordinates */
				int nd = (partition_nodes) ? (*partition_nodes)[j] : j;
				const double& X = reference_coords(nd,i);
				double& d = displacement(nd,i);
				double  x = X + d;

				/* shift displacements - back in the box */
				if (x > x_max) {
					double dist = x - x_max;
					int cell = int(dist/x_len); /* which periodic cell */
					d -= (cell+1)*x_len;
				}
				else if (x < x_min) {
					double dist = x_min - x;
					int cell = int(dist/x_len);
					d += (cell+1)*x_len;
				}

				/* collect nodes close to periodic boundaries */
				x = X + d;
				if (x - x_min < skin) {
					fPBCNodes.Append(nd);
					fPBCNodes_face.Append(x_lower);
				}
				if (x_max - x < skin) {
					fPBCNodes.Append(nd);
					fPBCNodes_face.Append(x_upper);
				}
			}
			
			/* create images of images */
			for (int j = 0; j < ngh; j++)
			{
				/* current coordinates */
				int nd = fPBCNodes[j];
				const double& X = reference_coords(nd,i);
				double& d = displacement(nd,i);
				double  x = X + d;
			
				/* close to the boundaries */
				if (x - x_min < skin) {
					fPBCNodes.Append(nd);
					fPBCNodes_face.Append(fPBCNodes_face[j] | x_lower); /* bitwise OR */
				}
				if (x_max - x < skin) {
					fPBCNodes.Append(nd);
					fPBCNodes_face.Append(fPBCNodes_face[j] | x_upper); /* bitwise OR */
				}			
			}
		}
		
		/* next direction */
		x_lower <<= 2;
		x_upper <<= 2;
	}

	/* configure ghost nodes */
	if (has_periodic)
	{
		/* communicate number of ghost nodes per partition */
		iArrayT ghost_count(fComm.Size());
		fComm.AllGather(fPBCNodes.Length(), ghost_count);
		ngn = ghost_count.Sum();

		/* local numbering of ghost nodes */
		int ghost_num_start = fNumRealNodes;
		for (int i = 0; i < fComm.Rank(); i++)
			ghost_num_start += ghost_count[i];

		fPBCNodes_ghost.Dimension(fPBCNodes.Length());
		int ghost_num = ghost_num_start;
		for (int i = 0; i < fPBCNodes_ghost.Length(); i++)
			fPBCNodes_ghost[i] = ghost_num++;

		/* resize all coordinate and field arrays */
		fNodeManager->ResizeNodes(fNumRealNodes + ngn);

		/* generate reference coordinates - write access to the coordinates from the
		 * model manager, not from the node manager. The node manager does not change
		 * the initial coordinates during CopyNodeToNode */
		const dArray2DT& coords = fModelManager.Coordinates();
		dArray2DT ghost_coords;
		if (fPBCNodes_ghost.Length() > 0) {
			ghost_coords.Alias(fPBCNodes_ghost.Length(), coords.MinorDim(), coords(ghost_num_start));

			/* copy coordinates */
			for (int i = 0; i < fPBCNodes_ghost.Length(); i++)
				ghost_coords.SetRow(i, coords(fPBCNodes[i]));
			
			/* loop over directions and shift coords */
			int x_lower = 1; /* ...000001 */
			int x_upper = 2; /* ...000010 */
			for (int j = 0; j < fIsPeriodic.Length(); j++)
			{
				if (fIsPeriodic[j])
				{
					for (int i = 0; i < fPBCNodes_ghost.Length(); i++)
					{
						/* periodic shift */
						int face = fPBCNodes_face[i];
						
						/* at lower bound */
						if (face & x_lower) /* bitwise AND */
							ghost_coords(i,j) += fPeriodicLength[j];
						else if (face & x_upper) /* bitwise AND */
							ghost_coords(i,j) -= fPeriodicLength[j];
					}
				}

				/* next direction */
				x_lower <<= 2;
				x_upper <<= 2;
			}
		}
		else ghost_coords.Set(0, coords.MinorDim(), NULL);

		/* exchange */
		if (ngn > 0) {
			AllGatherT all_gather(fComm, 0);
			all_gather.Initialize(ghost_coords);
			dArray2DT ghost_coords_all(ngn, coords.MinorDim(), coords(fNumRealNodes));
			all_gather.AllGather(ghost_coords_all);
		}

		/* persistent communications */
		for (int i = 0; i < fGhostCommunications.Length(); i++)
		{
			/* retrieve communications pointer */
			AllGatherT* all_gather = TB_DYNAMIC_CAST(AllGatherT*, fGhostCommunications[i]);
			if (!all_gather) ExceptionT::GeneralFail(caller, "could not recover ghost communication");

			/* (re-)set message size */
			all_gather->Initialize(fPBCNodes.Length()*fNumValues[i]);
		}
		
		/* reset the node-to-processor map */
		fProcessor.Resize(coords.MajorDim());
		int* n2p = fProcessor.Pointer(fNumRealNodes);
		for (int i = 0; i < ghost_count.Length(); i++)
			for (int j = 0; j < ghost_count[i]; j++)
				*n2p++ = -(i+1);

		/* collect partition nodes */
		CollectPartitionNodes(fProcessor, fComm.Rank(), fPartitionNodes, fExternalNodes);

		/* copy nodal information */
		fNodeManager->CopyNodeToNode(fPBCNodes, fPBCNodes_ghost);
	}
}

/* init data needed for reconfiguring across processors */
void CommManagerT::InitConfigure(iArray2DT& i_values, nVariArray2DT<int>& i_values_man, 
	dArray2DT& new_init_coords, nVariArray2DT<double>& new_init_coords_man,
	dArray2DT& new_curr_coords, nVariArray2DT<double>& new_curr_coords_man)
{
	/* current number of nodes owned by this partition */
	int npn = fPartitionNodes.Length();

	/* size of integer data per node */
	int r = int(sizeof(int)/sizeof(char));
	int nns = fModelManager.NumNodeSets();
	int nsd = fModelManager.NumDimensions();
	int i_offset = 1 + fNodalAttributes.NumRegistered();
	int i_size = i_offset + int(ceil(double(nns)/r)); /* {global ID, attributes, (char) {0|1)}} */

	/* initialize integer data per node */
	i_values_man.SetWard(10, i_values, i_size);
	i_values_man.SetMajorDimension(npn, false);

	/* write global ID and flags - these are length (total processor nodes) and the
	 * first npn are the processor nodes */
	i_values = 0; /* initialize */
	i_values.SetColumn(0, fNodeMap.Pointer()); /* 0: global ID */
	for (int i = 0; i < fNodalAttributes.NumRegistered(); i++)
		i_values.SetColumn(i+1, fNodalAttributes.Array(i).Pointer()); /* i: attribute_i */

	/* node set flags */
	const ArrayT<StringT>& ns_ID = fModelManager.NodeSetIDs();
	for (int i = 0; i < nns; i++) /* mark node set participation */
	{
		const iArrayT& ns = fModelManager.NodeSet(ns_ID[i]);
		for (int j = 0; j < ns.Length(); j++) {
			int node = ns[j];
			if (node < npn) /* processor nodes are 0...(npn-1) */{
				char* ns_flag = (char*)(i_values(node) + i_offset); /* field beyond int values */
				ns_flag[i] = 1; /* node is member */
			}
		}
	}

	/* size of double data per node */
	int pack_size = NodeManager().PackSize();
	int d_size = nsd + nsd + pack_size; /* {X, x, {kinematic data}} */

	/* work space for collecting processor coordinates */
	new_init_coords_man.SetWard(10, new_init_coords, nsd);
	new_init_coords_man.SetMajorDimension(npn, false);
	new_init_coords.BlockRowCopyAt(NodeManager().InitialCoordinates(), 0, npn);

	new_curr_coords_man.SetWard(10, new_curr_coords, nsd);
	new_curr_coords_man.SetMajorDimension(npn, false);
	new_curr_coords.BlockRowCopyAt(NodeManager().CurrentCoordinates(), 0, npn);

	/* initialize exchange buffers - 10% of atoms */
	fd_send_buffer_man.Dimension(npn/10, d_size);
	fd_recv_buffer_man.Dimension(npn/10, d_size);
	fi_send_buffer_man.Dimension(npn/10, i_size);
	fi_recv_buffer_man.Dimension(npn/10, i_size);

	/* initialize dimension of exchange node lists */
	for (int i = 0; i < fSendNodes.Length(); i++) {
		fSendNodes[i].SetHeadRoom(10);
		fSendNodes[i].Dimension(npn/10);
		fRecvNodes[i].SetHeadRoom(10);
		fRecvNodes[i].Dimension(npn/10);		
	}
}

/* finalize configuration */
void CommManagerT::CloseConfigure(iArray2DT& i_values, dArray2DT& new_init_coords)
{
	const char caller[] = "CommManagerT::CloseConfigure";

	/* dimensions */
	int npn = new_init_coords.MajorDim();
	int nns = fModelManager.NumNodeSets();
	
	/* reset model geometry */
	fModelManager.UpdateNodes(new_init_coords, false);
	fNodeManager->UpdateCurrentCoordinates();

	/* reset node map */
	fNodeMap.Dimension(i_values.MajorDim());
	i_values.ColumnCopy(0, fNodeMap.Pointer());

#ifdef CommManagerT_DEBUG_SPATIAL
iArrayT tmp;
tmp.Alias(fNodalAttributes.Array(0));
cout << "types: BEFORE\n" << tmp.wrap(10) << '\n';
#endif
	
	/* restore attributes */
	fNodalAttributes.Dimension(i_values.MajorDim(), false);
	for (int i = 0; i < fNodalAttributes.NumRegistered(); i++)
		i_values.ColumnCopy(i+1, fNodalAttributes.Array(i));

#ifdef CommManagerT_DEBUG_SPATIAL
tmp.Alias(fNodalAttributes.Array(0));
cout << "types: AFTER\n" << tmp.wrap(10) << '\n';
#endif

	/* reset partition nodes */
	CollectPartitionNodes(fProcessor, fComm.Rank(), fPartitionNodes, fExternalNodes);
	fPartitionNodes_inv.ClearMap();

	/* update partition */
	if (fPartition) {

		/* gather info */
		iArrayT part_processor(fPartitionNodes.Length(), fProcessor.Pointer());
		iArrayT part_node_map(fPartitionNodes.Length(), fNodeMap.Pointer());
		iArrayT tmp(fPartitionNodes.Length());
		tmp.SetValueToPosition();
		iArray2DT connects(tmp.Length(), 1, tmp.Pointer());
		ArrayT<const iArray2DT* > connects_1(1);
		connects_1[0] = &connects;
		ArrayT<const RaggedArray2DT<int>* > connects_2;

		/* update nodes */
		fPartition->Set(fComm.Size(), fPartition->ID(), part_processor, part_node_map, 
			connects_1, connects_2);
			
		/* update elements */
		const ArrayT<StringT>& block_ID = fPartition->BlockID();
		if (block_ID.Length() > 1)
			ExceptionT::GeneralFail(caller, "expecting only 1 element block %d", block_ID.Length());
		fPartition->InitElementBlocks(block_ID); /* wipe it */
	}

	/* reset node sets */
	AutoArrayT<int> ns_tmp;
	iArrayT ns;
	const ArrayT<StringT>& ns_ID = fModelManager.NodeSetIDs();
	int i_offset = 1 + fNodalAttributes.NumRegistered();
	for (int i = 0; i < nns; i++) /* mark node set participation */
	{
		/* count numbers of set nodes */
		int count = 0;
		for (int j = 0; j < npn; j++) {
			char* ns_flag = (char*)(i_values(j) + i_offset); /* field beyond int values */
			if (ns_flag[i] == 1) /* node is member */
				count++;
		}
	
		/* collect nodes */
		ns_tmp.Dimension(count);
		count = 0;
		for (int j = 0; j < npn; j++) {
			char* ns_flag = (char*)(i_values(j) + i_offset); /* field beyond int values */
			if (ns_flag[i] == 1) /* node is member */
				ns_tmp[count++] = j;
		}			

		/* reset model */
		ns.Alias(ns_tmp);
		fModelManager.UpdateNodeSet(ns_ID[i], ns, false);
	}
}

/* distribute nodes on spatial grid */
void CommManagerT::Distribute(iArray2DT& i_values, nVariArray2DT<int>& i_values_man, 
	dArray2DT& new_init_coords, nVariArray2DT<double>& new_init_coords_man,
	dArray2DT& new_curr_coords, nVariArray2DT<double>& new_curr_coords_man)
{
	/* dimensions */
	int pack_size = fNodeManager->PackSize();
	int npn = new_init_coords.MajorDim();
	int nsd = new_init_coords.MinorDim();

	/* exhange nodes that are out of bounds */
	dArrayT pack(pack_size), d_alias;
	MPI_Request request;
	int normal[2] = {-1, 1};
	for (int ii = 0; ii < fSwap.MajorDim(); ii++)
//	for (int i = 0 ; i < nsd; i++)
//		for (int j = 0; j < 2; j++) /* shift low and high */
		{
			/* swap info */
			int dir = fSwap(ii,0);
			int sgn = fSwap(ii,1);
		
			/* neighboring process ID's */
			int ID_s = fAdjacentCommID(dir,sgn);
			int ID_r = fAdjacentCommID(dir,1-sgn);
			
			/* has neighboring processes different from self */
			bool do_send = (ID_s != -1 && ID_s != Rank());
			bool do_recv = (ID_r != -1 && ID_r != Rank());

			/* message tags */
			int tag_s = Rank();
			int tag_r = ID_r;

			/* number of nodes to exchange */
			int n_s = 0;
			int n_r = 0;
		
			/* post receive */
			if (do_recv) fComm.PostReceive(n_r, ID_r, tag_r, request);
		
			/* collect outgoing */
			if (do_send) 
			{
				/* reset outgoing buffer size */
				fi_send_buffer_man.SetMajorDimension(0, false);
				fd_send_buffer_man.SetMajorDimension(0, false);
		
				/* loop over (current) processor nodes */
				double bound = fBounds(dir,sgn);
				int k = 0;
				while (k < npn)
				{
					double dxn = (new_curr_coords(k,dir) - bound)*normal[sgn];
					bool out = (sgn == 0) ? (dxn > 0.0) : (dxn >= 0.0); /* on upper boundary is "out" */
					if (out) /* out of bounds */
					{
						/* copy into send buffers */
						n_s++;
						fi_send_buffer_man.SetMajorDimension(n_s, true);
						fi_send_buffer.SetRow(n_s - 1, i_values(k));

						fd_send_buffer_man.SetMajorDimension(n_s, true);
						new_init_coords.RowCopy(k, fd_send_buffer(n_s - 1));
						new_curr_coords.RowCopy(k, fd_send_buffer(n_s - 1) + nsd);
						d_alias.Alias(pack_size, fd_send_buffer(n_s - 1) + 2*nsd);
						fNodeManager->Pack(k, d_alias);

						/* remove node from local database (shuffle down) */
						npn--;
						i_values.CopyRowFromRow(k, npn);
						new_init_coords.CopyRowFromRow(k, npn);
						new_curr_coords.CopyRowFromRow(k, npn);
						fNodeManager->Pack(npn, pack);
						fNodeManager->Unpack(k, pack);
					}
					else
						k++; /* next node */
				}
				
				/* send size */
				fComm.Send(n_s, ID_s, tag_s);

#ifdef CommManagerT_DEBUG_SPATIAL
cout << "s: sending " << n_s << " to " << ID_s << '\n';
#endif /* CommManagerT_DEBUG_SPATIAL */

				/* set flag */
				do_send = (n_s > 0); 
			}

#ifdef CommManagerT_DEBUG_SPATIAL
else
cout << "s: no neighbor in coordinate " << dir << " direction " << sgn << '\n';
#endif /* CommManagerT_DEBUG_SPATIAL */
			
			/* complete and post next receive */
			if (do_recv) 
			{
				/* complete receive */
				fComm.Wait(request);

#ifdef CommManagerT_DEBUG_SPATIAL
cout << "r: receiving " << n_r << " from " << ID_r << '\n';
#endif /* CommManagerT_DEBUG_SPATIAL */
			
				/* post receive for integer data */
				if (n_r > 0)
				{
					/* reset incoming buffer size */
					fi_recv_buffer_man.SetMajorDimension(n_r, false);

					/* post receive */
					fComm.PostReceive(fi_recv_buffer, ID_r, tag_r, request);
				}
				else /* set flag */
					do_recv = false;					
			}
#ifdef CommManagerT_DEBUG_SPATIAL
else
cout << "r: no neighbor in coordinate " << dir << " direction " << sgn << '\n';
#endif /* CommManagerT_DEBUG_SPATIAL */
			
			/* send integer data */
			if (do_send) fComm.Send(fi_send_buffer, ID_s, tag_s);

			/* complete and post next receive */
			if (do_recv) 
			{
				/* complete receive */
				fComm.Wait(request);

				/* reset incoming buffer size */
				fd_recv_buffer_man.SetMajorDimension(n_r, false);

				/* post receive */
				fComm.PostReceive(fd_recv_buffer, ID_r, tag_r, request);
			}

			/* send double data */
			if (do_send) fComm.Send(fd_send_buffer, ID_s, tag_s);

			/* re-dimension workspace */
			i_values_man.SetMajorDimension(npn + n_r, true);
			new_init_coords_man.SetMajorDimension(npn + n_r, true);
			new_curr_coords_man.SetMajorDimension(npn + n_r, true);
			fNodeManager->ResizeNodes(npn + n_r);

			/* process incoming */
			if (do_recv)
			{				
				/* complete receive */
				fComm.Wait(request);

				/* flags */
				i_values.BlockRowCopyAt(fi_recv_buffer, npn, n_r);

				/* coordinates and kinematic data */
				for (int k = 0; k < n_r; k++) {
					new_init_coords.SetRow(npn + k, fd_recv_buffer(k));
					new_curr_coords.SetRow(npn + k, fd_recv_buffer(k) + nsd);
					d_alias.Alias(pack_size, fd_recv_buffer(k) + 2*nsd);
					fNodeManager->Unpack(npn + k, d_alias);
				}
				
				/* update size */
				npn += n_r;
			}
		}		

	/* update the node-to-processor map */
	fProcessor.Dimension(new_curr_coords.MajorDim());
	fProcessor = Rank();
}

/* set border information */
void CommManagerT::SetExchange(iArray2DT& i_values, nVariArray2DT<int>& i_values_man, 
	dArray2DT& new_init_coords, nVariArray2DT<double>& new_init_coords_man,
	dArray2DT& new_curr_coords, nVariArray2DT<double>& new_curr_coords_man)
{
	/* dimensions */
	int pack_size = fNodeManager->PackSize();
	int npn = new_init_coords.MajorDim();
	int nsd = new_init_coords.MinorDim();

	/* reset buffer sizes */
	fi_send_buffer_man.SetMajorDimension(0, false);
	fd_send_buffer_man.SetMajorDimension(0, false);
	fi_recv_buffer_man.SetMajorDimension(0, false);
	fd_recv_buffer_man.SetMajorDimension(0, false);

	/* exhange nodes involved with swaps */
	dArrayT pack(pack_size), d_alias;
	MPI_Request request;
	int normal[2] = {-1, 1};
	for (int i = 0; i < fSwap.MajorDim(); i++)
	{
		/* swap info */
		int dir = fSwap(i,0);
		int sgn = fSwap(i,1);

		/* neighboring process ID's */
		int ID_s = fAdjacentCommID(dir,sgn);
		int ID_r = fAdjacentCommID(dir,1-sgn);
			
		/* has neighboring processes different from self */
		bool do_send = (ID_s != -1 && ID_s != Rank());
		bool do_recv = (ID_r != -1 && ID_r != Rank());

		/* message tags */
		int tag_s = Rank();
		int tag_r = ID_r;

		/* number of nodes to exchange */
		int n_s = 0;
		int n_r = 0;

		/* exchange lists */
		AutoArrayT<int>& send_nodes = fSendNodes[i];
		AutoArrayT<int>& recv_nodes = fRecvNodes[i];
		send_nodes.Dimension(0);
		recv_nodes.Dimension(0);

		/* post receive of incoming dimension */
		if (do_recv) fComm.PostReceive(n_r, ID_r, tag_r, request);

		/* collect nodes */
		if (do_send)
		{
			double bound = fBounds(dir,sgn);
			for (int k = 0; k < npn; k++)
			{
				double dxn = -(new_curr_coords(k,dir) - bound)*normal[sgn];
				bool out_low = (sgn == 0) ? (dxn >= 0.0) : (dxn > 0.0); /* on lower boundary is "sent" */
				if (out_low && dxn < fSkin)
				{
					/* copy into send buffers */
					n_s++;
					fi_send_buffer_man.SetMajorDimension(n_s, true);
					fi_send_buffer.SetRow(n_s - 1, i_values(k));

					fd_send_buffer_man.SetMajorDimension(n_s, true);
					new_init_coords.RowCopy(k, fd_send_buffer(n_s - 1));
					new_curr_coords.RowCopy(k, fd_send_buffer(n_s - 1) + nsd);
					d_alias.Alias(pack_size, fd_send_buffer(n_s - 1) + 2*nsd);
					fNodeManager->Pack(k, d_alias);
					
					/* add to list */
					send_nodes.Append(k);
				}
			}

			/* send size */
			fComm.Send(n_s, ID_s, tag_s);
				
			/* set flag */
			do_send = (n_s > 0);
		}

		/* receive dimension and post next receive */
		if (do_recv)
		{
			/* complete receive */
			fComm.Wait(request);
		
			/* post receive for the integer data */
			if (n_r > 0)
			{
				/* reset incoming buffer size */
				fi_recv_buffer_man.SetMajorDimension(n_r, false);

				/* post receive */
				fComm.PostReceive(fi_recv_buffer, ID_r, tag_r, request);
			}
			else /* set flag */
				do_recv = false;
		}

		/* send integer data */
		if (do_send) fComm.Send(fi_send_buffer, ID_s, tag_s);

		/* complete and post next receive */
		if (do_recv) 
		{
			/* complete receive */
			fComm.Wait(request);

			/* reset outgoing buffer size */
			fd_recv_buffer_man.SetMajorDimension(n_r, false);

			/* post receive */
			fComm.PostReceive(fd_recv_buffer, ID_r, tag_r, request);
		}

		/* send double data */
		if (do_send) fComm.Send(fd_send_buffer, ID_s, tag_s);

		/* process incoming */
		if (do_recv)
		{				
			/* complete receive */
			fComm.Wait(request);

			/* flags */
			i_values_man.SetMajorDimension(npn + n_r, true);
			i_values.BlockRowCopyAt(fi_recv_buffer, npn, n_r);
			
			/* coordinates and kinematic data */
			new_init_coords_man.SetMajorDimension(npn + n_r, true);
			new_curr_coords_man.SetMajorDimension(npn + n_r, true);
			fNodeManager->ResizeNodes(npn + n_r);
			recv_nodes.Dimension(n_r);
			for (int k = 0; k < n_r; k++) {
				new_init_coords.SetRow(npn + k, fd_recv_buffer(k));
				new_curr_coords.SetRow(npn + k, fd_recv_buffer(k) + nsd);
				d_alias.Alias(pack_size, fd_recv_buffer(k) + 2*nsd);
				fNodeManager->Unpack(npn + k, d_alias);
				
				/* add to incoming list */
				recv_nodes[k] = npn + k;
				
				/* update the node-to-processor map */
				fProcessor.Append(ID_r);
				/* NOTE: this isn't really the processor which owns the incoming node, but
				 *       instead is where this processor is getting information about the
				 *       node, which means that the same node can appear to be from different
				 *       processors. Trouble? */
			}
			
			/* nodes on processor */
			npn += n_r;
		}
	}
}

/* determine the coordinate bounds of this processor */
void CommManagerT::GetProcessorBounds(const dArray2DT& coords, dArray2DT& bounds, iArray2DT& adjacent_ID) const
{
	const char caller[] = "CommManagerT::GetProcessorBounds";
	if (!fPartition) ExceptionT::GeneralFail(caller, "partition not set");

	/* dimensions */
	int nsd = coords.MinorDim();
	bounds.Dimension(nsd, 2);
	adjacent_ID.Dimension(nsd, 2);

	/* spatial grid */
	SpatialGridT grid;
	grid.Dimension(fPartition->GridDimensions());

	/* determine global bounds */
	grid.SetBounds(fComm, coords, NULL);
	
	/* grid position of this processor */
	iArrayT grid_pos = fPartition->GridPosition();

	/* bounds of this grid */
	grid.GridBounds(grid_pos, bounds);

	/* ID's of surrounding processors */
	for (int i = 0; i < nsd; i++) {
		int& cell = grid_pos[i];
		cell--; /* low side */
		adjacent_ID(i,0) = grid.Grid2Processor(grid_pos);
		cell += 2; /* high side */
		adjacent_ID(i,1) = grid.Grid2Processor(grid_pos);
		cell--; /* restore */
	}
}
