/* $Id: CommManagerT.cpp,v 1.9.32.2 2004-11-08 02:16:06 d-farrell2 Exp $ */
#include "CommManagerT.h"
#include "CommunicatorT.h"
#include "ModelManagerT.h"
#include "NodeManagerT.h"
#include "PartitionT.h"
#include "InverseMapT.h"
#include "FieldT.h"
#include <float.h>

/* message types */
#include "AllGatherT.h"
#include "PointToPointT.h"

using namespace Tahoe;

CommManagerT::CommManagerT(CommunicatorT& comm, ModelManagerT& model_manager):
	fComm(comm),
	fModelManager(model_manager),
	fSize(fComm.Size()),
	fRank(fComm.Rank()),
	fPartition(NULL),
	fNodeManager(NULL),
	fFirstConfigure(true),
	fNumRealNodes(0)
{
	//TEMP - with inline coordinate information (serial) this map
	//       need to be set ASAP
	if (fComm.Size() == 1)
	{
		fProcessor.Dimension(fModelManager.NumNodes());
		fProcessor = fComm.Rank();
	}
}

CommManagerT::~CommManagerT(void)
{
	/* free any remaining communications */
	for (int i = 0; i < fCommunications.Length(); i++)
		delete fCommunications[i];

	/* free any remaining ghost communications */
	for (int i = 0; i < fGhostCommunications.Length(); i++)
		delete fGhostCommunications[i];
}

/* set or clear partition information */
void CommManagerT::SetPartition(PartitionT* partition)
{
	/* pointer to the partition info */
	fPartition = partition;
	
	/* fetch partition information */
	if (fPartition)
	{
		/* node map */
		fNodeMap.Alias(fPartition->NodeMap());
		
		/* nodes owned by this processor */
		fProcessor.Dimension(fNodeMap.Length());
		fPartition->ReturnProcessorMap(fProcessor);
		CollectPartitionNodes(fProcessor, fComm.Rank(), fPartitionNodes);
		
		/* external nodes */
		fExternalNodes.Alias(fPartition->Nodes_External());
		
		/* border nodes */
		fBorderNodes.Alias(fPartition->Nodes_Border());
	}
	else /* clear partition information */
	{
		/* node map */
		fNodeMap.Dimension(0);
		
		/* nodes owned by this processor */
		fPartitionNodes.Dimension(0);
		
		/* clear external nodes */
		fExternalNodes.Dimension(0);

		/* clear external nodes */
		fBorderNodes.Dimension(0);
	}
	
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

/* enforce the periodic boundary conditions */
void CommManagerT::EnforcePeriodicBoundaries(double skin)
{
	const char caller[] = "CommManagerT::EnforcePeriodicBoundaries";
	if (!fNodeManager) ExceptionT::GeneralFail(caller, "node manager node set");

	/* only implemented for atom decomposition (or serial) */
	if (fPartition && fPartition->DecompType() != PartitionT::kIndex) return;
	
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
			AllGatherT all_gather(fComm);
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

		/* create partition nodes list if */
		if (fPartitionNodes.Length() == 0)
			CollectPartitionNodes(fProcessor, fComm.Rank(), fPartitionNodes);

		/* copy nodal information */
		fNodeManager->CopyNodeToNode(fPBCNodes, fPBCNodes_ghost);
	}
}

/* configure local environment */
void CommManagerT::Configure(void)
{
	/* first time through */
	if (fFirstConfigure) {
		FirstConfigure();
		fFirstConfigure = false;
	}

	/* determine the current processor bounds */
	//GetBounds(current_coords, , fBounds);
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
		nArrayT<int> nd;
		nd.Alias(*nodes);
		non_const_this->fPartitionNodes_inv.SetMap(nd);
		
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

	PartitionT::DecompTypeT decomp_type = fPartition->DecompType();
	if (decomp_type == PartitionT::kGraph) /* non-blocking send-receive */
	{
		/* new point-to-point */
		PointToPointT* p2p = new PointToPointT(fComm, *fPartition);

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
		AllGatherT* all_gather = new AllGatherT(fComm);

		/* set message size */
		all_gather->Initialize(fPartition->NumPartitionNodes()*num_vals);
		
		/* store */
		fCommunications.Append(all_gather);

		/* new all gather for gh */
		AllGatherT* ghost_all_gather = new AllGatherT(fComm);

		/* set message size */
		ghost_all_gather->Initialize(fPBCNodes.Length()*num_vals);
		
		/* store */
		fGhostCommunications.Append(ghost_all_gather);
	}
	else if (decomp_type == PartitionT::kIndex) /* shifts */
	{
		ExceptionT::GeneralFail(caller, "not implemented yet");
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized decomp type: %d", decomp_type);

	/* ID is just the array position */
	return fCommunications.Length() - 1;
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
		fdExchange.Set(fNumRealNodes, values.MinorDim(), values(0));
		all_gather->AllGather(fdExchange);

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
			fdExchange.Set(fModelManager.NumNodes() - fNumRealNodes, values.MinorDim(), values(fNumRealNodes));
			all_gather->AllGather(fdExchange);
		}		
	}
	else if (decomp_type == PartitionT::kSpatial) /* shifts */
	{
		ExceptionT::GeneralFail(caller, "not implemented yet");
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
		fiExchange.Set(fNumRealNodes, values.MinorDim(), values(0));
		all_gather->AllGather(fiExchange);
		
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
			fiExchange.Set(fModelManager.NumNodes() - fNumRealNodes, values.MinorDim(), values(fNumRealNodes));
			all_gather->AllGather(fiExchange);
		}		
	}
	else if (decomp_type == PartitionT::kSpatial) /* shifts */
	{
		ExceptionT::GeneralFail(caller, "not implemented yet");
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized decomp type: %d", decomp_type);
}

/*************************************************************************
 * Private
 *************************************************************************/

/* collect partition nodes */
void CommManagerT::CollectPartitionNodes(const ArrayT<int>& n2p_map, int part, 
	AutoArrayT<int>& part_nodes) const
{
	int count = 0;
	const int* p = n2p_map.Pointer();
	int len = n2p_map.Length();
	for (int i = 0; i < len; i++)
		if (*p++ == part)
			count++;
			
	part_nodes.Dimension(count);
	int* pn = part_nodes.Pointer();
	p = n2p_map.Pointer();
	count = 0;
	for (int i = 0; i < len; i++)
		if (*p++ == part)
			*pn++ = i;
}

/* perform actions needed the first time CommManagerT::Configure is called. */
void CommManagerT::FirstConfigure(void)
{
	/* serial */
	if (!fPartition) {
		fNumRealNodes = fModelManager.NumNodes();
		fProcessor.Dimension(fNumRealNodes);
		fProcessor = fComm.Rank();
		return;
	}
	
	/* graph decomposition - nothing to do? */
	if (fPartition->DecompType() == PartitionT::kGraph)
		return;
	
	if (fPartition->DecompType() == PartitionT::kIndex)
	{
		/* change the numbering scope of partition data to global */
		fPartition->SetScope(PartitionT::kGlobal);
	
		/* number of nodes owned by this partition */
		int npn = fPartition->NumPartitionNodes();
		int nsd = fModelManager.NumDimensions();

		/* total number of nodes */
		AllGatherT all_gather(fComm);
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
			AllGatherT gather_node_set(fComm);
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
		CollectPartitionNodes(fProcessor, fComm.Rank(), fPartitionNodes);
		
		/* clear the inverse map */
		fPartitionNodes_inv.ClearMap();
		
		/* (potentially) all are adjacent to external nodes */
		fBorderNodes.Alias(fPartitionNodes);
		
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
			
		// get the first node owned by this partition (for start/end points of AddScaled/AddCombination
		if (fPartitionNodes.Length() > 0) {
			fPartStartNum = fPartitionNodes.First();
			fPartEndNum = fPartitionNodes.Last();
			fPartFieldStart = (nsd*fPartStartNum) - nsd; // start point in field arrays
			fPartFieldEnd = (nsd*fPartEndNum) - nsd; // end point in field arrays
			
			// check that (End - Start) + 1 is the same as the number of nodes in partition
			if ((fPartEndNum - fPartStartNum + 1) != npn )
			{
				cout << "\n problem with partition bounds \n" << endl; 
			}
		}
	}
}

/** determine the local coordinate bounds */
void CommManagerT::GetBounds(const dArray2DT& coords_all, const iArrayT& local, 
	dArray2DT& bounds) const
{
	/* dimensions */
	int nsd = coords_all.MinorDim();
	bounds.Dimension(nsd, 2);

	/* initialize */
	if (local.Length() == 0)
	{
		bounds = 0.0;
		return;
	}
	else
	{
		bounds.SetColumn(0, DBL_MAX);
		bounds.SetColumn(1, DBL_MIN);
	}

	/* resolve by spatial dimension */
	if (nsd == 1)
	{
		const int* p_i = local.Pointer();
		double& x_min = bounds(0,0);
		double& x_max = bounds(0,1);
		for (int i = 0; i < local.Length(); i++)
		{
			const double* x = coords_all(*p_i++);
			x_min = (*x < x_min) ? *x : x_min;
			x_max = (*x > x_max) ? *x : x_max;
		}
	}
	else if (nsd == 2)
	{
		const int* p_i = local.Pointer();
		double& x_min = bounds(0,0);
		double& x_max = bounds(0,1);
		double& y_min = bounds(1,0);
		double& y_max = bounds(1,1);
		for (int i = 0; i < local.Length(); i++)
		{
			const double* x = coords_all(*p_i++);
			x_min = (*x < x_min) ? *x : x_min;
			x_max = (*x > x_max) ? *x : x_max;
			x++;
			y_min = (*x < y_min) ? *x : y_min;
			y_max = (*x > y_max) ? *x : y_max;
		}
	}
	else if (nsd == 3)
	{
		const int* p_i = local.Pointer();
		double& x_min = bounds(0,0);
		double& x_max = bounds(0,1);
		double& y_min = bounds(1,0);
		double& y_max = bounds(1,1);
		double& z_min = bounds(2,0);
		double& z_max = bounds(2,1);
		for (int i = 0; i < local.Length(); i++)
		{
			const double* x = coords_all(*p_i++);
			x_min = (*x < x_min) ? *x : x_min;
			x_max = (*x > x_max) ? *x : x_max;
			x++;
			y_min = (*x < y_min) ? *x : y_min;
			y_max = (*x > y_max) ? *x : y_max;
			x++;
			z_min = (*x < z_min) ? *x : z_min;
			z_max = (*x > z_max) ? *x : z_max;
		}
	}
	else ExceptionT::OutOfRange("CommManagerT::SetBounds");
}
