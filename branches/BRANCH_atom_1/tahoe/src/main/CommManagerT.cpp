/* $Id: CommManagerT.cpp,v 1.1.2.10 2003-01-11 01:18:41 paklein Exp $ */
#include "CommManagerT.h"
#include "CommunicatorT.h"
#include "ModelManagerT.h"
#include "NodeManagerT.h"
#include "PartitionT.h"
#include "InverseMapT.h"
#include <float.h>

/* message types */
#include "AllGatherT.h"
#include "PointToPointT.h"

using namespace Tahoe;

CommManagerT::CommManagerT(CommunicatorT& comm, ModelManagerT& model_manager):
	fComm(comm),
	fModelManager(model_manager),
	fPartition(NULL),
	fNodeManager(NULL),
	fPeriodicBoundaries(0,2),
	fFirstConfigure(true)
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
}

/* set boundaries */
void CommManagerT::SetPeriodicBoundaries(int i, double x_i_min, double x_i_max)
{
	/* dimension dynamically */
	if (i >= fIsPeriodic.Length())
	{
		int len = i + 1;
		fIsPeriodic.Resize(len, false);
		fPeriodicBoundaries.Resize(len, 0.0);
		fPeriodicLength.Resize(len, 0.0);
	}

	if (x_i_min > x_i_max) ExceptionT::GeneralFail();
	fIsPeriodic[i] = true;
	fPeriodicBoundaries(i,0) = x_i_min;
	fPeriodicBoundaries(i,1) = x_i_max;
	fPeriodicLength[i] = x_i_max - x_i_min;

//TEMP
ExceptionT::Stop("CommManagerT::SetPeriodicBoundaries", "not supported yet");
}

/* unset boundaries */
void CommManagerT::ClearPeriodicBoundaries(int i) { fIsPeriodic[i] = false; }

/* enforce the periodic boundary conditions */
void CommManagerT::EnforcePeriodicBoundaries(dArray2DT& displacement, double skin)
{
	const char caller[] = "CommManagerT::EnforcePeriodicBoundaries";
	if (!fNodeManager) ExceptionT::GeneralFail(caller, "node manager node set");

	/* only implemented for atom decomposition (or serial) */
	if (fPartition && fPartition->DecompType() != PartitionT::kAtom) return;

	/* reference coordinates */
	const dArray2DT& reference_coords = fModelManager.Coordinates();
	if (reference_coords.MinorDim() > fIsPeriodic.Length()) ExceptionT::SizeMismatch(caller);

	/* nodes owned by this partition */
	const ArrayT<int>* partition_nodes = PartitionNodes();
	int nnd = (partition_nodes) ? partition_nodes->Length() : reference_coords.MajorDim();

	/* number of "real" nodes */
	int ngn = fComm.Sum(fPBCNodes.Length());
	int nrn = fModelManager.NumNodes() - ngn;

	/* reset ghost nodes */
	fPBCNodes.Dimension(0);
	fPBCNodes_ghost.Dimension(0);
	fPBCNodes_face.Dimension(0);

	/* loop over directions */
	for (int i = 0; i < fIsPeriodic.Length(); i++)
		if (fIsPeriodic[i])
		{		
			/* coordinate limits */
			double x_min = fPeriodicBoundaries(i,0);
			double x_max = fPeriodicBoundaries(i,1);
			double x_len = x_max - x_min;

			/* nodes owned by this partition */
			for (int j = 0; j < nnd; j++)
			{
				/* current coordinates */
				int nd = (partition_nodes) ? (*partition_nodes)[j] : j;
				double& X = reference_coords(nd,i);
				double& d = displacement(nd,i);
				double  x = X + d;
			
				/* shift displacements - back in the box */
				if (x > x_max) 
					d -= x_len;
				else if (x < x_min)
					d += x_len;
					
				/* collect nodes close to periodic boundaries */
				x = X + d;
				if (x - x_min < skin) {
					fPBCNodes.Append(nd);
					fPBCNodes_face.Append(-(i+1));
				}
				if (x_max - x < skin) {
					fPBCNodes.Append(nd);
					fPBCNodes_face.Append(i+1);
				}
			}
		}
		
	/* configure ghost nodes */
	if (fPBCNodes.Length() > 0)
	{
		/* communicate number of ghost nodes per partition */
		iArrayT ghost_count(fComm.Size());
		fComm.AllGather(fPBCNodes.Length(), ghost_count);

		/* local numbering of ghost nodes */
		int ghost_num_start = nrn;
		for (int i = 0; i < fComm.Rank(); i++)
			ghost_num_start += ghost_count[i];

		fPBCNodes_ghost.Dimension(fPBCNodes.Length());
		int ghost_num = ghost_num_start;
		for (int i = 0; i < fPBCNodes_ghost.Length(); i++)
			fPBCNodes_ghost[i] = ghost_num++;

		/* resize all coordinate and field arrays */
		fNodeManager->ResizeNodes(nrn + fPBCNodes_ghost.Length());

		/* copy nodal information */
		fNodeManager->CopyNodeToNode(fPBCNodes, fPBCNodes_ghost);

		/* generate reference coordinates - write access to the coordinates from the
		 * model manager, not from the node manager. The node manager does not change
		 * the initial coordinates during CopyNodeToNode */
		const dArray2DT& coords = fModelManager.Coordinates();
		dArray2DT ghost_coords(fPBCNodes_ghost.Length(), coords.MinorDim(), coords(ghost_num_start));
		for (int i = 0; i < fPBCNodes_ghost.Length(); i++)
		{
			/* copy coordinates */
			ghost_coords.SetRow(i, coords(fPBCNodes[i]));
		
			/* periodic shift */
			int face = fPBCNodes_face[i];
			int dim = int(fabs(double(face))) - 1;
			int shift = (face > 0) ? -1 : 1;
			ghost_coords(i,dim) += shift*fPeriodicLength[dim];
		}

		/* exchange */
		AllGatherT all_gather(fComm);
		all_gather.Initialize(ghost_coords);
		all_gather.AllGather(ghost_coords);

		/* persistent communications */
		for (int i = 0; i < fGhostCommunications.Length(); i++)
		{
			/* retrieve communications pointer */
			AllGatherT* all_gather = dynamic_cast<AllGatherT*>(fGhostCommunications[i]);
			if (!all_gather) ExceptionT::GeneralFail(caller, "could not recover ghost communication");

			/* (re-)set message size */
			all_gather->Initialize(fPBCNodes.Length()*fNumValues[i]);
		}
	}
}

/* configure local environment */
void CommManagerT::Configure(void)
{
	/* nothing to do */
	if (fComm.Size() == 1)
	{
		bool has_periodic = false;
		for (int i = 0; i < fIsPeriodic.Length(); i++)
			has_periodic = has_periodic || fIsPeriodic[i];
		if (!has_periodic) return;
	}
	if (!fPartition || fPartition->DecompType() != PartitionT::kAtom) return;

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
	else if (decomp_type == PartitionT::kAtom) /* all gather */
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
		all_gather->Initialize(fPBCNodes.Length()*num_vals);
		
		/* store */
		fGhostCommunications.Append(all_gather);
	}
	else if (decomp_type == PartitionT::kAtom) /* shifts */
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
		PointToPointT* p2p = dynamic_cast<PointToPointT*>(fCommunications[id]);
		if (!p2p) ExceptionT::GeneralFail(caller);

		/* do it */
		p2p->AllGather(values);
	}
	else if (decomp_type == PartitionT::kAtom) /* all gather */
	{
		/* retrieve pointer */
		AllGatherT* all_gather = dynamic_cast<AllGatherT*>(fCommunications[id]);
		if (!all_gather) ExceptionT::GeneralFail(caller);
		
		/* exchange */
		fdExchange.Set(fPartition->NumPartitionNodes(), values.MinorDim(), values.Pointer());
		all_gather->AllGather(fdExchange);

		/* ghost nodes */
		if (fPBCNodes.Length() > 0) {
		
			/* copy values to ghost nodes */
			int ngh = fPBCNodes.Length();
			for (int i = 0; i < ngh; i++)
				values.SetRow(fPBCNodes_ghost[i], values(fPBCNodes[i]));

			/* retrieve communications pointer */
			AllGatherT* all_gather = dynamic_cast<AllGatherT*>(fGhostCommunications[id]);
			if (!all_gather) ExceptionT::GeneralFail(caller);

			/* exchange */
			int nrn = fModelManager.NumNodes() - fPBCNodes.Length();
			fdExchange.Set(fPBCNodes.Length(), values.MinorDim(), values(nrn));
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
		PointToPointT* p2p = dynamic_cast<PointToPointT*>(fCommunications[id]);
		if (!p2p) ExceptionT::GeneralFail(caller);

		/* do it */
		p2p->AllGather(values);
	}
	else if (decomp_type == PartitionT::kAtom) /* all gather */
	{
		/* retrieve communications pointer */
		AllGatherT* all_gather = dynamic_cast<AllGatherT*>(fCommunications[id]);
		if (!all_gather) ExceptionT::GeneralFail(caller);
		
		/* exchange */
		fiExchange.Set(fPartition->NumPartitionNodes(), values.MinorDim(), values.Pointer());
		all_gather->AllGather(values);
		
		/* ghost nodes */
		if (fPBCNodes.Length() > 0) {
		
			/* copy values to ghost nodes */
			int ngh = fPBCNodes.Length();
			for (int i = 0; i < ngh; i++)
				values.SetRow(fPBCNodes_ghost[i], values(fPBCNodes[i]));

			/* retrieve communications pointer */
			AllGatherT* all_gather = dynamic_cast<AllGatherT*>(fGhostCommunications[id]);
			if (!all_gather) ExceptionT::GeneralFail(caller);

			/* exchange */
			int nrn = fModelManager.NumNodes() - fPBCNodes.Length();
			fiExchange.Set(fPBCNodes.Length(), values.MinorDim(), values(nrn));
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
	int* p = n2p_map.Pointer();
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
	/* nothing to do in serial */
	if (!fPartition) {
		fProcessor.Dimension(fModelManager.NumNodes());
		fProcessor = fComm.Rank();
		return;
	}
	
	//TEMP
	if (fPartition->DecompType() == PartitionT::kGraph)
		ExceptionT::GeneralFail("CommManagerT::FirstConfigure", "not yet implemented for kGraph");
	
	if (fPartition->DecompType() == PartitionT::kAtom)
	{
		/* change the numbering scope of partition data to global */
		fPartition->SetScope(PartitionT::kGlobal);
	
		/* number of nodes owned by this partition */
		int npn = fPartition->NumPartitionNodes();

		/* total number of nodes */
		int nsd = fModelManager.NumDimensions();
		AllGatherT all_gather(fComm);
		all_gather.Initialize(npn*nsd);
		int ntn = all_gather.Total()/nsd;

		/* collect global reference coordinate list */
		const dArray2DT& reference_coords = fModelManager.Coordinates();
		if (reference_coords.MajorDim() != ntn)
		{
			/* collect */
			dArray2DT coords_all(ntn, nsd);		
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
		fExternalNodes.Dimension(ntn - fPartitionNodes.Length());
		int lower, upper;
		lower = upper = ntn;
		if (fPartitionNodes.Length() > 0) {
			lower = fPartitionNodes.First();
			upper = fPartitionNodes.Last();
		}
		int dex = 0;
		for (int i = 0; i < lower; i++)
			fExternalNodes[dex++] = i;
		for (int i = upper+1; i < ntn; i++)
			fExternalNodes[dex++] = i;
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
		int* p_i = local.Pointer();
		double& x_min = bounds(0,0);
		double& x_max = bounds(0,1);
		for (int i = 0; i < local.Length(); i++)
		{
			double* x = coords_all(*p_i++);
			x_min = (*x < x_min) ? *x : x_min;
			x_max = (*x > x_max) ? *x : x_max;
		}
	}
	else if (nsd == 2)
	{
		int* p_i = local.Pointer();
		double& x_min = bounds(0,0);
		double& x_max = bounds(0,1);
		double& y_min = bounds(1,0);
		double& y_max = bounds(1,1);
		for (int i = 0; i < local.Length(); i++)
		{
			double* x = coords_all(*p_i++);
			x_min = (*x < x_min) ? *x : x_min;
			x_max = (*x > x_max) ? *x : x_max;
			x++;
			y_min = (*x < y_min) ? *x : y_min;
			y_max = (*x > y_max) ? *x : y_max;
		}
	}
	else if (nsd == 3)
	{
		int* p_i = local.Pointer();
		double& x_min = bounds(0,0);
		double& x_max = bounds(0,1);
		double& y_min = bounds(1,0);
		double& y_max = bounds(1,1);
		double& z_min = bounds(2,0);
		double& z_max = bounds(2,1);
		for (int i = 0; i < local.Length(); i++)
		{
			double* x = coords_all(*p_i++);
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
