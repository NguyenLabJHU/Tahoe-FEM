/* $Id: MappedPeriodicT.cpp,v 1.8 2004-06-17 07:41:57 paklein Exp $ */
/* created: paklein (04/07/1997) */

#include "MappedPeriodicT.h"

#include "NodeManagerT.h"
#include "FEManagerT.h"
#include "ifstreamT.h"
#include "BasicFieldT.h"

/* column indeces */

using namespace Tahoe;

const int kMaster = 0;
const int kSlave  = 1;

/* constructor */
MappedPeriodicT::MappedPeriodicT(NodeManagerT& node_manager, BasicFieldT& field):
	KBC_ControllerT(node_manager),
	fField(field),
	fnumLTf(-1),
	fSchedule(NULL),
	fFperturb(fNodeManager.NumSD()),
	fF(fNodeManager.NumSD()),
	fD_sm(fNodeManager.NumSD()),
	fDummySchedule(1.0)
{
	fF.Identity();
}

/* initialize data - called immediately after construction */
void MappedPeriodicT::Initialize(ifstreamT& in)
{
	/* schedule for fFperturb */
	in >> fnumLTf; fnumLTf--;
	if (fnumLTf < 0) throw ExceptionT::kBadInputValue;
	fSchedule = fNodeManager.Schedule(fnumLTf);	
	if (!fSchedule) throw ExceptionT::kBadInputValue;

	/* specified deformation gradient */
	in >> fFperturb;
	if (!in.good()) throw ExceptionT::kBadInputValue;

	/* list of mapped nodes */
	ArrayT<StringT> id_list;
	ReadNodes(in, id_list, fMappedNodeList);

	/* read master nodes */
	iArrayT tmp;
	ReadNodes(in, id_list, tmp);
	fSlaveMasterPairs.Dimension(tmp.Length(), 2);
	fSlaveMasterPairs.SetColumn(kMaster, tmp);

	/* read corresponding slave nodes */
	ReadNodes(in, id_list, tmp);
	if (tmp.Length() != fSlaveMasterPairs.MajorDim())
	{
		cout << "\n MappedPeriodicT::Initialize: length of master node list "
		     << fSlaveMasterPairs.MajorDim() << " does\n"
		     <<   "     not match the length of the slave node list "
		     << tmp.Length() << endl;
		throw ExceptionT::kBadInputValue;
	}
	fSlaveMasterPairs.SetColumn(kSlave, tmp);
	
	/* generate BC cards */
	int num_BC = fMappedNodeList.Length() + fSlaveMasterPairs.MajorDim();
	int nsd = fFperturb.Rows();
	fKBC_Cards.Dimension(num_BC*nsd);
	fMappedCards.Set(fMappedNodeList.Length()*nsd, fKBC_Cards.Pointer());
	fSlaveCards.Set(fSlaveMasterPairs.MajorDim()*nsd,
		fKBC_Cards.Pointer(fMappedCards.Length()));

	/* mapped nodes */
	int dex = 0;
	for (int i = 0; i < fMappedNodeList.Length(); i++)
		for (int j = 0; j < nsd; j++)
			fMappedCards[dex++].SetValues(fMappedNodeList[i], j, KBC_CardT::kDsp, 0, 0.0);	

	/* slave nodes */
	dex = 0;
	for (int ii = 0; ii < fSlaveMasterPairs.MajorDim(); ii++)
		for (int jj = 0; jj < nsd; jj++)
		{
			/* set values */
			fSlaveCards[dex].SetValues(fSlaveMasterPairs(ii, kSlave), jj, KBC_CardT::kDsp, 0, 0.0);
	
			/* dummy schedule */
			fSlaveCards[dex].SetSchedule(&fDummySchedule);
			dex++;
		}	
}

void MappedPeriodicT::WriteParameters(ostream& out) const
{
	/* inherited */
	KBC_ControllerT::WriteParameters(out);

	iArrayT tmp;
	out << "\n Mapping Parameters:\n";
	out << " Mapping load time function. . . . . . . . . . . = " << fnumLTf << '\n';
	out << " Mapping perturbation:\n";
	out << fFperturb << '\n';
	out << " Number of mapped nodes. . . . . . . . . . . . . = "
	    << fMappedNodeList.Length() << '\n';
	tmp.Alias(fMappedNodeList);
	tmp++;
	out << tmp.wrap(5) << '\n';
	tmp--;
	out << " Number of linked node pairs . . . . . . . . . . = "
	    << fSlaveMasterPairs.MajorDim() << '\n';
	tmp.Alias(fSlaveMasterPairs);
	tmp++;
	fSlaveMasterPairs.WriteNumbered(out);
	tmp--;
	out << '\n';
}

/* initial condition */
void MappedPeriodicT::InitialCondition(void)
{
	/* reference coordinates */
	const dArray2DT& init_coords = fNodeManager.InitialCoordinates();

	/* compute mapping */
	fF = fFperturb;
	fF.PlusIdentity();

	int nsd = fF.Rows();
	dArrayT	X;
	dArrayT d(nsd);
	int mappedcount = fMappedNodeList.Length();
	int dex = 0;
	for (int i = 0; i < mappedcount; i++)
	{
		int node = fMappedNodeList[i];
	
		/* fetch pointers */
		init_coords.RowAlias(node, X);
	
		/* map current coord */
		fF.Multx(X, d);
		
		/* set prescribed displacements */
		for (int j = 0; j < nsd; j++)
		{
			fMappedCards[dex].SetValues(node, j, KBC_CardT::kDsp, 0, d[j] - X[j]);
			fMappedCards[dex].SetSchedule(fSchedule);
			dex++;
		}
	}
}

/* set BC cards for current step */
void MappedPeriodicT::InitStep(void)
{
	/* inherited */
	KBC_ControllerT::InitStep();
	
	/* just compute mapping */
	if (fSlaveMasterPairs.MajorDim() == 0)
	{
		/* compute F - 1  = fFperturb */
		fF.SetToScaled(fSchedule->Value(), fFperturb);
		fF.PlusIdentity();
	}
	else /* apply mapping */
	{
		/* nodal information */
		const dArray2DT& init_coords = fNodeManager.InitialCoordinates();
		const dArray2DT& disp = fField[0];	

		/* compute F - 1  = fFperturb */
		fF.SetToScaled(fSchedule->Value(), fFperturb);

		int dex = 0;
		int nsd = fF.Rows();
		int num_pairs = fSlaveMasterPairs.MajorDim();	
		dArrayT	X_m, d_m;      //master data
		dArrayT X_s, d_s(nsd); //slave data
		for (int i = 0; i < num_pairs; i++)
		{
			int node_m = fSlaveMasterPairs(i, kMaster);
			int node_s = fSlaveMasterPairs(i, kSlave);

			/* fetch pointers */
			init_coords.RowAlias(node_m, X_m);
			disp.RowAlias(node_m, d_m);		
			init_coords.RowAlias(node_s, X_s);

			/* periodic displacements */		
			fD_sm.DiffOf(X_s, X_m);
			fF.Multx(fD_sm, d_s);
		
			/* set cards */
			for (int j = 0; j < nsd; j++)
				fSlaveCards[dex++].SetValues(node_s, j, KBC_CardT::kDsp, 0, d_s[j] + d_m[j]);
		}
	
		/* set mapping */
		fF.PlusIdentity();
	}
}

/* output current configuration */
void MappedPeriodicT::WriteOutput(ostream& out) const
{
	/* inherited */
	KBC_ControllerT::WriteOutput(out);

	/* write output */
	out << "\n Mapping:\n";
	out << fF << endl;
}
