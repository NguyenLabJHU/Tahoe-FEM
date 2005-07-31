/* $Id: KBC_ControllerT.cpp,v 1.1.1.1 2001-01-29 08:20:40 paklein Exp $ */
/* created: paklein (09/05/2000)                                          */

#include "KBC_ControllerT.h"

#include "NodeManagerT.h"
#include "FEManagerT.h"
#include "ModelFileT.h"
#include "ExodusT.h"
#include "fstreamT.h"

/* constructor */
KBC_ControllerT::KBC_ControllerT(NodeManagerT& node_manager):
	fNodeManager(node_manager)
{

}

/* destructor */
KBC_ControllerT::~KBC_ControllerT(void) { }

/* initialization */
void KBC_ControllerT::WriteParameters(ostream& out) const
{
#pragma unused(out)
}

void KBC_ControllerT::ReadRestart(istream& in)
{
#pragma unused(in)
}

void KBC_ControllerT::WriteRestart(ostream& out) const
{
#pragma unused(out)
}

/* initialize/finalize step */
void KBC_ControllerT::InitStep(void) { }
void KBC_ControllerT::CloseStep(void) { }
void KBC_ControllerT::Reset(void) { }

/* returns true if the internal force has been changed since
* the last time step */
GlobalT::RelaxCodeT KBC_ControllerT::RelaxSystem(void)
{
	return GlobalT::kNoRelax;
}

/* output current configuration */
void KBC_ControllerT::WriteOutput(ostream& out) const
{
#pragma unused(out)
}

/**********************************************************************
* Protected
**********************************************************************/

/* read nodes from stream */
void KBC_ControllerT::ReadNodes(ifstreamT& in, iArrayT& id_list,
	iArrayT& nodes) const
{
	/* top level */
	const FEManagerT& fe_man = fNodeManager.FEManager();

	/* resolve format */
	IOBaseT::FileTypeT type = fe_man.InputFormat();
	switch (type)
	{
		case IOBaseT::kTahoe:
		{
			ifstreamT tmp;
			ifstreamT& in2 = fe_man.OpenExternal(in, tmp, cout, true,
				"KBC_ControllerT::ReadNodes: could not open file");

			int num_nodes;
			in2 >> num_nodes;
			nodes.Allocate(num_nodes);
			in2 >> nodes;
			break;
		}
		case IOBaseT::kTahoeII:
		{
			/* number of node sets */
			int num_sets;
			in >> num_sets;
			if (num_sets > 0)
			{
				/* open database */
				ModelFileT model_file;
				model_file.OpenRead(fe_man.ModelFile());

				/* echo set ID's */
				id_list.Allocate(num_sets);
				in >> id_list;
				
				/* collect */
				if (model_file.GetNodeSets(id_list, nodes) !=
				    ModelFileT::kOK) throw eBadInputValue;
			}
			break;
		}
		case IOBaseT::kExodusII:
		{
			/* number of node sets */
			int num_sets;
			in >> num_sets;
			if (num_sets > 0)
			{
				/* echo set ID's */
				id_list.Allocate(num_sets);
				in >> id_list;

				/* open database */
				ExodusT database(cout);
				database.OpenRead(fe_man.ModelFile());
				
				/* read collect all nodes in sets */
				database.ReadNodeSets(id_list, nodes);
			}
			break;
		}
		default:

			cout << "\n KBC_ControllerT::ReadNodes: unsupported input format: ";
			cout << type << endl;
			throw eGeneralFail;
	}

	/* internal numbering */
	nodes--;
}
