/* $Id: KBC_ControllerT.cpp,v 1.1.1.1.6.2 2001-10-10 20:09:19 sawimme Exp $ */
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
	ModelManagerT* model = fe_man.ModelManager();

	/* resolve format */
	if (fe_man.InputFormat() == IOBaseT::kTahoe)
	  {
	    ifstreamT tmp;
	    ifstreamT& in2 = model->OpenExternal(in, tmp, cout, true,
		  "KBC_ControllerT::ReadNodes: could not open file");
	    
	    int num_nodes;
	    in2 >> num_nodes;
	    nodes.Allocate(num_nodes);
	    in2 >> nodes;
	  }
	else
	  {
	    /* number of node sets */
	    int num_sets;
	    in >> num_sets;
	    if (num_sets > 0)
	      {
		/* echo set ID's */
		id_list.Allocate(num_sets); // future: change to string
		ArrayT<StringT> id_name (num_sets);

		/* collect sets */
		iAutoArrayT temp;
		for (int i=0; i < num_sets; i++)
		  {
		    in >> id_name[i];
		    int index = model->NodeSetIndex (id_name[i]);
		    id_list[i] = index+1;
		    iArrayT tn = model->NodeSet (index);
		    temp.AppendUnique (tn);
		  }

		nodes.Allocate (temp.Length());
		nodes.CopyPart (0, temp, 0, temp.Length());
		nodes.SortAscending ();
	      }
	  }
}
