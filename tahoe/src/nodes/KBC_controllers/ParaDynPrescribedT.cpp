/* $Id: ParaDynPrescribedT.cpp,v 1.1.2.2 2003-09-18 22:38:29 cjkimme Exp $ */
#include "ParaDynPrescribedT.h"
#include "NodeManagerT.h"
#include "FEManagerT.h"
#include "ModelManagerT.h"
#include "ifstreamT.h"
#include "RandomNumberT.h"
#include "MessageT.h"
#include "CommunicatorT.h"
#include "iArrayT.h"
#include "dArray2DT.h"

using namespace Tahoe;

dArray2DT* ParaDynPrescribedT::coord_ptr = NULL;
dArray2DT* ParaDynPrescribedT::vel_ptr = NULL;

const double fkB = 0.00008617385;



/* constructor */
ParaDynPrescribedT::ParaDynPrescribedT(NodeManagerT& node_manager, BasicFieldT& field):
	KBC_ControllerT(node_manager),
	fField(field),
	fNodes(),
	fDummySchedule(1.0),
	n_index(0)
{
	// nein
}

ParaDynPrescribedT::~ParaDynPrescribedT(void)
{
  
}

void ParaDynPrescribedT::WriteParameters(ostream& out) const
{
	/* inherited */
	KBC_ControllerT::WriteParameters(out);

}

void ParaDynPrescribedT::Initialize(ifstreamT& in)
{

	if (fNodeManager.NumSD() != 3)
		ExceptionT::GeneralFail("ParaDynPrescribedT::Initialize","Must be in 3D\n");

	in >> disp_file_root;
	in >> disp_file_suffix;
	in >> vel_file_root;
	in >> vel_file_suffix;
	in >> n_start;
	in >> n_end;
	
	coords.Dimension(fNodeManager.NumNodes(),3);
	vels.Dimension(fNodeManager.NumNodes(),3);
	
	coords = fNodeManager.InitialCoordinates();
	vels = 0.;
	
	coord_ptr = &coords;
	vel_ptr = &vels;
	
	fNodes.Dimension(fNodeManager.NumNodes());
	fNodes.SetValueToPosition();
	
	SetBCCards();

}

void ParaDynPrescribedT::InitStep(void)
{

	if (n_index < n_end)
		n_index++;
	else
		cout << "ParaDynPrescribedT::InitStep : index out of range but carrying on anyway\n";

	StringT disp_file_name, vel_file_name;
	
	MakeFileName(disp_file_name, disp_file_root, disp_file_suffix, n_index);
	MakeFileName(vel_file_name, vel_file_root, vel_file_suffix, n_index);
	
	ReadArray(disp_file_name,coords,8);
	
	ReadArray(vel_file_name,vels,6);
	
	SetBCCards();
					
	/* inherited */
	KBC_ControllerT::InitStep();

}


void ParaDynPrescribedT::InitialCondition(void)
{
#if 0
	/* number of scaled nodes */
	int n_scaled = 0;
	int ndof = fField.NumDOF();
		
	/* get MPI stuff */
	CommunicatorT& communicator = fNodeManager.FEManager().Communicator();
	int nProcs = fNodeManager.Size();
	int thisProc = fNodeManager.Rank();
	const ArrayT<int>* pMap = fNodeManager.ProcessorMap();
	
	/* only prescribe nodes on this processor */
	iArrayT myNodes; 
	
	if (fNodeManager.Size() == 1 || !pMap)
	{
		myNodes.Set(fNodes.Length(), fNodes.Pointer());
		n_scaled = fNodes.Length();
	}
	else
	{
		for (int i = 0; i < myNodes.Length(); i++)
			if ((*pMap)[fNodes[i]] == thisProc)
			{
				myNodes[i] = fNodes[i];
				n_scaled++;
			}
			else
				myNodes[i] = -1;
	}
	
	/* generate BC cards */
	fKBC_Cards.Dimension(2*n_scaled*ndof);
	KBC_CardT* pcard_x = fKBC_Cards.Pointer();
	KBC_CardT* pcard_v = pcard_x + n_scaled*ndof;
	for (int i = 0; i < myNodes.Length(); i++)
	{
		if (myNodes[i] >= 0)
		{
			double* x_i = coords(myNodes[i]);
			double* v_i = vels(myNodes[i]);	
			
	    	for (int j = 0; j < ndof; j++)
			{	
			
				pcard_x->SetValues(myNodes[i], j, KBC_CardT::kDsp, 0, *x_i++);
				pcard_v->SetValues(myNodes[i], j, KBC_CardT::kVel, 0, *v_i++);
		

				/* dummy schedule */
				pcard_x->SetSchedule(&fDummySchedule);
				pcard_v->SetSchedule(&fDummySchedule);
				pcard_x++;
				pcard_v++;
			}
		} 
	}
#endif
}

/**********************************************************************
 * Protected
 **********************************************************************/

/* initialize the current step */
void ParaDynPrescribedT::SetBCCards(void)
{
	/* number of scaled nodes */
	int n_scaled = 0; 
	int ndof = fField.NumDOF();
	
	/* get MPI stuff */
	CommunicatorT& communicator = fNodeManager.FEManager().Communicator();
	int nProcs = fNodeManager.Size();
	int thisProc = fNodeManager.Rank();
	const ArrayT<int>* pMap = fNodeManager.ProcessorMap();	

	/* figure out which nodes to affect */
	iArrayT myNodes;

	if (fNodeManager.Size() == 1 || !pMap)
	{
		myNodes.Set(fNodes.Length(), fNodes.Pointer());
		n_scaled = fNodes.Length();
	}
	else
	{
		for (int i = 0; i < myNodes.Length(); i++)
			if ((*pMap)[fNodes[i]] == thisProc)
			{
				myNodes[i] = fNodes[i];
				n_scaled++;
			}
			else
				myNodes[i] = -1;
	}

	/* grab the positions */
	const dArray2DT& initCoords = fNodeManager.InitialCoordinates();

	/* generate BC cards */
	fKBC_Cards.Dimension(/*2*/n_scaled*ndof);
	KBC_CardT* pcard_x = fKBC_Cards.Pointer();
//	KBC_CardT* pcard_v = pcard_x + n_scaled*ndof;
	for (int i = 0; i < myNodes.Length(); i++)
	{
		int index = myNodes[i];
		if (index >= 0)
		{
			double* x_i = coords(index);
			double* X_i = initCoords(index);
//			double* v_i = vels(index);	
			
	    	for (int j = 0; j < ndof; j++)
			{	
			
				pcard_x->SetValues(index, j, KBC_CardT::kDsp, 0, *x_i++ - *X_i++);
//				pcard_v->SetValues(index, j, KBC_CardT::kVel, 0, *v_i++);

				/* dummy schedule */
				pcard_x->SetSchedule(&fDummySchedule);
//				pcard_v->SetSchedule(&fDummySchedule);
				pcard_x++;
//				pcard_v++;
			}
		} 
	}
	
}


void ParaDynPrescribedT::ReadArray(StringT& file_name, dArray2DT& info, int numGarbageLines)
{

	ifstreamT file_in(file_name);

	if (!file_in.is_open())
		ExceptionT::BadInputValue("ParaDynPrescribedT::ReadArray", "error opening file: %s", file_name.Pointer());

	StringT garbage;
	garbage.GetLineFromStream(file_in);
	
	int natoms;
	
	file_in >> natoms;
	if (natoms != fNodes.Length())
		ExceptionT::GeneralFail("ParaDynPrescribedT::ReadArray","Mismatch in number of atoms\n");

	for (int i = 0; i < numGarbageLines; i++)
		garbage.GetLineFromStream(file_in);
	
	for (int i = 0; i < natoms; i++)
	{
		int index;
		file_in >> index;
		
		int igarbage;
		file_in >> igarbage;
		
		index--;
		double* p = info(index);
		for (int j = 0; j < 3; j++)
			file_in >> *p++;
	}

}

void ParaDynPrescribedT::MakeFileName(StringT& fileName, StringT& rootName, StringT& suffixName, int index)
{
	fileName = rootName;
	fileName.Append(".");
	fileName.Append(index,3);
	fileName.Append(".");
	fileName.Append(suffixName);
}

/********************************************************************
 *   static access functions                                        *
 ********************************************************************/
 
dArray2DT* ParaDynPrescribedT::SendCoordinates(void)
{
	return coord_ptr;
}

dArray2DT* ParaDynPrescribedT::SendVelocities(void)
{
	return vel_ptr;
}


