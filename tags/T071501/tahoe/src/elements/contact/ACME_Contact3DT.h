/* $Id: ACME_Contact3DT.h,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (10/15/2000)                                          */

#ifndef _ACME_PENALTY_CONTACT3D_T_H_
#define _ACME_PENALTY_CONTACT3D_T_H_

/* base class */
#include "PenaltyContact3DT.h"

/* library support options */
#ifdef __ACME__

/* direct members */
#include "ContactSearch.h"

class ACME_Contact3DT: public PenaltyContact3DT
{
public:

	/* constructor */
	ACME_Contact3DT(FEManagerT& fe_manager);

	/* destructor */
	~ACME_Contact3DT(void);

	/* initialization after constructor */
	virtual void Initialize(void);

	/* restart functions */
	virtual void ReadRestart(istream& in);
	virtual void WriteRestart(ostream& out) const;

protected:

	/* initialization steps */
	virtual void EchoConnectivityData(ifstreamT& in, ostream& out);
	virtual void SetWorkSpace(void); /* constuct ContactSearch object */

	/* steps in setting contact configuration */
	virtual bool SetActiveInteractions(void); // "internal" data

private:

	/* translate native facet data to ACME format */
	void GenerateACMEConnectivities(iArrayT& connectivities);
	 	
protected:

	/* ACME search object */
	ContactSearch* fContactSearch;

	/* parameters */
	ArrayT<ContactSearch::ContactNode_Type> Node_Block_Types;
	iArrayT Number_Nodes_in_Blocks;
	iArrayT Node_Global_IDs;
	ArrayT<ContactSearch::ContactFace_Type> Face_Block_Types;
	iArrayT Number_Faces_in_Blocks;
	iArrayT Connectivity;
	iArrayT Nodal_Comm_Proc_IDs;
	iArrayT Number_Nodes_to_Partner;
	iArrayT Communication_Nodes;
	int mpi_communicator;

	/* contact search data */
	dArrayT fSearchData;
	
	/* interactions */
	AutoArrayT<int> fnode_block_ids;
	//AutoArrayT<int> fnode_indices_in_block;
	//AutoArrayT<int> fface_block_ids;
	//AutoArrayT<int> fface_indices_in_block;
	AutoArrayT<int> fface_proc;
	AutoArrayT<double> fdata;
};

#endif /* __ACME__ */
#endif /* _ACME_PENALTY_CONTACT3D_T_H_ */
