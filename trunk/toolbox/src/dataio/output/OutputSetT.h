/* $Id: OutputSetT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (03/07/2000)                                          */

#ifndef _OUTPUTSET_T_H_
#define _OUTPUTSET_T_H_

/* direct members */
#include "GeometryT.h"
#include "StringT.h"
#include "iArrayT.h"

/* forward declarations */
class iArray2DT;

class OutputSetT
{
public:

	/* constructor */
	OutputSetT(int ID, GeometryT::CodeT geometry_code,
		const iArray2DT& connectivities, const ArrayT<StringT>& n_labels,
		const ArrayT<StringT>& e_labels, bool changing);
	OutputSetT(const OutputSetT& source);

	/* print step counter */
	int PrintStep(void) const;
	void ResetPrintStep(void);
	void IncrementPrintStep(void);

	/* accessors */
	int ID(void) const;
	bool Changing(void) const;
	GeometryT::CodeT Geometry(void) const;
	const iArray2DT& Connectivities(void) const;
	const iArrayT& NodesUsed(void) const;
	const ArrayT<StringT>& NodeOutputLabels(void) const;
	const ArrayT<StringT>& ElementOutputLabels(void) const;

	/* dimensions */
	int NumNodes(void) const;
	int NumNodeValues(void) const;
	int NumElements(void) const;
	int NumElementValues(void) const;

private:

	/* not allowed */
	OutputSetT& operator=(OutputSetT&) { return *this; }

	/* set nodes used */
	void SetNodesUsed(void);

private:

	/* current print step */
	int fPrintStep;

	int  fID;
	bool fChanging;
	GeometryT::CodeT fGeometry;

	/* set nodes */
	const iArray2DT& fConnectivities;
	
	/* output labels */
	ArrayT<StringT> fNodeOutputLabels;
	ArrayT<StringT> fElementOutputLabels;
	
	/* derived */
	iArrayT fNodesUsed;
};

/* inlines */
inline int OutputSetT::PrintStep(void) const { return fPrintStep; }
inline void OutputSetT::ResetPrintStep(void) { fPrintStep = -1; }
inline void OutputSetT::IncrementPrintStep(void) { fPrintStep++; }

inline int OutputSetT::ID(void) const { return fID; }
inline bool OutputSetT::Changing(void) const { return fChanging; }
inline GeometryT::CodeT OutputSetT::Geometry(void) const { return fGeometry; }
inline const iArray2DT& OutputSetT::Connectivities(void) const
{
	return fConnectivities;
}

inline const iArrayT& OutputSetT::NodesUsed(void) const
{
	return fNodesUsed;
}

inline const ArrayT<StringT>& OutputSetT::NodeOutputLabels(void) const
{
	return fNodeOutputLabels;
}

inline const ArrayT<StringT>& OutputSetT::ElementOutputLabels(void) const
{
	return fElementOutputLabels;
}

/* dimensions */
inline int OutputSetT::NumNodeValues(void) const
{
	return fNodeOutputLabels.Length();
}

inline int OutputSetT::NumElementValues(void) const
{
	return fElementOutputLabels.Length();
}

#endif /* _OUTPUTSET_T_H_ */
