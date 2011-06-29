/* created saw 13 oct 2006           */

#ifndef _ABAQUSINP_T_H_
#define _ABAQUSINP_T_H_

/* direct members */
#include "ios_fwd_decl.h"
#include "StringT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "GeometryT.h"
#include "sArrayT.h"

namespace Tahoe {

/* forward declarations */
class iAutoArrayT;
class iArray2DT;

class AbaqusINPT
{
public:
	AbaqusINPT (void);
	
	bool OpenRead (const StringT& filename);
	
	int NumNodes (void) const;
	int NumElSets (void) const;
	int NumNodeSets (void) const;
	int NumSurfaces (void) const;

	void NodeIDs (iArrayT& ids) const;
	void Coordinates (dArray2DT& coords) const;
	
	void SetNames (ArrayT<StringT>& name, const char* setname, const char* settype) const;
	int NodeSetLength (const StringT& name) const;
	int ElSetLength (const StringT& name) const;
	void NodeSet (const StringT& name, iArrayT& set) const;
	void ElSet (const StringT& name, iArrayT& set) const;

	int NumElementNodesforSet (const StringT& name) const;
	bool Connectivity (const StringT& name, iArray2DT& conn) const;
	GeometryT::CodeT GeometryCode (const StringT& name);
	
	/* the *surface keyword would be where side set data would come from */
	//void SurfaceNames (ArrayT<StringT>& names) const;
	//void Surface (const StringT& name, iArray2DT& set) const;
	
private:
	void ScanFile (void);
	void GatherSimpleData (void);
	void ReadSetData (ifstream& in, iAutoArrayT& ids, sArrayT& names, StringT& line);
	void ReadGeneratedSetData (ifstreamT& in, iAutoArrayT& ids, StringT& line);
	bool NextKeyWord (ifstream& in, StringT& keyword) const;
	
	bool ExtractOptionName (const StringT& keyword, const char* settype, StringT& name) const;
	
	GeometryT::CodeT TranslateElementName (const char* name, int& numelenodes) const;
	GeometryT::CodeT TranslateContinuum (const char* name, int& numelnodes) const;
	GeometryT::CodeT TranslateSpring (const char* name, int& numelnodes) const;
	GeometryT::CodeT Translate2D (const char* name, int& numelnodes) const;
	GeometryT::CodeT Translate3D (const char* name, int& numelnodes) const;
	GeometryT::CodeT TranslateShell (const char* name, int& numelnodes) const;
	GeometryT::CodeT TranslateRigid (const char* name, int& numelnodes) const;
	
private:
	StringT fInputFile;
	int fNumNodes;
	int fNumElSets;
	int fNumNodeSets;
	int fNumSurfaces;
	
	sArrayT fElementSets;
	sArrayT fNodeSets;
	
	iArrayT fNodeIDs;
	ArrayT<iArrayT> fElementIDs;
	ArrayT<sArrayT> fSetsinElSets;
	iArrayT fElementSetNodes;
	ArrayT<iArrayT> fNodeSetIDs;
	ArrayT<sArrayT> fSetsinNSets;
};

inline int AbaqusINPT::NumNodes (void) const { return fNumNodes; }
inline int AbaqusINPT::NumElSets (void) const {	return fNumElSets; }
inline int AbaqusINPT::NumNodeSets (void) const { return fNumNodeSets; }
inline int AbaqusINPT::NumSurfaces (void) const { return fNumSurfaces; }

} //namespace Tahoe
#endif
