// file: Quad2Tri.h

// created: SAW 12/21/99

#ifndef _QUAD2TRI_H_
#define _QUAD2TRI_H_

#include "MakeCSE_ElementBaseT.h"
#include "CSEConstants.h"

namespace Tahoe {

class NodeManagerPrimitive;
class MakeCSEIOManager;
class ModelManagerT;
class dArray2DT;

class Quad2Tri : public MakeCSE_ElementBaseT
{
 public:

  Quad2Tri (ostream& fMainOut, NodeManagerPrimitive& NMP, CSEConstants::SplitMethodT method, const StringT& ID);

 protected:
  virtual void EchoConnectivity (ModelManagerT& theInput);
  virtual void EchoSideSets (ModelManagerT& model, MakeCSE_IOManager& theInput);

 private:

  void Translate (void);
  void Allocate (int numQuadNodes);
  int  ElementCentroid (int* quad, int numQuadNodes, const dArray2DT& coords);

  void XMethodNumbering (int& count, int newnode, int* quad);
  void SlashNumbering (int& count, int* quad);
  void BackSlashNumbering (int& count, int* quad);
  void StarNumbering (int& count, int newnode, int* quad);

  void XMethodSideSets (ArrayT<iArray2DT>& sidesets);
  void SlashSideSets (ArrayT<iArray2DT>& sidesets);
  void BackSlashSideSets (ArrayT<iArray2DT>& sidesets);
  void StarSideSets (ArrayT<iArray2DT>& sidesets);

 private:
  iArray2DT fConn;
  NodeManagerPrimitive* theNodes;
  CSEConstants::SplitMethodT fMethod;
};
}
#endif
