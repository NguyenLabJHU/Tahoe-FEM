// $Id: Quad2Tri.h,v 1.4 2002-10-08 20:51:51 paklein Exp $

// created: SAW 12/21/99

#ifndef _QUAD2TRI_H_
#define _QUAD2TRI_H_

#include "MakeCSE_ElementBaseT.h"

class NodeManagerPrimitive;
class MakeCSE_IOManager;

using namespace Tahoe;

class Quad2Tri : public MakeCSE_ElementBaseT
{
 public:

  enum Methods { kXMethod = 0,
		 kSlashMethod,
		 kBackSlashMethod,
		 kStarMethod };

  Quad2Tri (ostream& fMainOut, NodeManagerPrimitive& NMP, int method, int ID);

 protected:
  virtual void EchoConnectivity (MakeCSE_IOManager& theInput);
  virtual void EchoSideSets (MakeCSE_IOManager& theInput);

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
  int fMethod;
};

#endif
