#ifndef _ABAQUS_VARIABLES_T_H_
#define _ABAQUS_VARIABLES_T_H_

#include "StringT.h"

/** look up table for abaqus variables */
class AbaqusVariablesT
{
 public:
  enum PointT { kElementIntegration = 0,
		kElementSection,
		kElementWhole,
		kNodePoint };

  enum TypeT { kNotUsed = -1,
	       kNode = 0,
	       kElement,
	       kQuadrature };

  AbaqusVariablesT (void);
  void Set (const char* name, int key, bool nodenumberflag, AbaqusVariablesT::PointT point);
  void SetIOData (int dimension, AbaqusVariablesT::TypeT t, int index);

  const StringT& Name (void) const;
  int Key (void) const;
  bool NodeNumberFlag (void) const;
  PointT Point (void) const;
  int Dimension (void) const;
  TypeT Type (void) const;
  int IOIndex (void) const;

 private:
  StringT fName;
  int fKey;
  bool fNodeNumberFlag;
  PointT fPoint;

  int fDimension;
  TypeT fType;
  int fIOIndex;
};

inline AbaqusVariablesT::AbaqusVariablesT (void) :
  fName (""),
  fKey (-1),
  fNodeNumberFlag (false),
  fPoint (kNodePoint),
  fDimension (0),
  fType (kNotUsed),
  fIOIndex (-1)
{
}

inline void AbaqusVariablesT::Set (const char* name, int key, bool nodenumberflag, AbaqusVariablesT::PointT point)
{
  fName = name;
  fKey = key;
  fNodeNumberFlag = nodenumberflag;
  fPoint = point;
}

inline void AbaqusVariablesT::SetIOData (int dim, AbaqusVariablesT::TypeT t, int index)
{
  fDimension = dim;
  fType = t;
  fIOIndex = index;
}

inline const StringT& AbaqusVariablesT::Name (void) const { return fName; }
inline int AbaqusVariablesT::Key (void) const { return fKey; }
inline bool AbaqusVariablesT::NodeNumberFlag (void) const { return fNodeNumberFlag; }
inline AbaqusVariablesT::PointT AbaqusVariablesT::Point (void) const { return fPoint; }
inline int AbaqusVariablesT::Dimension (void) const { return fDimension; }
inline AbaqusVariablesT::TypeT AbaqusVariablesT::Type (void) const { return fType; }
inline int AbaqusVariablesT::IOIndex (void) const { return fIOIndex; }

#endif
