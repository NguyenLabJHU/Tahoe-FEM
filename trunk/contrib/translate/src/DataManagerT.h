/* created: sawimme August 2001 */

#ifndef _DATA_MANAGER_T_H_
#define _DATA_MANAGER_T_H_

#include "ModelManagerT.h"

class DataManagerT : public ModelManagerT
{
 public:
  DataManagerT (ostream& message);

  int NumTimeSteps (void);
  void TimeSteps (dArrayT& steps);

  int NumNodeVariables (void);
  int NumElementVariables (void);
  int NumQuadratureVariables (void);

  void NodeLabels (ArrayT<StringT>& labels);
  void ElementLabels (ArrayT<StringT>& labels);
  void QuadratureLabels (ArrayT<StringT>& labels);

  void AllNodeVariables (int stepindex, dArray2DT& values);
  void NodeVariables (int stepindex, StringT& elsetname, dArray2DT& values);
  void NodeSetVariables (int stepindex, StringT& nsetname, dArray2DT& values);

  void AllElementVariables (int stepindex, dArray2DT& values);
  void ElementVariables (int stepindex, StringT& elsetname, dArray2DT& values);

  void AllQuadratureVariables (int stepindex, dArray2DT& values);
  void QuadratureVariables (int stepindex, StringT& elsetname, dArray2DT& values);
};

inline DataManagerT::DataManagerT (ostream& message) :
  ModelManagerT (message) {}

inline int DataManagerT::NumTimeSteps (void)
{ return fInput->NumTimeSteps(); }

inline void DataManagerT::TimeSteps (dArrayT& steps)
{ fInput->ReadTimeSteps (steps); }

inline int DataManagerT::NumNodeVariables (void)
{ return fInput->NumNodeVariables (); }

inline int DataManagerT::NumElementVariables (void)
{ return fInput->NumElementVariables (); }

inline int DataManagerT::NumQuadratureVariables (void)
{ return fInput->NumQuadratureVariables (); }

inline void DataManagerT::NodeLabels (ArrayT<StringT>& labels)
{ fInput->ReadNodeLabels (labels); }

inline void DataManagerT::ElementLabels (ArrayT<StringT>& labels)
{ fInput->ReadElementLabels (labels); }

inline void DataManagerT::QuadratureLabels (ArrayT<StringT>& labels)
{ fInput->ReadQuadratureLabels (labels); }

inline void DataManagerT::AllNodeVariables (int stepindex, dArray2DT& values)
{ fInput->ReadAllNodeVariables (stepindex, values); }

inline void DataManagerT::NodeVariables (int stepindex, StringT& elsetname, dArray2DT& values)
{ fInput->ReadNodeVariables (stepindex, elsetname, values); }

inline void DataManagerT::NodeSetVariables (int stepindex, StringT& nsetname, dArray2DT& values)
{ fInput->ReadNodeSetVariables (stepindex, nsetname, values); }

inline void DataManagerT::AllElementVariables (int stepindex, dArray2DT& values)
{ fInput->ReadAllElementVariables (stepindex, values); }

inline void DataManagerT::ElementVariables (int stepindex, StringT& elsetname, dArray2DT& values)
{ fInput->ReadElementVariables (stepindex, elsetname, values); }

inline void DataManagerT::AllQuadratureVariables (int stepindex, dArray2DT& values)
{ fInput->ReadAllQuadratureVariables (stepindex, values); }

inline void DataManagerT::QuadratureVariables (int stepindex, StringT& elsetname, dArray2DT& values)
{ fInput->ReadQuadratureVariables (stepindex, elsetname, values); }

#endif
