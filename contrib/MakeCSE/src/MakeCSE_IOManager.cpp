/* $Id: MakeCSE_IOManager.cpp,v 1.3 2002-10-28 21:36:33 sawimme Exp $ */
#include "MakeCSE_IOManager.h"
#include "ExceptionT.h"
#include "ifstreamT.h"

using namespace Tahoe;

MakeCSE_IOManager::MakeCSE_IOManager (void)
{
}

CSEConstants::CSEMethodT MakeCSE_IOManager::int2CSEMethodT (int i)
{
  switch (i)
    {
    case CSEConstants::kFacet: return CSEConstants::kFacet;
    case CSEConstants::kZone: return CSEConstants::kZone;
    case CSEConstants::kBoundary: return CSEConstants::kBoundary;
    default:
      {
	cout << "\n Invalid Method \n";
	throw ExceptionT::kGeneralFail;
      }
    }
  return CSEConstants::kFacet;
}

CSEConstants::ZoneEdgeT MakeCSE_IOManager::int2ZoneEdgeT (int i)
{
  switch (i)
    {
    case CSEConstants::kSingleZE: return CSEConstants::kSingleZE;
    case CSEConstants::kDoubleZE: return CSEConstants::kDoubleZE;
    case CSEConstants::kMixSingZE: return CSEConstants::kMixSingZE;
    case CSEConstants::kMixDoubZE: return CSEConstants::kMixDoubZE;
    default:
      {
	cout << "\n Invalid Zone Edging \n";
	throw ExceptionT::kGeneralFail;
      }
    } 
  return CSEConstants::kSingleZE;
}

CSEConstants::NodeMapMethodT MakeCSE_IOManager::int2NodeMapMethodT (int i)
{
  switch (i)
    {
    case CSEConstants::kSurface1: return CSEConstants::kSurface1;
    case CSEConstants::kSurface2: return CSEConstants::kSurface2;
    case CSEConstants::kMap:      return CSEConstants::kMap;
    case CSEConstants::kSplit:    return CSEConstants::kSplit;
    default:        return CSEConstants::kMap;
    }
}

CSEConstants::RenumberMethodT MakeCSE_IOManager::int2RenumberMethodT (int i)
{
  switch (i)
    {
    case CSEConstants::kNoRenumber:    return	 CSEConstants::kNoRenumber;
    case CSEConstants::kRenumberAdded: return	 CSEConstants::kRenumberAdded;
    case CSEConstants::kRenumberAll:   return	 CSEConstants::kRenumberAll;
    default:
      {
	cout << "\n Invalid Renumbering Option " << i << "\n";
	throw ExceptionT::kGeneralFail;
      }
    }
  return CSEConstants::kNoRenumber;
}

CSEConstants::SplitMethodT MakeCSE_IOManager::int2SplitMethodT (int i)
{
  switch (i)
    {
    case CSEConstants::kXMethod:         return CSEConstants::kXMethod;
    case CSEConstants::kSlashMethod:     return CSEConstants::kSlashMethod;
    case CSEConstants::kBackSlashMethod: return CSEConstants::kBackSlashMethod;
    case CSEConstants::kStarMethod:      return CSEConstants::kStarMethod;
    default:               return CSEConstants::kXMethod;
    }
}
