#ifndef _CSE_CONSTANTS_H_
#define _CSE_CONSTANTS_H_

namespace Tahoe {

class CSEConstants
{
 public:
  enum SettingsT { kNotSet = -1 };
  
  enum EdgeTypeT { kNoNeighbor = kNotSet, 
		   kExteriorFacet = -5 };
  
  enum CSEMethodT { kFacet = 1, 
		    kZone, 
		    kBoundary };
  
  enum ZoneEdgeT { kSingleZE = 1,
		   kDoubleZE = 2,
		   kMixSingZE = 3,
		   kMixDoubZE = 4 };
  
  enum NodeMapMethodT { kSurface1 = 0,
			kSurface2,
			kMap,
			kSplit };
  
  enum RenumberMethodT { kNoRenumber = 0,
			 kRenumberAdded = 1,
			 kRenumberAll = 2 };

  enum SplitMethodT { kXMethod = 0,
		      kSlashMethod,
		      kBackSlashMethod,
		      kStarMethod };
};

}
#endif
