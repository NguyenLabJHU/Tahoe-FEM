/* created sawimme (05/17/2001)  */

#include "PatranT.h"
#include "ifstreamT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "ExceptionCodes.h"
#include "iAutoArrayT.h"

PatranT::PatranT (ostream &message_out) :
  fOut (message_out)
{
}

PatranT::~PatranT (void) 
{
}

bool PatranT::OpenRead (const StringT& filename)
{
  ifstreamT tmp (filename);
  if (!tmp.is_open()) 
    {
      fOut << "PatranT::OpenRead, unable to open file: "
	   << filename << "\n";
      return false;
    }

  file_name = filename;
  return true;
}

bool PatranT::OpenWrite (const StringT& filename)
{
  fOut << "PatranT::OpenWrite, not programed for writing\n";
  return false;
}

int PatranT::NumNodes (void) const
{
  int ID, IV, KC, num_nodes;
  ifstream in (file_name);
  if (!AdvanceTo (in, kSummary, ID, IV, KC)) 
    {
      fOut << "PatranT::NumNodes, no nodes found\n";
      return -1;
    }
  in >> num_nodes;
  return num_nodes;
}

int PatranT::NumElements (void) const
{
  int ID, IV, KC, num_nodes, num_elements;
  ifstream in (file_name);
  if (!AdvanceTo (in, kSummary, ID, IV, KC))  
    {
      fOut << "PatranT::NumElements, no elements found\n";
      return -1;
    }
  in >> num_nodes >> num_elements;
  return num_elements;
}

int PatranT::NumNamedComponents (void) const
{
  int ID, IV, KC, num = 0;
  ifstream in (file_name);
  while (AdvanceTo (in, kNamedComponents, ID, IV, KC)) 
    {
      num++;
      ClearPackets (in, KC+1);
    }
  return num;
}

int PatranT::NumDimensions (void) const
{
  int ID, IV, KC, num_dims;
  ifstream in (file_name);
  if (!AdvanceTo (in, kElement, ID, IV, KC))   
    {
      fOut << "PatranT::NumDimensions, no elements found\n";
      return -1;
    }
  switch (IV)
    {
    case kBarShape:
    case kTriShape:
    case kQuadShape:
      return 2;
      break;
    case kTetShape:
    case kWedgeShape:
    case kHexShape:
      return 3;
      break;
    }
  fOut << "\n PatranT::NumDimensions: Unknown element shape ID=" << ID
       << " shape= " << IV << "\n";
  return -1;
}

bool PatranT::NamedComponents (ArrayT<StringT>& names) const
{
  int ID, IV, KC, num = 0;
  ifstream in (file_name);
  while (AdvanceTo (in, kNamedComponents, ID, IV, KC)) 
    {
      if (num >= names.Length()) 
	{
	  fOut << "PatranT::NamedComponents, incorrect allocation\n";
	  return false;
	}
      ClearPackets (in, 1);
      in >> names[num++];
      ClearPackets (in, KC);
    }
  if (num != names.Length()) 
    {
      fOut << "PatranT::NamedComponents, incorrect amount\n";
      return false;
    }
  return true;
  
}

bool PatranT::ReadGlobalNodeMap (iArrayT& map) const
{
  int numnodes = NumNodes ();
  int ID, IV, KC, num = 0;
  map.Allocate (numnodes);
  ifstream in (file_name);
  while (AdvanceTo (in, kNode, ID, IV, KC))
    {
      if (num >= map.Length())
	{
	  fOut << "PatranT::ReadGlobalNodeMap, incorrect allocation\n";
	  return false;
	}
      map [num++] = ID;
      ClearPackets (in, KC+1);
    } 
  if (num != numnodes)  
    {
      fOut << "PatranT::ReadGlobalNodeMap num != numnodes\n"
	   << num << " " << numnodes << "\n";
      return false;
    }
  return true;
}

bool PatranT::ReadGlobalElementMap (iArrayT& map) const
{
  int numelems = NumElements ();
  int ID, IV, KC, num = 0;
  ifstream in (file_name);
  map.Allocate (numelems);
  while (AdvanceTo (in, kElement, ID, IV, KC))
    {
      if (num >= map.Length())
	{
	  fOut << "PatranT::ReadGlobalElementMap, incorrect allocation\n";
	  return false;
	}
      map [num++] = ID;
      ClearPackets (in, KC+1);
    } 
  if (num != numelems) 
    {
      fOut << "PatranT::ReadGlobalElementMap num != numelems\n"
	   << num << " " << numelems << "\n";
      return false;
    }
  return true;
}

bool PatranT::NumNodesInSet (StringT& title, int& num) const
{
  num = -1;
  iArrayT list;
  if (!ReadNamedComponent (title, list)) 
    {
      fOut << "PatranT::NumNodesInSet, unable to read named components\n";
      return false;
    }

  /* pull node IDs from component list */
  num = 0;
  int *it = list.Pointer();
  int *il = list.Pointer() + 1;
  for (int i=0; i < list.Length(); i++, it += 2, il += 2)
    if (*it == kNodeType)
      num++;

  return true;
}

bool PatranT::ReadCoordinates (dArray2DT& coords, int dof) const
{
  if (dof != coords.MinorDim())
    {
      fOut << "PatranT::ReadCoordinates, incorrect minor allocation\n";
      return false;
    }
  int ID, IV, KC;
  ifstream in (file_name);
  int count = 0;
  dArrayT temp (dof);
  while (count < coords.MajorDim())
    {
      if (!AdvanceTo (in, kNode, ID, IV, KC)) 
	{
	  fOut << "PatranT::ReadCoordinates, unable to find node, count = "
	       << count << "\n";
	  return false;
	}
      ClearPackets (in, 1);
      in >> temp;
      coords.SetRow (count, temp);
      ClearPackets (in, KC);
      count ++;
    }
  return true;
}

bool PatranT::ReadElementBlockDims (const StringT& title, int& num_elems, int& num_elem_nodes) const
{
  int ID, IV, KC, num_nodes, code;
  iArrayT elems;
  if (!ReadElementSet (title, code, elems))
    {
      fOut << "PatranT::ReadElementBlockDims, unable to read set\n";
      return false;
    }
  num_elems = elems.Length();
  ifstream in (file_name);
  while (AdvanceTo (in, kElement, ID, IV, KC))
    {
      if (elems.HasValue (ID))
	{
	  ClearPackets (in, 1);
	  in >> num_elem_nodes;
	  return true;
	}
      else
	ClearPackets (in, KC + 1);
    }
  fOut << "\n PatranT::ReadElementBlockDims, unable to find element for set "
       << title << "\n";
  return false;
}

bool PatranT::ReadConnectivity (const StringT& title, int& namedtype, iArray2DT& connects) const
{
  iArrayT elems;
  if (!ReadElementSet (title, namedtype, elems))
    {
      fOut << "PatranT::ReadConnectivity, unable to read set";
      return false;
    }
  
  int count = 0, num_nodes;
  int *pc = connects.Pointer();
  int ID, IV, KC;
  ifstream in (file_name);
  iArrayT temp (connects.MinorDim());
  while (count < connects.MajorDim())
    {
      if (!AdvanceTo (in, kElement, ID, IV, KC)) 
	{
	  fOut << "PatranT::ReadConnectivity, unable to find element";
	  return false;
	}
      if (elems.HasValue (ID))
	{
	  ClearPackets (in, 1);
	  KC--;

	  in >> num_nodes;
	  if (num_nodes != connects.MinorDim()) 
	    {
	      fOut << "\n PatranT::ReadConnectivity: Warning num_nodes_per_element mismatch for "
		   << title << " ID=" << ID << " has " << num_nodes 
		   << ", which doesn't match " << connects.MinorDim() << "\n";
	    }
	  ClearPackets (in, 1);
	  KC--;

	  in >> temp;
	  connects.SetRow (count, temp);
	  KC -= (num_nodes + 9)/10;
	  ClearPackets (in, KC);

	  count ++;
	}
      else
	ClearPackets (in, KC + 1);
    }
  return true;
}

bool PatranT::ReadAllElements (ArrayT<iArrayT>& connects, iArrayT& elementtypes) const
{
  int ID, IV, KC, num_nodes, count = 0;
  ifstream in (file_name);
  while (AdvanceTo (in, kElement, ID, IV, KC)) 
    {
      if (count >= connects.Length() ||
	  count >= elementtypes.Length())
	{
	  fOut << "\nPatranT::ReadAllElements, incorrect allocation\n";
	  throw eSizeMismatch;
	}

      ClearPackets (in, 1);
      KC--;

      in >> num_nodes;
      connects[count].Allocate (num_nodes);
      ClearPackets (in, 1);
      KC--;

      in >> connects[count];
      KC -= (num_nodes + 9)/10;
      ClearPackets (in, KC);

      elementtypes[count] = IV;
      count ++;
    }
  return true;
}

bool PatranT::ReadElementSet (const StringT& title, int& namedtype, iArrayT& elems) const
{
  iArrayT list;
  if (!ReadNamedComponent (title, list)) 
    {
      fOut << "PatranT::ReadElementSet, unable to read named component\n";
      return false;
    }

  /* pull element IDs from component list */
  iAutoArrayT set;
  namedtype = -1;
  int *it = list.Pointer();
  int *il = list.Pointer() + 1;
  int geocode = -1;
  for (int i=0; i < list.Length()/2; i++, it += 2, il += 2)
    if ((*it >   5 && *it < 19) ||
	(*it > 105 && *it < 119) ||
	(*it > 205 && *it < 219))
      {
	set.Append (*il);
	if (namedtype == -1)
	  namedtype = *it;
      }

  elems.Allocate (set.Length());
  elems.CopyPart (0, set, 0, set.Length());
  return true;
}

bool PatranT::ReadDistLoadSetDims (int setID, int& num_elems) const
{
  int ID, IV, KC;
  ifstream in (file_name);
  num_elems = 0;
  while (AdvanceTo (in, kDistLoads, ID, IV, KC))
    {
      if (setID == IV) num_elems++;
      ClearPackets (in, KC + 1);
    }
  return true;
}

bool PatranT::ReadDistLoadSet (int setID, iArray2DT& facets) const
{
  int ID, IV, KC;
  ifstream in (file_name);
  int num_elems = 0;
  while (num_elems < facets.MajorDim())
    {
      if (!AdvanceTo (in, kDistLoads, ID, IV, KC)) 
	{
	  fOut << "PatranT::ReadDistLoadSet, unable to find dist load\n";
	  return false;
	}

      if (setID == IV) 
	{
	  int itemp;
	  StringT temp;
	  ClearPackets (in, 1);
	  in >> temp >> itemp;
	  facets (num_elems, 0) = ID;
	  facets (num_elems, 1) = itemp;
	  num_elems++;
	  ClearPackets (in, KC);
	}
      else
	ClearPackets (in, KC + 1);
    }
  return true;
}

bool PatranT::ReadNodeSet (const StringT& title, iArrayT& nodes) const
{
  iArrayT list;
  if (!ReadNamedComponent (title, list)) 
    {
      fOut << "PatranT::ReadNodeSet, unable to read named components\n";
      return false;
    }

  /* pull node IDs from component list */
  iAutoArrayT set;
  int *it = list.Pointer();
  int *il = list.Pointer() + 1;
  for (int i=0; i < list.Length()/2; i++, it += 2, il += 2)
    if (*it == kNodeType)
      set.Append (*il);
    
  nodes.Allocate (set.Length());
  nodes.CopyPart (0, set, 0, set.Length());
  return true;
}

bool PatranT::ReadNodeSets (const ArrayT<StringT>& title, iArrayT& nodes) const
{
  iArrayT count (title.Length());
  for (int i=0; i < title.Length(); i++)
    if (!NumNodesInSet (title[i], count[i])) 
      {
	fOut << "PatranT::ReadNodeSets, unable to read num nodes in set\n";
	return false;
      }

  nodes.Allocate (count.Sum());
  iArrayT nt;
  for (int j=0, k=0; j < title.Length(); j++)
    {
      nt.Allocate (count[j]);
      if (!ReadNodeSet (title[j], nt)) 
	{
	  fOut << "PatranT::ReadNodeSets, unable to read node set "
	       << title[j] << "\n";
	  return false;
	}
      nodes.CopyPart (k, nt, 0, count[j]);
      k += count[j];
    }
}

/**************************************************************************
* Private
**************************************************************************/

bool PatranT::ReadNamedComponent (const StringT &title, iArrayT& list) const
{
  /* read list of elements */
  int ID, IV, KC;
  ifstream in (file_name);
  bool foundtitle = false;
  StringT name;
  while (!foundtitle)
    {
      if (!AdvanceTo (in, kNamedComponents, ID, IV, KC)) return false;
      ClearPackets (in, 1);
      in >> name;
      if (strncmp (name.Pointer(), title.Pointer(), title.Length()) == 0) 
	foundtitle = true;
      else 
	ClearPackets (in, KC);
    }
  list.Allocate (IV);
  in >> list;

  return true;
}

bool PatranT::AdvanceTo (ifstream &in, int target, int& ID, int &IV, int &KC) const
{
  int IT = 0;
  while (in >> IT)
    {
      in >> ID >> IV >> KC;
      
      if (IT == target) return true;

      /* skip to next entry */
      ClearPackets (in, KC + 1);
    }

  //fOut << "PatranT::AdvanceTo: Cannot find: " << target << '\n';
  return false;
}

void PatranT::ClearPackets (ifstream &in, int KC) const
{
  char line[255];
  for (int i=0; i < KC; i++)
    in.getline (line, 254);
}
