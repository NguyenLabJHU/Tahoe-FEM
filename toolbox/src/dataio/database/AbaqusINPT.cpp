/* created: S.A. Wimmer, Oct 2006    */
/* assumption: 
   1. One *NODE, optionally using *NSET
   2. *ELEMENT optionally using *ELSET
   3. GENERATE is only allowed in *NSET and *ELSET
   4. Make sure there is an /n at the last line 
   */
/* need to do to make complete
	1. Implement reading of *SURFACE for side set information.
	2. I don't think I accounted for all possible element types (higher order)
	*/

#include "AbaqusINPT.h"
#include "ifstreamT.h"
#include <ctype.h>
#include "iAutoArrayT.h"
#include "iArray2DT.h"
#include "GeometryT.h"

using namespace Tahoe;

/* constructor */
AbaqusINPT::AbaqusINPT (void) :
	fNumNodes (0),
	fNumElSets (0),
	fNumNodeSets (0)
{
}

bool AbaqusINPT::OpenRead (const StringT& filename)
{
	fInputFile = filename;
	fInputFile.ToNativePathName();
	ifstreamT in (fInputFile);
	if (!in.is_open())
	{
		cout << "\n AbaqusINPT::Open: error opening file: " << in.filename() << endl;
		return false;
	}
	
	ScanFile ();
	return true;
}

void AbaqusINPT::NodeIDs (iArrayT& ids) const
{
	ids = fNodeIDs;
}

void AbaqusINPT::Coordinates (dArray2DT& coords) const
{
	// technically *NODE can appear multiple times
	// will assume it appears once for now
	// but write code so multiple *NODE can be accounted for in future
	// will assume no GENERATE options
	// need to account for 2D among 3D
	ifstreamT in (fInputFile);
	StringT keyword;
	while (NextKeyWord (in, keyword))
	{
		if (strncmp (keyword.Pointer(), "NODE", 4) == 0) 
		{
			// comma deliminated
			char c;
			int id;
			double *p = coords.Pointer();
			for (int i = 0; i < coords.MajorDim(); i++)
			{
				in >> id >> c >> *p >> c >> *(p+1) >> c;
				
				// is there a third coordinate to read?
				// if there is no comma, no third coordinate?
				if (c != ',') 
				{
					in.putback (c);
					*(p+2) = 0;
				}
				else
					in >> *(p+2);

				//cout << id << " " << c << " ";
				//coords.PrintRow(i,cout);
				
				p += 3;
			}
			return;
		}
	}
}

void AbaqusINPT::SetNames (ArrayT<StringT>& names, const char* setname, const char* settype) const
{
	ifstreamT in (fInputFile);
	StringT keyword;
	int count = 0;
	StringT temp;
	int len1 = strlen (setname);
	int namelength = strlen (settype);
	while (NextKeyWord (in, keyword))
	{
		if (strncmp (keyword.Pointer(), settype, namelength) == 0)
		{
			if (ExtractOptionName (keyword, setname, temp))
				names[count++] = temp;
		}
		else if (strncmp (keyword.Pointer(), setname, len1) == 0)
		{
			StringT te = keyword;
			te.Drop (len1+1);
			if (ExtractOptionName (te, setname, temp))
				names[count++] = temp;
		}
	}
}

int AbaqusINPT::NodeSetLength (const StringT& name) const
{
	iArrayT set;
	NodeSet (name, set);
	return set.Length();
}

int AbaqusINPT::ElSetLength (const StringT& name) const
{
	iArrayT set;
	ElSet (name, set);
	return set.Length();
}

void AbaqusINPT::NodeSet (const StringT& name, iArrayT& set) const
{
	iAutoArrayT s;
	set.Free ();
	int index;
	if (fNodeSets.HasValue (name, index))
	{
		s = fNodeSetIDs[index];
		for (int j=0; j < fSetsinNSets[index].Length(); j++)
		{
			iArrayT temp;
			NodeSet (fSetsinNSets[index][j], temp);
			s.AppendUnique (temp);
		}
	}
	if (s.Length() > 0)
	{
		set.Dimension (s.Length());
		int *p = set.Pointer(), *q = s.Pointer();
		for (int h=0; h < set.Length(); h++)
			*p++ = *q++;
		set.SortAscending ();
	}
}

void AbaqusINPT::ElSet (const StringT& name, iArrayT& set) const
{
	iAutoArrayT s;
	set.Free ();
	int index;
	if (fElementSets.HasValue (name, index))
	{
		s = fElementIDs[index];
		for (int j=0; j < fSetsinElSets[index].Length(); j++)
		{
			iArrayT temp;
			ElSet (fSetsinElSets[index][j], temp);
			s.AppendUnique (temp);
		}
	}
	if (s.Length() > 0)
	{
		set.Dimension (s.Length());
		int *p = set.Pointer(), *q = s.Pointer();
		for (int h=0; h < set.Length(); h++)
			*p++ = *q++;
		set.SortAscending ();
	}
}

int AbaqusINPT::NumElementNodesforSet (const StringT& name) const
{
	int index;
	if (fElementSets.HasValue(name, index))
	{
		int nodes = fElementSetNodes[index];
		if (nodes > -1)
			return nodes;
		
		// see if set within set has this data
		for (int j=0; j < fSetsinElSets[index].Length(); j++)
		{
			nodes = NumElementNodesforSet (fSetsinElSets[index][j]);
			if (nodes > -1) 
				return nodes;
		}
		
		// deeper quering
		for (int e=0; e < fElementIDs[index].Length(); e++)
		{
			int target = fElementIDs[index][e];
			for (int k=0; k < fNumElSets; k++)
			{
				if (k != index && fElementIDs[k].HasValue (target))
				{
					nodes = fElementSetNodes[k];
					if (nodes > -1)
						return nodes;
				}
			}
		}
	}
	return 0;
}

bool AbaqusINPT::Connectivity (const StringT& name, iArray2DT& conn) const
{
	// check the name to see what kind of set it is (ELEMENT, or ELSET)
	int index;
	if (!fElementSets.HasValue(name, index)) return false;
	
	if (fSetsinElSets[index].Length() > 0)
	{
		int* r = conn.Pointer();
		for (int s=0; s < fSetsinElSets[index].Length(); s++)
		{
			iArrayT set;
			ElSet (fSetsinElSets[index][s], set);
			if (set.Length() < 1) return false;
			iArray2DT subconn (set.Length(), conn.MinorDim());
			if (!Connectivity (fSetsinElSets[index][s], subconn)) return false;
			int* sub = subconn.Pointer();
			for (int h=0; h < subconn.Length(); h++)
				*r++ = *sub++;
		}
		return true;
	}
	else
	{
		ifstreamT in (fInputFile);
		StringT keyword;
		bool good = NextKeyWord (in, keyword);
		while (good)
		{
			if (strncmp (keyword.Pointer(), "ELEMENT", 7) == 0)
			{
				StringT temp;
				if (ExtractOptionName (keyword, "ELSET", temp))
				{
					if (strncmp (temp.Pointer(), name.Pointer(), name.Length()) == 0)
					{
						int *p = conn.Pointer();
						int id;
						char c;
						for (int i=0; i < conn.MajorDim(); i++)
						{
							in >> id >> c;
							for (int j=0; j < conn.MinorDim() - 1; j++)
								in >> *p++ >> c;
							in >> *p++;
						}
					return true;
					}
				}
			}
			good = NextKeyWord (in, keyword);
		}
	}
	
	// need to also take into account ELSETS with element ids
	
	return false;
}

GeometryT::CodeT AbaqusINPT::GeometryCode (const StringT& name)
{
	ifstreamT in (fInputFile);
	StringT keyword, temp;
	int numelnodes;
	while (NextKeyWord (in, keyword))
		if (strncmp (keyword.Pointer(), "ELEMENT", 7) == 0)
			if (ExtractOptionName (keyword, "TYPE", temp))
				return TranslateElementName (temp, numelnodes);
	return GeometryT::kNone;
}

/*        private          */

void AbaqusINPT::ScanFile (void)
{
	fNumNodes = 0;
	fNumElSets = 0;
	fNumNodeSets = 0;
	fNumSurfaces = 0;
	ifstreamT in (fInputFile);
	StringT keyword;
	while (NextKeyWord (in, keyword))
	{
		if (strncmp (keyword.Pointer(), "NODE", 4) == 0 &&
		    strncmp (keyword.Pointer(), "NODE OUTPUT", 9) != 0 &&
		    strncmp (keyword.Pointer(), "NODE PRINT", 9) != 0 &&
		    strncmp (keyword.Pointer(), "NODE FILE", 9) != 0 &&
		    strncmp (keyword.Pointer(), "NODE RESPONSE", 9) != 0)
		{
			if (keyword.StringMatch ("NSET")) fNumNodeSets ++;

			iAutoArrayT ids;
			char n = in.next_char();
			int id;
			StringT line;
			while (n != '*')
			{
				in >> id;
				line.GetLineFromStream (in);
				ids.Append (id);
				n = in.next_char();
			}
				
			fNumNodes = ids.Length();
			fNodeIDs.Dimension (fNumNodes);
			int *nn = fNodeIDs.Pointer();
			int *f = ids.Pointer();
			for (int j=0; j < fNumNodes; j++)
				*nn++ = *f++;
		}
		else if (strncmp (keyword.Pointer(), "ELEMENT", 7) == 0 &&
		         strncmp (keyword.Pointer(), "ELEMENT MATRIX", 12) != 0 &&
		         strncmp (keyword.Pointer(), "ELEMENT OUTPUT", 12) != 0 &&
		         strncmp (keyword.Pointer(), "ELEMENT PROPERTIES", 12) != 0 &&
		         strncmp (keyword.Pointer(), "ELEMENT RESPONSE", 12) != 0)
		{	
			if (keyword.StringMatch ("ELSET")) fNumElSets ++;
		}
		else if (strncmp (keyword.Pointer (), "ELSET", 5) == 0) fNumElSets ++;
		
		else if (strncmp (keyword.Pointer (), "NSET", 4) == 0) fNumNodeSets ++;
		
		else if (strncmp (keyword.Pointer (), "SURFACE", 7) == 0 &&
		         strncmp (keyword.Pointer (), "SURFACE INTERACTION", 15) != 0 &&
		         strncmp (keyword.Pointer (), "SURFACE BEHAVIOR", 15) != 0 &&
		         strncmp (keyword.Pointer (), "SURFACE FLAW", 12) != 0 &&
		         strncmp (keyword.Pointer (), "SURFACE PROPERTY", 15) != 0 &&
		         strncmp (keyword.Pointer (), "SURFACE SECTION", 15) != 0) fNumSurfaces ++;
	}

	fElementSets.Dimension (fNumElSets);
	fElementSetNodes.Dimension (fNumElSets);
	fNodeSets.Dimension (fNumNodeSets);
	fElementIDs.Dimension (fNumElSets);
	fNodeSetIDs.Dimension (fNumNodeSets);
	fSetsinElSets.Dimension (fNumElSets);
	fSetsinNSets.Dimension (fNumNodeSets);
	
	cout << "\n  Number of NODES   : " << fNumNodes << endl;
	cout <<   "  Number of ELSETS  : " << fNumElSets << endl;
	cout <<   "  Number of NSETS   : " << fNumNodeSets << endl;
	cout <<   "  Number of SURFACES: " << fNumSurfaces << endl << endl;
	
	GatherSimpleData ();
}

void AbaqusINPT::GatherSimpleData (void)
{
	ifstreamT in (fInputFile);
	StringT keyword;
	bool good = NextKeyWord (in, keyword);
	int el = 0;
	int n = 0;
	while (good)
	{
		if (strncmp (keyword.Pointer(), "NODE", 4) == 0 &&
		    strncmp (keyword.Pointer(), "NODE OUTPUT", 9) != 0 &&
		    strncmp (keyword.Pointer(), "NODE PRINT", 9) != 0 &&
		    strncmp (keyword.Pointer(), "NODE FILE", 9) != 0 &&
		    strncmp (keyword.Pointer(), "NODE RESPONSE", 9) != 0)
		{
			StringT nna;
			if (ExtractOptionName (keyword, "NSET", nna))
			{
				fNodeSetIDs[n] = fNodeIDs;
				fSetsinNSets[n].Dimension (0);
				fNodeSets[n++] = nna;
			}
			good = NextKeyWord (in, keyword);
		}
		else if (strncmp (keyword.Pointer(), "ELEMENT", 7) == 0 &&
		         strncmp (keyword.Pointer(), "ELEMENT MATRIX", 12) != 0 &&
		         strncmp (keyword.Pointer(), "ELEMENT OUTPUT", 12) != 0 &&
		         strncmp (keyword.Pointer(), "ELEMENT PROPERTIES", 12) != 0 &&
		         strncmp (keyword.Pointer(), "ELEMENT RESPONSE", 12) != 0)
		{
			StringT na;
			int numelnodes = 0;
			if (ExtractOptionName (keyword, "TYPE", na))
				GeometryT::CodeT co = TranslateElementName (na, numelnodes);
		
			iAutoArrayT ids;
			char n = in.next_char();
			int id, nid;
			char commaspace;
			while (n != '*')
			{
				in >> id >> commaspace;
				for (int j=0; j < numelnodes - 1; j++)
					in >> nid >> commaspace;
				in >> nid;
				ids.Append (id);
				n = in.next_char();
			}
						
			if (ExtractOptionName (keyword, "ELSET", na))
			{
				fElementSets[el] = na;
				fElementSetNodes[el] = numelnodes;
				fElementIDs[el].Dimension (ids.Length());
				fSetsinElSets[el].Dimension (0);
				int *e = fElementIDs[el].Pointer();
				int *f = ids.Pointer();
				for (int j=0; j < ids.Length(); j++)
					*e++ = *f++;
				el++;
			}

			StringT line;
			line.GetLineFromStream(in);
			if (line.Length() < 1)
				good = false;
			else
			{
				line.Drop (1);
				keyword = line;
				keyword.ToUpper();
			}
		}
		else if (strncmp (keyword.Pointer(), "ELSET", 5) == 0)
		{
			StringT line;
			iAutoArrayT ids;
			sArrayT names;
			
			if (keyword.StringMatch ("GENERATE"))
				ReadGeneratedSetData (in, ids, line);
			else
				ReadSetData (in, ids, names, line);

			StringT na;
			if (ExtractOptionName (keyword, "ELSET", na))
			{
				fElementSets[el] = na;
				fElementSetNodes[el] = -1;
				fElementIDs[el].Dimension (ids.Length());
				fSetsinElSets[el] = names;
				for (int j=0; j < ids.Length(); j++)
					fElementIDs [el][j] = ids[j];
				
				el++;
			}

			if (line.Length() < 1)
				good = false;
			else
			{
				line.Drop (1);
				keyword = line;
				keyword.ToUpper();
			}
		}
		else if (strncmp (keyword.Pointer(), "NSET", 4) == 0)
		{
			StringT line;
			iAutoArrayT ids;
			sArrayT names;
			
			if (keyword.StringMatch ("GENERATE"))
				ReadGeneratedSetData (in, ids, line);
			else
				ReadSetData (in, ids, names, line);

			StringT na;
			if (ExtractOptionName (keyword, "NSET", na))
			{
				fNodeSetIDs[n].Dimension (ids.Length());
				for (int j=0; j < ids.Length(); j++)
					fNodeSetIDs[n][j] = ids[j];
				fSetsinNSets[n] = names;
				fNodeSets[n++] = na;
			}

			if (line.Length() < 1)
				good = false;
			else
			{
				line.Drop (1);
				keyword = line;
				keyword.ToUpper();
			}
		}
		else if (strncmp (keyword.Pointer(), "SURFACE", 7) == 0 &&
		         strncmp (keyword.Pointer (), "SURFACE INTERACTION", 15) != 0 &&
		         strncmp (keyword.Pointer (), "SURFACE BEHAVIOR", 15) != 0 &&
		         strncmp (keyword.Pointer (), "SURFACE FLAW", 12) != 0 &&
		         strncmp (keyword.Pointer (), "SURFACE PROPERTY", 15) != 0 &&
		         strncmp (keyword.Pointer (), "SURFACE SECTION", 15) != 0)
		{
			good = NextKeyWord (in, keyword);
			;// future
		}
		else
			good = NextKeyWord (in, keyword);

	}
	/*for (int e=0; e < fElementSets.Length(); e++)
	{
		cout << "ELSET " << fElementSets[e] << " length " << fElementIDs[e].Length() << " "
			 << "sets " << fSetsinElSets[e].Length() << " nodes " << fElementSetNodes[e];
		cout << " elsetlength " << ElSetLength(fElementSets[e]) << " elsetnodes "
			 << NumElementNodesforSet (fElementSets[e]) << endl;
	}
	cout << endl;
	for (int n=0; n < fNodeSets.Length(); n++)
	{
		cout << "NSET " << fNodeSets[n] << " length " << fNodeSetIDs[n].Length() << " sets "
			 << fSetsinNSets[n].Length();
		cout << " nsetlength " << NodeSetLength (fNodeSets[n]) << endl;
	}*/
}

void AbaqusINPT::ReadSetData(ifstream& in, iAutoArrayT& ids, sArrayT& names, StringT& line)
{
	// assume no GENERATE
	// need to account for an ELSET/NSET listed in an ELSET/NSET
	// need to account for comma delimination but no space after comma, no comma at end of line
	StringT temp;
	line.GetLineFromStream(in);
	int nc = 0;
	while (line[0] != '*' && in.good())
	{
		temp.Clear();
		for (int i=0; i < line.Length(); i++)
		{
			if (line[i] != ',' && line[i] != ' ' && line[i] != '\0' && 
				line[i] != '\n' && line[i] != '\t')
			{
				temp.Append (line[i]);
			}
			else 
			{
				if (temp.Length() > 0)
				{
					if (isalpha(temp[0]))
					{
						names.Resize (nc+1);
						names[nc] = temp;
						names[nc].ToUpper();
						nc++;
					}
					else
						ids.Append (atoi (temp.Pointer()));
					temp.Clear();
				}
			}
		}
		line.GetLineFromStream(in);
	}
}

void AbaqusINPT::ReadGeneratedSetData (ifstreamT& in, iAutoArrayT& ids, StringT& line)
{
	// GENERATE
	int first, last, increment;
	char comma;
	char n = in.next_char();
	while (n != '*')
	{
		in >> first >> comma >> last >> comma >> increment;
		n = in.next_char ();
		for (int i=first; i < last+increment; i += increment)
			ids.Append (i);
	}
	line.GetLineFromStream(in);
}

bool AbaqusINPT::NextKeyWord (ifstream& in, StringT& keyword) const
{
	bool found = false;
	StringT line;
	while (in.good() && !found)
	{
		line.GetLineFromStream (in);
		if (line[0] == '*')
		{
			line.Drop (1);
			if (line[0] != '*') // account for comment lines
			{
				keyword = line;
				keyword.ToUpper();
				return true;
			}
		}
	}
	return false;
}

bool AbaqusINPT::ExtractOptionName (const StringT& keyword, const char* settype, StringT& name) const
{
	const char* p = keyword.StringMatch (settype);
	if (p == NULL) return false;
	name = p;
	int e = name.FirstPositionOf ('=');
	name.Drop (e+1);
	name.DropLeadingSpace();
	while (name.FirstPositionOf(',') != -1) 
		name.Root (',');
	name.DropTrailingSpace();
	return true;
}

GeometryT::CodeT AbaqusINPT::TranslateElementName (const char* name, int& numelnodes) const
{
	numelnodes = -1;
	if (strncmp (name, "C", 1) == 0)
		return TranslateContinuum (name+1, numelnodes);
	else if (strncmp (name, "DC", 2) == 0)
		return TranslateContinuum (name+2, numelnodes);
	else if (strncmp (name, "AC", 2) == 0)
		return TranslateContinuum (name+2, numelnodes);
	else if (strncmp (name, "DCC", 3) == 0)
		return TranslateContinuum (name+3, numelnodes);
	else if (strncmp (name, "SPRING", 6) == 0)
		return TranslateSpring (name+6, numelnodes);
	else if (strncmp (name, "S", 1) == 0)
		return TranslateShell (name+1, numelnodes);
	else if (strncmp (name, "R", 1) == 0)
		return TranslateRigid (name+1, numelnodes);
	else
		return GeometryT::kNone;
}

GeometryT::CodeT AbaqusINPT::TranslateContinuum (const char* name, int& numelnodes) const
{
  if (strncmp (name, "PE", 2) == 0)
    return Translate2D (name+2, numelnodes);
  else if (strncmp (name, "PS", 2) == 0)
    return Translate2D (name+2, numelnodes);
  else if (strncmp (name, "2D", 2) == 0)
    return Translate2D (name+2, numelnodes);
  else if (strncmp (name, "GPE", 3) == 0)
    return Translate2D (name+3, numelnodes);
  else if (strncmp (name, "3D", 2) == 0)
    return Translate3D (name+2, numelnodes);
  return GeometryT::kNone;
}

GeometryT::CodeT AbaqusINPT::TranslateSpring (const char* name, int& numelnodes) const
{
	if (name[0] == '1')
	{
		numelnodes = 1;
		return GeometryT::kPoint;
	}
	else if (name[0] == '2' || name [0] == 'A')
	{
		numelnodes = 2;
		return GeometryT::kLine;
	}
	else
		return GeometryT::kNone;
}

GeometryT::CodeT AbaqusINPT::Translate2D (const char* name, int& numelnodes) const
{
	if (name[0] == '3')
	{
		numelnodes = 3;
		return GeometryT::kTriangle;
	}
	else if (name[0] == '6')
	{
		numelnodes = 6;
		return GeometryT::kTriangle;
	}
	else if (name[0] == '4')
	{
		numelnodes = 4;
		return GeometryT::kQuadrilateral;
	}
	else if (name[0] == '8')
	{
		numelnodes = 8;
		return GeometryT::kQuadrilateral;
	}
	else
		return GeometryT::kNone;
}

GeometryT::CodeT AbaqusINPT::Translate3D (const char* name, int& numelnodes) const
{
	if (name[0] == '4')
	{
		numelnodes = 4;
		return GeometryT::kTetrahedron;
	}
	else if (name[0] == '6')
	{
		numelnodes = 6;
		return GeometryT::kPentahedron;
	}
	else if (name[0] == '8')
	{
		numelnodes = 8;
		return GeometryT::kHexahedron;
	}
	else if (strncmp (name, "10", 2) == 0)
	{
		numelnodes = 10;
		return GeometryT::kTetrahedron;
	}
	else if (strncmp (name, "15", 2) == 0)
	{
		numelnodes = 15;
		return GeometryT::kPentahedron;
	}
	else if (strncmp (name, "20", 2) == 0)
	{
		numelnodes = 20;
		return GeometryT::kHexahedron;
	}
	else		
		return GeometryT::kNone;
}

GeometryT::CodeT AbaqusINPT::TranslateShell (const char* name, int& numelnodes) const
{
	if (name[0] == '3')
	{
		numelnodes = 3;
		return GeometryT::kTriangle;
	}
	else if (name[0] == '4')
	{
		numelnodes = 4;
		return GeometryT::kQuadrilateral;
	}
	else if (strncmp (name, "8R", 2) == 0)
	{
		numelnodes = 8;
		return GeometryT::kQuadrilateral;
	}
	else
		return GeometryT::kNone;
}

GeometryT::CodeT AbaqusINPT::TranslateRigid (const char* name, int& numelnodes) const
{
	if (strncmp (name, "3D", 2) == 0)
		return Translate2D (name+2, numelnodes);
	else if (strncmp (name, "2D", 2) == 0)
	{
		numelnodes = 1;
		return GeometryT::kPoint;
	}
	else if (strncmp (name, "AX", 2) == 0)
	{
		numelnodes = 1;
		return GeometryT::kPoint;
	}
	else
		return GeometryT::kNone;
}

