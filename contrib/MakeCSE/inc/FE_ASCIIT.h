// file: FE_ASCIIT.h

// created      : SAW (05/20/1999)
// last modified: PAK (11/08/1999)

#ifndef _FE_ASCII_T_H_
#define _FE_ASCII_T_H_

/* base class */
#include "OutputBaseT.h"

/* direct members */
#include <fstream.h>

class FE_ASCIIT: public OutputBaseT
{
  public:

	/* constructor */
	FE_ASCIIT(ostream& out, bool external, const ArrayT<StringT>& out_strings);

	/* destructor */
	~FE_ASCIIT(void);

	/* output functions */
	virtual void WriteGeometry(void);
	virtual void WriteOutput(double time, int ID, const dArray2DT& n_values, 
		const dArray2DT& e_values);

  private:

	/* open streams */
	void OpenStreams(void);

	/* set-by-set output */
	void WriteGeometryData(int ID);
	void WriteOutputData(int ID, const dArray2DT& n_values, 
		const dArray2DT& e_values);

	/* variable headers */
	void WriteNodeHeader(ostream& out, int num_output_nodes, 
		const ArrayT<StringT>& labels) const;
	void WriteElementHeader(ostream& out, int num_output_elems, 
		const ArrayT<StringT>& labels) const;

	/* writing data */
	void WriteNodeValues(ostream& out, const iArrayT& node_numbers, 
		const dArray2DT& values) const;
  	void WriteElementValues(ostream& out, const dArray2DT& values) const;

#if 0
	void PrintSingleTextFile(StringT& file, const iArray2DT& a) const;
	void PrintSingleTextFile(StringT& file, const iArrayT& a) const;
#endif

  private:

	ofstream fmovie;
	ofstream fgeo;
	bool fExternTahoeII;
};

#endif // _FE_ASCII_T_H_
