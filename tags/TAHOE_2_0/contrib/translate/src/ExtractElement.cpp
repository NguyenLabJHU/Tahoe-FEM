/* $Id: ExtractElement.cpp,v 1.2 2003-02-25 14:34:36 sawimme Exp $ */
#include "ExtractElement.h"

using namespace Tahoe;

ExtractElement::ExtractElement (ostream& out, istream& in, bool write):
  ExtractIOManager (out, in, write)
{

}

/**************** PROTECTED **********************/

void ExtractElement::Initialize (void)
{
	InitializeElementVariables ();
	if (fNumEV < 1) {
		fMessage << "\n No element variables found.";
		return;
	}

	/* select element ID */
	const ArrayT<StringT>& elem_ID = fModel.ElementGroupIDs();
	int elem_dex = 0;

	while (elem_dex < 1 || elem_dex > elem_ID.Length()) {
		fMessage << " Number of element groups: " << elem_ID.Length() << '\n';
		for (int i = 0; i < elem_ID.Length(); i++)
			fMessage << setw(kIntWidth) << i+1 << ": " << elem_ID[i] << '\n';
		fMessage << "\n Index of the element group: ";

		fIn >> elem_dex;
		if (fEcho) fEchoOut << elem_dex << "\n";
	}
	elem_dex--;
	fExtract_ID = elem_ID[elem_dex];

	iArrayT elements;
	SelectElements(fExtract_ID, elements, fItemIndex);
	fNumItems = elements.Length();
	if (fNumItems < 1) {
		fMessage << "\n No elements found.";
		return;
	}

	fItemNames.Allocate(fNumItems);
	for (int i = 0; i < fNumItems; i++)
		fItemNames[i].Append(elements[i]);
}

void ExtractElement::TranslateVariables (void)
{
	PrepFiles(fEVUsed, fElementLabels);

	int nel, nen;
	fModel.ElementGroupDimensions(fExtract_ID, nel, nen);

	fVarData.Allocate(nel, fNumEV);
	for (int t = 0; t < fNumTS; t++)
	{
		fModel.ElementVariables(fTimeIncs[t], fExtract_ID, fVarData);
		WriteVarData(fEVUsed, t);
	}
}

