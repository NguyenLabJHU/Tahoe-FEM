/* $Id: LocalizerT.cpp,v 1.11 2004-06-17 07:40:07 paklein Exp $ */
/* created: paklein (02/19/1998) */
#include "LocalizerT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ifstreamT.h"
#include "toolboxConstants.h"
#include "SolidMaterialT.h"
#include "SolidMatList2DT.h"
#include "SolidMatList3DT.h"
#include "ShapeFunctionT.h"
#include "EdgeFinderT.h"
#include "iAutoArrayT.h"
#include "ExodusT.h"

//for strain check below
#include "FSSolidMatT.h"

/* flag parameter to mark localized elements */

using namespace Tahoe;

const int kMonitorLocalized = -1;

/* strain check flags */
const int kStrainCheckNever =-1;
const int kStrainCheckRHS   = 0;
const int kStrainCheckLHS   = 1;
const int kStrainCheckEvery = 2;

/* crack tracking increment flags (values > 0 imply check at interval) */
const int kLocCheckNever   = 0;
const int kLocCheckAtPrint =-1;

/* constructor */
LocalizerT::LocalizerT(const ElementSupportT& support, const FieldT& field):
	UpdatedLagrangianT(support, field),
	fAvgStretch(NumSD())
{
	/* flags */
	ifstreamT& in = ElementSupport().Input();

	in >> fStrainCheckFlag;
	in >> fCriticalStretch;
	in >> fLocCheckInc;

	if (fStrainCheckFlag != kStrainCheckNever &&
	    fStrainCheckFlag != kStrainCheckRHS   &&
	    fStrainCheckFlag != kStrainCheckLHS   &&
	    fStrainCheckFlag != kStrainCheckEvery)
		throw ExceptionT::kBadInputValue;
	
	if (fStrainCheckFlag != kStrainCheckNever &&
	    fCriticalStretch < 1.0)
		throw ExceptionT::kBadInputValue;
	
	if (fLocCheckInc != kLocCheckNever   &&
	    fLocCheckInc != kLocCheckAtPrint &&
	    fLocCheckInc < 0)
		throw ExceptionT::kBadInputValue;
}

/* set work space */
void LocalizerT::Initialize(void)
{
	/* inherited */
	UpdatedLagrangianT::Initialize();

	/* dimension */
	fKloc.Dimension(NumElementNodes()*NumDOF());
	
	/* echo group-specific input data */
	EchoData(ElementSupport().Input(), ElementSupport().Output());

	/* dimension */	
	fElementMonitor.Resize(NumElements(), false);

//TEMP - needs rethinking
#if 0
	/* check material list for localizing materials */
	if (!fMaterialList->HasLocalizingMaterials())
		cout << "\n LocalizerT::Initialize: WARNING: no localizing materials" << endl;
#endif

	/* open output stream for localization data */
	if (fLocCheckInc != kLocCheckNever)
	{
		int groupnum = ElementSupport().ElementGroupNumber(this) + 1;
	
		StringT outfile;
		outfile.Root(ElementSupport().Input().filename());
		outfile.Append(".loc.elem", groupnum);
		fLocOut.open(outfile);

		outfile.Root(ElementSupport().Input().filename());
		outfile.Append(".TOC.elem", groupnum);
		fLocTOC.open(outfile);
	}		

	/* localization check work space */
	iArrayT eq_temp;
	switch (GeometryCode())
	{
		case GeometryT::kQuadrilateral:
		{
			int quad_eqs[5] = {2,4,5,6,7};
			eq_temp.Set(5,quad_eqs);
			break;
		}
		case GeometryT::kTriangle:
		{
			int tri_eqs[3] = {2,4,5};
			eq_temp.Set(3,tri_eqs);
			break;
		}
		case GeometryT::kHexahedron:
		{
			int hex_eqs[18] = {4,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
			eq_temp.Set(18,hex_eqs);
			break;
		}
		case GeometryT::kTetrahedron:
		{
			int tet_eqs[6] = {3,6,7,9,10,11};
			eq_temp.Set(6,tet_eqs);
			break;
		}
		default:

			throw ExceptionT::kGeneralFail;
	}
	fnoRBeqs = eq_temp;

	/* dimension eigenvalue solver */
	fEigSolver.Dimension(fnoRBeqs.Length());
	fEigs.Dimension(fnoRBeqs.Length());
}

/* finalize time increment */
void LocalizerT::CloseStep(void)
{
	/* inherited */
	UpdatedLagrangianT::CloseStep();
	
	/* set flag to OFF for elements marked last time */
	fElementMonitor.MarkedToOFF();
}

/* writing results */
void LocalizerT::WriteOutput(void)
{	
	/* inherited */
	UpdatedLagrangianT::WriteOutput();

	ostream& out = ElementSupport().Output();

	/* insert strain check info */
	if (fStrainCheckFlag != kStrainCheckNever)
	{
		out << "\n Elements marked as OFF:\n\n";
		fElementMonitor.PrintOFF(out);		
	}
	
	/* insert localization information */
	if (fLocCheckInc != kLocCheckNever)
	{
		/* search for localization */
		if (fLocCheckInc == kLocCheckAtPrint) Localization();
		
		/* print localized element list (at every call) */
		PrintLocalized(out);
		
		/* output search info to localization file */
		//if (mode == IOBaseT::kAtFinal)
		if (true)
		{
			fLocOut << "\n Localization search data:\n\n";
			fLocOut << " Number of check elements. . . . . . . . . . . . = ";
			fLocOut << fLocCheckList.Length() << '\n';
	
			if (fLocCheckList.Length() > 0)
			{
				/* header */
				fLocOut << '\n';
				fLocOut << setw(kIntWidth) << "no.";
				fLocOut << setw(kIntWidth) << "element";
				fLocOut << setw(kIntWidth) << "status" << '\n';
				
				for (int i = 0; i < fLocCheckList.Length(); i++)
				{
					int eltag = fLocCheckList[i];
					
					fLocOut << setw(kIntWidth) <<     i + 1;
					fLocOut << setw(kIntWidth) << eltag + 1;
					fLocOut << setw(kIntWidth) << fElementMonitor.Status(eltag) << '\n';
				}

				fLocOut << '\n';
			}
		}
	}
}

/* returns true if the internal force has been changed since
* the last time step */
GlobalT::RelaxCodeT LocalizerT::RelaxSystem(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = UpdatedLagrangianT::RelaxSystem();

	/* check localization at intervals */
	if (++fLocCheckCount == fLocCheckInc)
	{
		Localization();
		fLocCheckCount = 0;
	}

	/* element status has changed */
	if (fElementMonitor.Changed())
	{
		int numMARKED = fElementMonitor.Count(MonitorT::kMarked);
	
		if (numMARKED > 0)
		{
			/* to output file */
			ostream& out = ElementSupport().Output();
			out << "\n Marked elements:\n";
			fElementMonitor.PrintValued(cout, MonitorT::kMarked, 5);
			out << endl;
		
			return GlobalT::MaxPrecedence(relax, GlobalT::kRelax);
		}
		else
			/* reset changed flag */
			fElementMonitor.Reset(); //why this??
	}
	return GlobalT::MaxPrecedence(relax, GlobalT::kNoRelax);
}

/* initial condition/restart functions
*
* Set to initial conditions.  The restart functions
* should read/write any data that overrides the default
* values */
void LocalizerT::InitialCondition(void)
{
	/* inherited */
	UpdatedLagrangianT::InitialCondition();

	fLocCheckCount = 0;

	fElementMonitor.AllToON();
}

void LocalizerT::ReadRestart(istream& in)
{
	in >> fLocCheckCount;

	fElementMonitor.ReadStatus(in);
	
	/* rebuild the localization check list */
	if (fLocCheckInc == kLocCheckAtPrint || fLocCheckInc > 0)
	{
		/* set initial localization check list */
		cout << "\n Initializing localization check list: ";
		for (int j = 0; j < NumElements(); j++)
			if (fElementMonitor.Status(j) == kMonitorLocalized)
				fLocCheckList.Append(j);
		cout << fLocCheckList.Length() << endl;

		/* (re-)build element neighbor list */
		iArray2DT nodefacetmap;
		fShapes->NeighborNodeMap(nodefacetmap);
		cout << " Rebuilding element neighbor list: ";
		EdgeFinderT edger(fConnectivities, nodefacetmap);
		fNeighborList = edger.Neighbors();
		cout << "done" << endl;

		/* append neighbors to the check list */
		cout << " Finalizing localization check list: ";
		int init_list_length = fLocCheckList.Length();
		for (int i = 0; i < init_list_length; i++)
		{
			int elnum = fLocCheckList[i];

			int* p = fNeighborList(elnum);
			for (int i = 0; i < fNeighborList.MinorDim(); i++)
			{
				if (*p > -1) fLocCheckList.AppendUnique(*p);
				p++;
			}
		}	
		cout << fLocCheckList.Length() << endl;
	}
}

void LocalizerT::WriteRestart(ostream& out) const
{
	out << fLocCheckCount << '\n';

	fElementMonitor.WriteStatus(out);
}

/* resets to the last converged solution */
GlobalT::RelaxCodeT LocalizerT::ResetStep(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = UpdatedLagrangianT::ResetStep();

	/* reset flagged elements */
	fElementMonitor.MarkedToON();
	fElementMonitor.ResetFlags(kMonitorLocalized, MonitorT::kON);

	return relax;
}

/***********************************************************************
* Protected
***********************************************************************/

/* print element group data */
void LocalizerT::PrintControlData(ostream& out) const
{
	/* inherited */
	UpdatedLagrangianT::PrintControlData(out);
	
	/* control parameters */
	out << " Strain check flag . . . . . . . . . . . . . . . = " << fStrainCheckFlag << '\n';
	out << "    eq." << kStrainCheckNever << ", never (skip all strain checks)\n";
	out << "    eq." << kStrainCheckRHS   << ", while forming RHS\n";
	out << "    eq." << kStrainCheckLHS   << ", while forming LHS\n";
	out << "    eq." << kStrainCheckEvery << ", while forming LHS and RHS (ALL)\n";
	out << " Critical stretch. . . . . . . . . . . . . . . . = " << fCriticalStretch << '\n';
	out << " Crack tracking increment. . . . . . . . . . . . = " << fLocCheckInc << '\n';
	out << "    eq." << kLocCheckAtPrint << ", at kinematic print increments\n";
	out << "    eq." << kLocCheckNever   << ", never (skip all localization checks)\n";
	out << "    eq." << fLocCheckInc	 << ", at every nth step\n";
}

/* skip elements that are off (when?) */
bool LocalizerT::NextElement(void)
{
	/* inherited */
	bool notend = UpdatedLagrangianT::NextElement();
	
	/* skip OFF elements */
	if (fStrainCheckFlag != kStrainCheckNever)
	{
		while (notend &&
			fElementMonitor.Status(fElementCards.Position()) == MonitorT::kOFF)
			notend = UpdatedLagrangianT::NextElement();
	}
		
	return notend;
}

/* construct materials manager and read data */
void LocalizerT::ReadMaterialData(ifstreamT& in)
{
	/* inherited */
	UpdatedLagrangianT::ReadMaterialData(in);
	
//TEMP - needs rethinking
#if 0
	/* check for localizing materials */
	if (!fMaterialList->HasLocalizingMaterials())
	{
		const char message[] =
			"\n LocalizerT::ReadMaterialData: WARNING: no localizing materials";
	
		cout << message << endl;
		ElementSupport().Output() << message << endl;
	}
#endif	
}

/* form the element stiffness matrix */
void LocalizerT::FormStiffness(double constK)
{		
	/* inherited */
	UpdatedLagrangianT::FormStiffness(constK);
	
	/* disable the rest */
	return;
	
// CheckLocalization now write a TOC file to
// interpret the fLocOut file. Not implemented
// for stiffness yet	
	
	if (fLocCheckInc == kLocCheckNever) return;

	/* acoustical tensor check */
	int Qcheck = 0;
	if (CheckLocalization(fLocOut) > 0)
	{
		/* set flag */
		Qcheck = 1;

		fLocOut << "\n Acoustical tensor is NOT positive-definite:\n";
		fLocOut << " orientation of the normal: \n" << fNormal << endl;
	}
		
	/* initialize eigenvalue solver */
	fKloc = fLHS;
	fKloc.CopySymmetric();
	fKloc.CopyBlock(fnoRBeqs,fEigSolver);	
	fEigSolver.Diagonalize(fEigs);
			
	/* significant negative eigenvalue */
	int Kcheck = 0;
	if (fEigs.Min() < -1.0e-6)
	{		
		/* set flag */
		Kcheck = 1;

		/* sort */
		fLocOut << "\n eigenvalues:\n";
		fEigs.SortDescending();
		fLocOut << fEigs.wrap(6) << '\n';
		fLocOut << '\n';

		double neg_eig = fEigs.Min();

		/* int pt stiffness matrix */
		fLocOut << " K :\n\n";
		fKloc.WriteWithFormat(fLocOut,15,12,1);

		/* get material state */
		const dSymMatrixT& sij = fCurrMaterial->s_ij();
		const dMatrixT& cijkl = fCurrMaterial->c_ijkl();

		/* print out stress and moduli */
		fLocOut << "\n Cauchy stress:\n" << sij	<< '\n';
		fLocOut << "\n Material tangent moduli:\n" << cijkl << '\n';

		fKloc.CopyBlock(fnoRBeqs,fEigSolver);
		fEigSolver.Eigenvector(neg_eig,fEigs);
		
		/* full length */
		dArrayT fullvec(fRHS.Length());
		fullvec = 0.0;
		for (int i = 0; i < fEigs.Length(); i++)
			fullvec[fnoRBeqs[i]] = fEigs[i];
		
		fLocOut << "\n Eigenvector: \n";
		fLocOut << fullvec << '\n';

		fLocCurrCoords.ReturnTranspose(fullvec);
		fLocOut << "\n current coords: \n";
		fLocOut << fullvec << '\n';		
	}
		
	if (Kcheck || Qcheck) throw ExceptionT::kBadJacobianDet;
}

/***********************************************************************
* Private
***********************************************************************/

/* element specific input */
void LocalizerT::EchoData(ifstreamT& in, ostream& out)
{
	/* read initial localization check elements */
	int numlocelems;
	in >> numlocelems;
	if (numlocelems < 0) throw ExceptionT::kBadInputValue;
	
	/* indicates "smart" element checking */
	if (numlocelems > 0)
	{
		/* read element list */
		fLocCheckList.Dimension(numlocelems);
		for (int i = 0; i < numlocelems; i++)
		{
			int temp;
			in >> temp;
						
			/* range check */
			if (temp > NumElements())
			{
				cout << "\n LocalizerT::EchoSpecialData: element localization list";
				cout << " member " << temp << " is out of range" << endl;
				throw ExceptionT::kBadInputValue;
			}
			
			/* add to list */
			fLocCheckList[i] = (--temp);
		}

		/* build element neighbor list */
		iArray2DT nodefacetmap;
		fShapes->NeighborNodeMap(nodefacetmap);
		EdgeFinderT edger(fConnectivities, nodefacetmap);
		fNeighborList = edger.Neighbors();

		/* append first neighbors to the check list */
		iArrayT neighs;
		for (int j = 0; j < numlocelems; j++)
		{
			fNeighborList.RowAlias(fLocCheckList[j], neighs);
		
			for (int i = 0; i < neighs.Length(); i++)
				if (neighs[i] > -1)
					fLocCheckList.AppendUnique(neighs[i]);
		}
	}
	
	/* echo to output */
	out << "\n Localization checking:\n\n";
	out << " Initial number of check elements. . . . . . . . = ";
	if (numlocelems == 0) out << "<ALL>" << "\n\n";
		else out << fLocCheckList.Length() << "\n\n";

	/* output (corrected) element numbers */
	iArrayT temp(fLocCheckList.Length(), fLocCheckList.Pointer());
	temp++;
	out << temp.wrap(6) << '\n';
	temp--;
	out << '\n';
}

/* element localization check driver - loop over all/list */
void LocalizerT::Localization(void)
{
//TEMP - needs rethinking
#if 0
	/* check material behavior */
	if (!fMaterialList->HasLocalizingMaterials()) return;
#endif

	/* write the time */
	fLocOut << ElementSupport().Time() << '\n';

	/* check ALL */
	if (fLocCheckList.Length() == 0)
	{
		/* TOC entries */
		iArrayT numlocip(NumElements());
		numlocip = 0;
	
		int* pnumloc = numlocip.Pointer();
	 	Top();
		while ( NextElement() )
			*pnumloc++ = CheckLocalization(fLocOut);
			
		/* write TOC */
		int numOK = numlocip.Count(0);
		if (numOK < NumElements())
		{
			int numloc = NumElements() - numOK;
			fLocTOC << numloc << '\n';
			
			int* ploc = numlocip.Pointer();
			for (int i = 0; i < NumElements(); i++)
			{
				if (*ploc > 0)
					fLocTOC << i+1 << " " << *ploc << '\n';
					// write element number and the number of localized ip's
					
				ploc++;
			}
		}
		else
			fLocTOC << 0 << '\n';					
	}
	/* only elements in the list */
	else
	{
		/* TOC entries */
		iArrayT numlocip(fLocCheckList.Length());
		numlocip = 0;		
	
		for (int i = 0; i < fLocCheckList.Length(); i++)
		{
			int elnum = fLocCheckList[i];
	
			/* set current pointers */
			ElementCardT& element = fElementCards[elnum];

			/* cast is safe since class contructs materials list */
			ContinuumMaterialT* pcont_mat = (*fMaterialList)[element.MaterialNumber()];
			fCurrMaterial = (SolidMaterialT*) pcont_mat;
	
			/* check for localization */
			numlocip[i] = CheckLocalization(fLocOut);
			
			/* first time localized */
			if (numlocip[i] > 0 &&
			    fElementMonitor.Status(elnum) != kMonitorLocalized)
			{
				/* mark element */
				fElementMonitor.SetFlag(elnum,kMonitorLocalized);
			
				/* append neighbors to the check list */
				int* p = fNeighborList(elnum);
				for (int i = 0; i < fNeighborList.MinorDim(); i++)
				{
					if (*p > -1) fLocCheckList.AppendUnique(*p);					
					p++;
				}
			}
			
			/* expand TOC for added elements */
			if (fLocCheckList.Length() > numlocip.Length())
				numlocip.Resize(fLocCheckList.Length());
		}
		
		/* write TOC */
		int numOK = numlocip.Count(0);
		if (numOK < numlocip.Length())
		{
			int numloc = numlocip.Length() - numOK;
			fLocTOC << numloc << '\n';
			
			int* ploc = numlocip.Pointer();
			for (int i = 0; i < numlocip.Length(); i++)
			{
				if (*ploc > 0)
					fLocTOC << fLocCheckList[i]+1 << " " << *ploc << '\n';
					// write element number and the number of localized ip's
				
				ploc++;
			}
		}
		else
			fLocTOC << 0 << '\n';					
	}
}

/* check localization in current element the number of localized ip's
* and write localized ip coordinates to out */
int LocalizerT::CheckLocalization(ostream& out)
{
	/* compute global shape functions */
	SetGlobalShape();

	dArrayT vec(NumSD());
	int numloc = 0;
	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		/* set material state */
		const dSymMatrixT& sij = fCurrMaterial->s_ij();
		const dMatrixT& cijkl = fCurrMaterial->c_ijkl();
		
		/* localization check */
		if (fCurrMaterial->IsLocalized(fNormal))
		{
			numloc++;
			fShapes->IPCoords(vec);
			out << vec.no_wrap() << '\n';
		}
	}
	
	return numloc;
}	

/* output list of localized element numbers */
void LocalizerT::PrintLocalized(ostream& out)
{
	/* header */
	out << "\n Localized elements:\n\n";
	out << " Number of localized elements. . . . . . . . . . = ";
	out << fElementMonitor.Count(kMonitorLocalized) << "\n\n";

	/* output flagged */
	fElementMonitor.PrintValued(out, kMonitorLocalized, 5);


	//compile connectivity for the localized elements
	iAutoArrayT loc_elems;
	if (fLocCheckInc != kLocCheckNever)
		for (int i = 0; i < fLocCheckList.Length(); i++)
	    {
			int eltag = fLocCheckList[i];
			if (fElementMonitor.Status(eltag) == kMonitorLocalized)
				loc_elems.Append (eltag);
	    }

//NOTE - localized elements as separate "changing geometry" output
//       group not set yet
//	if (loc_elems.Length() > 0)
//	{
//		fOutputConn.Allocate (loc_elems.Length(), fConnectivities.MinorDim());
//		int* ploc_conn = fOutputConn.Pointer();
//	    for (int j = 0; j < loc_elems.Length(); j++)
//			for (int k = 0; k < fConnectivities.MinorDim(); k++)
//				*ploc_conn++ = fConnectivities (loc_elems[j], k);
//	}
}

/* check strain magnitude */
void LocalizerT::CheckStrain(void)
{
	/* work space */
	dSymMatrixT temp(NumSD());

	/* initialize */
	fAvgStretch = 0.0;

	/* check strain magnitude at every int pt */
	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		/* compute the stretch */
		temp.MultATA(DeformationGradient());
	
		/* accumulate over element */
		fAvgStretch += temp;
	}	
		
	/* element average */
	fAvgStretch	/= fShapes->NumIP();
	double avgstretch = sqrt(fabs(fAvgStretch.Invariant2()));
	if ( avgstretch > fCriticalStretch)
		fElementMonitor.Mark( fElementCards.Position() );
}
