/* $Id: FieldT.cpp,v 1.36.2.2 2004-11-15 04:15:03 d-farrell2 Exp $ */
#include "FieldT.h"

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "nIntegratorT.h"
#include "KBC_ControllerT.h"
#include "FBC_ControllerT.h"
#include "RaggedArray2DT.h"
#include "LinkedListT.h"
#include "LocalArrayT.h"
#include "FieldSupportT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "ModelManagerT.h"

using namespace Tahoe;

/* constructor */
FieldT::FieldT(const FieldSupportT& field_support):
	ParameterInterfaceT("field"),
	fFieldSupport(field_support),
	fGroup(-1),
	fIntegrator(NULL),	
	fnIntegrator(NULL),
	fEquationStart(0),
	fNumEquations(0)
{

}

/* configure the field */
void FieldT::Initialize(const StringT& name, int ndof, int order)
{
	/* initialize base class */
	BasicFieldT::Initialize(name, ndof, order);

	/* allocate history */
	fField_last.Dimension(order + 1);

	/* register arrays */
	for (int i = 0; i < fField_last.Length(); i++)
		RegisterArray2D(fField_last[i]);
	RegisterArray2D(fUpdate);
}

/* destructor */
FieldT::~FieldT(void)
{
	delete fIntegrator;

	for (int i = 0; i < fSourceOutput.Length(); i++)
		delete fSourceOutput[i];
		
	for (int i = 0; i < fKBC_Controllers.Length(); i++)
		delete fKBC_Controllers[i];

	for (int i = 0; i < fFBC_Controllers.Length(); i++)
		delete fFBC_Controllers[i];
}

void FieldT::RegisterLocal(LocalArrayT& array) const	
{
	const char caller[] = "FieldT::RegisterLocal";
	switch (array.Type())
	{
		case LocalArrayT::kDisp:
		{
			array.SetGlobal(fField[0]);
			break;
		}
		case LocalArrayT::kVel:
		{
			if (Order() < 1) ExceptionT::GeneralFail(caller, "only up to order %d: 1", Order());
			array.SetGlobal(fField[1]);
			break;			
		}
		case LocalArrayT::kAcc:
		{
			if (Order() < 2) ExceptionT::GeneralFail(caller, "only up to order %d: 2", Order());
			array.SetGlobal(fField[2]);
			break;			
		}
		case LocalArrayT::kLastDisp:
		{
			array.SetGlobal(fField_last[0]);
			break;
		}
		case LocalArrayT::kLastVel:
		{
			if (Order() < 1) ExceptionT::GeneralFail(caller, "only up to order %d: 1", Order());
			array.SetGlobal(fField_last[1]);
			break;			
		}
		case LocalArrayT::kLastAcc:
		{
			if (Order() < 2) ExceptionT::GeneralFail(caller, "only up to order %d: 2", Order());
			array.SetGlobal(fField_last[2]);
			break;			
		}
		default:
			ExceptionT::GeneralFail(caller, "unrecognized type: %d", array.Type());
	}
}

/* set all field values to 0.0 */
void FieldT::Clear(void)
{
	/* inherited */
	BasicFieldT::Clear();
	
	/* clear */
	for (int i = 0; i < fField_last.Length(); i++)
		fField_last[i] = 0.0;
}

/* set number of nodes */
void FieldT::Dimension(int nnd, bool copy_in)
{
	/* inherited */
	BasicFieldT::Dimension(nnd, copy_in);

	/* set dimensions with integrator */
	nIntegrator().Dimension(*this);
}

/* append the equation sets generated by the KBC_ControllerT's and
 * FBC_ControllerT's. */
void FieldT::EquationSets(AutoArrayT<const iArray2DT*>& eq_1, 
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	/* KBC controllers */
	for (int i = 0; i < fKBC_Controllers.Length(); i++)
		fKBC_Controllers[i]->Equations(eq_1);

	/* FBC controllers */
	for (int i = 0; i < fFBC_Controllers.Length(); i++)
		fFBC_Controllers[i]->Equations(eq_1, eq_2);
}

/* set the time step */
void FieldT::SetTimeStep(double dt) { Integrator().SetTimeStep(dt); }

/* append connectivities */
void FieldT::Connectivities(AutoArrayT<const iArray2DT*>& connects_1,
	AutoArrayT<const RaggedArray2DT<int>*>& connects_2,
	AutoArrayT<const iArray2DT*>& equivalents) const
{
	/* KBC controllers */
	for (int i = 0; i < fKBC_Controllers.Length(); i++)
		fKBC_Controllers[i]->Connectivities(connects_1, equivalents);

	/* FBC controllers */
	for (int i = 0; i < fFBC_Controllers.Length(); i++)
		fFBC_Controllers[i]->Connectivities(connects_1, connects_2, equivalents);
}

/* return the GlobalT::SystemTypeT for the  group */
GlobalT::SystemTypeT FieldT::SystemType(void) const
{
	/* normally no contribution to the stiffness */
	GlobalT::SystemTypeT type = GlobalT::kUndefined;

	/* FBC controllers */
	for (int i = 0; i < fFBC_Controllers.Length(); i++)
		type = GlobalT::MaxPrecedence(type, fFBC_Controllers[i]->TangentType());

	return type;
}

/* beginning of time series */
void FieldT::InitialCondition(void)
{
	/* initial fields */
	for (int i = 0; i < fField.Length(); i++)
		fField[i] = 0.0;

	/* apply initial cards */
	bool all_active = true;
	for (int i = 0; i < fIC.Length(); i++)
		if (!Apply_IC(fIC[i]))
			all_active = false;
	if (!all_active)
		cout << "\n FieldT::InitialCondition: initial conditions applied to prescribed\n" 
		     <<   "     equations in field \"" << FieldName() << "\" are being ignored" << endl;

	/* KBC controllers */
	for (int k = 0; k < fKBC_Controllers.Length(); k++)
		fKBC_Controllers[k]->InitialCondition();

	/* FBC controllers */
	for (int i = 0; i < fFBC_Controllers.Length(); i++)
		fFBC_Controllers[i]->InitialCondition();

	/* apply KBC cards */
	for (int i = 0; i < fKBC.Length(); i++)
		nIntegrator().ConsistentKBC(*this, fKBC[i]);
		
	/* apply KBC cards generated by KBC controllers */
	for (int j = 0; j < fKBC_Controllers.Length(); j++)
	{
		/* cards */
		const ArrayT<KBC_CardT>& cards = fKBC_Controllers[j]->KBC_Cards();
		/* apply KBC cards */
		for (int i = 0; i < cards.Length(); i++)
			nIntegrator().ConsistentKBC(*this, cards[i]);
	}

	/* initial history */
	fField_last = fField;
}

/* apply predictor to all degrees of freedom */
void FieldT::InitStep(void)
{
	/* integrator */
	nIntegratorT& integrator = nIntegrator();

	/* predictor to all DOF's */
	integrator.Predictor(*this);

	/* KBC controllers */
	for (int i = 0; i < fKBC_Controllers.Length(); i++)
		fKBC_Controllers[i]->InitStep();

	/* FBC controllers */
	for (int j = 0; j < fFBC_Controllers.Length(); j++)
		fFBC_Controllers[j]->InitStep();

	/* apply KBC cards */
	for (int i = 0; i < fKBC.Length(); i++)
		integrator.ConsistentKBC(*this, fKBC[i]);
		
	/* apply KBC cards generated by KBC controllers */
	for (int j = 0; j < fKBC_Controllers.Length(); j++)
	{
		/* cards */
		const ArrayT<KBC_CardT>& cards = fKBC_Controllers[j]->KBC_Cards();
		/* apply KBC cards */
		for (int i = 0; i < cards.Length(); i++)
			integrator.ConsistentKBC(*this, cards[i]);
	}
}
#pragma message("Roll up the redundancy after it works")
// apply predictor to all owned degrees of freedom 
void FieldT::InitStep(int fieldstart, int fieldend)
{
	/* integrator */
	nIntegratorT& integrator = nIntegrator();
///////// This here is what I want to do only for the owned nodes
	// predictor to all DOF's owned by this proc 
	integrator.Predictor(*this, fieldstart, fieldend);

	/* KBC controllers */
	for (int i = 0; i < fKBC_Controllers.Length(); i++)
		fKBC_Controllers[i]->InitStep();

	/* FBC controllers */
	for (int j = 0; j < fFBC_Controllers.Length(); j++)
		fFBC_Controllers[j]->InitStep();

	/* apply KBC cards */
	for (int i = 0; i < fKBC.Length(); i++)
		integrator.ConsistentKBC(*this, fKBC[i]);
		
	/* apply KBC cards generated by KBC controllers */
	for (int j = 0; j < fKBC_Controllers.Length(); j++)
	{
		/* cards */
		const ArrayT<KBC_CardT>& cards = fKBC_Controllers[j]->KBC_Cards();
		/* apply KBC cards */
		for (int i = 0; i < cards.Length(); i++)
			integrator.ConsistentKBC(*this, cards[i]);
	}
///////// 
}

/* assemble contributions to the residual */
void FieldT::FormRHS(void)
{
	/* has force bpundary conditions */
	if (fFBC.Length() > 0)
	{
		/* collect nodal forces */
		for (int i = 0; i < fFBC.Length(); i++)
			fFBCValues[i] = fFBC[i].CurrentValue();	
	
		/* assemble */
		fFieldSupport.AssembleRHS(Group(), fFBCValues, fFBCEqnos);
	}

	/* KBC controllers */
	for (int i = 0; i < fKBC_Controllers.Length(); i++)
		fKBC_Controllers[i]->FormRHS();
	
	/* FBC controllers */
	for (int i = 0; i < fFBC_Controllers.Length(); i++)
		fFBC_Controllers[i]->ApplyRHS();
}

/* assemble contributions to the tangent */
void FieldT::FormLHS(GlobalT::SystemTypeT sys_type)
{
	/* KBC controllers */
	for (int i = 0; i < fKBC_Controllers.Length(); i++)
		fKBC_Controllers[i]->FormLHS(sys_type);

	/* FBC controllers */
	for (int i = 0; i < fFBC_Controllers.Length(); i++)
		fFBC_Controllers[i]->ApplyLHS(sys_type);
}

/* update history */
void FieldT::CloseStep(void)
{
	/* update history */
	fField_last = fField;

	/* KBC controllers */
	for (int i = 0; i < fKBC_Controllers.Length(); i++)
		fKBC_Controllers[i]->CloseStep();

	/* FBC controllers */
	for (int i = 0; i < fFBC_Controllers.Length(); i++)
		fFBC_Controllers[i]->CloseStep();
}

/* overwrite the update values in the FieldT::Update array */
void FieldT::AssembleUpdate(const dArrayT& update)
{
	int *peq = fEqnos.Pointer();
	int len = fEqnos.Length();
	double *p = fUpdate.Pointer();
	for (int i = 0; i < len; i++)
	{
		/* local equation */
		int eq = *peq++ - fEquationStart;
		
		/* active dof */
		if (eq > -1 && eq < fNumEquations)
			*p = update[eq];
		else
			*p = 0.0;

		/* next */
		p++;
	}
#pragma message("FieldT -- Needs FBC controllers too?")	
	/* KBC controllers */
	for (int i = 0; i < fKBC_Controllers.Length(); i++)
		fKBC_Controllers[i]->Update(update);

}

/* update the active degrees of freedom */
void FieldT::ApplyUpdate(int fPartFieldStart /*= 0*/, int fPartFieldEnd /*= -1*/)
{
	/* corrector */
	nIntegrator().Corrector(*this, fUpdate,fPartFieldStart, fPartFieldEnd, 0);
}

/* copy nodal information */
void FieldT::CopyNodeToNode(const ArrayT<int>& source, const ArrayT<int>& target)
{
	/* copy data from source nodes to target nodes */
	for (int i = 0 ; i < fField.Length() ; i++ )
	{
		dArray2DT& field = fField[i];
		dArray2DT& field_last = fField_last[i];
		for (int j = 0; j < source.Length(); j++)
		{
			int from = source[j];
			int to = target[j];
			field.CopyRowFromRow(to,from);
			field_last.CopyRowFromRow(to,from);
		}
	}
}

/* check for relaxation */
GlobalT::RelaxCodeT FieldT::RelaxSystem(void)
{
	GlobalT::RelaxCodeT relax = GlobalT::kNoRelax;

	/* kinematics BC controllers */
	for (int i = 0; i < fKBC_Controllers.Length(); i++)
	{
		/* check controller */
		GlobalT::RelaxCodeT code = fKBC_Controllers[i]->RelaxSystem();

//NOTE - do the boundary condition really need to be reapplied here?

		/* reset BC */
		if (code != GlobalT::kNoRelax)
		{
			/* boundary condition cards generated by the controller */
			const ArrayT<KBC_CardT>& cards = fKBC_Controllers[i]->KBC_Cards();

			/* apply BC's */
			nIntegratorT& integrator = nIntegrator();
			for (int j = 0; j < cards.Length(); j++)
				integrator.ConsistentKBC(*this, cards[j]);
		}
		relax = GlobalT::MaxPrecedence(relax, code);
	}
	
	/* force BC controllers */
	for (int i = 0; i < fFBC_Controllers.Length(); i++)
		relax = GlobalT::MaxPrecedence(relax, fFBC_Controllers[i]->RelaxSystem());
	
	return relax;
}

/* reset displacements (and configuration to the last known solution) */
GlobalT::RelaxCodeT FieldT::ResetStep(void)
{
	/* reset field */
	fField = fField_last;

	/* KBC controllers */
	for (int j = 0; j < fKBC_Controllers.Length(); j++)
		fKBC_Controllers[j]->Reset();

	/* FBC controllers */
	for (int i = 0; i < fFBC_Controllers.Length(); i++)
		fFBC_Controllers[i]->Reset();
		
	return GlobalT::kNoRelax;
}

/* mark prescribed equations */
void FieldT::InitEquations(void)
{
	/* use the allocated space */
	fEqnos.Dimension(NumNodes(), NumDOF());
	fEqnos = FieldT::kInit;
	
	/* mark KBC nodes */
	for (int j = 0; j < fKBC.Length(); j++)
		SetBCCode(fKBC[j]);
	
	/* mark nodes used by KBC controllers */
	for (int j = 0; j < fKBC_Controllers.Length(); j++)
	{
		/* boundary condition cards generated by the controller */
		const ArrayT<KBC_CardT>& cards = fKBC_Controllers[j]->KBC_Cards();

		if (!fKBC_Controllers[j]->IsICController())
		{
			/* mark global equation numbers */
			for (int i = 0; i < cards.Length(); i++)
				SetBCCode(cards[i]);
		}
	}
}

/* set the equations array and the number of active equations */
void FieldT::FinalizeEquations(int eq_start, int num_eq)
{
	/* store parameters */
	fEquationStart = eq_start;
	fNumEquations = num_eq;

	/* set force boundary condition destinations */
	SetFBCEquations();

	/* run through KBC controllers */
	for (int j = 0; j < fKBC_Controllers.Length(); j++)
	  fKBC_Controllers[j]->SetEquations();
}

/* Collect the local element lists */
void FieldT::SetLocalEqnos(const iArray2DT& nodes, iArray2DT& eqnos) const
{
/* consistency checks */
#if __option (extended_errorcheck)
	if (nodes.MajorDim() != eqnos.MajorDim() ||
	    eqnos.MinorDim() < nodes.MinorDim()*NumDOF()) ExceptionT::GeneralFail("FieldT::SetLocalEqnos");
		//must have enough space (and maybe more)
#endif

	int numel = nodes.MajorDim();
	int nen   = nodes.MinorDim();
	int ndof  = NumDOF();
	for (int i = 0; i < numel; i++)
	{
		const int* pnodes = nodes(i);
		int* pien   = eqnos(i);
		for (int j = 0; j < nen; j++)
		{
			int nodenum = *pnodes++;
			for (int k = 0; k < ndof; k++)
				*pien++ = fEqnos(nodenum, k);
		}
	}
}

/* collect equation numbers */
void FieldT::SetLocalEqnos(ArrayT<const iArray2DT*> nodes, iArray2DT& eqnos) const
{
	const char caller[] = "FieldT::SetLocalEqnos";
	
	int row = 0;
	for (int i = 0; i < nodes.Length(); i++)
	{
		const iArray2DT& nd = *(nodes[i]);
	
		/* check */
		if (row + nd.MajorDim() > eqnos.MajorDim()) ExceptionT::OutOfRange(caller);

		/* single block */	
		iArray2DT eq(nd.MajorDim(), eqnos.MinorDim(), eqnos(row));
		SetLocalEqnos(nd, eq);
	
		/* next */
		row += nd.MajorDim();
	}
	
	/* check - must fill all of eqnos */
	if (row != eqnos.MajorDim()) ExceptionT::SizeMismatch(caller);
}

void FieldT::SetLocalEqnos(const RaggedArray2DT<int>& nodes,
	RaggedArray2DT<int>& eqnos, const iArrayT& which_dofs) const
{
/* consistency checks */
#if __option(extended_errorcheck)
	const char caller[] = "FieldT::SetLocalEqnos";
	if (nodes.MajorDim() != eqnos.MajorDim()) ExceptionT::SizeMismatch(caller);
	if (which_dofs.Length() != nodes.MajorDim()) ExceptionT::SizeMismatch(caller);
#endif
	
	int numel = nodes.MajorDim();
	for (int i = 0; i < numel; i++)
	{
#if __option(extended_errorcheck)
		if (eqnos.MinorDim(i) != nodes.MinorDim(i)) ExceptionT::SizeMismatch(caller);
		//must have enough space (and maybe more)
#endif
		int  nen    = nodes.MinorDim(i);
		const int* pnodes = nodes(i);
		int* pien   = eqnos(i);
		int k = which_dofs[i];
		for (int j = 0; j < nen; j++)
			*pien++ = fEqnos(*pnodes++, k);
	}
}

void FieldT::SetLocalEqnos(const RaggedArray2DT<int>& nodes,
	RaggedArray2DT<int>& eqnos) const
{
/* consistency checks */
#if __option(extended_errorcheck)
	const char caller[] = "FieldT::SetLocalEqnos";
	if (nodes.MajorDim() != eqnos.MajorDim()) ExceptionT::SizeMismatch(caller);
#endif
	
	int numel = nodes.MajorDim();
	for (int i = 0; i < numel; i++)
	{
#if __option(extended_errorcheck)
	   	if (eqnos.MinorDim(i) < nodes.MinorDim(i)*NumDOF()) ExceptionT::SizeMismatch(caller);
		//must have enough space (and maybe more)
#endif
		int  nen    = nodes.MinorDim(i);
		const int* pnodes = nodes(i);
		int* pien   = eqnos(i);
		int ndof    = NumDOF();
		
		for (int j = 0; j < nen; j++)
		{
			int nodenum = *pnodes++;
			for (int k = 0; k < ndof; k++)
				*pien++ = fEqnos(nodenum, k);
		}
	}
}

/* NB that nodes is not declared const since traversing the linked list
 * modifies pointers. 
 */
void FieldT::SetLocalEqnos(ArrayT< LinkedListT<int> >& nodes,
	RaggedArray2DT<int>& eqnos) const
{
/* consistency checks */
#if __option(extended_errorcheck)
	const char caller[] = "FieldT::SetLocalEqnos";
	if (nodes.Length() != eqnos.MajorDim()) ExceptionT::SizeMismatch(caller);
#endif
	
	int numel = nodes.Length();
	int ndof    = NumDOF();
	for (int i = 0; i < numel; i++)
	{
#if __option(extended_errorcheck)
		if (eqnos.MinorDim(i) < nodes[i].Length()*ndof) ExceptionT::SizeMismatch(caller);
		//must have enough space (and maybe more)
#endif
		int  nen    = eqnos.MinorDim(i)/ndof; // do this to avoid traversing the list twice
		int pnodes;
		int* pien   = eqnos(i);
		LinkedListT<int>& nodeList = nodes[i];
		nodeList.Top();
		while (nodeList.Next(pnodes))
		{
			for (int k = 0; k < ndof; k++)
				*pien++ = fEqnos(pnodes, k);
		}
	}
}

void FieldT::ReadRestart(ifstreamT& in, const ArrayT<int>* nodes)
{
	/* external file */
	StringT file;
	ifstreamT u_in;

	/* read fields */
	StringT deriv;
	for (int i = 0; i < fField.Length(); i++)
	{
		/* file name */
		StringT file = in.filename();
		file.Append(".", FieldName());
		file.Append(".", deriv, "u");

		/* read */
		u_in.open(file);
		if (u_in.is_open()) 
		{
			if (nodes) /* select nodes */
			{
				dArray2DT tmp(nodes->Length(), NumDOF());
				u_in >> tmp;
				fField[i].Assemble(*nodes, tmp);
			}
			else /* all nodes */
				u_in >> fField[i];

			u_in.close();
		} 
		else 
		{
			cout << "\n FieldT::ReadRestart: file not found: " 
                 << file << '\n' <<   "     assuming order " << i << " = 0.0" << endl;
			fField[i] = 0.0;
		}

		/* next */
		deriv.Append("D");
	}

	/* KBC controllers */
	for (int j = 0; j < fKBC_Controllers.Length(); j++)
		fKBC_Controllers[j]->ReadRestart(in);

	/* FBC controllers */
	for (int i = 0; i < fFBC_Controllers.Length(); i++)
		fFBC_Controllers[i]->ReadRestart(in);
}

void FieldT::WriteRestart(ofstreamT& out, const ArrayT<int>* nodes) const
{
	/* external file */
	ofstreamT u_out;
	u_out.precision(out.precision());

	/* write field */
	StringT deriv;
	for (int i = 0; i < fField.Length(); i++)
	{
		/* file name */
		StringT file = out.filename();
		file.Append(".", FieldName());
		file.Append(".", deriv, "u");

		/* write */
		u_out.open(file);
		if (nodes) /* select nodes */
		{
			dArray2DT tmp(nodes->Length(), NumDOF());
			tmp.RowCollect(*nodes, fField[i]);
			u_out << tmp << '\n';
		}
		else /* all nodes */
			u_out << fField[i] << '\n';
		u_out.close();
	
		/* next */
		deriv.Append("D");
	}

	/* KBC controllers */
	for (int j = 0; j < fKBC_Controllers.Length(); j++)
		fKBC_Controllers[j]->WriteRestart(out);

	/* FBC controllers */
	for (int i = 0; i < fFBC_Controllers.Length(); i++)
		fFBC_Controllers[i]->WriteRestart(out);
}

/* register results for output */
void FieldT::RegisterOutput(void)
{
	/* FBC controllers */
	for (int i = 0; i < fFBC_Controllers.Length(); i++)
		fFBC_Controllers[i]->RegisterOutput();
}

/* write output data */
void FieldT::WriteOutput(ostream& out) const
{
	/* KBC controllers */
	for (int j = 0; j < fKBC_Controllers.Length(); j++)
		fKBC_Controllers[j]->WriteOutput(out);

	/* FBC controllers */
	for (int i = 0; i < fFBC_Controllers.Length(); i++)
		fFBC_Controllers[i]->WriteOutput(out);
}

/* write field parameters to output stream */
void FieldT::WriteParameters(ostream& out) const
{
	out << "\n F i e l d : \"" << FieldName() << "\"\n\n";
	out << " Number of degrees of freedom. . . . . . . . . . = " << NumDOF() << '\n';
	for (int i = 0; i < fLabels.Length(); i++)
		out << '\t' << fLabels[i] << '\n';
	out << " Number of time derivatives. . . . . . . . . . . = " << Order() << '\n';
	out << " Group number. . . . . . . . . . . . . . . . . . = " << Group() << '\n';
	out.flush();
}

/* accumulate source terms */
void FieldT::RegisterSource(const StringT& ID, const dArray2DT& source) const
{
	/* search */
	int dex = SourceIndex(ID);
	
	/* NOTE: need this function to be const else ElementBaseT cannot call it */
	FieldT* non_const_this = const_cast<FieldT*>(this);

	/* add new */
	if (dex == -1) 
	{
		/* add ID to list */
		non_const_this->fID.Append(ID);
	
		/* source */
		dArray2DT* new_source = new dArray2DT;
		non_const_this->fSourceOutput.Append(new_source);
	}

	/* register sources */
	non_const_this->fSourceID.Append(ID);
	non_const_this->fSourceBlocks.Append(&source);
}

/* element source terms */
const dArray2DT* FieldT::Source(const StringT& ID) const
{
	/* search */
	int dex = SourceIndex(ID);

	/* accmulate and return */
	if (dex > -1) {
	
		/* return value */
		dArray2DT& source = *(fSourceOutput[dex]);
	
		/* accumulate */
		int count = 0;
		for (int i = 0; i < fSourceID.Length(); i++)
			if (fSourceID[i] == ID)
			{
				/* implied that all sources have the same dimension */
				if (count == 0)
					source = *(fSourceBlocks[i]);
				else
					source += *(fSourceBlocks[i]);
				count++;
			}
			
		return &source;
	}
	else
		return NULL;
}

/* describe the parameters needed by the interface */
void FieldT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	/* field name */
	list.AddParameter(ParameterT::Word, "field_name");

	/* solution group */
	ParameterT solver_group(ParameterT::Integer, "solution_group");
	solver_group.AddLimit(1, LimitT::LowerInclusive);
	solver_group.SetDefault(1);
	list.AddParameter(solver_group);
	
	/* integrator number */
	ParameterT integrator(ParameterT::Enumeration, "integrator");
	integrator.AddEnumeration(     "linear_static", IntegratorT::kLinearStatic);
	integrator.AddEnumeration(            "static", IntegratorT::kStatic);
	integrator.AddEnumeration(         "trapezoid", IntegratorT::kTrapezoid);
	integrator.AddEnumeration(        "linear_HHT", IntegratorT::kLinearHHT);
	integrator.AddEnumeration(     "nonlinear_HHT", IntegratorT::kNonlinearHHT);
	integrator.AddEnumeration("central_difference", IntegratorT::kExplicitCD);
	integrator.AddEnumeration(            "Verlet", IntegratorT::kVerlet);
	integrator.AddEnumeration(             "Gear6", IntegratorT::kGear6);
	integrator.SetDefault(IntegratorT::kStatic);
	list.AddParameter(integrator);
}

/* information about subordinate parameter lists */
void FieldT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* degrees of freedom - prescribe number or labels */
	sub_list.AddSub("ndof_specification", ParameterListT::Once, true);

	/* initial conditions */
	sub_list.AddSub("initial_condition", ParameterListT::Any);

	/* kinematic boundary conditions */
	sub_list.AddSub("kinematic_BC", ParameterListT::Any);
	
	/* force boundary conditions */
	sub_list.AddSub("force_BC", ParameterListT::Any);
	
	/* KBC controllers */
	sub_list.AddSub("KBC_controllers", ParameterListT::Any, true);
	
	/* FBC controllers */
	sub_list.AddSub("FBC_controllers", ParameterListT::Any, true);
}

/* return the description of the given inline subordinate parameter list */
void FieldT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "ndof_specification")
	{
		order = ParameterListT::Choice;

		/* just give an integer */
		sub_lists.AddSub("dof_count");		

		/* provide a list labels */
		sub_lists.AddSub("dof_labels");		
	}
	else if (name == "KBC_controllers")
	{
		order = ParameterListT::Choice;
		
		/* choices - KBC_ControllerT::Code must translate names */
		sub_lists.AddSub("K-field");	
		sub_lists.AddSub("bi-material_K-field");	
		sub_lists.AddSub("torsion");	
		sub_lists.AddSub("mapped_nodes");
		sub_lists.AddSub("scaled_velocity");
		sub_lists.AddSub("tied_nodes");
		sub_lists.AddSub("periodic_nodes");
		sub_lists.AddSub("conveyor");
	}
	else if (name == "FBC_controllers")
	{
		order = ParameterListT::Choice;
		
		/* choices - FBC_controllers::Code must translate names */
		sub_lists.AddSub("sphere_penalty");
		sub_lists.AddSub("sphere_augmented_Lagrangian");
		sub_lists.AddSub("sphere_penalty_meshfree");
		sub_lists.AddSub("wall_penalty");
		sub_lists.AddSub("wall_augmented_Lagrangian");
		sub_lists.AddSub("cylinder_penalty");
		sub_lists.AddSub("cylinder_augmented_Lagrangian");
		sub_lists.AddSub("augmented_Lagrangian_KBC_meshfree");
	}
	else /* inherited */
		ParameterInterfaceT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FieldT::NewSub(const StringT& name) const
{
	/* non-const this */
	FieldT* non_const_this = (FieldT*) this;

	/* (try to) translate to KBC code */
	KBC_ControllerT::CodeT KBC_code = KBC_ControllerT::Code(name);
	if (KBC_code != KBC_ControllerT::kNone)
		return fFieldSupport.NewKBC_Controller(*non_const_this, KBC_code);

	/* (try to) translate to FBC code */
	FBC_ControllerT::CodeT FBC_code = FBC_ControllerT::Code(name);
	if (FBC_code != FBC_ControllerT::kNone)
		return fFieldSupport.NewFBC_Controller(FBC_code);

	if (name == "dof_count")
	{
		ParameterContainerT* dof_count = new ParameterContainerT(name);

		ParameterT ndof(ParameterT::Integer, "ndof");
		ndof.SetDescription("number of unknown per node");
		ndof.AddLimit(0, LimitT::LowerInclusive);
		dof_count->AddParameter(ndof);

		return dof_count;
	}
	else if (name == "dof_labels")
	{
		StringListT* dof_labels = new StringListT("dof_labels");
		dof_labels->SetMinLength(1);
		return dof_labels;
	}
	else if (name == "initial_condition")
	{
		ParameterContainerT* ic = new ParameterContainerT(name);

		/* description */
		ic->SetDescription("apply to node set or all");

		/* define as node set or all */		
		ic->AddParameter(ParameterT::Word, "node_ID", ParameterListT::ZeroOrOnce);
		ParameterT all_nodes(ParameterT::Boolean, "all_nodes");
		all_nodes.SetDefault(true);
		ic->AddParameter(all_nodes, ParameterListT::ZeroOrOnce);

		ic->AddParameter(ParameterT::Integer, "dof");
		ParameterT IC_type(ParameterT::Enumeration, "type");
		IC_type.AddEnumeration("u", 0);
		IC_type.AddEnumeration("D_u", 1);
		IC_type.AddEnumeration("DD_u", 2);
		IC_type.AddEnumeration("D3_u", 3);
		IC_type.AddEnumeration("D4_u", 4);
		IC_type.SetDefault(0);
		ic->AddParameter(IC_type);
		ParameterT value(ParameterT::Double, "value");
		value.SetDefault(0.0);
		ic->AddParameter(value);
	
		return ic;	
	}
	else if (name == "kinematic_BC")
	{
		ParameterContainerT* kbc = new ParameterContainerT(name);
		
		kbc->AddParameter(ParameterT::Word, "node_ID");
		kbc->AddParameter(ParameterT::Integer, "dof");
		ParameterT BC_type(ParameterT::Enumeration, "type");
		BC_type.AddEnumeration("fixed", -1);
		BC_type.AddEnumeration("u", 0);
		BC_type.AddEnumeration("D_u", 1);
		BC_type.AddEnumeration("DD_u", 2);
		BC_type.AddEnumeration("D3_u", 3);
		BC_type.AddEnumeration("D4_u", 4);
		BC_type.SetDefault(-1);
		kbc->AddParameter(BC_type);
		ParameterT schedule(ParameterT::Integer, "schedule");
		schedule.SetDefault(0);
		kbc->AddParameter(schedule);
		ParameterT value(ParameterT::Double, "value");
		value.SetDefault(0.0);
		kbc->AddParameter(value);
	
		return kbc;
	}
	else if (name == "force_BC")
	{
		ParameterContainerT* fbc = new ParameterContainerT(name);
		
		fbc->AddParameter(ParameterT::Word, "node_ID");
		fbc->AddParameter(ParameterT::Integer, "dof");
		ParameterT schedule(ParameterT::Integer, "schedule");
		schedule.SetDefault(0);
		fbc->AddParameter(schedule);
		ParameterT value(ParameterT::Double, "value");
		value.SetDefault(0.0);
		fbc->AddParameter(value);
	
		return fbc;	
	}
	/* inherited */
	return ParameterInterfaceT::NewSub(name);
}

/* accept parameter list */
void FieldT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "FieldT::TakeParameterList";

	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* field name */
	const StringT& field_name = list.GetParameter("field_name");

	/* number of degrees of freedom per node */
	const ParameterListT& ndof_spec = list.GetListChoice(*this, "ndof_specification");
	ArrayT<StringT> labels;
	int ndof = 0;
	if (ndof_spec.Name() == "dof_count")
		ndof = ndof_spec.GetParameter("ndof");
	else if (ndof_spec.Name() == "dof_labels")
	{
		/* labels */
		const ArrayT<ParameterListT>& dof_labels = ndof_spec.Lists();
		ndof = dof_labels.Length();
		labels.Dimension(ndof);
		for (int i = 0; i < ndof; i++)
			labels[i] = dof_labels[i].GetParameter("value");
	}
	else
		ExceptionT::GeneralFail(caller, "not expecting \"%s\" for \"ndof_specification\" in \"%s\"",
			ndof_spec.Name().Pointer(), list.Name().Pointer());

	/* solution group */
	fGroup = list.GetParameter("solution_group");
	fGroup--;

	/* construct integrator */
	int integrator_type = list.GetParameter("integrator");
	fIntegrator = IntegratorT::New(integrator_type, true);
	
	/* cast to nodal interface */
	fnIntegrator = &(fIntegrator->nIntegrator());
	if (!fnIntegrator) ExceptionT::GeneralFail(caller);

	/* configure the field */
	int order = fIntegrator->Order();
	Initialize(field_name, ndof, order);
	if (labels.Length() > 0) SetLabels(labels);
	
	/* construct controllers and count numbers of IC, KBC, and FBC */
	int num_IC = 0;
	int num_KBC = 0;
	int num_FBC = 0;
	const ArrayT<ParameterListT>& subs = list.Lists();
	for (int i = 0; i < subs.Length(); i++)
	{
		const StringT& name = subs[i].Name();
		bool resolved = false;
	
		/* try KBC */
		if (!resolved) {
			KBC_ControllerT::CodeT KBC_code = KBC_ControllerT::Code(name);
			if (KBC_code != KBC_ControllerT::kNone) {
				KBC_ControllerT* KBC_controller = fFieldSupport.NewKBC_Controller(*this, KBC_code);
				
				/* store */
				AddKBCController(KBC_controller);

				/* initialize */
				KBC_controller->TakeParameterList(subs[i]);
				resolved = true;
			}
		}
			
		/* try FBC */
		if (!resolved) {
			FBC_ControllerT::CodeT FBC_code = FBC_ControllerT::Code(name);
			if (FBC_code != FBC_ControllerT::kNone) {
				FBC_ControllerT* FBC_controller = fFieldSupport.NewFBC_Controller(FBC_code);

				/* store */
				AddFBCController(FBC_controller);

				/* initialize */
				FBC_controller->SetField(*this);
				FBC_controller->TakeParameterList(subs[i]);
				resolved = true;
			}
		}
		
		/* look for other sublists */
		if (!resolved) {
			if (name == "initial_condition") {			

				/* node set or all nodes */				
				const ParameterT* node_ID = subs[i].Parameter("node_ID");
				const ParameterT* all_nodes = subs[i].Parameter("all_nodes");
				if ((node_ID && all_nodes) || (!node_ID && !all_nodes))
					ExceptionT::GeneralFail(caller, "expecting either \"node_ID\" or \"all_nodes\" in \"initial_condition\"");
				else if (node_ID) /* get size of the node set */ {
					const StringT& ID = *node_ID;
					ModelManagerT& model_manager = fFieldSupport.ModelManager();
					num_IC += model_manager.NodeSet(ID).Length();
				}
				else /* single card for all nodes */ {
					bool all = *all_nodes;
					if (all) num_IC++;
				}
			}
			else if (name == "kinematic_BC") {
				const StringT& node_ID = subs[i].GetParameter("node_ID");
				ModelManagerT& model_manager = fFieldSupport.ModelManager();
				num_KBC += model_manager.NodeSet(node_ID).Length();
			}
			else if (name == "force_BC") {
				const StringT& node_ID = subs[i].GetParameter("node_ID");
				ModelManagerT& model_manager = fFieldSupport.ModelManager();
				num_FBC += model_manager.NodeSet(node_ID).Length();
			}
		}
	}

	/* construct controllers and count numbers of IC, KBC, and FBC */
	if (num_IC > 0 || num_KBC > 0 || num_FBC > 0)
	{
		ModelManagerT& model_manager = fFieldSupport.ModelManager();
		fIC.Dimension(num_IC);
		fKBC.Dimension(num_KBC);
		fFBC.Dimension(num_FBC);
		num_IC = num_KBC = num_FBC = 0;
		for (int i = 0; i < subs.Length(); i++)
		{
			const ParameterListT& sub = subs[i];
			const StringT& name = sub.Name();
			if (name == "initial_condition") {

				/* extract values */
				int dof = sub.GetParameter("dof"); dof--;
				int order = sub.GetParameter("type");
				double value = sub.GetParameter("value");

				/* look for node set ID - exclusive choice checked above */
				const ParameterT* node_ID = sub.Parameter("node_ID");
				const ParameterT* all_nodes = sub.Parameter("all_nodes");
				if (node_ID) /* set cards */ {
					const StringT& ID = *node_ID;
					const iArrayT& set = model_manager.NodeSet(ID);
					for (int i = 0; i < set.Length(); i++)
						fIC[num_IC++].SetValues(set[i], dof, order, value);
				}
				else /* all nodes */ {
					bool all = *all_nodes;
					if (all) fIC[num_IC++].SetValues(-1, dof, order, value); /* -1 => "all" */
				}
			}
			else if (name == "kinematic_BC") {

				/* extract values */
				const StringT& node_ID = sub.GetParameter("node_ID");
				int dof = sub.GetParameter("dof"); dof--;
				int typ = sub.GetParameter("type");
				KBC_CardT::CodeT code = KBC_CardT::int2CodeT(typ + 1);
				int schedule_no = sub.GetParameter("schedule"); schedule_no--;
				double value = sub.GetParameter("value");

				/* get the schedule */
				const ScheduleT* schedule = (schedule_no > -1) ? fFieldSupport.Schedule(schedule_no) : NULL;

				/* set cards */
				const iArrayT& set = model_manager.NodeSet(node_ID);
				for (int i = 0; i < set.Length(); i++)
					fKBC[num_KBC++].SetValues(set[i], dof, code, schedule, value);
			}
			else if (name == "force_BC") {
		
				/* extract values */
				const StringT& node_ID = sub.GetParameter("node_ID");
				int node = atoi(node_ID);
				int dof = sub.GetParameter("dof"); dof--;
				int schedule_no = sub.GetParameter("schedule"); schedule_no--;
				double value = sub.GetParameter("value");

				/* get the schedule */
				const ScheduleT* schedule = (schedule_no > -1) ? fFieldSupport.Schedule(schedule_no) : NULL;
		
				/* set card */
				const NodeManagerT& node_man = fFieldSupport.NodeManager();
				const iArrayT& set = model_manager.NodeSet(node_ID);
				for (int i = 0; i < set.Length(); i++)
					fFBC[num_FBC++].SetValues(set[i], dof, schedule, value);
			}
		}
	}
}

/**********************************************************************
 * Private
 **********************************************************************/

/* return the index for the source of the given ID */
int FieldT::SourceIndex(const StringT& ID) const
{
	/* search */
	for (int i = 0; i < fID.Length(); i++)
		if (fID[i] == ID)
			return i;

	/* fail */
	return -1;
}

bool FieldT::Apply_IC(const IC_CardT& card)
{
	const char caller[] = "FieldT::Apply_IC";

	/* check order - do not allow initial conditions on highest order derivative
	 * of the field */
	if (card.Order() > Order()-1)
		ExceptionT::OutOfRange(caller, "order is out of range {0,%d}: %d", 
			Order()-1, card.Order());

	/* decode */
	int node     = card.Node();
	int dof      = card.DOF();
	double value = card.Value();

	/* set */
	bool active = true;
	dArray2DT& field = fField[card.Order()];
	if (node == -1)
	{
		/* must all be free */
		iArrayT tmp(fEqnos.MajorDim());
		fEqnos.ColumnCopy(dof,tmp);
		if (tmp.Min() < 1) active = false;

		/* set value */
		field.SetColumn(dof, value);
	}
	else
	{
		/* must be free DOF */
		if (node != -1 && EquationNumber(node,dof) < 1) active = false;

		/* set value */
		field(node, dof) = value;		
	}
	
	return active;
}

/* mark global equations with the specified BC */
void FieldT::SetBCCode(const KBC_CardT& card)
{
	/* mark equation */
	fEqnos(card.Node(), card.DOF()) = kPrescribed;
}

/* resolve destinations for force BC's */
void FieldT::SetFBCEquations(void)
{
	/* allocate */
	int nbc = fFBC.Length();
	fFBCValues.Dimension(nbc);
	fFBCEqnos.Dimension(nbc);

	/* force cards */
	int bad_count = 0;
	for (int i = 0; i < nbc; i++)
	{
		int nodenum, dofnum;
		fFBC[i].Destination(nodenum, dofnum);
		fFBCEqnos[i] = fEqnos(nodenum, dofnum);

		/* component is prescribed */
		if (fFBCEqnos[i] == kPrescribed) bad_count++;
	}
	
	/* warn */
	if (bad_count > 0)
		cout << "\n FieldT::SetFBCEquations: WARNING: found "
		     << bad_count << " prescribed forces\n"
		     <<   "     on equations with prescribed values" << endl;
}
