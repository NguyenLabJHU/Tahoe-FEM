/* $Id: LocalArrayT.cpp,v 1.19 2005-01-26 19:53:27 paklein Exp $ */
/* created: paklein (07/10/1996) */
#include "LocalArrayT.h"
#include "dArray2DT.h"
#include "iArrayT.h"

using namespace Tahoe;

namespace Tahoe {
/* array behavior */
DEFINE_TEMPLATE_STATIC const bool ArrayT<LocalArrayT::TypeT>::fByteCopy = true;
} /* namespace Tahoe */

/* cconstructors */
LocalArrayT::LocalArrayT(void):
	fType(kUnspecified),
	fNumNodes(0),
	fMinorDim(0),
	fGlobal(NULL)
{

}

LocalArrayT::LocalArrayT(TypeT type):
	fType(type),
	fNumNodes(0),
	fMinorDim(0),
	fGlobal(NULL)
{

}

LocalArrayT::LocalArrayT(TypeT type, int numnodes, int minordim):
	fType(type),
	fGlobal(NULL)
{
	Dimension(numnodes, minordim);
}

LocalArrayT::LocalArrayT(const LocalArrayT& source):
	fType(source.fType)
{
	*this = source;
}

/* copy data from an nArrayT */
void LocalArrayT::Copy(int numnodes, int minordim, const nArrayT<double>& source)
{
	const char caller[] = "LocalArrayT::Copy";
#if __option(extended_errorcheck)
	if (numnodes*minordim != source.Length()) ExceptionT::SizeMismatch(caller);
#endif

	/* dimensions */
	fNumNodes = numnodes;
	fMinorDim = minordim;
	if (fGlobal && fMinorDim != fGlobal->MinorDim()) ExceptionT::SizeMismatch(caller);

	/* inherited */
	nArrayT<double>::operator=(source);
}

/* assignment operator */
LocalArrayT& LocalArrayT::operator=(const LocalArrayT& RHS)
{
	/* inherited */
	dArrayT::operator=(RHS);

	/* set dimensions */
	fNumNodes = RHS.fNumNodes;
	fMinorDim = RHS.fMinorDim;

	/* source */
	fGlobal = RHS.fGlobal;
	
	return *this;
}

/* combining arrays - inserts all of source at start_node */
void LocalArrayT::BlockCopyAt(const LocalArrayT& source, int start_node)
{
#if __option (extended_errorcheck)
	const char caller[] = "LocalArrayT::BlockCopyAt";
	if (source.MinorDim() != MinorDim()) ExceptionT::SizeMismatch(caller);
	if (start_node < 0 ||
	    start_node + source.NumberOfNodes() > NumberOfNodes())
	    ExceptionT::OutOfRange(caller);
#endif

	int size = sizeof(double)*source.NumberOfNodes();
	for (int i = 0; i < fMinorDim; i++)
		memcpy((*this)(i) + start_node, source(i), size);
}

/* collect subset */
void LocalArrayT::Collect(const ArrayT<int>& nodes, const LocalArrayT& source)
{
#if __option(extended_errorcheck)
	if (nodes.Length() != fNumNodes || source.fMinorDim != fMinorDim)
		ExceptionT::SizeMismatch("LocalArrayT::Collect");
#endif

	for (int j = 0; j < fMinorDim; j++)
		for (int i = 0; i < nodes.Length(); i++)
			(*this)(i,j) = source(nodes[i],j); 
}

/* compute the array average value */
void LocalArrayT::Average(dArrayT& avg) const
{
	/* dimension */
	avg.Dimension(MinorDim());
	avg = 0.0;
	for (int i = 0; i < MinorDim(); i++)
	{
		const double* p = (*this)(i);
		double& s = avg[i];
		for (int j = 0; j < NumberOfNodes(); j++)
			s += *p++;
	}
	avg /= NumberOfNodes();
}

/* return the vector with transposed indexing */
void LocalArrayT::ReturnTranspose(nArrayT<double>& transpose) const
{
#if __option (extended_errorcheck)
	/* dimensions check */
	if(fLength != transpose.Length()) ExceptionT::SizeMismatch("LocalArrayT::ReturnTranspose");
#endif

	double* ptrans = transpose.Pointer();
	for (int i = 0; i < fNumNodes; i++)
	{
		double* pthis = fArray + i;
		for (int j = 0; j < fMinorDim; j++)
		{
			*ptrans++ = *pthis;
			pthis += fNumNodes;
		}
	}
}

void LocalArrayT::FromTranspose(const double* transpose)
{
	for (int i = 0; i < fNumNodes; i++)
	{
		double* pthis = fArray + i;	
		for (int j = 0; j < fMinorDim; j++)
		{
			*pthis = *transpose++;
			pthis += fNumNodes;
		}
	}
}

void LocalArrayT::AddScaledTranspose(double scale, const nArrayT<double>& transpose)
{
#if __option (extended_errorcheck)
	/* dimension check */
	if(fLength != transpose.Length()) ExceptionT::SizeMismatch("LocalArrayT::AddScaledTranspose");
#endif

	const double* ptrans = transpose.Pointer();
	for (int i = 0; i < fNumNodes; i++)
	{
		double* pthis = fArray + i;	
		for (int j = 0; j < fMinorDim; j++)
		{
			*pthis += scale*(*ptrans++);
			pthis += fNumNodes;
		}
	}
}

/* for registered arrays - preset source for SetLocal */
void LocalArrayT::SetGlobal(const dArray2DT& global)
{
#if __option(extended_errorcheck)
	if (global.MinorDim() != fMinorDim) ExceptionT::SizeMismatch("LocalArrayT::SetGlobal");
#endif
	
	fGlobal = &global;
}

void LocalArrayT::SetLocal(const ArrayT<int>& keys)
{
#if __option (extended_errorcheck)
	const char caller[] = "LocalArrayT::SetLocal";
	if (!fGlobal) ExceptionT::GeneralFail(caller);
	if (keys.Length() != fNumNodes) ExceptionT::SizeMismatch(caller);
#endif

	fGlobal->SetLocal(keys,*this);
}
