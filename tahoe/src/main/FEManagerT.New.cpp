/* $Id: FEManagerT.New.cpp,v 1.4 2005-04-06 15:43:39 paklein Exp $ */
#include "FEManagerT.h"

/* element configuration header */
#include "ElementsConfig.h"

#ifdef BRIDGING_ELEMENT
#include "MultiManagerT.h"
#include "FEManagerT_bridging.h"
#ifdef __DEVELOPMENT__
#include "DevelopmentElementsConfig.h"
#include "BridgingScaleManagerT.h"
#include "FEManagerT_THK.h"
#include "ThermomechanicalCouplingManagerT.h"
#endif
#endif

using namespace Tahoe;

/* factory method */
FEManagerT* FEManagerT::New(const StringT& name, const StringT& input_file, ofstreamT& output, 
	CommunicatorT& comm, const ArrayT<StringT>& argv, TaskT task)
{
	const char caller[] = "FEManagerT::New";

	if (name == "tahoe")
		return new FEManagerT(input_file, output, comm, argv, task);
	else if (name == "tahoe_bridging")
	{
#ifdef BRIDGING_ELEMENT
		return new FEManagerT_bridging(input_file, output, comm, argv, task);
#else
		ExceptionT::GeneralFail(caller, "\"%s\" requires BRIDGING_ELEMENT",
			name.Pointer());
		return NULL;
#endif
	}
	else if (name == "tahoe_multi")
	{
#ifdef BRIDGING_ELEMENT
		return new MultiManagerT(input_file, output, comm, argv, task);
#else
		ExceptionT::GeneralFail(caller, "\"%s\" requires BRIDGING_ELEMENT",
			name.Pointer());
		return NULL;
#endif
	}
	else if (name == "tahoe_bridging_scale")
	{
#if defined(BRIDGING_ELEMENT) && defined(BRIDGING_ELEMENT_DEV)
		return new BridgingScaleManagerT(input_file, output, comm, argv, task);
#else
		ExceptionT::GeneralFail(caller, "\"%s\" requires BRIDGING_ELEMENT and __DEVELOPMENT__",
			name.Pointer());
		return NULL;
#endif
	}
	else if (name == "tahoe_thermomechanical_coupling")
	{
#if defined(BRIDGING_ELEMENT) && defined(BRIDGING_ELEMENT_DEV)
		return new ThermomechanicalCouplingManagerT(input_file, output, comm, argv, task);
#else
		ExceptionT::GeneralFail(caller, "\"%s\" requires BRIDGING_ELEMENT and __DEVELOPMENT__",
			name.Pointer());
		return NULL;
#endif
	}
	else if (name == "tahoe_THK")
	{
#if defined(BRIDGING_ELEMENT) && defined(BRIDGING_ELEMENT_DEV)
		return new FEManagerT_THK(input_file, output, comm, argv, task);
#else
		ExceptionT::GeneralFail(caller, "\"%s\" requires BRIDGING_ELEMENT and __DEVELOPMENT__",
			name.Pointer());
		return NULL;
#endif
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized name \"%s\"",
			name.Pointer());

	return NULL;
}
