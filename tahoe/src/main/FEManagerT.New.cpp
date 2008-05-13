/* $Id: FEManagerT.New.cpp,v 1.11 2008-05-13 22:38:04 regueiro Exp $ */
#include "FEManagerT.h"

/* element configuration header */
#include "ElementsConfig.h"

#ifdef BRIDGING_ELEMENT
#include "MultiManagerT.h"
#include "FEManagerT_bridging.h"
#include "BridgingScaleManagerT.h"
#include "FEManagerT_THK.h"
#endif

#ifdef __DEVELOPMENT__
#include "DevelopmentElementsConfig.h"
#include "ThermomechanicalCouplingManagerT.h"
#endif

#ifdef DEM_COUPLING_DEV
#include "FEDEManagerT.h"
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
#if defined(BRIDGING_ELEMENT)
		return new BridgingScaleManagerT(input_file, output, comm, argv, task);
#else
		ExceptionT::GeneralFail(caller, "\"%s\" requires BRIDGING_ELEMENT",
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
#if defined(BRIDGING_ELEMENT)
		return new FEManagerT_THK(input_file, output, comm, argv, task);
#else
		ExceptionT::GeneralFail(caller, "\"%s\" requires BRIDGING_ELEMENT",
			name.Pointer());
		return NULL;
#endif
	}
	else if (name == "tahoe_DEM_coupling")
	{
#if defined(DEM_COUPLING_DEV)
	    return new FEDEManagerT(input_file, output, comm, argv, task);
#else
	    ExceptionT::GeneralFail(caller, "\"%s\" requires DEM_COUPLING_DEV",
				    name.Pointer());
	    return NULL;
#endif

	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized name \"%s\"",
			name.Pointer());

	return NULL;
}
