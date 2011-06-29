/* $Id: SurfacePotentialT_factory.cpp,v 1.10 2007-03-14 18:17:50 regueiro Exp $ */
#include "SurfacePotentialT.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentElementsConfig.h"
#endif

#include "XuNeedleman2DT.h"
#include "XuNeedleman3DT.h"
#include "TvergHutch2DT.h"
#include "TvergHutch3DT.h"
#include "TvergHutchIrrev3DT.h"
#include "ViscTvergHutch2DT.h"
#include "Tijssens2DT.h"
#include "RateDep2DT.h"
#include "YoonAllen2DT.h"
#include "YoonAllen3DT.h"
#include "SIMOD_2DT.h"
//#include "LinearDamage2DT.h"

#ifdef COHESIVE_SURFACE_ELEMENT_DEV
#include "InelasticDuctile_RP2DT.h"
#include "MR2DT.h"
#include "MR3DT.h"
#include "MR_RP2DT.h"
#include "MR_NodalRP2DT.h"
#endif

#include <string.h>

using namespace Tahoe;

/* factory method */
SurfacePotentialT* SurfacePotentialT::New(const char* name)
{
	if (strcmp(name, "Xu-Needleman_2D") == 0)
		return new XuNeedleman2DT;
	else if (strcmp(name, "Xu-Needleman_3D") == 0)
		return new XuNeedleman3DT;
	else if (strcmp(name, "Tvergaard-Hutchinson_2D") == 0)
		return new TvergHutch2DT;
	else if (strcmp(name, "Tvergaard-Hutchinson_3D") == 0)
		return new TvergHutch3DT;
	else if (strcmp(name, "Tvergaard-Hutchinson_Irreversible_3D") == 0)
	        return new TvergHutchIrrev3DT;
	else if (strcmp(name, "viscous_Tvergaard-Hutchinson_2D") == 0)
		return new ViscTvergHutch2DT;
	else if (strcmp(name, "Tijssens_2D") == 0)
		return new Tijssens2DT;
	else if (strcmp(name, "Tvergaard-Hutchinson_rate_dep_2D") == 0)
		return new RateDep2DT;
	else if (strcmp(name, "Yoon-Allen_2D") == 0)
		return new YoonAllen2DT;
	else if (strcmp(name, "Yoon-Allen_3D") == 0)
		return new YoonAllen3DT;

#ifdef COHESIVE_SURFACE_ELEMENT_DEV
	else if (strcmp(name, "rigid-inelastic_BCJ_2D") == 0)
		return new InelasticDuctile_RP2DT;
	else if (strcmp(name, "elastoplastic_MR_2D") == 0)
		return new MR2DT;
	else if (strcmp(name, "elastoplastic_MR_3D") == 0)
		return new MR3DT;
	else if (strcmp(name, "rigid-plastic_MR_RP2D") == 0)
		return new MR_RP2DT;
	else if (strcmp(name, "nodal-rigid-plastic_MR_RP2D") == 0)
		return new MR_NodalRP2DT;
#endif

#ifdef __SIMOD__
	else if (strcmp(name, "SIMOD_2D") == 0)
		return new SIMOD_2DT;
#endif

	else
		return NULL;
}
