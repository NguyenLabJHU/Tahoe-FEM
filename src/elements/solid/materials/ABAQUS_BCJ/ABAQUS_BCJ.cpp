/* $Id: ABAQUS_BCJ.cpp,v 1.1 2003-09-06 07:09:42 paklein Exp $ */
/* created: paklein (05/09/2000) */
#include "ABAQUS_BCJ.h"

#ifdef __F2C__

using namespace Tahoe;

/* function prototype */
extern "C" {
int bcjumat_(doublereal *stress, doublereal *statev, doublereal
	*ddsdde, doublereal *sse, doublereal *spd, doublereal *scd,
	doublereal *rpl, doublereal *ddsddt, doublereal *drplde, doublereal *
	drpldt, doublereal *stran, doublereal *dstran, doublereal *time,
	doublereal *dtime, doublereal *temp, doublereal *dtemp, doublereal *
	predef, doublereal *dpred, char *cmname, integer *ndi, integer *nshr,
	integer *ntens, integer *nstatv, doublereal *props, integer *nprops,
	doublereal *coords, doublereal *drot, doublereal *pnewdt, doublereal *
	celent, doublereal *dfgrd0, doublereal *dfgrd1, integer *noel,
	integer *npt, integer *layer, integer *kspt, integer *kstep, integer *
	kinc, ftnlen cmname_len);
}

/* constructor */
ABAQUS_BCJ::ABAQUS_BCJ(ifstreamT& in, const FSMatSupportT& support):
	ABAQUS_UMAT_BaseT(in, support)
{

}

/***********************************************************************
 * Private
 ***********************************************************************/

/* UMAT function wrapper */
void ABAQUS_BCJ::UMAT(doublereal *stress, doublereal *statev, doublereal
	*ddsdde, doublereal *sse, doublereal *spd, doublereal *scd,
	doublereal *rpl, doublereal *ddsddt, doublereal *drplde, doublereal *
	drpldt, doublereal *stran, doublereal *dstran, doublereal *time,
	doublereal *dtime, doublereal *temp, doublereal *dtemp, doublereal *
	predef, doublereal *dpred, char *cmname, integer *ndi, integer *nshr,
	integer *ntens, integer *nstatv, doublereal *props, integer *nprops,
	doublereal *coords, doublereal *drot, doublereal *pnewdt, doublereal *
	celent, doublereal *dfgrd0, doublereal *dfgrd1, integer *noel,
	integer *npt, integer *layer, integer *kspt, integer *kstep, integer *
	kinc, ftnlen cmname_len)
{
	/* call UMAT */
	bcjumat_(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde,
		drpldt, stran, dstran, time, dtime, temp, dtemp, predef, dpred,
		cmname, ndi, nshr, ntens, nstatv, props, nprops, coords, drot,
		pnewdt, celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep,
		kinc, cmname_len);
}

/* set material output */
void ABAQUS_BCJ::SetOutputVariables(iArrayT& variable_index,
	ArrayT<StringT>& output_labels)
{
	int num_output = 4;

	/* number of output */
#if 0
// these are for cycdmg7
	variable_index.Dimension(num_output);
	variable_index[0] = 6;
	variable_index[1] = 7;
	variable_index[2] = 8;
	variable_index[3] = 13;
#endif

#if 1
// these are for bcjumat
	variable_index.Dimension(num_output);
	variable_index[0] = 6;
	variable_index[1] = 7;
	variable_index[2] = 8;
	variable_index[3] = 9;
#endif

	/* labels */
	output_labels.Dimension(num_output);
	output_labels[0] = "kappa";
	output_labels[1] = "temp";
	output_labels[2] = "pl_strn";
	output_labels[3] = "damage";
}
#endif /* __F2C__ */
