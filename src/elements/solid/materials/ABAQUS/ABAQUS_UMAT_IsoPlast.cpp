/* $Id: ABAQUS_UMAT_IsoPlast.cpp,v 1.1 2005-08-04 07:09:17 paklein Exp $ */
#include "ABAQUS_UMAT_IsoPlast.h"
#include "fortran_names.h"

using namespace Tahoe;

/* function prototype */
extern "C" {
int FORTRAN_NAME(isoplast_umat)(doublereal *stress, doublereal *statev, doublereal
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
ABAQUS_UMAT_IsoPlast::ABAQUS_UMAT_IsoPlast(void):
	ParameterInterfaceT("ABAQUS_UMAT_IsoPlast")
{

}

/***********************************************************************
 * Private
 ***********************************************************************/

/* UMAT function wrapper */
void ABAQUS_UMAT_IsoPlast::UMAT(doublereal *stress, doublereal *statev, doublereal
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
	FORTRAN_NAME(isoplast_umat)(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde,
		drpldt, stran, dstran, time, dtime, temp, dtemp, predef, dpred,
		cmname, ndi, nshr, ntens, nstatv, props, nprops, coords, drot,
		pnewdt, celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep,
		kinc, cmname_len);
}

/* set material output */
void ABAQUS_UMAT_IsoPlast::SetOutputVariables(iArrayT& variable_index,
	ArrayT<StringT>& output_labels)
{
	int num_output = 1;

	/* number of output */
	variable_index.Dimension(num_output);
	variable_index[0] = 2*ntens;

	/* labels */
	output_labels.Dimension(num_output);
	output_labels[0] = "EQPLAS";
}
