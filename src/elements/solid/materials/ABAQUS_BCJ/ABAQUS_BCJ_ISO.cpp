/* $Id: ABAQUS_BCJ_ISO.cpp,v 1.3 2004-08-01 20:42:35 paklein Exp $ */
#include "ABAQUS_BCJ_ISO.h"

#ifdef __F2C__

using namespace Tahoe;

/* function prototype */
extern "C" {
int bcj_iso_(doublereal *stress, doublereal *statev, doublereal
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
ABAQUS_BCJ_ISO::ABAQUS_BCJ_ISO(void):
	ParameterInterfaceT("ABAQUS_UMAT_BCJ_iso-damage")
{

}

/* accept parameter list */
void ABAQUS_BCJ_ISO::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ABAQUS_UMAT_BaseT::TakeParameterList(list);

	/* set isotropic properties */
	double   Young = double(fProperties[3]);
	double Poisson = double(fProperties[4]);
	Set_E_nu(Young, Poisson);
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* UMAT function wrapper */
void ABAQUS_BCJ_ISO::UMAT(doublereal *stress, doublereal *statev, doublereal
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
	bcj_iso_(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde,
		drpldt, stran, dstran, time, dtime, temp, dtemp, predef, dpred,
		cmname, ndi, nshr, ntens, nstatv, props, nprops, coords, drot,
		pnewdt, celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep,
		kinc, cmname_len);
}

/* set material output */
void ABAQUS_BCJ_ISO::SetOutputVariables(iArrayT& variable_index,
	ArrayT<StringT>& output_labels)
{
	int num_output = 3;

	/* number of output */
	if (NumSD() == 2) {
		variable_index.Dimension(num_output);
		variable_index[0] =  9;
		variable_index[1] =  8;
		variable_index[2] = 10;
	}	
	else if (NumSD() == 3) {
		variable_index.Dimension(num_output);
		variable_index[0] = 13;
		variable_index[1] = 12;
		variable_index[2] = 14;
	}
	else
		ExceptionT::GeneralFail("ABAQUS_BCJ_ISO::SetOutputVariables");

	/* labels */
	output_labels.Dimension(num_output);
	output_labels[0] = "kappa";
	output_labels[1] = "pl_strn";
	output_labels[2] = "damage";
}

#endif /* __F2C__ */
