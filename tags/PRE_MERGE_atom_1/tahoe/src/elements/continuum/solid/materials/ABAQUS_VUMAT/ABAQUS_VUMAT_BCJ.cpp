/* $Id: ABAQUS_VUMAT_BCJ.cpp,v 1.8 2002-11-14 17:05:59 paklein Exp $ */
/* created: paklein (05/09/2000) */
#include "ABAQUS_VUMAT_BCJ.h"

#ifdef __F2C__

using namespace Tahoe;

/* function prototype */
extern "C" {
  //int umat_(doublereal *stress, doublereal *statev, doublereal
  //*ddsdde, doublereal *sse, doublereal *spd, doublereal *scd,
  //doublereal *rpl, doublereal *ddsddt, doublereal *drplde, doublereal *
  //drpldt, doublereal *stran, doublereal *dstran, doublereal *time,
  //doublereal *dtime, doublereal *temp, doublereal *dtemp, doublereal *
  //predef, doublereal *dpred, char *cmname, integer *ndi, integer *nshr,
  //integer *ntens, integer *nstatv, doublereal *props, integer *nprops,
  //doublereal *coords, doublereal *drot, doublereal *pnewdt, doublereal *
  //celent, doublereal *dfgrd0, doublereal *dfgrd1, integer *noel,
  //integer *npt, integer *layer, integer *kspt, integer *kstep, integer *
  //kinc, ftnlen cmname_len);
          int vumat_(integer *nblock, integer *ndi, integer *nshr, integer *nstatv,
	  integer *nfieldv, integer *nprops, integer *lanneal, doublereal *steptime, doublereal *totaltime,
	  doublereal *dtime, char *cmname, doublereal *coords, doublereal *celent,
	  doublereal *props, doublereal *density, doublereal *dstran, doublereal *relspininc,
	  doublereal *tempOld, doublereal *stretchold, doublereal *dfgrd0, doublereal *predef,
	  doublereal *stressold, doublereal *statevold, doublereal *enerInternOld, 
	  doublereal *enerInelasOld, doublereal *tempNew, doublereal *stretchnew,
	  doublereal *dfgrd1, doublereal *dpred, doublereal *stressnew, doublereal *statevnew,
	  doublereal *enerInternNew, doublereal *enerInelasNew);
}

#if 0
/* fortran function prototype (append underscode to fortran subroutine name) */
extern "C" {
int cycdmg7_([full argument list]);
}
#endif

/* constructor */
ABAQUS_VUMAT_BCJ::ABAQUS_VUMAT_BCJ(ifstreamT& in, const FDMatSupportT& support):
	ABAQUS_VUMAT_BaseT(in, support)
{

}

/***********************************************************************
* Private
***********************************************************************/

/* VUMAT function wrapper */
//void ABAQUS_VUMAT_BCJ::UMAT(doublereal *stress, doublereal *statev, doublereal
//*ddsdde, doublereal *sse, doublereal *spd, doublereal *scd,
//doublereal *rpl, doublereal *ddsddt, doublereal *drplde, doublereal *
//drpldt, doublereal *stran, doublereal *dstran, doublereal *time,
//doublereal *dtime, doublereal *temp, doublereal *dtemp, doublereal *
//predef, doublereal *dpred, char *cmname, integer *ndi, integer *nshr,
//integer *ntens, integer *nstatv, doublereal *props, integer *nprops,
//doublereal *coords, doublereal *drot, doublereal *pnewdt, doublereal *
//celent, doublereal *dfgrd0, doublereal *dfgrd1, integer *noel,
//integer *npt, integer *layer, integer *kspt, integer *kstep, integer *
//kinc, ftnlen cmname_len)
void ABAQUS_VUMAT_BCJ::VUMAT(integer *nblock, integer *ndi, integer *nshr, integer *nstatv,
	  integer *nfieldv, integer *nprops, integer *lanneal, doublereal *steptime, doublereal *totaltime,
	  doublereal *dtime, char *cmname, doublereal *coords, doublereal *celent,
	  doublereal *props, doublereal *density, doublereal *dstran, doublereal *relspininc,
	  doublereal *tempOld, doublereal *stretchold, doublereal *dfgrd0, doublereal *predef,
	  doublereal *stressold, doublereal *statevold, doublereal *enerInternOld, 
	  doublereal *enerInelasOld, doublereal *tempNew, doublereal *stretchnew,
	  doublereal *dfgrd1, doublereal *dpred, doublereal *stressnew, doublereal *statevnew,
	  doublereal *enerInternNew, doublereal *enerInelasNew)
{
	/* call VUMAT */
  //umat_(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde,
  //	drpldt, stran, dstran, time, dtime, temp, dtemp, predef, dpred,
  //	cmname, ndi, nshr, ntens, nstatv, props, nprops, coords, drot,
  //	pnewdt, celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep,
  //	kinc, cmname_len);
        vumat_(nblock, ndi, nshr, nstatv,
	  nfieldv, nprops, lanneal, steptime, totaltime,
	  dtime, cmname, coords, celent,
	  props, density, dstran, relspininc,
	  tempOld, stretchold, dfgrd0, predef,
	  stressold, statevold, enerInternOld, 
	  enerInelasOld, tempNew, stretchnew,
	  dfgrd1, dpred, stressnew, statevnew,
	  enerInternNew, enerInelasNew);
}

/* set material output */
void ABAQUS_VUMAT_BCJ::SetOutputVariables(iArrayT& variable_index,
	ArrayT<StringT>& output_labels)
{
	int num_output = 5;

	/* number of output */
	variable_index.Dimension(num_output);
	variable_index[0] = 6;
	variable_index[1] = 7;
	variable_index[2] = 11;
	variable_index[3] = 9;
	variable_index[4] = 8;

	/* labels */
	output_labels.Dimension(num_output);
	output_labels[0] = "kappa";
	output_labels[1] = "temp";
	output_labels[2] = "pl_strn";
	output_labels[3] = "damage";
	output_labels[4] = "pl_strn_rate";
}
#endif /* __F2C__ */
