/*****************************************************************************/
/*                                                                           */
/*   (C) Copyright IBM Corporation, 1997                                     */
/*   (C) Copyright Regents of the University of Minnesota, 1997              */
/*                                                                           */
/*   pspaces.h                                                               */
/*                                                                           */
/*   Written by Mahesh Joshi, U of MN.                                       */
/*                                                                           */
/*****************************************************************************/
/*                                                                           */
/* This code is meant to be used solely for educational, research, and       */
/* benchmarking purposes by non-profit institutions and US government        */
/* agencies only.  Use by any other organization requires prior written      */
/* permission from both IBM Corporation and the University of Minnesota.     */
/* The software may not be sold or redistributed.  One may make copies       */
/* of the software or modify it for their use provided that the copies,      */
/* modified or otherwise, are not sold or distributed, are used under the    */
/* same terms and conditions, and this notice and any part of the source     */
/* code that follows this notice are not separated.                          */
/*                                                                           */
/* As unestablished research software, this code is provided on an           */
/* ``as is'' basis without warranty of any kind, either expressed or         */
/* implied, including but not limited to implied warranties of               */
/* merchantability and fitness for a particular purpose.  IBM does not       */
/* warrant that the functions contained in this software will meet the       */
/* user's requirements or that the operation of its routines will be         */
/* uninterrupted or error-free.  Acceptance and use of this program          */
/* constitutes the user's understanding that he/she will have no recourse    */
/* to IBM for any actual or consequential damages, including, but not        */
/* limited to, lost profits or savings, arising out of the use or inability  */
/* to use these libraries.  Even if the user informs IBM of the possibility  */
/* of such damages, IBM expects the user to accept the risk of any such      */
/* harm, or the user shall not attempt to use these libraries for any        */
/* purpose.                                                                  */
/*                                                                           */
/* The downloading, compiling, or executing any part of this software        */
/* constitutes an implicit agreement to these terms.  These terms and        */
/* conditions are subject to change at any time without prior notice.        */
/*                                                                           */
/*****************************************************************************/
/* $Id: pspaces.h,v 1.1.1.1 2004-10-07 16:05:26 paklein Exp $ */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <mpi.h>

/* Following for xPSPACEz functions */

#define O_NOPTS 1
#define Y_NOPTS 3
#define F_NOPTS 4
#define T_NOPTS 1

#define FILLFACTOR 3.0
#define IFOPTS_SIZE 25

#define MAX_OPEN_PSPCOMMS	256

typedef struct ptrst{
  struct ptrst *Ylink;
  int howManyNfacts;
  int *rowdist,*order,*ranmasks;
  int *wrkint,*cinfo,*painds,*paptrs,*supinds;
  int *lptrs,*linds,*tptrs,*tinds,*sup,*tsind;
  int *lc,*iptrs,*ifopts;
  double *lvals;
  MPI_Comm *pmcomm;
} PTRS;

/* Following for ParMETIS */
/* Uncomment the Following #define if ParMetis was built with short. */
/** NOTE: CURRENTLY PSPASES DOES NOT WORK WITH short DATATYPE OF ParMETIS. */
#define IDXTYPE_INT

#ifdef IDXTYPE_INT
typedef int idxtype;
#define IDX_DATATYPE    MPI_INT
#else
typedef short idxtype;
#define IDX_DATATYPE    MPI_SHORT
#endif

struct KeyValueType {
  idxtype key;
  idxtype val;
};

typedef struct KeyValueType KeyValueType;

#define ikeysort ikeysort__

#ifdef UNDERSCORE

#define Fpspaceo pspaceo_ 
#define Fpspacey pspacey_ 
#define Fdpspacen dpspacen_ 
#define Fdpspacef dpspacef_ 
#define Fdpspacet dpspacet_
#define Fpspacec pspacec_
#define Fcheckb_ax checkb_ax_

#define porder porder_
#define emovea emovea_
#define pmovea pmovea_
#define gentree gentree_
#define parsymb parsymb_
#define parfact1 parfact1_
#define compmysan compmysan_
#define trisolve trisolve_
#define initrnd initrnd_
#define mydrand48n mydrand48n_
#define db_ax db_ax_
#define parometisf parometisf_
#define ikeysortf ikeysortf_
#define ygentree ygentree_
#define premovea premovea_
#define moveai moveai_
#define moveav moveav_
#define pparfact1 pparfact1_
#define eparfact1 eparfact1_

#elif defined(CAPS)

#define FCMATCH

#define Fpspaceo PSPACEO 
#define Fpspacey PSPACEY 
#define Fdpspacen DPSPACEN
#define Fdpspacef DPSPACEF
#define Fdpspacet DPSPACET
#define Fpspacec PSPACEC
#define Fcheckb_ax CHECKB_AX

#define porder PORDER
#define emovea EMOVEA
#define pmovea PMOVEA
#define gentree GENTREE
#define parsymb PARSYMB
#define parfact1 PARFACT1
#define compmysan COMPMYSAN
#define trisolve TRISOLVE
#define initrnd INITRND
#define mydrand48n MYDRAND48N
#define db_ax DB_AX
#define parometisf PAROMETISF
#define ikeysortf IKEYSORTF
#define ygentree YGENTREE
#define premovea PREMOVEA
#define moveai MOVEAI
#define moveav MOVEAV
#define pparfact1 PPARFACT1
#define eparfact1 EPARFACT1

#else

#define Fpspaceo pspaceo
#define Fpspacey pspacey
#define Fdpspacen dpspacen
#define Fdpspacef dpspacef
#define Fdpspacet dpspacet
#define Fpspacec pspacec
#define Fcheckb_ax checkb_ax

#endif

/* some macros */
#define max(a,b)	((a) > (b)) ? (a) : (b)
#define min(a,b)	((a) < (b)) ? (a) : (b)

/* prototypes for user callable routines.*/
void Fpspaceo(int *rowdista,int *aptrs,int *ainds,int *order,
              int *sizes,int *options,MPI_Comm *pcomm);

void Fpspacey(int *rowdista,int *aptrs,int *ainds,int *order,int *sizes,
              int *ioptions,double *doptions,long *pspcomm,MPI_Comm *pcomm);

void Fdpspacen(int *rowdista,int *aptrs,int *ainds,double *avals,
              long *pspcommin,long *pspcommout,MPI_Comm *pcomm);

void Fdpspacef(int *rowdista,int *aptrs,int *ainds,double *avals,
               int *ioptions,double *doptions,long *pspcomm,MPI_Comm *pcomm);

void Fdpspacet(int *rowdistb,int *pnrhs,double *b,int *pldb,double *x,
               int *pldx,int *options,long *pspcomm,MPI_Comm *pcomm);

void Fpspacec(long *pspcomm,int *option);

void Fcheckb_ax(int *rowdista,int *aptrs,int *ainds,double *avals,
                int *rowdistb,int *pnrhs,double *b,int *pldb,double *x,
                int *pldx,double *perr,MPI_Comm *pcomm);

