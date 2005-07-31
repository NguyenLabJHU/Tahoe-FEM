/* Include this file in C applications which use PSPASES calls. */

void PSPACEO(int *rowdista,int *aptrs,int *ainds,int *order,
             int *sizes,int *options,MPI_Comm *pcomm);

void PSPACEY(int *rowdista,int *aptrs,int *ainds,int *order,int *sizes,
             int *ioptions,double *doptions,long *pspcomm,MPI_Comm *pcomm);

void DPSPACEN(int *rowdista,int *aptrs,int *ainds,double *avals,
              long *pspcommin,long *pspcommout,MPI_Comm *pcomm);

void DPSPACEF(int *rowdista,int *aptrs,int *ainds,double *avals,
              int *ioptions,double *doptions,long *pspcomm,MPI_Comm *pcomm);

void DPSPACET(int *rowdistbx,int *pnrhs,double *b,int *pldb,double *x,
              int *pldx,int *options,long *pspcomm,MPI_Comm *pcomm);

void PSPACEC(long *pspcomm,int *option);

void CHECKB_AX(int *rowdista,int *aptrs,int *ainds,double *avals,
               int *rowdistb,int *pnrhs,double *b,int *pldb,double *x,
               int *pldx,double *perr,MPI_Comm *pcomm);

