/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
   utility file nr.h.  Do not confuse this file with the same-named
   file nr.h that is supplied in the 'misc' subdirectory.
   *That* file is the one from the book, and contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only ANSI C.               */

#ifndef _NR_H_
#define _NR_H_

#include <cstdio>
using namespace std;

namespace dem {

#ifndef _FCOMPLEX_DECLARE_T_
typedef struct FCOMPLEX {long double r,i;} fcomplex;
#define _FCOMPLEX_DECLARE_T_
#endif /* _FCOMPLEX_DECLARE_T_ */

#ifndef _ARITHCODE_DECLARE_T_
typedef struct {
	unsigned long *ilob,*iupb,*ncumfq,jdif,nc,minint,nch,ncum,nrad;
} arithcode;
#define _ARITHCODE_DECLARE_T_
#endif /* _ARITHCODE_DECLARE_T_ */

#ifndef _HUFFCODE_DECLARE_T_
typedef struct {
	unsigned long *icod,*ncod,*left,*right,nch,nodemax;
} huffcode;
#define _HUFFCODE_DECLARE_T_
#endif /* _HUFFCODE_DECLARE_T_ */

void addint(long double **uf, long double **uc, long double **res, int nf);
void airy(long double x, long double *ai, long double *bi, long double *aip, long double *bip);
void amebsa(long double **p, long double y[], int ndim, long double pb[],	long double *yb,
	long double ftol, long double (*funk)(long double []), int *iter, long double temptr);
void amoeba(long double **p, long double y[], int ndim, long double ftol,
	long double (*funk)(long double []), int *iter);
long double amotry(long double **p, long double y[], long double psum[], int ndim,
	long double (*funk)(long double []), int ihi, long double fac);
long double amotsa(long double **p, long double y[], long double psum[], int ndim, long double pb[],
	long double *yb, long double (*funk)(long double []), int ihi, long double *yhi, long double fac);
void anneal(long double x[], long double y[], int iorder[], int ncity);
long double anorm2(long double **a, int n);
void arcmak(unsigned long nfreq[], unsigned long nchh, unsigned long nradd,
	arithcode *acode);
void arcode(unsigned long *ich, unsigned char **codep, unsigned long *lcode,
	unsigned long *lcd, int isign, arithcode *acode);
void arcsum(unsigned long iin[], unsigned long iout[], unsigned long ja,
	int nwk, unsigned long nrad, unsigned long nc);
void asolve(unsigned long n, long double b[], long double x[], int itrnsp);
void atimes(unsigned long n, long double x[], long double r[], int itrnsp);
void avevar(long double data[], unsigned long n, long double *ave, long double *var);
void balanc(long double **a, int n);
void banbks(long double **a, unsigned long n, int m1, int m2, long double **al,
	unsigned long indx[], long double b[]);
void bandec(long double **a, unsigned long n, int m1, int m2, long double **al,
	unsigned long indx[], long double *d);
void banmul(long double **a, unsigned long n, int m1, int m2, long double x[], long double b[]);
void bcucof(long double y[], long double y1[], long double y2[], long double y12[], long double d1,
	long double d2, long double **c);
void bcuint(long double y[], long double y1[], long double y2[], long double y12[],
	long double x1l, long double x1u, long double x2l, long double x2u, long double x1,
	long double x2, long double *ansy, long double *ansy1, long double *ansy2);
void beschb(long double x, long double *gam1, long double *gam2, long double *gampl,
	long double *gammi);
long double bessi(int n, long double x);
long double bessi0(long double x);
long double bessi1(long double x);
void bessik(long double x, long double xnu, long double *ri, long double *rk, long double *rip,
	long double *rkp);
long double bessj(int n, long double x);
long double bessj0(long double x);
long double bessj1(long double x);
void bessjy(long double x, long double xnu, long double *rj, long double *ry, long double *rjp,
	long double *ryp);
long double bessk(int n, long double x);
long double bessk0(long double x);
long double bessk1(long double x);
long double bessy(int n, long double x);
long double bessy0(long double x);
long double bessy1(long double x);
long double beta(long double z, long double w);
long double betacf(long double a, long double b, long double x);
long double betai(long double a, long double b, long double x);
long double bico(int n, int k);
void bksub(int ne, int nb, int jf, int k1, int k2, long double ***c);
long double bnldev(long double pp, int n, long *idum);
long double brent(long double ax, long double bx, long double cx,
	long double (*f)(long double), long double tol, long double *xmin);
void broydn(long double x[], int n, int *check,
	void (*vecfunc)(int, long double [], long double []));
void bsstep(long double y[], long double dydx[], int nv, long double *xx, long double htry,
	long double eps, long double yscal[], long double *hdid, long double *hnext,
	void (*derivs)(long double, long double [], long double []));
void caldat(long julian, int *mm, int *id, int *iyyy);
void chder(long double a, long double b, long double c[], long double cder[], int n);
long double chebev(long double a, long double b, long double c[], int m, long double x);
void chebft(long double a, long double b, long double c[], int n, long double (*func)(long double));
void chebpc(long double c[], long double d[], int n);
void chint(long double a, long double b, long double c[], long double cint[], int n);
long double chixy(long double bang);
void choldc(long double **a, int n, long double p[]);
void cholsl(long double **a, int n, long double p[], long double b[], long double x[]);
void chsone(long double bins[], long double ebins[], int nbins, int knstrn,
	long double *df, long double *chsq, long double *prob);
void chstwo(long double bins1[], long double bins2[], int nbins, int knstrn,
	long double *df, long double *chsq, long double *prob);
void cisi(long double x, long double *ci, long double *si);
void cntab1(int **nn, int ni, int nj, long double *chisq,
	long double *df, long double *prob, long double *cramrv, long double *ccc);
void cntab2(int **nn, int ni, int nj, long double *h, long double *hx, long double *hy,
	long double *hygx, long double *hxgy, long double *uygx, long double *uxgy, long double *uxy);
void convlv(long double data[], unsigned long n, long double respns[], unsigned long m,
	int isign, long double ans[]);
void copy(long double **aout, long double **ain, int n);
void correl(long double data1[], long double data2[], unsigned long n, long double ans[]);
void cosft(long double y[], int n, int isign);
void cosft1(long double y[], int n);
void cosft2(long double y[], int n, int isign);
void covsrt(long double **covar, int ma, int ia[], int mfit);
void crank(unsigned long n, long double w[], long double *s);
void cyclic(long double a[], long double b[], long double c[], long double alpha, long double beta,
	long double r[], long double x[], unsigned long n);
void daub4(long double a[], unsigned long n, int isign);
long double dawson(long double x);
long double dbrent(long double ax, long double bx, long double cx,
	long double (*f)(long double), long double (*df)(long double), long double tol, long double *xmin);
void ddpoly(long double c[], int nc, long double x, long double pd[], int nd);
int decchk(char string[], int n, char *ch);
void derivs(long double x, long double y[], long double dydx[]);
long double df1dim(long double x);
void dfour1(long double data[], unsigned long nn, int isign);
void dfpmin(long double p[], int n, long double gtol, int *iter, long double *fret,
	long double (*func)(long double []), void (*dfunc)(long double [], long double []));
long double dfridr(long double (*func)(long double), long double x, long double h, long double *err);
void dftcor(long double w, long double delta, long double a, long double b, long double endpts[],
	long double *corre, long double *corim, long double *corfac);
void dftint(long double (*func)(long double), long double a, long double b, long double w,
	long double *cosint, long double *sinint);
void difeq(int k, int k1, int k2, int jsf, int is1, int isf,
	int indexv[], int ne, long double **s, long double **y);
void dlinmin(long double p[], long double xi[], int n, long double *fret,
	long double (*func)(long double []), void (*dfunc)(long double [], long double[]));
long double dpythag(long double a, long double b);
void drealft(long double data[], unsigned long n, int isign);
void dsprsax(long double sa[], unsigned long ija[], long double x[], long double b[],
	unsigned long n);
void dsprstx(long double sa[], unsigned long ija[], long double x[], long double b[],
	unsigned long n);
void dsvbksb(long double **u, long double w[], long double **v, int m, int n, long double b[],
	long double x[]);
void dsvdcmp(long double **a, int m, int n, long double w[], long double **v);
void eclass(int nf[], int n, int lista[], int listb[], int m);
void eclazz(int nf[], int n, int (*equiv)(int, int));
long double ei(long double x);
void eigsrt(long double d[], long double **v, int n);
long double elle(long double phi, long double ak);
long double ellf(long double phi, long double ak);
long double ellpi(long double phi, long double en, long double ak);
void elmhes(long double **a, int n);
long double erfcc(long double x);
long double erff(long double x);
long double erffc(long double x);
void eulsum(long double *sum, long double term, int jterm, long double wksp[]);
long double evlmem(long double fdt, long double d[], int m, long double xms);
long double expdev(long *idum);
long double expint(int n, long double x);
long double f1(long double x);
long double f1dim(long double x);
long double f2(long double y);
long double f3(long double z);
long double factln(int n);
long double factrl(int n);
void fasper(long double x[], long double y[], unsigned long n, long double ofac, long double hifac,
	long double wk1[], long double wk2[], unsigned long nwk, unsigned long *nout,
	unsigned long *jmax, long double *prob);
void fdjac(int n, long double x[], long double fvec[], long double **df,
	void (*vecfunc)(int, long double [], long double []));
void fgauss(long double x, long double a[], long double *y, long double dyda[], int na);
void fill0(long double **u, int n);
void fit(long double x[], long double y[], int ndata, long double sig[], int mwt,
	long double *a, long double *b, long double *siga, long double *sigb, long double *chi2, long double *q);
void fitexy(long double x[], long double y[], int ndat, long double sigx[], long double sigy[],
	long double *a, long double *b, long double *siga, long double *sigb, long double *chi2, long double *q);
void fixrts(long double d[], int m);
void fleg(long double x, long double pl[], int nl);
void flmoon(int n, int nph, long *jd, long double *frac);
long double fmin(long double x[]);
void four1(long double data[], unsigned long nn, int isign);
void fourew(FILE *file[5], int *na, int *nb, int *nc, int *nd);
void fourfs(FILE *file[5], unsigned long nn[], int ndim, int isign);
void fourn(long double data[], unsigned long nn[], int ndim, int isign);
void fpoly(long double x, long double p[], int np);
void fred2(int n, long double a, long double b, long double t[], long double f[], long double w[],
	long double (*g)(long double), long double (*ak)(long double, long double));
long double fredin(long double x, int n, long double a, long double b, long double t[], long double f[], long double w[],
	long double (*g)(long double), long double (*ak)(long double, long double));
void frenel(long double x, long double *s, long double *c);
void frprmn(long double p[], int n, long double ftol, int *iter, long double *fret,
	long double (*func)(long double []), void (*dfunc)(long double [], long double []));
void ftest(long double data1[], unsigned long n1, long double data2[], unsigned long n2,
	long double *f, long double *prob);
long double gamdev(int ia, long *idum);
long double gammln(long double xx);
long double gammp(long double a, long double x);
long double gammq(long double a, long double x);
long double gasdev(long *idum);
void gaucof(int n, long double a[], long double b[], long double amu0, long double x[], long double w[]);
void gauher(long double x[], long double w[], int n);
void gaujac(long double x[], long double w[], int n, long double alf, long double bet);
void gaulag(long double x[], long double w[], int n, long double alf);
void gauleg(long double x1, long double x2, long double x[], long double w[], int n);
void gaussj(long double **a, int n, long double **b, int m);
void gcf(long double *gammcf, long double a, long double x, long double *gln);
long double golden(long double ax, long double bx, long double cx, long double (*f)(long double), long double tol,
	long double *xmin);
void gser(long double *gamser, long double a, long double x, long double *gln);
void hpsel(unsigned long m, unsigned long n, long double arr[], long double heap[]);
void hpsort(unsigned long n, long double ra[]);
bool hqr(long double **a, int n, long double wr[], long double wi[]);
void hufapp(unsigned long index[], unsigned long nprob[], unsigned long n,
	unsigned long i);
void hufdec(unsigned long *ich, unsigned char *code, unsigned long lcode,
	unsigned long *nb, huffcode *hcode);
void hufenc(unsigned long ich, unsigned char **codep, unsigned long *lcode,
	unsigned long *nb, huffcode *hcode);
void hufmak(unsigned long nfreq[], unsigned long nchin, unsigned long *ilong,
	unsigned long *nlong, huffcode *hcode);
void hunt(long double xx[], unsigned long n, long double x, unsigned long *jlo);
void hypdrv(long double s, long double yy[], long double dyyds[]);
fcomplex hypgeo(fcomplex a, fcomplex b, fcomplex c, fcomplex z);
void hypser(fcomplex a, fcomplex b, fcomplex c, fcomplex z,
	fcomplex *series, fcomplex *deriv);
unsigned short icrc(unsigned short crc, unsigned char *bufptr,
	unsigned long len, short jinit, int jrev);
unsigned short icrc1(unsigned short crc, unsigned char onech);
unsigned long igray(unsigned long n, int is);
void iindexx(unsigned long n, long arr[], unsigned long indx[]);
void indexx(unsigned long n, long double arr[], unsigned long indx[]);
void interp(long double **uf, long double **uc, int nf);
int irbit1(unsigned long *iseed);
int irbit2(unsigned long *iseed);
void jacobi(long double **a, int n, long double d[], long double **v, int *nrot);
void jacobn(long double x, long double y[], long double dfdx[], long double **dfdy, int n);
long julday(int mm, int id, int iyyy);
void kendl1(long double data1[], long double data2[], unsigned long n, long double *tau, long double *z,
	long double *prob);
void kendl2(long double **tab, int i, int j, long double *tau, long double *z, long double *prob);
void kermom(long double w[], long double y, int m);
void ks2d1s(long double x1[], long double y1[], unsigned long n1,
	void (*quadvl)(long double, long double, long double *, long double *, long double *, long double *),
	long double *d1, long double *prob);
void ks2d2s(long double x1[], long double y1[], unsigned long n1, long double x2[], long double y2[],
	unsigned long n2, long double *d, long double *prob);
void ksone(long double data[], unsigned long n, long double (*func)(long double), long double *d,
	long double *prob);
void kstwo(long double data1[], unsigned long n1, long double data2[], unsigned long n2,
	long double *d, long double *prob);
void laguer(fcomplex a[], int m, fcomplex *x, int *its);
void lfit(long double x[], long double y[], long double sig[], int ndat, long double a[], int ia[],
	int ma, long double **covar, long double *chisq, void (*funcs)(long double, long double [], int));
void linbcg(unsigned long n, long double b[], long double x[], int itol, long double tol,
	 int itmax, int *iter, long double *err);
void linmin(long double p[], long double xi[], int n, long double *fret,
	long double (*func)(long double []));
void lnsrch(int n, long double xold[], long double fold, long double g[], long double p[], long double x[],
	 long double *f, long double stpmax, int *check, long double (*func)(long double []));
void load(long double x1, long double v[], long double y[]);
void load1(long double x1, long double v1[], long double y[]);
void load2(long double x2, long double v2[], long double y[]);
void locate(long double xx[], unsigned long n, long double x, unsigned long *j);
void lop(long double **out, long double **u, int n);
void lubksb(long double **a, int n, int *indx, long double b[]);
void ludcmp(long double **a, int n, int *indx, long double *d);
void machar(int *ibeta, int *it, int *irnd, int *ngrd,
	int *machep, int *negep, int *iexp, int *minexp, int *maxexp,
	long double *eps, long double *epsneg, long double *xmin, long double *xmax);
void matadd(long double **a, long double **b, long double **c, int n);
void matsub(long double **a, long double **b, long double **c, int n);
void medfit(long double x[], long double y[], int ndata, long double *a, long double *b, long double *abdev);
void memcof(long double data[], int n, int m, long double *xms, long double d[]);
int metrop(long double de, long double t);
void mgfas(long double **u, int n, int maxcyc);
void mglin(long double **u, int n, int ncycle);
long double midexp(long double (*funk)(long double), long double aa, long double bb, int n);
long double midinf(long double (*funk)(long double), long double aa, long double bb, int n);
long double midpnt(long double (*func)(long double), long double a, long double b, int n);
long double midsql(long double (*funk)(long double), long double aa, long double bb, int n);
long double midsqu(long double (*funk)(long double), long double aa, long double bb, int n);
void miser(long double (*func)(long double []), long double regn[], int ndim, unsigned long npts,
	long double dith, long double *ave, long double *var);
void mmid(long double y[], long double dydx[], int nvar, long double xs, long double htot,
	int nstep, long double yout[], void (*derivs)(long double, long double[], long double[]));
void mnbrak(long double *ax, long double *bx, long double *cx, long double *fa, long double *fb,
	long double *fc, long double (*func)(long double));
void mnewt(int ntrial, long double x[], int n, long double tolx, long double tolf);
void moment(long double data[], int n, long double *ave, long double *adev, long double *sdev,
	long double *var, long double *skew, long double *curt);
void mp2dfr(unsigned char a[], unsigned char s[], int n, int *m);
void mpadd(unsigned char w[], unsigned char u[], unsigned char v[], int n);
void mpdiv(unsigned char q[], unsigned char r[], unsigned char u[],
	unsigned char v[], int n, int m);
void mpinv(unsigned char u[], unsigned char v[], int n, int m);
void mplsh(unsigned char u[], int n);
void mpmov(unsigned char u[], unsigned char v[], int n);
void mpmul(unsigned char w[], unsigned char u[], unsigned char v[], int n,
	int m);
void mpneg(unsigned char u[], int n);
void mppi(int n);
void mprove(long double **a, long double **alud, int n, int indx[], long double b[],
	long double x[]);
void mpsad(unsigned char w[], unsigned char u[], int n, int iv);
void mpsdv(unsigned char w[], unsigned char u[], int n, int iv, int *ir);
void mpsmu(unsigned char w[], unsigned char u[], int n, int iv);
void mpsqrt(unsigned char w[], unsigned char u[], unsigned char v[], int n,
	int m);
void mpsub(int *is, unsigned char w[], unsigned char u[], unsigned char v[],
	int n);
void mrqcof(long double x[], long double y[], long double sig[], int ndata, long double a[],
	int ia[], int ma, long double **alpha, long double beta[], long double *chisq,
	void (*funcs)(long double, long double [], long double *, long double [], int));
void mrqmin(long double x[], long double y[], long double sig[], int ndata, long double a[],
	int ia[], int ma, long double **covar, long double **alpha, long double *chisq,
	void (*funcs)(long double, long double [], long double *, long double [], int), long double *alamda);
void newt(long double x[], int n, int *check,
	void (*vecfunc)(int, long double [], long double []));
void odeint(long double ystart[], int nvar, long double x1, long double x2,
	long double eps, long double h1, long double hmin, int *nok, int *nbad,
	void (*derivs)(long double, long double [], long double []),
	void (*rkqs)(long double [], long double [], int, long double *, long double, long double,
	long double [], long double *, long double *, void (*)(long double, long double [], long double [])));
void orthog(int n, long double anu[], long double alpha[], long double beta[], long double a[],
	long double b[]);
void pade(long double cof[], int n, long double *resid);
void pccheb(long double d[], long double c[], int n);
void pcshft(long double a, long double b, long double d[], int n);
void pearsn(long double x[], long double y[], unsigned long n, long double *r, long double *prob,
	long double *z);
void period(long double x[], long double y[], int n, long double ofac, long double hifac,
	long double px[], long double py[], int np, int *nout, int *jmax, long double *prob);
void piksr2(int n, long double arr[], long double brr[]);
void piksrt(int n, long double arr[]);
void pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k,
	long double ***c, long double **s);
long double plgndr(int l, int m, long double x);
long double poidev(long double xm, long *idum);
void polcoe(long double x[], long double y[], int n, long double cof[]);
void polcof(long double xa[], long double ya[], int n, long double cof[]);
void poldiv(long double u[], int n, long double v[], int nv, long double q[], long double r[]);
void polin2(long double x1a[], long double x2a[], long double **ya, int m, int n,
	long double x1, long double x2, long double *y, long double *dy);
void polint(long double xa[], long double ya[], int n, long double x, long double *y, long double *dy);
void powell(long double p[], long double **xi, int n, long double ftol, int *iter, long double *fret,
	long double (*func)(long double []));
void predic(long double data[], int ndata, long double d[], int m, long double future[], int nfut);
long double probks(long double alam);
void psdes(unsigned long *lword, unsigned long *irword);
void pwt(long double a[], unsigned long n, int isign);
void pwtset(int n);
long double pythag(long double a, long double b);
void pzextr(int iest, long double xest, long double yest[], long double yz[], long double dy[],
	int nv);
long double qgaus(long double (*func)(long double), long double a, long double b);
void qrdcmp(long double **a, int n, long double *c, long double *d, int *sing);
long double qromb(long double (*func)(long double), long double a, long double b);
long double qromo(long double (*func)(long double), long double a, long double b,
	long double (*choose)(long double (*)(long double), long double, long double, int));
void qroot(long double p[], int n, long double *b, long double *c, long double eps);
void qrsolv(long double **a, int n, long double c[], long double d[], long double b[]);
void qrupdt(long double **r, long double **qt, int n, long double u[], long double v[]);
long double qsimp(long double (*func)(long double), long double a, long double b);
long double qtrap(long double (*func)(long double), long double a, long double b);
long double quad3d(long double (*func)(long double, long double, long double), long double x1, long double x2);
void quadct(long double x, long double y, long double xx[], long double yy[], unsigned long nn,
	long double *fa, long double *fb, long double *fc, long double *fd);
void quadmx(long double **a, int n);
void quadvl(long double x, long double y, long double *fa, long double *fb, long double *fc, long double *fd);
long double ran0(long *idum);
long double ran1(long *idum);
long double ran2(long *idum);
long double ran3(long *idum);
long double ran4(long *idum);
void rank(unsigned long n, unsigned long indx[], unsigned long irank[]);
void ranpt(long double pt[], long double regn[], int n);
void ratint(long double xa[], long double ya[], int n, long double x, long double *y, long double *dy);
void ratlsq(long double (*fn)(long double), long double a, long double b, int mm, int kk,
	long double cof[], long double *dev);
long double ratval(long double x, long double cof[], int mm, int kk);
long double rc(long double x, long double y);
long double rd(long double x, long double y, long double z);
void realft(long double data[], unsigned long n, int isign);
void rebin(long double rc, int nd, long double r[], long double xin[], long double xi[]);
void red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
	int ic1, int jc1, int jcf, int kc, long double ***c, long double **s);
void relax(long double **u, long double **rhs, int n);
void relax2(long double **u, long double **rhs, int n);
void resid(long double **res, long double **u, long double **rhs, int n);
long double revcst(long double x[], long double y[], int iorder[], int ncity, int n[]);
void reverse(int iorder[], int ncity, int n[]);
long double rf(long double x, long double y, long double z);
long double rj(long double x, long double y, long double z, long double p);
void rk4(long double y[], long double dydx[], int n, long double x, long double h, long double yout[],
	void (*derivs)(long double, long double [], long double []));
void rkck(long double y[], long double dydx[], int n, long double x, long double h,
	long double yout[], long double yerr[], void (*derivs)(long double, long double [], long double []));
void rkdumb(long double vstart[], int nvar, long double x1, long double x2, int nstep,
	void (*derivs)(long double, long double [], long double []));
void rkqs(long double y[], long double dydx[], int n, long double *x,
	long double htry, long double eps, long double yscal[], long double *hdid, long double *hnext,
	void (*derivs)(long double, long double [], long double []));
void rlft3(long double ***data, long double **speq, unsigned long nn1,
	unsigned long nn2, unsigned long nn3, int isign);
long double rofunc(long double b);
void rotate(long double **r, long double **qt, int n, int i, long double a, long double b);
void rsolv(long double **a, int n, long double d[], long double b[]);
void rstrct(long double **uc, long double **uf, int nc);
long double rtbis(long double (*func)(long double), long double x1, long double x2, long double xacc);
long double rtflsp(long double (*func)(long double), long double x1, long double x2, long double xacc);
long double rtnewt(void (*funcd)(long double, long double *, long double *), long double x1, long double x2,
	long double xacc);
long double rtsafe(void (*funcd)(long double, long double *, long double *), long double x1, long double x2,
	long double xacc);
long double rtsec(long double (*func)(long double), long double x1, long double x2, long double xacc);
void rzextr(int iest, long double xest, long double yest[], long double yz[], long double dy[], int nv);
void savgol(long double c[], int np, int nl, int nr, int ld, int m);
void score(long double xf, long double y[], long double f[]);
void scrsho(long double (*fx)(long double));
long double select(unsigned long k, unsigned long n, long double arr[]);
long double selip(unsigned long k, unsigned long n, long double arr[]);
void shell(unsigned long n, long double a[]);
void shoot(int n, long double v[], long double f[]);
void shootf(int n, long double v[], long double f[]);
void simp1(long double **a, int mm, int ll[], int nll, int iabf, int *kp,
	long double *bmax);
void simp2(long double **a, int m, int n, int *ip, int kp);
void simp3(long double **a, int i1, int k1, int ip, int kp);
void simplx(long double **a, int m, int n, int m1, int m2, int m3, int *icase,
	int izrov[], int iposv[]);
void simpr(long double y[], long double dydx[], long double dfdx[], long double **dfdy,
	int n, long double xs, long double htot, int nstep, long double yout[],
	void (*derivs)(long double, long double [], long double []));
void sinft(long double y[], int n);
void slvsm2(long double **u, long double **rhs);
void slvsml(long double **u, long double **rhs);
void sncndn(long double uu, long double emmc, long double *sn, long double *cn, long double *dn);
long double snrm(unsigned long n, long double sx[], int itol);
void sobseq(int *n, long double x[]);
void solvde(int itmax, long double conv, long double slowc, long double scalv[],
	int indexv[], int ne, int nb, int m, long double **y, long double ***c, long double **s);
void sor(long double **a, long double **b, long double **c, long double **d, long double **e,
	long double **f, long double **u, int jmax, long double rjac);
void sort(unsigned long n, long double arr[]);
void sort2(unsigned long n, long double arr[], long double brr[]);
void sort3(unsigned long n, long double ra[], long double rb[], long double rc[]);
void spctrm(FILE *fp, long double p[], int m, int k, int ovrlap);
void spear(long double data1[], long double data2[], unsigned long n, long double *d, long double *zd,
	long double *probd, long double *rs, long double *probrs);
void sphbes(int n, long double x, long double *sj, long double *sy, long double *sjp, long double *syp);
void splie2(long double x1a[], long double x2a[], long double **ya, int m, int n, long double **y2a);
void splin2(long double x1a[], long double x2a[], long double **ya, long double **y2a, int m, int n,
	long double x1, long double x2, long double *y);
void spline(long double x[], long double y[], int n, long double yp1, long double ypn, long double y2[]);
void splint(long double xa[], long double ya[], long double y2a[], int n, long double x, long double *y);
void spread(long double y, long double yy[], unsigned long n, long double x, int m);
void sprsax(long double sa[], unsigned long ija[], long double x[], long double b[],
	unsigned long n);
void sprsin(long double **a, int n, long double thresh, unsigned long nmax, long double sa[],
	unsigned long ija[]);
void sprspm(long double sa[], unsigned long ija[], long double sb[], unsigned long ijb[],
	long double sc[], unsigned long ijc[]);
void sprstm(long double sa[], unsigned long ija[], long double sb[], unsigned long ijb[],
	long double thresh, unsigned long nmax, long double sc[], unsigned long ijc[]);
void sprstp(long double sa[], unsigned long ija[], long double sb[], unsigned long ijb[]);
void sprstx(long double sa[], unsigned long ija[], long double x[], long double b[],
	unsigned long n);
void stifbs(long double y[], long double dydx[], int nv, long double *xx,
	long double htry, long double eps, long double yscal[], long double *hdid, long double *hnext,
	void (*derivs)(long double, long double [], long double []));
void stiff(long double y[], long double dydx[], int n, long double *x,
	long double htry, long double eps, long double yscal[], long double *hdid, long double *hnext,
	void (*derivs)(long double, long double [], long double []));
void stoerm(long double y[], long double d2y[], int nv, long double xs,
	long double htot, int nstep, long double yout[],
	void (*derivs)(long double, long double [], long double []));
void svbksb(long double **u, long double w[], long double **v, int m, int n, long double b[],
	long double x[]);
void svdcmp(long double **a, int m, int n, long double w[], long double **v);
void svdfit(long double x[], long double y[], long double sig[], int ndata, long double a[],
	int ma, long double **u, long double **v, long double w[], long double *chisq,
	void (*funcs)(long double, long double [], int));
void svdvar(long double **v, int ma, long double w[], long double **cvm);
void toeplz(long double r[], long double x[], long double y[], int n);
void tptest(long double data1[], long double data2[], unsigned long n, long double *t, long double *prob);
void tqli(long double d[], long double e[], int n, long double **z);
long double trapzd(long double (*func)(long double), long double a, long double b, int n);
void tred2(long double **a, int n, long double d[], long double e[]);
void tridag(long double a[], long double b[], long double c[], long double r[], long double u[],
	unsigned long n);
long double trncst(long double x[], long double y[], int iorder[], int ncity, int n[]);
void trnspt(int iorder[], int ncity, int n[]);
void ttest(long double data1[], unsigned long n1, long double data2[], unsigned long n2,
	long double *t, long double *prob);
void tutest(long double data1[], unsigned long n1, long double data2[], unsigned long n2,
	long double *t, long double *prob);
void twofft(long double data1[], long double data2[], long double fft1[], long double fft2[],
	unsigned long n);
void vander(long double x[], long double w[], long double q[], int n);
void vegas(long double regn[], int ndim, long double (*fxn)(long double [], long double), int init,
	unsigned long ncall, int itmx, int nprn, long double *tgral, long double *sd,
	long double *chi2a);
void voltra(int n, int m, long double t0, long double h, long double *t, long double **f,
	long double (*g)(int, long double), long double (*ak)(int, int, long double, long double));
void wt1(long double a[], unsigned long n, int isign,
	void (*wtstep)(long double [], unsigned long, int));
void wtn(long double a[], unsigned long nn[], int ndim, int isign,
	void (*wtstep)(long double [], unsigned long, int));
void wwghts(long double wghts[], int n, long double h,
	void (*kermom)(long double [], long double ,int));
int zbrac(long double (*func)(long double), long double *x1, long double *x2);
void zbrak(long double (*fx)(long double), long double x1, long double x2, int n, long double xb1[],
	long double xb2[], int *nb);
long double zbrent(long double (*func)(long double), long double x1, long double x2, long double tol);
bool zrhqr(long double a[], int m, long double rtr[], long double rti[]);
long double zriddr(long double (*func)(long double), long double x1, long double x2, long double xacc);
void zroots(fcomplex a[], int m, fcomplex roots[], int polish);

} // namespace dem ends

#endif /* _NR_H_ */
